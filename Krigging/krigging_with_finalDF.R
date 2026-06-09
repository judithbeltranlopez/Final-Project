library(tidyverse)
library(sf)
library(sp)
library(gstat)
library(readxl)
library(readr)
library(stringr)
library(janitor)

# -----------------------------
# 1) Paths and parameters
# -----------------------------
school_excel_path <- "data/Data_on_Individual_Schools_Mainstream_2024_25.xlsx"
school_sheet      <- 2
counties_shp_path <- "data/gadm41_IRL_shp/gadm41_IRL_2.shp"
small_area_path   <- "data/Small_Area_National_Statistical_Boundaries_2022.geojson"
saps_path         <- "data/SAPS_2022_Small_Area_UR_171024.csv"
output_dir        <- "outputs"

target_epsg <- 2157

# Kriging parameters
variogram_width  <- 5000
variogram_cutoff <- 200000
vgm_start <- gstat::vgm(
  psill  = 0.8,
  model  = "Exp",
  range  = 80000,
  nugget = 0.1
)

nmax_kriging <- 150
n_side_sa    <- 3

# -----------------------------
# 2) Helper functions
# -----------------------------
to_num <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  readr::parse_number(as.character(x))
}

safe_divide <- function(num, den) {
  if_else(is.finite(den) & den > 0, num / den, NA_real_)
}

check_required_cols <- function(data, cols, data_name = "data") {
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in ", data_name, ": ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
}

choose_sa_id <- function(sa) {
  id_candidates <- c(
    "SA_PUB2022", "GEOGID", "SA_GUID", "SAID", "SA_ID",
    "SMALL_AREA", "OBJECTID", "FID"
  )
  
  sa_id <- id_candidates[id_candidates %in% names(sa)][1]
  
  if (is.na(sa_id)) {
    sa <- sa %>% mutate(sa_uid = row_number())
    sa_id <- "sa_uid"
    warning(
      "No standard Small Area ID column found. Created sa_uid, but SAPS joins may fail."
    )
  }
  
  list(sa = sa, sa_id = sa_id)
}

# -----------------------------
# 3) Prepare school kriging data
# -----------------------------
prepare_school_kriging_data <- function(
    excel_path,
    sheet = 2,
    shp_path,
    region_name = NULL,
    region_field = "NAME_1",
    target_epsg = 2157,
    min_points = 50
) {
  
  schools_df <- read_excel(excel_path, sheet = sheet) %>%
    rename(
      longitude = `School Longitude`,
      latitude  = `School Latitude`
    ) %>%
    filter(!is.na(latitude), !is.na(longitude)) %>%
    mutate(
      enrolment = suppressWarnings(as.numeric(`Enrolment per Return`)),
      log_enrol = log(enrolment)
    ) %>%
    filter(is.finite(enrolment), enrolment > 0, is.finite(log_enrol))
  
  ireland_counties <- sf::st_read(shp_path, quiet = TRUE)
  
  if (is.null(region_name)) {
    region_sf <- ireland_counties %>%
      st_transform(target_epsg) %>%
      st_make_valid() %>%
      st_union() %>%
      st_as_sf() %>%
      mutate(id = 1)
  } else {
    region_sf <- ireland_counties %>%
      filter(.data[[region_field]] == region_name) %>%
      st_transform(target_epsg) %>%
      st_make_valid() %>%
      st_union() %>%
      st_as_sf() %>%
      mutate(id = 1)
  }
  
  schools_sf <- st_as_sf(
    schools_df,
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  ) %>%
    st_transform(target_epsg)
  
  schools_region_sf <- schools_sf[region_sf, ] %>%
    filter(is.finite(log_enrol))
  
  schools_sp <- as(schools_region_sf, "Spatial")
  
  df_sp <- as.data.frame(schools_sp)
  xy <- coordinates(schools_sp)
  df_sp$x <- xy[, 1]
  df_sp$y <- xy[, 2]
  
  schools_agg_df <- df_sp %>%
    group_by(x, y) %>%
    summarise(log_enrol = mean(log_enrol, na.rm = TRUE), .groups = "drop")
  
  coordinates(schools_agg_df) <- ~ x + y
  proj4string(schools_agg_df) <- proj4string(schools_sp)
  
  if (nrow(as.data.frame(schools_agg_df)) <= min_points) {
    stop("Not enough points after filtering/aggregation. Increase area or lower min_points.")
  }
  
  list(
    region_sf      = region_sf,
    schools_agg_sp = schools_agg_df,
    schools_agg_sf = st_as_sf(schools_agg_df)
  )
}

# -----------------------------
# 4) Create 3 x 3 prediction points inside each Small Area
# -----------------------------
make_sa_3x3_points <- function(sa, sa_id, n_side = 3) {
  sa <- st_make_valid(sa)
  
  pts_all <- purrr::map_dfr(seq_len(nrow(sa)), function(i) {
    poly <- sa[i, ]
    
    g <- st_make_grid(
      poly,
      n = c(n_side, n_side),
      what = "polygons"
    )
    
    g <- st_sf(
      local_cell_id = seq_along(g),
      geometry = g,
      crs = st_crs(sa)
    )
    
    id_cols <- unique(c(sa_id, "SA_PUB2022"))
    
    clipped <- suppressWarnings(
      st_intersection(
        g,
        poly %>% select(all_of(id_cols))
      )
    )
    
    clipped <- clipped %>%
      st_make_valid() %>%
      filter(!st_is_empty(geometry)) %>%
      mutate(area_w = as.numeric(st_area(geometry))) %>%
      filter(area_w > 0)
    
    if (nrow(clipped) == 0) return(NULL)
    
    pts <- st_point_on_surface(clipped)
    
    pts %>%
      mutate(
        sa_sample_id = paste0(.data[[sa_id]], "_", local_cell_id)
      )
  })
  
  pts_all %>%
    mutate(sample_row = row_number())
}

# -----------------------------
# 5) Krige at Small-Area 3 x 3 points and aggregate to one value per Small Area
# -----------------------------
run_sa_3x3_kriging <- function(
    schools_agg_sp,
    prediction_points_sf,
    prediction_points_sp,
    sa,
    sa_id,
    width,
    cutoff,
    vgm_start,
    nmax = 150,
    nmin = 1,
    maxdist = 1e6
) {
  v <- variogram(
    log_enrol ~ 1,
    data = schools_agg_sp,
    width = width,
    cutoff = cutoff
  )
  
  fv <- fit.variogram(v, vgm_start)
  
  kriged_sp <- krige(
    log_enrol ~ 1,
    locations = schools_agg_sp,
    newdata   = prediction_points_sp,
    model     = fv,
    nmax      = nmax,
    nmin      = nmin,
    maxdist   = maxdist
  )
  
  kriged_sf <- st_as_sf(kriged_sp) %>%
    mutate(
      !!sa_id := prediction_points_sf[[sa_id]],
      SA_PUB2022 = prediction_points_sf$SA_PUB2022,
      local_cell_id = prediction_points_sf$local_cell_id,
      sa_sample_id = prediction_points_sf$sa_sample_id,
      area_w = prediction_points_sf$area_w,
      enrol_pred = exp(var1.pred)
    )
  
  join_cols <- unique(c(sa_id, "SA_PUB2022"))
  
  sa_agg <- kriged_sf %>%
    st_drop_geometry() %>%
    group_by(across(all_of(join_cols))) %>%
    summarise(
      krig_value = mean(enrol_pred,na.rm = TRUE),
      krig_var   = mean(var1.var,na.rm = TRUE),
      n_samples  = n(),
      .groups = "drop"
    )
  
  sa_with_krig <- sa %>%
    left_join(sa_agg, by = join_cols)
  
  list(
    variogram = v,
    variogram_model = fv,
    kriged_points_sf = kriged_sf,
    sa_with_krig = sa_with_krig
  )
}

# -----------------------------
# 6) Prepare SAPS variables to add to final dataset
# -----------------------------
prepare_saps_variables <- function(saps_path) {
  stopifnot(file.exists(saps_path))
  
  saps_raw <- read_csv(saps_path, show_col_types = FALSE) %>%
    mutate(SA_PUB2022 = str_trim(as.character(GEOGID)))
  
  target_ages <- 3:18
  target_age_cols <- paste0("T1_1AGE", target_ages, "T")
  
  age_demo_cols <- c(
    "T1_1AGE0T", "T1_1AGE1T", "T1_1AGE2T", "T1_1AGE3T", "T1_1AGE4T",
    "T1_1AGE5T", "T1_1AGE6T", "T1_1AGE7T", "T1_1AGE8T", "T1_1AGE9T",
    "T1_1AGE10T", "T1_1AGE11T", "T1_1AGE12T", "T1_1AGE13T", "T1_1AGE14T",
    "T1_1AGE15T", "T1_1AGE16T", "T1_1AGE17T", "T1_1AGE18T", "T1_1AGE19T",
    "T1_1AGE20_24T", "T1_1AGE25_29T", "T1_1AGE30_34T", "T1_1AGE35_39T",
    "T1_1AGE40_44T", "T1_1AGE45_49T", "T1_1AGE50_54T", "T1_1AGE55_59T",
    "T1_1AGE60_64T", "T1_1AGE65_69T", "T1_1AGE70_74T", "T1_1AGE75_79T",
    "T1_1AGE80_84T", "T1_1AGEGE_85T"
  )
  
  car_cols <- c("T6_5_NCH", "T6_5_OCH", "T6_5_NGCH", "T6_5_ECH", "T6_5_CCH")
  
  edu_cols <- c(
    "T10_4_TT", "T10_4_NFT", "T10_4_PT", "T10_4_LST", "T10_4_UST",
    "T10_4_TVT", "T10_4_HCT", "T10_4_HDPQT", "T10_4_PDT", "T10_4_DT",
    "T10_2_SAST"
  )
  
  active_cols <- c(
    "T11_1_TSCCC", "T11_1_WMFHSCCC", "T11_1_NSSCCC",
    "T11_1_FSCCC", "T11_1_BISCCC"
  )
  
  required_cols <- c(
    "GEOGID", "UR_Category_Desc", "T1_1AGETT",
    target_age_cols, age_demo_cols, car_cols, edu_cols, active_cols
  ) %>%
    unique()
  
  check_required_cols(saps_raw, required_cols, "SAPS")
  
  numeric_cols <- c(
    "T1_1AGETT", target_age_cols, age_demo_cols, car_cols, edu_cols, active_cols
  ) %>%
    unique()
  
  saps <- saps_raw %>%
    mutate(across(all_of(numeric_cols), to_num))
  
  saps_vars <- saps %>%
    mutate(
      category = str_extract(UR_Category_Desc, "(?<=\\.).*"),
      category = str_trim(category),
      urban_rural = case_when(
        category %in% c(
          "Cities",
          "Satellite Urban Towns",
          "Independent urban towns"
        ) ~ "Urban",
        category %in% c(
          "Rural areas with high urban influence",
          "Rural areas with moderate urban influence",
          "Highly rural/remote areas"
        ) ~ "Rural",
        TRUE ~ NA_character_
      ),
      pop_3_18 = rowSums(across(all_of(target_age_cols)), na.rm = TRUE),
      denom_active = T11_1_TSCCC - T11_1_WMFHSCCC - T11_1_NSSCCC,
      perc_active = if_else(
        is.finite(denom_active) & denom_active > 0,
        100 * (T11_1_FSCCC + T11_1_BISCCC) / denom_active,
        NA_real_
      ),
      total_pop = T1_1AGETT,
      pct_0_12 = safe_divide(
        100 * rowSums(across(c(
          T1_1AGE0T, T1_1AGE1T, T1_1AGE2T, T1_1AGE3T, T1_1AGE4T,
          T1_1AGE5T, T1_1AGE6T, T1_1AGE7T, T1_1AGE8T, T1_1AGE9T,
          T1_1AGE10T, T1_1AGE11T, T1_1AGE12T
        )), na.rm = TRUE),
        total_pop
      ),
      pct_13_21 = safe_divide(
        100 * rowSums(across(c(
          T1_1AGE13T, T1_1AGE14T, T1_1AGE15T, T1_1AGE16T,
          T1_1AGE17T, T1_1AGE18T, T1_1AGE19T, T1_1AGE20_24T
        )), na.rm = TRUE),
        total_pop
      ),
      pct_22_40 = safe_divide(
        100 * rowSums(across(c(
          T1_1AGE25_29T, T1_1AGE30_34T, T1_1AGE35_39T, T1_1AGE40_44T
        )), na.rm = TRUE),
        total_pop
      ),
      pct_41_70 = safe_divide(
        100 * rowSums(across(c(
          T1_1AGE45_49T, T1_1AGE50_54T, T1_1AGE55_59T,
          T1_1AGE60_64T, T1_1AGE65_69T, T1_1AGE70_74T
        )), na.rm = TRUE),
        total_pop
      ),
      pct_71plus = safe_divide(
        100 * rowSums(across(c(
          T1_1AGE75_79T, T1_1AGE80_84T, T1_1AGEGE_85T
        )), na.rm = TRUE),
        total_pop
      ),
      total_hh = T6_5_NCH + T6_5_OCH + T6_5_NGCH + T6_5_ECH + T6_5_CCH,
      pct_0car = safe_divide(100 * T6_5_NCH, total_hh),
      pct_1car = safe_divide(100 * T6_5_OCH, total_hh),
      pct_2plus = safe_divide(100 * (T6_5_NGCH + T6_5_ECH + T6_5_CCH), total_hh),
      total_edu = T10_4_TT,
      pct_none_primary = safe_divide(100 * (T10_4_NFT + T10_4_PT), total_edu),
      pct_secondary = safe_divide(100 * (T10_4_LST + T10_4_UST + T10_4_TVT), total_edu),
      pct_tertiary = safe_divide(100 * (T10_4_HCT + T10_4_HDPQT + T10_4_PDT + T10_4_DT), total_edu),
      pct_still_in_edu = safe_divide(100 * T10_2_SAST, total_pop)
    ) %>%
    filter(total_pop <= 950) %>%
    select(
      SA_PUB2022,
      category,
      urban_rural,
      pop_3_18,
      perc_active,
      total_pop,
      pct_0_12,
      pct_13_21,
      pct_22_40,
      pct_41_70,
      pct_71plus,
      total_hh,
      pct_0car,
      pct_1car,
      pct_2plus,
      total_edu,
      pct_none_primary,
      pct_secondary,
      pct_tertiary,
      pct_still_in_edu
    ) %>%
    distinct(SA_PUB2022, .keep_all = TRUE)
  
  saps_vars
}

# -----------------------------
# 7) Run everything
# -----------------------------
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 7.1 School data
dat_irl <- prepare_school_kriging_data(
  excel_path  = school_excel_path,
  sheet       = school_sheet,
  shp_path    = counties_shp_path,
  region_name = NULL,
  target_epsg = target_epsg,
  min_points  = 50
)

# 7.2 Small Areas
sa <- st_read(small_area_path, quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(st_crs(dat_irl$region_sf))

sa_id_info <- choose_sa_id(sa)
sa <- sa_id_info$sa
sa_id <- sa_id_info$sa_id

sa <- sa %>%
  mutate(SA_PUB2022 = str_trim(as.character(.data[[sa_id]])))

# 7.3 Prediction points: 3 x 3 inside each Small Area
sa_3x3_pts_sf <- make_sa_3x3_points(
  sa = sa,
  sa_id = sa_id,
  n_side = n_side_sa
)

sa_3x3_pts_sp <- as(sa_3x3_pts_sf, "Spatial")

# 7.4 Kriging directly at those Small-Area points
res_sa <- run_sa_3x3_kriging(
  schools_agg_sp       = dat_irl$schools_agg_sp,
  prediction_points_sf = sa_3x3_pts_sf,
  prediction_points_sp = sa_3x3_pts_sp,
  sa                   = sa,
  sa_id                = sa_id,
  width                = variogram_width,
  cutoff               = variogram_cutoff,
  vgm_start            = vgm_start,
  nmax                 = nmax_kriging
)

# 7.5 SAPS variables
saps_vars <- prepare_saps_variables(saps_path)

# 7.6 Final Small-Area dataset with kriging + other variables
final_sf <- res_sa$sa_with_krig %>%
  mutate(SA_PUB2022 = str_trim(as.character(SA_PUB2022))) %>%
  left_join(saps_vars, by = "SA_PUB2022") %>%
  st_transform(target_epsg) %>%
  mutate(
    area_km2 = as.numeric(st_area(geometry)) / 1e6,
    
    # Raw kriging value kept for reference
    krig_value_raw = krig_value,
    
    # Normalised kriging value by population aged 3 to 18
    krig_value_per_child_3_18 = safe_divide(krig_value_raw, pop_3_18),
    krig_value_per_1000_3_18 = 1000 * krig_value_per_child_3_18,
    
    # If you want krig_value itself to be the normalised variable:
    krig_value = krig_value_per_1000_3_18,
    
    # Optional density variables
    places_density_km2 = safe_divide(krig_value_raw, area_km2),
    
    # Standardised versions
    krig_value_z = as.numeric(scale(krig_value)),
    krig_value_raw_z = as.numeric(scale(krig_value_raw))
  )

final_df <- final_sf %>%
  st_drop_geometry() %>%
  filter(
    !is.na(urban_rural),
    !is.na(perc_active),
    is.finite(krig_value),
    krig_value > 0,
    is.finite(pop_3_18),
    pop_3_18 > 0
  ) %>%
  select(
    SA_PUB2022,
    category,
    urban_rural,
    perc_active,
    
    # Normalised kriging value: predicted value per 1000 children aged 3-18
    krig_value,
    
    # Raw kriging value before normalisation
    krig_value_raw,
    
    krig_var,
    n_samples,
    area_km2,
    pop_3_18,
    places_density_km2,
    krig_value_per_child_3_18,
    krig_value_per_1000_3_18,
    krig_value_z,
    krig_value_raw_z,
    total_pop,
    pct_0_12,
    pct_13_21,
    pct_22_40,
    pct_41_70,
    pct_71plus,
    total_hh,
    pct_0car,
    pct_1car,
    pct_2plus,
    total_edu,
    pct_none_primary,
    pct_secondary,
    pct_tertiary,
    pct_still_in_edu
  )

final_df_urban <- final_df %>%
  filter(urban_rural == "Urban")

# -----------------------------
# 8) Save outputs
# -----------------------------
saveRDS(final_df, file.path(output_dir, "final_df.rds"))
write_csv(final_df, file.path(output_dir, "final_df.csv"))

saveRDS(final_df_urban, file.path(output_dir, "final_df_urban.rds"))
write_csv(final_df_urban, file.path(output_dir, "final_df_urban.csv"))

saveRDS(final_sf, file.path(output_dir, "final_sf_with_geometry.rds"))

saveRDS(
  list(
    dat_irl = dat_irl,
    sa_id = sa_id,
    variogram = res_sa$variogram,
    variogram_model = res_sa$variogram_model,
    kriged_points_sf = res_sa$kriged_points_sf,
    sa_with_krig = res_sa$sa_with_krig,
    final_sf = final_sf,
    final_df = final_df,
    final_df_urban = final_df_urban
  ),
  file.path(output_dir, "combined_kriging_objects.rds")
)

# -----------------------------
# 9) Checks
# -----------------------------
cat("\nSaved outputs in:", output_dir, "\n")
cat("Final dataset rows:", nrow(final_df), "\n")
cat("Urban-only rows:", nrow(final_df_urban), "\n")

print(
  final_df %>%
    summarise(
      n = n(),
      missing_krig = sum(is.na(krig_value)),
      missing_perc_active = sum(is.na(perc_active)),
      min_krig = min(krig_value, na.rm = TRUE),
      max_krig = max(krig_value, na.rm = TRUE),
      min_krig_norm = min(krig_density_per_1000_3_18, na.rm = TRUE),
      max_krig_norm = max(krig_density_per_1000_3_18, na.rm = TRUE)
    )
)















