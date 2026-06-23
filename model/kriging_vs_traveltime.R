
#comparing the two accessibility measures.

library(tidyverse)
library(sf)
library(leaflet)
library(scales)


# Load the saved data

df_urban_model <- readRDS("D:/Final-Project/model/data2/df_urban_model.rds")
obj <- readRDS("D:/Final-Project/outputs/combined_kriging_objects.rds")

df_urban_model <- st_drop_geometry(df_urban_model)

final_sf   <- obj$final_sf
schools_sf <- obj$dat_irl$schools_agg_sf


#comparison table and rank the two measures

compare_df <- df_urban_model %>%
  filter(
    is.finite(krig_value_per_100_3_18),
    is.finite(access_walk_per_100),
    is.finite(pop_3_18)
  ) %>%
  mutate(
    rank_krig     = rank(krig_value_per_100_3_18),
    rank_walk     = rank(access_walk_per_100),
    rank_diff     = rank_krig - rank_walk,
    abs_rank_diff = abs(rank_diff)
  )

cat("Number of CSAs:", nrow(compare_df), "\n")


#Scatter plot of the two measures
# point size = number of children, so we can see if a big value

p_scatter <- ggplot(compare_df,
                    aes(x = krig_value_per_100_3_18,
                        y = access_walk_per_100,
                        size = pop_3_18,
                        colour = rank_diff)) +
  geom_point(alpha = 0.55) +
  scale_x_continuous(
    limits = c(0, quantile(compare_df$krig_value_per_100_3_18, 0.99))
  ) +
  scale_y_log10(
    labels = label_number(big.mark = ","),
    limits = c(NA, quantile(compare_df$access_walk_per_100, 0.99))
  ) +
  scale_colour_gradient2(
    name = "Kriging - walk\nrank difference",
    low = "#2166AC", mid = "grey85", high = "#B2182B", midpoint = 0
  ) +
  scale_size_continuous(name = "Children\naged 3-18", range = c(0.3, 4)) +
  labs(
    x = "Kriged school-place density (per 100 children aged 3-18)",
    y = "Walking accessibility in 30 min (per 100 children, log scale)",
    title = "Kriging vs travel-time accessibility, urban CSAs"
  ) +
  theme_minimal()

print(p_scatter)


#The 20 CSAs where the two measures disagree the most

divergent <- compare_df %>%
  arrange(desc(abs_rank_diff)) %>%
  slice_head(n = 20) %>%
  select(SA_PUB2022, krig_value_per_100_3_18, access_walk_per_100,
         rank_krig, rank_walk, rank_diff, pop_3_18)

cat("\nTop 20 most divergent CSAs:\n")
print(divergent)

cat("\nHow many of these have fewer than 20 children:",
    sum(divergent$pop_3_18 < 20), "\n")


# Map of those divergent CSAs with schools on top

divergent_sf <- final_sf %>%
  filter(SA_PUB2022 %in% divergent$SA_PUB2022) %>%
  left_join(divergent, by = "SA_PUB2022") %>%
  st_transform(4326)

schools_ll <- st_transform(schools_sf, 4326)
schools_in <- schools_ll[divergent_sf, ]

pal <- colorNumeric("RdBu", domain = divergent_sf$rank_diff)

map_divergent <- leaflet(divergent_sf) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(
    fillColor = ~pal(rank_diff),
    fillOpacity = 0.6, weight = 1, color = "grey30",
    label = ~paste0("CSA ", SA_PUB2022)
  ) %>%
  addCircleMarkers(
    data = schools_in,
    radius = 3, color = "black", fillColor = "yellow",
    fillOpacity = 0.9, stroke = TRUE, weight = 1
  ) %>%
  addLegend(pal = pal, values = divergent_sf$rank_diff,
            title = "Kriging - walk<br>rank difference")

map_divergent


# 6. South-west Dublin anomaly check
# Neil flagged a very high raw kriging value near SW Dublin.
# Look at the 15 highest raw kriging values.

top_raw <- final_sf %>%
  st_drop_geometry() %>%
  filter(is.finite(krig_value_raw)) %>%
  arrange(desc(krig_value_raw)) %>%
  slice_head(n = 15) %>%
  select(SA_PUB2022, krig_value_raw, krig_value_per_100_3_18,
         pop_3_18, n_samples)

cat("\nTop 15 CSAs by raw kriging value:\n")
print(top_raw)

# map those hotspots with schools, to check for stacked schools
hotspot_sf <- final_sf %>%
  filter(SA_PUB2022 %in% top_raw$SA_PUB2022) %>%
  st_transform(4326)

schools_hot <- schools_ll[hotspot_sf, ]

map_hotspot <- leaflet(hotspot_sf) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(fillColor = "red", fillOpacity = 0.4,
              weight = 1, color = "grey30",
              label = ~paste0("CSA ", SA_PUB2022)) %>%
  addCircleMarkers(
    data = schools_hot,
    radius = 4, color = "black", fillColor = "yellow",
    fillOpacity = 0.9, stroke = TRUE, weight = 1
  )

map_hotspot


# 7. Which county are the hotspots in?

counties_sf <- st_read("D:/Final-Project/data/gadm41_IRL_shp/gadm41_IRL_2.shp",
                       quiet = TRUE) %>%
  st_make_valid()

hotspot_named <- final_sf %>%
  filter(SA_PUB2022 %in% top_raw$SA_PUB2022) %>%
  st_transform(st_crs(counties_sf)) %>%
  st_join(select(counties_sf, NAME_1, NAME_2)) %>%
  st_drop_geometry() %>%
  select(SA_PUB2022, NAME_2, krig_value_raw, pop_3_18, n_samples)

print(hotspot_named)

# how high is ~500 compared to the rest of the country?
summary(final_sf$krig_value_raw)
quantile(final_sf$krig_value_raw, c(0.5, 0.9, 0.95, 0.99), na.rm = TRUE)