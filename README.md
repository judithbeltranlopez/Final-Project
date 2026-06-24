# Distance, Mobility, and Mode Choice: A Spatial Analysis of School Travel Behaviour in Ireland

MSc Health Data Science final project, University of Galway.

**Authors:** Sananthaa Jagadesan Sethamaraikannan, Judith Beltrán López
**Supervisor:** Neil O'Leary


## Overview

This project studies whether school accessibility predicts active travel (walking
and cycling to school) across Census Small Areas (CSAs) in Ireland. We compare
two different measures of accessibility to school places:

1. **Kriging-based density** — a smooth interpolated surface of school-place
   density, based on straight-line distance with a distance-decay
   function over a wide area.
2. **Travel-time accessibility** — the number of school places reachable on foot
   within 30 minutes, based on the actual road network, with a hard 30-minute
   cut-off and no decay.

Both measures are normalised per 100 children aged 3–18, and their effect on the
proportion of active travel is modelled using Bayesian beta regression.


## How to run (pipeline order)

The analysis runs in the following order. Each stage saves outputs that the next
stage reads.

### 1. Data preparation and exploratory analysis — `EDA/`
- `EDA/data_prep.Rmd` — assembles the working dataset.
- `EDA/eda_census.Rmd` — exploratory analysis of the census variables.
- `EDA/eda_schools.Rmd` — exploratory analysis of the schools dataset.

### 2. Census variables and kriging density — `Krigging/`
- `Krigging/krigging_with_finalDF.R` — generates the Census-derived active travel
  proportion and the demographic variables, **and** the kriging density surface.
- `Krigging/Krigging_schools.Rmd` — kriging workflow and diagnostic maps.

### 3. Travel-time accessibility — `Travel Time - Gridpoints/`
- `Travel Time - Gridpoints/R5R Accessibilty Analysis.Rmd` — the current
  travel-time script. Builds a 3×3 grid of points inside each Small Area
  (~160,000 points), computes walking accessibility to school places with `r5r`,
  and aggregates back to one value per Small Area.
  - Network built from the Geofabrik OpenStreetMap extract for Ireland and
    Northern Ireland (`ireland-and-northern-ireland-latest.osm.pbf`).
  - Walking, 30-minute cut-off, step decay (all places within the cut-off count
    equally).
  - *(The `Travel Time - Centroids/` folder contains earlier exploratory
    versions of this analysis and is superseded by the grid-points script above.)*

### 4. Modelling — `model/`
- `model/model.Rmd` — Bayesian beta regression models (one for each accessibility
  measure), with ILR-transformed compositional covariates, spline terms, and full
  convergence diagnostics (Rhat, ESS, divergent transitions). Saves fitted models
  to `model/models/`.

### 5. Comparison and severance analysis — `Comparitive Analysis - Krig Vs Travel Time/`
- `kriging_vs_traveltime.Rmd` — compares the two accessibility measures, finds the
  Small Areas where they diverge most in rank, and investigates these via
  interactive maps with school locations overlaid. Includes the severance analysis
  (areas where schools are close in straight-line distance but not reachable on
  foot because of barriers such as motorways, rivers and water).


## Report and outputs

- **Report:** `Final_Report/FInal_Report/Final_Report.Rmd` (and compiled
  `Final_Report.pdf`). 
- **Poster:** `poster/`(Poster presentation files)
- **Slides:** `slides/`(quarto of the presentation)
- **Saved models and datasets:** `model/models/`, `outputs/`


## Data

Input data lives in `DATA/`. Key inputs:
- `Small_Area_National_Statistical_Boundaries_2022.geojson` — Census Small Area boundaries.
- `schools_sf.rds` — school locations and enrolment.
- Census SAPS files (`SAPS_2022_*.csv`) — demographic variables.
- `gadm41_IRL_shp/` — county and province boundaries.

The OpenStreetMap network file (`ireland-and-northern-ireland-latest.osm.pbf`) and
the r5r build artefacts (`network.dat`, `*.mapdb`) are large and not central to
the analysis; if missing, the `.pbf` can be downloaded from
[Geofabrik](https://download.geofabrik.de/europe/ireland-and-northern-ireland.html)
and rebuilt by re-running the travel-time script.


## Key methods

- **Accessibility:** kriging (gstat) and network travel time (r5r).
- **Compositional covariates:** isometric log-ratio (ILR) transformation.
- **Model:** Bayesian beta regression (brms), weighted by the square root of the
  response denominator, with spline terms for the accessibility measures.

## Notes on filtering

CSAs are excluded where total population > 900 (typically non-residential, e.g.
hospitals, hotels) and where the population aged 3–18 is < 5 (to avoid unstable
normalised values). The number and proportion of CSAs excluded are reported in
the report.

## Overview

This project investigates the relationship between school place accessibility and active travel, specifically walking and cycling, across Ireland. Using spatial analysis methods applied to Census 2022 Small Area data, the project explores whether children living in areas with greater school place accessibility are more likely to walk or cycle to school.(yet to be updated)

## Authors

Judith Beltrán López — @judithbeltranlopez  
Sananthaa Jagadesan Senthamaraikannan — @Sananthaasenthamaraikannan  

MSc Health Data Science, University of Galway, 2025–2026
