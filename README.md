
## **ARSET — Introduction to Thermal Remote Sensing and Applications in Urban Heat Island Mapping**

Part 2 | Savannah Cooley, NASA Ames Research Center / Bay Area Environmental Research Institute
# Demo 1: Nighttime Cooling Patterns During the LA September 2024 Heatwave

---

## Overview

This script compares nighttime land surface temperature (LST) between a natural green space and an adjacent impervious urban surface during the September 2024 Los Angeles heatwave, using two ECOSTRESS acquisitions approximately 23 hours apart. The central questions are:

1. Which surface type cools more effectively overnight?
2. Does the thermal benefit of urban green space persist — or erode — after a full day of accumulated extreme heat?

The analysis extracts pixel-level LST from user-defined areas of interest (AOIs), computes summary statistics, and produces two publication-ready figures: a distribution plot and a thermal gap trend.

---

## Prerequisites

### R packages

Install any missing packages before running:

```r
install.packages(c("terra", "tidyverse", "tidyterra", "lubridate", "patchwork"))
```

| Package | Purpose |
|---------|---------|
| `terra` | Raster I/O, masking, pixel extraction |
| `tidyverse` | Data wrangling and ggplot2 plotting |
| `tidyterra` | ggplot2 geoms for terra rasters and vectors (`geom_spatraster`, `geom_spatvector`) |
| `lubridate` | Date handling |
| `patchwork` | Combining ggplot panels |

### R version

Tested on R ≥ 4.2. The `terra` package requires R ≥ 4.1.

---

## Data Requirements

### ECOSTRESS L2T LSTE Collection 2 GeoTIFFs

Three files are required **per acquisition** (six files total):

| Band | Filename pattern | Description |
|------|-----------------|-------------|
| LST | `ECO_L2T_LSTE.002_LST_*_aid0001_11N.tif` | Land surface temperature (float32, Kelvin) |
| QC | `ECO_L2T_LSTE.002_QC_*_aid0001_11N.tif` | Quality control flags |
| Cloud | `ECO_L2T_LSTE.002_cloud_*_aid0001_11N.tif` | Cloud mask |

The script uses two nighttime acquisitions over the LA region:

| Label | UTC timestamp | Local time (PDT) |
|-------|--------------|-----------------|
| Night 1 | 20240904T061421 | Sep 3, 11:14 PM |
| Night 2 | 20240905T052621 | Sep 4, 10:26 PM |

### Downloading the data

Data used in the example area of interest (Los Angeles County, CA) is available here: https://doi.org/10.5281/zenodo.20090795  

Data for other areas of interest can be accessed via **NASA AppEEARS** (Earthdata login required):

1. Go to [appeears.earthdatacloud.nasa.gov](https://appeears.earthdatacloud.nasa.gov)
2. Submit an Area Sample request for product `ECO_L2T_LSTE.002`
3. Select bands: `LST`, `QC`, `cloud`
4. Draw or upload a bounding box covering your study region
5. Set the date range to include your target acquisitions
6. Download the GeoTIFF outputs

> A step-by-step AppEEARS tutorial is available here: https://www.earthdata.nasa.gov/learn/trainings/new-sensor-highlight-ecostress 

### Notes on the L2T GeoTIFF format

- LST values are stored as **float32 in Kelvin** — no 0.02 scale factor is needed (that applies only to the HDF5 swath product). The script subtracts 273.15 to convert to °C.
- Despite `_11N` in the filename (a UTM zone reference used internally by the ECOSTRESS tiling grid), the GeoTIFFs are delivered in **WGS84 geographic coordinates (lon/lat)**. The script explicitly sets this CRS to avoid reprojection issues.
- QC masking uses **bits 0–1** (Mandatory QA) extracted via `qc %% 4`. Pixels with bit value ≥ 2 are rejected.
- Cloud masking is applied separately as required by the Collection 2 User Guide (v4.2): any pixel where `cloud > 0` is masked out.

---

## Setup

### 1. Set the data directory

Edit the `data_dir` path in Section 2 to point to the folder containing your six ECOSTRESS GeoTIFFs:

```r
data_dir <- "/path/to/your/ecostress/data"
```

### 2. Choose your AOIs

The script ships with two pre-configured study area pairs in **Section 4**. Uncomment the block you want to use and comment out the other.

#### Default — Chatsworth (original demo)

```r
aoi_preserve   <- ext(-118.64167, -118.61652, 34.22849, 34.23602) %>%
  vect(crs = wgs84_wkt)
aoi_preserve_name <- "Chatsworth Nature Preserve"

aoi_commercial <- ext(-118.60405, -118.57015, 34.23226, 34.24024) %>%
  vect(crs = wgs84_wkt)
aoi_commercial_name <- "Commercial Area"

study_region   <- ext(-118.65, -118.55, 34.22, 34.25) %>%
  vect(crs = wgs84_wkt)
```

#### Alternative — Hansen Dam (homework exercise)

```r
aoi_preserve   <- ext(-118.38833, -118.37166, 34.26197, 34.26767) %>%
  vect(crs = wgs84_wkt)
aoi_preserve_name <- "Hansen Dam Wildlife Preserve"

aoi_commercial <- ext(-118.40428, -118.39430, 34.25404, 34.25793) %>%
  vect(crs = wgs84_wkt)
aoi_commercial_name <- "Commercial Area"

study_region   <- ext(-118.42, -118.36, 34.25, 34.28) %>%
  vect(crs = wgs84_wkt)
```

#### Using your own AOIs

You can substitute any AOI pair by following the same pattern. Coordinates must be in **WGS84 decimal degrees** in the order `ext(xmin, xmax, ymin, ymax)` (i.e., west, east, south, north). Update `aoi_preserve_name` and `aoi_commercial_name` with descriptive labels — these strings populate figure legends and titles automatically.

Make sure `study_region` is large enough to contain both AOIs with a small buffer.

---

## Script Structure

| Section | Description |
|---------|-------------|
| 1 | Load libraries |
| 2 | Define file paths and acquisition metadata |
| 3 | Load and quality-filter LST rasters (QC + cloud masking, K → °C) |
| 4 | Define CRS and AOI bounding boxes |
| 5 | Produce stacked nighttime LST maps with dynamic AOI outlines |
| 6 | Extract pixel samples by AOI and compute summary statistics (`df_summary`) |
| 7 | Figure 1 — Violin + boxplot + jitter distribution plot by AOI and night |
| 8 | Figure 2 — Absolute mean LST trend and thermal gap (commercial − preserve) |

---

## Outputs

### Console

`df_summary` is printed after Section 6 and contains:

| Column | Description |
|--------|-------------|
| `AOI` | Area of interest label (`Nature Preserve` or `Commercial Area`) |
| `Night` | Acquisition label |
| `mean_lst` | Mean LST in °C across all valid pixels in the AOI |
| `sd_lst` | Standard deviation of LST (°C) |
| `n_pixels` | Number of valid (unmasked) 70 m pixels in the AOI |

### Figures

**Section 5 — Map panel**
Two stacked ECOSTRESS LST maps (Night 1 above, Night 2 below) showing the study region with AOI bounding boxes outlined in green (preserve) and orange (commercial). A shared colorbar appears on the right of the bottom panel. Legend labels and the plot title update automatically from `aoi_preserve_name`.

**Figure 1 (Section 7) — LST distribution**
Violin plot with overlaid boxplot and jittered raw pixels, faceted by night. Shows the full spread of LST values across all valid ECOSTRESS pixels within each AOI, enabling visual comparison of both central tendency and variability.

**Figure 2 (Section 8) — Thermal persistence**
Two-panel figure:
- *Panel A* — Line plot of mean LST (± 1 SD) per AOI across both nights, showing whether each surface type warmed from Night 1 to Night 2.
- *Panel B* — Bar chart of the thermal gap (commercial − preserve mean LST) per night, with the temperature difference annotated. A persistent or growing gap indicates the green space is maintaining its cooling function under sustained heat stress.

---

## Key Concepts

**Why nighttime?**
Daytime LST reflects both absorbed solar energy and surface emissivity, making it difficult to isolate land cover effects. At night, solar input ceases and the thermal signal reflects how well each surface type has shed its stored heat — making surface cover contrasts cleaner and more interpretable.

**Why two nights?**
A single nighttime snapshot could reflect transient conditions. Comparing Night 1 (early in the heatwave) with Night 2 (after a full day of extreme heat accumulation) tests whether the cooling benefit of green space is robust or degrades as the event intensifies.

**The thermal gap**
The difference in mean nighttime LST between the commercial area and the natural preserve (commercial − preserve) quantifies the local urban heat island effect. If this gap narrows on Night 2, accumulated heat stress is eroding the preserve's buffering capacity. If it holds or widens, the green space is functioning as a thermal refuge even under sustained extreme conditions.

---

## Troubleshooting

**`Error: [rast] file does not exist`**
Check that `data_dir` points to the correct folder and that all six GeoTIFF filenames match the patterns in Section 2 exactly (including the UTC timestamp strings).

**All LST values masked to NA after QC filtering**
Verify that the QC and cloud files correspond to the same acquisition as the LST file. Mismatched timestamps will cause alignment failures. Also confirm that the GeoTIFFs were downloaded as Collection 2 (`ECO_L2T_LSTE.002`) — Collection 1 files use a different QC bit structure.

**`df_summary` shows 0 pixels for one or both AOIs**
The AOI extent may fall outside the ECOSTRESS swath footprint for these dates, or all pixels within the AOI may have been masked by cloud or QC flags. Check the map output from Section 5 to confirm the AOI bounding boxes overlap with valid (non-grey) pixels.

**Legend labels still show `"preserve"` / `"commercial"` instead of place names**
Make sure `aoi_preserve_name` and `aoi_commercial_name` are defined as character strings *before* Section 5 runs, and that neither is `NULL` or an empty string.

---

## Citation

NASA ECOSTRESS data accessed via NASA AppEEARS.

*ECOSTRESS L2T Land Surface Temperature and Emissivity, Collection 2.*
NASA EOSDIS Land Processes Distributed Active Archive Center (LP DAAC), USGS Earth Resources Observation and Science (EROS) Center, Sioux Falls, South Dakota.

Script author: Savannah Cooley, NASA Ames Research Center / Bay Area Environmental Research Institute.
