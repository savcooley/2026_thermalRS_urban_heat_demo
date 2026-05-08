# ARSET: Introduction to Thermal Remote Sensing
# Part 2 Demo 1: Nighttime Cooling Patterns During the LA September 2024 Heatwave
# Author: Savannah Cooley
# ============================================================================
# OVERVIEW
# ============================================================================
# During a multi-day heatwave, not all surfaces cool down overnight equally.
# Green spaces and natural areas release stored heat through evapotranspiration
# and longwave radiation, returning to cooler temperatures by nightfall.
# Impervious urban surfaces (asphalt, concrete) absorb heat during the day
# and re-radiate it slowly, keeping nighttime temperatures dangerously elevated.
#
# This demo uses two ECOSTRESS nighttime acquisitions ~23 hours apart to ask:
#   1. Which surfaces cool most effectively overnight?
#   2. Do the thermal benefits of urban green space persist after a full day
#      of extreme heat accumulation — i.e., does the heat mitigation "hold"
#      as the heatwave intensifies?
#
# Acquisitions used (both nighttime, PDT):
#   Image 1: Sep 3, 11:14 PM PDT (20240904T061421 UTC) — Night 1 of heatwave
#   Image 2: Sep 4, 10:26 PM PDT (20240905T052621 UTC) — Night 2, ~23 hrs later
#
# Study areas (Chatsworth, northwestern LA County):
#   AOI 1: Chatsworth Nature Preserve — natural open space, native vegetation
#   AOI 2: Adjacent commercial/parking area — impervious urban surface

# ============================================================================
# 1. LOAD LIBRARIES
# ============================================================================
library(terra)
library(tidyverse)
library(tidyterra)
library(lubridate)
library(patchwork)   # For combining ggplot panels side-by-side

# ============================================================================
# 2. FILE PATHS AND LOAD DATA
# ============================================================================
repo_dir <- "/Users/sscooley/Documents/GitHub/ARSET/2026/"

# Define the Zenodo URL and the local destination path for the zip file
zip_url <- "https://zenodo.org/records/20090796/files/NASAARSET/2026_thermalRS_urban_heat_demo-init.zip?download=1"
dest_file <- paste0(repo_dir, "2026_thermalRS_urban_heat_demo-init.zip")

# Download the file using curl
download.file(url = zip_url, destfile = dest_file, method = "curl")

# Unzip the downloaded file directly into the repo directory
unzip(zipfile = dest_file, exdir = repo_dir)

# Define the data directory based on the unzipped contents
data_dir <- paste0(repo_dir, "/2026_thermalRS_urban_heat_demo/ECOSTRESS_L2_data")

# Night 1: Sep 3 11:14 PM PDT (23 hrs into heatwave onset)
night1 <- list(
  lst   = file.path(data_dir, "ECO_L2T_LSTE.002_LST_20240904T061421_aid0001_11N.tif"),
  qc    = file.path(data_dir, "ECO_L2T_LSTE.002_QC_20240904T061421_aid0001_11N.tif"),
  cloud = file.path(data_dir, "ECO_L2T_LSTE.002_cloud_20240904T061421_aid0001_11N.tif"),
  label = "Night 1: Sep 3, 11:14 PM PDT",
  date  = "2024-09-03"
)

# Night 2: Sep 4 10:26 PM PDT (after one full day of extreme heat)
night2 <- list(
  lst   = file.path(data_dir, "ECO_L2T_LSTE.002_LST_20240905T052621_aid0001_11N.tif"),
  qc    = file.path(data_dir, "ECO_L2T_LSTE.002_QC_20240905T052621_aid0001_11N.tif"),
  cloud = file.path(data_dir, "ECO_L2T_LSTE.002_cloud_20240905T052621_aid0001_11N.tif"),
  label = "Night 2: Sep 4, 10:26 PM PDT",
  date  = "2024-09-04"
)

# ============================================================================
# 3. PROCESSING FUNCTION
# ============================================================================
# Note on scale factor: The L2T GeoTIFF product stores LST as float32 already
# in Kelvin — no 0.02 scale factor is needed (that applies only to the HDF5
# swath product). We subtract 273.15 to convert K → °C.
# QC bits 0–1 (Mandatory QA) are extracted via modulo 4.
# Cloud mask applied separately — required in Collection 2 (User Guide v4.2).

process_ecostress <- function(lst_path, qc_path, cloud_path) {
  lst   <- rast(lst_path)[[1]]   # Band 1 = LST; Band 2 = LST_err — drop it
  qc    <- rast(qc_path)
  cloud <- rast(cloud_path)
  
  lst_c     <- lst - 273.15                              # K → °C
  qc_bits01 <- qc %% 4                                  # Extract bits 0–1 only
  lst_c     <- mask(lst_c, qc_bits01 >= 2, maskvalues = TRUE)
  lst_c     <- mask(lst_c, cloud > 0,      maskvalues = TRUE)
  
  return(lst_c)
}

# Process both nights
lst_n1 <- process_ecostress(night1$lst, night1$qc, night1$cloud)
lst_n2 <- process_ecostress(night2$lst, night2$qc, night2$cloud)

# Sanity check
message(sprintf("Night 1 LST range: %.1f°C to %.1f°C",
                global(lst_n1, "min", na.rm = TRUE)[[1]],
                global(lst_n1, "max", na.rm = TRUE)[[1]]))
message(sprintf("Night 2 LST range: %.1f°C to %.1f°C",
                global(lst_n2, "min", na.rm = TRUE)[[1]],
                global(lst_n2, "max", na.rm = TRUE)[[1]]))

# ============================================================================
# 4. DEFINE CRS AND AOIs
# ============================================================================
# The L2T GeoTIFF files are delivered in geographic coordinates (lon/lat WGS84)
# despite the "_11N" in the filename — that suffix refers to the UTM zone used
# internally by the ECOSTRESS tiling grid, not the storage CRS of the GeoTIFF.
# We keep everything in WGS84 to avoid any proj.db reprojection dependency.

wgs84_wkt <- 'GEOGCS["WGS 84",DATUM["WGS_1984",
  SPHEROID["WGS 84",6378137,298.257223563]],
  PRIMEM["Greenwich",0],
  UNIT["degree",0.0174532925199433]]'

# Reset raster CRS labels to WGS84 (lon/lat matches actual pixel coordinates)
crs(lst_n1) <- wgs84_wkt
crs(lst_n2) <- wgs84_wkt

# ----------------------------------------------------------------------------
# DEFAULT AOIs — Chatsworth (original demo)
# To switch study areas, comment out this block and uncomment the Hansen Dam
# block below (or substitute your own coordinates).
# ----------------------------------------------------------------------------
aoi_preserve   <- ext(-118.64167, -118.61652, 34.22849, 34.23602) %>%
  vect(crs = wgs84_wkt)
aoi_preserve_name <- "Chatsworth Nature Preserve"

aoi_commercial <- ext(-118.60405, -118.57015, 34.23226, 34.24024) %>%
  vect(crs = wgs84_wkt)
aoi_commercial_name <- "Commercial Area"

study_region   <- ext(-118.65, -118.55, 34.22, 34.25) %>%
  vect(crs = wgs84_wkt)

# ----------------------------------------------------------------------------
# ALTERNATIVE AOIs — Hansen Dam (homework exercise)
# Uncomment this block and comment out the Chatsworth block above to switch.
# ----------------------------------------------------------------------------
# aoi_preserve   <- ext(-118.38833, -118.37166, 34.26197, 34.26767) %>%
#   vect(crs = wgs84_wkt)
# aoi_preserve_name <- "Hansen Dam Wildlife Preserve"
# 
# aoi_commercial <- ext(-118.41, -118.39430, 34.25404, 34.26) %>%
#   vect(crs = wgs84_wkt)
# aoi_commercial_name <- "Commercial Area"
# 
# study_region   <- ext(-118.42, -118.36, 34.25, 34.28) %>%
#   vect(crs = wgs84_wkt)

# Verify — should all print TRUE (geographic = lonlat)
message("Raster is lonlat: ",       is.lonlat(lst_n1))
message("Preserve is lonlat: ",     is.lonlat(aoi_preserve))
message("study_region is lonlat: ", is.lonlat(study_region))

# ============================================================================
# 5. SIDE-BY-SIDE MAP: Stacked panels, shared colorbar, dynamic AOI legend
# ============================================================================

# Crop rasters to study region for faster plotting
lst_n1_crop <- crop(lst_n1, study_region)
lst_n2_crop <- crop(lst_n2, study_region)

# Shared color scale limits across both panels
shared_min <- min(
  global(lst_n1_crop, "min", na.rm = TRUE)[[1]],
  global(lst_n2_crop, "min", na.rm = TRUE)[[1]]
)
shared_max <- max(
  global(lst_n1_crop, "max", na.rm = TRUE)[[1]],
  global(lst_n2_crop, "max", na.rm = TRUE)[[1]]
)

# ----------------------------------------------------------------------------
# Map function
# ----------------------------------------------------------------------------
# AOI outlines are mapped to internal keys ("preserve", "commercial") inside
# aes() so that ggplot recognises them as a discrete color scale. The display
# labels — drawn from aoi_preserve_name and aoi_commercial_name — are applied
# via the `labels` argument of scale_color_manual(). This means the legend
# updates automatically whenever the AOI name variables change, with no edits
# needed inside the function itself.
# ----------------------------------------------------------------------------
make_lst_map <- function(lst_raster, title_label, show_legend = FALSE) {
  ggplot() +
    geom_spatraster(data = lst_raster) +
    # Use fixed internal keys in aes() so ggplot builds the discrete scale;
    # display labels are resolved from the name variables in scale_color_manual
    geom_spatvector(data = aoi_preserve,
                    aes(color = "preserve"),
                    fill = NA, linewidth = 1) +
    geom_spatvector(data = aoi_commercial,
                    aes(color = "commercial"),
                    fill = NA, linewidth = 1) +
    scale_color_manual(
      name   = NULL,
      values = c("preserve"   = "forestgreen",
                 "commercial" = "darkorange"),
      # labels maps internal keys → display strings from the name variables,
      # so switching AOIs only requires changing aoi_preserve_name /
      # aoi_commercial_name in Section 4 above.
      labels = c("preserve"   = aoi_preserve_name,
                 "commercial" = aoi_commercial_name),
      guide  = guide_legend(
        override.aes = list(fill = NA, linewidth = 1.2),
        order = 2
      )
    ) +
    scale_fill_whitebox_c(
      palette  = "muted",
      limits   = c(shared_min, shared_max),
      labels   = function(x) paste0(round(x, 1), "°C"),
      name     = "LST (°C)",
      na.value = "grey85",
      guide    = if (show_legend) guide_colorbar(order = 1) else "none"
    ) +
    labs(title = title_label) +
    theme_minimal() +
    theme(
      plot.title      = element_text(face = "bold", size = 11),
      legend.position = "right"
    )
}

# Only Night 2 (bottom panel) carries the colorbar; Night 1 suppresses it.
# Both panels show the AOI outline legend.
map_n1 <- make_lst_map(lst_n1_crop, night1$label, show_legend = FALSE)
map_n2 <- make_lst_map(lst_n2_crop, night2$label, show_legend = TRUE)

# Stack vertically with patchwork; collect shared legends.
# plot_annotation title uses aoi_preserve_name so it updates automatically.
map_n1 / map_n2 +
  plot_layout(guides = "collect") +
  plot_annotation(
    title   = paste0("Nighttime LST: ", aoi_preserve_name, "During LA Heatwave"),
    caption = "NASA ECOSTRESS L2T LSTE 70m | Collection 2"
  ) &
  theme(legend.position = "right")

# ============================================================================
# 6. EXTRACT PIXEL SAMPLES BY AOI
# ============================================================================
# Extract all valid pixels within each AOI for both nights.
# We build a tidy long-format table for downstream plotting.

extract_aoi <- function(raster, aoi_vect, aoi_name, night_label) {
  vals <- terra::extract(raster, aoi_vect, ID = FALSE)
  colnames(vals) <- "LST_C"
  vals %>%
    drop_na() %>%
    mutate(
      AOI   = aoi_name,
      Night = night_label
    )
}

df_aoi <- bind_rows(
  extract_aoi(lst_n1, aoi_preserve,   "Nature Preserve", night1$label),
  extract_aoi(lst_n1, aoi_commercial, "Commercial Area", night1$label),
  extract_aoi(lst_n2, aoi_preserve,   "Nature Preserve", night2$label),
  extract_aoi(lst_n2, aoi_commercial, "Commercial Area", night2$label)
) %>%
  mutate(
    AOI   = factor(AOI,   levels = c("Nature Preserve", "Commercial Area")),
    Night = factor(Night, levels = c(night1$label, night2$label))
  )

# Summary statistics for annotation
df_summary <- df_aoi %>%
  group_by(AOI, Night) %>%
  summarise(
    mean_lst = mean(LST_C),
    sd_lst   = sd(LST_C),
    n_pixels = n(),
    .groups  = "drop"
  )

print(df_summary)

# ============================================================================
# 7. FIGURE 1: Violin + Boxplot + Jittered raw pixels
# ============================================================================
ggplot(df_aoi, aes(x = AOI, y = LST_C, fill = AOI, color = AOI)) +
  geom_violin(alpha = 0.35, trim = TRUE, color = NA) +
  # Raw pixels underneath the boxplot — jittered to avoid overplotting
  geom_jitter(width = 0.15, size = 0.15, alpha = 0.25, shape = 16) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.85,
               color = "grey20", fill = "white") +
  facet_wrap(~Night) +
  scale_fill_manual(values  = c("Nature Preserve" = "forestgreen",
                                "Commercial Area"  = "darkorange")) +
  scale_color_manual(values = c("Nature Preserve" = "forestgreen",
                                "Commercial Area"  = "darkorange")) +
  labs(
    title    = "Nighttime Surface Temperature Distribution",
    x        = NULL,
    y        = "Land Surface Temperature (°C)",
    caption  = "Each point = one 70m ECOSTRESS pixel | Heatwave: Sep 4–9, 2024"
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "none"
  )

# ============================================================================
# 8. FIGURE 2: ΔT Across Nights — Does the Cooling Benefit Hold?
# ============================================================================
# Core question: after a full day of extreme heat on Sep 4, do we still see
# the same thermal gap between the preserve and the commercial area on Night 2?
#
# We calculate mean LST per AOI per night and plot:
#   (a) The absolute temperature for each AOI across both nights — shows
#       whether each land cover type warmed from Night 1 to Night 2.
#   (b) The thermal gap (commercial − preserve) on each night — if this
#       shrinks, the cooling benefit of the preserve is being eroded by
#       accumulated heat stress; if it persists or grows, the preserve is
#       maintaining its function even under sustained extreme conditions.

# Panel A: Absolute nighttime mean LST by AOI across both nights
p_abs <- ggplot(df_summary, aes(x = Night, y = mean_lst,
                                color = AOI, group = AOI)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_lst - sd_lst,
                    ymax = mean_lst + sd_lst),
                width = 0.1, linewidth = 0.8, alpha = 0.6) +
  scale_color_manual(values = c("Nature Preserve" = "forestgreen",
                                "Commercial Area"  = "darkorange")) +
  labs(
    title  = "A. Nighttime Mean LST",
    x      = NULL,
    y      = "Mean LST (°C)",
    color  = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(face = "bold"),
    axis.text.x     = element_text(size = 8),
    legend.position = "bottom"
  )

# Panel B: Thermal gap (commercial − preserve) per night
df_gap <- df_summary %>%
  select(AOI, Night, mean_lst) %>%
  pivot_wider(names_from = AOI, values_from = mean_lst) %>%
  mutate(gap = `Commercial Area` - `Nature Preserve`)

p_gap <- ggplot(df_gap, aes(x = Night, y = gap, group = 1)) +
  geom_col(fill = "steelblue", alpha = 0.75, width = 0.4) +
  geom_text(aes(label = paste0("+", round(gap, 1), "°C")),
            vjust = -0.5, fontface = "bold", size = 4.5) +
  ylim(0, max(df_gap$gap) * 1.25) +
  labs(
    title    = "B. Thermal Gap",
    subtitle = "Commercial − Nature Preserve",
    x        = NULL,
    y        = "Temperature Difference (°C)",
    caption  = "Larger gap = stronger urban heat island effect"
  ) +
  theme_minimal() +
  theme(
    plot.title  = element_text(face = "bold"),
    axis.text.x = element_text(size = 8)
  )

# Combine panels
(p_abs | p_gap) +
  plot_annotation(
    title    = "Urban Heat Island Persistence During the LA Heatwave",
    subtitle = "Does green space maintain its cooling benefit after a day of extreme heat?",
    caption  = "NASA ECOSTRESS L2T LSTE 70m | Collection 2 | Sep 2024"
  )