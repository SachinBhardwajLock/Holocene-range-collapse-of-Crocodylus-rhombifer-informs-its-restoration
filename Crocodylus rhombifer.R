##############################################################################
#  SPECIES DISTRIBUTION MODEL — Crocodylus rhombifer (Cuban Crocodile)
#
#  Climate data:
#    Present  — WorldClim 2.1 at 2.5 arc-min (~5 km), VIF-selected
#    Past     — CHELSA-TraCE21k via pastclim (7–4 kya cell-wise median)
#    Future   — CMIP6 (MIROC6 SSP245 & SSP585, 2061–2080)
#
#  Topography: GMTED-derived 5 km layers (VIF-selected from full set)
#  Geology:    Caribbean lithological map (geo6bg.shp), dummy-encoded
#
#  Occurrence data:
#    Present  — CSV with climate.period, latitude, longitude
#    Past     — CSV with climate.period, latitude, longitude
#
#  Study region: West Indies — Cuba east to Puerto Rico
#    (Bahamas, Turks & Caicos, Cayman Is, Cuba, Jamaica,
#     Hispaniola, Puerto Rico and surrounding islands)
#    Excludes: Florida, mainland Central & South America
#
#  Models fitted:
#    1. Modern-only   (present occurrences × present climate)
#    2. Past-only     (past occurrences × TraCE21k 7–4 kya median paleoclimate)
#    3. Combined      (pooled occurrences × time-matched climate)
#       Combined model uses Past PB MaxSSS threshold for habitat classification
##############################################################################


# ============================================================
# 0) PACKAGES
# ============================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(install.packages.compile.from.source = "never")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

pacman::p_load(
  maxnet,           # MaxEnt modelling
  rlang,            # tidyverse core
  blockCV,          # spatial cross-validation
  dismo,            # SDM utilities (evaluate, randomPoints)
  ecospat,          # Boyce index
  precrec,          # AUC calculation
  raster,           # raster operations (maxnet compatibility)
  terra,            # modern raster/vector handling (pastclim needs this)
  sf,               # simple features
  tidyverse,        # data wrangling + ggplot2
  spThin,           # spatial thinning
  car,              # VIF
  usdm,             # VIF for rasters
  rnaturalearth,    # country polygons
  rnaturalearthdata,# high-res country data
  openxlsx          # write formatted xlsx (Appendix S4)
)

# pastclim must be installed from GitHub if not already available
if (!requireNamespace("pastclim", quietly = TRUE)) {
  pacman::p_load(remotes)
  remotes::install_github("EvolEcolGroup/pastclim")
}
library(pastclim)

# rnaturalearthhires is needed for scale = 10 maps; install silently if missing
if (!requireNamespace("rnaturalearthhires", quietly = TRUE)) {
  install.packages("rnaturalearthhires",
                   repos = "https://ropensci.r-universe.dev",
                   type  = "source")
}

# Parallel processing (adjust workers to your machine)
library(future)
parallel::detectCores()
plan(multisession, workers = 12)


# ============================================================
# 1) STUDY REGION — West Indies (Cuba to Puerto Rico)
# ============================================================
# Include all major islands from the Bahamas/Cuba east to Puerto Rico.
# Exclude Florida and all mainland Central/South America.

# Get island nation polygons at highest resolution
# NOTE: rnaturalearth uses "The Bahamas" not "Bahamas", and names
# can vary between scales. We try multiple approaches.

wi_countries <- c("Cuba", "Jamaica", "Haiti", "Dominican Republic",
                  "Puerto Rico", "The Bahamas", "Cayman Islands",
                  "Turks and Caicos Islands")

wi_sf <- tryCatch(
  ne_countries(country = wi_countries, scale = 10, returnclass = "sf"),
  error = function(e) {
    # Fallback: try alternate names
    cat("Retrying with alternate country names...\n")
    wi_countries_alt <- c("Cuba", "Jamaica", "Haiti", "Dominican Rep.",
                          "Puerto Rico", "Bahamas", "Cayman Is.",
                          "Turks and Caicos Is.")
    ne_countries(country = wi_countries_alt, scale = 10, returnclass = "sf")
  }
)

# Verify we got the Bahamas (critical — many historic records are Bahamian)
wi_names <- wi_sf$name
cat("Countries loaded:", paste(wi_names, collapse = ", "), "\n")

has_bahamas <- any(grepl("bahamas", wi_names, ignore.case = TRUE))
if (!has_bahamas) {
  warning("Bahamas NOT found in country polygons! ",
          "Historic records from the Bahamas will be outside the mask.\n",
          "Try: ne_countries(scale = 10) and check available names.")
} else {
  cat("Bahamas confirmed in study region.\n")
}

# Merge all polygons into one multipolygon for masking
wi_union <- st_union(wi_sf)
wi_union <- st_sf(geometry = wi_union, crs = 4326)

# Bounding box for the study region
wi_bbox   <- st_bbox(wi_union)
wi_extent <- extent(wi_bbox["xmin"], wi_bbox["xmax"],
                    wi_bbox["ymin"], wi_bbox["ymax"])

# Slightly buffer the extent for coastal cells
wi_extent <- extent(
  wi_extent@xmin - 0.5, wi_extent@xmax + 0.5,
  wi_extent@ymin - 0.5, wi_extent@ymax + 0.5
)

# Helper: crop and mask a raster stack to the West Indies islands
crop_to_wi <- function(r) {
  r <- raster::crop(r, wi_extent)
  r <- raster::mask(r, as_Spatial(wi_union))
  r
}

cat("West Indies study region extent:", as.vector(wi_extent), "\n")
cat("Countries included:", paste(wi_countries, collapse = ", "), "\n")


# ============================================================
# 2) HELPER FUNCTIONS
# ============================================================

align_env_to_model <- function(env_stack, model_vars) {
  missing <- setdiff(model_vars, names(env_stack))
  extra   <- setdiff(names(env_stack), model_vars)
  
  if (length(missing) > 0) {
    stop("ENV STACK is missing predictors required by model: ",
         paste(missing, collapse = ", "))
  }
  if (length(extra) > 0) {
    message("Dropping extra predictors not used by model: ",
            paste(extra, collapse = ", "))
  }
  
  env_stack <- raster::subset(env_stack, model_vars)
  stopifnot(all(names(env_stack) == model_vars))
  env_stack
}

# Build future climate env stack from CMIP6 multi-band GeoTIFF
# static_layers: time-invariant layers (VIF-selected topo + geology dummies)
build_future_env <- function(tif_path, bio_indices, bio_names,
                             static_layers, ref_env) {
  fut_all <- raster::stack(tif_path)
  fut <- raster::subset(fut_all, bio_indices)
  names(fut) <- bio_names
  fut <- raster::resample(fut, static_layers)
  env <- raster::stack(fut, static_layers)
  env <- raster::crop(env, raster::extent(ref_env))
  env <- raster::mask(env, as_Spatial(wi_union))
  env
}

# Suitability bin boundaries from a MaxSSS threshold
make_bins <- function(maxss) {
  stopifnot(is.finite(maxss), maxss > 0, maxss < 1)
  list(
    maxss          = maxss,
    thr_unsuitable = maxss,
    thr_poor       = maxss + (1 - maxss) / 3,
    thr_moderate   = maxss + 2 * (1 - maxss) / 3,
    thr_high       = 1
  )
}

# Calculate habitat area (km2) by suitability bin
calc_habitat_areas <- function(rast, bins) {
  cell_areas <- raster::area(rast)
  vals  <- raster::getValues(rast)
  areas <- raster::getValues(cell_areas)
  
  valid <- is.finite(vals) & is.finite(areas)
  vals  <- vals[valid]
  areas <- areas[valid]
  
  data.frame(
    Unsuitable = sum(areas[vals <  bins$thr_unsuitable]),
    Low        = sum(areas[vals >= bins$thr_unsuitable & vals < bins$thr_poor]),
    Moderate   = sum(areas[vals >= bins$thr_poor       & vals < bins$thr_moderate]),
    High       = sum(areas[vals >= bins$thr_moderate]),
    MaxSSS     = bins$maxss
  )
}

# --- Create output directories early so all sections can write ---
dir.create("output",              showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables",       showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures",      showWarnings = FALSE, recursive = TRUE)
dir.create("output/habitat_maps", showWarnings = FALSE, recursive = TRUE)


# ============================================================
# 3) LOAD & CLEAN OCCURRENCE DATA
# ============================================================

# --- 3a) Present records ---
modern.records <- read.csv("data/distribution/Crocodylus_rhombifer_modern.csv")

cat("Modern records loaded:", nrow(modern.records), "\n")
plot(modern.records$longitude, modern.records$latitude,
     main = "Modern C. rhombifer records", pch = 19)

# --- 3b) Historic records ---
historic.records <- read.csv("data/distribution/Crocodylus_rhombifer_historic.csv")

cat("Historic records loaded:", nrow(historic.records), "\n")
plot(historic.records$longitude, historic.records$latitude,
     main = "Historic C. rhombifer records", pch = 19)


# ============================================================
# 4) LOAD CURRENT CLIMATE DATA — WorldClim 2.1 (all 19 bioclim)
# ============================================================
# Uses WorldClim 2.1 at 2.5 arc-min (~5 km) resolution.

wc_dir <- "data/climate/wc2.1_2.5m"

bio_stack <- raster::stack(
  lapply(1:19, function(i) {
    f <- file.path(wc_dir, paste0("wc2.1_2.5m_bio_", i, ".tif"))
    if (!file.exists(f)) stop("Cannot find: ", f)
    raster(f)
  })
)
names(bio_stack) <- paste0("bio_", 1:19)

cat("Loaded WorldClim 2.1 bioclim stack\n")


# ============================================================
# 5) LOAD TOPOGRAPHIC DATA (GMTED-derived, 5 km)
# ============================================================

topo_dir <- "data/topography/GMTED"

load_topo <- function(varname, suffix = "GMTEDmd") {
  f <- file.path(topo_dir, paste0(varname, "_5KMmd_", suffix, ".tif"))
  if (!file.exists(f)) stop("Cannot find topo file: ", f)
  raster(f)
}

elevation_md   <- load_topo("elevation")
slope_md       <- load_topo("slope")
aspect_cos_md  <- load_topo("aspectcosine")
tpi_md         <- load_topo("tpi")
vrm_md         <- load_topo("vrm")

topography <- stack(
  elevation_md, slope_md, aspect_cos_md, tpi_md, vrm_md)

cat("Loaded", nlayers(topography), "topographic layers\n")


# ============================================================
# 5b) LOAD & RASTERIZE GEOLOGY (simplified to 5 classes)
# ============================================================
# 33 raw lithological classes are collapsed into 5 ecologically
# meaningful groups based on substrate type:
#   1. Alluvium      — lowland, wetland-forming substrate
#   2. Limestone     — carbonate/marine strata, karst terrain
#   3. Volcanic      — volcanic rocks, flows, tuff, pyroclastics
#   4. Plutonic      — intrusive/plutonic rocks, gabbro
#   5. Metamorphic   — metamorphic, ultramafic, eclogite
# "Water" and "Undetermined" are excluded (mapped to NA).

geo_sf <- st_read("data/geology/geo6bg.shp", quiet = TRUE)

# Crop to West Indies study region
geo_wi <- st_intersection(geo_sf, wi_union)

# Define the grouping lookup
geo_group_map <- c(
  # Alluvium
  "Quaternary alluvium"                                                          = "alluvium",
  
  # Limestone / carbonate / marine strata
  "Cretaceous marine strata"                                                     = "limestone",
  "Eocene and(or) Paleocene marine strata"                                       = "limestone",
  "Post-Eocene marine strata"                                                    = "limestone",
  "Quaternary and Tertiary marine limestone, sandstone, and shale"               = "limestone",
  "Tertiary marine strata"                                                       = "limestone",
  "Upper Cretaceous marine strata"                                               = "limestone",
  "Jurassic marine and continental strata"                                       = "limestone",
  
  # Volcanic
  "Cretaceous andesitic to silicic volcanic rocks"                               = "volcanic",
  "Cretaceous sedimentary and volcanic rocks"                                    = "volcanic",
  "Cretaceous volcanic rocks"                                                    = "volcanic",
  "Eocene and Paleocene volcanic flows and associated pyroclastic and volcanoge" = "volcanic",
  "Lower Cretaceous flows, breccia, and tuff, partly of marine deposition"       = "volcanic",
  "Mesozoic sedimentary and volcanic rocks"                                      = "volcanic",
  "Quaternary and Tertiary volcanic edifices, flows, tuff, silicic pyroclastric" = "volcanic",
  "Quaternary volcanic edifices, flows, and pyroclastic deposits"                = "volcanic",
  "Tertiary and Cretaceous volcanic rocks"                                       = "volcanic",
  "Tertiary volcanic rocks"                                                      = "volcanic",
  "Volcanic rocks"                                                               = "volcanic",
  "Tertiary and Cretaceous complex of deformed sedimentary rocks"                = "volcanic",
  
  # Plutonic / intrusive
  "Cretaceous plutons, mostly intermediate to silicic"                           = "plutonic",
  "Gabbro and related rocks"                                                     = "plutonic",
  "Intrusive rocks, undivided, mostly intermediate to silicic"                   = "plutonic",
  "Tertiary and Cretaceous plutons, mostly intermediate to silicic"              = "plutonic",
  "Tertiary plutons, mostly intermediate to silicic"                             = "plutonic",
  
  # Metamorphic
  "Eclogite"                                                                     = "metamorphic",
  "Jurassic complex of sedimentary, igneous, and metamorphic rocks"              = "metamorphic",
  "Mesozoic amphibolites and associated metasedimentary rocks"                   = "metamorphic",
  "Mesozoic metasedimentary and metaigneous rocks, low to intermediate metamorp" = "metamorphic",
  "Mesozoic metavolcanic and associated metasedimentary rocks"                   = "metamorphic",
  "Ultramafic rocks"                                                             = "metamorphic"
)

# Apply grouping
geo_wi$geo_group <- geo_group_map[geo_wi$DESCRPTN]

# Drop Water, Undetermined, and any unmatched classes
geo_wi <- geo_wi[!is.na(geo_wi$geo_group), ]

geo_group_classes <- sort(unique(geo_wi$geo_group))
cat("Simplified geology classes:", paste(geo_group_classes, collapse = ", "), "\n")

# Assign integer codes
geo_wi$geo_code <- match(geo_wi$geo_group, geo_group_classes)

# Reference raster at topo resolution, cropped to WI
ref_raster <- crop_to_wi(topography[[1]])

# Rasterize geology
geo_raster <- rasterize(as_Spatial(geo_wi), ref_raster,
                        field = "geo_code", fun = "first")

# Create dummy binary rasters for each simplified class
geo_dummy_names <- paste0("geo_", geo_group_classes)

geo_dummies <- raster::stack(lapply(seq_along(geo_group_classes), function(i) {
  dummy <- geo_raster == i
  dummy <- raster::mask(dummy, ref_raster)
  names(dummy) <- geo_dummy_names[i]
  dummy
}))

# Drop any class with <10 cells (unlikely after grouping, but safe)
geo_cell_counts <- sapply(1:nlayers(geo_dummies), function(i) {
  sum(values(geo_dummies[[i]]) == 1, na.rm = TRUE)
})
keep_geo <- geo_cell_counts >= 10
if (any(!keep_geo)) {
  cat("Dropping rare geology classes (<10 cells):",
      names(geo_dummies)[!keep_geo], "\n")
  geo_dummies <- raster::subset(geo_dummies, which(keep_geo))
}

cat("Retained", nlayers(geo_dummies), "geology dummy layers:",
    names(geo_dummies), "\n")


# ============================================================
# 6) BUILD PRESENT ENVIRONMENTAL STACK + VIF SELECTION
# ============================================================

# Resample climate to match topography resolution
climate_resampled <- resample(bio_stack, topography)

# Stack all 19 bioclim + 5 topo (continuous variables only for VIF)
env_continuous <- stack(climate_resampled, topography)

# Crop and mask to West Indies
env_continuous_wi <- crop_to_wi(env_continuous)

cat("Continuous predictors before VIF:", nlayers(env_continuous_wi), "\n")

# --- VIF-based variable selection (continuous variables only) ---
# Geology dummies are categorical/binary and excluded from VIF.
# Extract raster values as a data frame (sample if too large for memory)
env_vals <- as.data.frame(values(env_continuous_wi))
env_vals <- na.omit(env_vals)
cat("Non-NA cells for VIF:", nrow(env_vals), "\n")

# If very large, take a random sample for speed
if (nrow(env_vals) > 50000) {
  set.seed(42)
  env_vals <- env_vals[sample(nrow(env_vals), 50000), ]
  cat("Sampled 50000 cells for VIF computation\n")
}

vif_results <- usdm::vifstep(env_vals, th = 3)
cat("\n--- VIF results (threshold = 3) ---\n")
print(vif_results)

# Retain only selected continuous variables
selected_continuous <- vif_results@results$Variables
cat("\nRetained continuous predictors after VIF:", selected_continuous, "\n")

# Build the final present environmental stack:
#   VIF-selected continuous variables + geology dummies
env_vif <- raster::subset(env_continuous_wi, selected_continuous)

# Resample geology dummies to match
geo_dummies_matched <- resample(geo_dummies, env_vif)
geo_dummies_matched <- raster::mask(geo_dummies_matched, env_vif[[1]])

# Round dummy values back to 0/1 after resampling interpolation
# IMPORTANT: round() strips layer names, so we save and restore them
geo_names_saved <- names(geo_dummies_matched)
for (i in 1:nlayers(geo_dummies_matched)) {
  geo_dummies_matched[[i]] <- round(geo_dummies_matched[[i]])
}
names(geo_dummies_matched) <- geo_names_saved

cat("Geology dummies after processing:", names(geo_dummies_matched), "\n")

env_present <- stack(env_vif, geo_dummies_matched)

cat("Final present env stack:", names(env_present), "\n")

# Define predictor column names for use throughout the script
predictor_cols <- names(env_present)

# Build the time-invariant static layers (topo + geology) for reuse
# in past and future env stacks. These are already correctly matched.
bio_predictors  <- predictor_cols[grepl("^bio_", predictor_cols)]
geo_predictors  <- predictor_cols[grepl("^geo_", predictor_cols)]
topo_predictors <- predictor_cols[!grepl("^bio_|^geo_", predictor_cols)]

static_wi <- raster::subset(env_present, c(topo_predictors, geo_predictors))

cat("Bio predictors:", bio_predictors, "\n")
cat("Topo predictors:", topo_predictors, "\n")
cat("Geology predictors:", geo_predictors, "\n")
cat("Static layers (topo + geo):", names(static_wi), "\n")


# ============================================================
# 7) LOAD PALEOCLIMATE — CHELSA-TraCE21k VIA pastclim
# ============================================================
# Using CHELSA-TraCE21k for the 7–4 kya period.
# We compute a cell-wise median over 500-year intervals from 7000–4000 BP.
# Historic record climate values are extracted from this median surface.

# Set pastclim data path (adjust to your local path)
# pastclim::set_data_path("data/climate/past/pastclim_data")

# Map VIF-selected bioclim predictor names to pastclim variable names
pastclim_vars <- gsub("^bio_", "bio", bio_predictors)
pastclim_vars <- gsub("bio(\\d)$", "bio0\\1", pastclim_vars)

cat("Requesting pastclim variables:", pastclim_vars, "\n")

# Download any missing variables
cat("Downloading CHELSA TraCE21k variables (if not already cached)...\n")
pastclim::download_dataset(
  dataset       = "CHELSA_trace21k_1.0_0.5m_vsi",
  bio_variables = pastclim_vars
)


# --- Build the 7–4 kya MEDIAN raster stack ---
# Fetch region_slice at 500-year intervals from 7000 to 4000 BP,
# then compute the cell-wise median across all 7 time slices.

past_time_steps <- seq(-7000, -4000, by = 500)
cat("Fetching TraCE21k region slices for median:", past_time_steps, "\n")
cat("This may take several minutes (streaming", length(past_time_steps),
    "slices via VSI)...\n")

past_slice_list <- list()
for (t in past_time_steps) {
  cat("  Fetching time_bp =", t, "...\n")
  r_terra <- pastclim::region_slice(
    time_bp       = t,
    bio_variables = pastclim_vars,
    dataset       = "CHELSA_trace21k_1.0_0.5m_vsi",
    ext           = terra::ext(wi_extent[1], wi_extent[2],
                               wi_extent[3], wi_extent[4])
  )
  past_slice_list[[as.character(t)]] <- r_terra
}

# Compute cell-wise median for each bioclim variable across all time slices
cat("Computing cell-wise median across", length(past_time_steps), "time slices...\n")

# Stack all slices per variable, then take median
median_layers <- list()
for (v in seq_along(pastclim_vars)) {
  vname <- pastclim_vars[v]
  # Extract the v-th layer from each time slice as a terra SpatRaster
  var_stack <- terra::rast(lapply(past_slice_list, function(r) r[[v]]))
  # Cell-wise median
  var_median <- terra::app(var_stack, fun = median, na.rm = TRUE)
  names(var_median) <- vname
  median_layers[[v]] <- var_median
}

# Combine into a single SpatRaster
past_median_terra <- terra::rast(median_layers)

# Convert terra SpatRaster to raster stack for maxnet compatibility
past_rasters <- raster::stack(raster::raster(past_median_terra))
for (i in 2:terra::nlyr(past_median_terra)) {
  past_rasters <- raster::addLayer(past_rasters, raster::raster(past_median_terra[[i]]))
}

# Fix names: pastclim uses bio01, bio02 ... → convert to bio_1, bio_2 ...
names(past_rasters) <- gsub("^bio0?(\\d+)$", "bio_\\1", names(past_rasters))

cat("Past median raster stack (7–4 kya):", names(past_rasters), "\n")

# Resample past climate to match static layers, stack together
past_climate_rasters <- resample(past_rasters, static_wi)
past_env <- stack(past_climate_rasters, static_wi)
past_env <- raster::mask(past_env, as_Spatial(wi_union))

# Ensure only VIF-selected + geology predictors are kept
past_env <- raster::subset(past_env, predictor_cols)

cat("Past (7–4 kya median) env stack:", names(past_env), "\n")
stopifnot(all(names(past_env) == predictor_cols))


# ============================================================
# 8) SPATIAL THINNING
# ============================================================

# --- 8a) Modern records ---
modern.records$species <- "rhombifer"

set.seed(42)
rhombifer_thin_modern <- thin(
  loc.data  = as.data.frame(modern.records),
  lat.col   = "latitude",
  long.col  = "longitude",
  spec.col  = "species",
  thin.par  = 5,    # 5 km minimum distance
  reps      = 100,
  write.files = FALSE,
  locs.thinned.list.return = TRUE
)

best_thin <- rhombifer_thin_modern[[which.max(sapply(rhombifer_thin_modern, nrow))]]
modern_pres <- modern.records %>%
  inner_join(
    data.frame(longitude = best_thin$Longitude, latitude = best_thin$Latitude),
    by = c("longitude", "latitude")
  )

modern_thin <- modern_pres %>%
  transmute(
    species = 1,
    climate.period = climate.period,
    Longitude = longitude,
    Latitude  = latitude
  ) %>%
  as.data.frame()

cat("Modern records after distance thinning:", nrow(modern_thin), "\n")

# Reduce to one record per unique raster grid cell to prevent
# pseudo-replication from spatially clustered records (e.g., Zapata Swamp).
modern_cells <- raster::cellFromXY(env_present,
                                   modern_thin[, c("Longitude", "Latitude")])
modern_thin <- modern_thin[!duplicated(modern_cells), ]

cat("Modern records after grid-cell deduplication:", nrow(modern_thin), "\n")

# --- 8b) Historic records ---
historic.records$species <- "rhombifer"

set.seed(42)
rhombifer_thin_historic <- thin(
  loc.data  = as.data.frame(historic.records),
  lat.col   = "latitude",
  long.col  = "longitude",
  spec.col  = "species",
  thin.par  = 5,    # 5 km minimum distance
  reps      = 100,
  write.files = FALSE,
  locs.thinned.list.return = TRUE
)

best_thin_hist <- rhombifer_thin_historic[[which.max(sapply(rhombifer_thin_historic, nrow))]]
historic_pres <- historic.records %>%
  inner_join(
    data.frame(longitude = best_thin_hist$Longitude,
               latitude  = best_thin_hist$Latitude),
    by = c("longitude", "latitude")
  )

historic_thin <- historic_pres %>%
  transmute(
    species = 1,
    climate.period = climate.period,
    Longitude = longitude,
    Latitude  = latitude
  ) %>%
  as.data.frame()

cat("Historic records after thinning:", nrow(historic_thin), "\n")


# ============================================================
# 9) PRESENT-ONLY HSM
# ============================================================

climate_mask <- env_present
climate_crop <- climate_mask

# Convert occurrences to spatial
dat_current <- modern_thin %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  as_Spatial()

plot(climate_mask[[1]], main = "West Indies — Present Env + Modern Occurrences")
points(dat_current@coords[, 1], dat_current@coords[, 2], col = "red", pch = 19, cex = 1.5)

# Background points
set.seed(42)
bg_pts <- randomPoints(climate_mask, n = 10000) %>%
  as.data.frame() %>%
  dplyr::select(Longitude = x, Latitude = y) %>%
  mutate(species = 0) %>%
  dplyr::select(species, Longitude, Latitude)

dat_current_df <- as.data.frame(dat_current) %>%
  dplyr::select(species, Longitude = coords.x1, Latitude = coords.x2)

plot(climate_mask[[1]], main = "Present Occurrences (red) & Background (blue)")
points(dat_current_df$Longitude, dat_current_df$Latitude, col = "red", pch = 19, cex = 1.5)
points(bg_pts$Longitude, bg_pts$Latitude, pch = ".", col = "blue", cex = 1.5)

# Combine presence + background
pb_data <- rbind(dat_current_df, bg_pts)
pb_data <- st_as_sf(pb_data, coords = c("Longitude", "Latitude"), crs = 4326)
pb_data <- st_transform(pb_data, crs = st_crs(climate_mask))

# Spatial cross-validation
n_pres_modern <- sum(dat_current_df$species == 1)
cv_results <- NULL

if (n_pres_modern >= 10) {
  tryCatch({
    set.seed(42)
    cv_results <- blockCV::cv_spatial(
      x = pb_data,
      raster = climate_crop,
      k = 5,
      selection = "systematic"
    )
    blockCV::cv_plot(cv_results)
  }, error = function(e) {
    cat("WARNING: Spatial CV failed:", conditionMessage(e), "\n")
    cv_results <<- NULL
  })
} else {
  cat("NOTE: Only", n_pres_modern,
      "presence records — skipping spatial cross-validation.\n")
}

# Extract environmental data
pb_coords  <- st_coordinates(pb_data)
pb_data_df <- as.data.frame(pb_data)
pb_data_df$Longitude <- pb_coords[, 1]
pb_data_df$Latitude  <- pb_coords[, 2]
pb_data_df$ID <- 1:nrow(pb_data_df)

pb_climate_df <- raster::extract(climate_mask, pb_data, df = TRUE)
pb_all <- merge(pb_data_df, pb_climate_df, by = "ID")
pb_all_clean <- na.omit(pb_all)

stopifnot(all(predictor_cols %in% names(pb_all_clean)))

pb <- pb_all_clean$species
pb_climate <- as.data.frame(pb_all_clean[, predictor_cols])

pb_data_clean <- st_as_sf(
  pb_all_clean,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)

stopifnot(length(pb) == nrow(pb_climate))
stopifnot(all(colnames(pb_climate) == predictor_cols))

# Spatial CV model evaluation (if CV succeeded and enough presences)
if (!is.null(cv_results) && n_pres_modern >= 10) {
  folds <- cv_results$folds
  if (is.null(folds) || length(folds) == 0) folds <- cv_results$folds_list
  
  if (!is.null(folds) && length(folds) > 0) {
    AUCs <- CBIs <- ORMTPs <- OR10s <- rep(NA_real_, length(folds))
    
    for (k in seq_along(folds)) {
      trainSet <- unlist(folds[[k]][1])
      testSet  <- unlist(folds[[k]][2])
      
      trainSet <- trainSet[!is.na(trainSet)]
      testSet  <- testSet[!is.na(testSet)]
      
      # Filter to valid index range (pb may be shorter than original pb_data)
      trainSet <- trainSet[trainSet <= length(pb)]
      testSet  <- testSet[testSet <= length(pb)]
      
      if (length(trainSet) == 0 || length(testSet) == 0) next
      if (sum(pb[trainSet] == 1, na.rm = TRUE) == 0 ||
          sum(pb[trainSet] == 0, na.rm = TRUE) == 0) next
      
      mx <- maxnet(
        pb[trainSet],
        pb_climate[trainSet, ],
        maxnet.formula(pb[trainSet], pb_climate[trainSet, ], classes = "lq")
      )
      
      prediction <- predict(mx, pb_climate, type = "cloglog")
      
      test_lab <- pb[testSet]
      if (length(unique(test_lab)) > 1) {
        test_pred <- prediction[testSet]
        AUCs[k] <- precrec::auc(precrec::evalmod(test_pred, test_lab))[1, 4]
      }
      
      p_vals    <- prediction[pb == 1]
      ORMTPs[k] <- min(p_vals, na.rm = TRUE)
      OR10s[k]  <- quantile(p_vals, 0.10, na.rm = TRUE)
      
      try({
        CBIs[k] <- ecospat.boyce(prediction, p_vals, PEplot = FALSE)$cor
      }, silent = TRUE)
    }
    
    evaluation_stats_modern <- tibble(
      AUC_mean   = mean(AUCs,   na.rm = TRUE),
      AUC_sd     = sd(AUCs,     na.rm = TRUE),
      CBI_mean   = mean(CBIs,   na.rm = TRUE),
      CBI_sd     = sd(CBIs,     na.rm = TRUE),
      ORMTP_mean = mean(ORMTPs, na.rm = TRUE),
      OR10_mean  = mean(OR10s,  na.rm = TRUE)
    )
    print(evaluation_stats_modern)
  }
}

# Fit final modern MaxEnt model
me <- maxnet(
  p    = pb,
  data = pb_climate,
  f    = maxnet.formula(pb, pb_climate, classes = "lq")
)

pred_all <- predict(me, pb_climate, type = "cloglog")

eval_all <- dismo::evaluate(
  p = pred_all[pb == 1],
  a = pred_all[pb == 0]
)

maxss_modern <- eval_all@t[which.max(eval_all@TPR + eval_all@TNR)]
cat("Modern model MaxSSS threshold =", round(maxss_modern, 4), "\n")

train_vars_modern <- colnames(pb_climate)

# Predict to present West Indies
final_env <- env_present
final_env <- align_env_to_model(final_env, train_vars_modern)
final_suitability <- raster::predict(final_env, me, type = "cloglog")

plot(final_suitability, main = "C. rhombifer — Present Suitability (Modern Model)")

png("output/final_suitability_wi_modern.png", width = 2000, height = 1600, res = 200)
plot(final_suitability, main = "C. rhombifer — Present Suitability (Modern Model)")
dev.off()

# Save modern PB objects
pb_all_clean_modern  <- pb_all_clean
pb_data_clean_modern <- pb_data_clean


# ============================================================
# 10) COMBINED MODERN + PAST HSM
# ============================================================

ENV_present    <- env_present
past_env_wi    <- past_env

stopifnot(all(names(ENV_present) == names(past_env_wi)))

# Helper to build PB table
build_pb_table <- function(occ_df, env_stack, n_bg = 10000, period = NA_character_) {
  stopifnot(all(c("Longitude", "Latitude") %in% names(occ_df)))
  
  pres_sf <- occ_df %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(env_stack))
  pres_sf$species <- 1
  
  set.seed(42)
  bg_xy <- dismo::randomPoints(env_stack, n = n_bg) %>%
    as.data.frame() %>%
    dplyr::rename(Longitude = x, Latitude = y)
  bg_sf <- st_as_sf(bg_xy, coords = c("Longitude", "Latitude"),
                    crs = st_crs(env_stack))
  bg_sf$species <- 0
  
  all_sf <- rbind(pres_sf[, "species"], bg_sf[, "species"])
  all_sp <- as(all_sf, "Spatial")
  vals   <- raster::extract(env_stack, all_sp, df = TRUE)
  
  X <- vals[, -1, drop = FALSE]
  y <- all_sf$species
  
  keep <- complete.cases(X) & is.finite(y)
  X <- X[keep, , drop = FALSE]
  y <- y[keep]
  
  list(
    X = as.data.frame(X),
    y = as.numeric(y),
    period = period,
    points_sf = all_sf[keep, ]
  )
}

# Prepare occurrence data
modern_pres_df <- modern_thin %>% dplyr::select(Longitude, Latitude)
historic_pres_df <- historic_thin %>% dplyr::select(Longitude, Latitude)

# Build PB tables for each period
mod_obj  <- build_pb_table(modern_pres_df,   ENV_present, n_bg = 10000, period = "Present")
past_obj <- build_pb_table(historic_pres_df, past_env_wi, n_bg = 10000, period = "Past")

dat_modern   <- bind_cols(species = mod_obj$y,  mod_obj$X)  %>%
  mutate(period = "Present") %>% as.data.frame()
dat_historic <- bind_cols(species = past_obj$y, past_obj$X) %>%
  mutate(period = "Past") %>% as.data.frame()

# Combine
dat_all <- bind_rows(dat_modern, dat_historic)
y <- dat_all$species
predictor_cols_combined <- setdiff(names(dat_all), c("species", "period"))
X <- dat_all[, predictor_cols_combined, drop = FALSE]

stopifnot(nrow(X) == length(y))

# Fit combined model
me_combined <- maxnet::maxnet(
  p = y,
  data = X,
  f = maxnet::maxnet.formula(p = y, data = X, classes = "lq")
)

train_vars_combined <- colnames(X)

# Predict combined model to present
ENV_present_pred <- raster::subset(ENV_present, train_vars_combined)
pred_present_from_combined <- raster::predict(
  ENV_present_pred, me_combined, type = "cloglog"
)

plot(pred_present_from_combined,
     main = "C. rhombifer — Present Suitability (Combined Model)", zlim = c(0, 1))

png("output/combined_model_present_wi.png", width = 2000, height = 1600, res = 200)
plot(pred_present_from_combined,
     main = "C. rhombifer — Present Suitability (Combined Model)", zlim = c(0, 1))
dev.off()


# ============================================================
# 11) PAST-ONLY HSM (7–4 kya median model)
# ============================================================

historic_sf <- historic_thin %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

plot(past_env_wi[[1]], main = "Past (7–4 kya median) Environment + Historic Occurrences")
points(st_coordinates(historic_sf)[, 1], st_coordinates(historic_sf)[, 2],
       col = "red", pch = 19, cex = 1.5)

set.seed(42)
bg_pts_past <- randomPoints(past_env_wi, n = 10000) %>%
  as.data.frame() %>%
  dplyr::select(Longitude = x, Latitude = y) %>%
  mutate(species = 0)

historic_coords <- st_coordinates(historic_sf)
historic_df <- data.frame(
  Longitude = historic_coords[, "X"],
  Latitude  = historic_coords[, "Y"],
  species   = 1
)

pb_data_past <- rbind(historic_df, bg_pts_past)
pb_data_past <- st_as_sf(pb_data_past, coords = c("Longitude", "Latitude"), crs = 4326)

pb_coords_past  <- st_coordinates(pb_data_past)
pb_data_past_df <- as.data.frame(pb_data_past)
pb_data_past_df$Longitude <- pb_coords_past[, 1]
pb_data_past_df$Latitude  <- pb_coords_past[, 2]
pb_data_past_df$ID <- 1:nrow(pb_data_past_df)

pb_env_df <- raster::extract(past_env_wi, pb_data_past, df = TRUE)
pb_all_past <- merge(pb_data_past_df, pb_env_df, by = "ID")
pb_all_clean_past <- na.omit(pb_all_past)

pb_past     <- pb_all_clean_past$species
pb_env_past <- as.data.frame(pb_all_clean_past[, predictor_cols])

pb_data_clean_past <- st_as_sf(
  pb_all_clean_past,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)

stopifnot(length(pb_past) == nrow(pb_env_past))
stopifnot(all(colnames(pb_env_past) == predictor_cols))

# Fit past MaxEnt model
me_past <- maxnet(
  p    = pb_past,
  data = pb_env_past,
  f    = maxnet.formula(p = pb_past, data = pb_env_past, classes = "lq")
)

# Predict past suitability
past_env_aligned <- align_env_to_model(past_env_wi, colnames(pb_env_past))
suitability_past <- raster::predict(past_env_aligned, me_past, type = "cloglog")

plot(suitability_past, main = "C. rhombifer — 7–4 kya Median Suitability (Past-Only Model)")

png("output/past_7_4kya_suitability_wi.png", width = 2000, height = 1600, res = 200)
plot(suitability_past, main = "C. rhombifer — 7–4 kya Median Suitability (Past-Only Model)")
dev.off()


# ============================================================
# 11b) PAST-ONLY MODEL → PRESENT-DAY PROJECTION
# ============================================================
# Project the past-only model onto present-day environmental conditions.
# The prediction itself does not depend on thresholds; the Zapata
# suitability analysis (which uses thr_past) follows in section 13b
# after thresholds are defined.

past_env_present <- align_env_to_model(env_present, colnames(pb_env_past))
pred_present_from_past <- raster::predict(past_env_present, me_past, type = "cloglog")

plot(pred_present_from_past,
     main = "C. rhombifer — Present Suitability (Past-Only Model)", zlim = c(0, 1))
png("output/past_model_present_projection_wi.png", width = 2000, height = 1600, res = 200)
plot(pred_present_from_past,
     main = "C. rhombifer — Present Suitability (Past-Only Model)", zlim = c(0, 1))
dev.off()


# ============================================================
# 12) FUTURE PROJECTIONS — CMIP6 (MIROC6 SSP245 & SSP585)
# ============================================================

bio_band_indices <- as.integer(gsub("bio_", "", bio_predictors))
bio_band_names   <- bio_predictors

# static_wi already defined in section 6 (VIF-selected topo + geology)

# === Future projections using MODERN model (me) ===

# --- SSP245 ---
future_env_ssp245 <- build_future_env(
  tif_path      = "data/climate/cmip6/wc2.1_2.5m_bioc_MIROC6_ssp245_2061-2080.tif",
  bio_indices   = bio_band_indices,
  bio_names     = bio_band_names,
  static_layers = static_wi,
  ref_env       = final_env
)
future_env_ssp245 <- align_env_to_model(future_env_ssp245, train_vars_modern)
future_pred_ssp245 <- raster::predict(future_env_ssp245, me, type = "cloglog")

plot(future_pred_ssp245,
     main = "C. rhombifer — SSP245 2061–2080 (Modern Model)", zlim = c(0, 1))
png("output/future_ssp245_modern_wi.png", width = 2000, height = 1600, res = 200)
plot(future_pred_ssp245,
     main = "C. rhombifer — SSP245 2061–2080 (Modern Model)", zlim = c(0, 1))
dev.off()

# --- SSP585 ---
future_env_ssp585 <- build_future_env(
  tif_path      = "data/climate/cmip6/wc2.1_2.5m_bioc_MIROC6_ssp585_2061-2080.tif",
  bio_indices   = bio_band_indices,
  bio_names     = bio_band_names,
  static_layers = static_wi,
  ref_env       = final_env
)
future_env_ssp585 <- align_env_to_model(future_env_ssp585, train_vars_modern)
future_pred_ssp585 <- raster::predict(future_env_ssp585, me, type = "cloglog")

plot(future_pred_ssp585,
     main = "C. rhombifer — SSP585 2061–2080 (Modern Model)", zlim = c(0, 1))
png("output/future_ssp585_modern_wi.png", width = 2000, height = 1600, res = 200)
plot(future_pred_ssp585,
     main = "C. rhombifer — SSP585 2061–2080 (Modern Model)", zlim = c(0, 1))
dev.off()

# === Future projections using COMBINED model (me_combined) ===

# --- SSP245 ---
future_env_ssp245_comb <- build_future_env(
  tif_path      = "data/climate/cmip6/wc2.1_2.5m_bioc_MIROC6_ssp245_2061-2080.tif",
  bio_indices   = bio_band_indices,
  bio_names     = bio_band_names,
  static_layers = static_wi,
  ref_env       = ENV_present
)
future_env_ssp245_comb <- align_env_to_model(future_env_ssp245_comb, train_vars_combined)
future_pred_ssp245_comb <- raster::predict(future_env_ssp245_comb, me_combined, type = "cloglog")

plot(future_pred_ssp245_comb,
     main = "C. rhombifer — SSP245 2061–2080 (Combined Model)", zlim = c(0, 1))
png("output/future_ssp245_combined_wi.png", width = 2000, height = 1600, res = 200)
plot(future_pred_ssp245_comb,
     main = "C. rhombifer — SSP245 2061–2080 (Combined Model)", zlim = c(0, 1))
dev.off()

# --- SSP585 ---
future_env_ssp585_comb <- build_future_env(
  tif_path      = "data/climate/cmip6/wc2.1_2.5m_bioc_MIROC6_ssp585_2061-2080.tif",
  bio_indices   = bio_band_indices,
  bio_names     = bio_band_names,
  static_layers = static_wi,
  ref_env       = ENV_present
)
future_env_ssp585_comb <- align_env_to_model(future_env_ssp585_comb, train_vars_combined)
future_pred_ssp585_comb <- raster::predict(future_env_ssp585_comb, me_combined, type = "cloglog")

plot(future_pred_ssp585_comb,
     main = "C. rhombifer — SSP585 2061–2080 (Combined Model)", zlim = c(0, 1))
png("output/future_ssp585_combined_wi.png", width = 2000, height = 1600, res = 200)
plot(future_pred_ssp585_comb,
     main = "C. rhombifer — SSP585 2061–2080 (Combined Model)", zlim = c(0, 1))
dev.off()


# ============================================================
# 13) MODEL COMPARISONS
# ============================================================

eval_maxnet_on_pb <- function(model, X, y, model_name = "model") {
  X <- as.data.frame(X)
  y <- as.numeric(y)
  
  keep <- complete.cases(X) & is.finite(y)
  X <- X[keep, , drop = FALSE]
  y <- y[keep]
  
  if (length(unique(y)) < 2) {
    warning(model_name, ": only one class present; cannot evaluate.")
    return(tibble(
      Model = model_name, N = nrow(X), N_pres = sum(y == 1),
      N_bg = sum(y == 0), AUC = NA_real_, CBI = NA_real_,
      ORMTP = NA_real_, OR10 = NA_real_, MSS = NA_real_
    ))
  }
  
  pred <- predict(model, X, type = "cloglog")
  pres_pred <- pred[y == 1]
  back_pred <- pred[y == 0]
  
  if (length(pres_pred) < 2 || length(back_pred) < 2) {
    warning(model_name, ": insufficient presences/background.")
    return(tibble(
      Model = model_name, N = nrow(X), N_pres = sum(y == 1),
      N_bg = sum(y == 0), AUC = NA_real_, CBI = NA_real_,
      ORMTP = min(pres_pred, na.rm = TRUE),
      OR10  = as.numeric(quantile(pres_pred, 0.10, na.rm = TRUE, names = FALSE)),
      MSS = NA_real_
    ))
  }
  
  ev  <- dismo::evaluate(p = pres_pred, a = back_pred)
  mss <- ev@t[which.max(ev@TPR + ev@TNR)]
  
  ormtp <- min(pres_pred, na.rm = TRUE)
  or10  <- as.numeric(quantile(pres_pred, 0.10, na.rm = TRUE, names = FALSE))
  
  cbi <- NA_real_
  try({
    boy <- ecospat::ecospat.boyce(
      fit = pred, obs = pres_pred,
      nclass = 0, window.w = "default", res = 100, PEplot = FALSE
    )
    cbi <- as.numeric(boy$cor)
  }, silent = TRUE)
  
  tibble(
    Model = model_name, N = nrow(X), N_pres = sum(y == 1),
    N_bg = sum(y == 0), AUC = as.numeric(ev@auc), CBI = cbi,
    ORMTP = as.numeric(ormtp), OR10 = as.numeric(or10),
    MSS = as.numeric(mss)
  )
}

# Build evaluation PB data
X_modern_eval  <- pb_climate
y_modern_eval  <- pb

X_past_eval    <- pb_env_past
y_past_eval    <- pb_past

X_comb_mod  <- dat_modern  %>% dplyr::select(-species, -period)
y_comb_mod  <- dat_modern$species
X_comb_past <- dat_historic %>% dplyr::select(-species, -period)
y_comb_past <- dat_historic$species

X_comb_all <- bind_rows(X_comb_mod, X_comb_past)
y_comb_all <- c(y_comb_mod, y_comb_past)

eval_results <- bind_rows(
  eval_maxnet_on_pb(me,         X_modern_eval, y_modern_eval, "Modern (present-only)"),
  eval_maxnet_on_pb(me_past,    X_past_eval,   y_past_eval,   "Past-only (7–4 kya median)"),
  eval_maxnet_on_pb(me_combined, X_comb_mod,   y_comb_mod,    "Combined (present PB)"),
  eval_maxnet_on_pb(me_combined, X_comb_past,  y_comb_past,   "Combined (past PB)"),
  eval_maxnet_on_pb(me_combined, X_comb_all,   y_comb_all,    "Combined (pooled PB)")
)

print(eval_results)
write.csv(eval_results, "output/evaluation_metrics_by_model.csv", row.names = FALSE)


# ============================================================
# 14) MaxSSS-BASED HABITAT BINS + AREA SUMMARIES
# ============================================================

get_mss <- function(label) {
  thr <- eval_results %>%
    dplyr::filter(Model == label) %>%
    dplyr::pull(MSS)
  if (length(thr) != 1 || !is.finite(thr)) {
    stop("Could not get MSS for: ", label)
  }
  thr
}

thr_modern   <- get_mss("Modern (present-only)")
thr_past     <- get_mss("Past-only (7–4 kya median)")
# Use the past-only model's MaxSSS threshold for the combined model.
# The combined model's own thresholds are problematic: the pooled PB
# threshold (0.758) is inflated by geographic separability of Zapata
# records, and the past PB threshold (0.146) is too permissive.
# The past-only model's threshold reflects discrimination learned
# from the full historic range and provides a more ecologically
# defensible suitability cutoff.
thr_combined <- thr_past

cat("\n--- MaxSSS thresholds ---\n")
cat("Modern:",   round(thr_modern,   4), "\n")
cat("Past:",     round(thr_past,     4), "\n")
cat("Combined (using past-only model threshold):", round(thr_combined, 4), "\n\n")


# ============================================================
# 13b) ZAPATA SUITABILITY ANALYSIS (Past-only model → present)
# ============================================================
# Evaluate the past-only model's present-day prediction at modern
# Zapata Peninsula records. Requires thr_past (defined above).

zapata_suit_past_model <- raster::extract(
  pred_present_from_past,
  modern_thin[, c("Longitude", "Latitude")]
)

cat("\n--- Past-only model: present-day suitability at modern (Zapata) records ---\n")
cat("N modern records:", length(zapata_suit_past_model), "\n")
cat("Mean suitability:", round(mean(zapata_suit_past_model, na.rm = TRUE), 4), "\n")
cat("Median suitability:", round(median(zapata_suit_past_model, na.rm = TRUE), 4), "\n")
cat("Range:", round(min(zapata_suit_past_model, na.rm = TRUE), 4), "–",
    round(max(zapata_suit_past_model, na.rm = TRUE), 4), "\n")
cat("Past-only MaxSSS threshold:", round(thr_past, 4), "\n")
cat("Records above threshold:", sum(zapata_suit_past_model >= thr_past, na.rm = TRUE),
    "of", sum(!is.na(zapata_suit_past_model)), "\n")

# Percentile rank of Zapata within the distribution of suitable cells
past_pres_vals <- raster::getValues(pred_present_from_past)
past_pres_vals <- past_pres_vals[is.finite(past_pres_vals) & past_pres_vals >= thr_past]
zapata_median <- median(zapata_suit_past_model, na.rm = TRUE)
zapata_pctile <- round(100 * mean(past_pres_vals <= zapata_median), 1)
cat("Percentile rank of Zapata median among suitable cells:", zapata_pctile, "%\n")

# --- Area of suitable habitat (past-only model → present day) ---
# Using the same binning approach as section 14 for consistency
bins_past_pres <- make_bins(thr_past)

# Calculate areas for both time periods using the helper function
areas_past_present <- calc_habitat_areas(pred_present_from_past, bins_past_pres)
areas_past_hol     <- calc_habitat_areas(suitability_past, bins_past_pres)

# Total study area (all valid cells)
past_pres_all_vals  <- raster::getValues(pred_present_from_past)
past_pres_cell_area <- raster::getValues(raster::area(pred_present_from_past))
past_pres_valid     <- is.finite(past_pres_all_vals) & is.finite(past_pres_cell_area)
total_study_area    <- sum(past_pres_cell_area[past_pres_valid])

cat("\n--- Past-only model: MaxSSS bin thresholds ---\n")
cat("MaxSSS (suitable threshold):", round(bins_past_pres$maxss, 4), "\n")
cat("Low suitability range:      ", round(bins_past_pres$maxss, 4), "–",
    round(bins_past_pres$thr_poor, 4), "\n")
cat("Moderate suitability range: ", round(bins_past_pres$thr_poor, 4), "–",
    round(bins_past_pres$thr_moderate, 4), "\n")
cat("High suitability range:     ", round(bins_past_pres$thr_moderate, 4), "– 1.0\n")

# --- Present-day areas ---
suit_present <- areas_past_present$Low + areas_past_present$Moderate +
  areas_past_present$High

cat("\n--- Past-only model → present-day habitat area by bin ---\n")
cat(sprintf("  Unsuitable:  %8.0f km²  (%5.1f%%)\n",
            areas_past_present$Unsuitable,
            100 * areas_past_present$Unsuitable / total_study_area))
cat(sprintf("  Low:         %8.0f km²  (%5.1f%%)\n",
            areas_past_present$Low,
            100 * areas_past_present$Low / total_study_area))
cat(sprintf("  Moderate:    %8.0f km²  (%5.1f%%)\n",
            areas_past_present$Moderate,
            100 * areas_past_present$Moderate / total_study_area))
cat(sprintf("  High:        %8.0f km²  (%5.1f%%)\n",
            areas_past_present$High,
            100 * areas_past_present$High / total_study_area))
cat(sprintf("  ─────────────────────────────\n"))
cat(sprintf("  Total suitable: %5.0f km²  (%5.1f%%)\n",
            suit_present, 100 * suit_present / total_study_area))
cat(sprintf("  Total study area: %.0f km²\n", total_study_area))

# --- Holocene (7–4 kya) areas ---
suit_hol <- areas_past_hol$Low + areas_past_hol$Moderate + areas_past_hol$High

past_hol_all_vals  <- raster::getValues(suitability_past)
past_hol_cell_area <- raster::getValues(raster::area(suitability_past))
past_hol_valid     <- is.finite(past_hol_all_vals) & is.finite(past_hol_cell_area)
total_study_area_hol <- sum(past_hol_cell_area[past_hol_valid])

cat("\n--- Past-only model → Holocene (7–4 kya) habitat area by bin ---\n")
cat(sprintf("  Unsuitable:  %8.0f km²  (%5.1f%%)\n",
            areas_past_hol$Unsuitable,
            100 * areas_past_hol$Unsuitable / total_study_area_hol))
cat(sprintf("  Low:         %8.0f km²  (%5.1f%%)\n",
            areas_past_hol$Low,
            100 * areas_past_hol$Low / total_study_area_hol))
cat(sprintf("  Moderate:    %8.0f km²  (%5.1f%%)\n",
            areas_past_hol$Moderate,
            100 * areas_past_hol$Moderate / total_study_area_hol))
cat(sprintf("  High:        %8.0f km²  (%5.1f%%)\n",
            areas_past_hol$High,
            100 * areas_past_hol$High / total_study_area_hol))
cat(sprintf("  ─────────────────────────────\n"))
cat(sprintf("  Total suitable: %5.0f km²  (%5.1f%%)\n",
            suit_hol, 100 * suit_hol / total_study_area_hol))

# --- Niche displacement comparison ---
cat("\n--- Niche displacement: Holocene vs present (past-only model) ---\n")
cat(sprintf("                   Holocene    Present    Change\n"))
cat(sprintf("  Low:           %7.0f km²  %7.0f km²  %+.1f%%\n",
            areas_past_hol$Low, areas_past_present$Low,
            ifelse(areas_past_hol$Low > 0,
                   100 * (areas_past_present$Low - areas_past_hol$Low) / areas_past_hol$Low,
                   NA)))
cat(sprintf("  Moderate:      %7.0f km²  %7.0f km²  %+.1f%%\n",
            areas_past_hol$Moderate, areas_past_present$Moderate,
            ifelse(areas_past_hol$Moderate > 0,
                   100 * (areas_past_present$Moderate - areas_past_hol$Moderate) / areas_past_hol$Moderate,
                   NA)))
cat(sprintf("  High:          %7.0f km²  %7.0f km²  %+.1f%%\n",
            areas_past_hol$High, areas_past_present$High,
            ifelse(areas_past_hol$High > 0,
                   100 * (areas_past_present$High - areas_past_hol$High) / areas_past_hol$High,
                   NA)))
cat(sprintf("  Total suitable:%7.0f km²  %7.0f km²  %+.1f%%\n",
            suit_hol, suit_present,
            100 * (suit_present - suit_hol) / suit_hol))

# Save comparison table to CSV
past_model_area_comparison <- data.frame(
  Period       = c("Holocene (7–4 kya)", "Present day"),
  Unsuitable   = round(c(areas_past_hol$Unsuitable, areas_past_present$Unsuitable)),
  Low          = round(c(areas_past_hol$Low, areas_past_present$Low)),
  Moderate     = round(c(areas_past_hol$Moderate, areas_past_present$Moderate)),
  High         = round(c(areas_past_hol$High, areas_past_present$High)),
  Total_suitable = round(c(suit_hol, suit_present)),
  Pct_of_study = round(c(100 * suit_hol / total_study_area_hol,
                         100 * suit_present / total_study_area), 1)
)
write.csv(past_model_area_comparison,
          "output/tables/past_model_area_by_bin.csv", row.names = FALSE)
cat("\nSaved: output/tables/past_model_area_by_bin.csv\n")

# Zoomed map of Cuba / Zapata region
cuba_sf <- wi_sf[wi_sf$name %in% c("Cuba"), ]
if (nrow(cuba_sf) == 0) {
  cuba_sf <- wi_sf[grepl("Cuba", wi_sf$name, ignore.case = TRUE), ]
}

if (nrow(cuba_sf) > 0) {
  cuba_ext <- raster::extent(as_Spatial(cuba_sf))
  pred_cuba <- raster::crop(pred_present_from_past, cuba_ext)
  pred_cuba <- raster::mask(pred_cuba, as_Spatial(cuba_sf))
  
  png("output/past_model_present_cuba_zoom.png", width = 2400, height = 1200, res = 200)
  par(mar = c(4, 4, 3, 5))
  
  plot(st_geometry(cuba_sf), col = "grey90", border = "black",
       main = "C. rhombifer — Present Suitability at Cuba (Past-Only Model)",
       axes = TRUE)
  plot(pred_cuba, zlim = c(0, 1), col = rev(terrain.colors(50)),
       add = TRUE, legend = TRUE)
  plot(st_geometry(cuba_sf), add = TRUE, border = "black", lwd = 0.8)
  
  points(modern_thin$Longitude, modern_thin$Latitude,
         pch = 21, bg = "red", col = "black", cex = 1.2)
  legend("bottomleft",
         legend = c("Modern records", paste0("MaxSSS = ", round(thr_past, 3))),
         pch = c(21, NA), pt.bg = c("red", NA),
         col = c("black", NA), cex = 0.8, bg = "white")
  
  dev.off()
  cat("Saved: output/past_model_present_cuba_zoom.png\n")
}

# Save Zapata suitability analysis
zapata_analysis <- data.frame(
  metric = c("N modern records",
             "Mean suitability (past model, present climate)",
             "Median suitability (past model, present climate)",
             "Min suitability",
             "Max suitability",
             "Past-only MaxSSS threshold",
             "Records above threshold",
             "Percentile rank of Zapata median among suitable cells"),
  value = c(length(zapata_suit_past_model),
            round(mean(zapata_suit_past_model, na.rm = TRUE), 4),
            round(median(zapata_suit_past_model, na.rm = TRUE), 4),
            round(min(zapata_suit_past_model, na.rm = TRUE), 4),
            round(max(zapata_suit_past_model, na.rm = TRUE), 4),
            round(thr_past, 4),
            paste0(sum(zapata_suit_past_model >= thr_past, na.rm = TRUE),
                   " of ", sum(!is.na(zapata_suit_past_model))),
            paste0(zapata_pctile, "%"))
)
write.csv(zapata_analysis, "output/tables/zapata_suitability_past_model.csv",
          row.names = FALSE)
cat("Saved: output/tables/zapata_suitability_past_model.csv\n\n")


models <- list(
  Modern_present          = list(r = final_suitability,          thr = thr_modern),
  Past_7_4kya             = list(r = suitability_past,           thr = thr_past),
  Past_present            = list(r = pred_present_from_past,     thr = thr_past),
  Combined_present        = list(r = pred_present_from_combined, thr = thr_combined),
  
  Future_SSP245_Modern    = list(r = future_pred_ssp245,         thr = thr_modern),
  Future_SSP585_Modern    = list(r = future_pred_ssp585,         thr = thr_modern),
  
  Future_SSP245_Combined  = list(r = future_pred_ssp245_comb,    thr = thr_combined),
  Future_SSP585_Combined  = list(r = future_pred_ssp585_comb,    thr = thr_combined)
)

area_results <- bind_rows(lapply(names(models), function(nm) {
  cat("Processing area bins for:", nm, "\n")
  obj  <- models[[nm]]
  bins <- make_bins(obj$thr)
  out  <- calc_habitat_areas(obj$r, bins)
  out$Model <- nm
  out
})) %>%
  dplyr::mutate(Total_suitable = Low + Moderate + High) %>%
  dplyr::select(Model, `MaxSSS Threshold` = MaxSSS,
                Unsuitable, Low, Moderate, High, `Total Suitable` = Total_suitable)

cat("\n--- Appendix S4: Area (km²) of suitable habitat by MaxSSS bin ---\n")
print(area_results, row.names = FALSE)

write.csv(area_results, "output/area_results_MaxSSS_bins.csv", row.names = FALSE)

# Write formatted Appendix S4 xlsx
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Area Results")
  openxlsx::addWorksheet(wb, "MaxSSS Habitat Suitability Maps")
  
  # Write data starting at row 1
  openxlsx::writeData(wb, "Area Results", area_results, startRow = 1)
  
  # Caption row beneath the data
  caption_row <- nrow(area_results) + 2
  openxlsx::writeData(wb, "Area Results",
                      "Appendix S4. Area (km2) of suitable habitat predicted by species distribution models using MaxSSS thresholds",
                      startCol = 1, startRow = caption_row)
  
  # Header formatting
  header_style <- openxlsx::createStyle(textDecoration = "bold", halign = "center")
  openxlsx::addStyle(wb, "Area Results", header_style,
                     rows = 1, cols = 1:ncol(area_results), gridExpand = TRUE)
  
  # Number formatting for area columns (no decimals would lose precision; keep 2dp)
  num_style <- openxlsx::createStyle(numFmt = "#,##0.00")
  openxlsx::addStyle(wb, "Area Results", num_style,
                     rows = 2:(nrow(area_results) + 1),
                     cols = 3:ncol(area_results), gridExpand = TRUE)
  
  # Auto column widths
  openxlsx::setColWidths(wb, "Area Results",
                         cols = 1:ncol(area_results), widths = "auto")
  
  openxlsx::saveWorkbook(wb, "output/Appendix_S4.xlsx", overwrite = TRUE)
  cat("Saved: output/Appendix_S4.xlsx\n")
} else {
  cat("Note: openxlsx not available — Appendix S4 saved as CSV only.\n")
}


# ============================================================
# 14b) CLASSIFIED HABITAT SUITABILITY MAPS (MaxSSS-binned)
# ============================================================
# Reclassify each prediction surface into 4 habitat classes using
# the model-specific MaxSSS threshold, then save as individual PNGs.

# Colour scheme: grey = unsuitable, yellow = low, orange = moderate, red = high
habitat_cols  <- c("grey85", "#FFD700", "#FF8C00", "#CC0000")
habitat_labs  <- c("Unsuitable", "Low", "Moderate", "High")

classify_habitat <- function(rast, thr) {
  bins <- make_bins(thr)
  # Reclassify: 0 = unsuitable, 1 = low, 2 = moderate, 3 = high
  rcl <- matrix(c(
    -Inf,              bins$thr_unsuitable, 0,
    bins$thr_unsuitable, bins$thr_poor,     1,
    bins$thr_poor,       bins$thr_moderate, 2,
    bins$thr_moderate,   Inf,               3
  ), ncol = 3, byrow = TRUE)
  raster::reclassify(rast, rcl)
}

plot_habitat_map <- function(classified_rast, title, outfile,
                             width = 2400, height = 1600, res = 200) {
  png(outfile, width = width, height = height, res = res)
  
  # Get raster extent for consistent plot bounds
  e <- raster::extent(classified_rast)
  
  # Plot country outlines FIRST to establish coordinate system
  plot(st_geometry(wi_sf),
       xlim = c(e@xmin, e@xmax),
       ylim = c(e@ymin, e@ymax),
       col  = NA, border = "grey50", lwd = 0.3,
       main = title, axes = TRUE)
  
  # Add classified raster on top
  plot(classified_rast,
       col    = habitat_cols,
       breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5),
       legend = FALSE,
       add    = TRUE)
  
  # Re-draw outlines on top of raster for visibility
  plot(st_geometry(wi_sf), add = TRUE, border = "black", lwd = 0.5)
  
  # Add legend
  legend("bottomleft",
         legend = habitat_labs,
         fill   = habitat_cols,
         title  = "Habitat suitability",
         cex    = 0.8,
         bg     = "white")
  
  dev.off()
  cat("Saved:", outfile, "\n")
}

# Plot to screen + save to file for each model
for (nm in names(models)) {
  obj <- models[[nm]]
  classified <- classify_habitat(obj$r, obj$thr)
  
  # Clean title
  title <- paste("C. rhombifer —", gsub("_", " ", nm))
  
  # Save to file
  outfile <- file.path("output/habitat_maps",
                       paste0("habitat_", nm, ".png"))
  plot_habitat_map(classified, title, outfile)
  
  # Also plot to screen
  e <- raster::extent(classified)
  plot(st_geometry(wi_sf),
       xlim = c(e@xmin, e@xmax),
       ylim = c(e@ymin, e@ymax),
       col  = NA, border = "grey50", lwd = 0.3,
       main = title, axes = TRUE)
  plot(classified,
       col    = habitat_cols,
       breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5),
       legend = FALSE,
       add    = TRUE)
  plot(st_geometry(wi_sf), add = TRUE, border = "black", lwd = 0.5)
  legend("bottomleft", legend = habitat_labs, fill = habitat_cols,
         title = "Habitat suitability", cex = 0.8, bg = "white")
}

cat("\nAll habitat maps saved to output/habitat_maps/\n")
list.files("output/habitat_maps", pattern = "\\.png$")


# ============================================================
# 15) ENVIRONMENTAL VARIABLE RANGES (suitable vs unsuitable)
# ============================================================

model_env_list <- list(
  Modern_present         = list(pred = final_suitability,          env = final_env,              thr = thr_modern),
  Past_7_4kya            = list(pred = suitability_past,           env = past_env_aligned,       thr = thr_past),
  Combined_present       = list(pred = pred_present_from_combined, env = ENV_present_pred,       thr = thr_combined),
  
  Future_SSP245_Modern   = list(pred = future_pred_ssp245,         env = future_env_ssp245,      thr = thr_modern),
  Future_SSP585_Modern   = list(pred = future_pred_ssp585,         env = future_env_ssp585,      thr = thr_modern),
  
  Future_SSP245_Combined = list(pred = future_pred_ssp245_comb,    env = future_env_ssp245_comb, thr = thr_combined),
  Future_SSP585_Combined = list(pred = future_pred_ssp585_comb,    env = future_env_ssp585_comb, thr = thr_combined)
)

extract_ranges <- function(pred_raster, env_stack, thr, which_side = c("suitable", "unsuitable")) {
  which_side <- match.arg(which_side)
  
  pv <- raster::getValues(pred_raster)
  if (which_side == "suitable")   idx <- which(is.finite(pv) & pv >= thr)
  if (which_side == "unsuitable") idx <- which(is.finite(pv) & pv <  thr)
  
  if (length(idx) == 0) return(NULL)
  
  out <- vector("list", raster::nlayers(env_stack))
  for (i in seq_len(raster::nlayers(env_stack))) {
    lyr <- env_stack[[i]]
    vnm <- names(lyr)
    lv  <- raster::getValues(lyr)[idx]
    lv  <- lv[is.finite(lv)]
    if (length(lv) == 0) next
    out[[i]] <- data.frame(Variable = vnm, Min = min(lv), Max = max(lv),
                           stringsAsFactors = FALSE)
  }
  bind_rows(out)
}

range_results_suitable <- bind_rows(lapply(names(model_env_list), function(nm) {
  cat("Suitable ranges for:", nm, "\n")
  obj <- model_env_list[[nm]]
  rr  <- extract_ranges(obj$pred, obj$env, obj$thr, "suitable")
  if (is.null(rr) || nrow(rr) == 0) return(NULL)
  rr$Model <- nm; rr$MaxSSS <- obj$thr; rr
})) %>% dplyr::select(Model, MaxSSS, Variable, Min, Max)

print(range_results_suitable)
write.csv(range_results_suitable, "output/env_ranges_suitable.csv", row.names = FALSE)

range_results_unsuitable <- bind_rows(lapply(names(model_env_list), function(nm) {
  cat("Unsuitable ranges for:", nm, "\n")
  obj <- model_env_list[[nm]]
  rr  <- extract_ranges(obj$pred, obj$env, obj$thr, "unsuitable")
  if (is.null(rr) || nrow(rr) == 0) return(NULL)
  rr$Model <- nm; rr$MaxSSS <- obj$thr; rr
})) %>% dplyr::select(Model, MaxSSS, Variable, Min, Max)

print(range_results_unsuitable)
write.csv(range_results_unsuitable, "output/env_ranges_unsuitable.csv", row.names = FALSE)


# ============================================================
# 16) PERMUTATION IMPORTANCE
# ============================================================

align_X_to_maxnet <- function(model, X) {
  X <- as.data.frame(X)
  req <- names(model$varmin)
  if (is.null(req) || length(req) == 0) {
    stop("Could not determine required predictors from model$varmin.")
  }
  missing <- setdiff(req, names(X))
  extra   <- setdiff(names(X), req)
  if (length(missing) > 0) stop("X missing: ", paste(missing, collapse = ", "))
  if (length(extra) > 0)   message("Dropping: ", paste(extra, collapse = ", "))
  X[, req, drop = FALSE]
}

perm_importance_auc_pct <- function(model, X, y, nperm = 10, seed = 1) {
  y <- as.integer(as.character(y))
  y[!(y %in% c(0L, 1L))] <- NA_integer_
  
  X <- align_X_to_maxnet(model, X)
  
  keep <- complete.cases(X) & !is.na(y)
  X <- X[keep, , drop = FALSE]
  y <- y[keep]
  
  stopifnot(nrow(X) == length(y))
  if (length(unique(y)) < 2) stop("Both classes (0/1) required for AUC.")
  
  base_pred <- predict(model, X, type = "cloglog")
  base_auc  <- dismo::evaluate(p = base_pred[y == 1], a = base_pred[y == 0])@auc
  
  vars  <- colnames(X)
  drops <- numeric(length(vars))
  
  set.seed(seed)
  for (i in seq_along(vars)) {
    v <- vars[i]
    auc_perm <- numeric(nperm)
    for (p in seq_len(nperm)) {
      Xp <- X
      Xp[[v]] <- sample(Xp[[v]])
      pred_p <- predict(model, Xp, type = "cloglog")
      auc_perm[p] <- dismo::evaluate(p = pred_p[y == 1], a = pred_p[y == 0])@auc
    }
    drops[i] <- mean(base_auc - auc_perm, na.rm = TRUE)
  }
  
  tibble(Variable = vars, AUC_base = base_auc, AUC_drop = drops) %>%
    mutate(
      Contribution_pct = if (sum(AUC_drop, na.rm = TRUE) > 0) {
        100 * AUC_drop / sum(AUC_drop, na.rm = TRUE)
      } else { 0 }
    ) %>%
    arrange(desc(Contribution_pct))
}

model_objects <- list(
  Modern   = list(model = me,         X = pb_climate,  y = pb),
  Past     = list(model = me_past,    X = pb_env_past, y = pb_past),
  Combined = list(model = me_combined, X = X,          y = y)
)

varimp_all <- bind_rows(lapply(names(model_objects), function(nm) {
  cat("Permutation importance:", nm, "\n")
  obj <- model_objects[[nm]]
  imp <- perm_importance_auc_pct(
    model = obj$model, X = obj$X, y = obj$y,
    nperm = 10, seed = 1
  )
  imp %>% mutate(Model = nm) %>%
    dplyr::select(Model, Variable, AUC_base, AUC_drop, Contribution_pct)
}))

print(n = 30, varimp_all)
write.csv(varimp_all, "output/variable_contributions_permutation.csv", row.names = FALSE)


# ============================================================
# 17) RESPONSE CURVES
# ============================================================

outdir <- file.path(getwd(), "output", "response_curves")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

safe_name <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

save_response_curves_model <- function(model, X, model_name, outdir,
                                       ncol = 2, n = 200,
                                       qlow = 0.01, qhigh = 0.99) {
  X <- as.data.frame(X)
  is_num <- vapply(X, is.numeric, logical(1))
  X <- X[, is_num, drop = FALSE]
  if (ncol(X) == 0) stop("No numeric predictors for: ", model_name)
  
  req <- names(model$varmin)
  if (!is.null(req) && length(req) > 0) {
    missing <- setdiff(req, names(X))
    if (length(missing) > 0)
      stop("X missing predictors for ", model_name, ": ", paste(missing, collapse = ", "))
    X <- X[, req, drop = FALSE]
  }
  
  typical <- vapply(X, function(z) {
    z <- z[is.finite(z)]
    if (length(z) == 0) NA_real_ else median(z, na.rm = TRUE)
  }, numeric(1))
  
  # Only plot continuous variables — exclude geology dummies (binary 0/1)
  # Geo dummies are held at their median value during prediction but
  # don't get their own response curve panel.
  vars_actual <- names(typical)[is.finite(typical)]
  vars_actual <- vars_actual[!grepl("^geo_", vars_actual)]
  if (length(vars_actual) == 0) stop("No plottable predictors in: ", model_name)
  
  typical_df <- as.data.frame(as.list(typical))[1, , drop = FALSE]
  
  n_panels <- length(vars_actual)
  nrow_panels <- ceiling(n_panels / ncol)
  
  outfile <- file.path(outdir, paste0("response_curves_", safe_name(model_name), ".png"))
  
  png(outfile,
      width  = 2400,
      height = 1200 + 450 * max(0, (nrow_panels - 1)),
      res    = 200)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit({ par(oldpar); dev.off() }, add = TRUE)
  
  par(
    mfrow = c(nrow_panels, ncol),
    mar   = c(7.5, 7.5, 3, 2),
    oma   = c(0, 0, 4, 0),
    cex.lab  = 1.15,
    cex.axis = 0.9
  )
  
  for (i in seq_along(vars_actual)) {
    v <- vars_actual[i]
    xv <- X[[v]]
    xv <- xv[is.finite(xv)]
    
    if (length(xv) < 2) {
      plot.new(); title(main = v)
      text(0.5, 0.5, "Not enough finite values", cex = 0.95)
      next
    }
    
    lo <- as.numeric(quantile(xv, qlow, na.rm = TRUE, names = FALSE))
    hi <- as.numeric(quantile(xv, qhigh, na.rm = TRUE, names = FALSE))
    
    if (!is.finite(lo) || !is.finite(hi) || lo == hi) {
      plot.new(); title(main = v)
      text(0.5, 0.5, "Insufficient variation", cex = 0.95)
      next
    }
    
    xseq <- seq(lo, hi, length.out = n)
    newdata <- typical_df[rep(1, n), , drop = FALSE]
    newdata[[v]] <- xseq
    
    yhat <- tryCatch(
      predict(model, newdata, type = "cloglog"),
      error = function(e) rep(NA_real_, length(xseq))
    )
    
    keep <- is.finite(xseq) & is.finite(yhat)
    if (sum(keep) < 2) {
      plot.new(); title(main = v)
      text(0.5, 0.5, "Prediction failed", cex = 0.95)
      next
    }
    
    plot(xseq[keep], yhat[keep],
         type = "l", ylim = c(0, 1),
         xlab = v, ylab = "Predicted suitability (cloglog)")
    title(main = v, cex.main = 1.05)
  }
  
  mtext(
    paste("C. rhombifer —", model_name, "SDM response curves (Maxnet cloglog)"),
    outer = TRUE, cex = 1.35
  )
  
  cat("Saved:", outfile, "\n")
}

save_response_curves_model(me,         pb_climate,  "Modern",   outdir)
save_response_curves_model(me_past,    pb_env_past, "Past",     outdir)
save_response_curves_model(me_combined, X,          "Combined", outdir)


# ============================================================
# 18) MODERN DISTRIBUTION SUITABILITY ANALYSIS
# ============================================================

R_potential_present <- pred_present_from_combined
thr <- thr_combined

stopifnot(is.finite(thr), thr >= 0, thr <= 1)

modern_pts <- modern_thin %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(R_potential_present))

s_mod <- raster::extract(R_potential_present, modern_pts)
s_mod <- s_mod[is.finite(s_mod)]

pot_vals <- raster::getValues(R_potential_present)
pot_vals <- pot_vals[is.finite(pot_vals) & pot_vals >= thr]

if (length(pot_vals) == 0) stop("No cells >= MaxSSS threshold; thr=", thr)
if (length(s_mod) == 0)    stop("No finite suitability at modern points.")

pct_rank <- sapply(s_mod, function(x) mean(pot_vals <= x))

summary_tbl <- tibble(
  metric = c("MaxSSS threshold (combined)",
             "Potential niche median (cells>=thr)",
             "Modern points median (combined→present)",
             "Modern points mean percentile rank",
             "Prop modern points in bottom 25% of niche",
             "Prop modern points in top 25% of niche"),
  value = c(thr,
            median(pot_vals),
            median(s_mod),
            mean(pct_rank),
            mean(pct_rank <= 0.25),
            mean(pct_rank >= 0.75))
)

print(summary_tbl)
cat("\nPercentile ranks (modern points within potential niche):\n")
print(summary(pct_rank))


# ============================================================
# 18b) CONTIGUOUS SUITABLE HABITAT PATCH ANALYSIS
# ============================================================
# Identify contiguous patches of suitable habitat (cells >= MaxSSS),
# calculate their areas, and determine which patch the modern
# C. rhombifer population occupies.

cat("\n=== CONTIGUOUS PATCH ANALYSIS ===\n")
cat("Using combined model projected onto present climate\n")
cat("MaxSSS threshold:", round(thr, 4), "\n\n")

# Binarize: 1 = suitable (>= threshold), NA = unsuitable
suitable_binary <- R_potential_present >= thr
suitable_binary[suitable_binary == 0] <- NA

# Identify contiguous patches (requires igraph for queen's case adjacency)
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
patches <- raster::clump(suitable_binary, directions = 8)  # 8-connected (queen's case)

n_patches <- max(values(patches), na.rm = TRUE)
cat("Total contiguous patches of suitable habitat:", n_patches, "\n")

# Calculate area of each patch (km²)
cell_areas <- raster::area(patches)  # km² per cell for lon/lat rasters
patch_ids  <- raster::getValues(patches)
cell_km2   <- raster::getValues(cell_areas)

patch_area_df <- data.frame(patch_id = patch_ids, area_km2 = cell_km2) %>%
  filter(!is.na(patch_id)) %>%
  group_by(patch_id) %>%
  summarise(
    area_km2 = sum(area_km2, na.rm = TRUE),
    n_cells  = n(),
    .groups  = "drop"
  ) %>%
  arrange(desc(area_km2)) %>%
  mutate(rank = row_number())

cat("\nTop 10 largest patches:\n")
print(head(patch_area_df, 10))

# Identify which patch the modern records occupy
modern_patch_ids <- raster::extract(patches, modern_pts)
modern_patch_ids <- modern_patch_ids[!is.na(modern_patch_ids)]

if (length(modern_patch_ids) == 0) {
  cat("\nWARNING: No modern records fall within suitable habitat patches.\n")
  cat("This may indicate records fall in cells just below the threshold.\n")
} else {
  # Most common patch (should be unanimous for clustered records)
  occupied_patch <- as.integer(names(sort(table(modern_patch_ids), decreasing = TRUE))[1])
  n_in_patch <- sum(modern_patch_ids == occupied_patch)
  
  occupied_area <- patch_area_df %>%
    filter(patch_id == occupied_patch)
  
  cat("\n--- Modern population patch ---\n")
  cat("Patch ID:", occupied_patch, "\n")
  cat("Modern records in this patch:", n_in_patch, "of", length(modern_patch_ids), "\n")
  cat("Patch area:", format(round(occupied_area$area_km2, 1), big.mark = ","), "km²\n")
  cat("Patch rank by area:", occupied_area$rank, "of", n_patches, "\n")
  
  # How does the occupied patch compare?
  if (occupied_area$rank == 1) {
    cat("The modern population occupies the LARGEST contiguous area of suitable habitat.\n")
  } else {
    largest <- patch_area_df %>% filter(rank == 1)
    cat("Largest patch:", format(round(largest$area_km2, 1), big.mark = ","), "km²\n")
  }
  
  # Other large patches (potential reintroduction targets)
  cat("\n--- Other contiguous patches > 1,000 km² ---\n")
  large_patches <- patch_area_df %>%
    filter(area_km2 >= 1000, patch_id != occupied_patch)
  
  if (nrow(large_patches) > 0) {
    print(large_patches)
  } else {
    cat("No other patches exceed 1,000 km²\n")
  }
}

# Save patch analysis
write.csv(patch_area_df, "output/contiguous_patches_combined_present.csv",
          row.names = FALSE)
cat("\nPatch analysis saved to output/contiguous_patches_combined_present.csv\n")

# --- Map of top 10 largest patches ---
top10_ids <- patch_area_df %>% head(10) %>% pull(patch_id)

# Create a raster showing only the top 10 patches, coded 1–10 by rank
patch_top10 <- patches
patch_vals  <- raster::getValues(patch_top10)

# Map patch IDs to their rank (1 = largest)
rank_lookup <- setNames(1:10, top10_ids)
new_vals <- rank_lookup[as.character(patch_vals)]
patch_top10 <- raster::setValues(patch_top10, as.integer(new_vals))

# Colour palette: distinct colours, with the occupied patch highlighted
patch_colours <- c("#CC0000",   # 1 - dark red
                   "#FF6600",   # 2 - orange
                   "#FFD700",   # 3 - gold
                   "#66CC00",   # 4 - green
                   "#009999",   # 5 - teal
                   "#3366CC",   # 6 - blue
                   "#6633CC",   # 7 - purple
                   "#CC3399",   # 8 - magenta
                   "#996633",   # 9 - brown
                   "#999999")   # 10 - grey

# Build legend labels with area
top10_labels <- patch_area_df %>%
  head(10) %>%
  mutate(label = paste0("#", rank, ": ",
                        format(round(area_km2), big.mark = ","), " km²")) %>%
  pull(label)

# Mark the occupied patch in the legend
if (exists("occupied_patch")) {
  occ_rank <- rank_lookup[as.character(occupied_patch)]
  if (!is.na(occ_rank)) {
    top10_labels[occ_rank] <- paste0(top10_labels[occ_rank], " *")
  }
}

# Plot function
plot_patches <- function(outfile = NULL) {
  if (!is.null(outfile)) {
    png(outfile, width = 2800, height = 1800, res = 200)
  }
  
  e <- raster::extent(patch_top10)
  
  # Country outlines first
  plot(st_geometry(wi_sf),
       xlim = c(e@xmin, e@xmax), ylim = c(e@ymin, e@ymax),
       col = "grey95", border = "grey60", lwd = 0.3,
       main = "C. rhombifer — Top 10 contiguous suitable habitat patches\n(Combined model, present climate)",
       axes = TRUE, cex.main = 1.0)
  
  # Add patches
  plot(patch_top10,
       col    = patch_colours,
       breaks = seq(0.5, 10.5, by = 1),
       legend = FALSE,
       add    = TRUE)
  
  # Re-draw outlines
  plot(st_geometry(wi_sf), add = TRUE, border = "black", lwd = 0.5)
  
  # Add modern occurrence points
  points(modern_thin$Longitude, modern_thin$Latitude,
         pch = 21, bg = "white", col = "black", cex = 1.2, lwd = 1.5)
  
  # Legend
  legend("bottomleft",
         legend = top10_labels,
         fill   = patch_colours,
         title  = "Patch rank (* = occupied)",
         cex    = 0.65,
         bg     = "white",
         ncol   = 1)
  
  if (!is.null(outfile)) {
    dev.off()
    cat("Saved:", outfile, "\n")
  }
}

# Save to file
plot_patches("output/habitat_maps/top10_patches_combined_present.png")

# Plot to screen
plot_patches()


# ============================================================
# 18c) THRESHOLD SENSITIVITY ANALYSIS
# ============================================================
# Compare habitat area estimates for the Combined model under
# three alternative MaxSSS thresholds:
#   1. Past-only model's MaxSSS (primary analysis threshold)
#   2. Combined model evaluated on pooled PB (own threshold)
#   3. Combined model evaluated on past PB only
# This addresses reviewer concerns about applying one model's
# threshold to another model.
# ============================================================

cat("\n\n========== THRESHOLD SENSITIVITY ANALYSIS ==========\n")

# Extract sensitivity thresholds from eval_results
# (get_mss defined in section 14 above)
thr_primary     <- thr_past
thr_pooled_pb   <- get_mss("Combined (pooled PB)")
thr_past_pb     <- get_mss("Combined (past PB)")

cat("Threshold 1 (Past-only model MaxSSS, primary):     ", round(thr_primary,   4), "\n")
cat("Threshold 2 (Combined model, pooled PB MaxSSS):    ", round(thr_pooled_pb,  4), "\n")
cat("Threshold 3 (Combined model, past PB only MaxSSS): ", round(thr_past_pb,    4), "\n\n")

# Surfaces to evaluate: combined model present + two futures
sensitivity_surfaces <- list(
  Combined_present       = pred_present_from_combined,
  Combined_SSP245        = future_pred_ssp245_comb,
  Combined_SSP585        = future_pred_ssp585_comb
)

threshold_labels <- c(
  "Past-only model MaxSSS (primary)" = thr_primary,
  "Combined pooled PB MaxSSS"        = thr_pooled_pb,
  "Combined past PB MaxSSS"          = thr_past_pb
)

# Calculate area under each threshold x surface combination
sensitivity_rows <- list()
for (surf_name in names(sensitivity_surfaces)) {
  r <- sensitivity_surfaces[[surf_name]]
  cell_areas <- raster::area(r)
  vals  <- raster::getValues(r)
  areas <- raster::getValues(cell_areas)
  valid <- is.finite(vals) & is.finite(areas)
  vals  <- vals[valid]
  areas <- areas[valid]
  
  for (thr_label in names(threshold_labels)) {
    thr_val <- threshold_labels[[thr_label]]
    
    suitable_area <- sum(areas[vals >= thr_val])
    total_area    <- sum(areas)
    pct_suitable  <- round(100 * suitable_area / total_area, 1)
    
    # Also compute binned areas
    bins <- make_bins(thr_val)
    low_area  <- sum(areas[vals >= bins$thr_unsuitable & vals < bins$thr_poor])
    mod_area  <- sum(areas[vals >= bins$thr_poor       & vals < bins$thr_moderate])
    high_area <- sum(areas[vals >= bins$thr_moderate])
    
    sensitivity_rows[[length(sensitivity_rows) + 1]] <- data.frame(
      Surface          = surf_name,
      Threshold_Label  = thr_label,
      Threshold_Value  = round(thr_val, 4),
      Total_Suitable_km2 = round(suitable_area, 1),
      Pct_Suitable     = pct_suitable,
      Low_km2          = round(low_area, 1),
      Moderate_km2     = round(mod_area, 1),
      High_km2         = round(high_area, 1),
      stringsAsFactors = FALSE
    )
  }
}

sensitivity_df <- do.call(rbind, sensitivity_rows)
print(sensitivity_df)

write.csv(sensitivity_df, "output/tables/threshold_sensitivity_analysis.csv",
          row.names = FALSE)
cat("\nSaved: output/tables/threshold_sensitivity_analysis.csv\n")

# --- Percentage change relative to primary threshold ---
cat("\n--- Percentage change in total suitable area vs primary threshold ---\n")
for (surf_name in names(sensitivity_surfaces)) {
  primary_area <- sensitivity_df$Total_Suitable_km2[
    sensitivity_df$Surface == surf_name &
      sensitivity_df$Threshold_Label == "Past-only model MaxSSS (primary)"]
  
  for (thr_label in names(threshold_labels)) {
    if (thr_label == "Past-only model MaxSSS (primary)") next
    alt_area <- sensitivity_df$Total_Suitable_km2[
      sensitivity_df$Surface == surf_name &
        sensitivity_df$Threshold_Label == thr_label]
    pct_change <- round(100 * (alt_area - primary_area) / primary_area, 1)
    cat(sprintf("  %s | %s: %+.1f%%\n", surf_name, thr_label, pct_change))
  }
}

# --- Summary figure: grouped bar chart ---
cat("\nGenerating sensitivity analysis figure...\n")

png("output/figures/threshold_sensitivity_barplot.png",
    width = 2400, height = 1200, res = 200)

# Reshape for plotting
sens_plot <- sensitivity_df
sens_plot$Surface <- factor(sens_plot$Surface,
                            levels = c("Combined_present", "Combined_SSP245", "Combined_SSP585"),
                            labels = c("Present", "SSP2-4.5\n(2061-2080)", "SSP5-8.5\n(2061-2080)"))
sens_plot$Threshold_Label <- factor(sens_plot$Threshold_Label,
                                    levels = c("Past-only model MaxSSS (primary)",
                                               "Combined pooled PB MaxSSS",
                                               "Combined past PB MaxSSS"))

# Grouped barplot
bar_data <- tapply(sens_plot$Total_Suitable_km2,
                   list(sens_plot$Threshold_Label, sens_plot$Surface),
                   identity)

bar_cols <- c("#2166AC", "#B2182B", "#F4A582")
bp <- barplot(bar_data, beside = TRUE,
              col = bar_cols,
              main = "Threshold Sensitivity Analysis — Combined Model",
              ylab = expression("Total suitable area (km"^2*")"),
              las = 1, ylim = c(0, max(sensitivity_df$Total_Suitable_km2) * 1.15),
              border = NA, cex.main = 1.1)

# Add value labels above bars
text(bp, bar_data, labels = format(round(bar_data), big.mark = ","),
     pos = 3, cex = 0.65, font = 1)

legend("topright",
       legend = c(paste0("Past-only MaxSSS (", round(thr_primary, 3), ")"),
                  paste0("Combined pooled PB (", round(thr_pooled_pb, 3), ")"),
                  paste0("Combined past PB (", round(thr_past_pb, 3), ")")),
       fill = bar_cols, border = NA, bty = "n", cex = 0.75)

dev.off()
cat("Saved: output/figures/threshold_sensitivity_barplot.png\n")


# ============================================================
# 18d) HUMAN FOOTPRINT INDEX OVERLAY (Combined model, present)
# ============================================================
# Overlay the Human Footprint Index (HFP-100 v1.2; Gassert et al., 2023)
# on the combined model's present-day suitability surface to visualise the
# relationship between habitat quality and anthropogenic pressure.
#
# HFP-100 v1.2: 100m resolution, scores 0–50 (0 = no pressure, 50 = severe)
# Download from: https://source.coop/vizzuality/hfp-100
# Raw values are x1000 (uint16); 65535 = no data
# Place the GeoTIFF in data/anthropogenic/
# ============================================================

hfi_path <- "data/anthropogenic/hfp_2020_100m_v1-2_cog.tif"

if (file.exists(hfi_path)) {
  cat("\n\n========== HUMAN FOOTPRINT INDEX OVERLAY ==========\n")
  cat("Loading HFP-100 v1.2 (2020) and clipping to study area...\n")
  
  # Use terra for efficient windowed reading of the large COG file
  hfi_terra <- terra::rast(hfi_path)
  cat("HFP CRS:", terra::crs(hfi_terra, describe = TRUE)$name, "\n")
  
  # Create study area extent in the HFP's native CRS for efficient cropping
  wi_vect <- terra::vect(as_Spatial(wi_union))
  wi_reproj <- terra::project(wi_vect, hfi_terra)
  
  # Crop in native CRS first (fast), then reproject to WGS84
  hfi_cropped <- terra::crop(hfi_terra, terra::ext(wi_reproj) + 0.5)
  
  # Reproject to WGS84 to match SDM data
  hfi_wgs84 <- terra::project(hfi_cropped, "EPSG:4326")
  
  # Mask to study area (now both in WGS84)
  hfi_wi_terra <- terra::mask(hfi_wgs84, terra::vect(as_Spatial(wi_union)))
  
  # Convert from uint16 (x1000) to actual HFP score (0–50)
  hfi_wi_terra[hfi_wi_terra >= 65000] <- NA
  hfi_wi_terra <- hfi_wi_terra / 1000
  
  cat("HFP range in study area:", round(terra::minmax(hfi_wi_terra)[1], 2),
      "to", round(terra::minmax(hfi_wi_terra)[2], 2), "\n")
  
  # Convert to raster and resample to match SDM grid
  hfi_wi <- raster::raster(hfi_wi_terra)
  hfi_wi <- raster::resample(hfi_wi, pred_present_from_combined)
  
  # --- Side-by-side: combined present suitability vs HFP ---
  png("output/figures/suitability_vs_HFP_sidebyside.png",
      width = 3200, height = 1600, res = 200)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 5))
  
  plot(pred_present_from_combined,
       main = "Present Suitability\n(Combined Model)",
       zlim = c(0, 1), col = rev(terrain.colors(50)))
  plot(st_geometry(wi_sf), add = TRUE, border = "black", lwd = 0.5)
  
  plot(hfi_wi,
       main = "Human Footprint Index\n(HFP-100 v1.2, 2020)",
       col = colorRampPalette(c("darkgreen", "yellow", "orange", "red"))(50))
  plot(st_geometry(wi_sf), add = TRUE, border = "black", lwd = 0.5)
  
  dev.off()
  cat("Saved: output/figures/suitability_vs_HFP_sidebyside.png\n")
  
  # --- Suitable habitat in low human footprint areas (HFP < 6) ---
  # Di Marco et al. (2018) threshold
  hfp_threshold <- 6
  suit_low_hfp <- pred_present_from_combined
  
  # Mask out BOTH unsuitable habitat AND high-HFP areas
  suit_low_hfp[pred_present_from_combined < thr_combined] <- NA
  suit_low_hfp[hfi_wi >= hfp_threshold | is.na(hfi_wi)] <- NA
  
  # Calculate area of suitable + low-HFP habitat
  suit_low_hfp_binary <- suit_low_hfp >= thr_combined
  suit_low_hfp_binary[is.na(suit_low_hfp)] <- NA
  low_hfp_area <- raster::cellStats(
    raster::area(suit_low_hfp_binary) * suit_low_hfp_binary, sum, na.rm = TRUE
  )
  
  cat("\n--- Suitable habitat with HFP <", hfp_threshold, "---\n")
  cat("Area:", round(low_hfp_area), "km2\n")
  
  # Total suitable area for comparison
  suit_binary <- pred_present_from_combined >= thr_combined
  suit_binary[is.na(pred_present_from_combined)] <- NA
  total_suitable_area <- raster::cellStats(
    raster::area(suit_binary) * suit_binary, sum, na.rm = TRUE
  )
  cat("Total suitable area:", round(total_suitable_area), "km2\n")
  cat("Proportion with low HFP:", round(100 * low_hfp_area / total_suitable_area, 1), "%\n")
  
  # Map
  png("output/figures/suitable_habitat_low_HFP.png",
      width = 2400, height = 1600, res = 200)
  par(mar = c(4, 4, 3, 5))
  
  plot(st_geometry(wi_sf), col = "grey90", border = "grey50",
       main = paste0("C. rhombifer — Suitable habitat with low human pressure\n",
                     "(suitability >= ", round(thr_combined, 3),
                     " & HFP < ", hfp_threshold, ")"),
       axes = TRUE)
  plot(suit_low_hfp, col = rev(terrain.colors(50)), zlim = c(0, 1),
       add = TRUE, legend = TRUE)
  plot(st_geometry(wi_sf), add = TRUE, border = "black", lwd = 0.5)
  
  legend("topleft",
         legend = paste0("Low-HFP suitable: ", format(round(low_hfp_area), big.mark = ","),
                         " km2"),
         cex = 0.7, bg = "white")
  
  dev.off()
  cat("Saved: output/figures/suitable_habitat_low_HFP.png\n")
  
  # --- Quantitative comparison: HFP at modern range vs all suitable habitat ---
  hfi_at_modern <- raster::extract(hfi_wi, modern_thin[, c("Longitude", "Latitude")])
  suit_vals <- raster::getValues(pred_present_from_combined)
  hfi_vals <- raster::getValues(hfi_wi)
  valid <- is.finite(suit_vals) & is.finite(hfi_vals) & suit_vals >= thr_combined
  hfi_at_suitable <- hfi_vals[valid]
  
  cat("\n--- HFP at modern Zapata records ---\n")
  cat("Mean HFP at modern range:", round(mean(hfi_at_modern, na.rm = TRUE), 2), "\n")
  cat("Median HFP at modern range:", round(median(hfi_at_modern, na.rm = TRUE), 2), "\n")
  cat("Mean HFP across all suitable habitat:", round(mean(hfi_at_suitable, na.rm = TRUE), 2), "\n")
  cat("Median HFP across all suitable habitat:", round(median(hfi_at_suitable, na.rm = TRUE), 2), "\n")
  
  hfi_mean_modern <- mean(hfi_at_modern, na.rm = TRUE)
  hfi_mean_suitable <- mean(hfi_at_suitable, na.rm = TRUE)
  
  if (is.finite(hfi_mean_modern) && is.finite(hfi_mean_suitable)) {
    if (hfi_mean_modern < hfi_mean_suitable) {
      cat("=> Modern range has LOWER human pressure than average suitable habitat,\n")
      cat("   consistent with persistence in a low-disturbance refugium.\n")
    } else {
      cat("=> Modern range has HIGHER human pressure than average suitable habitat.\n")
    }
  } else {
    cat("WARNING: HFP extraction returned NAs — check CRS alignment.\n")
    print(hfi_at_modern)
  }
  
  # --- Cuba-focused HFP + suitability overlay ---
  if (nrow(cuba_sf) > 0) {
    cuba_ext <- raster::extent(as_Spatial(cuba_sf))
    hfi_cuba <- raster::crop(hfi_wi, cuba_ext)
    hfi_cuba <- raster::mask(hfi_cuba, as_Spatial(cuba_sf))
    suit_cuba <- raster::crop(pred_present_from_combined, cuba_ext)
    suit_cuba <- raster::mask(suit_cuba, as_Spatial(cuba_sf))
    
    png("output/figures/suitability_vs_HFP_cuba_zoom.png",
        width = 3200, height = 1200, res = 200)
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 5))
    
    plot(st_geometry(cuba_sf), col = "grey90", border = "black",
         main = "Present Suitability at Cuba\n(Combined Model)", axes = TRUE)
    plot(suit_cuba, zlim = c(0, 1), col = rev(terrain.colors(50)),
         add = TRUE, legend = TRUE)
    plot(st_geometry(cuba_sf), add = TRUE, border = "black", lwd = 0.8)
    points(modern_thin$Longitude, modern_thin$Latitude,
           pch = 21, bg = "red", col = "black", cex = 1.0)
    
    plot(st_geometry(cuba_sf), col = "grey90", border = "black",
         main = "Human Footprint Index at Cuba\n(HFP-100 v1.2, 2020)", axes = TRUE)
    plot(hfi_cuba,
         col = colorRampPalette(c("darkgreen", "yellow", "orange", "red"))(50),
         add = TRUE, legend = TRUE)
    plot(st_geometry(cuba_sf), add = TRUE, border = "black", lwd = 0.8)
    points(modern_thin$Longitude, modern_thin$Latitude,
           pch = 21, bg = "cyan", col = "black", cex = 1.0)
    
    dev.off()
    cat("Saved: output/figures/suitability_vs_HFP_cuba_zoom.png\n")
  }
  
  # Save HFP comparison table
  hfp_comparison <- data.frame(
    metric = c("Mean HFP at modern range (Zapata)",
               "Median HFP at modern range (Zapata)",
               "Mean HFP across all suitable habitat",
               "Median HFP across all suitable habitat",
               "HFP score range (study area)",
               "Total suitable area (km2)",
               "Suitable area with HFP < 6 (km2)",
               "% suitable with low HFP"),
    value = c(round(mean(hfi_at_modern, na.rm = TRUE), 2),
              round(median(hfi_at_modern, na.rm = TRUE), 2),
              round(mean(hfi_at_suitable, na.rm = TRUE), 2),
              round(median(hfi_at_suitable, na.rm = TRUE), 2),
              paste0(round(raster::cellStats(hfi_wi, min, na.rm = TRUE), 2), "–",
                     round(raster::cellStats(hfi_wi, max, na.rm = TRUE), 2)),
              round(total_suitable_area),
              round(low_hfp_area),
              round(100 * low_hfp_area / total_suitable_area, 1))
  )
  write.csv(hfp_comparison, "output/tables/HFP_comparison_modern_vs_suitable.csv",
            row.names = FALSE)
  cat("Saved: output/tables/HFP_comparison_modern_vs_suitable.csv\n")
  
} else {
  cat("\nHFP file not found at:", hfi_path, "\n")
  cat("Download HFP-100 v1.2 (2020) from: https://source.coop/vizzuality/hfp-100\n")
  cat("Place in data/anthropogenic/ and rerun this section.\n")
}


# ============================================================
# 19) SESSION INFO
# ============================================================
cat("\n\n========== SESSION INFO ==========\n")
sessionInfo()
cat("Script completed successfully.\n")