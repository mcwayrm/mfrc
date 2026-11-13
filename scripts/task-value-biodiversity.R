# First-pass biodiversity (nongame habitat quality) valuation for Minnesota forests.
#
# Quantity: breeding-bird community diversity from the USGS North American
# Breeding Bird Survey (BBS). Route-year counts are cleaned (RunType = 1,
# RPID = 101, roadside 50-stop routes, diurnal landbirds only, no
# unidentified/hybrid codes), summarized as effective species (Hill q = 1,
# exp of Shannon H), and interpolated with a spatiotemporal GAM
# (log effN ~ te(x, y, year) + observer random effect + first-year indicator)
# onto the 300 m forest grid. The habitat-quality change for each pixel is the
# percent change of predicted effective species relative to the pixel's
# 1966-1970 baseline mean.
#
# Price: ForCAST stated-preference (choice experiment) WTP of $1.15 per
# household per month for a 10% increase in nongame species habitat quality
# (https://mnforcast.org/ecosystem-valuation, Table 1), scaled linearly,
# applied to Census household counts, and allocated across forest pixels so
# pixel values sum to the statewide annual WTP.

required_packages <- c("terra", "sf", "tmap", "jsonlite", "data.table", "mgcv")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
}
if (utils::packageVersion("tmap") < "4.0") stop("Biodiversity maps require tmap >= 4.0.")
source_paths <- vapply(sys.frames(), function(frame) {
  if (is.character(frame$ofile) && length(frame$ofile) == 1L) frame$ofile else NA_character_
}, character(1))
source_paths <- source_paths[!is.na(source_paths)]
if (length(source_paths)) {
  script_path <- tail(source_paths, 1L)
} else {
  script_argument <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_argument) != 1L) stop("Run with Rscript or source() using the script path.")
  script_path <- sub("^--file=", "", script_argument)
}
root_dir <- dirname(dirname(normalizePath(script_path, mustWork = TRUE)))

effective_species <- function(counts) {
  if (!length(counts) || any(!is.finite(counts) | counts < 0) || sum(counts) <= 0) {
    stop("Effective species requires finite nonnegative counts with a positive total.")
  }
  shares <- counts[counts > 0] / sum(counts)
  exp(-sum(shares * log(shares)))
}

download_bbs_raw <- function(raw_dir) {
  release <- "https://www.sciencebase.gov/catalog/file/get/6a0b0b0ab66b0188da36aedd?f="
  files <- c(
    "Routes.csv" = "__disk__85%2Ff3%2F30%2F85f33018ad86e77e43031f08e198e53d73532eb8",
    "Weather.csv" = "__disk__99%2F9b%2Faa%2F999baa5fdb09c96980a32e4d599892ae563b28c8",
    "SpeciesList.csv" = "__disk__67%2F26%2F5b%2F67265bbc0d9312830188fb7e683ed344f19f5452",
    "States.zip" = "__disk__fb%2Ff3%2Fb6%2Ffbf3b68cd779b281a77c7b10b8289c8bf6b77497"
  )
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(raw_dir)) stop("Could not create BBS raw-data directory: ", raw_dir)
  options(timeout = max(3600, getOption("timeout")))
  for (name in names(files)) {
    path <- file.path(raw_dir, name)
    if (!file.exists(path)) {
      message("Downloading BBS 2026 release file: ", name)
      utils::download.file(paste0(release, files[[name]]), path, mode = "wb", quiet = TRUE)
    }
    if (!file.exists(path) || file.info(path)$size == 0) stop("Missing or empty BBS file: ", path)
  }
  states_dir <- file.path(raw_dir, "States")
  if (!dir.exists(states_dir)) {
    utils::unzip(file.path(raw_dir, "States.zip"), exdir = states_dir)
  }
  state_files <- list.files(states_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  if (!length(state_files)) stop("No extracted BBS state count files under: ", states_dir)
  state_files
}

read_route_diversity <- function(raw_dir, boundary, buffer_m = 100000) {
  state_files <- download_bbs_raw(raw_dir)
  # Minnesota plus adjacent states/provinces support predictions near borders.
  wanted <- c("Minneso", "Wiscons", "Iowa", "NDakota", "SDakota", "Manitob", "Ontario")
  selected_files <- state_files[sub("\\.csv$", "", basename(state_files)) %in% wanted]
  if (length(selected_files) != length(wanted)) {
    stop("Missing BBS state count files: ",
         paste(setdiff(wanted, sub("\\.csv$", "", basename(selected_files))), collapse = ", "))
  }
  routes <- data.table::fread(file.path(raw_dir, "Routes.csv"), showProgress = FALSE)
  route_columns <- c("CountryNum", "StateNum", "Route", "Latitude", "Longitude", "RouteTypeID")
  if (!all(route_columns %in% names(routes))) stop("Routes.csv is missing required columns.")
  routes <- routes[RouteTypeID == 1 & is.finite(Latitude) & is.finite(Longitude) &
                     Latitude >= -90 & Latitude <= 90 & Longitude >= -180 & Longitude <= 180]
  if (anyDuplicated(routes[, c("CountryNum", "StateNum", "Route"), with = FALSE])) {
    stop("Routes.csv contains duplicated route identifiers.")
  }
  points <- sf::st_transform(sf::st_as_sf(
    as.data.frame(routes), coords = c("Longitude", "Latitude"), crs = 4326
  ), 26915)
  zone <- sf::st_buffer(sf::st_union(boundary), buffer_m)
  keep <- lengths(sf::st_intersects(points, zone)) > 0L
  coordinates <- sf::st_coordinates(points)
  routes[, `:=`(x = coordinates[, 1], y = coordinates[, 2])]
  routes <- routes[keep]
  if (nrow(routes) < 50L) stop("Fewer than 50 roadside BBS routes near Minnesota.")
  weather <- data.table::fread(
    file.path(raw_dir, "Weather.csv"),
    select = c("RouteDataID", "CountryNum", "StateNum", "Route", "RPID", "Year",
               "ObsN", "RunType"),
    showProgress = FALSE
  )
  # Standard-protocol runs that met BBS data-quality criteria.
  weather <- weather[RPID == 101 & RunType == 1]
  if (anyNA(weather) || anyDuplicated(weather$RouteDataID)) {
    stop("Filtered BBS Weather records must be complete with unique run identifiers.")
  }
  weather <- merge(weather, routes[, .(CountryNum, StateNum, Route, x, y)],
                   by = c("CountryNum", "StateNum", "Route"), sort = FALSE)
  # One run per route-year: keep the lowest run identifier deterministically.
  data.table::setorder(weather, CountryNum, StateNum, Route, Year, RouteDataID)
  weather <- weather[, .SD[1L], by = .(CountryNum, StateNum, Route, Year)]
  species <- data.table::fread(file.path(raw_dir, "SpeciesList.csv"), showProgress = FALSE,
                               encoding = "Latin-1")
  if (!all(c("AOU", "English_Common_Name", "Order") %in% names(species))) {
    stop("SpeciesList.csv is missing required columns.")
  }
  # Diurnal landbird orders the BBS samples adequately; excludes waterbirds,
  # shorebirds, raptors, and nocturnal species per the valuation design.
  landbird_orders <- c("Passeriformes", "Piciformes", "Columbiformes", "Cuculiformes",
                       "Apodiformes", "Coraciiformes", "Galliformes", "Trogoniformes",
                       "Psittaciformes")
  kept_species <- species[Order %in% landbird_orders &
                            !grepl("unid\\.|hybrid", English_Common_Name, ignore.case = TRUE)]
  if (nrow(kept_species) < 200L) stop("Unexpectedly small landbird species list.")
  counts <- data.table::rbindlist(lapply(selected_files, function(path) {
    data.table::fread(path, select = c("RouteDataID", "RPID", "AOU", "SpeciesTotal"),
                      showProgress = FALSE)
  }))
  counts <- counts[RPID == 101 & RouteDataID %in% weather$RouteDataID &
                     AOU %in% kept_species$AOU]
  if (!nrow(counts) || any(!is.finite(counts$SpeciesTotal) | counts$SpeciesTotal <= 0)) {
    stop("BBS counts must be positive route totals for retained runs and species.")
  }
  totals <- counts[, .(count = sum(SpeciesTotal)), by = .(RouteDataID, AOU)]
  diversity <- totals[, .(
    effective_species = effective_species(count),
    richness = data.table::uniqueN(AOU),
    individuals = sum(count)
  ), by = RouteDataID]
  samples <- merge(weather, diversity, by = "RouteDataID", sort = FALSE)
  samples[, route_id := paste(CountryNum, StateNum, Route, sep = "-")]
  samples[, first_year_observer := as.integer(Year == min(Year)), by = .(route_id, ObsN)]
  samples[, `:=`(observer = factor(ObsN), year = as.numeric(Year),
                 log_effective_species = log(effective_species))]
  if (any(!is.finite(samples$log_effective_species))) stop("Nonfinite route diversity values.")
  if (nrow(samples) < 1500L) {
    stop("Fewer than 1500 quality route-years near Minnesota; refusing to fit the surface.")
  }
  message("Prepared ", nrow(samples), " BBS route-years from ",
          data.table::uniqueN(samples$route_id), " routes (",
          min(samples$Year), "-", max(samples$Year), ").")
  samples
}

fit_diversity_model <- function(samples) {
  fit <- mgcv::bam(
    log_effective_species ~ first_year_observer +
      te(x, y, year, d = c(2, 1), bs = c("tp", "cr"), k = c(30, 10)) +
      s(observer, bs = "re"),
    data = samples, method = "fREML", discrete = TRUE
  )
  if (any(!is.finite(stats::coef(fit)))) stop("Diversity GAM produced nonfinite coefficients.")
  fitted_correlation <- stats::cor(stats::fitted(fit), samples$log_effective_species)
  if (!is.finite(fitted_correlation) || fitted_correlation <= 0.3) {
    stop("Diversity GAM fits the route-year data too poorly (correlation <= 0.3).")
  }
  message("Fitted spatiotemporal diversity GAM; deviance explained: ",
          round(summary(fit)$dev.expl * 100, 1), "%.")
  fit
}

# Exact population-level prediction (intercept + tensor smooth only). The
# observer random effect and first-year indicator contribute zero at the
# population level, so building their prediction columns is skipped; results
# match predict(fit, ..., exclude = "s(observer)") with first_year_observer = 0.
predict_log_diversity <- function(fit, xy, year, with_se = FALSE, chunk_size = 50000L) {
  labels <- vapply(fit$smooth, function(smooth) smooth$label, character(1))
  tensor_index <- which(startsWith(labels, "te("))
  if (length(tensor_index) != 1L) stop("Expected exactly one tensor smooth in the GAM.")
  smooth <- fit$smooth[[tensor_index]]
  columns <- c(1L, smooth$first.para:smooth$last.para)
  beta <- stats::coef(fit)[columns]
  covariance <- fit$Vp[columns, columns]
  n <- nrow(xy)
  estimates <- numeric(n)
  standard_errors <- if (with_se) numeric(n) else NULL
  for (start in seq.int(1L, n, by = chunk_size)) {
    rows <- start:min(start + chunk_size - 1L, n)
    design <- cbind(1, mgcv::PredictMat(smooth, data.frame(
      x = xy[rows, 1], y = xy[rows, 2], year = year
    )))
    estimates[rows] <- drop(design %*% beta)
    if (with_se) standard_errors[rows] <- sqrt(rowSums((design %*% covariance) * design))
  }
  if (any(!is.finite(estimates)) || (with_se && any(!is.finite(standard_errors) |
                                                    standard_errors < 0))) {
    stop("Nonfinite diversity predictions for year ", year, ".")
  }
  list(fit = estimates, se = standard_errors)
}

run_biodiversity <- function(root_dir) {
  years <- c(2005L, 2010L, 2015L, 2020L, 2025L)
  baseline_years <- 1966:1970
  wtp_per_household_per_month_per_10pct <- 1.15
  # Annual USD per household for each one percent habitat-quality change.
  annual_price_per_pct <- wtp_per_household_per_month_per_10pct * 12 / 10
  households <- c(2020144, 2091548, 2147262, 2253990, 2391870)
  household_sources <- c(
    "ACS 2005 1-year, table B11001 (summary file sequence 46)",
    "ACS 2010 1-year, table B11001",
    "ACS 2015 1-year, table B11001",
    "2020 Decennial Census DHC, table P16 (households)",
    "ACS 2025 1-year, table B11001"
  )
  input_dir <- file.path(root_dir, "data", "forest-layer")
  raw_dir <- file.path(root_dir, "data", "BBS", "raw")
  output_dir <- file.path(root_dir, "data", "GEP")
  quantity_dir <- file.path(root_dir, "data", "forest-biodiversity")
  map_dir <- file.path(root_dir, "outputs", "maps")
  figure_dir <- file.path(root_dir, "outputs", "figures")
  table_dir <- file.path(root_dir, "outputs", "tables")
  template_paths <- file.path(input_dir, sprintf("forest-layer-%d.tif", years))
  metadata_paths <- sub("\\.tif$", ".json", template_paths)
  quantity_paths <- file.path(quantity_dir, sprintf("forest-biodiversity-%d.tif", years))
  value_paths <- file.path(output_dir, sprintf("gep-biodiversity-%d.tif", years))
  sample_path <- file.path(quantity_dir, "bbs-route-diversity.csv")
  price_path <- file.path(output_dir, "biodiversity-wtp-prices.csv")
  summary_path <- file.path(table_dir, "biodiversity-time-series.csv")
  coverage_path <- file.path(table_dir, "biodiversity-route-coverage.csv")
  boundary_path <- file.path(root_dir, "data", "david-wilson", "GIS", "bdry_counpy2.shp")
  inputs <- c(template_paths, metadata_paths, boundary_path)
  if (any(!file.exists(inputs))) {
    stop("Missing inputs: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
  }
  outputs <- c(quantity_paths, sub("\\.tif$", ".json", quantity_paths), value_paths,
               sub("\\.tif$", ".json", value_paths), sample_path, price_path,
               summary_path, coverage_path)
  if (any(file.exists(outputs))) {
    stop("Refusing to overwrite biodiversity outputs: ",
         paste(outputs[file.exists(outputs)], collapse = ", "))
  }
  metadata <- lapply(metadata_paths, jsonlite::fromJSON)
  templates <- lapply(seq_along(years), function(index) {
    raster <- terra::rast(template_paths[index])
    info <- metadata[[index]]
    if (!terra::same.crs(raster, "EPSG:26915") || any(terra::res(raster) != 300) ||
        info$target_year != years[index] ||
        info$mask_year != if (years[index] == 2025L) 2022L else years[index]) {
      stop("Unexpected forest template or metadata: ", template_paths[index])
    }
    raster[[1]]
  })
  for (index in seq_along(years)[-1]) {
    terra::compareGeom(templates[[1]], templates[[index]], stopOnError = TRUE)
  }
  boundary <- sf::st_read(boundary_path, quiet = TRUE)
  if (is.na(sf::st_crs(boundary)) || any(!sf::st_is_valid(boundary))) {
    stop("Invalid Minnesota boundary.")
  }
  boundary <- sf::st_sf(geometry = sf::st_union(sf::st_transform(boundary, 26915)))
  samples <- read_route_diversity(raw_dir, boundary)
  if (!all(baseline_years %in% samples$Year) || max(samples$Year) < max(years)) {
    stop("BBS route-years do not span the baseline and target years.")
  }
  fit <- fit_diversity_model(samples)
  for (directory in c(output_dir, quantity_dir, map_dir, figure_dir, table_dir)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(directory)) stop("Could not create directory: ", directory)
  }
  warning("Breeding-season landbird diversity from roadside BBS routes is a partial proxy for ",
          "nongame habitat quality. The $1.15/household/month WTP from an 11-county ForCAST ",
          "survey is scaled linearly and transferred statewide without ownership attribution; ",
          "validate before GEP inference.")
  data.table::fwrite(samples[, .(route_id, year = Year, x, y, observer = ObsN,
                                 first_year_observer, effective_species, richness,
                                 individuals)], sample_path)
  coverage <- samples[, .(route_years = .N, routes = data.table::uniqueN(route_id),
                          mean_effective_species = mean(effective_species)),
                      keyby = .(year = Year)]
  data.table::fwrite(coverage, coverage_path)
  price_table <- data.frame(
    target_year = years,
    wtp_usd_per_household_per_month_per_10pct = wtp_per_household_per_month_per_10pct,
    annual_usd_per_household_per_pct = annual_price_per_pct,
    households = households,
    household_source = household_sources,
    price_source = "https://mnforcast.org/ecosystem-valuation (Table 1, nongame +10%)",
    assumption = paste("Linear WTP scaling per percent; symmetric for losses;",
                       "11-county survey transferred to all Minnesota households;",
                       "nominal survey-period USD (price year not stated by source)")
  )
  utils::write.csv(price_table, price_path, row.names = FALSE)
  mask_cells <- lapply(templates, function(template) {
    which(!is.na(terra::values(template, mat = FALSE)))
  })
  for (index in seq_along(years)) {
    if (length(mask_cells[[index]]) != metadata[[index]]$forest_pixels) {
      stop("Template forest count mismatch: ", years[index])
    }
  }
  union_cells <- sort(unique(unlist(mask_cells)))
  xy <- terra::xyFromCell(templates[[1]], union_cells)
  message("Predicting baseline log diversity for ", length(union_cells), " forest pixels.")
  baseline_log <- rowMeans(vapply(baseline_years, function(year) {
    predict_log_diversity(fit, xy, year)$fit
  }, numeric(length(union_cells))))
  cell_area_ha <- prod(terra::res(templates[[1]])) / 10000
  summaries <- vector("list", length(years))
  values_by_year <- vector("list", length(years))
  for (index in seq_along(years)) {
    year <- years[index]
    message("Predicting ", year, " diversity surface.")
    prediction <- predict_log_diversity(fit, xy, year, with_se = TRUE)
    rows <- match(mask_cells[[index]], union_cells)
    effective <- exp(prediction$fit[rows])
    pct_change <- 100 * (exp(prediction$fit[rows] - baseline_log[rows]) - 1)
    standard_error <- prediction$se[rows]
    cells <- mask_cells[[index]]
    forest_area_ha <- length(cells) * cell_area_ha
    # Pixel values decompose statewide WTP: summing value density times pixel
    # area reproduces households * annual price * (area-weighted mean pct / 1%).
    value_per_ha <- households[index] * annual_price_per_pct * pct_change / forest_area_ha
    quantity <- terra::rast(templates[[index]], nlyrs = 3)
    terra::values(quantity) <- NA_real_
    quantity[[1]][cells] <- effective
    quantity[[2]][cells] <- pct_change
    quantity[[3]][cells] <- standard_error
    names(quantity) <- c("effective_species_hill_q1",
                         "pct_change_from_1966_1970_baseline",
                         "log_scale_prediction_se")
    value <- terra::rast(templates[[index]])
    terra::values(value) <- NA_real_
    value[cells] <- value_per_ha
    names(value) <- "annual_nongame_habitat_value_USD_per_ha_per_year"
    for (item in list(list(raster = quantity, path = quantity_paths[index]),
                      list(raster = value, path = value_paths[index]))) {
      terra::writeRaster(item$raster, item$path, overwrite = FALSE, datatype = "FLT8S",
                         NAflag = -9999, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES"))
      saved <- terra::rast(item$path)
      terra::compareGeom(item$raster, saved, stopOnError = TRUE)
      for (band in seq_len(terra::nlyr(saved))) {
        actual <- terra::values(saved[[band]], mat = FALSE)
        expected <- terra::values(item$raster[[band]], mat = FALSE)
        if (!identical(is.na(actual), is.na(expected)) ||
            any(abs(actual[!is.na(actual)] - expected[!is.na(expected)]) > 1e-8)) {
          stop("Written quantity/value raster failed verification: ", item$path)
        }
      }
    }
    info <- list(
      target_year = year, baseline_years = range(baseline_years),
      quantity = "Predicted effective landbird species (Hill q = 1) and percent change",
      quantity_units = "effective species; percent change from 1966-1970 pixel baseline",
      value_units = "nominal survey-period USD per hectare per year",
      formula = paste("households * ($1.15 * 12 / 10%) * pixel_pct_change /",
                      "forest_area_ha; pixels sum to statewide annual WTP"),
      wtp_usd_per_household_per_month_per_10pct = wtp_per_household_per_month_per_10pct,
      households = households[index], household_source = household_sources[index],
      price_source = "https://mnforcast.org/ecosystem-valuation (Table 1)",
      bbs_release = paste("2026 Release - North American Breeding Bird Survey Dataset",
                          "(1966-2025); https://www.sciencebase.gov/catalog/item/",
                          "6a0b0b0ab66b0188da36aedd"),
      data_filters = paste("RunType = 1; RPID = 101; RouteTypeID = 1; diurnal landbird",
                           "orders only; unidentified and hybrid codes dropped;",
                           "lowest RouteDataID per route-year"),
      model = paste("bam(log effN ~ first-year observer + te(x, y, year, d = c(2, 1),",
                    "bs = c('tp', 'cr'), k = c(30, 10)) + s(observer, 're'), fREML)"),
      model_route_years = nrow(samples),
      model_routes = data.table::uniqueN(samples$route_id),
      deviance_explained = summary(fit)$dev.expl,
      prediction = "Observer effects excluded; first-year indicator fixed at zero",
      mask_template = template_paths[index], esa_mask_year = metadata[[index]]$mask_year,
      forest_pixels = length(cells), forest_area_ha = forest_area_ha,
      mean_pct_change = mean(pct_change),
      total_annual_value_USD = sum(value_per_ha) * cell_area_ha,
      interpretation = paste(
        "Breeding-season roadside landbird diversity proxies nongame habitat quality;",
        "it under-represents interior habitat, rare species, and other taxa.",
        "WTP is scaled linearly per percent and symmetrically for losses, although the",
        "source study reports non-linear and loss-averse preferences.",
        "An 11-county survey WTP is transferred to all Minnesota households.",
        "Statewide WTP is allocated only to forest pixels in proportion to their",
        "predicted change; negative pixel changes yield negative values.",
        "Predictions are smoothed; pixels far from routes are weakly informed.",
        "The 2020 BBS season was largely cancelled, so 2020 is interpolated.",
        "SE layer is the log-scale GAM prediction error for the target year only,",
        "ignoring baseline uncertainty; treat downstream use as a generated regressor."
      )
    )
    for (path in c(sub("\\.tif$", ".json", quantity_paths[index]),
                   sub("\\.tif$", ".json", value_paths[index]))) {
      jsonlite::write_json(info, path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
    }
    summaries[[index]] <- data.frame(
      target_year = year, esa_mask_year = metadata[[index]]$mask_year,
      forest_cells = length(cells), forest_area_ha = forest_area_ha,
      mean_effective_species = mean(effective),
      mean_pct_change_from_baseline = mean(pct_change),
      minimum_pct_change = min(pct_change), maximum_pct_change = max(pct_change),
      households = households[index],
      annual_usd_per_household_per_pct = annual_price_per_pct,
      total_annual_value_USD = sum(value_per_ha) * cell_area_ha,
      total_annual_value_fixed_2005_households_USD =
        sum(households[1] * annual_price_per_pct * pct_change / forest_area_ha) * cell_area_ha,
      mean_annual_value_USD_per_ha = mean(value_per_ha),
      scope = "Statewide BBS landbird-diversity interpolation; forest pixels only"
    )
    values_by_year[[index]] <- value
    message("Saved ", year, " biodiversity quantity and value rasters; mean change ",
            round(mean(pct_change), 2), "% from the 1966-1970 baseline.")
  }
  summary <- do.call(rbind, summaries)
  utils::write.csv(summary, summary_path, row.names = FALSE)
  annual_series <- statewide_annual_series(fit, xy, seq(min(samples$Year), max(samples$Year)))
  visualize_biodiversity(values_by_year, summary, samples, annual_series, boundary,
                         baseline_years, map_dir, figure_dir)
  invisible(summary)
}

statewide_annual_series <- function(fit, xy, years, sample_size = 4000L) {
  set.seed(42)
  rows <- sort(sample.int(nrow(xy), min(sample_size, nrow(xy))))
  subgrid <- xy[rows, , drop = FALSE]
  means <- vapply(years, function(year) {
    mean(exp(predict_log_diversity(fit, subgrid, year)$fit))
  }, numeric(1))
  data.frame(year = years, mean_effective_species = means)
}

visualize_biodiversity <- function(values, summary, samples, annual_series, boundary,
                                   baseline_years, map_dir, figure_dir) {
  old_mode <- tmap::tmap_mode("plot")
  on.exit(tmap::tmap_mode(old_mode), add = TRUE)
  first_year <- summary$target_year[1]
  last_year <- tail(summary$target_year, 1)
  last <- values[[length(values)]]
  change <- last - values[[1]]
  maps <- list()
  maps[[paste0("biodiversity-value-", last_year)]] <- list(
    raster = last, title = paste0("A. Nongame habitat-quality value: ", last_year),
    legend = "Annual value\nUSD / ha / year",
    note = paste0("BBS landbird diversity vs 1966-1970 baseline; ESA mask: ",
                  tail(summary$esa_mask_year, 1),
                  ".\nWTP: $1.15/household/month per 10% change (ForCAST Table 1).")
  )
  maps[[paste0("biodiversity-value-change-", first_year, "-", last_year)]] <- list(
    raster = change,
    title = paste0("B. Annual value difference: ", first_year, "-", last_year),
    legend = "Value difference\nUSD / ha / year",
    note = paste0("Shared forest pixels only; absolute difference, not percentage.",
                  "\nIncludes household growth between ", first_year, " and ", last_year, ".")
  )
  for (name in names(maps)) {
    item <- maps[[name]]
    vals <- terra::values(item$raster, mat = FALSE)
    bound <- max(abs(vals), na.rm = TRUE)
    if (!is.finite(bound) || bound <= 0) stop("No nonzero values to visualize: ", name)
    map <- tmap::tm_shape(boundary, crs = 26915, is.main = TRUE) +
      tmap::tm_polygons(fill = "#eeeeee", col = NA) +
      tmap::tm_shape(item$raster) +
      tmap::tm_raster(
        col.scale = tmap::tm_scale_continuous(
          limits = c(-bound, bound), midpoint = 0,
          values = c("#b2182b", "#f7f7f7", "#2166ac"), value.na = NA
        ),
        col.legend = tmap::tm_legend(title = item$legend, na.show = FALSE)
      ) +
      tmap::tm_shape(boundary) + tmap::tm_borders(col = "#555555", lwd = 0.7) +
      tmap::tm_title(item$title, size = 0.85, position = tmap::tm_pos_out("center", "top")) +
      tmap::tm_credits(item$note, size = 0.65, position = tmap::tm_pos_out("left", "bottom")) +
      tmap::tm_layout(frame = FALSE, meta.margins = c(0.12, 0, 0.07, 0.25),
                      legend.position = tmap::tm_pos_out("right", "center"),
                      legend.title.size = 0.8, legend.text.size = 0.7)
    for (extension in c("png", "pdf")) {
      path <- file.path(map_dir, paste0(name, ".", extension))
      tmap::tmap_save(map, path, width = 7, height = 7, units = "in", dpi = 300)
      if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) {
        stop("Map export failed: ", path)
      }
    }
  }
  for (extension in c("png", "pdf")) {
    path <- file.path(figure_dir, paste0("biodiversity-time-series.", extension))
    if (extension == "png") {
      grDevices::png(path, width = 8, height = 8, units = "in", res = 300, type = "cairo")
    } else {
      grDevices::pdf(path, width = 8, height = 8)
    }
    tryCatch({
      graphics::par(mfrow = c(2, 1), mar = c(4.5, 5, 2.5, 1), oma = c(4, 0, 2, 0), las = 1)
      graphics::plot(samples$Year, samples$effective_species, pch = 16, cex = 0.3,
                     col = grDevices::adjustcolor("#777777", alpha.f = 0.25),
                     xlab = "", ylab = "Effective species (Hill q = 1)",
                     main = "BBS route-year landbird diversity near Minnesota")
      graphics::rect(min(baseline_years), graphics::par("usr")[3], max(baseline_years),
                     graphics::par("usr")[4],
                     col = grDevices::adjustcolor("#fee391", alpha.f = 0.4), border = NA)
      graphics::lines(annual_series$year, annual_series$mean_effective_species,
                      col = "#238443", lwd = 2.5)
      graphics::legend("topleft",
                       c("Route-year observations", "Model statewide forest mean",
                         "Baseline window"),
                       pch = c(16, NA, 15), lty = c(NA, 1, NA), lwd = c(NA, 2.5, NA),
                       col = c("#777777", "#238443",
                               grDevices::adjustcolor("#fee391", alpha.f = 0.8)),
                       bty = "n", cex = 0.8)
      annual <- summary$total_annual_value_USD / 1e6
      fixed <- summary$total_annual_value_fixed_2005_households_USD / 1e6
      graphics::plot(summary$target_year, annual, type = "o", pch = 16, col = "#2166ac",
                     lwd = 2, xaxt = "n", ylim = range(c(annual, fixed, 0)),
                     xlab = "Target year", ylab = "Million USD / year",
                     main = "Statewide annual WTP for nongame habitat-quality change")
      graphics::axis(1, at = summary$target_year)
      graphics::abline(h = 0, col = "#aaaaaa", lty = 3)
      graphics::lines(summary$target_year, fixed, type = "o", pch = 1, lty = 2,
                      col = "#666666")
      graphics::legend("topleft", c("Target-year households", "Fixed 2005 households"),
                       lty = c(1, 2), col = c("#2166ac", "#666666"), bty = "n", cex = 0.8)
      graphics::mtext("C. Exploratory statewide biodiversity valuation", side = 3,
                      outer = TRUE, font = 2)
      graphics::mtext(
        paste("Roadside breeding-season landbird diversity vs the 1966-1970 baseline;",
              "WTP $1.15/household/month per 10% change.\nSmoothed interpolation of BBS",
              "routes, not a census of biodiversity or family-forest-only GEP."),
        side = 1, outer = TRUE, line = 0.5, cex = 0.75
      )
    }, finally = grDevices::dev.off())
    if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) {
      stop("Time-series export failed: ", path)
    }
  }
}

run_biodiversity(root_dir)
