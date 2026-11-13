# Historical timber yield and stumpage-value figures for the ecosystem-service report.

required_packages <- c("terra", "sf", "tmap", "data.table", "rFIA", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
}

source_paths <- vapply(sys.frames(), function(frame) {
  if (is.character(frame$ofile) && length(frame$ofile) == 1L) frame$ofile else NA_character_
}, character(1))
source_paths <- source_paths[!is.na(source_paths)]
if (length(source_paths)) {
  script_path <- tail(source_paths, 1L)
} else {
  script_argument <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_argument) != 1L) stop("Run this script with Rscript.")
  script_path <- sub("^--file=", "", script_argument)
}
script_path <- normalizePath(script_path, mustWork = TRUE)
root_dir <- dirname(dirname(script_path))
source(file.path(dirname(script_path), "task-value-timber.R"))

historical_prices <- function(year) {
  if (year == 2021L) {
    prices <- stumpage_prices_2021()[, .(product, species_group, price_per_unit)]
  } else if (year == 2005L) {
    pulpwood <- data.table::data.table(
      group = c("Aspen", "Balm", "Birch", "Ash", "Elm", "Oak", "Basswood",
                "Maple", "Mixed/Other Hdwds", "Balsam fir", "Spruce Species",
                "Tamarack", "Cedar", "Jack Pine", "Red/White Pine"),
      volume = c(732923.7, 17927.2, 110195.2, 11446.8, 33.5, 19953.9,
                 14809.0, 54741.4, 27395.4, 76520.4, 119210.7, 59419.5,
                 2062.3, 88264.7, 29107.4),
      price = c(59.70, 45.25, 20.57, 5.43, 0.60, 16.28, 10.64, 13.30,
                22.98, 33.54, 43.39, 9.84, 5.50, 30.66, 29.06)
    )
    pulp_bolts <- data.table::data.table(
      group = pulpwood$group,
      volume = c(865293.3, 20978.1, 122630.2, 15144.6, 132.2, 30171.9,
                 19410.3, 58057.4, 33992.7, 87167.2, 135529.3, 62232.8,
                 3373.5, 165235.6, 72988.1),
      price = c(65.14, 47.09, 24.99, 21.76, 32.47, 42.24, 23.89, 16.30,
                23.69, 41.38, 48.03, 10.07, 7.47, 50.81, 55.17)
    )
    cord <- rbind(pulpwood, pulp_bolts)[, .(
      price_per_unit = weighted.mean(price, volume),
      total_volume = sum(volume)
    ), by = group]
    data.table::setnames(cord, "group", "species_group")
    spruce_price <- cord[species_group == "Spruce Species", price_per_unit]
    pine_price <- cord[species_group == "Red/White Pine", price_per_unit]
    hardwood_groups <- c("Aspen", "Balm", "Birch", "Ash", "Elm", "Oak",
      "Basswood", "Maple", "Mixed/Other Hdwds")
    softwood_groups <- c("Balsam fir", "Spruce Species", "Tamarack", "Cedar",
      "Jack Pine", "Red/White Pine")
    composite_prices <- data.table::data.table(
      species_group = c("W Spruce", "B Spruce", "Red Pine", "White Pine",
        "Hdwd Cd Composite", "Sftwd Cd Composite"),
      price_per_unit = c(
        spruce_price, spruce_price, pine_price, pine_price,
        weighted.mean(cord[species_group %in% hardwood_groups, price_per_unit],
          cord[species_group %in% hardwood_groups, total_volume]),
        weighted.mean(cord[species_group %in% softwood_groups, price_per_unit],
          cord[species_group %in% softwood_groups, total_volume])
      )
    )
    cord <- data.table::rbindlist(list(
      cord[, .(species_group, price_per_unit)],
      composite_prices
    ), use.names = TRUE)
    cord[, product := "cord"]
    saw <- data.table::data.table(
      species_group = c("Ash", "Aspen", "Balm", "Birch", "Elm", "Oak",
        "Basswood", "Maple", "Mixed/Other Hdwds", "Balsam fir",
        "Spruce Species", "Tamarack", "Cedar", "Jack Pine", "Red/White Pine",
        "MBF Composite"),
      price_per_unit = c(144.62, 190.44, 115.89, 128.30, 86.52, 185.90,
        133.10, 131.53, 66.62, 244.43, 204.73, 37.19, 26.38, 184.79,
        181.21, 180.73),
      product = "sawtimber"
    )
    prices <- rbind(cord, saw)
  } else {
    cord_by_year <- list(
      `2010` = c(Ash=8.60, Aspen=25.16, Balm=21.22, Birch=10.37, Oak=19.51,
        Basswood=10.45, Maple=9.36, `Mixed/Other Hdwds`=12.29,
        `Balsam fir`=17.38, `W Spruce`=22.43, `B Spruce`=22.43,
        Tamarack=5.03, Cedar=6.19, `Jack Pine`=27.42, `Red Pine`=29.42,
        `White Pine`=29.42, `Hdwd Cd Composite`=20.39,
        `Sftwd Cd Composite`=25.87),
      `2015` = c(Ash=10.88, Aspen=36.64, Balm=27.68, Birch=15.60, Oak=26.06,
        Basswood=14.19, Maple=13.50, `Mixed/Other Hdwds`=14.08,
        `Balsam fir`=15.84, `W Spruce`=22.12, `B Spruce`=23.79,
        Tamarack=8.09, Cedar=6.41, `Jack Pine`=29.85, `Red Pine`=38.75,
        `White Pine`=38.66, `Hdwd Cd Composite`=28.47,
        `Sftwd Cd Composite`=29.56),
      `2020` = c(Ash=8.55, Aspen=30.14, Balm=25.22, Birch=12.64, Oak=27.81,
        Basswood=11.96, Maple=11.82, `Mixed/Other Hdwds`=11.73,
        `Balsam fir`=11.81, `W Spruce`=20.66, `B Spruce`=20.84,
        Tamarack=6.72, Cedar=6.93, `Jack Pine`=24.85, `Red Pine`=31.62,
        `White Pine`=27.34, `Hdwd Cd Composite`=25.52,
        `Sftwd Cd Composite`=21.61)
    )
    saw_by_year <- list(
      `2010` = c(Ash=56.27, Birch=38.92, Elm=45.08, Oak=243.09,
        Basswood=63.47, Maple=219.83, `Mixed/Other Hdwds`=29.80,
        `W Spruce`=102.15, `Jack Pine`=105.00, `Red Pine`=123.36,
        `White Pine`=123.36, `MBF Composite`=115.16),
      `2015` = c(Ash=54.17, Birch=42.84, Elm=42.50, Oak=265.54,
        Basswood=59.24, Maple=126.70, `Mixed/Other Hdwds`=47.87,
        `W Spruce`=74.68, `Jack Pine`=89.56, `Red Pine`=177.18,
        `White Pine`=88.92, `MBF Composite`=160.13),
      `2020` = c(Ash=61.14, Birch=51.69, Elm=72.91, Oak=161.13,
        Basswood=75.34, Maple=94.29, `Mixed/Other Hdwds`=67.78,
        `W Spruce`=82.53, `Jack Pine`=105.86, `Red Pine`=128.10,
        `White Pine`=109.90, `MBF Composite`=125.95)
    )
    cord <- data.table::data.table(
      species_group = names(cord_by_year[[as.character(year)]]),
      price_per_unit = unname(cord_by_year[[as.character(year)]]),
      product = "cord"
    )
    saw <- data.table::data.table(
      species_group = names(saw_by_year[[as.character(year)]]),
      price_per_unit = unname(saw_by_year[[as.character(year)]]),
      product = "sawtimber"
    )
    prices <- data.table::rbindlist(list(cord, saw), fill = TRUE)
  }
  source_urls <- c(
    `2005` = "https://files.dnr.state.mn.us/forestry/timber_sales/stumpage/aveprice05.pdf",
    `2010` = "https://files.dnr.state.mn.us/forestry/timber_sales/stumpage/stumpage-review-report-2015.pdf",
    `2015` = "https://files.dnr.state.mn.us/forestry/timber_sales/stumpage/stumpage-review-report-2015.pdf",
    `2020` = "https://files.dnr.state.mn.us/forestry/timber_sales/stumpage/stumpage-review-report-2020.pdf",
    `2021` = "https://files.dnr.state.mn.us/forestry/timber_sales/stumpage/stumpage-review-report-2021.pdf"
  )
  table_sources <- c(
    `2005` = "All-agencies tables 1-2; cord prices are volume-weighted pulpwood and pulp-and-bolts averages",
    `2010` = "2015 report historical tables 2-3, 2010 column",
    `2015` = "Tables 2-3",
    `2020` = "Tables 2-3",
    `2021` = "Tables 2-3"
  )
  prices[, `:=`(
    price_year = as.integer(year), currency = "nominal USD",
    source_url = unname(source_urls[as.character(year)]),
    source_table = unname(table_sources[as.character(year)])
  )]
  prices[]
}

target_years <- c(2005L, 2010L, 2015L, 2020L, 2025L)
price_years <- c(2005L, 2010L, 2015L, 2020L, 2021L)
fia_dir <- file.path(root_dir, "data", "FIA")
forest_dir <- file.path(root_dir, "data", "forest-layer")
gep_dir <- file.path(root_dir, "data", "GEP")
map_dir <- file.path(root_dir, "outputs", "maps")
figure_dir <- file.path(root_dir, "outputs", "figures")
table_dir <- file.path(root_dir, "outputs", "tables")
boundary_path <- file.path(root_dir, "data", "david-wilson", "GIS", "bdry_counpy2.shp")
for (directory in c(gep_dir, map_dir, figure_dir, table_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(directory)) stop("Could not create output directory: ", directory)
}
boundary <- sf::st_read(boundary_path, quiet = TRUE)
if (is.na(sf::st_crs(boundary)) || any(!sf::st_is_valid(boundary))) {
  stop("Invalid Minnesota boundary geometry or CRS.")
}
boundary <- sf::st_sf(geometry = sf::st_union(sf::st_transform(boundary, 26915)))

output_paths <- c(
  file.path(gep_dir, "gep-timber-2005.tif"),
  file.path(root_dir, "data", "forest-timber", "timber-stumpage-prices-2005-2021.csv"),
  file.path(map_dir, "timber-value-percent-change-2005-2025.png"),
  file.path(map_dir, "timber-value-percent-change-2005-2025.pdf"),
  file.path(figure_dir, "timber-valuation-time-series-2005-2025.png"),
  file.path(figure_dir, "timber-valuation-time-series-2005-2025.pdf"),
  file.path(table_dir, "timber-valuation-time-series-2005-2025.csv")
)
if (any(file.exists(output_paths))) {
  stop("Refusing to overwrite timber history outputs: ",
       paste(output_paths[file.exists(output_paths)], collapse = ", "))
}

evaluations <- data.table::fread(file.path(fia_dir, "MN_POP_EVAL.csv"),
  select = c("CN", "STATECD", "EVALID", "START_INVYR", "END_INVYR"),
  showProgress = FALSE)
eval_types <- data.table::fread(file.path(fia_dir, "MN_POP_EVAL_TYP.csv"),
  select = c("EVAL_CN", "EVAL_TYP"), showProgress = FALSE)
evals <- merge(evaluations[STATECD == 27L], unique(eval_types[, .(EVAL_CN, EVAL_TYP)]),
  by.x = "CN", by.y = "EVAL_CN", all.x = TRUE, sort = FALSE)
required_types <- c("EXPGROW", "EXPMORT", "EXPREMV")
for (i in seq_along(target_years)) {
  target_evalid <- unique(evals[END_INVYR == target_years[i] &
    EVAL_TYP %in% required_types, EVALID])
  if (length(target_evalid) != 1L ||
      !all(required_types %in% evals[EVALID %in% target_evalid, EVAL_TYP])) {
    stop("No unique complete FIA timber evaluation for target year ", target_years[i])
  }
}

plots <- data.table::fread(file.path(fia_dir, "MN_PLOT.csv"),
  select = c("CN", "STATECD", "UNITCD", "COUNTYCD", "PLOT", "MEASYEAR",
    "MEASMON", "MEASDAY", "INVYR", "PLOT_STATUS_CD", "LAT", "LON"),
  colClasses = list(character = c("CN", "STATECD", "UNITCD", "COUNTYCD", "PLOT")),
  showProgress = FALSE)
conditions <- data.table::fread(file.path(fia_dir, "MN_COND.csv"),
  select = c("PLT_CN", "COND_STATUS_CD", "CONDPROP_UNADJ"),
  colClasses = list(character = "PLT_CN"), showProgress = FALSE)
forest_share <- conditions[COND_STATUS_CD == 1L, .(
  forest_share = sum(CONDPROP_UNADJ, na.rm = TRUE)
), by = PLT_CN]
plots[, forest_share := forest_share$forest_share[match(CN, forest_share$PLT_CN)]]
for (column in c("MEASYEAR", "MEASMON", "MEASDAY", "INVYR", "PLOT_STATUS_CD")) {
  data.table::set(plots, j = column, value = as.integer(plots[[column]]))
}
for (column in c("LAT", "LON")) data.table::set(plots, j = column, value = as.numeric(plots[[column]]))
for (column in c("MEASMON", "MEASDAY")) data.table::set(plots, which(is.na(plots[[column]])),
  j = column, value = 0L)

rate_cache_path <- file.path(Sys.getenv("TEMP"), "mfrc-timber-history-rates-cache.rds")
cache_is_fresh <- file.exists(rate_cache_path) &&
  difftime(Sys.time(), file.info(rate_cache_path)$mtime, units = "hours") < 24
if (cache_is_fresh) {
  cached_rates <- readRDS(rate_cache_path)
  net <- cached_rates$net
  saw <- cached_rates$saw
  message("Reusing cached FIA annualized timber rates.")
} else {
  message("Loading FIA tables for historical timber snapshots.")
  fia <- rFIA::readFIA(dir = fia_dir, states = "MN", inMemory = FALSE)
  message("Calculating historical growing-stock volume rates.")
  net <- data.table::as.data.table(rFIA::growMort(
    fia, byPlot = TRUE, bySpecies = TRUE, landType = "forest",
    treeType = "gs", method = "annual", stateVar = "NETVOL"
  ))
  message("Calculated historical net cubic-foot rates.")
  saw <- data.table::as.data.table(rFIA::growMort(
    fia, byPlot = TRUE, bySpecies = TRUE, landType = "forest",
    treeType = "gs", method = "annual", stateVar = "SAWVOL_BF"
  ))
  message("Calculated historical sawtimber rates.")
  saveRDS(list(net = net, saw = saw), rate_cache_path, compress = FALSE)
}
net[, `:=`(plot_id = format_fia_id(PLT_CN), YEAR = as.integer(YEAR))]
saw[, `:=`(plot_id = format_fia_id(PLT_CN), YEAR = as.integer(YEAR))]
prices_used <- data.table::rbindlist(lapply(price_years, historical_prices), fill = TRUE)

build_plot_rates <- function(target_year, price_year) {
  start_year <- target_year - 4L
  n <- net[YEAR >= start_year & YEAR <= target_year]
  s <- saw[YEAR >= start_year & YEAR <= target_year]
  ids <- unique(c(n$plot_id, s$plot_id))
  visits <- plots[MEASYEAR >= start_year & MEASYEAR <= target_year & CN %in% ids]
  if (!nrow(visits)) stop("No FIA plot visits for ", target_year)
  keys <- c("STATECD", "UNITCD", "COUNTYCD", "PLOT")
  data.table::setorderv(visits, c(keys, "MEASYEAR", "MEASMON", "MEASDAY", "INVYR", "CN"))
  selected <- visits[, .SD[.N], by = keys]
  selected <- selected[PLOT_STATUS_CD == 1L & is.finite(forest_share) &
    forest_share > 0 & is.finite(LAT) & is.finite(LON) &
    LAT >= -90 & LAT <= 90 & LON >= -180 & LON <= 180]
  ids <- selected$CN
  n <- n[plot_id %in% ids]
  s <- s[plot_id %in% ids]
  message(target_year, " selected ", length(ids), " plots; rate records: ",
    nrow(n), " net and ", nrow(s), " sawtimber.")
  if (!nrow(n) || !nrow(s)) stop("Missing annualized growth or sawtimber records for ", target_year)
  nr <- c("PLT_CN", "YEAR", "SPCD", "COMMON_NAME", "RECR_NETVOL_ACRE",
    "MORT_NETVOL_ACRE", "REMV_NETVOL_ACRE", "GROW_NETVOL_ACRE")
  sr <- c("PLT_CN", "YEAR", "SPCD", "COMMON_NAME", "RECR_SAWVOL_BF_ACRE",
    "MORT_SAWVOL_BF_ACRE", "REMV_SAWVOL_BF_ACRE", "GROW_SAWVOL_BF_ACRE")
  if (!all(nr %in% names(n)) || !all(sr %in% names(s))) {
    stop("Unexpected rFIA field set.")
  }
  n <- n[, .(plot_id, YEAR, SPCD = as.integer(SPCD), COMMON_NAME,
    growth_cf_acre = GROW_NETVOL_ACRE, recruitment_cf_acre = RECR_NETVOL_ACRE,
    mortality_cf_acre = MORT_NETVOL_ACRE, removals_cf_acre = REMV_NETVOL_ACRE)]
  s <- s[, .(plot_id, YEAR, SPCD = as.integer(SPCD), COMMON_NAME,
    growth_bf_acre = GROW_SAWVOL_BF_ACRE, recruitment_bf_acre = RECR_SAWVOL_BF_ACRE,
    mortality_bf_acre = MORT_SAWVOL_BF_ACRE, removals_bf_acre = REMV_SAWVOL_BF_ACRE)]
  z <- merge(n, s, by = c("plot_id", "YEAR", "SPCD"), all = TRUE, suffixes = c("", "_saw"))
  for (field in c("growth_cf_acre", "recruitment_cf_acre", "mortality_cf_acre",
    "removals_cf_acre", "growth_bf_acre", "recruitment_bf_acre",
    "mortality_bf_acre", "removals_bf_acre")) {
    data.table::set(z, which(is.na(z[[field]])), field, 0)
    if (any(!is.finite(z[[field]]))) stop("Nonfinite FIA timber rate: ", field)
  }
  z[, COMMON_NAME := data.table::fifelse(
    is.na(COMMON_NAME) | COMMON_NAME == "", COMMON_NAME_saw, COMMON_NAME
  )]
  z <- z[plot_id %in% ids]
  z[, `:=`(
    sustainable_cords_ha = pmax(
      pmax(growth_cf_acre + recruitment_cf_acre - mortality_cf_acre, 0) / 80 -
        pmax(growth_bf_acre + recruitment_bf_acre - mortality_bf_acre, 0) /
          1000 * 2,
      0
    ) * 2.4710538147,
    sustainable_mbf_ha = pmax(growth_bf_acre + recruitment_bf_acre -
      mortality_bf_acre, 0) / 1000 * 2.4710538147
  )]
  g <- species_price_groups(z$SPCD, z$COMMON_NAME)
  z[, `:=`(
    species_group = g$species_group,
    cord_group = g$cord_group,
    saw_group = g$saw_group
  )]
  if (target_year == 2005L) {
    z[cord_group %in% c("W Spruce", "B Spruce"), cord_group := "Spruce Species"]
    z[cord_group %in% c("Red Pine", "White Pine"), cord_group := "Red/White Pine"]
    z[species_group %in% c("W Spruce", "B Spruce"), saw_group := "Spruce Species"]
    z[species_group %in% c("Red Pine", "White Pine"), saw_group := "Red/White Pine"]
    z[species_group %in% c("Balsam fir", "Aspen", "Balm"), saw_group := species_group]
  }
  price_rows <- which(prices_used[["price_year"]] == price_year)
  p <- prices_used[price_rows]
  if (!nrow(p)) stop("No DNR stumpage prices found for price year ", price_year)
  if (anyDuplicated(p[, .(product, species_group)])) {
    stop("Duplicate DNR product/species price keys for price year ", price_year)
  }
  cp <- p[product == "cord", .(cord_group = species_group, cord_price = price_per_unit)]
  sp <- p[product == "sawtimber", .(saw_group = species_group, saw_price = price_per_unit)]
  z <- merge(z, cp, by = "cord_group", all.x = TRUE)
  z <- merge(z, sp, by = "saw_group", all.x = TRUE)
  if (anyNA(z$cord_price) || anyNA(z$saw_price)) {
    stop("Missing annual DNR price crosswalk for target year ", target_year)
  }
  z[, `:=`(
    value_nominal_usd_ha = sustainable_cords_ha * cord_price +
      sustainable_mbf_ha * saw_price
  )]
  rates <- z[, .(
    sustainable_cords_ha = sum(sustainable_cords_ha),
    sustainable_mbf_ha = sum(sustainable_mbf_ha),
    value_nominal_usd_ha = sum(value_nominal_usd_ha)
  ), by = plot_id]
  selected <- merge(selected, rates, by.x = "CN", by.y = "plot_id", all = FALSE)
  coordinates <- sf::st_transform(
    sf::st_as_sf(selected, coords = c("LON", "LAT"), crs = 4326, remove = FALSE), 26915
  )
  inside <- lengths(sf::st_intersects(coordinates, boundary)) > 0L
  result <- data.table::as.data.table(sf::st_drop_geometry(coordinates[inside, ]))
  if (nrow(result) < 32L) stop("Fewer than 32 valid Minnesota FIA plots in ", target_year)
  xy <- sf::st_coordinates(coordinates[inside, ])
  result[, `:=`(x = xy[, 1], y = xy[, 2])]
  result
}

series <- vector("list", length(target_years))
rates_by_year <- vector("list", length(target_years))
for (i in seq_along(target_years)) {
  year <- target_years[i]
  price_year <- price_years[i]
  rates <- build_plot_rates(year, price_year)
  rates_by_year[[i]] <- rates
  forest_path <- file.path(forest_dir, sprintf("forest-layer-%d.tif", year))
  metadata_path <- sub("\\.tif$", ".json", forest_path)
  if (!file.exists(forest_path) || !file.exists(metadata_path)) {
    stop("Missing historical forest layer or metadata for ", year)
  }
  forest <- terra::rast(forest_path)[[1]]
  metadata <- jsonlite::fromJSON(metadata_path)
  if (metadata$target_year != year || !terra::same.crs(forest, "EPSG:26915") ||
      any(terra::res(forest) != 300)) stop("Unexpected forest layer for ", year)
  area_ha <- sum(!is.na(terra::values(forest, mat = FALSE))) * 9
  series[[i]] <- data.frame(
    target_year = year, measurement_window_start = year - 4L,
    measurement_window_end = year, price_year = price_year,
    price_currency = "year-specific nominal USD", FIA_plot_visits = nrow(rates),
    mapped_forest_area_ha = area_ha,
    mean_annual_cords_per_ha = mean(rates$sustainable_cords_ha),
    total_annual_cords_extrapolated = mean(rates$sustainable_cords_ha) * area_ha,
    mean_annual_sawtimber_Mbf_per_ha = mean(rates$sustainable_mbf_ha),
    total_annual_sawtimber_Mbf_extrapolated = mean(rates$sustainable_mbf_ha) * area_ha,
    mean_annual_stumpage_value_nominal_USD_per_ha = mean(rates$value_nominal_usd_ha),
    total_annual_stumpage_value_nominal_USD_extrapolated =
      mean(rates$value_nominal_usd_ha) * area_ha,
    extrapolation = "Unweighted FIA plot mean times mapped forest area; descriptive, not design-based"
  )
  message("Prepared historical timber snapshot ", year, " from ", nrow(rates), " FIA plot visits.")
}
series <- do.call(rbind, series)

message("Fitting 2005 spatial value surface for the percent-change map.")
first <- rates_by_year[[1]]
first_sites <- sf::st_as_sf(
  first[, .(x, y, annual_rate = mean(value_nominal_usd_ha)), by = .(x, y)],
  coords = c("x", "y"), crs = 26915
)
model <- fit_timber_variogram(first_sites)$model
forest_2005 <- terra::rast(file.path(forest_dir, "forest-layer-2005.tif"))[[1]]
cells_2005 <- which(!is.na(terra::values(forest_2005, mat = FALSE)))
prediction_2005 <- krige_timber_rate(first_sites, forest_2005, cells_2005, model)
value_2005 <- terra::rast(forest_2005)
terra::values(value_2005) <- NA_real_
value_2005[cells_2005] <- prediction_2005$rate
value_2025 <- terra::rast(file.path(gep_dir, "gep-timber-2025.tif"))[[1]]
terra::compareGeom(value_2005, value_2025, stopOnError = TRUE)
change <- terra::ifel(
  !is.na(value_2005) & !is.na(value_2025) & value_2005 > 0,
  100 * (value_2025 - value_2005) / value_2005, NA
)
if (!is.finite(terra::global(change, "min", na.rm = TRUE)[1, 1])) {
  stop("No shared positive-baseline timber value pixels for the change map.")
}
terra::writeRaster(value_2005, output_paths[1], overwrite = FALSE,
  datatype = "FLT8S", NAflag = -9999,
  gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES"))

data.table::fwrite(prices_used, output_paths[2])
utils::write.csv(series, output_paths[7], row.names = FALSE)
names(change) <- "nominal_stumpage_value_percent_change_2005_2025"
tmap::tmap_mode("plot")
map <- tmap::tm_shape(boundary, crs = 26915, is.main = TRUE) +
  tmap::tm_polygons(fill = "#eeeeee", col = NA) +
  tmap::tm_shape(change) +
  tmap::tm_raster(
    col.scale = tmap::tm_scale_intervals(
      breaks = c(-Inf, -50, -10, 10, 50, 100, 200, Inf), midpoint = NA,
      labels = c("< -50%", "-50 to -10%", "-10 to +10%", "+10 to +50%",
        "+50 to +100%", "+100 to +200%", ">= +200%"),
      values = c("#b2182b", "#ef8a62", "#f7f7f7", "#d1e5f0",
        "#92c5de", "#4393c3", "#2166ac"), value.na = NA
    ),
    col.legend = tmap::tm_legend(title = "Value change\n% of 2005 value", na.show = FALSE)
  ) +
  tmap::tm_shape(boundary) + tmap::tm_borders(col = "#555555", lwd = 0.7) +
  tmap::tm_title("Annual timber stumpage-value change: 2005-2025",
    size = 0.85, position = tmap::tm_pos_out("center", "top")) +
  tmap::tm_credits(
    "Species/product stumpage prices by FIA snapshot year (nominal USD).\nShared forest pixels; public FIA coordinates; exploratory kriging.",
    size = 0.65, position = tmap::tm_pos_out("left", "bottom")
  ) +
  tmap::tm_layout(frame = FALSE, meta.margins = c(0.12, 0, 0.07, 0.25),
    legend.position = tmap::tm_pos_out("right", "center"),
    legend.title.size = 0.8, legend.text.size = 0.7)
for (extension in c("png", "pdf")) {
  path <- file.path(map_dir, paste0("timber-value-percent-change-2005-2025.", extension))
  tmap::tmap_save(map, path, width = 7, height = 7, units = "in", dpi = 300)
  if (!file.exists(path) || file.info(path)$size <= 0) stop("Map export failed: ", path)
}

draw_series <- function() {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(3, 1), mar = c(4, 5, 2.5, 1), oma = c(2, 0, 1, 0), las = 1)
  graphics::plot(series$target_year, series$total_annual_cords_extrapolated / 1e6,
    type = "o", pch = 16, col = "#238443", lwd = 2,
    xlab = "", ylab = "Million cords / year", xaxt = "n",
    main = "Annual sustainable cord-product quantity")
  graphics::axis(1, at = series$target_year)
  graphics::plot(series$target_year, series$total_annual_sawtimber_Mbf_extrapolated / 1e6,
    type = "o", pch = 16, col = "#2166ac", lwd = 2,
    xlab = "", ylab = "Million Mbf / year", xaxt = "n",
    main = "Annual sustainable sawtimber quantity")
  graphics::axis(1, at = series$target_year)
  graphics::plot(series$target_year,
    series$total_annual_stumpage_value_nominal_USD_extrapolated / 1e6,
    type = "o", pch = 16, col = "#b35806", lwd = 2,
    xlab = "FIA five-year window endpoint; value in year's nominal USD",
    ylab = "Million USD / year", xaxt = "n",
    main = "Gross annual stumpage-value proxy")
  graphics::axis(1, at = series$target_year)
}
for (extension in c("png", "pdf")) {
  path <- file.path(figure_dir, paste0("timber-valuation-time-series-2005-2025.", extension))
  if (extension == "png") grDevices::png(path, width = 2100, height = 2400, res = 300)
  else grDevices::pdf(path, width = 7, height = 8)
  draw_series()
  grDevices::dev.off()
  if (!file.exists(path) || file.info(path)$size <= 0) stop("Series export failed: ", path)
}
if (file.exists(rate_cache_path)) unlink(rate_cache_path)
message("Saved timber percentage-change map, five-point series, and source price table.")
