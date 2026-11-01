required_packages <- c("rFIA", "data.table")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install required packages before running: install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_argument) != 1L) {
  stop("Run this script with Rscript scripts/task-pull-fia-data.R.")
}
script_path <- normalizePath(sub("^--file=", "", script_argument), mustWork = TRUE)
root_dir <- dirname(dirname(script_path))
data_dir <- file.path(root_dir, "data", "FIA")
metadata_dir <- file.path(data_dir, "metadata")
requested_years <- 2000:2025

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(data_dir)) {
  stop("Could not create FIA output directory: ", data_dir)
}
options(timeout = max(3600, getOption("timeout")))

# Keep all historical measurements and population tables needed for estimation.
rFIA::getFIA(
  states = "MN",
  dir = data_dir,
  common = FALSE,
  load = FALSE,
  nCores = 1
)

csv_files <- list.files(
  data_dir, pattern = "^MN_.*\\.csv$", full.names = TRUE, ignore.case = TRUE
)
read_table <- function(table) {
  path <- csv_files[
    toupper(basename(csv_files)) == paste0("MN_", table, ".CSV")
  ]
  if (length(path) != 1L || file.info(path)$size == 0) {
    stop("Missing or empty downloaded table: MN_", table, ".csv")
  }
  data.table::fread(path, colClasses = "character", showProgress = FALSE)
}

evaluations <- read_table("POP_EVAL")
evaluation_types <- read_table("POP_EVAL_TYP")
required_columns <- c("CN", "STATECD", "EVALID", "START_INVYR", "END_INVYR")
if (!all(required_columns %in% names(evaluations))) {
  stop("POP_EVAL is missing required evaluation-year columns.")
}
if (!all(c("EVAL_CN", "EVAL_TYP") %in% names(evaluation_types))) {
  stop("POP_EVAL_TYP is missing required evaluation-type columns.")
}
if (any(is.na(evaluations$STATECD) | evaluations$STATECD != "27")) {
  stop("Downloaded evaluations contain missing or non-Minnesota state codes.")
}

evaluation_years <- as.integer(evaluations$END_INVYR)
if (anyNA(evaluation_years)) {
  stop("POP_EVAL contains missing or invalid END_INVYR values.")
}
selected_evaluations <- evaluations[evaluation_years %in% requested_years]
if (!nrow(selected_evaluations)) {
  stop("No Minnesota evaluations available for 2000-2025.")
}
selected_evaluations <- merge(
  selected_evaluations,
  unique(evaluation_types[, c("EVAL_CN", "EVAL_TYP"), with = FALSE]),
  by.x = "CN", by.y = "EVAL_CN", all.x = TRUE, sort = FALSE
)
if (anyNA(selected_evaluations$EVAL_TYP)) {
  stop("Selected evaluations have no matching evaluation types.")
}
data.table::setorderv(selected_evaluations, c("END_INVYR", "EVALID", "EVAL_TYP"))

coverage <- data.table::data.table(
  year = requested_years,
  evaluation_count = vapply(
    requested_years,
    function(year) length(unique(evaluations$EVALID[evaluation_years == year])),
    integer(1)
  )
)
coverage[, available := evaluation_count > 0L]
dir.create(metadata_dir, showWarnings = FALSE)
if (!dir.exists(metadata_dir)) {
  stop("Could not create FIA metadata directory: ", metadata_dir)
}
data.table::fwrite(
  selected_evaluations, file.path(metadata_dir, "MN_evaluations_2000_2025.csv")
)
data.table::fwrite(
  coverage, file.path(metadata_dir, "MN_year_coverage_2000_2025.csv")
)
file_sizes <- file.info(csv_files)$size
if (anyNA(file_sizes) || any(file_sizes == 0)) {
  stop("The downloaded database contains missing or empty CSV files.")
}
data.table::fwrite(
  data.table::data.table(file = basename(csv_files), bytes = file_sizes),
  file.path(metadata_dir, "MN_download_files.csv")
)
data.table::fwrite(
  data.table::data.table(
    state = "MN",
    first_requested_year = min(requested_years),
    last_requested_year = max(requested_years),
    year_field = "POP_EVAL.END_INVYR",
    source_url = "https://apps.fs.usda.gov/fia/datamart/CSV/MN_CSV.zip",
    downloaded_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    rFIA_version = as.character(utils::packageVersion("rFIA"))
  ),
  file.path(metadata_dir, "MN_download_metadata.csv")
)

missing_years <- coverage[available == FALSE, year]
if (length(missing_years)) {
  warning(
    "No published Minnesota evaluations for requested years: ",
    paste(missing_years, collapse = ", "),
    ". See metadata/MN_year_coverage_2000_2025.csv. No years were fabricated."
  )
}
message(
  "Saved ", length(csv_files), " Minnesota FIA tables to ", data_dir,
  ". Evaluation coverage: ", sum(coverage$available), "/",
  length(requested_years), " requested years."
)