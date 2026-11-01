# README: Minnesota Forest Resources Council Project
Ryan McWay

# Overview

This project contains resources and scripts related to the Minnesota
Forest Resources Council (MFRC). The purpose of this project is to
assess the monetary and non-monetary value of Minnesota family forests.
The main approach for this project involves collecting and analyzing
data on Minnesota family forests to determine their economic and
ecological contributions within the concept of ecosystem services.

# Structure of the Project

The project is organized into the following directory format, which
contains data, scripts, and logs that facilitate data collection,
analysis, and reporting. Key components include:

- `data/`: Contains raw and processed data related to Minnesota family
  forests.
- `scripts/`: Includes scripts for data processing, analysis, and
  visualization.
- `logs/`: Stores logs for issues and questions encountered during the
  project.
- `README.qmd`: Provides an overview and documentation for the project.

# Usage

To use this project, follow these steps:

1.  Clone the repository to your local machine.
2.  Ensure that all required dependencies for the scripts are installed.
3.  Place any new raw data files in the `data/` directory.
4.  Within the `scripts/` directory, run the `scripts/run-pipeline.sh`
    script to execute the data processing and analysis workflow.
5.  While in development, I am using `scripts/plan-analysis.qmd` to plan
    and document the analysis steps.
6.  Record any issues or questions in the `logs/` directory that arise
    during the execution of the scripts.

# Downloading Minnesota FIA data

Install the R packages once with `install.packages(c("rFIA", "data.table"))`,
then run `Rscript scripts/task-pull-fia-data.R` from the repository root.
The script downloads the complete Minnesota FIADB CSV database using
`rFIA::getFIA()` and saves it in `data/FIA/`.

The requested 2000-2025 range refers to annual estimation years
(`POP_EVAL.END_INVYR`), not plot measurement years. Historical supporting
records are retained rather than filtered by `INVYR`. Evaluation IDs, types,
year availability, and download metadata are saved under `data/FIA/metadata/`.
The script warns if any requested years have no published evaluations.
Rerunning the script refreshes the downloaded tables and metadata.
The pipeline also runs this download step before rendering the logs and
records its runtime and status. Both `Rscript` and `quarto` must be on `PATH`.

For analysis, load the database using
`db <- rFIA::readFIA("data/FIA", states = "MN", inMemory = FALSE)`.
Read `data/FIA/metadata/MN_evaluations_2000_2025.csv` and pass its unique
`EVALID` values to
`rFIA::clipFIA(db, mostRecent = FALSE, evalid = evaluation_ids)`.
Availability of specific estimates depends on the evaluation types published
for each year; an available year does not imply every estimate type is available.

# Creating exploratory forest layers

The pipeline runs `scripts/task-create-forest-layer.py` after downloading FIA
data and requires `python` on `PATH`. Install its dependencies with
`python -m pip install -r scripts/requirements-forest-layer.txt`, then run
`python scripts/task-create-forest-layer.py`. No raster generation was run
during implementation. Existing forest-layer outputs must be reviewed and
moved or removed before rerunning this pipeline step.

It creates `data/forest-layer/forest-layer-[year].tif` for 2005, 2010, 2015,
2020, and 2025. Each uses the newest physical-plot observation in the trailing
five measurement years (`MEASYEAR`), not annual evaluation samples.
Forest-condition `BALIVE` is weighted by `CONDPROP_UNADJ` and divided by the
summed forest proportion, yielding basal area in square feet per forest acre.
Nonforest/nonsampled plots and plots with invalid forest attributes are excluded;
coincident public coordinates are averaged.

The script fits an exponential variogram using up to 1,000 reproducibly sampled
sites, then performs local ordinary kriging using the nearest 32 of all eligible
sites. Outputs use a common 300-metre EPSG:26915 grid clipped to the existing
Minnesota county boundary. Only ESA forest classes 50, 60, 61, 62, 70, 71, 72,
80, 81, 82, 90, and 100 are predicted; mosaics and flooded forest are excluded.
Categorical masks are reprojected with nearest-neighbor resampling.

The default ESA source is
`G:\Shared drives\NatCapTEEMs\Files\base_data\lulc\esa`.
Masks match the target years except that 2025 explicitly uses
`lulc_esa_2022.tif`. Override input locations using `--fia-dir`, `--esa-dir`,
`--boundary`, and `--output-dir`. Missing inputs fail explicitly, and existing
outputs are not overwritten.

Band 1 is nonnegative basal area; band 2 is latent kriging variance in squared
basal-area units. Nonforest, out-of-state, and unavailable mask pixels are
NoData (`-9999`), not zero. JSON sidecars record sampling, mask years, variogram
parameters, zero-clipping counts, and maximum nearest-sample distance.

These are exploratory five-year composites, not design-based rFIA estimates
or exact annual maps. Public FIA coordinates are fuzzed/swapped, and 300-metre
pixels do not imply 300-metre accuracy. All masked forest pixels are predicted,
including extrapolation beyond sampled locations. Variograms and spatial
cross-validation require review before inference; band 2 does not include
land-cover, temporal, or coordinate uncertainty.

# Blame

If you encounter any issues or have questions regarding the project,
please refer to the logs in the `logs/` directory or contact the project
maintainer – Ryan McWay (mcway005@umn.edu). We wish to acknowledge that
generative AI tools were used to assist in the creation and editing of
this project.
