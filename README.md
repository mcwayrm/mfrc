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

The pipeline runs the implemented forest-layer, ownership, carbon, timber,
biodiversity, FEMA, and ESVD valuation and visualization tasks before rendering the
logs. Air- and water-quality task files are planning notes and are not run.
Install each task's required packages and provide its external input data first.
Several tasks refuse to overwrite existing outputs, so a rerun stops unless those
outputs are moved or removed deliberately; the pipeline does not delete them.

The ESVD benchmark requires the forest and all-ecoregion CSV exports in
`data/EVSD/`. It uses 2020 international dollars and is an illustrative,
non-additive benefit-transfer benchmark, not a Minnesota-specific valuation.
Its selection rules, service summaries, maps, and limitations are documented in
the ESVD section of `papers/report.tex`.

# Downloading Minnesota FIA data

Install the R packages once with `install.packages(c("rFIA", "data.table"))`,
then run `Rscript scripts/task-pull-fia-data.R` from the repository root.
Alternatively, in an R console, run
`source("f:\\mfrc\\scripts\\task-pull-fia-data.R", encoding = "UTF-8")`
(adjust the absolute path for your checkout). Both methods locate the output
directory relative to the script, regardless of the working directory.
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

If Conda stalls while solving dependencies, activate the intended environment
and run `python -m pip install --only-binary=:all: -r scripts/requirements-forest-layer.txt`
from the repository root. This uses prebuilt wheels and leaves already-satisfied
dependencies unchanged. If pip reports an `ssl.SSLError` while loading the
Windows certificate store, add `--use-deprecated=legacy-certs` to use pip's
CA bundle instead; HTTPS certificate verification remains enabled.
Check imports without generating outputs with
`python scripts/task-create-forest-layer.py --help`.

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

# Visualizing FIA and the forest layer

Install R packages `terra`, `sf`, `tmap` (version 4 or newer), `data.table`,
and `jsonlite`, then run `Rscript scripts/task-visualize-forest-layer.R`.
The script also supports `source()` with its absolute file path, regardless
of the working directory. It reads all five basal-area rasters and their JSON
sidecars, the FIA plot table, and the Minnesota county boundary.

Four individual maps and `forest-layer-four-panel` are saved as PNG (300 dpi)
and PDF in `outputs/maps/`:

- **A:** Distinct measurement years per physical FIA plot during 2000-2025,
  keyed by `STATECD`, `UNITCD`, `COUNTYCD`, and `PLOT`. All plot statuses are
  included. The latest public coordinates are shown without jitter; invalid
  or out-of-state locations are excluded with warnings.
- **B:** 2025 kriged basal-area density in square feet per forest acre,
  using 2021-2025 measurements and the 2022 ESA mask.
- **C:** The exact 2022 ESA forest footprint preserved by the 2025 model.
  The producer predicts every eligible forest pixel, so band 1's valid cells
  preserve the mask without introducing differences from a second reprojection.
  Gray means outside the forest mask, including unavailable land cover; it
  is not a separate estimate of nonforest area. Forest classes are
  50, 60, 61, 62, 70, 71, 72, 80, 81, 82, 90, and 100.
- **D:** `100 * (2025 basal area - 2005 basal area) / 2005 basal area`.
  Only pixels forested in both saved masks with positive 2005 basal area
  are compared. Zero baselines and nonoverlapping forest pixels are excluded,
  not assigned zero change. This is density change, not forest-cover change.

Three CSVs in `outputs/tables/` summarize measurement frequency, annual mapped
forest area and basal-area distributions, and shared-pixel change/exclusions.
Areas use the projected 300-metre grid (9 hectares per cell); density summaries
are equal-cell statistics over each year's forest mask, not design-based FIA
estimates. Mean pixel percentage change is not percentage change in mean
density, and annual means use different forest footprints. These maps are
exploratory five-year composites, not exact annual ecosystem-service flows.
Small baselines amplify percentage changes. Rerunning replaces only the named
visualization outputs, leaving input data unchanged.

Useful report companions would be an absolute basal-area change map (less
sensitive to small baselines), a shared-forest change histogram, and a kriging
uncertainty or nearest-sample-distance map. Kriging variance alone does not
capture coordinate, land-cover, or temporal uncertainty.

# Valuing live aboveground carbon storage

The forest-layer script now also produces `forest-carbon-stock-[year].tif`
and JSON sidecars in `data/forest-layer/`, while retaining the basal-area
outputs. The carbon layers use the same sampling windows, ESA masks, grid,
and kriging method, with independently fitted carbon variograms. Live forest
trees (`STATUSCD=1` and forest conditions) supply
`sum(TREE.CARBON_AG * TREE.TPA_UNADJ)` in pounds C per plot acre. Divide by
the summed forest-condition proportion once, then multiply by
`0.00045359237 / 0.40468564224` for metric tonnes C per forest hectare.
Forest plots without live-tree records contribute zero; invalid live-tree
carbon or expansion values exclude the plot. Only live aboveground tree
carbon is represented, not soil, deadwood, litter, roots, or understory.

After review, install `terra`, `sf`, `tmap` (>= 4.0), and `jsonlite` in R and run
`Rscript scripts/task-value-carbon-storage.R`, or run it as part of the pipeline.
It also supports `source()` with its absolute file path. It requires the carbon rasters and
their metadata; basal-area rasters are not accepted as carbon quantities.

Outputs are `data/GEP/gep-carbon-stock-[year].tif` for 2005, 2010, 2015,
2020, and 2025, in **2007 USD per hectare**, computed as
`carbon tonnes C/ha * (44/12) * SCC USD/metric tonne CO2`.
The grid and NoData forest mask are preserved; no pixel-area multiplier
is used for per-hectare values.

The [EPA archived table](https://19january2017snapshot.epa.gov/climatechange/social-cost-carbon_.html)
supplies the 3% average SCC prices: $36 (2015), $42 (2020), and $46 (2025).
The user specified $21.40 for both 2005 and 2010, referencing the
[2010 TSD](https://www.epa.gov/sites/default/files/2016-12/documents/scc_tsd_2010.pdf).
Applying that value to 2005 is an explicit assumption, not a published 2005
entry from the linked snapshot table. Prices remain in 2007 dollars, with
no inflation adjustment. A price CSV and per-year JSON sidecars record the
formula, sources, assumptions, and upstream interpolation metadata.
Existing GEP outputs are not overwritten.

The script finishes by creating two maps in `outputs/maps/`: a 2025 heat map
of stock value per hectare (`carbon-stock-value-hotspots-2025`) and a
2005-2025 percentage value-change map
(`carbon-stock-value-percent-change-2005-2025`). The heat map shows high-value
locations, not statistically tested hotspots. Change is computed only on
shared forest pixels with positive 2005 value; zero baselines and nonoverlapping
forest are excluded. Value change includes SCC price changes as well as carbon
changes; it is not a map of sequestration or carbon-stock change alone.

`outputs/figures/carbon-stock-value-time-series` shows aligned quantity and
total-value panels, plus a fixed-$21.40 SCC comparison that isolates quantity
and area changes from price changes. All three visualizations are saved as PNG
(300 pixels per inch at the intended print size) and PDF.
`outputs/tables/carbon-stock-value-time-series.csv` contains the five annual
totals, areas, prices, and fixed-price comparison;
`carbon-stock-value-change-summary.csv` records shared-pixel change and exclusions.
Totals multiply per-hectare quantities and values by 9 hectares per projected
300-metre cell. Each year uses its own forest footprint, so temporal differences
include mapped-area changes. The 2025 layer uses ESA 2022. These are statewide
mapped-forest totals, not yet family-ownership estimates or design-based FIA
estimates. Existing valuation rasters are protected; visualization files are
replaced when the visualization function is run.

This is an SCC-valued **storage stock proxy**, not annual sequestration GEP.
SCC describes marginal damages from emissions; valuing the entire stock
implicitly assumes avoided release and does not establish that counterfactual.
Do not sum repeated stock maps as annual benefits. The public FIA coordinate,
interpolation, and five-year-composite limitations still apply.

# Valuing annualized live-carbon change

Run `Rscript scripts/task-value-carbon-sequestration.R`, or `source()` the
absolute script path in R. Required packages are `terra`, `sf`, `tmap` (>= 4),
`jsonlite`, `data.table`, and `gstat`. This is a separate matched-FIA analysis,
not subtraction of the independently kriged carbon-stock maps.

The script links visits using `PLOT.PREV_PLT_CN`, verifies physical plot keys,
and requires sampled status and identical `DESIGNCD` at both visits. Both
visits must have only forest conditions with forest shares within `1e-6` of
one; mixed forest/nonforest plots and land-use transitions are excluded.
Stocks use live forest-tree `CARBON_AG * TPA_UNADJ`, normalized by each visit's
forest share and converted to metric tonnes C/ha. A fully forested sampled
plot with no live trees has zero live-tree stock. Invalid stocks are excluded.
The rate is `(current C - previous C) / PLOT.REMPER`, using the actual FIA
remeasurement interval in years, not a fixed five-year divisor.

For targets 2005, 2010, 2015, 2020, and 2025, the newest endpoint in the
trailing five measurement years is selected per physical plot before
eligibility checks. Its predecessor may fall outside that window and the
interval may exceed five years. Early records do not support a comparable
2005 estimate under these rules: missing predecessors and design changes
prevent a valid comparison. No 2005 quantity or valuation is fabricated.
`outputs/tables/carbon-sequestration-year-coverage.csv` records availability,
and `carbon-sequestration-pair-exclusions.csv` records sequential exclusions.
Currently, supported targets are 2010, 2015, 2020, and 2025.

Each supported target fits an exponential variogram to up to 1,000
reproducibly selected matched sites (seed 42), with 20 equal-distance bins
and at least 30 pairs per bin. It uses bounded, pair-count-weighted normalized
least squares, following the stock-layer fitter's objective and bounds, with
an R optimizer. `gstat` supplies ordinary kriging and point prediction variance
(rather than the stock script's latent-surface variance). Predictions use the nearest 32
of all eligible sites. Coincident coordinates are averaged. Negative observed
rates and negative predictions are retained, not clipped. Fully forested
matched-plot rates are interpolated to all valid forest pixels of the saved
carbon-stock grid, including pixels with no nearby matched observations.
This conditional sampling and extrapolation require validation.

Two-band quantity rasters (`forest-carbon-sequestration-[year].tif`, rate
and prediction variance), matched-pair CSVs, and metadata are saved in
`data/forest-sequestration/`. Single-band valued rates
(`gep-carbon-sequestration-[year].tif`) and metadata are saved in `data/GEP/`.
Values are `rate_tonnes_C_per_ha_per_year * (44/12) * SCC`, in
**2007 USD/ha/year**, using the storage task's 3% average SCC schedule.
An endpoint SCC is applied to the annualized interval rate; this is not
year-by-year pricing of emissions. Per-hectare rates are multiplied by
9 hectares per cell only when forming statewide totals.

The script creates a signed annual-value map and an **absolute** value-change
map between the earliest and latest supported targets (currently 2010-2025),
plus aligned quantity/value time-series panels with a fixed-SCC comparison.
Percentage change is avoided because baseline sequestration can be zero or
negative. PNG/PDF outputs are in `outputs/maps/` and `outputs/figures/`;
`outputs/tables/carbon-sequestration-time-series.csv` includes signed totals,
positive accumulation, negative losses, areas, prices, and interval ranges.
Input stock layers and storage valuations remain unchanged. Existing quantity,
valuation, sample, and summary outputs are protected against overwriting.

These are **annualized net live-aboveground pool-change proxies**, not a
complete atmospheric sequestration balance: mortality and removals reduce
the live pool but may transfer carbon to deadwood or harvested products
rather than immediately to the atmosphere. Other pools, land-use transitions,
population expansion, family-ownership attribution, and spatial
cross-validation are not included. Do not report these as definitive
family-forest GEP or add the full stored-carbon value to annual flow values.

# Blame

If you encounter any issues or have questions regarding the project,
please refer to the logs in the `logs/` directory or contact the project
maintainer – Ryan McWay (mcway005@umn.edu). We wish to acknowledge that
generative AI tools were used to assist in the creation and editing of
this project.
