# DRI-RAD Eddy Covariance Processing Pipeline

This repository processes logger and EddyPro data for NICE NET (DRI) eddy
covariance sites: ECDP, EDVG, EDVP, ERVA, ERVP, ECSM, and ESTW. It runs meteorological
QAQC, multi-site meteorological gap filling, EC and LI-710 QAQC, REddyProc MDS
gap filling, ET calculations, and diagnostic figures.

An AmeriFlux review export is created automatically after each site's Level 4
step, including runs with no new L4 observations. It contains the complete saved
L4 period. `dry_run = TRUE` only previews this step; runs excluding L4 do not
export. Use `export_ameriflux = FALSE` to disable it explicitly.
The default review export includes existing EC `LE_PI_F`, `H_PI_F`, and `FC_PI_F`
(L4 QC, storage fallback and MDS; no EBR correction). These retain local PI names
and require submission mapping with the AmeriFlux Data Team. Other PI/corrected
columns are opt-in via `include_pi = TRUE`. Base meteorology
is from L1 QC before L2 filling. Bulk/tipping-bucket precipitation use
`P_1_1_1`/`P_1_1_2`; `P_F_1_1_1` is the separate L2 filled precipitation series.
EC LE/H use `_1_1_1`; LI-710 LE/H use `_1_1_2`, with their L4 MDS-filled results
separately labeled `LE_F_1_1_2`/`H_F_1_1_2`. Existing EC storage-inclusive PI
results are not relabeled as turbulent filled fluxes. Scientific calculations
are unchanged. Confirm sensor indices against the submitted site metadata.
To export existing Level 4 results independently, use
`source("run_level5_ameriflux.R")` followed by
`run_level5_ameriflux(base_dir = data_root)`. All seven sites are selected by
default. The exporter preserves the input files. See
[Level 5 instructions and variable mapping](AMERIFLUX_LEVEL5.md), including
the requested flux QC policy and the distinction between local `_PI` review
columns and final AmeriFlux submission labels.

EddyPro random uncertainties (`LE_rand_err`, `H_rand_err`, `FC_rand_err`,
`TAU_rand_err`) are retained in L3 and carried unchanged into L4. Existing older
rows without those fields read as NA. Reprocess L3_EC and L4, or the full pipeline,
to populate them from the original CSVs; incremental processing alone does not
backfill historical uncertainty. These are not filled-product uncertainty estimates.

## Software and package installation

Use R 4.4 or later and a recent RStudio Desktop. The pipeline requires:

| Package | Use |
| --- | --- |
| tidyverse | Data import, transformation, and plotting |
| lubridate | Date and timestamp handling |
| zoo | Short-gap interpolation |
| cowplot | Multi-panel figures |
| openeddy | EC despiking and QAQC utilities |
| bigleaf | Micrometeorological calculations |
| REddyProc | MDS gap filling |
| ggpmisc | Regression annotations |
| gridExtra | Figure layout |
| Metrics | Model metrics |

Open RStudio and run this installation block once. Internet access is required.

```r
cran_packages <- c(
  "tidyverse", "lubridate", "zoo", "cowplot", "bigleaf",
  "REddyProc", "ggpmisc", "gridExtra", "Metrics", "remotes"
)

missing_cran <- setdiff(cran_packages, rownames(installed.packages()))
if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

# openeddy is installed from GitHub because it is not distributed through CRAN.
if (!requireNamespace("openeddy", quietly = TRUE)) {
  remotes::install_github("lsigut/openeddy")
}
```

On Windows, install the Rtools release matching your R version if installation
reports that build tools are missing. Restart RStudio and rerun the block.

Verify the installation:

```r
required_packages <- c(
  "tidyverse", "lubridate", "zoo", "cowplot", "openeddy",
  "bigleaf", "REddyProc", "ggpmisc", "gridExtra", "Metrics"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) == 0) {
  message("All required packages are installed.")
} else {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}
```

Package installation is normally required only once. Do not put
`install.packages()` inside the processing scripts.

## Project layout

```text
DRI-RAD-EC_Proc_operational.Rproj   RStudio project
run_pipeline.R                     Main runner and options
run_all.R                          Immediate all-site run
run_example.R                      All-site and single-site examples
R/                                 Scientific stages and shared I/O
meta_files/file_directions_new.csv Site paths
meta_files/manual_flag.csv         Manual QAQC flags
meta_files/meta_G.csv              Soil heat-flux metadata
tests/                             Synthetic tests
docs/                              Scientific workflow documentation
```

Open `DRI-RAD-EC_Proc_operational.Rproj`. Do not source individual
`R/Level_*_new.R` files directly; they receive paths, sites, periods, and staging
information from `run_pipeline.R`.

## Configure the data root

Each collaborator can use a different data location. Set it in the command that
calls `run_pipeline()`; the runner source code does not need to be edited.

```r
# DRI network drive
data_root <- "Z:/NWI/Task_3_Water_Use_ET_and_Meteoroligical_Monitoring/networks/eddy_stations"

# Example on another computer
# data_root <- "D:/eddy_stations"
```

Use forward slashes. Paths in `meta_files/file_directions_new.csv` are resolved
relative to `base_dir`, unless they are already absolute drive-letter or UNC
paths.

EddyPro files are discovered below each configured site folder:

```text
<base_dir>/EddyPro_Output_Files_YK/<site_id>/YYYY/*full_output*.csv
<base_dir>/EddyPro_Output_Files_YK/<site_id>/YYYY_01/*full_output*.csv
<base_dir>/EddyPro_Output_Files_YK/<site_id>/YYYY_02/*full_output*.csv
```

The reader detects field names and units, converts `-9999` to missing values,
creates `TIMESTAMP` from `date` and `time`, and uses the latest run when exports
overlap. `dir_eddypro` in `meta_files/file_directions_new.csv` must point to the
site folder (for example `/EddyPro_Output_Files_YK/ESTW`). All matching CSVs below
that folder are read recursively; an `output` subfolder is not required. These
CSVs now supply the entire EC history. The legacy `dir_EC` column is unused by
this runner; no Excel EC summary is read or used as a fallback.

For duplicate observation timestamps, the latest run timestamp in the filename
(`YYYY-MM-DDTHHMMSS`) wins, with modification time and path breaking ties. Files
without a run timestamp use modification time as their priority. The complete
winning row is retained, including missing values; older valid values do not
replace missing values from a newer run. Year folder names do not limit data
coverage. Keep only the intended full-output runs below each configured folder.

### Reprocess all seven sites after rerunning EddyPro

```r
source("run_pipeline.R")
run_pipeline(base_dir = data_root, start = NULL, end = NULL,
             reprocess = TRUE, dry_run = TRUE)
# After reviewing the plan:
run_pipeline(base_dir = data_root, start = NULL, end = NULL,
             reprocess = TRUE, dry_run = FALSE)
```

The default targets are ECDP, EDVG, EDVP, ERVA, ERVP, ECSM, and ESTW.
With `reprocess = TRUE`, existing timestamps are processed again. A NULL start
uses each stage's earliest available input; a NULL end stops at each site's
latest EddyPro observation. All L1 outputs are prepared before L2 begins.
An L2 reference site with no observations overlapping the calculation window
is omitted for that window, so a later deployment does not block older data.
Missing reference output files still require an upstream L1 run. Existing
output/archive behavior is unchanged; rows outside the reprocessed window are
retained. This does not delete historical output files to start over.

## First run: preview, then process

Loading the runner does not start processing:

```r
source("run_pipeline.R")
```

Preview all seven sites without writing output:

```r
run_pipeline(base_dir = data_root, dry_run = TRUE)
```

Review the pending and calculation periods printed in the console, then run:

```r
run_pipeline(base_dir = data_root)
```

The default order is Level 1 for all selected sites, then Level 2 for all sites,
then Level 3 EC, Level 3 LI-710, and Level 4. Sites run sequentially within each
stage. This order ensures that Level 2 uses current multi-site Level 1 data.

`source("run_all.R")` starts the all-site run immediately. Edit `data_root` near
the top of that file first. New users should preview through `run_pipeline()`
before using this shortcut.

## Select a period

`start` and `end` are inclusive half-hour interval-end timestamps in local
standard time (without DST). R uses a UTC clock label without shifting these
local clock readings; it does not mean the observations occurred in UTC. Use
September 1 at midnight to include August's final half hour.

```r
run_pipeline(
  base_dir = data_root,
  start = "2026-08-01 00:30:00",
  end = "2026-09-01 00:00:00",
  dry_run = TRUE
)
```

Remove `dry_run = TRUE` after checking the preview. Processed rows are skipped.
With `end = NULL` (the default), the latest EddyPro observation timestamp defines
the common processing end **for each site**, including runs selecting only some
stages. This uses the observation clock, not the export filename date or the last
nonmissing flux. EddyPro input is therefore required for automatic-end runs.
Sites may have different ends from one another.

Logger and LI-710 time grids extend to that end. Missing logger measurements and
LI-710 observations remain NA; absent LI-710 samples also have NA QC flags. An
interval absent from the LI-710 logger is treated as an absent observation. Level
4 retains its existing gap-filling rules, so filled estimates can still be
available where observations are missing. An explicit `end` retains the previous
source-coverage behavior instead of applying this automatic EddyPro limit.

All stages updated in an automatic-end run share one site-specific figure folder
covering the earliest pending timestamp across selected stages through the
EddyPro end. Individual figures still plot only that stage's pending data.

Existing cumulative output rows are preserved, even if an older output extends
beyond the current EddyPro end; this change does not truncate historical files.
If source values are later corrected at timestamps already saved (including NA
rows), update those rows with `reprocess = TRUE`; optionally bound the affected period.

## Run one site

```r
run_pipeline(
  sites = "ECSM",
  base_dir = data_root,
  start = "2026-08-01 00:30:00",
  end = "2026-09-01 00:00:00",
  dry_run = TRUE
)
```

Only ECSM raw data and outputs are processed. Level 2 reads the other six sites'
existing Level 1 outputs as meteorological predictors; it does not read their
raw data or update them. Run all sites first when the reference Level 1 products
need updating.

The reference pool can be restricted, although this changes gap-filling inputs:

```r
run_pipeline(
  sites = "ECSM",
  reference_sites = c("ECSM", "ECDP"),
  base_dir = data_root
)
```

`run_example.R` includes editable all-site and single-site examples. Only the
first uncommented `run_pipeline()` call runs.

## Level 4 MDS context

REddyProc requires at least 90 days of regular half-hour data. The pipeline uses
100 days by default. If a new or selected period is shorter, Level 4 adds enough
preceding observations to make a 100-day calculation window. A 10-day new period
therefore uses 90 preceding days.

```r
run_pipeline(base_dir = data_root, min_mds_days = 100)
```

Only pending/requested rows are merged into the cumulative output. Context-period calculations do not replace
historical results. Saved Level 1, Level 2 and Level 3 EC products must cover the context; LI710 is optional. If they do not, the error reports the missing period;
process the listed upstream stages for that earlier period and rerun Level 4.

`context_days` may add even more preceding context. Level 4 uses whichever start
is earlier: explicit context or automatic MDS context.

## LI710 current files and backups

For each selected site, the runner discovers its configured LI710 file and
matching `LI710.dat`, `LI710.*.backup`, `LI710.dat.backup`, and `LI710_*.backup`
files in the configured source directories. It reads their observation clocks,
not dates in filenames. Backups are ordered by modification time; the active
`.dat` file wins exact timestamp overlaps. Duplicate samples are removed before
half-hour averaging. Missing optional `data_qc` uses the legacy default of zero.

A missing current `.dat` is allowed when a matching backup exists. With no LI710
source coverage, L3 LI710 is skipped. Level 4 requires L1, L2 and L3 EC only:
available L3 LI710 is left-joined, and absent LI710 channels remain NA. LI710 MDS
is skipped for channels with no usable observations; EC processing continues.

If backup data arrive for timestamps already saved as missing, use
`reprocess = TRUE` for L3 LI710 and L4, optionally restricting start/end, to update those rows.
Adding a source file alone does not overwrite an existing processed timestamp.

## Diagnostic figures

Each stage writes its named figures into a directory for the pending period:

```text
<site output>/figures/period_20260716T1700_to_20260801T0000/
```

Stages with the same pending bounds share the directory. These figures exclude
MDS training context. Reprocessing the same period replaces matching filenames
in that directory. No automatic plots PDF is generated.

Level 4 also renders the same six diagnostic figures from the merged cumulative
output, directly under `figures/`, with a `_full_period.jpeg` suffix. Only these
full-period files are overwritten on subsequent L4 updates. Historical QAQC and
MDS are not rerun to draw cumulative figures. Existing figures, `figure_archive`
and legacy `runs` directories are left untouched for manual organization.

Daily energy-balance scatter and time-series panels (EC and LI-710, corrected
and uncorrected) use only days with 48 distinct half-hourly timestamps and finite
values for all four plotted components: NETRAD_PI_F, G_PI_F and the corresponding
LE/H pair. Dates use interval starts (TIMESTAMP minus 30 minutes). Incomplete
days remain missing on the plotting calendar rather than being connected by a
line. Existing half-hourly and mean-diurnal energy panels are unchanged.

Water-balance diagnostics are saved separately as
`<site>_Level_4_Figure_2a_water_balance_EC.jpeg` and
`<site>_Level_4_Figure_2b_water_balance_LI710.jpeg`. Each uses the first through
last timestamp where that instrument has both filled and corrected ET, within
the period supplied to the figure function. The two full-period files therefore
have independent instrument coverage. An absent instrument produces a labeled
placeholder. Old combined water-balance figures are not deleted.

Daily ET requires complete paired uncorrected/corrected days. Cumulative ET and
precipitation use the same days with all 48 ET, corrected ET and precipitation
values present. They resume after missing days as **complete-day subtotals**, not
estimates of flux during the missing intervals; the figure reports included and
total day counts. Missing dates break ribbons/lines. Ribbon bounds use finite
paired minima/maxima, including when corrected ET is below uncorrected ET.
SSM retains the 0.08 m soil-layer assumption and is referenced to the first
available daily soil-moisture value in that instrument's period.

## Outputs and reprocessing

Each stage keeps one current cumulative CSV named
`<site>_<stage-pattern>_<YYYY-MM-DD>.csv`. The date is the execution date.
Level 4 retains **all previous CSV versions and their period audits** in
`archive/`. Same-day archive collisions receive a unique suffix. Levels 1-3 keep
only their immediately preceding cumulative version, as before.
A companion `*.csv.period.csv` records coverage, calculation context, and
pending-period flux/QC/gap-filling statistics with explicit denominators.

See [Output storage and quality diagnostics](OUTPUT_STORAGE.md) for details.

A failed stage does not publish a completed CSV. A `.pipeline-lock` directory
prevents simultaneous writes for one site. After a forcibly terminated session,
remove a stale lock only after confirming that no other run is active.

Existing timestamps count as processed, even if scientific values are `NA`.
Corrected source data require `reprocess = TRUE`. Bounds are optional:

- `start = NULL`: begin at each selected stage's earliest available interval.
- `end = NULL`: stop at each site's latest EddyPro observation, as in a normal run.
- Both omitted: reprocess all available intervals through that EddyPro end.

Existing rows in the selected interval are replaced; saved history outside that
interval is retained. Stage input/deployment limits and the Level 4 minimum MDS
context requirement still apply. Use `dry_run = TRUE` to inspect the resolved range.

```r
# Reprocess from the beginning through each site's EddyPro end.
run_pipeline(base_dir = data_root, reprocess = TRUE)

# Reprocess from a specified start through each site's EddyPro end.
run_pipeline(base_dir = data_root, start = "2026-08-01 00:30:00", reprocess = TRUE)
```

Explicit bounds can restrict the correction to a particular period:

```r
run_pipeline(
  sites = "ECSM",
  base_dir = data_root,
  start = "2026-08-01 00:30:00",
  end = "2026-09-01 00:00:00",
  reprocess = TRUE
)
```

Keep the full stage chain for upstream corrections so downstream outputs are
refreshed over the same period.

## Short Level 2 processing windows

`min_mds_days` applies to Level 4 only. In Level 2, the optional Leuning soil
heat-flux estimate uses its original 30-day smoothing when at least 30 daily
estimates are available. Shorter windows retain interpolated daily estimates.
When no damping depth can be estimated, Leuning results remain NA without
stopping processing. The final calorimetric `G_PI` calculation is unchanged.

## Troubleshooting

- L4 keeps the actual half-hourly start, including a mid-day deployment. It does
  not align to 00:30 or insert leading NA rows. Any additional MDS context must
  already exist upstream. Missing interior upstream rows and genuinely insufficient
  history still stop processing. L4 planning errors identify the site.
- `min_mds_days = 100` is a requested minimum, not fabricated pre-deployment data.
  A site with 90-99 days of common upstream coverage can be run separately with
  `min_mds_days = 90`; sites below 90 days require additional observations.

- `base_dir is not an accessible data root`: confirm the path and mount the
  network drive in Windows first.
- `Missing L1/L2/L3 output`: run the full pipeline or the missing upstream stages.
- `Insufficient Level 4 context`: process the earlier interval identified by the
  error through Level 1, Level 2, Level 3 EC, and Level 3 LI-710.
- Package installation failure: restart RStudio, verify internet/CRAN/GitHub
  access, and install matching Rtools if R requests build tools.

## LE MDS predictions, quality flags, and VBR-Ts

L4 saves all available MDS FillAll predictions, without restricting them to
quality code 1. This applies to EC, LI-710, VBR-Ts and residual-energy LE only.
Failed predictions remain missing. FillAll predictions can exist at measured
timestamps; they are not the same as observed-plus-gap-filled series.

| MDS prediction | Prediction quality | Observed-plus-filled series | Series quality |
| --- | --- | --- | --- |
| LE_MDS | LE_MDS_QC | LE_PI_F | LE_PI_F_QC |
| LE_710_MDS | LE_710_MDS_QC | LE_710_PI_F | LE_710_PI_F_QC |
| LE_VBR_Ts_MDS | LE_VBR_Ts_MDS_QC | LE_VBR_Ts_F | LE_VBR_Ts_F_QC |
| LE_RE_MDS | LE_RE_MDS_QC | LE_RE_F | LE_RE_F_QC |

Prediction flags preserve REddyProc `*_fall_qc`; series flags preserve `*_fqc`.
Codes 1, 2 and 3 denote good, medium and poor MDS quality; 0 in the series flag
denotes an original input observation. Missing flags remain missing. These flags
describe MDS quality, not the quality of an energy-balance correction. No H or
FC MDS quality columns are added. Existing filled fluxes, EBR corrections and
retained ET products are unchanged. These are L4 columns; L5 mapping is unchanged.

Only `VBR_Ts = 1004.834 * sigma_Ts / (Lv * sigma_q)` is calculated, with existing
LE_VBR_Ts and ET_VBR_Ts names retained. The alternative with the temperature
subtraction, and its VBR/LE_VBR/ET_VBR output family, are removed. Diagnostic
VBR panels now use VBR-Ts.

When merging old L4 history, retired columns are dropped in memory and missing
quality columns are padded with NA. Saved files and archives are not edited on
read. Previously discarded predictions cannot be recovered without reprocessing;
old quality flags are not inferred. To regenerate the entire saved period:

```r
run_pipeline(base_dir = data_root, stages = "L4", reprocess = TRUE,
             start = NULL, end = NULL)
```

This uses existing upstream products and still requires their full MDS context.
Incremental runs retain historical prediction values and unknown flags outside
the newly processed interval. L4 publication continues to archive previous files.

## Tests

Tests use temporary synthetic data and do not write to the network data root.
From the project directory:

```r
system2(file.path(R.home("bin"), "Rscript"), "tests/test_vbr_mds.R")
system2(file.path(R.home("bin"), "Rscript"), "tests/test_layout.R")
system2(file.path(R.home("bin"), "Rscript"), "tests/test_pipeline.R")
system2(file.path(R.home("bin"), "Rscript"), "tests/test_publication.R")
system2(file.path(R.home("bin"), "Rscript"), "tests/test_figure_scope.R")
system2(file.path(R.home("bin"), "Rscript"), "tests/test_level2_damping.R")
system2(file.path(R.home("bin"), "Rscript"), "tests/test_li710_despike.R")
system2(file.path(R.home("bin"), "Rscript"), "tests/smoke_li710.R")
system2(file.path(R.home("bin"), "Rscript"), "tests/test_li710_sources.R")
system2(file.path(R.home("bin"), "Rscript"), "tests/test_mds_context.R")
system2(file.path(R.home("bin"), "Rscript"), "tests/smoke_level3.R")
system2(file.path(R.home("bin"), "Rscript"), "tests/smoke_level4.R")
```

The layout test checks paths and metadata. The pipeline test checks incremental,
all-site, and single-site operation, merging, locks, and recovery. The context
test checks 100-day window logic. Smoke tests execute the real Level 3 EC and
Level 4 scripts using synthetic data, including MDS and figures.

See `VALIDATION.md` for the validation record. Real-data scientific review is
still required before routine production use.
