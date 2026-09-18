# DRI-RAD Eddy Covariance Processing Pipeline

This repository processes logger and EddyPro data for NICE NET (DRI) eddy
covariance sites: ECDP, EDVG, EDVP, ERVA, ERVP, and ECSM. It runs meteorological
QAQC, multi-site meteorological gap filling, EC and LI-710 QAQC, REddyProc MDS
gap filling, ET calculations, and diagnostic figures.

To prepare an AmeriFlux review file from existing Level 4 results, use
`source("run_level5_ameriflux.R")` followed by
`run_level5_ameriflux(base_dir = data_root)`. All six sites are selected by
default. This separate export step preserves the input files. See
[Level 5 instructions and variable mapping](AMERIFLUX_LEVEL5.md), including
the requested flux QC policy and the distinction between local `_PI` review
columns and final AmeriFlux submission labels.

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
<base_dir>/EddyPro_Output_Files/<site>/<month_year>/output/
    eddypro_<site>_full_output_<run_timestamp>_adv.csv
```

The reader detects field names and units, converts `-9999` to missing values,
creates `TIMESTAMP` from `date` and `time`, and uses the latest run when exports
overlap.

## First run: preview, then process

Loading the runner does not start processing:

```r
source("run_pipeline.R")
```

Preview all six sites without writing output:

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
rows), update those rows with an explicitly bounded `reprocess = TRUE` run.

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

Only ECSM raw data and outputs are processed. Level 2 reads the other five sites'
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

If backup data arrive for timestamps already saved as missing, use explicit
start/end bounds and `reprocess = TRUE` for L3 LI710 and L4 to update those rows.
Adding a source file alone does not overwrite an existing processed timestamp.

## Diagnostic figures

Each stage writes its named figures into a directory for the pending period:

```text
<site output>/figures/period_20260716T1700_to_20260801T0000/
```

Stages with the same pending bounds share the directory. These figures exclude
MDS training context. Reprocessing the same period replaces matching filenames
in that directory. No automatic plots PDF is generated.

Level 4 also renders the same five diagnostic figures from the merged cumulative
output, directly under `figures/`, with a `_full_period.jpeg` suffix. Only these
full-period files are overwritten on subsequent L4 updates. Historical QAQC and
MDS are not rerun to draw cumulative figures. Existing figures, `figure_archive`
and legacy `runs` directories are left untouched for manual organization.

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
Corrected source data require explicit bounds and `reprocess = TRUE`:

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

- `base_dir is not an accessible data root`: confirm the path and mount the
  network drive in Windows first.
- `Missing L1/L2/L3 output`: run the full pipeline or the missing upstream stages.
- `Insufficient Level 4 context`: process the earlier interval identified by the
  error through Level 1, Level 2, Level 3 EC, and Level 3 LI-710.
- Package installation failure: restart RStudio, verify internet/CRAN/GitHub
  access, and install matching Rtools if R requests build tools.

## VBR and VBR_Ts MDS comparison outputs

`LE_VBR_MDS` and `LE_VBR_Ts_MDS` retain their respective MDS estimates only
where that estimate's `*_fall_qc` equals 1. VBR without the `_Ts` suffix is
retained for testing only and was not used for the study analysis; VBR_Ts was used.

The earlier export block overwrote the VBR estimate with the Ts estimate and
left the Ts comparison column unmasked. This is now corrected. The separate
filled series `LE_VBR_Ts_F` and `ET_VBR_Ts_F` are unchanged by this export fix.
Existing CSVs are not rewritten automatically: use explicit period bounds and
`reprocess = TRUE` for Level 4 to regenerate already processed rows.

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
