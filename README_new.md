# Incremental EC processing

This operational layout stores stages under `R/` and metadata under `meta_files/`. Use `run_pipeline.R` and the
`Level_*_new.R` scripts together. The changes concern input discovery, site and
period selection, and output publication. Scientific QAQC rules, thresholds,
regressions, and gap-filling calculations remain in the copied stage scripts.

## Quick start

Open the R project and run:

```r
source("run_pipeline.R")

# Preview August processing without writing output.
run_pipeline(
  start = "2026-08-01 00:30:00",
  end = "2026-09-01 00:00:00",
  dry_run = TRUE
)

# Run the same period. Existing output rows are skipped automatically.
run_pipeline(
  start = "2026-08-01 00:30:00",
  end = "2026-09-01 00:00:00"
)

# Later runs discover new observations without editing each Level script.
run_pipeline(sites = c("ECDP", "EDVG", "EDVP", "ERVA", "ERVP", "ECSM"))
```

`run_pipeline()` defaults to all six main-network sites. `run_all.R` starts processing directly; `run_example.R` starts in preview mode. Do not
source a new Level script directly; it requires the runner's period and staging
configuration. A dry run projects downstream coverage from planned upstream
rows; it does not execute or validate the scientific calculations.

## Site paths

`meta_files/file_directions_new.csv` retains the original logger, sensor, and output paths
and adds `dir_eddypro`. The legacy `dir_EC` Excel path is retained for reference
but is not read by the new EC stage.

| Site | Folder under EddyPro_Output_Files |
| --- | --- |
| ECDP | Carson_Desert_EC |
| ECSM | Columbus |
| EDVP | Desert_Valley |
| EDVG | DValfalfa |
| ERVA | Rvalfalfa |
| ERVP | RVphreat |

For example, the ECSM reader discovers:

```text
<base_dir>/EddyPro_Output_Files/Columbus/Aug_2026/output/
    eddypro_Columbus_full_output_2026-08-11T162337_adv.csv
```

`base_dir` defaults to the existing Z: network location. Override it in
`run_pipeline(base_dir = "D:/my_copy", ...)` for another storage root. A leading
single slash in configuration means relative to `base_dir`, matching the
original CSV convention. Drive-letter and UNC paths are accepted as absolute.

The reader detects the CSV header and units row, uses field names instead of
the former 116-column Excel slice, converts `-9999` to missing values, and builds
`TIMESTAMP` from `date` plus `time`. Optional random-error and storage fields
may be absent and become NA. Missing required scientific fields cause a clear
error. Files must be under an `output` folder and match
`eddypro_*full_output*.csv`. This excludes biomet and other unrelated exports.

Overlapping rows use the latest run timestamp in the filename, then file
modification time and path as tie-breakers. Month folder names are not treated
as observation limits because an export may span multiple months. Source CSVs
are read to determine their actual coverage; only the selected period is passed
to scientific processing. There is no persistent raw-file cache in this version.

## Period and incremental behavior

Each stage inspects TIMESTAMP rows from its existing four-header-row output
CSVs, including legacy snapshots. Save dates in filenames are not processing
boundaries. Overlapping historical snapshots are merged by modification time,
with newer rows taking precedence. Schema mismatches stop the run for explicit
migration rather than silently discarding columns.

The runner identifies timestamps available upstream but absent from that stage's
output. This supports both appending and missing interior rows. It calculates
the span from the first pending row through the last pending row, then publishes
only the pending rows. Already-processed rows inside that span can participate
in calculations but retain their previous saved values.

Logger Level 1 may advance ahead of EddyPro. EC Level 3 waits for the intersection
of L1/L2 coverage within the actual EddyPro time bounds. Missing observations
inside those bounds remain available to the original EasyFlux fallback. Level 4
waits for the common coverage of both logger stages, EC QAQC, and LI710 QAQC.
An absent future EddyPro month therefore does not get marked as completed EC.

Existing rows count as processed even if their flux is NA. Later corrections to
already-processed source observations, including late EddyPro values replacing
an earlier EasyFlux fallback, require explicit reprocessing:

```r
run_pipeline(
  start = "2026-08-01 00:30:00",
  end = "2026-09-01 00:00:00",
  reprocess = TRUE
)
```

Both bounds are required for reprocessing. Keep the default full stage chain
when correcting upstream data so downstream results are refreshed too.

Times use the original scripts' fixed UTC clock label, without DST conversion
or an automatic logger/EddyPro offset. Confirm that both sources use the same
clock convention. Bounds are inclusive, with half-hour interval-end timestamps.
Use a final midnight to include a complete day; a date without a time is not
accepted. The preserved REddyProc calculations expect suitable regular time
series and adequate training data; a tiny or incomplete-day interval can fail.

## Calculation context and references

`context_days = 0` is the default: process the pending span only. Setting, for
example, `context_days = 30` includes preceding observations in calculations but
still writes only pending rows. The immediately preceding cumulative rain-gauge
sample is retained to calculate the first precipitation increment.

Changing the input period can change statistical results even though the code
for despiking, precipitation corrections, regressions, and MDS remains the same.
Thirty days is an example context length, not a scientifically validated default.
Use a sufficiently long interval/context for the site and review the diagnostics;
this implementation does not promise numerical equivalence to a full-history run.

Level 2 retains all six main-network sites as its default reference pool. Running
one target site reads the other sites' existing Level 1 outputs for the same
calculation window. It only writes the selected target's Level 2 result. Run all
six sites to update every reference site's L1 first. If only a subset is available,
select it explicitly:

```r
run_pipeline(sites = "ECSM", reference_sites = c("ECSM", "ECDP"))
```

Reducing the reference pool affects cross-site gap filling. A missing reference
output or no reference observations in the calculation window stops processing.

## Outputs and recovery

The runner stages each site's inputs, calculations, and figures in a temporary
directory. Once the script succeeds and its final table includes every pending
timestamp, it merges the new rows with history and publishes a uniquely named
cumulative CSV. Existing files are retained. The new CSV keeps the same four
metadata/header rows and column names, with proper CSV quoting.

The committed final CSV is the authoritative checkpoint. Its accompanying
`*.period.csv` records the stage, coverage, new rows, total rows, context start,
and reprocessing mode. Diagnostic figures and intermediate exports are retained
under `runs/<version>/`. Audit CSVs are excluded from input discovery.

If a scientific stage fails, its final output is not published; completed earlier
stages remain reusable on the next run. A `.pipeline-lock` directory prevents
concurrent writes for the same site. Normal errors release it. After a killed R
process, remove a stale lock only after confirming that no other run is active.

Do not mix the old standalone readers with new versioned filenames: the old
readers assume a date-only suffix. Continue using the new runner for downstream
processing. Disk use grows because cumulative versions and diagnostics are kept;
archival/deletion is intentionally left to the operator.

## Files and scope

| File | Responsibility |
| --- | --- |
| run_pipeline.R | Site/stage ordering, periods, locks, staging and commits |
| run_example.R | Editable execution example |
| R/pipeline_io.R | EddyPro import, timestamp handling, history and CSV I/O |
| file_directions_new.csv | Updated site configuration |
| Level_1_met_QAQC_new.R | Existing logger QAQC with period-aware I/O |
| Level_2_met_PI_vars_new.R | Existing PI/gap filling with selected targets/references |
| Level_3_EC_data_QAQC_new.R | Monthly EddyPro CSV input and existing EC QAQC |
| Level_3_LI710_QAQC_new.R | Existing LI710 QAQC with period-aware I/O |
| Level_4_gapfilling_ET_calculation_new.R | Existing gap filling/ET with period-aware I/O |

The runner supports the six main-network sites. OSU, UofIdaho, S2, and RRR
variation scripts are not included here and are not routed through the main
network algorithm. The previous Level 5 file was an unfinished selection draft
with no export implementation; it is not included here and is not an executable pipeline
stage. The main Level 4 code assumes LI710 fields, so the runner checks for LI710
configuration instead of changing that algorithm.

## Validation

From the project root:

```sh
Rscript tests/test_pipeline.R
Rscript tests/smoke_level3.R
```

See `VALIDATION.md` for the tested scope and remaining validation on real data.
