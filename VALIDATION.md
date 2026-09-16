# Operational layout validation

## September 15: context and diagnostic updates

- test_mds_context.R checks exact 100-day padding for a 10-day pending interval,
  longer intervals, extra context, missing history and figure-period selection.
- test_pipeline.R additionally checks single-site execution with deliberately
  invalid raw paths for all other sites. Only L2 reads their saved reference
  outputs and those output files remain unchanged.
- smoke_level4.R executes the real Level 4 calculations and MDS with 4,800 input
  rows, including missing initial fluxes. It checks that only 480 requested rows
  are published, five diagnostic JPEGs are present, and a second run is skipped.
- All user-restored Level 4 figure blocks are retained. The MDS calculations use
  the complete context window, while a separate table controls figure coverage.

The installed REddyProc fCheckHHTimeSeries function explicitly requires at least
90 days. The new default is 100 days; adequate usable predictors still depend on
the observations. Historical context must exist in the saved upstream products.

The earlier layout-validation record follows.

- Restored the missing runner and execution examples.
- Confirmed all five stage paths resolve under R/.
- Loaded the actual manual-flag and soil metadata from meta_files/ while running
  from another working directory.
- Checked the default six-site selection and all 30 site/stage combinations with
  deterministic test calculations. Every L2 target received all six L1 references.
- Checked no-new-data reruns, partial history merging, explicit replacement,
  locking, and failure recovery using the existing synthetic regression tests.
- Executed the actual Level 3 EC script with 1,488 synthetic half-hour observations,
  including QAQC, diagnostic generation, output publication and a skipped rerun.

The test environment uses R 4.4.2. Package build-version and synthetic regression
warnings may be emitted. Full scientific Level 1 through Level 4 validation on
real observations is not claimed: access to the Z: data root was denied in the
execution environment. No operational network outputs were modified.

At the initial layout revision, changes to stage files were limited to metadata
paths. This later revision also adjusts Level 4 input-window selection and figure
handling, with shared I/O and runner updates as described above. Metadata files
and Level 1 through Level 3 scientific code remain unchanged.

## Output retention revision (2026-09-15)

- `test_publication.R`: same-day rotation, one previous cumulative version,
  repeated-period figure preservation, legacy run migration, injected publication
  failure with rollback, and explicit diagnostic denominators passed.
- `test_pipeline.R`: default all-site and single-site orchestration, incremental
  updates, reprocessing, locks, and stage-failure recovery passed.
- `test_layout.R` and `test_mds_context.R`: paths, metadata and 100-day context
  passed. Figures now always use pending timestamps.
- `smoke_level3.R`: actual Level 3 EC QAQC and diagnostic publication passed.
- `smoke_level4.R`: actual REddyProc MDS on a 100-day window, 480-row (10-day)
  publication, all five JPEGs, and raw/filling diagnostic count assertions passed.

QAQC and gap-filling equations/thresholds are unchanged. Stage edits restrict
plot inputs to pending timestamps, repair the Level 2 staging figure path, and
capture MDS quality before it is discarded. Tests use synthetic local data;
operational network files were not processed or rotated during this revision.
R 4.4.2 emitted locale/package-version and synthetic-data plotting warnings.
- `test_figure_scope.R`: actual Level 2 plotting function/call generated its PDF in the absolute staging folder with only pending rows.

## Level 2 short-window repair (2026-09-15)

The soil Leuning reconstruction previously applied a 30-day rolling mean even
when the incremental input contained fewer than 30 daily estimates. This could
make the damping depths and phase shifts entirely missing and fail in lead(n=NA).
Windows shorter than 30 days now retain the interpolated daily damping estimates;
longer windows keep the original 30-day smoothing. One finite estimate can fill
edges; no finite estimate leaves the reconstruction missing. Invalid phase shifts
produce NA instead of aborting. Calorimetric G_PI remains unchanged.

`tests/test_level2_damping.R` executes the real reconstruction loop with 1, 16,
29, 30, 31 and 60 days, all-missing sun peaks, a single valid daily peak and a
midnight-only interval. It verifies unchanged timestamps/row counts/G_PI and
parity with the former smoothing on a supported 60-day input. All cases passed
with synthetic data. Missing-data plot warnings are expected in these fixtures.

## LI710 empty-subset repair (2026-09-15)

The installed openeddy despikeLF implementation fails in desp_loop when a
physical-range-filtered day/night subset is empty (var_minus assignment), or
when its available second differences do not exceed nVals. The LI710 adapter
now guards those two cases in a private function environment without changing
the installed package. Unassessed flags remain NA; existing range/spike flags
are retained. Supported subsets use the original calculation and thresholds.
Other errors still propagate. A period with no good-quality observations emits
a diagnostic status PDF rather than aborting figure publication.

- test_li710_despike.R reproduced the reported var_minus error and passed empty,
  all-rejected, daytime-only, sparse, and ordinary-data cases. Ordinary results
  exactly matched openeddy; invalid timestamp validation remained active.
- smoke_li710.R ran the real LI710 stage on 735 all-missing flux rows, with an
  empty legacy-file slice. It published all timestamps, missing QC results,
  two diagnostic status PDFs, and skipped a no-new-data rerun.

Tests used local synthetic observations; the actual ECDP observations were not
inspected, so the specific reason its subset became empty is not established.

## EC empty-subset protection (2026-09-15)

Level 3 EC LE/H/FC now use the shared pipeline_despike adapter. The LI710 name
remains as a backward-compatible alias. QC thresholds and supported-subset
calculations are unchanged. EC diagnostics write a status PDF when no values
pass all QC tests, as LI710 does. This does not fill missing flux observations.

The actual EC stage smoke test passed normal data, all missing fluxes and a
1,289-interval continuous gap (about 27 days) within a 1,488-row period. Each
case published the correct row count and three QC PDFs, retained missing
observations/flags, and skipped the subsequent no-new-data run. Shared helper
parity/validation tests also passed. Tests used synthetic local observations.

To repeat the additional EC smoke cases from the project root:

```r
Sys.setenv(EC_SMOKE_MISSING = "long_gap")
source("tests/smoke_level3.R")
Sys.setenv(EC_SMOKE_MISSING = "all")
source("tests/smoke_level3.R")
Sys.unsetenv("EC_SMOKE_MISSING")
```

## LI710 discovery, optional input and figure layout (2026-09-16)

- Added pipeline_li710.R to discover the configured stream's current file and
  matching dot/underscore backup variants. Paths are normalized, exact sample
  duplicates are resolved before half-hour averaging, and current .dat wins
  overlaps. The same discovered observations drive planning and L3 ingestion.
- Level 4 planning/context now requires L1/L2/L3 EC only. Missing LI710 inputs
  receive a stable NA schema; channels with no finite observations skip MDS.
  Numeric typing is preserved for entirely missing channels and daily figures.
- Extracted the original five L4 figures into pipeline_level4_figures.R. Figures
  render from the pending table and from the merged cumulative table, without
  rerunning historical QAQC/MDS. Removed the old ERVP plot-only start-date clamp
  so cumulative figures include the full saved period.
- Removed automatic_plots generation. Each stage publishes its explicit figures
  into figures/period_<first>_to_<last>. Full L4 figures use a _full_period suffix
  directly in figures and are overwritten on updates. Existing figure files,
  archives and legacy runs are untouched.
- L4 archives preserve every previous CSV and audit, with collision suffixes for
  same-day reruns. L1-L3 retain their one-previous-version policy. Rendering is
  completed before rotating published outputs, with rollback for copy failures.

Validation used R 4.4.2 and synthetic local data only:
- test_li710_sources.R: discovery, backup-only/no-file cases, site isolation,
  normalized path deduplication, overlap precedence and period slicing passed.
- smoke_li710.R: real stage passed with active or backup-only source files and
  735 all-missing flux rows; output and two status PDFs were published.
- smoke_level4.R: real 100-day MDS calculation passed with normal LI710, absent
  LI710, and only 2 days of LI710 with an absent H channel. Each saved 480 pending
  rows plus five period/five cumulative figures. No automatic PDF was produced.
- test_publication.R: period directories, cumulative figure rows, existing-figure
  preservation, same-day L4 archive retention and injected-failure rollback passed.
- test_pipeline.R: all-site/single-site orchestration and missing-source LI710 skip
  with successful L4 continuation passed. Layout, MDS context and L2 figure tests
  passed; actual EC smoke test passed with period-directory publication.

Optional smoke variants:
```r
Sys.setenv(EC_LI710_BACKUP_ONLY = "1")
source("tests/smoke_li710.R")
Sys.unsetenv("EC_LI710_BACKUP_ONLY")
Sys.setenv(EC_NO_LI710 = "1")
source("tests/smoke_level4.R")
Sys.unsetenv("EC_NO_LI710")
Sys.setenv(EC_LI710_PARTIAL = "1")
source("tests/smoke_level4.R")
Sys.unsetenv("EC_LI710_PARTIAL")
```
The prior validation sections describe earlier layouts; this section supersedes
their figure-archive and mandatory-LI710 behavior. Operational NAS data were not
processed or reorganized by these tests. Package build and synthetic regression/
missing-data plotting warnings remain possible.
