# Operational EC processing

Open `DRI-RAD-EC_Proc_operational.Rproj`, then run:

```r
source("run_pipeline.R")
run_pipeline(dry_run = TRUE)  # Preview all six sites.
run_pipeline()                # Process all six sites.
```

Alternatively, source `run_all.R` to start processing directly. Opening the
project or sourcing `run_pipeline.R` alone does not start processing.

The default sites are ECDP, EDVG, EDVP, ERVA, ERVP, and ECSM. All selected sites
complete Level 1 before any site's Level 2 starts. Level 2 keeps all six sites
as references by default, including when an explicit target subset is selected.
Sites execute sequentially within each level. Already processed rows are skipped.

```r
# Optional period: inclusive half-hour interval-end timestamps.
run_pipeline(start = "2026-08-01 00:30:00", end = "2026-09-01 00:00:00")

# Optional target subset; other sites' existing Level 1 outputs remain references.
run_pipeline(sites = "ECSM")
```

## Directory layout

```text
run_pipeline.R        Runner and default settings
run_all.R             Direct execution for all six sites
run_example.R         Preview example
R/                    Stage scripts and shared I/O functions
meta_files/           file_directions_new.csv, manual_flag.csv, meta_G.csv
tests/                Synthetic regression and layout tests
docs/                 Existing scientific pipeline documentation
```

Edit `meta_files/file_directions_new.csv` for source and output paths. The default
`base_dir` remains the supplied Z: network location; pass another `base_dir` to
use a local data copy. The stage scripts load metadata from the project directory,
independently of the current working directory. Use `source()` with the absolute
path to `run_pipeline.R` when starting outside the project directory.

QAQC and gap-filling calculation blocks are unchanged from the supplied archive.
Only runner restoration, layout references, defaults, and tests were adjusted.
The input period still affects statistical training. `context_days = 0` retains
the existing behavior; choose additional context explicitly if needed. Diagnostic
figures are still saved per successful site/stage run under `runs/<version>/figures`.

The archive's `.Rproj.user` session cache is excluded from this distribution.
Its scientific code, metadata, and existing docs are retained. See `README_new.md`
for the detailed incremental-processing behavior.

## Validation

```sh
Rscript tests/test_layout.R
Rscript tests/test_pipeline.R
Rscript tests/smoke_level3.R
```

See `VALIDATION.md` for tested behavior and limitations.
