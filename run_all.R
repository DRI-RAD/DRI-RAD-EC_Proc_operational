# Run all six main-network sites, completing each level before the next level.
# Source this file in RStudio, or execute Rscript /path/to/run_all.R.
.run_sources <- Filter(function(x) !is.null(x), lapply(sys.frames(), function(x) x$ofile))
.run_file <- if (length(.run_sources)) tail(.run_sources, 1)[[1]] else {
  arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(arg)) sub("^--file=", "", arg[1]) else "run_all.R"
}
source(file.path(dirname(normalizePath(.run_file, mustWork = TRUE)), "run_pipeline.R"))

# NULL bounds discover unprocessed observations automatically.
# For a selected period, use inclusive interval-end timestamps, for example:
# start = "2026-08-01 00:30:00", end = "2026-09-01 00:00:00".
# Set dry_run = TRUE to inspect the plan without processing or writing output.
run_pipeline(
  start = NULL,
  end = NULL,
  context_days = 0,
  dry_run = FALSE
)
