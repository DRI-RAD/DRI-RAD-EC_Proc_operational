# Preview the default six-site run for a complete month from the project root.
source("run_pipeline.R")
run_pipeline(
  start = "2026-08-01 00:30:00",
  end = "2026-09-01 00:00:00",
  dry_run = TRUE
)
