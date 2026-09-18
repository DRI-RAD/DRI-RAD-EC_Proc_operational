# Preview examples. Only the first example runs; uncomment one alternative.
source("run_pipeline.R")
data_root <- "Z:/NWI/Task_3_Water_Use_ET_and_Meteoroligical_Monitoring/networks/eddy_stations"
# Another computer can use, for example, data_root <- "D:/eddy_stations".

run_pipeline(
  sites = "EDVP",
  base_dir = data_root,
  start = "2026-06-20 00:30:00",
  end = "2026-08-01 00:00:00",
  min_mds_days = 100,
  figure_period = "new",  # Figures show only the pending period; MDS context is excluded.
  dry_run = FALSE,
  reprocess = TRUE
)

# All six sites. Every site's L1 finishes before L2 uses cross-site references.
run_pipeline(
  base_dir = data_root,
  start = "2026-08-01 00:00:00",
  end = NULL,
  min_mds_days = 100,
  dry_run = FALSE,
  reprocess = TRUE
)

# Single target site. Only L2 reads other sites' EXISTING Level 1 results.
# Other sites' raw logger/EddyPro files are not read and their outputs are not updated.
# run_pipeline(
#   sites = "ECSM",
#   base_dir = data_root,
#   start = "2026-08-01 00:30:00",
#   end = "2026-09-01 00:00:00",
#   reference_sites = c("ECDP", "EDVG", "EDVP", "ERVA", "ERVP", "ECSM"),
#   min_mds_days = 100,
#   figure_period = "new",  # Figures show only the pending period; MDS context is excluded.
#   dry_run = TRUE
# )

# If the necessary historical QAQC outputs do not yet exist, prepare them first.
# The error message reports the missing period; change these example bounds.
# run_pipeline(sites = "ECSM", base_dir = data_root,
#              start = "2026-05-24 00:30:00", end = "2026-09-01 00:00:00",
#              stages = c("L1", "L2", "L3_EC", "L3_LI710"))
# Then rerun the desired NEW period above. L4 publishes only its pending rows.
