# Jasoni: open DRI-RAD-EC_Proc_operational.Rproj first, then run this file.
# Only the main call below is active. Alternative examples are commented out.
source("R/run_pipeline.R")  # Also loads run_level5_ameriflux().

data_root <- "Y:/NWI/Task_3_Water_Use_ET_and_Meteoroligical_Monitoring/networks/eddy_stations"
all_sites <- c("ECDP", "EDVG", "EDVP", "ERVA", "ERVP", "ECSM", "ESTW")

# Default: NULL uses each site's dir_output in meta_files/file_directions_new.csv,
# relative to data_root (normally data_root/Post_Processed_Data/<site folder>).
output_root <- NULL
# To use another location, replace NULL above, for example:
# output_root <- "D:/EC_processed"
# A custom root stores/reads history under output_root/<site ID>.
# Non-target reference sites still use their configured default L1 locations.

# MAIN RUN: all sites, only new periods; L1 -> L2 -> L3 EC/LI710 -> L4 -> AmeriFlux.
run_pipeline(
  base_dir = data_root,
  sites = all_sites,
  output_dir = output_root,
  reprocess = FALSE,        # Skip timestamps already saved in each stage.
  start = NULL,            # Automatically select the pending period.
  end = NULL,              # End at each site's latest EddyPro observation.
  min_mds_days = 90,       # At least 90; use 100 for additional MDS context.
  dry_run = FALSE,         # TRUE previews the plan without writing products.
  export_ameriflux = TRUE,
  ameriflux_plots = TRUE
)

# ALTERNATIVES: comment out the main call before enabling an example.
# 1. Preview all sites without processing or saving.
# run_pipeline(base_dir = data_root, sites = all_sites,
#              output_dir = output_root, reprocess = FALSE, dry_run = TRUE)

# 2. Single site. L2 can read other sites' existing L1 for met gap filling.
# run_pipeline(base_dir = data_root, sites = "ECSM",
#              output_dir = output_root, reprocess = FALSE,
#              reference_sites = all_sites, min_mds_days = 90)

# 3. Restrict NEW processing to an inclusive local-standard half-hour period.
# Existing timestamps are still skipped because reprocess = FALSE.
# run_pipeline(base_dir = data_root, sites = all_sites,
#              output_dir = output_root, reprocess = FALSE,
#              start = "2026-11-01 00:30:00", end = "2026-11-20 00:00:00",
#              min_mds_days = 90)

# 4. Run only pending L4 from existing upstream products.
# Allowed stages: L1, L2, L3_EC, L3_LI710, L4 (there is no L5 stage).
# run_pipeline(base_dir = data_root, sites = all_sites, stages = "L4",
#              output_dir = output_root, reprocess = FALSE, min_mds_days = 90)

# 5. Regenerate ONLY AmeriFlux files/figures from existing complete L4.
# No QAQC/MDS is rerun. No reprocess argument is needed.
# processed_dir points to the custom L1-L4 root, or NULL for configured paths.
# Here output_dir is a DIFFERENT option: it controls only the AmeriFlux destination.
# NULL writes into each site's existing output folder under AmeriFlux/.
# run_level5_ameriflux(base_dir = data_root, sites = all_sites,
#                     processed_dir = output_root, output_dir = NULL,
#                     make_plots = TRUE, dry_run = FALSE)

# Further options for run_pipeline():
# stages = c("L1", "L2")  # Run selected stages only.
# context_days = 0        # Extra preceding calculation context; not extra output.
# figure_period = "new"   # Period figures cover pending rows (supported value).
# export_ameriflux = FALSE # Skip automatic export after L4.
# ameriflux_plots = FALSE  # Export CSV without AmeriFlux diagnostic plots.
# To replace previously saved results deliberately, change reprocess to TRUE.
# With TRUE, NULL bounds select all available history through the EddyPro end.
