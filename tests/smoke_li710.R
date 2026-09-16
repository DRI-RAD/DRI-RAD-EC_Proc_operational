# Actual LI710 ingestion -> QC -> figures -> publication with no usable flux.
source("run_pipeline.R")
root <- tempfile("li710-smoke-"); dir.create(root)
t <- seq(pipeline_time("2026-07-16 17:00"), by = 1800, length.out = 735)
units_for <- function(x) as.data.frame(setNames(lapply(names(x), function(n) c("unit", "avg")), names(x)))
out <- file.path(root, "output")
met <- data.frame(TIMESTAMP = t, PotRad = 500*pmax(0, sin(seq_along(t)/48*2*pi)))
pipeline_write_table(met, units_for(met), file.path(out, "ECDP_Level_1_QAQCed_logger_data_fixture.csv"), "fixture")
# Numeric missing-value sentinels exercise the real ingestion path.
raw <- data.frame(TIMESTAMP = t, LE_710 = -9999, H_710 = -9999, diag = 0, flow = 200, tilt = 0, data_qc = 0)
raw_file <- file.path(root, "Site_LI710.dat")
pipeline_write_table(raw, units_for(raw), raw_file, "fixture")
old_file <- file.path(root, "old.dat")
old <- raw[1:10, ]; old$TIMESTAMP <- old$TIMESTAMP - 100*86400
pipeline_write_table(old, units_for(old), old_file, "fixture")
cfg <- data.frame(site = "ECDP", dir_output = out, dir_LI710 = raw_file,
                  dir_LI710_old = old_file, dir_eddypro = NA_character_)
config <- file.path(root, "config.csv"); readr::write_csv(cfg, config)
if (Sys.getenv("EC_LI710_BACKUP_ONLY") == "1") {
  stopifnot(file.rename(raw_file, file.path(root, "Site_LI710.202608.backup")))
}
args <- list(sites = "ECDP", stages = "L3_LI710", base_dir = root, config_file = config)
result <- do.call(run_pipeline, args)
data <- pipeline_read_table(result[[1]])$data
stopifnot(nrow(data) == 735, all(is.na(data$LE_710)), all(is.na(data$H_710)),
          all(is.na(data$LE_710_QC_despike)), all(is.na(data$H_710_QC_despike)))
figures <- list.files(file.path(out, "figures"), pattern = "Level_3_.*pdf$", full.names = TRUE, recursive = TRUE)
stopifnot(length(figures) == 2, all(file.info(figures)$size > 1000))
stopifnot(do.call(run_pipeline, args)[[1]] == "skipped")
unlink(root, recursive = TRUE)
cat("PASS: actual LI710 all-missing flux, empty legacy window, 735 output rows, two status PDFs, and skipped rerun.\n")
