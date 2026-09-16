# Run from the project root: Rscript tests/test_pipeline.R
# Temporary fixtures exercise I/O and orchestration without modifying NAS data.
source("run_pipeline.R")
test_root <- tempfile("pipeline-tests-")
dir.create(test_root)
check_error <- function(expr, text = NULL) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"))
  if (!is.null(text)) stopifnot(grepl(text, conditionMessage(error), fixed = TRUE))
}
clock <- pipeline_time(c("2026-08-01 00:30", "2026-08-01 01:00", "2026-08-01 01:30", "2026-08-01 02:00"))
units <- data.frame(TIMESTAMP = c("TS", NA), value = c("W m-2", "avg"))
data <- data.frame(TIMESTAMP = clock, value = c(1, NA, 3, 4))

# Timestamp precision, path handling, CSV quoting, NA values and metadata rows.
stopifnot(format(clock[1], "%H:%M:%S", tz = "UTC") == "00:30:00")
check_error(pipeline_time("invalid"), "Invalid timestamp")
stopifnot(pipeline_path("Z:/base", "/site/a.csv") == "Z:/base/site/a.csv",
          pipeline_path("Z:/base", "C:/site/a.csv") == "C:/site/a.csv")
file <- file.path(test_root, "table.csv")
pipeline_write_table(data, units, file, 'Description, with "quotes"')
roundtrip <- pipeline_read_table(file)
stopifnot(identical(as.numeric(roundtrip$data$TIMESTAMP), as.numeric(clock)),
          identical(roundtrip$data$value, data$value), nrow(roundtrip$units) == 2)
check_error(pipeline_write_table(data, units, file, "do not overwrite"), "Refusing to overwrite")

# Completed rows, interior holes, bounds, overlap context, and no-new-data runs.
history <- list(data = data[c(1, 3), ], units = units)
plan <- pipeline_plan(clock, history)
stopifnot(identical(as.numeric(plan$pending), as.numeric(clock[c(2, 4)])))
stopifnot(is.null(pipeline_plan(clock, list(data = data))))
plan <- pipeline_plan(clock, history, start = clock[2], end = clock[3], context_days = 2)
stopifnot(identical(as.numeric(plan$pending), as.numeric(clock[2])),
          plan$read_start == clock[2] - 2 * 86400)
stopifnot(nrow(pipeline_slice(data, list(read_start = clock[3], end = clock[4]), predecessor = TRUE)) == 3)

# EddyPro category/header/units rows, reordered fields, optional fields, missing
# values, overlapping month folders, and latest-run precedence.
ep_fields <- c("date", "time", "LE", "qc_LE", "H", "qc_H", "co2_flux", "qc_co2_flux", "Tau", "qc_Tau",
  "co2_mole_fraction", "co2_molar_density", "co2_var", "h2o_mole_fraction", "h2o_molar_density", "h2o_var",
  "ts_var", "u_var", "v_var", "w_var", "sonic_temperature", "wind_speed", "max_wind_speed", "wind_dir",
  "u*", "(z-d)/L", "L", "x_90%", "x_70%", "x_50%", "x_30%", "x_10%")
make_ep <- function(file, value, fields = ep_fields) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  frame <- as.data.frame(setNames(lapply(fields, function(x) rep("1", 2)), fields), check.names = FALSE)
  frame$date <- c("2026-08-01", "2026-08-01")
  frame$time <- c("00:30", "01:00")
  frame$LE <- as.character(value)
  writeLines("EddyPro category row", file)
  write.table(matrix(fields, nrow = 1), file, append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
  write.table(matrix(rep("unit", length(fields)), nrow = 1), file, append = TRUE,
              sep = ",", row.names = FALSE, col.names = FALSE)
  readr::write_csv(frame, file, append = TRUE, col_names = FALSE)
}
ep_dir <- file.path(test_root, "Columbus")
ep1 <- file.path(ep_dir, "Aug_2026", "output", "eddypro_Columbus_full_output_2026-08-11T162337_adv.csv")
ep2 <- file.path(ep_dir, "Sep_2026", "output", "eddypro_Columbus_full_output_2026-09-01T120000_adv.csv")
make_ep(ep1, c(10, 20))
make_ep(ep2, c(30, -9999), rev(ep_fields))
ep <- pipeline_eddypro(ep_dir)
stopifnot(nrow(ep$data) == 2, ep$data$LE[1] == 30, is.na(ep$data$LE[2]),
          all(is.na(ep$data$LE_strg)), identical(as.numeric(ep$data$TIMESTAMP), as.numeric(clock[1:2])))
stopifnot(nrow(pipeline_eddypro(ep_dir, list(read_start = clock[2], end = clock[2]))$data) == 1)
bad <- file.path(test_root, "bad.csv")
make_ep(bad, c(1, 2), setdiff(ep_fields, "H"))
check_error(pipeline_eddypro_file(bad), "Missing EddyPro fields")

# Merge partial legacy snapshots using actual timestamps rather than save dates.
out <- file.path(test_root, "history")
dir.create(out)
pipeline_write_table(data[1:2, ], units, file.path(out, "ECSM_Level_1_QAQCed_logger_data_2026-01-01.csv"), "old")
pipeline_write_table(data[3, , drop = FALSE], units, file.path(out, "ECSM_Level_1_QAQCed_logger_data_2026-01-02.csv"), "partial")
hist <- pipeline_history(out, pipeline_patterns[["L1"]])
stopifnot(nrow(hist$data) == 3)
staging <- file.path(test_root, "staging")
dir.create(staging)
options(ec.pipeline = list(window = pipeline_plan(clock, hist)))
pipeline_write(data, units, file.path(staging, "ECSM_Level_1_QAQCed_logger_data_2026-01-03.csv"))
target <- pipeline_publish(staging, out, pipeline_patterns[["L1"]], getOption("ec.pipeline")$window, FALSE)
hist <- pipeline_history(out, pipeline_patterns[["L1"]])
stopifnot(nrow(hist$data) == 4, identical(hist$data$value, data$value),
          is.null(pipeline_plan(clock, hist)), file.exists(paste0(target, ".period.csv")))

# Runner integration uses deterministic stand-ins for expensive science.
# These exercise period propagation, publication, checkpoints, locks and retry.
real_figures <- pipeline_level4_figures
pipeline_level4_figures <- function(data, site, folder) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  writeLines(as.character(nrow(data)), file.path(folder, paste0(site, "_Level_4_fixture.pdf")))
}
real_mds_window <- pipeline_mds_window
# Tiny orchestration fixtures intentionally bypass only the scientific length
# guard. test_mds_context.R and smoke_level4.R exercise the real guard separately.
pipeline_mds_window <- function(window, available, min_days) window
real_root <- .pipeline_root
fixture <- file.path(test_root, "runner")
dir.create(fixture)
dir.create(file.path(fixture, "R"))
.pipeline_root <- fixture
raw <- file.path(fixture, "raw.dat")
pipeline_write_table(data, units, raw, "TOA5 fixture")
raw_li <- file.path(fixture, "Site_LI710.dat")
li <- data.frame(TIMESTAMP = clock, LE_710 = 1, H_710 = 1, diag = 0, flow = 200, tilt = 0)
u_li <- as.data.frame(setNames(lapply(names(li), function(n) c("unit", "avg")), names(li)))
pipeline_write_table(li, u_li, raw_li, "fixture")
cfg <- data.frame(site = "ECSM", dir_met = raw, dir_output = file.path(fixture, "output"),
                  dir_eddypro = ep_dir, dir_LI710 = raw_li, dir_LI710_old = NA_character_)
readr::write_csv(cfg, file.path(fixture, "config.csv"))
for (stage in names(pipeline_stages)) {
  script <- c('opt <- getOption("ec.pipeline")',
              'x <- data.frame(TIMESTAMP = opt$window$pending, value = 5)',
              'u <- data.frame(TIMESTAMP = c("TS", NA), value = c("W m-2", "avg"))',
              paste0('pipeline_write(x, u, file.path(opt$dirs$dir_output[opt$dirs$site == "ECSM"], "ECSM_',
                     pipeline_patterns[[stage]], '_fixture.csv"))'))
  writeLines(script, file.path(fixture, pipeline_stages[[stage]]))
}
args <- list(sites = "ECSM", reference_sites = "ECSM", base_dir = fixture, config_file = file.path(fixture, "config.csv"))
preview <- do.call(run_pipeline, c(args, list(dry_run = TRUE)))
stopifnot(length(preview) == 5, !dir.exists(cfg$dir_output))
first <- do.call(run_pipeline, args)
stopifnot(length(first) == 5, all(file.exists(unlist(first))))
second <- do.call(run_pipeline, args)
stopifnot(all(unlist(second) == "skipped"))
stopifnot(nrow(pipeline_history(cfg$dir_output, pipeline_patterns[["L1"]])$data) == 4,
          nrow(pipeline_history(cfg$dir_output, pipeline_patterns[["L4"]])$data) == 2)
before <- pipeline_history(cfg$dir_output, pipeline_patterns[["L1"]])$data
writeLines('stop("Injected stage failure")', file.path(fixture, pipeline_stages[["L1"]]))
check_error(do.call(run_pipeline, c(args, list(stages = "L1", start = clock[1], end = clock[2], reprocess = TRUE))),
            "Injected stage failure")
stopifnot(isTRUE(all.equal(pipeline_history(cfg$dir_output, pipeline_patterns[["L1"]])$data, before)),
          !dir.exists(file.path(cfg$dir_output, ".pipeline-lock")))
# Explicit replacement changes only the requested timestamps and resumes safely.
script <- c('opt <- getOption("ec.pipeline")',
            'x <- data.frame(TIMESTAMP = opt$window$pending, value = 99)',
            'u <- data.frame(TIMESTAMP = c("TS", NA), value = c("W m-2", "avg"))',
            'pipeline_write(x, u, file.path(opt$dirs$dir_output[opt$dirs$site == "ECSM"], "ECSM_Level_1_QAQCed_logger_data_fixture.csv"))')
writeLines(script, file.path(fixture, pipeline_stages[["L1"]]))
do.call(run_pipeline, c(args, list(stages = "L1", start = clock[1], end = clock[2], reprocess = TRUE)))
after <- pipeline_history(cfg$dir_output, pipeline_patterns[["L1"]])$data
stopifnot(identical(after$value, c(99, 99, 5, 5)))
dir.create(file.path(cfg$dir_output, ".pipeline-lock"))
check_error(do.call(run_pipeline, c(args, list(stages = "L1"))), "Site is locked")
unlink(file.path(cfg$dir_output, ".pipeline-lock"), recursive = TRUE)
# Verify the default six-site run and the stage-major reference dependency.
network <- c("ECDP", "EDVG", "EDVP", "ERVA", "ERVP", "ECSM")
network_cfg <- cfg[rep(1, length(network)), ]
network_cfg$site <- network
network_cfg$dir_output <- file.path(fixture, "network", network)
readr::write_csv(network_cfg, file.path(fixture, "network.csv"))
for (stage in names(pipeline_stages)) {
  script <- c('opt <- getOption("ec.pipeline")',
    'site <- opt$sites',
    'folder <- opt$dirs$dir_output[opt$dirs$site == site]',
    if (stage == "L2") c(
      'stopifnot(length(opt$reference_sites) == 6)',
      'for (ref in opt$reference_sites) {',
      '  ref_folder <- opt$dirs$dir_output[opt$dirs$site == ref]',
      '  stopifnot(length(pipeline_files(ref_folder, "Level_1_QAQCed_logger_data")) == 1)',
      '}'),
    'x <- data.frame(TIMESTAMP = opt$window$pending, value = 5)',
    'u <- data.frame(TIMESTAMP = c("TS", NA), value = c("W m-2", "avg"))',
    paste0('pipeline_write(x, u, file.path(folder, paste0(site, "_', pipeline_patterns[[stage]], '_fixture.csv")))'))
  writeLines(script, file.path(fixture, pipeline_stages[[stage]]))
}
network_result <- run_pipeline(base_dir = fixture, config_file = file.path(fixture, "network.csv"))
stopifnot(length(network_result) == 30, all(file.exists(unlist(network_result))))
stopifnot(identical(names(network_result)[1:6], paste(network, "L1", sep = ":")))
stopifnot(all(unlist(run_pipeline(base_dir = fixture, config_file = file.path(fixture, "network.csv"))) == "skipped"))
# A single target may load other sites only as L2 references. Invalid raw paths
# for unselected sites must never be accessed and their outputs must not change.
other_outputs <- unlist(lapply(network_cfg$dir_output[network_cfg$site != "ECSM"], list.files, full.names = TRUE))
before_other <- tools::md5sum(other_outputs[!dir.exists(other_outputs)])
network_cfg$dir_met[network_cfg$site != "ECSM"] <- "missing-reference-raw.dat"
network_cfg$dir_eddypro[network_cfg$site != "ECSM"] <- "missing-reference-eddypro"
readr::write_csv(network_cfg, file.path(fixture, "network.csv"))
single <- run_pipeline(sites = "ECSM", base_dir = fixture,
  config_file = file.path(fixture, "network.csv"), start = clock[1], end = clock[2], reprocess = TRUE)
stopifnot(length(single) == 5, all(startsWith(names(single), "ECSM:")),
          identical(before_other, tools::md5sum(names(before_other))))
# No current/backup LI710 files must skip L3 LI710 without blocking L4.
network_cfg$dir_LI710[network_cfg$site == "ECSM"] <- file.path(fixture, "no_li710", "Site_LI710.dat")
network_cfg$dir_met[network_cfg$site == "ECSM"] <- file.path(fixture, "no_li710", "met.dat")
readr::write_csv(network_cfg, file.path(fixture, "network.csv"))
optional <- run_pipeline(sites = "ECSM", stages = c("L3_LI710", "L4"),
  base_dir = fixture, config_file = file.path(fixture, "network.csv"),
  start = clock[1], end = clock[2], reprocess = TRUE)
stopifnot(optional[["ECSM:L3_LI710"]] == "skipped", file.exists(optional[["ECSM:L4"]]))
pipeline_level4_figures <- real_figures
pipeline_mds_window <- real_mds_window
.pipeline_root <- real_root
options(ec.pipeline = NULL)
for (file in unname(pipeline_stages)) parse(file.path(real_root, file))
unlink(test_root, recursive = TRUE)
cat("PASS: ingestion, periods, merges, dry run, incremental rerun, failure recovery, and R syntax.\n")
