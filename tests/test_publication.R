# Run from the repository root. All rotation tests use temporary local files.
source("run_pipeline.R")
root <- tempfile("publication-tests-"); dir.create(root)
out <- file.path(root, "out"); dir.create(out)
stage <- file.path(root, "stage"); dir.create(file.path(stage, "figures"), recursive = TRUE)
t <- pipeline_time(c("2026-08-01 00:30", "2026-08-01 01:00", "2026-08-01 01:30"))
u <- data.frame(TIMESTAMP = c("TS", NA), value = c("unit", "avg"))
p <- pipeline_patterns[["L1"]]
publish <- function(index, value, reprocess = FALSE) {
  unlink(list.files(stage, pattern = "csv$", full.names = TRUE))
  x <- data.frame(TIMESTAMP = t[index], value = value)
  pipeline_write_table(x, u, file.path(stage, paste0("ECSM_", p, "_stage.csv")), "fixture")
  w <- list(pending = t[index], read_start = min(t[index]), end = max(t[index]))
  options(ec.pipeline = list(window = w))
  pipeline_publish(stage, out, p, w, reprocess)
}
writeLines("first figure", file.path(stage, "figures", "ECSM_Level_1_test.pdf"))
first <- publish(1, 10)
stopifnot(basename(first) == paste0("ECSM_", p, "_", Sys.Date(), ".csv"))
writeLines("second figure", file.path(stage, "figures", "ECSM_Level_1_test.pdf"))
second <- publish(2, 20)
stopifnot(identical(first, second), length(pipeline_files(out, p)) == 1,
  length(pipeline_files(file.path(out, "archive"), p)) == 1,
  pipeline_read_table(pipeline_files(file.path(out, "archive"), p))$data$value == 10,
  nrow(pipeline_read_table(second)$data) == 2)
# Reprocessing replaces matching period figures; no figure archive is created.
publish(2, 25, TRUE)
stopifnot(identical(pipeline_read_table(second)$data$value, c(10, 25)),
  identical(pipeline_read_table(pipeline_files(file.path(out, "archive"), p))$data$value, c(10, 20)),
  !dir.exists(file.path(out, "figures", "figure_archive")),
  length(list.dirs(file.path(out, "figures"), recursive = FALSE)) == 2)
# Old figures and legacy runs are left untouched for manual organization.
legacy <- file.path(out, "runs", "old", "figures"); dir.create(legacy, recursive = TRUE)
legacy_files <- c(file.path(legacy, "old.pdf"), file.path(out, "figures", "old.jpeg"),
                  file.path(out, "figures", "ECSM_Level_4_fixture.pdf"))
for (f in legacy_files) writeLines("legacy figure", f)
before_legacy <- tools::md5sum(legacy_files)
publish(3, 30)
stopifnot(identical(before_legacy, tools::md5sum(legacy_files)))
# Inject a publication failure AFTER rotation. Previous tables and figures return.
before <- pipeline_read_table(second)$data
real_copy <- file.copy
file.copy <- function(from, to, ...) {
  if (basename(from) == "current.csv") return(FALSE)
  real_copy(from, to, ...)
}
err <- tryCatch(publish(3, 99, TRUE), error = identity)
rm(file.copy)
stopifnot(inherits(err, "error"), isTRUE(all.equal(before, pipeline_read_table(second)$data)),
  length(pipeline_files(out, p)) == 1,
  !length(list.files(out, pattern = "^\\.publish-", all.files = TRUE)))
# All L4 data versions survive, including same-day runs and legacy root files.
p <- pipeline_patterns[["L4"]]
real_figures <- pipeline_level4_figures
pipeline_level4_figures <- function(data, site, folder) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  writeLines(as.character(nrow(data)), file.path(folder, paste0(site, "_Level_4_fixture.pdf")))
}
# Staging contains only this stage's figures.
unlink(list.files(file.path(stage, "figures"), full.names = TRUE))
l4 <- publish(1, 1)
publish(2, 2)
archived <- pipeline_files(file.path(out, "archive"), p)
stopifnot(length(archived) == 1)
first_hash <- tools::md5sum(archived)
publish(3, 3)
publish(3, 4, TRUE)
stopifnot(length(pipeline_files(file.path(out, "archive"), p)) == 3,
  identical(first_hash, tools::md5sum(names(first_hash))),
  length(pipeline_files(out, p)) == 1,
  readLines(file.path(out, "figures", "ECSM_Level_4_fixture_full_period.pdf")) == "3")
period <- file.path(out, "figures", paste0("period_", pipeline_period_label(list(pending = t[3]))))
stopifnot(readLines(file.path(period, "ECSM_Level_4_fixture.pdf")) == "1",
          identical(before_legacy, tools::md5sum(legacy_files)))
# Multiple legacy root versions must each survive byte-for-byte in L4 archive.
legacy_csv <- file.path(out, paste0("ECSM_", p, "_2020-01-01.csv"))
stopifnot(file.copy(l4, legacy_csv))
legacy_hash <- unname(tools::md5sum(legacy_csv))
publish(3, 5, TRUE)
stopifnot(length(pipeline_files(file.path(out, "archive"), p)) == 5,
  identical(legacy_hash, unname(tools::md5sum(file.path(out, "archive", basename(legacy_csv))))))
# Failure after rotation must restore ALL bytes, including overwritten figures.
tracked <- list.files(out, recursive = TRUE, full.names = TRUE)
tracked <- tracked[!dir.exists(tracked)]
checksums <- tools::md5sum(tracked)
file.copy <- function(from, to, ...) {
  if (basename(from) == "current.csv") return(FALSE)
  real_copy(from, to, ...)
}
err <- tryCatch(publish(3, 99, TRUE), error = identity)
rm(file.copy)
stopifnot(inherits(err, "error"), identical(checksums, tools::md5sum(tracked)),
          length(pipeline_files(file.path(out, "archive"), p)) == 5)
pipeline_level4_figures <- real_figures
# QC percentages distinguish flag codes from actual removal, and zero gaps
# produce NA rather than an invented good-quality percentage.
x <- data.frame(LE = c(10, 20, NA, 40), LE_QC_SSITC = c(0, 2, NA, 1),
                LE_filtered2 = c(10, NA, NA, 40), LE_PI_F = c(10, 22, 33, 40))
d <- pipeline_flux_diagnostics(x)
stopifnot(d$LE_raw_observed_pct == 75, d$LE_qc_removed_n == 1,
  d$LE_qc_removed_denominator_n == 3, d$LE_gapfilled_pct == 50,
  d$LE_gapfill_success_pct == 100, d$LE_QC_SSITC_code_2_n == 1)
options(ec.pipeline = list(window = list(pending = t)))
pipeline_capture_mds(data.frame(TIMESTAMP = t, LE = c(1, NA, NA)),
  data.frame(LE_f = c(1, 2, 3), LE_fqc = c(0, 1, 2)), c(LE = "LE"))
stopifnot(getOption("ec.pipeline")$diagnostics$LE_good_gapfilled_pct == 50)
options(ec.pipeline = NULL)
unlink(root, recursive = TRUE)
cat("PASS: same-day rotation, L1 previous version, all L4 versions, period and cumulative figures, legacy preservation, rollback, and diagnostic denominators.\n")
