source("run_pipeline.R")
root <- tempfile("li710-sources-"); dir.create(root)
t <- seq(pipeline_time("2026-07-01 00:00"), by = 1800, length.out = 4)
units_for <- function(x) as.data.frame(setNames(lapply(names(x), function(n) c("unit", "avg")), names(x)))
write_raw <- function(name, indexes, le, qc = TRUE) {
  x <- data.frame(TIMESTAMP = t[indexes], LE_710 = le, H_710 = 2, diag = 0, flow = 200, tilt = 0)
  if (qc) x$data_qc <- 1
  pipeline_write_table(x, units_for(x), file.path(root, name), "fixture")
}
write_raw("Site_LI710.dat.backup", 1:2, c(10, 20), FALSE)
write_raw("Site_LI710.202607.backup", 2:3, c(25, 30))
write_raw("Site_LI710.dat", 3:4, c(300, 400))
write_raw("Another_LI710.dat", 1, 999)
row <- data.frame(dir_LI710 = file.path(root, "Site_LI710.dat"), dir_LI710_old = NA_character_)
files <- pipeline_li710_files(row, root)
stopifnot(length(files) == 3, basename(tail(files, 1)) == "Site_LI710.dat")
x <- pipeline_li710_read(files)$data
stopifnot(nrow(x) == 4, identical(x$LE_710, c(10, 25, 300, 400)), x$data_qc[1] == 0)
slice <- pipeline_li710_read(files, list(read_start = t[2], end = t[3]))$data
stopifnot(nrow(slice) == 2)
unlink(file.path(root, "Site_LI710.dat"))
stopifnot(length(pipeline_li710_files(row, root)) == 2,
          nrow(pipeline_li710_read(pipeline_li710_files(row, root))$data) == 3)
unlink(root, recursive = TRUE)
stopifnot(length(pipeline_li710_files(row, root)) == 0,
          nrow(pipeline_li710_read(character())$data) == 0)
cat("PASS: LI710 dat/backup discovery, site isolation, overlap precedence, optional QC, period slicing, backup-only and absent inputs.\n")
