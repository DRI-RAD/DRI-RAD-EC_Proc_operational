# Run from the repository root. Synthetic CSVs never access production inputs.
source("run_pipeline.R")
root <- tempfile("eddypro-years-")
dir.create(root)
fields <- c("date", "time", "LE", "qc_LE", "H", "qc_H", "co2_flux", "qc_co2_flux", "Tau", "qc_Tau",
  "co2_mole_fraction", "co2_molar_density", "co2_var", "h2o_mole_fraction", "h2o_molar_density", "h2o_var",
  "ts_var", "u_var", "v_var", "w_var", "sonic_temperature", "wind_speed", "max_wind_speed", "wind_dir",
  "u*", "(z-d)/L", "L", "x_90%", "x_70%", "x_50%", "x_30%", "x_10%")
write_ep <- function(path, date, value) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  x <- as.data.frame(setNames(rep(list("1"), length(fields)), fields), check.names = FALSE)
  x$date <- date; x$time <- "00:30"; x$LE <- as.character(value)
  writeLines(c("category", paste(fields, collapse = ","), paste(rep("unit", length(fields)), collapse = ",")), path)
  readr::write_csv(x, path, append = TRUE, col_names = FALSE)
}
for (site in eval(formals(run_pipeline)$sites)) {
  folder <- file.path(root, "EddyPro_Output_Files_YK", site)
  old <- file.path(folder, "2025", "site_full_output_2026-10-01T120000.csv")
  first <- file.path(folder, "2026_01", "eddypro_site_full_output_2026-10-02T120000_adv.csv")
  second <- file.path(folder, "2026_02", "custom_full_output_2026-10-03T120000.csv")
  write_ep(old, "2025-12-31", 10)
  write_ep(first, "2026-01-01", 20)
  write_ep(second, "2026-01-01", -9999)
  # Copy dates cannot override a newer explicit run timestamp.
  Sys.setFileTime(first, Sys.time() + 86400)
  writeLines("not a full output", file.path(folder, "2026_02", "biomet.csv"))
  x <- pipeline_eddypro(folder)
  stopifnot(length(x$files) == 3, nrow(x$data) == 2, x$data$LE[1] == 10,
            is.na(x$data$LE[2]), !anyDuplicated(x$data$TIMESTAMP),
            pipeline_eddypro_end(x) == pipeline_time("2026-01-01 00:30"))
  w <- list(read_start = pipeline_time("2026-01-01 00:30"), end = pipeline_time("2026-01-01 00:30"))
  stopifnot(nrow(pipeline_eddypro(folder, w)$data) == 1)
}
# Files without a run date use modification time, then path for ties.
folder <- file.path(root, "undated")
a <- file.path(folder, "2026_01", "custom_full_output_a.csv")
b <- file.path(folder, "2026_02", "custom_full_output_b.csv")
write_ep(a, "2026-01-01", 1); write_ep(b, "2026-01-01", 2)
Sys.setFileTime(a, as.POSIXct("2026-10-01", tz = "UTC"))
Sys.setFileTime(b, as.POSIXct("2026-10-02", tz = "UTC"))
stopifnot(pipeline_eddypro(folder)$data$LE == 2)
unlink(root, recursive = TRUE)
cat("PASS: seven sites, year/run folders, custom filenames, overlap/NA precedence, window and end.\n")
