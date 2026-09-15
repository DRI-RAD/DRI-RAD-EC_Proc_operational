# Exercise the actual Level 3 EC adapter and its unchanged QAQC using synthetic
# half-hour observations. Run from the project root: Rscript tests/smoke_level3.R
source("run_pipeline.R")
root <- tempfile("level3-smoke-")
dir.create(root)
t <- seq(pipeline_time("2026-08-01 00:30"), pipeline_time("2026-09-01 00:00"), by = 1800)
n <- length(t)
solar <- pmax(0, 600 * sin(2 * pi * ((seq_len(n) %% 48) / 48 - 0.25)))
make_units <- function(data) as.data.frame(setNames(lapply(names(data), function(name) c(if(name == "TIMESTAMP") "TS" else "unit", NA_character_)), names(data)), check.names = FALSE)
out <- file.path(root, "out")
met1 <- data.frame(TIMESTAMP = t, CO2_sig_strgth_Min = 0.9, H2O_sig_strgth_Min = 0.9, PotRad = solar)
met2 <- data.frame(TIMESTAMP = t, SW_IN_PI_F = solar, TA_PI_F = 20, RH_PI_F = 40, VPD_PI_F = 10)
pipeline_write_table(met1, make_units(met1), file.path(out, "ECSM_Level_1_QAQCed_logger_data_2026-09-01.csv"), "synthetic")
pipeline_write_table(met2, make_units(met2), file.path(out, "ECSM_Level_2_PI_vars_logger_data_2026-09-01.csv"), "synthetic")
easy_names <- c("LE", "LE_SSITC_TEST", "H", "H_SSITC_TEST", "FC", "FC_SSITC_TEST", "TAU", "TAU_SSITC_TEST",
 "CO2", "CO2_SIGMA", "H2O", "H2O_SIGMA", "T_SONIC", "T_SONIC_SIGMA", "WS", "WS_MAX", "WD", "USTAR", "ZL", "MO_LENGTH",
 "U_SIGMA", "V_SIGMA", "W_SIGMA", "FETCH_MAX", "FETCH_90", "FETCH_55", "FETCH_40")
easy <- as.data.frame(setNames(lapply(easy_names, function(x) rep(1, n)), easy_names))
easy <- cbind(data.frame(TIMESTAMP = t), easy)
easy$LE <- solar / 2 + sin(seq_len(n))
easy$H <- solar / 3 + sin(seq_len(n))
easy$FC <- -solar / 100 + sin(seq_len(n))
raw <- file.path(root, "easy.dat")
pipeline_write_table(easy, make_units(easy), raw, "TOA5 synthetic")
fields <- c("LE", "qc_LE", "H", "qc_H", "co2_flux", "qc_co2_flux", "Tau", "qc_Tau",
 "co2_mole_fraction", "co2_molar_density", "co2_var", "h2o_mole_fraction", "h2o_molar_density", "h2o_var",
 "ts_var", "u_var", "v_var", "w_var", "sonic_temperature", "wind_speed", "max_wind_speed", "wind_dir",
 "u*", "(z-d)/L", "L", "x_90%", "x_70%", "x_50%", "x_30%", "x_10%")
ep <- as.data.frame(setNames(lapply(fields, function(x) rep(1, n)), fields), check.names = FALSE)
ep$LE <- easy$LE
ep$H <- easy$H
ep$co2_flux <- easy$FC
ep$date <- format(t, "%Y-%m-%d", tz = "UTC")
ep$time <- format(t, "%H:%M", tz = "UTC")
epdir <- file.path(root, "Columbus", "Aug_2026", "output")
dir.create(epdir, recursive = TRUE)
epfile <- file.path(epdir, "eddypro_Columbus_full_output_2026-09-01T120000_adv.csv")
writeLines("category row", epfile)
write.table(matrix(names(ep), nrow = 1), epfile, sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
write.table(matrix(rep("unit", ncol(ep)), nrow = 1), epfile, sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
readr::write_csv(ep, epfile, append = TRUE, col_names = FALSE)
cfg <- data.frame(site = "ECSM", dir_output = out, dir_eddypro = file.path(root, "Columbus"),
                  dir_Easyflux = raw, dir_LI710 = NA_character_, dir_LI710_old = NA_character_)
config <- file.path(root, "config.csv")
readr::write_csv(cfg, config)
run_pipeline("ECSM", stages = "L3_EC", config_file = config, base_dir = root)
result <- pipeline_history(out, pipeline_patterns[["L3_EC"]])$data
stopifnot(nrow(result) == n, all(result$processing == "EddyPro"),
          all(c("LE_QC_despike", "FC_QC_longrun") %in% names(result)))
stopifnot(run_pipeline("ECSM", stages = "L3_EC", config_file = config, base_dir = root)[[1]] == "skipped")
unlink(root, recursive = TRUE)
cat("PASS: actual Level 3 CSV ingestion, QAQC, output publication, and no-new-data rerun.\n")
