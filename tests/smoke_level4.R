# Run the real Level 4 MDS and all five restored diagnostic figures.
# Uses temporary synthetic upstream products; never reads the network data root.
source("run_pipeline.R")
root <- tempfile("level4-smoke-")
dir.create(root)
t <- seq(pipeline_time("2026-05-01 00:30"), by = 1800, length.out = 110*48)
n <- length(t)
set.seed(42)
sun <- pmax(0, sin(2*pi*((seq_len(n) %% 48)/48 - 0.25)))
season <- 1 + 0.15*sin(seq_len(n)/(48*12))
rad <- 650 * sun * season
make_units <- function(df) as.data.frame(setNames(lapply(names(df), function(x) c(if(x == "TIMESTAMP") "TS" else "unit", NA_character_)), names(df)), check.names = FALSE)
drop <- c("T_DP_1_1_1", "T_DP_1_1_2", "T_DP_1_1_3", "T_nr", "R_LW_in_meas", "R_LW_out_meas", "sun_azimuth",
          "sun_elevation", "hour_angle", "sun_declination", "air_mass_coeff", "daytime", "WS_05103_rslt",
          "WS_05103_5s_TMx", "WS_05103_5s_Min", "WS_05103_5s_TMn")
l1 <- data.frame(TIMESTAMP = t, PotRad = rad, NETRAD = rad*0.7-30,
                 SWC_1_1_1 = 20 + sin(seq_len(n)/200), SWC_1_1_2 = 22 + cos(seq_len(n)/250))
for (name in drop) l1[[name]] <- 0
l2 <- data.frame(TIMESTAMP = t, TA_PI_F = 15+10*sun + sin(seq_len(n)/100), RH_PI_F = 60-20*sun,
                 VPD_PI_F = 5+10*sun, SW_IN_PI_F = rad, G_PI = rad*0.08,
                 P_PI_F = ifelse(seq_len(n) %% 432 == 0, 2, 0))
ec <- data.frame(TIMESTAMP = t, LE = rad*0.35+rnorm(n, 0, 3), H = rad*0.18+rnorm(n, 0, 2),
                 FC = 2-10*sun+rnorm(n, 0, 0.1), USTAR = 0.25+0.2*sun,
                 H2O = 12, H2O_SIGMA = 0.2+0.1*sun, T_SONIC_SIGMA = 0.3+0.2*sun,
                 SLE = NA_real_, SH = NA_real_, SC = NA_real_, processing = "EddyPro")
for (var in c("LE", "H", "FC")) for (flag in c("despike", "sig_str", "fetch", "longrun", "SSITC", "ustar"))
  ec[[paste(var, "QC", flag, sep = "_")]] <- 0
li <- data.frame(TIMESTAMP = t, LE_710 = ec$LE*0.95+rnorm(n), H_710 = ec$H*1.03+rnorm(n))
for (var in c("LE_710", "H_710")) for (flag in c("diag", "despike", "longrun", "flow"))
  li[[paste(var, "QC", flag, sep = "_")]] <- 0
# Missing fluxes at the beginning must not truncate any MDS input below 100 days.
ec$LE[1:(12*48)] <- NA_real_
li$LE_710[1:(12*48)] <- NA_real_
ec$LE[(105*48+3):(105*48+6)] <- NA_real_
out <- file.path(root, "outputs")
for (stage in if (Sys.getenv("EC_NO_LI710") == "1") c("L1", "L2", "L3_EC") else c("L1", "L2", "L3_EC", "L3_LI710")) {
  df <- switch(stage, L1 = l1, L2 = l2, L3_EC = ec, L3_LI710 = li)
  if (stage == "L3_LI710" && Sys.getenv("EC_LI710_PARTIAL") == "1") {
    df <- tail(df, 96)
    df$H_710 <- NA_real_
  }
  pipeline_write_table(df, make_units(df), file.path(out, paste0("ECSM_", pipeline_patterns[[stage]], "_fixture.csv")), "synthetic")
}
cfg <- data.frame(site = "ECSM", dir_output = out, dir_eddypro = NA_character_, dir_LI710 = if (Sys.getenv("EC_NO_LI710") == "1") NA_character_ else "configured")
config <- file.path(root, "config.csv")
readr::write_csv(cfg, config)
args <- list(sites = "ECSM", stages = "L4", base_dir = root, config_file = config,
             start = t[100*48+1], end = tail(t,1))
preview <- do.call(run_pipeline, c(args, list(dry_run = TRUE)))
stopifnot(preview[[1]]$calculation_days == 100, length(preview[[1]]$pending) == 480)
run <- do.call(run_pipeline, args)
product <- pipeline_read_table(run[[1]])$data
stopifnot(nrow(product) == 480, min(product$TIMESTAMP) == args$start,
          sum(!is.na(product$LE_PI_F)) > 470)
figures <- list.files(file.path(out, "figures"), pattern = "Level_4_Figure.*jpeg$", recursive = TRUE, full.names = TRUE)
stopifnot(length(figures) == 10, all(file.info(figures)$size > 10000))
stopifnot(length(list.files(file.path(out, "figures"), pattern = "Level_4_Figure.*jpeg$")) == 5,
          length(list.dirs(file.path(out, "figures"), recursive = FALSE)) == 1,
          !any(grepl("automatic_plots", list.files(root, recursive = TRUE))))
if (Sys.getenv("EC_NO_LI710") == "1") stopifnot(all(is.na(product$LE_710_PI_F)), all(is.na(product$ET_710_F)))
period <- readr::read_csv(paste0(run[[1]], ".period.csv"), show_col_types = FALSE)
stopifnot(period$LE_raw_observed_n == 476, period$LE_gapfilled_n == 4, period$LE_good_gapfilled_denominator_n == 4, period$new_rows == 480, period$calculation_days == 100, period$figure_period == "new")
stopifnot(do.call(run_pipeline, args)[[1]] == "skipped")
# Optional local visual QA artifacts, outside the repository under test.
artifacts <- Sys.getenv("EC_TEST_ARTIFACT_DIR")
if (nzchar(artifacts)) {
  dir.create(artifacts, recursive = TRUE, showWarnings = FALSE)
  file.copy(figures, artifacts)
}
unlink(root, recursive = TRUE)
cat("PASS: real Level 4 MDS on 100 days, 10-day-only publication, five period plus five cumulative figures, and skipped rerun.\n")
