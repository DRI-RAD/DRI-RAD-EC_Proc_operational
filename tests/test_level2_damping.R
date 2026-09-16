# Run from the repository root: Rscript tests/test_level2_damping.R
# Execute the actual soil reconstruction loop with synthetic data, not NAS data.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(lubridate); library(ggplot2); library(zoo) })
stage_path <- Sys.getenv("EC_LEVEL2_TEST_SCRIPT", "R/Level_2_met_PI_vars_new.R")
code <- parse(stage_path)
helper_names <- c("get_peak_hour", "level2_fill_damping", "level2_smooth_damping", "level2_lead_or_na")
for (expr in code) if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
                      as.character(expr[[2]])[1] %in% helper_names) eval(expr)
loops <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("for")) &&
                  any(grepl("level2_smooth_damping", deparse(x), fixed = TRUE)), as.list(code))
stopifnot(length(loops) == 1)

# Well-supported windows retain exactly the former smoothing results.
x <- 0.15 + 0.02*sin(seq_len(60)/5)
x[c(1, 8, 60)] <- NA
former <- zoo::na.approx(x, na.rm = FALSE, rule = 2)
former <- zoo::rollmean(former, k = 30, fill = NA, align = "center")
former <- zoo::na.approx(former, na.rm = FALSE, rule = 2)
stopifnot(isTRUE(all.equal(level2_smooth_damping(x), former)),
          identical(level2_lead_or_na(1:10, 2), dplyr::lead(1:10, 2L)))
for (n in c(NA_real_, NaN, Inf, -1, 1.5, 100))
  stopifnot(all(is.na(level2_lead_or_na(1:10, n))))
stopifnot(all(is.na(level2_smooth_damping(rep(NA_real_, 16)))),
          identical(level2_smooth_damping(c(NA, 0.2, NA)), rep(0.2, 3)))

root <- tempfile("damping-tests-"); dir.create(root)
pdf(file.path(root, "plots.pdf"))
all_sites <- "ECSM"
dirs <- data.frame(site = "ECSM")
meta_G <- data.frame(site = "ECSM", Sun_plate_depth_1 = 0.08,
                     Shade_plate_depth_1 = 0.08, G_sun_weight = 0.5)
pipeline_figure_data <- identity
run_case <- function(days, missing_sun = FALSE, one_peak = FALSE, midnight_only = FALSE) {
  t <- seq(as.POSIXct("2026-07-16", tz = "UTC"), by = 1800, length.out = days*48)
  if (midnight_only) t <- t[1]
  h <- hour(t) + minute(t)/60
  flux <- 30 + 50*cos((h - 14)/24*2*pi)
  input <- data.frame(TIMESTAMP = t, G_plate_sun = flux, G_plate_shade = flux*0.8,
                      G_calorimetric = flux+5, G_plate = flux)
  if (missing_sun) input$G_plate_sun <- NA_real_
  if (one_peak) input$G_plate_sun[as.Date(t) != as.Date(t[1])] <- NA_real_
  site_list <- list(ECSM = input)
  eval(loops[[1]])
  result <- site_list$ECSM
  stopifnot(nrow(result) == nrow(input), identical(as.numeric(result$TIMESTAMP), as.numeric(input$TIMESTAMP)),
            identical(result$G_PI, input$G_calorimetric))
  if (missing_sun || midnight_only) stopifnot(all(is.na(result$G_sun_leuning)))
  else stopifnot(any(is.finite(result$G_sun_leuning)))
  if (!midnight_only) stopifnot(any(is.finite(result$G_shade_leuning)))
}
for (days in c(1, 16, 29, 30, 31, 60)) run_case(days)
run_case(16, missing_sun = TRUE)
run_case(16, one_peak = TRUE)
run_case(1, midnight_only = TRUE)
dev.off()
unlink(root, recursive = TRUE)
cat("PASS: actual Level 2 harmonic loop on 1/16/29/30/31/60 days, missing peaks, one peak, midnight-only input; long-window parity and unchanged G_PI.\n")
