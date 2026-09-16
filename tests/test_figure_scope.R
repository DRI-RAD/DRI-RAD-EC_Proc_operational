# Exercise the actual Level 2 figure call and function without its science loop.
source("run_pipeline.R")
suppressPackageStartupMessages({ library(tidyverse); library(cowplot) })
root <- tempfile("figure-scope-"); dir.create(file.path(root, "figures"), recursive = TRUE)
exprs <- parse(pipeline_stages[["L2"]])
for (expr in exprs) if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
                       identical(expr[[2]], as.name("plot_gapfilled_pdf"))) eval(expr)
calls <- list()
walk <- function(x) {
  if (!is.call(x) && !is.expression(x)) return(invisible(NULL))
  if (is.call(x) && identical(x[[1]], as.name("plot_gapfilled_pdf"))) calls[[length(calls)+1L]] <<- x
  for (child in as.list(x)) {
    if (missing(child)) next
    if (is.call(child) || is.expression(child)) walk(child)
  }
}
walk(exprs)
stopifnot(length(calls) == 1L)
t <- seq(pipeline_time("2026-08-01 00:30"), by = 1800, length.out = 96)
df_PI <- data.frame(TIMESTAMP = t)
for (v in c("TA_PI_F", "RH_PI_F", "SW_IN_PI_F", "LW_IN_PI_F", "P_PI_F")) {
  df_PI[[v]] <- sin(seq_along(t)/10)
  df_PI[[paste0(v, "_method")]] <- rep(0:3, length.out = length(t))
}
options(ec.pipeline = list(window = list(pending = tail(t, 48))))
base_dir <- "Z:/unused-root"
dirs_use <- data.frame(dir_output = normalizePath(root, winslash = "/"))
i <- "ECSM"
real_plot <- plot_gapfilled_pdf
plot_gapfilled_pdf <- function(df, ...) {
  stopifnot(nrow(df) == 48, min(df$TIMESTAMP) == t[49])
  real_plot(df, ...)
}
eval(calls[[1]])
f <- file.path(root, "figures", "ECSM_Level_2_Figure_gapfilled_PI.pdf")
stopifnot(file.exists(f), file.info(f)$size > 1000)
options(ec.pipeline = NULL)
unlink(root, recursive = TRUE)
cat("PASS: actual Level 2 plot call, pending-only figure inputs, and absolute staging path.\n")
