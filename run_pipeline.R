# Source this file once, then call run_pipeline(). All settings are arguments;
# the runner never clears the user's workspace or installs packages.
.pipeline_sources <- Filter(function(x) !is.null(x), lapply(sys.frames(), function(x) x$ofile))
.pipeline_script <- if (length(.pipeline_sources)) tail(.pipeline_sources, 1)[[1]] else NA_character_
if (!length(.pipeline_script) || is.na(.pipeline_script)) {
  .pipeline_arg <- grep("^--file=", commandArgs(), value = TRUE)
  .pipeline_script <- if (length(.pipeline_arg)) sub("^--file=", "", .pipeline_arg[1]) else "run_pipeline.R"
}
.pipeline_root <- dirname(normalizePath(.pipeline_script, winslash = "/", mustWork = TRUE))
source(file.path(.pipeline_root, "R", "pipeline_io.R"))

pipeline_stages <- c(
  L1 = "R/Level_1_met_QAQC_new.R",
  L2 = "R/Level_2_met_PI_vars_new.R",
  L3_EC = "R/Level_3_EC_data_QAQC_new.R",
  L3_LI710 = "R/Level_3_LI710_QAQC_new.R",
  L4 = "R/Level_4_gapfilling_ET_calculation_new.R"
)
pipeline_patterns <- c(
  L1 = "Level_1_QAQCed_logger_data",
  L2 = "Level_2_PI_vars_logger_data",
  L3_EC = "Level_3_EC_with_qaqc",
  L3_LI710 = "Level_3_li710_with_qaqc",
  L4 = "Level_4_post_processed_data"
)

pipeline_publish <- function(staging, output, pattern, window, reprocess) {
  produced <- pipeline_files(staging, pattern)
  # Upstream files in staging have a different stage pattern.
  if (length(produced) != 1L) stop("Expected one final stage output: ", pattern)
  fresh <- pipeline_read_table(produced)
  if (!all(window$pending %in% fresh$data$TIMESTAMP))
    stop("Stage did not produce every planned timestamp; output was not committed.")
  history <- pipeline_history(output, pattern)
  data <- fresh$data
  if (!is.null(history)) {
    if (!setequal(names(history$data), names(data))) stop("Historical output schema differs: ", pattern)
    old <- history$data
    if (reprocess) old <- old[!old$TIMESTAMP %in% data$TIMESTAMP, ]
    else if (any(data$TIMESTAMP %in% old$TIMESTAMP)) stop("Unexpected overlap during append.")
    data <- dplyr::bind_rows(old, data)
  }
  data <- data[order(data$TIMESTAMP), ]
  version <- basename(tempfile(paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "-")))
  site <- sub(paste0("_", pattern, ".*"), "", basename(produced), ignore.case = TRUE)
  target <- file.path(output, paste0(site, "_", pattern, "_", version, ".csv"))
  pipeline_write_table(data, fresh$units, target, "Cumulative processed output; prior rows preserved outside the requested update.")
  # The final CSV is the checkpoint. Audit metadata is informative, never the
  # authority for skipping data after interruption or an incomplete write.
  audit <- data.frame(stage = pattern, output = basename(target),
                      first_timestamp = min(data$TIMESTAMP), last_timestamp = max(data$TIMESTAMP),
                      new_first = min(fresh$data$TIMESTAMP), new_last = max(fresh$data$TIMESTAMP),
                      new_rows = nrow(fresh$data), total_rows = nrow(data),
                      context_start = window$read_start, reprocess = reprocess)
  readr::write_csv(audit, paste0(target, ".period.csv"))
  # Keep diagnostics per run so old figures are not silently replaced.
  diagnostics <- file.path(output, "runs", version)
  dir.create(diagnostics, recursive = TRUE, showWarnings = FALSE)
  file.copy(list.files(staging, full.names = TRUE), diagnostics, recursive = TRUE)
  message("Committed ", pattern, ": ", nrow(fresh$data), " rows -> ", target)
  target
}

run_pipeline <- function(sites = c("ECDP", "EDVG", "EDVP", "ERVA", "ERVP", "ECSM"), start = NULL, end = NULL,
                         base_dir = "Z:/NWI/Task_3_Water_Use_ET_and_Meteoroligical_Monitoring/networks/eddy_stations",
                         config_file = file.path(.pipeline_root, "meta_files", "file_directions_new.csv"),
                         stages = names(pipeline_stages), reference_sites = NULL,
                         context_days = 0, reprocess = FALSE, dry_run = FALSE) {
  # start/end are inclusive half-hour observation timestamps, not file dates.
  # context_days adds preceding observations for boundary calculations only;
  # already-processed values are never rewritten unless reprocess = TRUE.
  if (!is.numeric(context_days) || length(context_days) != 1L ||
      !is.finite(context_days) || context_days < 0) stop("context_days must be nonnegative.")
  if (reprocess && (is.null(start) || is.null(end))) stop("Reprocessing requires explicit start and end.")
  if (!is.null(start)) start <- pipeline_time(start)
  if (!is.null(end)) end <- pipeline_time(end)
  if ((!is.null(start) && length(start) != 1L) || (!is.null(end) && length(end) != 1L))
    stop("start and end must each contain one timestamp.")
  if (!is.null(start) && !is.null(end) && start > end) stop("start is after end.")
  if (!length(stages) || any(!stages %in% names(pipeline_stages))) stop("Unknown or empty stage selection.")
  stages <- names(pipeline_stages)[names(pipeline_stages) %in% stages]
  required <- c("tidyverse", "lubridate", "zoo", "cowplot", "openeddy", "bigleaf", "REddyProc", "ggpmisc", "gridExtra", "Metrics")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages first: ", paste(missing, collapse = ", "))
  dirs <- readr::read_csv(config_file, show_col_types = FALSE)
  supported <- c("ECDP", "EDVG", "EDVP", "ERVA", "ERVP", "ECSM")
  if (!length(sites) || any(!sites %in% supported) || any(!sites %in% dirs$site))
    stop("The new runner supports the six main-network sites; variation scripts require separate adapters.")
  if (anyDuplicated(dirs$site)) stop("Duplicate site identifiers in configuration.")
  if (!"dir_eddypro" %in% names(dirs)) stop("Configuration must contain dir_eddypro.")
  sites <- unique(sites)
  # Preserve the original Level 2 cross-site reference pool by default.
  if (is.null(reference_sites)) reference_sites <- if ("L2" %in% stages) supported else sites
  reference_sites <- unique(c(sites, reference_sites))
  if (any(!reference_sites %in% supported) || any(!reference_sites %in% dirs$site)) stop("Unknown reference site.")
  oldwd <- getwd()
  oldopt <- options(ec.pipeline = NULL)
  on.exit({ setwd(oldwd); options(oldopt) }, add = TRUE)
  setwd(.pipeline_root)
  output_for <- function(site) pipeline_path(base_dir, dirs$dir_output[match(site, dirs$site)])
  history_for <- function(site, stage) pipeline_history(output_for(site), pipeline_patterns[[stage]])
  projected <- list()
  times_for <- function(site, stage) {
    key <- paste(site, stage, sep = ":")
    if (dry_run && !is.null(projected[[key]])) return(projected[[key]])
    x <- history_for(site, stage)
    if (is.null(x)) stop("Missing ", stage, " output for ", site, ". Run upstream stages first.")
    x$data$TIMESTAMP
  }
  # Load raw clocks once. Reading historical source rows is necessary to discover
  # coverage, but scientific processing is restricted to the selected window.
  logger <- eddy <- li710 <- list()
  for (site in sites) {
    row <- dirs[dirs$site == site, ]
    if ("L1" %in% stages) logger[[site]] <- pipeline_raw_times(pipeline_path(base_dir, row$dir_met))
    if ("L3_EC" %in% stages) eddy[[site]] <- pipeline_eddypro(pipeline_path(base_dir, row$dir_eddypro))
    if (any(c("L3_LI710", "L4") %in% stages) && is.na(row$dir_LI710))
      stop("The existing LI710/Level 4 algorithms require LI710 input for ", site)
    if ("L3_LI710" %in% stages) {
      t <- pipeline_raw_times(pipeline_path(base_dir, row$dir_LI710))
      if (!is.na(row$dir_LI710_old)) t <- c(t, pipeline_raw_times(pipeline_path(base_dir, row$dir_LI710_old)))
      li710[[site]] <- sort(unique(lubridate::floor_date(t, "30 minutes")))
    }
  }
  results <- list()
  # Stage-major order makes every selected site's new L1 available to L2.
  for (stage in stages) for (site in sites) {
    row <- dirs[dirs$site == site, ]
    output <- output_for(site)
    # A per-site directory lock prevents concurrent read/merge/write races.
    if (!dry_run) {
      dir.create(output, recursive = TRUE, showWarnings = FALSE)
      lock <- file.path(output, ".pipeline-lock")
      if (!dir.create(lock, showWarnings = FALSE)) stop("Site is locked by another run: ", lock)
    }
    tryCatch({
      available <- switch(stage,
        L1 = seq(min(logger[[site]]), max(logger[[site]]), by = 1800),
        L2 = times_for(site, "L1"),
        L3_EC = intersect(as.numeric(times_for(site, "L1")), as.numeric(times_for(site, "L2"))),
        L3_LI710 = times_for(site, "L1"),
        L4 = Reduce(intersect, lapply(c("L1", "L2", "L3_EC", "L3_LI710"),
                                      function(s) as.numeric(times_for(site, s))))
      )
      if (is.numeric(available) && !inherits(available, "POSIXt"))
        available <- as.POSIXct(available, origin = "1970-01-01", tz = "UTC")
      # Match the original Level 1 deployment cutoff without altering its QAQC.
      if (stage == "L1" && site == "EDVG")
        available <- available[available >= pipeline_time("2023-09-26 11:30:00")]
      # Stop at the flux source's actual coverage. Internal missing observations
      # remain in the regular grid for the original fallback/gap-filling logic.
      if (stage %in% c("L3_EC", "L3_LI710")) {
        raw_times <- if (stage == "L3_EC") eddy[[site]]$data$TIMESTAMP else li710[[site]]
        available <- available[available >= min(raw_times) & available <= max(raw_times)]
      }
      window <- pipeline_plan(available, history_for(site, stage), start, end, context_days, reprocess)
      key <- paste(site, stage, sep = ":")
      if (is.null(window)) {
        message(key, ": no unprocessed timestamps.")
        results[[key]] <- "skipped"
      } else {
        message(key, ": ", length(window$pending), " pending rows; ", window$start, " through ", window$end,
                "; calculation starts ", window$read_start)
        if (dry_run) {
          results[[key]] <- window
          prior <- history_for(site, stage)
          projected[[key]] <- sort(unique(c(prior$data$TIMESTAMP, window$pending)))
        } else {
          staging_root <- tempfile("ec-pipeline-")
          dir.create(staging_root)
          tryCatch({
            stage_dirs <- dirs
            refs <- if (stage == "L2") reference_sites else site
            upstream <- switch(stage, L1 = character(), L2 = "L1", L3_EC = c("L1", "L2"),
                               L3_LI710 = "L1", L4 = c("L1", "L2", "L3_EC", "L3_LI710"))
            for (ref in refs) {
              folder <- file.path(staging_root, ref)
              dir.create(file.path(folder, "figures"), recursive = TRUE)
              stage_dirs$dir_output[stage_dirs$site == ref] <- folder
              for (up in upstream) {
                item <- history_for(ref, up)
                if (is.null(item)) stop("Missing ", up, " reference input for ", ref)
                item$data <- pipeline_slice(item$data, window)
                if (!nrow(item$data)) stop("Reference input has no data in this window: ", ref)
                pipeline_write_table(item$data, item$units,
                    file.path(folder, paste0(ref, "_", pipeline_patterns[[up]], "_input.csv")),
                    "Merged upstream input for the current calculation window.")
              }
            }
            options(ec.pipeline = list(dirs = stage_dirs, base_dir = base_dir,
                sites = site, reference_sites = refs, window = window, eddypro = eddy[[site]]))
            env <- new.env(parent = environment(run_pipeline))
            # Route automatic R plots into staging, preserving the interactive device.
            grDevices::pdf(file.path(staging_root, site, "figures", "automatic_plots.pdf"))
            plot_device <- grDevices::dev.cur()
            tryCatch(sys.source(file.path(.pipeline_root, pipeline_stages[[stage]]), envir = env),
                     finally = { if (plot_device %in% grDevices::dev.list()) grDevices::dev.off(plot_device) })
            results[[key]] <- pipeline_publish(file.path(staging_root, site), output,
                                               pipeline_patterns[[stage]], window, reprocess)
          }, finally = unlink(staging_root, recursive = TRUE))
        }
      }
    }, finally = { if (!dry_run) unlink(lock, recursive = TRUE) })
  }
  invisible(results)
}
