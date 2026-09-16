# Publication is stage-specific: one current cumulative table. L4 retains all
# previous tables; L1-L3 retain one. Figures use period and cumulative locations.
pipeline_period_label <- function(window) {
  paste(format(range(window$pending), "%Y%m%dT%H%M", tz = "UTC"), collapse = "_to_")
}

# Summaries always describe pending timestamps, never the MDS training context.
# Each percentage includes a count and a documented denominator; unavailable
# measurements are NA, including percentages with a zero denominator.
pipeline_flux_diagnostics <- function(data) {
  out <- list(diagnostic_rows = nrow(data), diagnostic_scope = "pending timestamps")
  add <- function(name, mask, denominator) {
    count <- sum(mask, na.rm = TRUE)
    out[[paste0(name, "_n")]] <<- count
    out[[paste0(name, "_denominator_n")]] <<- denominator
    out[[paste0(name, "_pct")]] <<- if (denominator > 0) 100 * count / denominator else NA_real_
  }
  for (flux in intersect(c("LE", "H", "FC", "LE_710", "H_710"), names(data))) {
    raw <- is.finite(data[[flux]])
    add(paste0(flux, "_raw_observed"), raw, nrow(data))
    flags <- grep(paste0("^", flux, "_QC_"), names(data), value = TRUE)
    for (flag in flags) {
      value <- data[[flag]]
      add(paste0(flag, "_missing"), is.na(value), nrow(data))
      for (code in sort(unique(value[!is.na(value)])))
        add(paste0(flag, "_code_", code), !is.na(value) & value == code, sum(!is.na(value)))
    }
    filtered <- paste0(flux, if (grepl("710$", flux)) "_filtered" else "_filtered2")
    filled <- paste0(flux, "_PI_F")
    if (filtered %in% names(data)) {
      accepted <- is.finite(data[[filtered]])
      add(paste0(flux, "_qc_removed"), raw & !accepted, sum(raw))
      if (filled %in% names(data)) {
        gap <- !accepted
        filled_gap <- gap & is.finite(data[[filled]])
        add(paste0(flux, "_gapfilled"), filled_gap, nrow(data))
        add(paste0(flux, "_gapfill_success"), filled_gap, sum(gap))
        add(paste0(flux, "_remaining_gap"), !is.finite(data[[filled]]), nrow(data))
      }
    }
  }
  out
}

# Capture MDS quality directly from REddyProc before the scripts discard QC
# columns. FillAll predictions at observed timestamps must not count as filling.
pipeline_capture_mds <- function(input, result, mapping) {
  opt <- getOption("ec.pipeline")
  keep <- input$TIMESTAMP %in% opt$window$pending
  stats <- opt$diagnostics
  for (variable in names(mapping)) {
    flux <- mapping[[variable]]
    gap <- !is.finite(input[[variable]][keep])
    filled <- is.finite(result[[paste0(variable, "_f")]][keep]) & gap
    quality <- result[[paste0(variable, "_fqc")]][keep]
    if (length(quality) != sum(keep)) stop("Missing MDS quality column for ", variable)
    good <- filled & !is.na(quality) & quality == 1
    stats[[paste0(flux, "_good_gapfilled_n")]] <- sum(good)
    stats[[paste0(flux, "_good_gapfilled_denominator_n")]] <- sum(filled)
    stats[[paste0(flux, "_good_gapfilled_pct")]] <- if (sum(filled)) 100 * sum(good)/sum(filled) else NA_real_
  }
  opt$diagnostics <- stats
  options(ec.pipeline = opt)
  invisible(NULL)
}

pipeline_publish <- function(staging, output, pattern, window, reprocess) {
  produced <- pipeline_files(staging, pattern)
  if (length(produced) != 1L) stop("Expected one final stage output: ", pattern)
  fresh <- pipeline_read_table(produced)
  if (!setequal(window$pending, fresh$data$TIMESTAMP)) stop("Stage output does not match pending timestamps.")
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
  site <- sub(paste0("_", pattern, ".*"), "", basename(produced), ignore.case = TRUE)
  stage <- names(pipeline_patterns)[match(pattern, pipeline_patterns)]
  filename <- paste0(site, "_", pattern, "_", Sys.Date(), ".csv")
  target <- file.path(output, filename)
  archive <- file.path(output, "archive")
  figures <- file.path(output, "figures")
  period_folder <- file.path(figures, paste0("period_", pipeline_period_label(window)))
  for (folder in c(archive, figures, period_folder)) dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  # Never delete recovery files left by a forcibly interrupted publication.
  if (length(list.files(output, pattern = "^\\.publish-", all.files = TRUE)))
    stop("Interrupted publication found in ", output, ". Recover the .publish-* backup before retrying.")
  tx <- tempfile(".publish-", tmpdir = output)
  dir.create(tx)
  moved <- list(); installed <- character(); committed <- FALSE
  on.exit({
    if (!committed) {
      unlink(installed)
      restored <- vapply(rev(moved), function(x) file.rename(x$backup, x$original), logical(1))
      if (!all(restored)) {
        warning("Recovery incomplete; preserve and inspect ", tx)
      } else unlink(tx, recursive = TRUE)
    } else unlink(tx, recursive = TRUE)
  }, add = TRUE)
  stash <- function(paths) {
    for (path in unique(paths[file.exists(paths)])) {
      backup <- file.path(tx, paste0("backup-", length(moved) + 1L))
      entry <- list(original = path, backup = backup)
      # Persist the mapping before moving a file, for manual crash recovery.
      saveRDS(c(moved, list(entry)), file.path(tx, "recovery.rds"))
      if (!file.rename(path, backup)) stop("Cannot rotate ", path)
      moved[[length(moved) + 1L]] <<- entry
    }
  }
  install <- function(from, to) {
    if (file.exists(to)) stop("Publication destination exists: ", to)
    installed <<- c(installed, to)
    saveRDS(installed, file.path(tx, "installed.rds"))
    if (!file.copy(from, to)) stop("Cannot publish ", to)
  }
  prepared <- file.path(tx, "current.csv")
  pipeline_write_table(data, fresh$units, prepared, "Cumulative processed output; prior rows preserved outside the requested update.")
  audit <- data.frame(stage = pattern, output = filename,
    first_timestamp = min(data$TIMESTAMP), last_timestamp = max(data$TIMESTAMP),
    new_first = min(fresh$data$TIMESTAMP), new_last = max(fresh$data$TIMESTAMP),
    new_rows = nrow(fresh$data), total_rows = nrow(data), context_start = window$read_start,
    calculation_days = as.numeric(difftime(window$end, window$read_start, units = "days")) + 1/48,
    figure_period = "new", figure_period_label = pipeline_period_label(window),
    reprocess = reprocess)
  stats <- pipeline_flux_diagnostics(fresh$data)
  extra <- getOption("ec.pipeline")$diagnostics
  if (length(extra)) stats[names(extra)] <- extra
  for (name in names(stats)) audit[[name]] <- stats[[name]]
  readr::write_csv(audit, file.path(tx, "audit.csv"))
  # Finish rendering before touching any published CSVs or figures.
  if (stage == "L4") {
    pipeline_level4_figures(fresh$data, site, file.path(staging, "figures"))
    pipeline_level4_figures(data, site, file.path(tx, "full_figures"))
  }
  current <- pipeline_files(output, pattern)
  previous <- pipeline_files(archive, pattern)
  if (stage == "L4") {
    # Preserve every original L4 CSV and audit, including multiple legacy root
    # versions. A collision suffix keeps repeated same-day runs distinct.
    for (file in current) {
      dest <- file.path(archive, basename(file))
      if (file.exists(dest) || file.exists(paste0(dest, ".period.csv"))) {
        suffix <- basename(tempfile(paste0(format(Sys.time(), "%H%M%S"), "-")))
        dest <- file.path(archive, paste0(tools::file_path_sans_ext(basename(file)), "_", suffix, ".csv"))
      }
      install(file, dest)
      audit_source <- paste0(file, ".period.csv")
      if (file.exists(audit_source)) install(audit_source, paste0(dest, ".period.csv"))
    }
    stash(c(current, paste0(current, ".period.csv")))
  } else {
  # The first migration consolidates legacy partial snapshots into one previous
  # cumulative table. Subsequent rotations preserve the exact previous file.
  if (length(current)) {
    latest <- tail(current, 1)
    prior <- file.path(tx, "previous.csv")
    if (length(current) == 1L) {
      if (!file.copy(latest, prior)) stop("Cannot prepare previous output")
    } else pipeline_write_table(history$data, history$units, prior, "Consolidated legacy output before this update.")
    prior_audit <- paste0(latest, ".period.csv")
    if (length(current) == 1L && file.exists(prior_audit)) {
      if (!file.copy(prior_audit, file.path(tx, "previous.period.csv"))) stop("Cannot prepare previous audit")
    }
  }
  stash(c(current, paste0(current, ".period.csv"), previous, paste0(previous, ".period.csv")))
  if (length(current)) {
    install(prior, file.path(archive, basename(latest)))
    if (file.exists(file.path(tx, "previous.period.csv")))
      install(file.path(tx, "previous.period.csv"), file.path(archive, paste0(basename(latest), ".period.csv")))
  }
  }
  # Existing figures and legacy runs are never moved, archived or pruned.
  # Reprocessing replaces only matching figures in the same period directory.
  new_figures <- list.files(file.path(staging, "figures"), full.names = TRUE)
  for (file in new_figures) {
    dest <- file.path(period_folder, basename(file))
    stash(dest)
    install(file, dest)
  }
  full_figures <- list.files(file.path(tx, "full_figures"), full.names = TRUE)
  for (file in full_figures) {
    # A dedicated suffix protects old, manually organized diagnostic files.
    name <- paste0(tools::file_path_sans_ext(basename(file)), "_full_period.", tools::file_ext(file))
    dest <- file.path(figures, name)
    stash(dest)
    install(file, dest)
  }
  install(file.path(tx, "audit.csv"), paste0(target, ".period.csv"))
  install(prepared, target)
  committed <- TRUE


  message("Diagnostics: ", period_folder)
  message("Committed ", pattern, ": ", nrow(fresh$data), " rows -> ", target)
  target
}
