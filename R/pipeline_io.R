# Shared ingestion, period planning, and versioned output for the new pipeline.
# Scientific QAQC and gap-filling calculations remain in the Level_*_new.R files.

# openeddy's internal desp_loop cannot handle an empty day/night subset and
# stops when too few second differences are available. Guard only those cases
# in a private function environment; never modify the installed package.
pipeline_despike <- function(x, var, ...) {
  implementation <- openeddy::despikeLF
  package_env <- environment(implementation)
  original_loop <- get("desp_loop", envir = package_env)
  skipped <- integer()
  guarded_loop <- function(SD_sub, date, nVals, z, c, plot = FALSE) {
    values <- SD_sub$var
    n <- length(values)
    differences <- if (n >= 3L)
      (values[2:(n-1)] - values[1:(n-2)]) - (values[3:n] - values[2:(n-1)]) else numeric()
    if (sum(!is.na(differences)) <= nVals) {
      skipped <<- union(skipped, SD_sub$Index)
      # Spike remains NA for unassessed observations. Already detected spikes
      # and physical-range flags are preserved by the outer openeddy function.
      if (plot) return(list(SD = SD_sub, plots = list()))
      return(SD_sub)
    }
    original_loop(SD_sub, date, nVals, z, c, plot)
  }
  local_env <- new.env(parent = package_env)
  local_env$desp_loop <- guarded_loop
  environment(implementation) <- local_env
  result <- implementation(x, var, ...)
  if (length(skipped)) warning(
    var, ": insufficient observations for statistical despiking in one or more day/night subsets (",
    length(skipped), " rows). Unassessed flags remain NA; physical-range flags are retained.",
    call. = FALSE)
  result
}

# Backward-compatible name for the LI710 adapter.
pipeline_li710_despike <- pipeline_despike

# Assess each existing storage channel independently on the actual input clock.
pipeline_storage_qc <- function(data) {
  rows <- seq_len(nrow(data))
  for (variable in intersect(c("SLE", "SH", "SC"), names(data))) {
    if (!any(is.finite(data[[variable]][rows]))) next
    x <- as.data.frame(data[rows, , drop = FALSE])
    x$timestamp <- x$TIMESTAMP
    if (any(!is.finite(x$PotRad))) stop(
      "Storage QAQC: PotRad missing on actual upstream rows (first: ",
      x$TIMESTAMP[which(!is.finite(x$PotRad))[1]],
      "). Check L1 potential radiation.", call. = FALSE)
    flag <- pipeline_despike(x, variable, iter = 10, light = "PotRad", z = 10)
    rejected <- !is.na(flag) & flag >= 1
    data[[variable]][rows[rejected]] <- NA_real_
  }
  data
}

pipeline_time <- function(x) {
  if (inherits(x, "POSIXt")) return(as.POSIXct(x, tz = "UTC"))
  # UTC is a fixed clock label, matching the original logger scripts. No DST shift.
  x <- as.character(x)
  x <- sub("(Z| UTC)$", "", sub("T", " ", x, fixed = TRUE))
  # EddyPro normally writes HH:MM; logger exports may include seconds.
  x[grepl("^\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}$", x)] <-
    paste0(x[grepl("^\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}$", x)], ":00")
  out <- as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  if (anyNA(out)) stop("Invalid timestamp. Use YYYY-MM-DD HH:MM:SS.")
  out
}

pipeline_path <- function(root, path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) return(NA_character_)
  if (grepl("^[A-Za-z]:[/\\\\]|^[/\\\\]{2}", path)) return(path)
  file.path(root, sub("^[/\\\\]+", "", path))
}

pipeline_read_table <- function(file) {
  # Preserve the existing four-row format: description, names, units, aggregation.
  units <- readr::read_csv(file, skip = 1, n_max = 2,
                          col_types = readr::cols(.default = "c"),
                          na = c("", "NA"), name_repair = "minimal")
  if (!"TIMESTAMP" %in% names(units) || anyDuplicated(names(units)))
    stop("Missing or duplicate columns in ", file)
  data <- readr::read_csv(file, skip = 4, col_names = names(units),
                         guess_max = 100000, show_col_types = FALSE,
                         na = c("", "NA", "NaN"), name_repair = "minimal")
  data$TIMESTAMP <- pipeline_time(data$TIMESTAMP)
  if (anyDuplicated(data$TIMESTAMP)) stop("Duplicate timestamps in ", file)
  list(data = data, units = units)
}

pipeline_files <- function(path, pattern) {
  files <- list.files(path, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  files <- files[grepl("\\.csv$", files, ignore.case = TRUE) & !grepl("\\.period\\.csv$", files)]
  # Merge historical snapshots as well as new snapshots; dates in filenames are
  # save dates, not observation coverage. Newest modified rows win on overlap.
  files[order(file.info(files)$mtime, files)]
}

pipeline_uncertainty_schema <- function(item, pattern) {
  # Narrow additive migration for the four EddyPro uncertainty fields only.
  # Old saved observations have unknown uncertainty: never reconstruct it.
  if (!tolower(pattern) %in% c("level_3_ec_with_qaqc", "level_4_post_processed_data")) return(item)
  for (v in c("LE", "H", "FC", "TAU")) {
    name <- paste0(v, "_rand_err")
    if (!name %in% names(item$data)) {
      item$data[[name]] <- rep(NA_real_, nrow(item$data))
      unit <- if (v == "FC") "umolCO2 m-2 s-1" else if (v == "TAU") "kg m-1 s-2" else "W m-2"
      item$units[[name]] <- c(unit, rep(NA_character_, nrow(item$units) - 1L))
    }
  }
  item
}

pipeline_le_quality_columns <- function() {
  c("LE_PI_F_QC", "LE_MDS_QC", "LE_710_PI_F_QC", "LE_710_MDS_QC",
    "LE_VBR_Ts_F_QC", "LE_VBR_Ts_MDS_QC", "LE_RE_F_QC", "LE_RE_MDS_QC")
}

pipeline_output_schema <- function(item, pattern) {
  item <- pipeline_uncertainty_schema(item, pattern)
  if (tolower(pattern) != "level_4_post_processed_data") return(item)
  # Retire only the alternative VBR method. Old source/archive files are not
  # rewritten on read. Never reinterpret alternative-method values as VBR-Ts.
  retired <- c("VBR", "LE_VBR", "LE_VBR_F", "LE_VBR_MDS", "ET_VBR", "ET_VBR_F")
  item$data <- item$data[, setdiff(names(item$data), retired), drop = FALSE]
  item$units <- item$units[, setdiff(names(item$units), retired), drop = FALSE]
  # Historical QC and discarded lower-quality predictions cannot be recovered
  # from saved L4 values. Preserve values; require L4 reprocessing for backfill.
  for (name in pipeline_le_quality_columns()) {
    if (!name %in% names(item$data)) {
      item$data[[name]] <- rep(NA_real_, nrow(item$data))
      item$units[[name]] <- c("1", rep("quality flag", nrow(item$units) - 1L))
    }
  }
  item
}

pipeline_history <- function(path, pattern) {
  files <- pipeline_files(path, pattern)
  if (!length(files)) return(NULL)
  result <- NULL
  for (file in files) {
    item <- pipeline_output_schema(pipeline_read_table(file), pattern)
    if (!is.null(result)) {
      if (!setequal(names(result$data), names(item$data)))
        stop("Output schema changed in ", file, "; migrate the history explicitly.")
      result$data <- dplyr::bind_rows(result$data, item$data)
      result$data <- result$data[!duplicated(result$data$TIMESTAMP, fromLast = TRUE), ]
      result$units <- item$units
    } else result <- item
  }
  result$data <- result$data[order(result$data$TIMESTAMP), ]
  result$units <- result$units[, names(result$data), drop = FALSE]
  result
}

pipeline_latest <- function(files) {
  # Used by the new stage scripts. The runner provides a merged input snapshot,
  # so even legacy partial exports are available as one consistent table.
  if (!length(files)) stop("Required upstream output is missing.")
  files[order(file.info(files)$mtime, files, decreasing = TRUE)][1]
}

pipeline_slice <- function(data, window = getOption("ec.pipeline")$window,
                           predecessor = FALSE) {
  if (!"TIMESTAMP" %in% names(data)) stop("Input has no TIMESTAMP column.")
  data$TIMESTAMP <- pipeline_time(data$TIMESTAMP)
  data <- data[order(data$TIMESTAMP), ]
  keep <- data$TIMESTAMP >= window$read_start & data$TIMESTAMP <= window$end
  # A cumulative precipitation counter needs its immediately preceding sample.
  if (predecessor && any(data$TIMESTAMP < window$read_start))
    keep[max(which(data$TIMESTAMP < window$read_start))] <- TRUE
  data[keep, , drop = FALSE]
}

pipeline_raw_times <- function(file) {
  if (is.na(file) || !file.exists(file)) stop("Missing raw input: ", file)
  header <- readr::read_csv(file, skip = 1, n_max = 0, show_col_types = FALSE)
  x <- readr::read_csv(file, skip = 4, col_names = names(header),
                       col_types = readr::cols_only(TIMESTAMP = "c"))
  if (!nrow(x)) stop("Empty raw input: ", file)
  sort(unique(pipeline_time(x$TIMESTAMP)))
}

pipeline_eddypro_file <- function(file) {
  # EddyPro full_output CSV: optional category row, field-name row, units row,
  # then data. Detect the header instead of assuming a fixed column count.
  lines <- readLines(file, n = 12L, warn = FALSE)
  header <- which(vapply(lines, function(line) {
    fields <- strsplit(gsub('"', '', line), ",", fixed = TRUE)[[1]]
    all(c("date", "time", "LE", "qc_LE") %in% trimws(fields))
  }, logical(1)))
  if (length(header) != 1L) stop("Cannot identify EddyPro header: ", file)
  units <- readr::read_csv(file, skip = header - 1L, n_max = 1L,
                          col_types = readr::cols(.default = "c"),
                          name_repair = "minimal")
  if (anyDuplicated(names(units))) stop("Duplicate EddyPro fields: ", file)
  data <- readr::read_csv(file, skip = header + 1L, col_names = names(units),
                         col_types = readr::cols(.default = "c"),
                         na = c("", "NA", "-9999", "NaN"), name_repair = "minimal")
  if (!nrow(data)) stop("Empty EddyPro file: ", file)
  data$TIMESTAMP <- pipeline_time(paste(data$date, data$time))
  numeric_fields <- setdiff(names(data), c("TIMESTAMP", "filename", "date", "time"))
  # Convert only fields consumed by the existing QAQC code. Flag strings and
  # optional text metadata elsewhere in the file need not be numeric.
  required <- c("LE", "qc_LE", "H", "qc_H", "co2_flux", "qc_co2_flux", "Tau", "qc_Tau",
                "co2_mole_fraction", "co2_molar_density", "co2_var",
                "h2o_mole_fraction", "h2o_molar_density", "h2o_var", "ts_var",
                "u_var", "v_var", "w_var", "sonic_temperature", "wind_speed",
                "max_wind_speed", "wind_dir", "u*", "(z-d)/L", "L",
                "x_90%", "x_70%", "x_50%", "x_30%", "x_10%")
  optional <- c("rand_err_LE", "rand_err_H", "rand_err_co2_flux", "rand_err_Tau",
                "LE_strg", "H_strg", "co2_strg")
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Missing EddyPro fields in ", file, ": ", paste(missing, collapse = ", "))
  for (name in intersect(c(required, optional), numeric_fields)) {
    value <- suppressWarnings(as.numeric(data[[name]]))
    if (any(!is.na(data[[name]]) & is.na(value))) stop("Non-numeric EddyPro field ", name, " in ", file)
    data[[name]] <- value
  }
  for (name in setdiff(optional, names(data))) {
    data[[name]] <- NA_real_
    units[[name]] <- NA_character_
  }
  units$TIMESTAMP <- "TS"
  list(data = data, units = units)
}

pipeline_eddypro <- function(site_dir, window = NULL) {
  if (is.na(site_dir) || !dir.exists(site_dir)) stop("Configure an existing EddyPro site folder: ", site_dir)
  # dir_eddypro is the complete site archive, including YYYY and YYYY_NN runs.
  # Match the product name, not a mandatory EddyPro prefix or output subfolder.
  files <- list.files(site_dir, pattern = "full_output.*\\.csv$",
                      recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) stop("No EddyPro full_output CSV below dir_eddypro: ", site_dir)
  # Prefer the latest EddyPro run timestamp in filenames, then mtime and path.
  # Folder years/run suffixes are not observation bounds. If a filename has no
  # standard run timestamp, use its modification time as the fallback priority.
  stamp <- regmatches(basename(files), regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{6}", basename(files)))
  run <- as.numeric(file.info(files)$mtime)
  has_stamp <- grepl("[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{6}", basename(files))
  parsed <- as.POSIXct(stamp, format = "%Y-%m-%dT%H%M%S", tz = "UTC")
  if (anyNA(parsed)) stop("Invalid EddyPro run timestamp in a filename below: ", site_dir)
  run[has_stamp] <- as.numeric(parsed)
  files <- files[order(run, file.info(files)$mtime, files)]
  result <- NULL
  for (file in files) {
    item <- pipeline_eddypro_file(file)
    if (!is.null(window)) item$data <- pipeline_slice(item$data, window)
    if (is.null(result)) result <- item else {
      common <- intersect(names(result$units), names(item$units))
      same <- is.na(result$units[1, common]) | is.na(item$units[1, common]) |
        result$units[1, common] == item$units[1, common]
      if (any(!same)) stop("Inconsistent EddyPro units: ", file)
      result$data <- dplyr::bind_rows(result$data, item$data)
    }
  }
  result$data <- result$data[!duplicated(result$data$TIMESTAMP, fromLast = TRUE), ]
  result$data <- result$data[order(result$data$TIMESTAMP), ]
  result$files <- files
  result
}

# The latest observation, not the filename/run date or latest nonmissing flux,
# defines the automatic cap. Never silently fall back to a longer logger clock.
pipeline_eddypro_end <- function(input, site = "site") {
  times <- input$data$TIMESTAMP
  if (is.null(times) || !length(times) || anyNA(times))
    stop("No valid EddyPro observation timestamps for ", site, "; specify end explicitly or supply EddyPro output.")
  max(pipeline_time(times))
}

pipeline_grid_to_end <- function(times, end) {
  times <- pipeline_time(times)
  if (!length(times) || min(times) > end) return(times[FALSE])
  seq(min(times), end, by = 1800)
}

pipeline_plan <- function(available, history, start = NULL, end = NULL,
                          context_days = 0, reprocess = FALSE) {
  available <- sort(unique(pipeline_time(available)))
  if (!length(available)) return(NULL)
  if (!is.null(start)) available <- available[available >= pipeline_time(start)]
  if (!is.null(end)) available <- available[available <= pipeline_time(end)]
  if (!reprocess && !is.null(history))
    available <- available[!available %in% history$data$TIMESTAMP]
  if (!length(available)) return(NULL)
  list(start = min(available), end = max(available), pending = available,
       read_start = min(available) - context_days * 86400)
}

# Expand Level 4 backwards to at least min_days of half-hour observations.
# Historical QAQC inputs provide training context; pending timestamps never change.
pipeline_mds_window <- function(window, available, min_days = 100) {
  if (is.null(window)) return(NULL)
  required_start <- window$end - min_days * 86400 + 1800
  window$read_start <- min(window$read_start, required_start)
  # Keep the requested half-hour start. Do not round to midnight or create
  # pre-deployment rows; any earlier MDS context must exist upstream.
  expected <- seq(window$read_start, window$end, by = 1800)
  missing <- expected[!expected %in% available]
  first <- min(available)
  if (length(missing)) stop(
    "Insufficient Level 4 context: need common L1/L2/L3_EC timestamps from ",
    window$read_start, " through ", window$end, ". Missing ", length(missing),
    " rows (first: ", missing[1], "). Process the missing upstream period first with ",
    "stages = c('L1', 'L2', 'L3_EC'), then rerun Level 4. ",
    "Available common coverage: ", first, " through ", max(available),
    "; requested minimum: ", min_days, " days. If the site was deployed later, ",
    "earlier observations do not exist; choose an appropriate min_mds_days (at least 90) or wait for more data.")
  window$calculation_days <- length(expected) / 48
  window
}

pipeline_validate_mds <- function(timestamps, min_days = getOption("ec.pipeline")$min_mds_days) {
  t <- pipeline_time(timestamps)
  if (length(t) < min_days * 48 || anyNA(t) ||
      any(as.numeric(t) %% 1800 != 0) || any(diff(as.numeric(t)) != 1800))
    stop("MDS input must contain at least ", min_days, " days on a regular half-hour grid.")
  invisible(TRUE)
}

# Use a separate plot table so diagnostics cannot alter the calculation table.
pipeline_figure_data <- function(data) {
  opt <- getOption("ec.pipeline")
  # Figures always use the saved period, excluding preceding calculation context.
  data[data$TIMESTAMP %in% opt$window$pending, , drop = FALSE]
}

pipeline_figure_path <- function(name) {
  folder <- getOption("ec.pipeline")$figure_dir
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  file.path(folder, name)
}

pipeline_write_table <- function(data, units, file, description) {
  if (!identical(names(data), names(units))) stop("Output data/units column mismatch: ", file)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  # Write to a sibling temporary file and rename only after a successful readback.
  tmp <- tempfile(".pending-", tmpdir = dirname(file))
  on.exit(unlink(tmp), add = TRUE)
  write.table(matrix(c(description, rep("", ncol(data) - 1L)), nrow = 1), tmp,
              sep = ",", row.names = FALSE, col.names = FALSE, quote = TRUE, na = "NA")
  write.table(matrix(names(data), nrow = 1), tmp, sep = ",", append = TRUE,
              row.names = FALSE, col.names = FALSE, quote = TRUE)
  write.table(units, tmp, sep = ",", append = TRUE, row.names = FALSE,
              col.names = FALSE, quote = TRUE, na = "NA")
  readr::write_csv(data, tmp, append = TRUE, col_names = FALSE, na = "NA")
  checked <- pipeline_read_table(tmp)
  if (nrow(checked$data) != nrow(data)) stop("Output readback failed: ", file)
  if (file.exists(file)) stop("Refusing to overwrite an existing version: ", file)
  if (!file.rename(tmp, file)) stop("Cannot publish output: ", file)
  invisible(file)
}

pipeline_write <- function(df, units_df, file, description = "Data description") {
  opt <- getOption("ec.pipeline")
  df$TIMESTAMP <- pipeline_time(df$TIMESTAMP)
  df <- df[df$TIMESTAMP %in% opt$window$pending, , drop = FALSE]
  if (!nrow(df)) stop("The stage produced no rows in the pending period.")
  if (anyDuplicated(df$TIMESTAMP)) stop("Stage produced duplicate timestamps.")
  # Each stage writes to a staging directory. Publication happens only when the
  # whole script completes, so an intermediate EC export cannot advance state.
  pipeline_write_table(df, units_df, file, description)
}
