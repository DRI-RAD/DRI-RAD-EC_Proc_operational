# Shared ingestion, period planning, and versioned output for the new pipeline.
# Scientific QAQC and gap-filling calculations remain in the Level_*_new.R files.

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

pipeline_history <- function(path, pattern) {
  files <- pipeline_files(path, pattern)
  if (!length(files)) return(NULL)
  result <- NULL
  for (file in files) {
    item <- pipeline_read_table(file)
    if (!is.null(result)) {
      if (!setequal(names(result$data), names(item$data)))
        stop("Output schema changed in ", file, "; migrate the history explicitly.")
      result$data <- dplyr::bind_rows(result$data, item$data)
      result$data <- result$data[!duplicated(result$data$TIMESTAMP, fromLast = TRUE), ]
      result$units <- item$units
    } else result <- item
  }
  result$data <- result$data[order(result$data$TIMESTAMP), ]
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
  files <- list.files(site_dir, pattern = "^eddypro_.*full_output.*\\.csv$",
                      recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  files <- files[tolower(basename(dirname(files))) == "output"]
  if (!length(files)) stop("No EddyPro full_output CSV in month/output folders: ", site_dir)
  # Prefer the latest EddyPro run timestamp in filenames, then mtime and path.
  # Folder month names are not trusted as observation bounds: reruns may overlap.
  run <- sub(".*full_output_([0-9-]+T[0-9]+).*", "\\1", basename(files))
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
