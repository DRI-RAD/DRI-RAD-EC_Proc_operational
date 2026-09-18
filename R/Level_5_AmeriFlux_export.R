# AmeriFlux BASE-In export helpers. These functions never modify Level 1-4 data.
# Source through run_level5_ameriflux.R. UTC is only the pipeline's fixed-clock
# container: local-standard timestamps are preserved without a timezone shift.

amf_unit_key <- function(x) {
  x <- tolower(trimws(enc2utf8(as.character(x))))
  # R may serialize the micro symbol this way under a non-Unicode locale.
  x <- gsub("<u+00b5>", "u", x, fixed = TRUE)
  x <- gsub("<u+03bc>", "u", x, fixed = TRUE)
  x <- gsub("\u00b5|\u03bc", "u", x, perl = TRUE)
  x <- gsub("\u00b0", "deg", x)
  x <- gsub("+1", "", x, fixed = TRUE)
  gsub("[[:space:]\\[\\]\\^*+]", "", x, perl = TRUE)
}

amf_convert <- function(x, from, to, variable, difference = FALSE) {
  # Refuse unknown unit conversions; never infer units from the value magnitude.
  if (!is.numeric(x) && !is.logical(x)) stop(variable, ": nonnumeric observations.")
  x <- as.numeric(x)
  x[!is.finite(x) | x %in% c(-9999, -6999)] <- NA_real_
  a <- amf_unit_key(from); b <- amf_unit_key(to)
  if (length(a) != 1L || is.na(a)) stop(variable, ": missing input unit.")
  aliases <- list(
    "degc" = c("c", "degc", "celsius"),
    "wm-2" = c("wm-2", "w/m2", "w/m-2"),
    "ms-1" = c("ms-1", "m/s"),
    "umolco2m-2s-1" = c("umolco2m-2s-1", "umolm-2s-1", "umols-1m-2", "umol/m2/s"),
    "umolco2mol-1" = c("umolco2mol-1", "umolmol-1", "umol/mol", "ppm"),
    "mmolh2omol-1" = c("mmolh2omol-1", "mmolmol-1", "mmol/mol"),
    "umolphotonm-2s-1" = c("umolphotonm-2s-1", "umolm-2s-1", "umols-1m-2", "umol/m2/s"),
    "decimaldegrees" = c("decimaldegrees", "deg", "degrees", "degree"),
    "1" = c("1", "-", "#", "adimensional", "dimensionless", "unitless"),
    "kgm-1s-2" = c("kgm-1s-2", "kgm-1s-2", "nm-2", "n/m2", "pa"))
  if (identical(a, b) || a %in% aliases[[b]]) return(x)
  if (b == "degc" && a %in% c("k", "kelvin"))
    return(if (difference) x else x - 273.15)
  if (b == "kpa" && a %in% c("pa", "hpa", "mbar"))
    return(x * if (a == "pa") 0.001 else 0.1)
  if (b == "hpa" && a %in% c("pa", "kpa", "mbar"))
    return(x * switch(a, pa = 0.01, kpa = 10, mbar = 1))
  if (b == "%" && a %in% c("fraction", "m3m-3", "m3/m3", "ratio")) return(x * 100)
  if (b == "%" && a %in% c("percent", "percentage")) return(x)
  if (b == "mm" && a %in% c("m", "cm")) return(x * if (a == "m") 1000 else 10)
  stop(variable, ": unsupported unit '", from, "' -> '", to, "'. Supply an explicit unit override after review.")
}

amf_mapping <- function(columns) {
  # source, destination, and target unit; only unambiguous existing names are
  # included automatically. Instrument-specific nonstandard names are reported.
  out <- data.frame(source = character(), target = character(), unit = character(),
                    note = character(), stringsAsFactors = FALSE)
  add <- function(source, target = source, unit, note = "Level 4 retained observation") {
    out <<- rbind(out, data.frame(source, target, unit, note, stringsAsFactors = FALSE))
  }
  for (v in c("NETRAD", "SW_IN", "SW_OUT", "LW_IN", "LW_OUT")) add(v, unit = "W m-2")
  add("PPFD_IN", unit = "umolPhoton m-2 s-1")
  add("ALB", unit = "%")
  add("PA", unit = "kPa")
  add("P", unit = "mm")
  for (v in c("TA", "RH", "TS", "SWC")) {
    sensors <- grep(paste0("^", v, "(_[1-9][0-9]*_[1-9][0-9]*_[1-9][0-9]*)?$"), columns, value = TRUE)
    for (s in sensors) add(s, unit = if (v %in% c("TA", "TS")) "deg C" else "%")
  }
  for (v in c("WS", "WS_MAX", "USTAR", "U_SIGMA", "V_SIGMA", "W_SIGMA")) add(v, unit = "m s-1")
  add("WD", unit = "Decimal degrees")
  add("ZL", unit = "1")
  add("MO_LENGTH", unit = "m")
  add("TAU", unit = "kg m-1 s-2")
  for (v in c("CO2", "CO2_SIGMA")) add(v, unit = "umolCO2 mol-1")
  for (v in c("H2O", "H2O_SIGMA")) add(v, unit = "mmolH2O mol-1")
  for (v in c("T_SONIC", "T_SONIC_SIGMA")) add(v, unit = "deg C")
  for (v in c("FETCH_90", "FETCH_70")) add(v, unit = "m")
  for (v in c("LE_SSITC_TEST", "H_SSITC_TEST", "FC_SSITC_TEST", "TAU_SSITC_TEST")) add(v, unit = "1")
  add("G_plate", "G", "W m-2", "Mean soil heat flux at plate depth; soil storage exported separately as SG")
  add("SG", unit = "W m-2", note = "Soil heat storage above plates; G + SG gives surface calorimetric heat flux")
  add("G_PI", "G_PI", "W m-2", "PI calorimetric surface soil heat flux, includes soil storage; do not add SG again")
  # ALB is explicitly calculated as percent in L1.
  out[out$source %in% columns, , drop = FALSE]
}

amf_prepare <- function(item, flux_qc = "despike_fetch", include_storage = TRUE, local_standard_time = TRUE,
                        start = NULL, end = NULL, unit_overrides = character(),
                        extra_mapping = NULL) {
  if (!isTRUE(local_standard_time)) stop("Confirm local standard time (no DST) and interval-end timestamps first.")
  flux_qc <- match.arg(flux_qc, "despike_fetch")
  if (!is.logical(include_storage) || length(include_storage) != 1L || is.na(include_storage))
    stop("include_storage must be TRUE or FALSE.")
  x <- as.data.frame(item$data)
  for (bound in list(start, end)) if (!is.null(bound) && length(bound) != 1L)
    stop("start/end must each be a single interval-end timestamp.")
  if (!is.null(start) && !is.null(end) && pipeline_time(start) > pipeline_time(end)) stop("start is after end.")
  if (!"TIMESTAMP" %in% names(x) || !nrow(x)) stop("No Level 4 observations.")
  x$TIMESTAMP <- pipeline_time(x$TIMESTAMP)
  if (anyDuplicated(x$TIMESTAMP)) stop("Duplicate timestamps in Level 4.")
  if (any(as.numeric(x$TIMESTAMP) %% 1800 != 0)) stop("Timestamps must lie on half-hour boundaries.")
  x <- x[order(x$TIMESTAMP), , drop = FALSE]
  if (!is.null(start)) x <- x[x$TIMESTAMP >= pipeline_time(start), , drop = FALSE]
  if (!is.null(end)) x <- x[x$TIMESTAMP <= pipeline_time(end), , drop = FALSE]
  if (!nrow(x)) stop("No observations in selected interval.")
  clock <- seq(min(x$TIMESTAMP), max(x$TIMESTAMP), by = 1800)
  x <- x[match(clock, x$TIMESTAMP), , drop = FALSE]
  x$TIMESTAMP <- clock
  result <- data.frame(TIMESTAMP_START = format(clock - 1800, "%Y%m%d%H%M", tz = "UTC"),
                       TIMESTAMP_END = format(clock, "%Y%m%d%H%M", tz = "UTC"))
  audit <- data.frame(source = character(), target = character(), input_unit = character(),
                     output_unit = character(), status = character(), note = character())
  record <- function(source, target, from, to, status, note) {
    audit <<- rbind(audit, data.frame(source, target, input_unit = from,
                                    output_unit = to, status, note))
  }
  input_unit <- function(v) {
    if (v %in% names(unit_overrides)) return(unit_overrides[[v]])
    if (!v %in% names(item$units)) return(NA_character_)
    as.character(item$units[[v]][1])
  }
  # Apply ONLY the PI-requested despike/physical-range and fetch masks to the
  # retained raw turbulent flux. filtered1 already includes unwanted extra masks.
  for (v in c("LE", "H", "FC")) {
    src <- v
    if (!src %in% names(x)) stop("Required raw turbulent flux column missing: ", src)
    to <- if (v == "FC") "umolCO2 m-2 s-1" else "W m-2"
    values <- amf_convert(x[[src]], input_unit(src), to, src)
    for (flag in paste0(v, c("_QC_despike", "_QC_fetch"))) {
      if (!flag %in% names(x)) stop("Required QC column missing: ", flag)
      if (!is.numeric(x[[flag]]) && !is.logical(x[[flag]])) stop("Nonnumeric QC flag: ", flag)
      # Unassessed NA flags do not imply rejection, matching existing despike
      # handling. They are counted below so missing QC is visible to the PI.
      values[!is.na(x[[flag]]) & x[[flag]] >= 1] <- NA_real_
    }
    result[[v]] <- values
    record(src, v, input_unit(src), to, "exported", paste(flux_qc,
      "only; no SSITC, signal-strength, long-run or USTAR mask; no storage addition/filling/EBR; unassessed despike flags:",
      sum(is.na(x[[paste0(v, "_QC_despike")]])), "; unassessed fetch flags:", sum(is.na(x[[paste0(v, "_QC_fetch")]]))))
  }
  mapping <- amf_mapping(names(x))
  # Preserve the PI's processed estimates separately in the local review file.
  # These labels require a later submission mapping; they are not BASE-In labels.
  pi_sources <- c("LE_PI_F", "H_PI_F", "FC_PI_F", "NETRAD_PI_F", "G_PI_F",
                  "TA_PI_F", "RH_PI_F", "PA_PI_F", "VPD_PI_F", "SW_IN_PI_F", "LW_IN_PI_F", "P_PI_F",
                  "LE_PI_CORR", "H_PI_CORR")
  pi_units <- c("W m-2", "W m-2", "umolCO2 m-2 s-1", "W m-2", "W m-2",
                "deg C", "%", "kPa", "hPa", "W m-2", "W m-2", "mm", "W m-2", "W m-2")
  for (j in which(pi_sources %in% names(x))) mapping <- rbind(mapping,
    data.frame(source = pi_sources[j], target = pi_sources[j], unit = pi_units[j],
      note = if (pi_sources[j] %in% c("LE_PI_F", "H_PI_F", "FC_PI_F", "LE_PI_CORR", "H_PI_CORR"))
        "PI processed result; includes Level 4 storage fallback and gap filling; CORR also includes EBR correction"
      else "PI processed/gap-filled result; separate from observations"))
  if (include_storage) {
    for (v in intersect(c("SLE", "SH", "SC"), names(x))) {
      mapping <- rbind(mapping, data.frame(source = v, target = v,
        unit = if (v == "SC") "umolCO2 m-2 s-1" else "W m-2",
        note = "Separate storage estimate retained from Level 4; document EddyPro storage method"))
    }
  }
  if (!is.null(extra_mapping)) {
    required <- c("source", "target", "unit", "note")
    if (!all(required %in% names(extra_mapping))) stop("extra_mapping needs source, target, unit, note.")
    mapping <- rbind(mapping[!mapping$target %in% extra_mapping$target, ], extra_mapping[, required])
  }
  if (anyDuplicated(mapping$target) || any(mapping$target %in% names(result))) stop("Duplicate export variable.")
  if (any(!grepl("^[A-Z][A-Z0-9_]*$", mapping$target))) stop("Invalid export variable name.")
  for (i in seq_len(nrow(mapping))) {
    m <- mapping[i, ]; from <- input_unit(m$source)
    if (grepl("_SSITC_TEST$", m$source)) from <- "1" # L3 explicitly normalizes the 0/1/2 test scale.
    # Verified code provenance: ALB is recomputed as SW_OUT/SW_IN*100 in L1.
    if (m$source == "ALB" && !m$source %in% names(unit_overrides)) from <- "%"
    # L2 copies PA to PA_PI_F before filling, but writes an incorrect '%' header.
    if (m$source == "PA_PI_F" && identical(from, "%") && !m$source %in% names(unit_overrides)) {
      from <- input_unit("PA")
      m$note <- paste(m$note, "Legacy percent header ignored; unit inherited from PA by L2 calculation provenance")
    }
    if (!m$source %in% names(x)) {
      record(m$source, m$target, from, m$unit, "omitted", "Source column absent")
      next
    }
    converted <- tryCatch({
      # L3 retains EddyPro sonic_temperature in K but inherits the EasyFlux
      # header for the mixed column. Convert by processing provenance, not size.
      if (m$source == "T_SONIC" && "processing" %in% names(x)) {
        a <- rep(NA_real_, nrow(x))
        ep <- !is.na(x$processing) & x$processing == "EddyPro"
        ef <- !is.na(x$processing) & x$processing == "EasyFlux"
        if (any(is.finite(x[[m$source]]) & !(ep | ef))) stop("Unknown sonic-temperature processing provenance")
        a[ep] <- amf_convert(x[[m$source]][ep], "K", m$unit, m$source)
        if (any(ef)) a[ef] <- amf_convert(x[[m$source]][ef], from, m$unit, m$source)
        a
      } else amf_convert(x[[m$source]], from, m$unit, m$source, m$source == "T_SONIC_SIGMA")
    }, error = identity)
    if (inherits(converted, "error")) {
      record(m$source, m$target, from, m$unit, "omitted", conditionMessage(converted))
    } else {
      result[[m$target]] <- converted
      note <- m$note
      if (m$source == "T_SONIC" && "processing" %in% names(x))
        note <- paste(note, "EddyPro rows K -> deg C; EasyFlux rows use source unit")
      record(m$source, m$target, from, m$unit, "exported", note)
    }
  }
  used <- unique(c("TIMESTAMP", audit$source))
  for (v in setdiff(names(x), used)) record(v, "", input_unit(v), "", "not_mapped",
    "Not automatically exported: derived/filled/diagnostic or requires instrument mapping")
  list(data = result, audit = audit,
       summary = data.frame(variable = names(result)[-(1:2)],
         observations = vapply(result[-(1:2)], function(v) sum(is.finite(v)), integer(1)),
         missing = vapply(result[-(1:2)], function(v) sum(!is.finite(v)), integer(1))))
}

amf_write <- function(prepared, site, output_dir) {
  if (length(site) != 1L || is.na(site) || !grepl("^[A-Za-z0-9_-]+$", site)) stop("Invalid internal site code.")
  d <- prepared$data
  stem <- paste0(site, "_Level_5_ameriflux_ready")
  # A separate export directory prevents Level 4 discovery/rotation touching L5.
  destination <- file.path(output_dir, paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_", site))
  if (dir.exists(destination)) destination <- tempfile(paste0(basename(destination), "_"), tmpdir = output_dir)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  staged <- tempfile(".ameriflux-", tmpdir = output_dir)
  dir.create(staged)
  on.exit(unlink(staged, recursive = TRUE), add = TRUE)
  file <- file.path(staged, paste0(stem, ".csv"))
  for (v in names(d)[-(1:2)]) d[[v]][!is.finite(d[[v]])] <- NA_real_
  write.table(d, file, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = FALSE, na = "-9999", eol = "\n", fileEncoding = "ASCII")
  dir.create(file.path(staged, "metadata"))
  write.csv(prepared$audit, file.path(staged, "metadata", "variable_mapping.csv"), row.names = FALSE, na = "")
  write.csv(prepared$summary, file.path(staged, "metadata", "coverage.csv"), row.names = FALSE)
  writeLines(c("LOCAL LEVEL 5 REVIEW FILE - NOT A FINAL AMERIFLUX UPLOAD",
    "Internal site name and _PI labels retained at PI request.",
    "Before upload: use the registered CC-xxx ID and BASE-In variable labels.",
    "TIMESTAMP_START/END: local standard time, no DST; 30-minute intervals.",
    "Main LE/H/FC: despike/physical range and fetch masks only; storage separate.",
    "PI flux products retain existing Level 4 processing, including storage and filling.",
    "Review metadata/variable_mapping.csv for omitted/unsupported variables.",
    "Only the data CSV is a data product; metadata files are not upload inputs."),
    file.path(staged, "README.txt"))
  # Round-trip verify timestamps, sentinel formatting, and a single header row.
  check <- read.csv(file, colClasses = "character", check.names = FALSE)
  if (!identical(names(check), names(d)) || nrow(check) != nrow(d) ||
      !identical(check$TIMESTAMP_START, d$TIMESTAMP_START) ||
      !identical(check$TIMESTAMP_END, d$TIMESTAMP_END)) stop("Export round-trip validation failed.")
  if (!file.rename(staged, destination)) stop("Could not publish export directory.")
  file.path(destination, paste0(stem, ".csv"))
}
