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

# Keep each physical-variable family together; positional indices sort numerically.
# Names and values never change, only their published order.
amf_column_order <- function(columns, processed = character()) {
  vars <- setdiff(columns, c("TIMESTAMP_START", "TIMESTAMP_END"))
  family <- sub("_[0-9]+_[0-9]+_[0-9]+$", "", vars)
  family <- sub("(_PI)?_(F|CORR)$|_PI$", "", family)
  family <- sub("_(SIGMA|SSITC_TEST|MAX)$", "", family)
  family[grepl("^FETCH_", family)] <- "FETCH"
  product <- ifelse(grepl("_CORR($|_)", vars),3L,ifelse(grepl("_F($|_)",vars),2L,ifelse(grepl("_PI($|_)",vars),1L,0L)))
  tower <- vars %in% processed | product > 0L
  core <- sub("_[0-9]+_[0-9]+_[0-9]+$", "", vars)
  indices <- t(vapply(vars,function(v) {
    if(!grepl("_[0-9]+_[0-9]+_[0-9]+$",v)) return(c(0,0,0))
    as.numeric(tail(strsplit(v,"_",fixed=TRUE)[[1]],3))
  },numeric(3)))
  c(intersect(c("TIMESTAMP_START","TIMESTAMP_END"),columns),
    vars[order(tower,match(family,unique(family)),product,core,indices[,1],indices[,2],indices[,3])])
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
  add("PPFD_IN", unit = "umolPhoton m-2 s-1", note = "L1 QAQC including recurring PPFD deficit screening against independent SW_IN; no gap filling or scaling of the submitted PPFD units")
  add("ALB", unit = "%")
  add("PA", unit = "kPa")
  add("P_Bulk", "P_1_1_1", "mm", "L1 bulk precipitation (Pluvio or corrected SG400); no L2 gap filling")
  add("P", "P_1_1_2", "mm", "L1 tipping-bucket precipitation; no L2 gap filling")
  add("P_PI_F", "P_PI_F_1_1_1", "mm", "L2 filled precipitation: selected source, regression filling and final zero fallback; source can be tipping bucket where bulk is unavailable")
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
  for (v in c("LE_SSITC_TEST", "H_SSITC_TEST", "FC_SSITC_TEST", "TAU_SSITC_TEST"))
    add(v, target = if (v %in% c("LE_SSITC_TEST", "H_SSITC_TEST")) paste0(v, "_1_1_1") else v, unit = "1")
  add("SG", "SG_PI", unit = "W m-2", note = "Soil heat storage above plates; G + SG gives surface calorimetric heat flux")
  add("G_PI", "G_PI", "W m-2", "PI calorimetric surface soil heat flux, includes soil storage; do not add SG again")
  for (s in grep("^G_plate_[0-9]+_[0-9]+_[0-9]+$", columns, value = TRUE))
    add(s, sub("^G_plate_", "G_", s), "W m-2",
        "Individual L1 QAQC plate-depth soil heat flux; no above-plate storage addition or MDS filling; indices retained from source sensor")
  # Explicit instrument mappings: model names belong in metadata, not headers.
  # Reserve stable replicate slots; never renumber according to missing columns.
  for (i in 1:3) {
    add(paste0("TS_CS65X_1_1_", i), paste0("TS_1_1_", i + 3L), "deg C",
        "L1 QAQC soil temperature from CS65x; original source recorded; no gap filling")
    add(paste0("cs65x_ec_1_1_", i), paste0("SECA_1_1_", i), "dS m-1",
        "L1 CS65x bulk soil apparent electrical conductivity; not pore-water conductivity")
  }
  add("WS_05103_mean", "WS_1_1_2", "m s-1", "L1 RM Young 05103 mean wind speed; no gap filling")
  add("WD_05103_mean", "WD_1_1_2", "Decimal degrees", "L1 RM Young 05103 mean wind direction")
  add("WD_05103_StdDev", "WD_SIGMA_1_1_2", "Decimal degrees", "L1 RM Young 05103 wind-direction standard deviation")
  add("WS_05103_5s_Max", "WS_MAX_1_1_2", "m s-1", "L1 RM Young 05103 maximum five-second wind speed within averaging interval")
  add("TA_710", "TA_1_1_4", "deg C", "L1 retained LI-710 air temperature; no gap filling")
  add("RH_710", "RH_1_1_4", "%", "L1 retained LI-710 relative humidity; no gap filling")
  add("PA_710", "PA_1_1_2", "kPa", "L1 retained LI-710 pressure; no gap filling")
  out$target[out$source %in% c("WS", "WD", "WS_MAX", "PA")] <-
    paste0(out$source[out$source %in% c("WS", "WD", "WS_MAX", "PA")], "_1_1_1")
  # ALB is explicitly calculated as percent in L1.
  out[out$source %in% columns, , drop = FALSE]
}

amf_prepare <- function(item, flux_qc = "despike", include_storage = TRUE, local_standard_time = TRUE,
                        start = NULL, end = NULL, unit_overrides = character(),
                        extra_mapping = NULL, include_pi = TRUE, site = NULL) {
  if (!is.logical(include_pi) || length(include_pi) != 1L || is.na(include_pi))
    stop("include_pi must be TRUE or FALSE.")
  if (!isTRUE(local_standard_time)) stop("Confirm local standard time (no DST) and interval-end timestamps first.")
  flux_qc <- match.arg(flux_qc, "despike")
  if (!is.null(site) && (!is.character(site) || length(site) != 1L || is.na(site))) stop("site must be NULL or one site ID.")
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
  # Apply ONLY the PI-requested despike/physical-range masks to the
  # retained raw turbulent flux. filtered1 already includes unwanted extra masks.
  for (v in c("LE", "H", "FC")) {
    src <- v
    if (!src %in% names(x)) stop("Required raw turbulent flux column missing: ", src)
    to <- if (v == "FC") "umolCO2 m-2 s-1" else "W m-2"
    values <- amf_convert(x[[src]], input_unit(src), to, src)
    for (flag in paste0(v, "_QC_despike")) {
      if (!flag %in% names(x)) stop("Required QC column missing: ", flag)
      if (!is.numeric(x[[flag]]) && !is.logical(x[[flag]])) stop("Nonnumeric QC flag: ", flag)
      # Unassessed NA flags do not imply rejection, matching existing despike
      # handling. They are counted below so missing QC is visible to the PI.
      values[!is.na(x[[flag]]) & x[[flag]] >= 1] <- NA_real_
    }
    target <- if (v %in% c("LE", "H")) paste0(v, "_1_1_1") else v
    result[[target]] <- values
    record(src, target, input_unit(src), to, "exported", paste(flux_qc,
      "only; no fetch, SSITC, signal-strength, long-run or USTAR mask; no storage addition/filling/EBR; unassessed despike flags:",
      sum(is.na(x[[paste0(v, "_QC_despike")]]))))
  }
  # AmeriFlux FETCH_FILTER has the reverse polarity of internal *_QC_fetch:
  # 0 discard, 1 keep. Supply only at sites where footprint screening applies.
  if (!is.null(site) && site %in% c("EDVG", "ERVA")) {
    needed <- c("LE_QC_fetch", "H_QC_fetch", "FC_QC_fetch", "WD", "FETCH_70")
    if (!all(needed %in% names(x))) stop("Missing fetch metadata for ", site, ": ", paste(setdiff(needed,names(x)),collapse=", "))
    flag <- x$LE_QC_fetch
    if (any(!is.na(flag) & !flag %in% c(0,1,2))) stop("Unexpected internal fetch flag values for ", site)
    if (!identical(flag,x$H_QC_fetch) || !identical(flag,x$FC_QC_fetch)) stop("Flux fetch flags disagree for ", site)
    flag <- if ("FETCH_FILTER" %in% names(x)) x$FETCH_FILTER else rep(NA_real_,nrow(x))
    if (any(!is.na(flag) & !flag %in% c(0,1))) stop("Unexpected FETCH_FILTER values for ",site)
    # Historical L3 coerces unassessed fetch flags to zero. Missing footprint
    # or direction must not therefore be presented as an assessed pass.
    flag[!is.finite(x$WD) | !is.finite(x$FETCH_70) | x$FETCH_70 < 0 | x$WD < 0 | x$WD > 360] <- NA_real_
    result$FETCH_FILTER <- flag
    record("FETCH_FILTER", "FETCH_FILTER", "1", "1", "exported",
      "EDVG/ERVA footprint test; 0=discard, 1=keep, missing=unassessed. Internal H/FC fetch flags agree. Flag only: base flux is not fetch-filtered; final PI products retain existing fetch QC.")
  }
  mapping <- amf_mapping(names(x))
  # LI710 uses the existing unfilled L4 result of L3 diagnostic, physical-range/
  # despike, long-run and flow QC. Never substitute LE_710_PI_F or MDS estimates.
  for (v in c("LE", "H")) {
    src <- paste0(v, "_710_filtered")
    if (src %in% names(x)) mapping <- rbind(mapping, data.frame(
      source = src, target = paste0(v, "_1_1_2"), unit = "W m-2",
      note = "LI710 unfilled flux: L3 diag/despike/physical-range/long-run/flow flags applied in L4; no MDS or storage addition"))
    filled <- paste0(v, "_710_PI_F")
    if (filled %in% names(x)) mapping <- rbind(mapping, data.frame(
      source = filled, target = paste0(v, "_PI_F_1_1_2"), unit = "W m-2",
      note = "LI710 L4 REddyProc MDS-filled series after diagnostic/despike/physical-range/long-run/flow QC; no storage or EBR correction"))
  }
  if (!include_pi) mapping <- mapping[mapping$source != "G_PI", , drop = FALSE]
  # Retain existing EC MDS products with explicit storage provenance. Data Team
  # correspondence permits submission of tower-filled/corrected products with comments.
  # include_pi controls the OTHER local PI products, not these two columns.
  # FC_PI_F remains in L4 but is not submitted: its filling does not use the
  # standard ustar-filtered carbon-flux workflow. Retain observed FC and SC.
  ec_filled <- c("LE_PI_F", "H_PI_F")
  # Keep filled and corrected estimates distinct from unfilled turbulent flux.
  pi_sources <- c("LE_PI_F", "H_PI_F", "NETRAD_PI_F", 
                  "TA_PI_F", "RH_PI_F", "PA_PI_F", "VPD_PI_F", "SW_IN_PI_F", "LW_IN_PI_F", "P_PI_F",
                  "LE_PI_CORR", "H_PI_CORR", "LE_710_PI_CORR", "H_710_PI_CORR")
  pi_units <- c("W m-2", "W m-2", "W m-2",
                "deg C", "%", "kPa", "hPa", "W m-2", "W m-2", "mm", "W m-2", "W m-2", "W m-2", "W m-2")
  pi_targets <- c("LE_PI_F_1_1_1", "H_PI_F_1_1_1", "NETRAD_PI_F",
                  "TA_PI_F", "RH_PI_F", "PA_PI_F", "VPD_PI_F", "SW_IN_PI_F", "LW_IN_PI_F", "P_PI_F_1_1_1",
                  "LE_PI_CORR_1_1_1", "H_PI_CORR_1_1_1", "LE_PI_CORR_1_1_2", "H_PI_CORR_1_1_2")
  # P_PI_F already has an explicit PI-labeled mapping above.
  pi_sources_keep <- pi_sources != "P_PI_F"
  for (j in which(pi_sources_keep & (include_pi | pi_sources %in% ec_filled) & pi_sources %in% names(x))) mapping <- rbind(mapping,
    data.frame(source = pi_sources[j], target = pi_targets[j], unit = pi_units[j],
      note = if (grepl("710_PI_CORR", pi_sources[j]))
        "L4 LI710 QAQC + MDS + energy-balance correction; no EC storage added; custom CORR label documented for submission"
      else if (pi_sources[j] %in% c("LE_PI_F", "H_PI_F", "LE_PI_CORR", "H_PI_CORR"))
        "L4 EC full QAQC + available storage (turbulent fallback when storage absent) + MDS; CORR additionally applies EBR correction. Do not add storage again."
      else if (pi_sources[j] == "NETRAD_PI_F")
        "L4 MDS-filled net radiation; separate from retained NETRAD observations"
      else "L2 variable-specific filled meteorology; separate from L1 QAQC observations; VPD derives from filled TA/RH"))
  if (include_pi && "PPFD_IN_PI_F" %in% names(x)) mapping <- rbind(mapping,
    data.frame(source = "PPFD_IN_PI_F", target = "PPFD_IN_PI_F", unit = "umolPhoton m-2 s-1",
      note = "L2 PPFD: short interpolation, same-site observed SW regression, other-site PPFD regression, then local mean diurnal filling; photon flux units, not SW-scaled irradiance"))
  if (include_storage) {
    for (v in intersect(c("SLE", "SH", "SC"), names(x))) {
      mapping <- rbind(mapping, data.frame(source = v, target = if (v %in% c("SLE", "SH")) paste0(v, "_1_1_1") else v,
        unit = if (v == "SC") "umolCO2 m-2 s-1" else "W m-2",
        note = "Separate storage estimate retained from Level 4; document EddyPro storage method"))
    }
  }
  if (!is.null(extra_mapping)) {
    required <- c("source", "target", "unit", "note")
    if (!all(required %in% names(extra_mapping))) stop("extra_mapping needs source, target, unit, note.")
    mapping <- rbind(mapping[!mapping$target %in% extra_mapping$target, ], extra_mapping[, required])
  }
  # Allowlist only: no automatic passthrough of intermediate L4 variables.
  # Internal QC, MDS predictions, VBR/RE, instrument diagnostics and soil-method
  # components remain in L4 but are not submission columns.
  if (anyDuplicated(mapping$target) || any(mapping$target %in% names(result))) stop("Duplicate export variable.")
  if (any(!grepl("^[A-Z][A-Z0-9_]*$", mapping$target))) stop("Invalid export variable name.")
  for (i in seq_len(nrow(mapping))) {
    m <- mapping[i, ]; from <- input_unit(m$source)
    if (grepl("_SSITC_TEST$", m$source)) from <- "1" # L3 explicitly normalizes the 0/1/2 test scale.
    # Verified code provenance: ALB is recomputed as SW_OUT/SW_IN*100 in L1.
    if (m$source == "ALB" && !m$source %in% names(unit_overrides)) from <- "%"
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
    "Not exported: text provenance, optional product disabled, or requires explicit mapping")
  audit$source_stage <- vapply(audit$source, function(v) {
    if (grepl("_MDS|_filtered|_strg_added|_CORR|^ET_|^VBR|^LE_RE|^EBR|^(LE|H|FC).*_PI_F",v)) return("L4")
    if (grepl("_PI_F|^G_PI$|^G_(plate$|plate_sun|plate_shade|sun|shade|calorimetric|leuning)|^SG",v)) return("L2 (L4 MDS for NETRAD_PI_F/G_PI_F)")
    if (grepl("_QC_|_SSITC_TEST|^LE$|^H$|^FC$|^TAU|^SLE$|^SH$|^SC$|rand_err|^FETCH|^USTAR|^ZL$|^MO_LENGTH|^T_SONIC|^CO2$|^H2O$|SIGMA|^WS$|^WD$|^WS_MAX$|^(LE|H)_710|^(diag|flow|tilt|data_qc)_710",v)) return("L3; storage additionally screened in L4")
    "L1; retained through L2-L4"
  }, character(1))
  audit$source_scripts <- vapply(seq_len(nrow(audit)), function(i) {
    stage <- audit$source_stage[i]; src <- audit$source[i]
    if (src == "PPFD_IN") return("R/Level_1_met_QAQC_new.R; R/pipeline_ppfd_shadow.R")
    if (stage == "L4") return("R/Level_4_gapfilling_ET_calculation_new.R (retained L1-L3 inputs)")
    if (startsWith(stage,"L2")) return("R/Level_2_met_PI_vars_new.R; R/Level_4_gapfilling_ET_calculation_new.R for energy-driver MDS")
    if (startsWith(stage,"L3")) return(if (grepl("710",src)) "R/Level_3_LI710_QAQC_new.R; L4 flag application" else "R/Level_3_EC_data_QAQC_new.R; L4 storage screening; L5 base-flux mask")
    "R/Level_1_met_QAQC_new.R"
  }, character(1))
  # _PI is assigned by the network, not by the submitting tower team. Keep
  # original L4 source names in the audit, and move all tower products to the end.
  strip_pi <- function(v) gsub("_PI(?=_|$)", "", v, perl = TRUE)
  processed <- strip_pi(names(result)[grepl("_PI($|_)", names(result))])
  export_names <- strip_pi(names(result))
  if (anyDuplicated(export_names)) stop("Duplicate export variable after removing _PI.")
  names(result) <- export_names
  audit$target <- strip_pi(audit$target)
  result <- result[,amf_column_order(names(result), processed),drop=FALSE]
  # Mapping and coverage tables use the exact same order as the CSV/HTML.
  audit <- audit[order(ifelse(audit$status=="exported",0L,1L),
    match(audit$target,names(result)),seq_len(nrow(audit))),,drop=FALSE]
  rownames(audit) <- NULL
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
  comments <- c("AMERIFLUX SUBMISSION COMMENTS",
    "Curated observations plus explicitly requested tower-team final products; see metadata/variable_mapping.csv for every source and unit.",
    "The network assigns _PI; submitted headers omit it. _F denotes tower gap filling. CORR is a custom tower EBC label for Data Team review, not a standardized qualifier. Tower products follow all observation columns.",
    "TIMESTAMP_START/END: local standard time, no DST; 30-minute intervals. Use registered site ID in final upload filename.",
    "EC LE/H: _1_1_1; collocated equal-height LI710 LE/H: _1_1_2. Base EC applies despike/physical range only; fetch is a separate flag; no storage addition. Base LI710 uses L4 unfilled QC results.",
    "LE_F/H_F: existing final MDS-filled series. EC filled products include available storage; LI710 does not. Never add storage again.",
    "FC_F is excluded because L4 FC_PI_F does not use the standard ustar-filtered carbon-flux filling workflow. Observed FC and separate SC remain available.",
    "LE_CORR/H_CORR: existing filled and energy-balance-corrected series, separately indexed for EC and LI710. CORR denotes tower correction, not ONEFlux processing.",
    "G_i_j_k: individual plates. G: unfilled selected L2 surface heat flux including soil storage. SG: weighted L2 above-plate storage. Never add SG to G again. No sun/shade or alternative method components exported.",
    "P_1_1_1: bulk; P_1_1_2: tipping bucket; P_F_1_1_1: L2 selected and filled precipitation, including zero fallback.",
    "Meteorological observations precede filling; *_F are final filled products. NETRAD_F includes L4 MDS filling.",
    "CS65x temperature uses TS_1_1_4 through TS_1_1_6; conductivity uses SECA_1_1_1 through SECA_1_1_3. LI710 TA/RH use replicate 4 and PA replicate 2. 05103 wind uses replicate 2. Confirm positional indices against site sensor metadata.",
    "Instrument source names are recorded in metadata, not used as export variable base names. No internal QC/method codes, raw LI710 flux duplicates, MDS predictions, storage-added intermediates, VBR/RE or signal-strength columns are exported.")
  writeLines(comments, file.path(staged, "README.txt"))
  writeLines(comments, file.path(staged, "submission_comments.txt"))
  # Round-trip verify timestamps, sentinel formatting, and a single header row.
  check <- read.csv(file, colClasses = "character", check.names = FALSE)
  if (!identical(names(check), names(d)) || nrow(check) != nrow(d) ||
      !identical(check$TIMESTAMP_START, d$TIMESTAMP_START) ||
      !identical(check$TIMESTAMP_END, d$TIMESTAMP_END)) stop("Export round-trip validation failed.")
  if (!file.rename(staged, destination)) stop("Could not publish export directory.")
  file.path(destination, paste0(stem, ".csv"))
}
