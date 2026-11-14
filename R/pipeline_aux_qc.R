# Upstream QC for auxiliary sensors. Preserve input values and internal reason
# flags in L1/L3/L4; L5 exports only the screened standard quantity.
# Flags are additive bits: 1 invalid/range, 2 constant run, 4 wind 5-IQR.
# NA flag means the input was already missing. No interpolation is performed.
pipeline_aux_qc <- function(data, units, stage = c("L1", "L3")) {
  stage <- match.arg(stage)
  key <- function(u) gsub("[[:space:]\\[\\]\\^*+]", "", tolower(u), perl = TRUE)
  screen <- function(v, bounds = NULL, scale = 1, offset = 0, run = NULL, wind = FALSE) {
    if (!v %in% names(data)) return(invisible(NULL))
    if (!v %in% names(units)) stop("Missing auxiliary QC unit: ", v)
    raw <- data[[v]]
    if (!is.numeric(raw) && !is.logical(raw)) stop("Nonnumeric auxiliary sensor: ", v)
    y <- as.numeric(raw); flags <- rep(0L, length(y))
    invalid <- !is.finite(y) | y %in% c(-9999, -6999)
    if (!is.null(bounds)) invalid <- invalid | (!is.na(y) &
      (y * scale + offset < bounds[1] | y * scale + offset > bounds[2]))
    flags[invalid & !is.na(raw)] <- 1L
    y[invalid] <- NA_real_
    if (!is.null(run) && length(y)) {
      # A missing half-hour or value breaks a run; only exact repeats count.
      dt <- c(NA_real_, diff(as.numeric(data$TIMESTAMP)))
      same <- c(FALSE, !is.na(y[-1]) & !is.na(y[-length(y)]) &
        y[-1] == y[-length(y)]) & !is.na(dt) & dt == 1800
      group <- cumsum(!same)
      n <- tabulate(group)[group]
      reject <- !is.na(y) & n >= run
      flags[reject] <- flags[reject] + 2L; y[reject] <- NA_real_
    }
    if (wind) {
      z <- y[is.finite(y)]
      # Match L1's WS_05103_mean Q1-5*IQR / Q3+5*IQR envelope.
      # A zero-IQR or tiny sample cannot define a useful outlier envelope.
      if (length(z) >= 4L) {
        q <- as.numeric(quantile(z, c(.25, .75))); spread <- q[2] - q[1]
        if (spread > 0) {
          reject <- !is.na(y) & (y < q[1] - 5 * spread | y > q[2] + 5 * spread)
          flags[reject] <- flags[reject] + 4L; y[reject] <- NA_real_
        }
      }
    }
    flags[is.na(raw)] <- NA_integer_
    data[[paste0(v, "_raw")]] <<- raw
    data[[v]] <<- y
    data[[paste0(v, "_QC_aux")]] <<- flags
    units[[paste0(v, "_raw")]] <<- units[[v]]
    units[[paste0(v, "_QC_aux")]] <<- c("1", rep("derived", nrow(units)-1L))
  }
  if (stage == "L1") {
    for (v in intersect(c("TA_710", "RH_710", "PA_710"), names(data))) {
      u <- key(as.character(units[[v]][1]))
      if (v == "TA_710") {
        if (!u %in% c("degc", "c", "k", "kelvin")) stop("Unknown TA_710 unit: ", u)
        screen(v, c(-60, 60), offset = if (u %in% c("k", "kelvin")) -273.15 else 0, run = 16L)
      } else if (v == "RH_710") {
        if (!u %in% c("%", "percent", "fraction")) stop("Unknown RH_710 unit: ", u)
        # Same 0--110% tolerance as the existing L1 RH sensors.
        screen(v, c(0, 110), scale = if (u == "fraction") 100 else 1, run = 16L)
      } else {
        factors <- c(kpa=1, pa=.001, hpa=.1, mbar=.1)
        if (!u %in% names(factors)) stop("Unknown PA_710 unit: ", u)
        screen(v, c(50, 110), scale = factors[[u]], run = 16L)
      }
    }
    screen("WS_05103_5s_Max", c(0, Inf), wind = TRUE)
  } else {
    screen("WS_MAX", c(0, Inf), wind = TRUE)
    # Minimal sanity checks only. No statistical/upper-range/long-run or flux
    # flags applied to gas means or standard deviations; preserve real extremes.
    for (v in c("CO2", "H2O", "CO2_SIGMA", "H2O_SIGMA")) screen(v, c(0, Inf))
  }
  list(data = data, units = units)
}
