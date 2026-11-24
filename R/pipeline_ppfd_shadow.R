# SW_IN is the independent reference. Flag a repeated same-clock deficit, not
# individual cloudy observations. No irradiance correction is applied to SW_IN.
pipeline_ppfd_shadow <- function(data, min_days = 5L, min_sw = 100,
                                 deficit_fraction = .20, deficit_wm2 = 80,
                                 min_pairs = 100L, min_fit_days = 5L, context = NULL) {
  required <- c("TIMESTAMP", "SW_IN", "PPFD_IN")
  if (!all(required %in% names(data))) stop("PPFD shadow QC requires TIMESTAMP, SW_IN and PPFD_IN.")
  # Incremental runs may use preceding saved L1 observations for this check only.
  # Other L1 QC/calculation windows remain unchanged; return only current rows.
  if (!is.null(context) && nrow(data) && nrow(context)) {
    if (!all(required %in% names(context))) stop("Missing PPFD context fields.")
    context <- context[context$TIMESTAMP < min(data$TIMESTAMP) &
      context$TIMESTAMP >= min(data$TIMESTAMP)-14*86400,,drop=FALSE]
    if (nrow(context)) {
      if ("PPFD_IN_pre_shadow" %in% names(context)) {
        known <- !is.na(context$PPFD_IN_pre_shadow)
        context$PPFD_IN[known] <- context$PPFD_IN_pre_shadow[known]
      }
      all_data <- rbind(context[,required],data[,required])
      all_data <- all_data[order(all_data$TIMESTAMP),]
      result <- pipeline_ppfd_shadow(all_data,min_days,min_sw,deficit_fraction,deficit_wm2,min_pairs,min_fit_days)
      ix <- match(data$TIMESTAMP,result$data$TIMESTAMP)
      for(v in c("PPFD_IN","PPFD_IN_pre_shadow","PPFD_IN_SW_scaled","PPFD_IN_QC_shadow")) data[[v]]<-result$data[[v]][ix]
      result$data<-data;result$candidate<-result$candidate[ix];result$removed<-result$removed[ix]
      return(result)
    }
  }
  for (v in list(min_days, min_fit_days, min_pairs))
    if (length(v)!=1 || !is.numeric(v) || !is.finite(v) || v<1 || v!=floor(v)) stop("PPFD sample/day counts must be positive integers.")
  for (v in list(min_sw, deficit_wm2))
    if(length(v)!=1 || !is.numeric(v) || !is.finite(v) || v<0) stop("PPFD irradiance thresholds must be nonnegative.")
  if(length(deficit_fraction)!=1 || !is.numeric(deficit_fraction) || !is.finite(deficit_fraction) || deficit_fraction<=0 || deficit_fraction>=1) stop("deficit_fraction must lie between 0 and 1.")
  t <- data$TIMESTAMP; sw <- data$SW_IN; ppfd <- data$PPFD_IN
  if (!inherits(t,"POSIXt") || anyNA(t) || anyDuplicated(t) || is.unsorted(t) || any(as.numeric(t)%%1800!=0))
    stop("PPFD shadow QC requires sorted unique half-hour timestamps.")
  # The pipeline stores local standard clock values in a UTC container. Use
  # interval midpoints so midnight interval ends belong to the preceding day.
  mid <- t - 900
  day <- as.integer(as.Date(mid,tz="UTC")); slot <- format(mid,"%H:%M",tz="UTC")
  fit <- is.finite(sw) & sw >= min_sw & is.finite(ppfd) & ppfd > 0
  slope <- NA_real_; status <- "insufficient reference pairs"
  if (sum(fit)>=min_pairs && length(unique(day[fit]))>=min_fit_days) {
    # Robust preselection avoids letting a localized shadow drive the scale.
    # Then fit SW = slope * PPFD through the origin using reference inliers.
    ratio <- sw[fit]/ppfd[fit]; center <- median(ratio)
    tolerance <- max(3*mad(ratio), .20*center)
    idx <- which(fit); fit[idx[abs(ratio-center)>tolerance]] <- FALSE
    if (sum(fit)>=min_pairs && length(unique(day[fit]))>=min_fit_days) {
      slope <- sum(ppfd[fit]*sw[fit])/sum(ppfd[fit]^2)
      if (!is.finite(slope) || slope<=0) { slope<-NA_real_;status<-"invalid regression slope" } else status<-"assessed"
    }
  }
  scaled <- ppfd*slope; scaled[!is.finite(scaled)] <- NA_real_
  assessed <- is.finite(sw) & sw>=min_sw & is.finite(scaled) & ppfd>=0
  deficit <- sw-scaled
  candidate <- assessed & deficit>=deficit_wm2 & deficit>=deficit_fraction*sw
  removed <- rep(FALSE,nrow(data))
  for (clock in unique(slot)) {
    idx<-which(slot==clock)
    # Missing/invalid days, nondeficit days and missing clock records break runs.
    ok<-candidate[idx]
    group<-cumsum(c(TRUE,diff(day[idx])!=1 | !head(ok,-1) | !tail(ok,-1)))
    for(g in unique(group)) {
      rows<-idx[group==g]
      if(length(rows)>=min_days && all(candidate[rows])) removed[rows]<-TRUE
    }
  }
  data$PPFD_IN_pre_shadow <- ppfd
  data$PPFD_IN_SW_scaled <- scaled
  data$PPFD_IN_QC_shadow <- ifelse(assessed,as.numeric(removed),NA_real_)
  data$PPFD_IN[removed] <- NA_real_
  list(data=data, slope=slope, fit_pairs=sum(fit), status=status,
       candidate=candidate, removed=removed)
}
