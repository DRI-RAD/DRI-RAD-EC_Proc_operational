# Aggregate only complete local-standard days of interval-end observations.
# Keep missing dates as NA rows so plotted lines do not bridge absent days.
pipeline_figure_daily <- function(data, columns, fun = mean) {
  out <- data.frame(DATE = as.Date(character()), complete = logical())
  for (name in columns) out[[name]] <- numeric()
  if (!nrow(data)) return(out)
  dates <- as.Date(data$TIMESTAMP - 1800, tz = "UTC")
  days <- seq(min(dates), max(dates), by = "day")
  out <- data.frame(DATE = days, complete = FALSE)
  for (name in columns) out[[name]] <- NA_real_
  groups <- split(seq_len(nrow(data)), as.character(dates))
  for (i in seq_along(days)) {
    rows <- groups[[as.character(days[i])]]
    times <- as.numeric(data$TIMESTAMP[rows])
    clock_ok <- length(rows) == 48L && length(unique(times)) == 48L &&
      all(is.finite(times) & times %% 1800 == 0)
    valid <- clock_ok && all(columns %in% names(data)) &&
      all(vapply(data[rows, columns, drop = FALSE], function(x) all(is.finite(x)), logical(1)))
    if (!valid) next
    out$complete[i] <- TRUE
    for (name in columns) out[[name]][i] <- fun(data[[name]][rows])
  }
  out
}

# Each instrument has its own time extent. Filled and corrected ET must both
# exist to define a comparable ribbon. Interior missing dates remain gaps.
pipeline_water_daily <- function(data, et, corrected) {
  columns <- c(et, corrected, "P_PI_F")
  if (!all(columns %in% names(data))) return(data.frame())
  valid <- is.finite(data[[et]]) & is.finite(data[[corrected]])
  if (!any(valid)) return(data.frame())
  bounds <- range(data$TIMESTAMP[valid])
  data <- data[data$TIMESTAMP >= bounds[1] & data$TIMESTAMP <= bounds[2], , drop = FALSE]
  # ET ribbon uses paired complete days. P is also retained independently for
  # the daily time series, but cumulative P/ET use exactly the same dates.
  daily <- pipeline_figure_daily(data, c(et, corrected), sum)
  names(daily)[match(c(et, corrected), names(daily))] <- c("ET", "ET_corrected")
  p <- pipeline_figure_daily(data, "P_PI_F", sum)
  daily$P <- p$P_PI_F
  daily$ET_min <- pmin(daily$ET, daily$ET_corrected)
  daily$ET_max <- pmax(daily$ET, daily$ET_corrected)
  daily$paired <- daily$complete & p$complete
  # This is explicitly a complete-day subtotal, NOT an estimate across gaps.
  subtotal <- function(x) {
    value <- cumsum(ifelse(daily$paired, x, 0))
    value[!daily$paired] <- NA_real_
    value
  }
  daily$ET_cum <- subtotal(daily$ET)
  daily$ET_cum_corr <- subtotal(daily$ET_corrected)
  daily$P_cum <- subtotal(daily$P)
  daily$cum_min <- pmin(daily$ET_cum, daily$ET_cum_corr)
  daily$cum_max <- pmax(daily$ET_cum, daily$ET_cum_corr)
  # Explicit groups stop ribbons from spanning missing dates after NA removal.
  daily$et_segment <- cumsum(!daily$complete)
  daily$cum_segment <- cumsum(!daily$paired)
  daily$SSM <- NA_real_
  sensors <- intersect(c("SWC_1_1_1", "SWC_1_1_2"), names(data))
  if (length(sensors)) {
    dates <- as.Date(data$TIMESTAMP - 1800, tz = "UTC")
    groups <- split(seq_len(nrow(data)), as.character(dates))
    swc <- vapply(as.character(daily$DATE), function(day) {
      means <- vapply(data[groups[[day]], sensors, drop = FALSE], function(x) {
        x <- x[is.finite(x)]; if (length(x)) mean(x) else NA_real_
      }, numeric(1))
      if (any(is.finite(means))) mean(means[is.finite(means)]) else NA_real_
    }, numeric(1))
    if (any(is.finite(swc))) daily$SSM <- (swc - swc[which(is.finite(swc))[1]]) * 0.08 * 1000 / 100
  }
  daily
}

pipeline_water_figure <- function(data, site_id, device, et, corrected, file) {
  daily <- pipeline_water_daily(data, et, corrected)
  if (!nrow(daily) || !any(daily$complete | is.finite(daily$P))) {
    message <- if (!nrow(daily)) "No paired ET observations" else
      paste("No complete daily ET or precipitation:", paste(range(daily$DATE), collapse = " to "))
    plot <- ggplot() + theme_void() +
      annotate("text", x = 0, y = 0, label = paste(site_id, device, message, sep = "\n"))
  } else {
    period <- paste(format(range(daily$DATE), "%Y-%m-%d"), collapse = " to ")
    coverage <- paste0("Complete paired ET/P days: ", sum(daily$paired), "/", nrow(daily),
      ". Cumulative curves sum these days only; missing days are excluded.")
    offset <- 5; scale_factor <- 0.1
    a <- ggplot(daily, aes(DATE)) + theme_bw() +
      geom_line(aes(y = offset - P * scale_factor), color = "#0072B2", na.rm = TRUE) +
      geom_point(aes(y = offset - P * scale_factor), color = "#0072B2", size = 0.7, na.rm = TRUE) +
      geom_ribbon(aes(ymin = ET_min, ymax = ET_max, group = et_segment),
                  fill = "grey65", alpha = 0.5, na.rm = TRUE) +
      geom_line(aes(y = ET), color = "grey40", na.rm = TRUE) +
      geom_line(aes(y = ET_corrected), color = "black", na.rm = TRUE) +
      geom_point(aes(y = ET), color = "grey40", size = 0.7, na.rm = TRUE) +
      geom_point(aes(y = ET_corrected), color = "black", size = 0.7, na.rm = TRUE) +
      scale_y_continuous(name = "ET (mm/day)",
        sec.axis = sec_axis(~ (offset - .)/scale_factor, name = "Precipitation (mm/day)")) +
      labs(title = paste(site_id, device, "daily ET and precipitation"), subtitle = period,
        x = NULL, caption = "Gray: uncorrected ET; black: EBR-corrected ET. Daily totals require 48 valid intervals.")
    b <- ggplot(daily, aes(DATE)) + theme_bw() +
      geom_ribbon(aes(ymin = cum_min, ymax = cum_max, group = cum_segment),
                  fill = "grey65", alpha = 0.5, na.rm = TRUE) +
      geom_line(aes(y = ET_cum, color = "ET"), na.rm = TRUE) +
      geom_line(aes(y = ET_cum_corr, color = "ET corrected"), na.rm = TRUE) +
      geom_line(aes(y = P_cum, color = "P"), na.rm = TRUE) +
      geom_line(aes(y = SSM, color = "SSM"), na.rm = TRUE) +
      geom_point(aes(y = ET_cum, color = "ET"), size = 0.7, na.rm = TRUE) +
      geom_point(aes(y = ET_cum_corr, color = "ET corrected"), size = 0.7, na.rm = TRUE) +
      geom_point(aes(y = P_cum, color = "P"), size = 0.7, na.rm = TRUE) +
      geom_point(aes(y = SSM, color = "SSM"), size = 0.7, na.rm = TRUE) +
      scale_color_manual(values = c("ET" = "grey40", "ET corrected" = "black",
                                   "P" = "#0072B2", "SSM" = "#D55E00"), name = NULL) +
      labs(title = paste(site_id, device, "cumulative ET and precipitation"), subtitle = period,
        x = NULL, y = "Cumulative flux / storage change (mm)",
        caption = paste(coverage, "SSM: change from first available soil moisture (0.08 m layer).", sep = "\n"))
    # Label the last finite value of each series rather than an all-NA last row.
    for (name in c("ET_cum", "ET_cum_corr", "P_cum", "SSM")) {
      rows <- which(is.finite(daily[[name]]))
      if (!length(rows)) next
      label <- data.frame(DATE = daily$DATE[tail(rows, 1)], value = daily[[name]][tail(rows, 1)])
      b <- b + geom_text(data = label, aes(DATE, value, label = round(value, 1)),
                         inherit.aes = FALSE, hjust = -0.1, size = 3)
    }
    plot <- plot_grid(a, b, ncol = 1)
  }
  ggsave(file, plot = plot, width = 11, height = 7, dpi = 150)
  invisible(daily)
}

# Render six L4 diagnostics from either pending or cumulative rows.
# No QAQC or gap-filling is rerun for cumulative figures.
pipeline_level4_figures <- function(df_plot, site_id, figure_dir) {
  if (!nrow(df_plot)) return(invisible(NULL))
  # CSV readers infer entirely missing numeric columns as logical. Keep these
  # channels in the original numeric daily summaries, including EC-only runs.
  missing_numeric <- vapply(df_plot, function(x) is.logical(x) && all(is.na(x)), logical(1))
  df_plot[missing_numeric] <- lapply(df_plot[missing_numeric], as.numeric)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  pipeline_figure_path <- function(name) file.path(figure_dir, name)

for (name in setdiff(c("SWC_1_1_1", "SWC_1_1_2"), names(df_plot))) df_plot[[name]] <- NA_real_
# ec
ebr_val <- df_plot %>%
	mutate(
		RnG = NETRAD - G_PI,
		LEH = LE_strg_added + H_strg_added
	) %>%
	filter(!is.na(RnG) & !is.na(LEH)) %>%
	summarise(EBR_EC = sum(LEH) / sum(RnG)) %>%
	pull(EBR_EC) %>% round(3)

a <- df_plot %>%
	ggplot(aes(NETRAD - G_PI,LE_strg_added + H_strg_added)) +
	theme_bw() +
	geom_point(alpha = 0.2) +
	geom_abline(lty = 2) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "LE + H (W m-2)",
			 x = "Rn - G (W m-2)",
			 title = paste0(site_id,": Eddy covariance (half hourly)")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	annotate("text",x = -Inf, y = Inf,label = paste0("EBR = ", ebr_val),	hjust = -0.2,vjust = 4, size = 3.5)

b <- pipeline_figure_daily(df_plot, c("NETRAD_PI_F", "G_PI_F", "LE_PI_F", "H_PI_F")) %>%
	ggplot(aes(NETRAD_PI_F - G_PI_F,LE_PI_F + H_PI_F)) +
	theme_bw() +
	geom_point() +
	geom_abline(lty = 2) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "LE + H (W m-2)",
			 x = "Rn - G (W m-2)",
			 title = paste0(site_id,": Eddy covariance (daily)")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0)

c <- df_plot %>%
	mutate(
		HOUR = hour(TIMESTAMP) + minute(TIMESTAMP) / 60,
		`Rn - G` = NETRAD - G_PI,
		`LE + H` = LE_strg_added + H_strg_added
	) %>%
	filter(!is.na(`Rn - G`) & !is.na(`LE + H`)) %>%
	group_by(HOUR) %>%
	summarise(
		`Rn - G` = mean(`Rn - G`, na.rm = TRUE),
		`LE + H` = mean(`LE + H`, na.rm = TRUE),
		LE = mean(LE_strg_added, na.rm = TRUE),
		H = mean(H_strg_added, na.rm = TRUE),
		.groups = "drop"
	) %>%
	pivot_longer(
		cols = c(`Rn - G`, `LE + H`,LE,H),
		names_to = "var",
		values_to = "value"
	) %>%
	ggplot(aes(HOUR, value, color = var, shape = var)) +
	theme_bw() +
	geom_hline(yintercept = 0) +
	geom_line() +
	geom_point(size = 2) +
	scale_color_manual(
		name = NULL,
		breaks = c("Rn - G", "LE + H","LE", "H"),
		values = c(
			"Rn - G" = "#000000",        # black
			"LE + H" = "#999999",        # gray
			"LE" = "#0072B2",        # blue
			"H" = "#D55E00" # vermillion
		)
	) +
	scale_shape_manual(
		name = NULL,
		breaks = c("Rn - G", "LE + H","LE", "H"),
		values = c(
			"Rn - G" = 1,   # open circle
			"LE + H" = 16,  # filled circle
			"LE" = 18,
			"H" = 17  # triangle
		)
	) +
	labs(y = "Energy flux (W m-2)",
			 x = "",
			 title = paste0(site_id, ": EC mean diurnal cycle")
	) +
	theme(legend.position = "bottom",
				axis.title.x = element_blank())


d <- pipeline_figure_daily(df_plot, c("NETRAD_PI_F", "G_PI_F", "LE_PI_F", "H_PI_F")) %>%
	ggplot(aes(DATE)) +
	theme_bw() +
	geom_hline(yintercept = 0, lty = 2) +
	geom_line(aes(y = NETRAD_PI_F, color = "Rn")) +
	geom_line(aes(y = G_PI_F, color = "G")) +
	geom_line(aes(y = LE_PI_F, color = "LE")) +
	geom_line(aes(y = H_PI_F, color = "H")) +
	scale_color_manual(
		name = NULL,
		breaks = c("Rn", "G","LE", "H"),
		values = c(
			"Rn" = "#000000",        # black
			"G" = "#999999",        # gray
			"LE" = "#0072B2",        # blue
			"H" = "#D55E00" # vermillion
		)
	) +
	labs(y = "Energy flux (W m-2)",
			 x = "",
			 title = paste0(site_id, ": EC daily mean energy flux (EBR uncorrected)")
	)

p1 <- plot_grid(plot_grid(a,b,c,ncol = 3),d,ncol = 1)


# 710

ebr_val <- df_plot %>%
	mutate(
		RnG = NETRAD - G_PI,
		LEH = LE_710_filtered + H_710_filtered
	) %>%
	filter(!is.na(RnG) & !is.na(LEH)) %>%
	summarise(EBR_EC = sum(LEH) / sum(RnG)) %>%
	pull(EBR_EC) %>% round(3)

a <- df_plot %>%
	ggplot(aes(NETRAD - G_PI,LE_710_filtered + H_710_filtered)) +
	theme_bw() +
	geom_point(alpha = 0.2) +
	geom_abline(lty = 2) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "LE + H (W m-2)",
			 x = "Rn - G (W m-2)",
			 title = paste0(site_id,": LI710 (half hourly)")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	annotate("text",x = -Inf, y = Inf,label = paste0("EBR = ", ebr_val),	hjust = -0.2,vjust = 4, size = 3.5)


b <- pipeline_figure_daily(df_plot, c("NETRAD_PI_F", "G_PI_F", "LE_710_PI_F", "H_710_PI_F")) %>%
	ggplot(aes(NETRAD_PI_F - G_PI_F,LE_710_PI_F + H_710_PI_F)) +
	theme_bw() +
	geom_point() +
	geom_abline(lty = 2) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "LE + H (W m-2)",
			 x = "Rn - G (W m-2)",
			 title = paste0(site_id,": LI710 (daily)")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0)


c <- df_plot %>%
	mutate(
		HOUR = hour(TIMESTAMP) + minute(TIMESTAMP) / 60,
		`Rn - G` = NETRAD - G_PI,
		`LE + H` = LE_710_filtered + H_710_filtered
	) %>%
	filter(!is.na(`Rn - G`) & !is.na(`LE + H`)) %>%
	group_by(HOUR) %>%
	summarise(
		`Rn - G` = mean(`Rn - G`, na.rm = TRUE),
		`LE + H` = mean(`LE + H`, na.rm = TRUE),
		LE = mean(LE_710_filtered, na.rm = TRUE),
		H = mean(H_710_filtered, na.rm = TRUE),
		.groups = "drop"
	) %>%
	pivot_longer(
		cols = c(`Rn - G`, `LE + H`,LE,H),
		names_to = "var",
		values_to = "value"
	) %>%
	ggplot(aes(HOUR, value, color = var, shape = var)) +
	theme_bw() +
	geom_hline(yintercept = 0) +
	geom_line() +
	geom_point(size = 2) +
	scale_color_manual(
		name = NULL,
		breaks = c("Rn - G", "LE + H","LE", "H"),
		values = c(
			"Rn - G" = "#000000",        # black
			"LE + H" = "#999999",        # gray
			"LE" = "#0072B2",        # blue
			"H" = "#D55E00" # vermillion
		)
	) +
	scale_shape_manual(
		name = NULL,
		breaks = c("Rn - G", "LE + H","LE", "H"),
		values = c(
			"Rn - G" = 1,   # open circle
			"LE + H" = 16,  # filled circle
			"LE" = 18,
			"H" = 17  # triangle
		)
	) +
	labs(y = "Energy flux (W m-2)",
			 x = "",
			 title = paste0(site_id, ": 710 mean diurnal cycle")
	) +
	theme(legend.position = "bottom",
				axis.title.x = element_blank())


d <- pipeline_figure_daily(df_plot, c("NETRAD_PI_F", "G_PI_F", "LE_710_PI_F", "H_710_PI_F")) %>%
	ggplot(aes(DATE)) +
	theme_bw() +
	geom_hline(yintercept = 0, lty = 2) +
	geom_line(aes(y = NETRAD_PI_F, color = "Rn")) +
	geom_line(aes(y = G_PI_F, color = "G")) +
	geom_line(aes(y = LE_710_PI_F, color = "LE")) +
	geom_line(aes(y = H_710_PI_F, color = "H")) +
	scale_color_manual(
		name = NULL,
		breaks = c("Rn", "G","LE", "H"),
		values = c(
			"Rn" = "#000000",        # black
			"G" = "#999999",        # gray
			"LE" = "#0072B2",        # blue
			"H" = "#D55E00" # vermillion
		)
	) +
	labs(y = "Energy flux (W m-2)",
			 x = "",
			 title = paste0(site_id, ": 710 daily mean energy flux (EBR uncorrected)")
	)

p2 <- plot_grid(plot_grid(a,b,c,ncol = 3),d,ncol = 1)

ggsave(plot = plot_grid(p1,p2,ncol = 1),
			 filename = pipeline_figure_path(paste0(site_id,"_Level_4_Figure_1a_energy_balance.jpeg")),
			 width = 11, height = 13, dpi = 150
)



################################################################################
#### Diagnostic figure 1b: Energy balance (corrected)
################################################################################
# ec
ebr_val <- df_plot %>%
	mutate(
		RnG = NETRAD - G_PI,
		LEH = LE_PI_CORR + H_PI_CORR
	) %>%
	filter(!is.na(RnG) & !is.na(LEH)) %>%
	summarise(EBR_EC = sum(LEH) / sum(RnG)) %>%
	pull(EBR_EC) %>% round(3)

a <- df_plot %>%
	ggplot(aes(NETRAD - G_PI,LE_PI_CORR + H_PI_CORR)) +
	theme_bw() +
	geom_point(alpha = 0.2) +
	geom_abline(lty = 2) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "LE + H (W m-2)",
			 x = "Rn - G (W m-2)",
			 title = paste0(site_id,": corrected EC (half hourly)")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	annotate("text",x = -Inf, y = Inf,label = paste0("EBR = ", ebr_val),	hjust = -0.2,vjust = 4, size = 3.5)

b <- pipeline_figure_daily(df_plot, c("NETRAD_PI_F", "G_PI_F", "LE_PI_CORR", "H_PI_CORR")) %>%
	ggplot(aes(NETRAD_PI_F - G_PI_F,LE_PI_CORR + H_PI_CORR)) +
	theme_bw() +
	geom_point() +
	geom_abline(lty = 2) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "LE + H (W m-2)",
			 x = "Rn - G (W m-2)",
			 title = paste0(site_id,": corrected EC (daily)")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0)

c <- df_plot %>%
	mutate(
		HOUR = hour(TIMESTAMP) + minute(TIMESTAMP) / 60,
		`Rn - G` = NETRAD - G_PI,
		`LE + H` = LE_PI_CORR + H_PI_CORR
	) %>%
	filter(!is.na(`Rn - G`) & !is.na(`LE + H`)) %>%
	group_by(HOUR) %>%
	summarise(
		`Rn - G` = mean(`Rn - G`, na.rm = TRUE),
		`LE + H` = mean(`LE + H`, na.rm = TRUE),
		LE = mean(LE_PI_CORR, na.rm = TRUE),
		H = mean(H_PI_CORR, na.rm = TRUE),
		.groups = "drop"
	) %>%
	pivot_longer(
		cols = c(`Rn - G`, `LE + H`,LE,H),
		names_to = "var",
		values_to = "value"
	) %>%
	ggplot(aes(HOUR, value, color = var, shape = var)) +
	theme_bw() +
	geom_hline(yintercept = 0) +
	geom_line() +
	geom_point(size = 2) +
	scale_color_manual(
		name = NULL,
		breaks = c("Rn - G", "LE + H","LE", "H"),
		values = c(
			"Rn - G" = "#000000",        # black
			"LE + H" = "#999999",        # gray
			"LE" = "#0072B2",        # blue
			"H" = "#D55E00" # vermillion
		)
	) +
	scale_shape_manual(
		name = NULL,
		breaks = c("Rn - G", "LE + H","LE", "H"),
		values = c(
			"Rn - G" = 1,   # open circle
			"LE + H" = 16,  # filled circle
			"LE" = 18,
			"H" = 17  # triangle
		)
	) +
	labs(y = "Energy flux (W m-2)",
			 x = "",
			 title = paste0(site_id, ": corrected EC mdc")
	) +
	theme(legend.position = "bottom",
				axis.title.x = element_blank())


d <- pipeline_figure_daily(df_plot, c("NETRAD_PI_F", "G_PI_F", "LE_PI_CORR", "H_PI_CORR")) %>%
	ggplot(aes(DATE)) +
	theme_bw() +
	geom_hline(yintercept = 0, lty = 2) +
	geom_line(aes(y = NETRAD_PI_F, color = "Rn")) +
	geom_line(aes(y = G_PI_F, color = "G")) +
	geom_line(aes(y = LE_PI_CORR, color = "LE")) +
	geom_line(aes(y = H_PI_CORR, color = "H")) +
	scale_color_manual(
		name = NULL,
		breaks = c("Rn", "G","LE", "H"),
		values = c(
			"Rn" = "#000000",        # black
			"G" = "#999999",        # gray
			"LE" = "#0072B2",        # blue
			"H" = "#D55E00" # vermillion
		)
	) +
	labs(y = "Energy flux (W m-2)",
			 x = "",
			 title = paste0(site_id, ": EC daily mean energy flux (EBR corrected)")
	)

p1 <- plot_grid(plot_grid(a,b,c,ncol = 3),d,ncol = 1)


# 710

ebr_val <- df_plot %>%
	mutate(
		RnG = NETRAD - G_PI,
		LEH = LE_710_PI_CORR + H_710_PI_CORR
	) %>%
	filter(!is.na(RnG) & !is.na(LEH)) %>%
	summarise(EBR_EC = sum(LEH) / sum(RnG)) %>%
	pull(EBR_EC) %>% round(3)

a <- df_plot %>%
	ggplot(aes(NETRAD - G_PI,LE_710_PI_CORR + H_710_PI_CORR)) +
	theme_bw() +
	geom_point(alpha = 0.2) +
	geom_abline(lty = 2) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "LE + H (W m-2)",
			 x = "Rn - G (W m-2)",
			 title = paste0(site_id,": corrected LI710 (half hourly)")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	annotate("text",x = -Inf, y = Inf,label = paste0("EBR = ", ebr_val),	hjust = -0.2,vjust = 4, size = 3.5)


b <- pipeline_figure_daily(df_plot, c("NETRAD_PI_F", "G_PI_F", "LE_710_PI_CORR", "H_710_PI_CORR")) %>%
	ggplot(aes(NETRAD_PI_F - G_PI_F,LE_710_PI_CORR + H_710_PI_CORR)) +
	theme_bw() +
	geom_point() +
	geom_abline(lty = 2) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "LE + H (W m-2)",
			 x = "Rn - G (W m-2)",
			 title = paste0(site_id,": corrected LI710 (daily)")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0)


c <- df_plot %>%
	mutate(
		HOUR = hour(TIMESTAMP) + minute(TIMESTAMP) / 60,
		`Rn - G` = NETRAD - G_PI,
		`LE + H` = LE_710_PI_CORR + H_710_PI_CORR
	) %>%
	filter(!is.na(`Rn - G`) & !is.na(`LE + H`)) %>%
	group_by(HOUR) %>%
	summarise(
		`Rn - G` = mean(`Rn - G`, na.rm = TRUE),
		`LE + H` = mean(`LE + H`, na.rm = TRUE),
		LE = mean(LE_710_PI_CORR, na.rm = TRUE),
		H = mean(H_710_PI_CORR, na.rm = TRUE),
		.groups = "drop"
	) %>%
	pivot_longer(
		cols = c(`Rn - G`, `LE + H`,LE,H),
		names_to = "var",
		values_to = "value"
	) %>%
	ggplot(aes(HOUR, value, color = var, shape = var)) +
	theme_bw() +
	geom_hline(yintercept = 0) +
	geom_line() +
	geom_point(size = 2) +
	scale_color_manual(
		name = NULL,
		breaks = c("Rn - G", "LE + H","LE", "H"),
		values = c(
			"Rn - G" = "#000000",        # black
			"LE + H" = "#999999",        # gray
			"LE" = "#0072B2",        # blue
			"H" = "#D55E00" # vermillion
		)
	) +
	scale_shape_manual(
		name = NULL,
		breaks = c("Rn - G", "LE + H","LE", "H"),
		values = c(
			"Rn - G" = 1,   # open circle
			"LE + H" = 16,  # filled circle
			"LE" = 18,
			"H" = 17  # triangle
		)
	) +
	labs(y = "Energy flux (W m-2)",
			 x = "",
			 title = paste0(site_id, ": corrected 710 mdc")
	) +
	theme(legend.position = "bottom",
				axis.title.x = element_blank())


d <- pipeline_figure_daily(df_plot, c("NETRAD_PI_F", "G_PI_F", "LE_710_PI_CORR", "H_710_PI_CORR")) %>%
	ggplot(aes(DATE)) +
	theme_bw() +
	geom_hline(yintercept = 0, lty = 2) +
	geom_line(aes(y = NETRAD_PI_F, color = "Rn")) +
	geom_line(aes(y = G_PI_F, color = "G")) +
	geom_line(aes(y = LE_710_PI_CORR, color = "LE")) +
	geom_line(aes(y = H_710_PI_CORR, color = "H")) +
	scale_color_manual(
		name = NULL,
		breaks = c("Rn", "G","LE", "H"),
		values = c(
			"Rn" = "#000000",        # black
			"G" = "#999999",        # gray
			"LE" = "#0072B2",        # blue
			"H" = "#D55E00" # vermillion
		)
	) +
	labs(y = "Energy flux (W m-2)",
			 x = "",
			 title = paste0(site_id, ": 710 daily mean energy flux (EBR corrected)")
	)

p2 <- plot_grid(plot_grid(a,b,c,ncol = 3),d,ncol = 1)

ggsave(plot = plot_grid(p1,p2,ncol = 1),
			 filename = pipeline_figure_path(paste0(site_id,"_Level_4_Figure_1b_energy_balance_corrected.jpeg")),
			 width = 11, height = 13, dpi = 150
)

################################################################################
#### Diagnostic figure 2: water balance
################################################################################
# EC and LI-710 use independent data extents in both period and full-period output.
pipeline_water_figure(df_plot, site_id, "EC", "ET_EC_F", "ET_EC_EBC",
  pipeline_figure_path(paste0(site_id, "_Level_4_Figure_2a_water_balance_EC.jpeg")))
pipeline_water_figure(df_plot, site_id, "LI-710", "ET_710_F", "ET_710_EBC",
  pipeline_figure_path(paste0(site_id, "_Level_4_Figure_2b_water_balance_LI710.jpeg")))

################################################################################
#### Diagnostic figure 3: ET methods intercomparison
################################################################################
a <- df_plot %>%
	ggplot(aes(ET_EC,ET_710)) +
	theme_bw() +
	geom_abline(lty = 2) +
	geom_point(alpha = 0.2) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "710 ET (mm/30min)",
			 x = "EC ET (mm/30min)",
			 title = paste0(site_id,": EC vs 710")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)

b <- df_plot %>%
	filter(!is.na(ET_710)) %>%
	ggplot(aes(ET_EC,ET_710_EBC)) +
	theme_bw() +
	geom_abline(lty = 2) +
	geom_point(alpha = 0.2) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "corrected 710 ET (mm/30min)",
			 x = "EC ET (mm/30min)",
			 title = paste0(site_id,": EC vs EBC 710")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)

c <- df_plot %>%
	ggplot(aes(ET_EC,ET_VBR_Ts)) +
	theme_bw() +
	geom_abline(lty = 2) +
	geom_point(alpha = 0.2) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "VBR-Ts ET (mm/30min)",
			 x = "EC ET (mm/30min)",
			 title = paste0(site_id,": EC vs VBR-Ts")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)

d <- df_plot %>%
	ggplot(aes(ET_EC,ET_RE)) +
	theme_bw() +
	geom_abline(lty = 2) +
	geom_point(alpha = 0.2) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "RE ET (mm/30min)",
			 x = "EC ET (mm/30min)",
			 title = paste0(site_id,": EC vs RE")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)

e <- df_plot %>%
	filter(!is.na(ET_EC)) %>%
	ggplot(aes(ET_EC_EBC,ET_710)) +
	theme_bw() +
	geom_abline(lty = 2) +
	geom_point(alpha = 0.2) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "710 ET (mm/30min)",
			 x = "corrected EC ET (mm/30min)",
			 title = paste0(site_id,": EBC EC vs 710")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)

f <- df_plot %>%
	filter(!is.na(ET_EC) & !is.na(ET_710)) %>%
	ggplot(aes(ET_EC_EBC,ET_710_EBC)) +
	theme_bw() +
	geom_abline(lty = 2) +
	geom_point(alpha = 0.2) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "corrected 710 ET (mm/30min)",
			 x = "corrected EC ET (mm/30min)",
			 title = paste0(site_id,": EBC EC vs EBC 710")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)

g <- df_plot %>%
	filter(!is.na(ET_EC)) %>%
	ggplot(aes(ET_EC_EBC,ET_VBR_Ts)) +
	theme_bw() +
	geom_abline(lty = 2) +
	geom_point(alpha = 0.2) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "VBR-Ts ET (mm/30min)",
			 x = "corrected EC ET (mm/30min)",
			 title = paste0(site_id,": EBC EC vs VBR-Ts")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)

h <- df_plot %>%
	filter(!is.na(ET_EC)) %>%
	ggplot(aes(ET_EC_EBC,ET_RE)) +
	theme_bw() +
	geom_abline(lty = 2) +
	geom_point(alpha = 0.2) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "RE ET (mm/30min)",
			 x = "corrected EC ET (mm/30min)",
			 title = paste0(site_id,": EBC EC vs RE")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)

p1 <- plot_grid(a,b,c,d,e,f,g,h,ncol = 4)


a <- df_plot %>%
	mutate(date = date(TIMESTAMP - minutes(30))) %>%
	group_by(date) %>%
	summarise_if(is.numeric,sum) %>%
	ggplot(aes(ET_EC_F,ET_710_F)) +
	theme_bw() +
	geom_abline(lty = 2) +
	geom_point(alpha = 1) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "710 ET (mm/day)",
			 x = "EC ET (mm/day)",
			 title = paste0(site_id,": EC vs 710")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)

b <- df_plot %>%
	mutate(date = date(TIMESTAMP - minutes(30))) %>%
	group_by(date) %>%
	summarise_if(is.numeric,sum) %>%
	ggplot(aes(ET_EC_F,ET_710_EBC)) +
	theme_bw() +
	geom_abline(lty = 2) +
	geom_point(alpha = 1) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "corrected 710 ET (mm/day)",
			 x = "EC ET (mm/day)",
			 title = paste0(site_id,": EC vs EBC 710")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)

c <- df_plot %>%
	mutate(date = date(TIMESTAMP - minutes(30))) %>%
	group_by(date) %>%
	summarise_if(is.numeric,sum) %>%
	ggplot(aes(ET_EC_F,ET_VBR_Ts_F)) +
	theme_bw() +
	geom_abline(lty = 2) +
	geom_point(alpha = 1) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "VBR-Ts ET (mm/day)",
			 x = "EC ET (mm/day)",
			 title = paste0(site_id,": EC vs VBR-Ts")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)

d <- df_plot %>%
	mutate(date = date(TIMESTAMP - minutes(30))) %>%
	group_by(date) %>%
	summarise_if(is.numeric,sum) %>%
	ggplot(aes(ET_EC_F,ET_RE_F)) +
	theme_bw() +
	geom_abline(lty = 2) +
	geom_point(alpha = 1) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "RE ET (mm/day)",
			 x = "EC ET (mm/day)",
			 title = paste0(site_id,": EC vs RE")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)

e <- df_plot %>%
	mutate(date = date(TIMESTAMP - minutes(30))) %>%
	group_by(date) %>%
	summarise_if(is.numeric,sum) %>%
	ggplot(aes(ET_EC_EBC,ET_710_F)) +
	theme_bw() +
	geom_abline(lty = 2) +
	geom_point(alpha = 1) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "710 ET (mm/day)",
			 x = "corrected EC ET (mm/day)",
			 title = paste0(site_id,": EBC EC vs 710")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)

f <- df_plot %>%
	mutate(date = date(TIMESTAMP - minutes(30))) %>%
	group_by(date) %>%
	summarise_if(is.numeric,sum) %>%
	ggplot(aes(ET_EC_EBC,ET_710_EBC)) +
	theme_bw() +
	geom_abline(lty = 2) +
	geom_point(alpha = 1) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "corrected 710 ET (mm/day)",
			 x = "corrected EC ET (mm/day)",
			 title = paste0(site_id,": EBC EC vs EBC 710")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)

g <- df_plot %>%
	mutate(date = date(TIMESTAMP - minutes(30))) %>%
	group_by(date) %>%
	summarise_if(is.numeric,sum) %>%
	ggplot(aes(ET_EC_EBC,ET_VBR_Ts_F)) +
	theme_bw() +
	geom_abline(lty = 2) +
	geom_point(alpha = 1) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "VBR-Ts ET (mm/day)",
			 x = "corrected EC ET (mm/day)",
			 title = paste0(site_id,": EBC EC vs VBR-Ts")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)

h <- df_plot %>%
	mutate(date = date(TIMESTAMP - minutes(30))) %>%
	group_by(date) %>%
	summarise_if(is.numeric,sum) %>%
	ggplot(aes(ET_EC_EBC,ET_RE_F)) +
	theme_bw() +
	geom_abline(lty = 2) +
	geom_point(alpha = 1) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "RE ET (mm/day)",
			 x = "corrected EC ET (mm/day)",
			 title = paste0(site_id,": EBC EC vs RE")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)

p2 <- plot_grid(a,b,c,d,e,f,g,h,ncol = 4)

ggsave(plot = p1,
			 filename = pipeline_figure_path(paste0(site_id,"_Level_4_Figure_3_ET_intercomparison.jpeg")),
			 width = 11, height = 6.5, dpi = 150
)

ggsave(plot = p2,
       filename = pipeline_figure_path(paste0(site_id,"_Level_4_Figure_4_ET_intercomparison_daily.jpeg")),
       width = 11, height = 6.5, dpi = 150
)


 invisible(NULL)
}
