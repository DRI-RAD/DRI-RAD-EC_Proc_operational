# Render the original five L4 diagnostics from either pending or cumulative rows.
# No QAQC or gap-filling is rerun for cumulative figures.
pipeline_level4_figures <- function(df_plot, site_id, figure_dir) {
  if (!nrow(df_plot)) return(invisible(NULL))
  # CSV readers infer entirely missing numeric columns as logical. Keep these
  # channels in the original numeric daily summaries, including EC-only runs.
  missing_numeric <- vapply(df_plot, function(x) is.logical(x) && all(is.na(x)), logical(1))
  df_plot[missing_numeric] <- lapply(df_plot[missing_numeric], as.numeric)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  pipeline_figure_path <- function(name) file.path(figure_dir, name)
  EC_period <- range(df_plot$TIMESTAMP)
  LI710_period <- EC_period

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

b <- df_plot %>%
	mutate(DATE = date(TIMESTAMP - minutes(30))) %>%
	group_by(DATE) %>%
	summarise(across(where(is.numeric), ~ mean(.x))) %>%
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


d <- df_plot %>%
	filter(TIMESTAMP >= EC_period[1] & TIMESTAMP <= EC_period[2]) %>%
	mutate(DATE = date(TIMESTAMP - minutes(30))) %>%
	group_by(DATE) %>%
	summarise(across(where(is.numeric), ~ mean(.x, na.rm = FALSE))) %>%
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


b <- df_plot %>%
	mutate(DATE = date(TIMESTAMP - minutes(30))) %>%
	group_by(DATE) %>%
	summarise(across(where(is.numeric), ~ mean(.x))) %>%
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


d <- df_plot %>%
	filter(TIMESTAMP >= EC_period[1] & TIMESTAMP <= EC_period[2]) %>%
	mutate(DATE = date(TIMESTAMP - minutes(30))) %>%
	group_by(DATE) %>%
	summarise(across(where(is.numeric), ~ mean(.x, na.rm = FALSE))) %>%
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

b <- df_plot %>%
	mutate(DATE = date(TIMESTAMP - minutes(30))) %>%
	group_by(DATE) %>%
	summarise(across(where(is.numeric), ~ mean(.x))) %>%
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


d <- df_plot %>%
	filter(TIMESTAMP >= EC_period[1] & TIMESTAMP <= EC_period[2]) %>%
	mutate(DATE = date(TIMESTAMP - minutes(30))) %>%
	group_by(DATE) %>%
	summarise(across(where(is.numeric), ~ mean(.x, na.rm = FALSE))) %>%
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


b <- df_plot %>%
	mutate(DATE = date(TIMESTAMP - minutes(30))) %>%
	group_by(DATE) %>%
	summarise(across(where(is.numeric), ~ mean(.x))) %>%
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


d <- df_plot %>%
	filter(TIMESTAMP >= EC_period[1] & TIMESTAMP <= EC_period[2]) %>%
	mutate(DATE = date(TIMESTAMP - minutes(30))) %>%
	group_by(DATE) %>%
	summarise(across(where(is.numeric), ~ mean(.x, na.rm = FALSE))) %>%
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
offset <- 5
scale_factor <- 0.1


a <- df_plot %>%
	filter(TIMESTAMP >= EC_period[1] & TIMESTAMP <= EC_period[2]) %>%
	mutate(date = date(TIMESTAMP - minutes(30))) %>%
	group_by(date) %>%
	summarise_if(is.numeric,sum) %>%
	rowwise() %>%
	mutate(ET_max = max(ET_EC_F,ET_EC_EBC,na.rm=T),
				 ET_min = min(ET_EC_F,ET_EC_EBC,na.rm=T)
				 ) %>%
	ggplot(aes(x = date)) +
	geom_line(aes(y = offset - P_PI_F * scale_factor),
						color = "blue") +
	geom_ribbon(aes(ymax = ET_max, ymin = ET_min)) +
	scale_y_continuous(
		name = expression("EC ET (mm "*d^-1*")"),
		sec.axis = sec_axis(~ (offset - .)/scale_factor,
												name = "Precipitation (mm)")
	) +
	theme_bw() +
	labs(title = paste0(site_id,": Daily EC ET and Precipitation"))

b <- df_plot %>%
	filter(TIMESTAMP >= EC_period[1] & TIMESTAMP <= EC_period[2]) %>%
	mutate(date = date(TIMESTAMP - minutes(30))) %>%
	group_by(date) %>%
	summarise(P_PI_F = sum(P_PI_F),
						ET_EC_F = sum(ET_EC_F),
						ET_EC_EBC = sum(ET_EC_EBC),
						SWC_1_1_1 = mean(SWC_1_1_1,na.rm = TRUE),
						SWC_1_1_2 = mean(SWC_1_1_2,na.rm = TRUE),
						.groups = "drop"
	) %>%
	mutate(SWC = rowMeans(across(c(SWC_1_1_1, SWC_1_1_2)), na.rm = TRUE),
				 SSM = (SWC - first(SWC))*0.08 * 1000 / 100
				 ) %>%
	mutate(ET_cum = cumsum(ET_EC_F),
				 ET_cum_corr = cumsum(ET_EC_EBC),
				 Pr_cum = cumsum(P_PI_F)
				 ) %>%
	ggplot(aes(date)) +
	geom_ribbon(aes(ymax = ET_cum_corr, ymin = ET_cum, color = "EC ET"),alpha = 0.5) +
	geom_line(aes(y = Pr_cum,	color = "P")) +
	geom_line(aes(y = SSM,color = "SSM")) +
	theme_bw() +
  geom_text(
    data = ~ dplyr::slice_tail(.x, n = 1),
    aes(y = Pr_cum, label = round(Pr_cum,1)),
    hjust = -0.1, color = "#0072B2"
  ) +
  geom_text(
    data = ~ dplyr::slice_tail(.x, n = 1),
    aes(y = ET_cum_corr, label = round(ET_cum_corr,1)),
    hjust = -0.1, color = "#000000"
  ) +
  geom_text(
    data = ~ dplyr::slice_tail(.x, n = 1),
    aes(y = ET_cum, label = round(ET_cum,1)),
    hjust = -0.1, color = "#000000"
  ) +
  geom_text(
    data = ~ dplyr::slice_tail(.x, n = 1),
    aes(y = SSM, label = round(SSM,1)),
    hjust = -0.1, color = "#D55E00"
  ) +
	scale_color_manual(
		breaks = c("EC ET", "P","SSM"),
		values = c(
			"EC ET" = "#000000",        # black
			"P" = "#0072B2",        # blue
			"SSM" = "#D55E00" # vermillion
		)
	) +
	labs(title = paste0(site_id,": Cumulative EC ET vs precipitation"),
			 y = expression("Cumulative flux (mm)"))

c <- df_plot %>%
	filter(TIMESTAMP >= LI710_period[1] & TIMESTAMP <= LI710_period[2]) %>%
	filter(TIMESTAMP >= EC_period[1] & TIMESTAMP <= EC_period[2]) %>%
	mutate(date = date(TIMESTAMP - minutes(30))) %>%
	group_by(date) %>%
	summarise_if(is.numeric,sum) %>%
	rowwise() %>%
	mutate(ET_max = max(ET_710_F,ET_710_EBC,na.rm=T),
				 ET_min = min(ET_710_F,ET_710_EBC,na.rm=T)
	) %>%
	ggplot(aes(x = date)) +
	geom_line(aes(y = offset - P_PI_F * scale_factor),
						color = "blue") +
	geom_ribbon(aes(ymax = ET_max, ymin = ET_min)) +
	scale_y_continuous(
		name = expression("710 ET (mm "*d^-1*")"),
		sec.axis = sec_axis(~ (offset - .)/scale_factor,
												name = "Precipitation (mm)")
	) +
	theme_bw() +
	labs(title = paste0(site_id,": Daily 710 ET and Precipitation"))

d <- df_plot %>%
	filter(TIMESTAMP >= LI710_period[1] & TIMESTAMP <= LI710_period[2]) %>%
	filter(TIMESTAMP >= EC_period[1] & TIMESTAMP <= EC_period[2]) %>%
	mutate(date = date(TIMESTAMP - minutes(30))) %>%
	group_by(date) %>%
	summarise(P_PI_F = sum(P_PI_F),
						ET_EC_F = sum(ET_710_F),
						ET_EC_EBC = sum(ET_710_EBC),
						SWC_1_1_1 = mean(SWC_1_1_1,na.rm = TRUE),
						SWC_1_1_2 = mean(SWC_1_1_2,na.rm = TRUE),
						.groups = "drop"
	) %>%
	mutate(SWC = rowMeans(across(c(SWC_1_1_1, SWC_1_1_2)), na.rm = TRUE),
				 SSM = (SWC - first(SWC))*0.08 * 1000 / 100
	) %>%
	mutate(ET_cum = cumsum(ET_EC_F),
				 ET_cum_corr = cumsum(ET_EC_EBC),
				 Pr_cum = cumsum(P_PI_F)
	) %>%
	ggplot(aes(date)) +
	geom_ribbon(aes(ymax = ET_cum_corr, ymin = ET_cum, color = "710 ET"),alpha = 0.5) +
	geom_line(aes(y = Pr_cum,	color = "P")) +
	geom_line(aes(y = SSM,color = "SSM")) +
	theme_bw() +
  geom_text(
    data = ~ dplyr::slice_tail(.x, n = 1),
    aes(y = Pr_cum, label = round(Pr_cum,1)),
    hjust = -0.1, color = "#0072B2"
  ) +
  geom_text(
    data = ~ dplyr::slice_tail(.x, n = 1),
    aes(y = ET_cum_corr, label = round(ET_cum_corr,1)),
    hjust = -0.1, color = "#000000"
  ) +
  geom_text(
    data = ~ dplyr::slice_tail(.x, n = 1),
    aes(y = ET_cum, label = round(ET_cum,1)),
    hjust = -0.1, color = "#000000"
  ) +
  geom_text(
    data = ~ dplyr::slice_tail(.x, n = 1),
    aes(y = SSM, label = round(SSM,1)),
    hjust = -0.1, color = "#D55E00"
  ) +
	scale_color_manual(
		breaks = c("710 ET", "P","SSM"),
		values = c(
			"710 ET" = "#000000",        # black
			"P" = "#0072B2",        # blue
			"SSM" = "#D55E00" # vermillion
		)
	) +
	labs(title = paste0(site_id,": Cumulative 710 ET vs precipitation"),
			 y = expression("Cumulative flux (mm)"))

ggsave(plot = plot_grid(a,b,c,d,ncol = 1),
			 filename = pipeline_figure_path(paste0(site_id,"_Level_4_Figure_2_water_balance.jpeg")),
			 width = 11, height = 13, dpi = 150
)


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
	ggplot(aes(ET_EC,ET_VBR)) +
	theme_bw() +
	geom_abline(lty = 2) +
	geom_point(alpha = 0.2) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "VBR ET (mm/30min)",
			 x = "EC ET (mm/30min)",
			 title = paste0(site_id,": EC vs VBR")) +
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
	ggplot(aes(ET_EC_EBC,ET_VBR)) +
	theme_bw() +
	geom_abline(lty = 2) +
	geom_point(alpha = 0.2) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "VBR ET (mm/30min)",
			 x = "corrected EC ET (mm/30min)",
			 title = paste0(site_id,": EBC EC vs VBR")) +
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
	ggplot(aes(ET_EC_F,ET_VBR_F)) +
	theme_bw() +
	geom_abline(lty = 2) +
	geom_point(alpha = 1) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "VBR ET (mm/day)",
			 x = "EC ET (mm/day)",
			 title = paste0(site_id,": EC vs VBR")) +
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
			 title = paste0(site_id,"EBC EC vs EBC 710")) +
	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)

g <- df_plot %>%
	mutate(date = date(TIMESTAMP - minutes(30))) %>%
	group_by(date) %>%
	summarise_if(is.numeric,sum) %>%
	ggplot(aes(ET_EC_EBC,ET_VBR_F)) +
	theme_bw() +
	geom_abline(lty = 2) +
	geom_point(alpha = 1) +
	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
	geom_smooth(method = "lm", formula = "y ~ x + 0") +
	labs(y = "VBR ET (mm/day)",
			 x = "corrected EC ET (mm/day)",
			 title = paste0(site_id,": EBC EC vs VBR")) +
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
