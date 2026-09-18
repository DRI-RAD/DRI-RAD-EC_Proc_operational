# Run this stage through run_pipeline.R so period and output state are explicit.
if (is.null(getOption("ec.pipeline"))) stop("Use run_pipeline() from run_pipeline.R.")
pipeline_options <- getOption("ec.pipeline")

################################################################################
#### Load library and meta data
################################################################################
library(tidyverse)
library(lubridate)
library(zoo)
library(cowplot)
library(openeddy)
library(bigleaf)

# load direction meta
dirs <- pipeline_options$dirs

# load soil heat flux meta 
meta_G <- read_csv(file.path(.pipeline_root, "meta_files", "meta_G.csv"))
meta_G$Sun_plate_valid_from_1 <- as.POSIXct(meta_G$Sun_plate_valid_from_1, format = "%m/%d/%Y %H:%M", tz = "UTC")
meta_G$Sun_plate_valid_from_2 <- as.POSIXct(meta_G$Sun_plate_valid_from_2, format = "%m/%d/%Y %H:%M", tz = "UTC")
meta_G$Shade_plate_valid_from_1 <- as.POSIXct(meta_G$Shade_plate_valid_from_1, format = "%m/%d/%Y %H:%M", tz = "UTC")

# NAS direction for data
base_dir <- pipeline_options$base_dir

################################################################################
#### Load recent Level 1 data
################################################################################
#---------------------------------------
# site list 
#---------------------------------------
all_sites <- pipeline_options$sites

#---------------------------------------
# function: load one site
#---------------------------------------
load_latest_site <- function(site_id, dirs, base_dir) {
	
	dirs_use <- dirs %>% filter(site == site_id)
	
	files <- list.files(
		path = pipeline_path(base_dir, dirs_use$dir_output),
		pattern = "level_1_QAQCed_logger_data",
		full.names = TRUE,
		ignore.case = TRUE
	)
	
	
	latest_file <- pipeline_latest(files)
	
	df_unit <- read_csv(latest_file, skip = 1, n_max = 2)
	df <- read_csv(latest_file, skip = 4, col_names = FALSE)
	colnames(df) <- colnames(df_unit)
	# Restrict input to the requested calculation window before processing.
	df <- pipeline_slice(df)
	
	return(df)
}

#---------------------------------------
# load all sites
#---------------------------------------
site_list <- list()

for (i in seq_along(pipeline_options$reference_sites)) {
	site_list[[i]] <- load_latest_site(pipeline_options$reference_sites[i], dirs, base_dir)
}

#---------------------------------------
# naming
#---------------------------------------
names(site_list) <- pipeline_options$reference_sites

#---------------------------------------
# extract to envir
#---------------------------------------
# list2env(site_list, envir = .GlobalEnv)

################################################################################
#### long gap filling function 1: same site data
################################################################################
fill_long_gap_reg <- function(df, y_fill, y_ref, x_vars,
															maxgap = 4,
															window = 240,
															min_train = 50,
															force_zero = FALSE,
															min_pred_coverage = 0,
															prefer_full_gap = TRUE,
															min_predictors = 1) {
	
	out <- df[[y_fill]]
	y0  <- df[[y_ref]]
	
	is_na <- is.na(y0)
	if (!any(is_na)) return(out)
	
	r <- rle(is_na)
	ends <- cumsum(r$lengths)
	starts <- ends - r$lengths + 1
	
	gap_tbl <- data.frame(
		is_gap = r$values,
		start = starts,
		end = ends,
		len = r$lengths
	) %>%
		dplyr::filter(is_gap, len > maxgap)
	
	if (nrow(gap_tbl) == 0) return(out)
	
	for (j in seq_len(nrow(gap_tbl))) {
		
		g_start <- gap_tbl$start[j]
		g_end   <- gap_tbl$end[j]
		g_idx   <- g_start:g_end
		
		w_start <- max(1, g_start - window)
		w_end   <- min(nrow(df), g_end + window)
		
		train_idx <- setdiff(w_start:w_end, g_idx)
		
		x_avail <- x_vars[x_vars %in% names(df)]
		if (length(x_avail) == 0) next
		
		# training base
		train_df <- df[train_idx, c(y_ref, x_avail), drop = FALSE]
		names(train_df)[1] <- "y"
		
		#------------------------------------------
		# 1) training sufficiency + variation
		#------------------------------------------
		x_keep <- x_avail[sapply(x_avail, function(xn) {
			z <- train_df[[xn]]
			sum(!is.na(z)) >= min_train &&
				length(unique(z[!is.na(z)])) > 1
		})]
		
		if (length(x_keep) == 0) next
		
		#------------------------------------------
		# 2) prediction-period coverage
		#------------------------------------------
		pred_df_full <- df[g_idx, x_keep, drop = FALSE]
		
		x_keep <- x_keep[sapply(x_keep, function(xn) {
			z <- pred_df_full[[xn]]
			mean(!is.na(z)) >= min_pred_coverage
		})]
		
		if (length(x_keep) == 0) next
		
		#------------------------------------------
		# 3) full-gap predictor ordering
		#------------------------------------------
		if (prefer_full_gap) {
			full_gap_x <- x_keep[sapply(x_keep, function(xn) all(!is.na(pred_df_full[[xn]])))]
			partial_x  <- setdiff(x_keep, full_gap_x)
			x_keep <- c(full_gap_x, partial_x)
		}
		
		#------------------------------------------
		# 4) available predictor pattern
		#------------------------------------------
		pred_mat <- df[g_idx, x_keep, drop = FALSE]
		
		avail_list <- apply(pred_mat, 1, function(z) x_keep[!is.na(z)])
		
		pattern_key <- vapply(avail_list, function(v) {
			if (length(v) < min_predictors) return(NA_character_)
			paste(v, collapse = "|")
		}, character(1))
		
		unique_patterns <- unique(stats::na.omit(pattern_key))
		if (length(unique_patterns) == 0) next
		
		pred_val <- rep(NA_real_, length(g_idx))
		
		#------------------------------------------
		# 5) fit one model per available-predictor pattern
		#------------------------------------------
		for (pat in unique_patterns) {
			
			vars_use <- strsplit(pat, "\\|")[[1]]
			if (length(vars_use) < min_predictors) next
			
			train_use <- train_df[, c("y", vars_use), drop = FALSE]
			train_use <- train_use[complete.cases(train_use), , drop = FALSE]
			
			if (nrow(train_use) < min_train) next
			
			form <- if (force_zero) {
				as.formula(paste("y ~ 0 +", paste(vars_use, collapse = " + ")))
			} else {
				as.formula(paste("y ~", paste(vars_use, collapse = " + ")))
			}
			
			fit <- try(lm(form, data = train_use), silent = TRUE)
			if (inherits(fit, "try-error")) next
			
			idx_pat <- which(pattern_key == pat)
			if (length(idx_pat) == 0) next
			
			newdata_pat <- pred_mat[idx_pat, vars_use, drop = FALSE]
			ok_pat <- complete.cases(newdata_pat)
			if (!any(ok_pat)) next
			
			pred_tmp <- rep(NA_real_, nrow(newdata_pat))
			pred_tmp[ok_pat] <- predict(fit, newdata = newdata_pat[ok_pat, , drop = FALSE])
			
			pred_val[idx_pat] <- pred_tmp
		}
		
		fillable <- is.na(out[g_idx]) & !is.na(pred_val)
		out[g_idx][fillable] <- pred_val[fillable]
	}
	
	return(out)
}


################################################################################
#### long gap filling function 2: different site data
################################################################################
fill_long_gap_reg_other_sites <- function(df_target,
																					site_name,
																					site_list,
																					y_fill,
																					y_ref,
																					x_var,
																					maxgap = 4,
																					window = 240,
																					min_train = 50,
																					force_zero = FALSE,
																					min_pred_coverage = 0,
																					prefer_full_gap = TRUE,
																					min_predictors = 1) {
	
	out <- df_target[[y_fill]]
	y0  <- df_target[[y_ref]]
	
	# identify long gaps from original reference variable
	is_na <- is.na(y0)
	if (!any(is_na)) return(out)
	
	r <- rle(is_na)
	ends <- cumsum(r$lengths)
	starts <- ends - r$lengths + 1
	
	gap_tbl <- data.frame(
		is_gap = r$values,
		start = starts,
		end = ends,
		len = r$lengths
	) %>%
		dplyr::filter(is_gap, len > maxgap)
	
	if (nrow(gap_tbl) == 0) return(out)
	
	# candidate sites = all except target site
	other_sites <- setdiff(names(site_list), site_name)
	
	for (j in seq_len(nrow(gap_tbl))) {
		
		g_start <- gap_tbl$start[j]
		g_end   <- gap_tbl$end[j]
		g_idx   <- g_start:g_end
		
		w_start <- max(1, g_start - window)
		w_end   <- min(nrow(df_target), g_end + window)
		
		# target training data outside the gap but within the local window
		target_train <- df_target[w_start:w_end, c("TIMESTAMP", y_ref), drop = FALSE]
		target_train$in_gap <- seq_len(nrow(target_train)) %in% (g_idx - w_start + 1)
		target_train <- target_train %>%
			dplyr::filter(!in_gap) %>%
			dplyr::select(-in_gap)
		
		# build merged training/prediction table
		merged_train <- target_train
		merged_pred  <- df_target[g_idx, c("TIMESTAMP"), drop = FALSE]
		
		predictor_names <- c()
		
		for (s in other_sites) {
			df_other <- site_list[[s]]
			
			if (!all(c("TIMESTAMP", x_var) %in% names(df_other))) next
			
			new_name <- paste0(x_var, "_", s)
			
			tmp <- df_other %>%
				dplyr::select(TIMESTAMP, !!rlang::sym(x_var)) %>%
				dplyr::rename(!!new_name := !!rlang::sym(x_var))
			
			merged_train <- merged_train %>% dplyr::left_join(tmp, by = "TIMESTAMP")
			merged_pred  <- merged_pred  %>% dplyr::left_join(tmp, by = "TIMESTAMP")
			
			predictor_names <- c(predictor_names, new_name)
		}
		
		if (length(predictor_names) == 0) next
		
		train_df <- merged_train %>%
			dplyr::rename(y = !!rlang::sym(y_ref))
		
		x_keep <- predictor_names[predictor_names %in% names(train_df)]
		if (length(x_keep) == 0) next
		
		# 1) training sufficiency + variation
		x_keep <- x_keep[sapply(x_keep, function(xn) {
			z <- train_df[[xn]]
			sum(!is.na(z)) >= min_train && length(unique(z[!is.na(z)])) > 1
		})]
		if (length(x_keep) == 0) next
		
		# 2) prediction-period coverage filter
		x_keep <- x_keep[sapply(x_keep, function(xn) {
			z <- merged_pred[[xn]]
			mean(!is.na(z)) >= min_pred_coverage
		})]
		if (length(x_keep) == 0) next
		
		# 3) full-gap predictors first (optional)
		if (prefer_full_gap) {
			full_gap_x <- x_keep[sapply(x_keep, function(xn) all(!is.na(merged_pred[[xn]])))]
			partial_x  <- setdiff(x_keep, full_gap_x)
			x_keep <- c(full_gap_x, partial_x)
		}
		
		# available predictors pattern for each prediction row
		pred_mat <- merged_pred[, x_keep, drop = FALSE]
		avail_list <- apply(pred_mat, 1, function(z) x_keep[!is.na(z)])
		
		# encode predictor set as string
		pattern_key <- vapply(avail_list, function(v) {
			if (length(v) < min_predictors) return(NA_character_)
			paste(v, collapse = "|")
		}, character(1))
		
		unique_patterns <- unique(stats::na.omit(pattern_key))
		if (length(unique_patterns) == 0) next
		
		pred_val <- rep(NA_real_, nrow(merged_pred))
		
		# fit one model per available-predictor pattern
		for (pat in unique_patterns) {
			vars_use <- strsplit(pat, "\\|")[[1]]
			
			if (length(vars_use) < min_predictors) next
			
			# training rows complete for y and vars_use
			train_use <- train_df[, c("y", vars_use), drop = FALSE]
			train_use <- train_use[complete.cases(train_use), , drop = FALSE]
			
			if (nrow(train_use) < min_train) next
			
			# formula
			form <- if (force_zero) {
				as.formula(paste("y ~ 0 +", paste(vars_use, collapse = " + ")))
			} else {
				as.formula(paste("y ~", paste(vars_use, collapse = " + ")))
			}
			
			fit <- try(lm(form, data = train_use), silent = TRUE)
			if (inherits(fit, "try-error")) next
			
			idx_pat <- which(pattern_key == pat)
			if (length(idx_pat) == 0) next
			
			newdata_pat <- merged_pred[idx_pat, vars_use, drop = FALSE]
			ok_pat <- complete.cases(newdata_pat)
			if (!any(ok_pat)) next
			
			pred_tmp <- rep(NA_real_, nrow(newdata_pat))
			pred_tmp[ok_pat] <- predict(fit, newdata = newdata_pat[ok_pat, , drop = FALSE])
			
			pred_val[idx_pat] <- pred_tmp
		}
		
		fillable <- is.na(out[g_idx]) & !is.na(pred_val)
		out[g_idx][fillable] <- pred_val[fillable]
	}
	
	return(out)
}

################################################################################
#### long gap filling function 3: mean diurnal cycle
################################################################################
fill_mdc_local <- function(df, y_fill, timestamp = "TIMESTAMP",
													 window_days = 5,
													 nmin = 3) {
	
	out <- df[[y_fill]]
	time_vec <- df[[timestamp]]
	
	# time-of-day slot for 30-min data
	hhmm <- format(time_vec, "%H:%M")
	
	na_idx <- which(is.na(out))
	if (length(na_idx) == 0) return(out)
	
	for (i in na_idx) {
		t0 <- time_vec[i]
		slot0 <- hhmm[i]
		
		# local ±window_days
		use_idx <- which(
			!is.na(out) &
				hhmm == slot0 &
				time_vec >= (t0 - days(window_days)) &
				time_vec <= (t0 + days(window_days))
		)
		
		if (length(use_idx) >= nmin) {
			out[i] <- mean(out[use_idx], na.rm = TRUE)
		}
	}
	
	return(out)
}


################################################################################
#### function: visualize gap-filled PI variables
################################################################################
plot_gapfilled_pdf <- function(df,
															 site_id,
															 file_out = "gapfilled_PI.pdf",
															 chunk_size = 48 * 365 / 2,
															 n_panel_per_page = 4,
															 max_n_per_panel = 3000,
															 draw_points = TRUE) {
	
	stopifnot("TIMESTAMP" %in% names(df))
	
	vars_use <- c("TA_PI_F", "RH_PI_F", "SW_IN_PI_F", "LW_IN_PI_F", "P_PI_F")
	
	method_key <- c(
		"0" = "Original",
		"1" = "Short interpolation",
		"2" = "Same-site regression",
		"3" = "Other-site regression",
		"4" = "Final fill / MDC"
	)
	
	df <- df %>%
		filter(!is.na(TIMESTAMP)) %>%
		arrange(TIMESTAMP)
	
	n <- nrow(df)
	starts <- seq(1, n, by = chunk_size)
	
	#--------------------------------------------
	# collect plots
	#--------------------------------------------
	all_plots <- list()
	k <- 1
	
	for (s in starts) {
		
		e <- min(s + chunk_size - 1, n)
		df_sub <- df[s:e, ]
		
		if (nrow(df_sub) == 0) next
		
		# downsample
		if (nrow(df_sub) > max_n_per_panel) {
			step_n <- ceiling(nrow(df_sub) / max_n_per_panel)
			df_sub <- df_sub[seq(1, nrow(df_sub), by = step_n), ]
		}
		
		for (v in vars_use) {
			
			method_v <- paste0(v, "_method")
			
			dat_plot <- df_sub %>%
				transmute(
					TIMESTAMP,
					value = .data[[v]],
					method = factor(
						as.character(.data[[method_v]]),
						levels = names(method_key),
						labels = unname(method_key)
					)
				)
			
			p <- ggplot(dat_plot, aes(x = TIMESTAMP)) +
				geom_line(aes(y = value), color = "black", linewidth = 0.3, na.rm = TRUE)
			
			if (draw_points) {
				p <- p +
					geom_point(
						data = dat_plot %>% filter(method != "Original"),
						aes(y = value, color = method),
						size = 0.4,
						alpha = 0.9,
						na.rm = TRUE
					)
			}
			
			p <- p +
				scale_color_manual(
					values = c(
						"Original" = "black",
						"Short interpolation" = "dodgerblue3",
						"Same-site regression" = "orange2",
						"Other-site regression" = "purple3",
						"Final fill / MDC" = "red3"
					),
					drop = FALSE
				) +
				labs(
					title = paste0(site_id, " | ", v, "  (rows ", s, "-", e, ")"),
					subtitle = paste0(
						format(min(dat_plot$TIMESTAMP, na.rm = TRUE), "%Y-%m-%d %H:%M"),
						" to ",
						format(max(dat_plot$TIMESTAMP, na.rm = TRUE), "%Y-%m-%d %H:%M")
					),
					x = NULL,
					y = v
				) +
				theme_bw() +
				theme(
					plot.title = element_text(face = "bold", size = 10),
					plot.subtitle = element_text(size = 8),
					axis.text.x = element_text(size = 7),
					axis.text.y = element_text(size = 7),
					axis.title.y = element_text(size = 8),
					panel.grid.minor = element_blank(),
					#legend.position = "none",
					plot.margin = margin(3, 6, 3, 3)
				)
			
			all_plots[[k]] <- p
			k <- k + 1
		}
	}
	
	#--------------------------------------------
	# split pages (QC style)
	#--------------------------------------------
	page_id <- ceiling(seq_along(all_plots) / n_panel_per_page)
	pages <- split(all_plots, page_id)
	
	pdf(file_out, width = 12, height = 15, onefile = TRUE, useDingbats = FALSE)
	
	for (i in seq_along(pages)) {
		page_plot <- plot_grid(
			plotlist = pages[[i]],
			ncol = 1,
			align = "v",
			axis = "lr"
		)
		
		print(page_plot)
	}
	
	dev.off()
	
	message("Saved: ", file_out)
}

################################################################################
#### PI variables & gap filling
################################################################################
# To do (if there is needs): PPFD_PI_F

for(i in all_sites) {
	df <- site_list[[i]]
	
	# Select HMP temperature and humidity as represented data. 
	df_PI <- df %>%
		mutate(SW_IN_2 = PPFD_IN*0.48,
					 TA_PI_F = TA_1_1_3,
					 RH_PI_F = RH_1_1_3,
					 PA_PI_F = PA,
					 SW_IN_PI_F = SW_IN,
					 LW_IN_PI_F = LW_IN,
					 
					 TA_PI_F_method = ifelse(!is.na(TA_1_1_3), 0, NA),
					 RH_PI_F_method =  ifelse(!is.na(RH_1_1_3), 0, NA),
					 PA_PI_F_method =  ifelse(!is.na(PA), 0, NA),
					 SW_IN_PI_F_method =  ifelse(!is.na(SW_IN), 0, NA),
					 LW_IN_PI_F_method =  ifelse(!is.na(LW_IN), 0, NA),
					 
					 ) %>%
	#--------------------------------------------------
	# Gap filling step 1: short gaps (≤ 4 hours for met, 2 hours for SW)
	#--------------------------------------------------
		mutate(TA_PI_F = zoo::na.approx(TA_PI_F, maxgap = 8, na.rm = FALSE),
					 RH_PI_F = zoo::na.approx(RH_PI_F, maxgap = 8, na.rm = FALSE),
					 PA_PI_F = zoo::na.approx(PA_PI_F, maxgap = 8, na.rm = FALSE),
					 SW_IN_PI_F = zoo::na.approx(SW_IN_PI_F, maxgap = 4, na.rm = FALSE),
					 LW_IN_PI_F = zoo::na.approx(LW_IN_PI_F, maxgap = 8, na.rm = FALSE),
					 
					 TA_PI_F_method = ifelse(is.na(TA_PI_F_method) & !is.na(TA_PI_F), 1, TA_PI_F_method),
					 RH_PI_F_method = ifelse(is.na(RH_PI_F_method) & !is.na(RH_PI_F), 1, RH_PI_F_method),
					 PA_PI_F_method = ifelse(is.na(PA_PI_F_method) & !is.na(PA_PI_F), 1, PA_PI_F_method),
					 SW_IN_PI_F_method = ifelse(is.na(SW_IN_PI_F_method) & !is.na(SW_IN_PI_F), 1, SW_IN_PI_F_method),
					 LW_IN_PI_F_method = ifelse(is.na(LW_IN_PI_F_method) & !is.na(LW_IN_PI_F), 1, LW_IN_PI_F_method)
					 
					 ) 
	
	# precipitation (ERVA doesn't have bulk precipitation guage)
	# EDVG should use corrected P_Bulk while other sites don't require correction.
	if(i == "ERVA"){
		df_PI <- df_PI %>%
			mutate(P_PI_F = P,
						 P_PI_F_method = ifelse(!is.na(P_PI_F), 0, NA)
						 )
	} else{
		df_PI <- df_PI %>%
			mutate(P_PI_F = P_Bulk,
						 P_PI_F_method = ifelse(!is.na(P_PI_F), 0, NA)
						 )
	}
	
	
	#--------------------------------------------------
	# Gap filling step 2: long gaps using same-site sensors
	#--------------------------------------------------
	ta_xvars <- intersect(c("TA_1_1_1", "TA_1_1_2","TA_710"), names(df_PI))
	if (length(ta_xvars) > 0) {
		df_PI$TA_PI_F <- fill_long_gap_reg(
			df = df_PI,
			y_fill = "TA_PI_F",
			y_ref = "TA_1_1_3",
			x_vars = ta_xvars,
			maxgap = 8,
			window = 240,
			min_train = 50,
			force_zero = FALSE
		)
	}
	
	
	RH_xvars <- intersect(c("RH_1_1_1", "RH_1_1_2","RH_710"), names(df_PI))
	if (length(RH_xvars) > 0) {
		df_PI$RH_PI_F <- fill_long_gap_reg(
			df = df_PI,
			y_fill = "RH_PI_F",
			y_ref = "RH_1_1_3",
			x_vars = RH_xvars,
			maxgap = 8,
			window = 240,
			min_train = 50,
			force_zero = FALSE
		)
	}
	
	PA_xvars <- intersect(c("PA_710"), names(df_PI))
	if (length(PA_xvars) > 0) {
		df_PI$PA_PI_F <- fill_long_gap_reg(
			df = df_PI,
			y_fill = "PA_PI_F",
			y_ref = "PA",
			x_vars = PA_xvars,
			maxgap = 8,
			window = 240,
			min_train = 50,
			force_zero = FALSE
		)
	}
	
	SW_IN_xvars <- intersect(c("SW_IN_2"), names(df_PI))
	if (length(SW_IN_xvars) > 0) {
		df_PI$SW_IN_PI_F <- fill_long_gap_reg(
			df = df_PI,
			y_fill = "SW_IN_PI_F",
			y_ref = "SW_IN",
			x_vars = SW_IN_xvars,
			maxgap = 4,
			window = 240,
			min_train = 50,
			force_zero = TRUE
		)
	}
	
	SW_IN_xvars <- intersect(c("SW_IN_2"), names(df_PI))
	if (length(SW_IN_xvars) > 0) {
		df_PI$SW_IN_PI_F <- fill_long_gap_reg(
			df = df_PI,
			y_fill = "SW_IN_PI_F",
			y_ref = "SW_IN",
			x_vars = SW_IN_xvars,
			maxgap = 4,
			window = 240,
			min_train = 50,
			force_zero = TRUE
		)
	}
	
	
	P_xvars <- intersect(c("P"), names(df_PI))
	if (length(SW_IN_xvars) > 0) {
		df_PI$P_PI_F <- fill_long_gap_reg(
			df = df_PI,
			y_fill = "P_PI_F",
			y_ref = "P_PI_F",
			x_vars = P_xvars,
			maxgap = 0,
			window = 4800,
			min_train = 50,
			force_zero = TRUE
		)
	}
	
	
	df_PI <- df_PI %>%
	mutate(TA_PI_F_method = ifelse(is.na(TA_PI_F_method) & !is.na(TA_PI_F), 2, TA_PI_F_method),
				 RH_PI_F_method = ifelse(is.na(RH_PI_F_method) & !is.na(RH_PI_F), 2, RH_PI_F_method),
				 PA_PI_F_method = ifelse(is.na(PA_PI_F_method) & !is.na(PA_PI_F), 2, PA_PI_F_method),
				 SW_IN_PI_F_method = ifelse(is.na(SW_IN_PI_F_method) & !is.na(SW_IN_PI_F), 2, SW_IN_PI_F_method),
				 P_PI_F_method = ifelse(is.na(P_PI_F_method) & !is.na(P_PI_F), 2, P_PI_F_method)
	) 
	
	
	#--------------------------------------------------
	# Gap filling step 3: long gaps using different site data
	#--------------------------------------------------
	df_PI$TA_PI_F <- fill_long_gap_reg_other_sites(
		df_target = df_PI,
		site_name = i,
		site_list = site_list,
		y_fill = "TA_PI_F",
		y_ref = "TA_1_1_3",
		x_var = "TA_1_1_3",
		maxgap = 8,
		window = 480,
		min_train = 50,
		force_zero = FALSE
	)
	
	df_PI$RH_PI_F <- fill_long_gap_reg_other_sites(
		df_target = df_PI,
		site_name = i,
		site_list = site_list,
		y_fill = "RH_PI_F",
		y_ref = "RH_1_1_3",
		x_var = "RH_1_1_3",
		maxgap = 8,
		window = 480,
		min_train = 50,
		force_zero = FALSE
	)
	
	df_PI$PA_PI_F <- fill_long_gap_reg_other_sites(
		df_target = df_PI,
		site_name = i,
		site_list = site_list,
		y_fill = "PA_PI_F",
		y_ref = "PA",
		x_var = "PA",
		maxgap = 8,
		window = 480,
		min_train = 50,
		force_zero = FALSE
	)
	
	df_PI$SW_IN_PI_F <- fill_long_gap_reg_other_sites(
		df_target = df_PI,
		site_name = i,
		site_list = site_list,
		y_fill = "SW_IN_PI_F",
		y_ref = "SW_IN",
		x_var = "SW_IN",
		maxgap = 4,
		window = 480,
		min_train = 50,
		force_zero = TRUE
	)
	
	df_PI$LW_IN_PI_F <- fill_long_gap_reg_other_sites(
		df_target = df_PI,
		site_name = i,
		site_list = site_list,
		y_fill = "LW_IN_PI_F",
		y_ref = "LW_IN",
		x_var = "LW_IN",
		maxgap = 8,
		window = 480,
		min_train = 50,
		force_zero = FALSE
	)
	
	df_PI$P_PI_F <- fill_long_gap_reg_other_sites(
		df_target = df_PI,
		site_name = i,
		site_list = site_list,
		y_fill = "P_PI_F",
		y_ref = "P_PI_F",
		x_var = "P_Bulk",
		maxgap = 0,
		window = 48*60, #+- two month
		min_train = 50,
		force_zero = TRUE
	)
	
	df_PI <- df_PI %>%
		mutate(TA_PI_F_method = ifelse(is.na(TA_PI_F_method) & !is.na(TA_PI_F), 3, TA_PI_F_method),
					 RH_PI_F_method = ifelse(is.na(RH_PI_F_method) & !is.na(RH_PI_F), 3, RH_PI_F_method),
					 PA_PI_F_method = ifelse(is.na(PA_PI_F_method) & !is.na(PA_PI_F), 3, PA_PI_F_method),
					 SW_IN_PI_F_method = ifelse(is.na(SW_IN_PI_F_method) & !is.na(SW_IN_PI_F), 3, SW_IN_PI_F_method),
					 LW_IN_PI_F_method = ifelse(is.na(LW_IN_PI_F_method) & !is.na(LW_IN_PI_F), 3, LW_IN_PI_F_method),
					 P_PI_F_method = ifelse(is.na(P_PI_F_method) & !is.na(P_PI_F), 3, P_PI_F_method),
		) 
	
	
	#--------------------------------------------------
	# Gap filling step for LW: modelled LW from same site
	#--------------------------------------------------
	
	# todo: ASCE considering cloud fraction....
	df_PI <- df_PI %>%
		mutate(esat = Esat.slope(TA_PI_F, formula = "Allen_1998")$Esat,
					 e = RH_PI_F/100 *esat,
					 LW_IN_Abramowitz = 2.84*(TA_PI_F + 273.15) + 0.031*e*1000 - 522.5,
					 LW_IN_ASCE = 5.670367e-08 * (1 - 0.34 + 0.14*sqrt(e))*(TA_PI_F + 273.15)^4
		)
	LW_IN_xvars <- intersect(c("LW_IN_Abramowitz","LW_IN_ASCE"), names(df_PI))
	if (length(LW_IN_xvars) > 0) {
		df_PI$LW_IN_PI_F <- fill_long_gap_reg(
			df = df_PI,
			y_fill = "LW_IN_PI_F",
			y_ref = "LW_IN",
			x_vars = LW_IN_xvars,
			maxgap = 8,
			window = 240,
			min_train = 50,
			force_zero = FALSE
		)
	}
	
	df_PI <- df_PI %>%
		mutate(LW_IN_PI_F_method = ifelse(is.na(LW_IN_PI_F_method) & !is.na(LW_IN_PI_F), 2, LW_IN_PI_F_method)) 
	
	#--------------------------------------------------
	# Gap filling step 3: mean diurnal cycle
	#--------------------------------------------------
	df_PI$SW_IN_PI_F <- fill_mdc_local(
		df = df_PI,
		y_fill = "SW_IN_PI_F",
		timestamp = "TIMESTAMP",
		window_days = 5,
		nmin = 3
	)
	
	df_PI$LW_IN_PI_F <- fill_mdc_local(
		df = df_PI,
		y_fill = "LW_IN_PI_F",
		timestamp = "TIMESTAMP",
		window_days = 5,
		nmin = 3
	)
	
	df_PI$PA_PI_F <- fill_mdc_local(
		df = df_PI,
		y_fill = "PA_PI_F",
		timestamp = "TIMESTAMP",
		window_days = 5,
		nmin = 3
	)
	
	df_PI$TA_PI_F <- fill_mdc_local(
		df = df_PI,
		y_fill = "TA_PI_F",
		timestamp = "TIMESTAMP",
		window_days = 5,
		nmin = 3
	)
	
	df_PI$RH_PI_F <- fill_mdc_local(
		df = df_PI,
		y_fill = "RH_PI_F",
		timestamp = "TIMESTAMP",
		window_days = 5,
		nmin = 3
	)
	
	df_PI <- df_PI %>%
		mutate(TA_PI_F_method = ifelse(is.na(TA_PI_F_method) & !is.na(TA_PI_F), 4, TA_PI_F_method),
					 RH_PI_F_method = ifelse(is.na(RH_PI_F_method) & !is.na(RH_PI_F), 4, RH_PI_F_method),
					 PA_PI_F_method = ifelse(is.na(PA_PI_F_method) & !is.na(PA_PI_F), 4, PA_PI_F_method),
					 SW_IN_PI_F_method = ifelse(is.na(SW_IN_PI_F_method) & !is.na(SW_IN_PI_F), 4, SW_IN_PI_F_method),
					 LW_IN_PI_F_method = ifelse(is.na(LW_IN_PI_F_method) & !is.na(LW_IN_PI_F), 4, LW_IN_PI_F_method),
					 
					 P_PI_F = ifelse(is.na(P_PI_F), 0,P_PI_F),
					 P_PI_F_method = ifelse(is.na(P_PI_F_method) & !is.na(P_PI_F), 4, P_PI_F_method),
		)
	
	
	print(paste0(i," TA Gap: ",sum(is.na(df_PI$TA_PI_F))))
	print(paste0(i," RH Gap: ",sum(is.na(df_PI$RH_PI_F))))
	print(paste0(i," PA Gap: ",sum(is.na(df_PI$PA_PI_F))))
	print(paste0(i," SW Gap: ",sum(is.na(df_PI$SW_IN_PI_F))))
	print(paste0(i," LW Gap: ",sum(is.na(df_PI$LW_IN_PI_F))))
	print(paste0(i," P Gap: ",sum(is.na(df_PI$P_PI_F))))
	
	dirs_use <- dirs %>% filter(site == i)
	
	#--------------------------------------------------
	# plot
	#--------------------------------------------------
	plot_gapfilled_pdf(
		df = pipeline_figure_data(df_PI),
		site_id = i,
		file_out = paste0(
			pipeline_path(base_dir, dirs_use$dir_output),
			"/figures/",
			i,
			"_Level_2_Figure_gapfilled_PI.pdf"
		),
		chunk_size = 48 * 365 / 4,
		max_n_per_panel = 4000
	)
	
	#--------------------------------------------------
	# vpd
	#--------------------------------------------------
	df_PI$VPD_PI_F <- rH.to.VPD(rH = df_PI$RH_PI_F/100, 
														Tair = df_PI$TA_PI_F, 
														Esat.formula = "Allen_1998") * 10 #kPa -> hPa (Ameriflux unit)
	
	# save
	site_list[[i]] <- df_PI
}


################################################################################
#### Soil surface soil heat flux (G0)
################################################################################
# SG_PI
# G_plate_PI
# G0_PI 

for(i in all_sites) {
	df <- site_list[[i]]
	
	dirs_use <- dirs %>% filter(site == i)
	meta1 <- meta_G %>% filter(site == i)
	
	#--------------------------------------------------
	# Sun plate 1 -> G_plate_sun_1
	#--------------------------------------------------
	if (!all(is.na(meta1[, c("Sun_plate_1_1_1", "Sun_plate_1_1_2")]))) {
		
		v1 <- meta1$Sun_plate_1_1_1
		v2 <- meta1$Sun_plate_1_1_2
		t0 <- meta1$Sun_plate_valid_from_1
		
		vars_use <- c(v1, v2)
		vars_use <- vars_use[!is.na(vars_use)]
		vars_use <- vars_use[vars_use %in% names(df)]
		
		if (length(vars_use) > 0) {
			df$G_plate_sun_1 <- rowMeans(df[, vars_use, drop = FALSE], na.rm = TRUE)
			df$G_plate_sun_1[rowSums(!is.na(df[, vars_use, drop = FALSE])) == 0] <- NA
			
			if (!is.na(t0)) {
				df$G_plate_sun_1[df$TIMESTAMP < t0] <- NA
			}
		}
		
	}
	
	#--------------------------------------------------
	# Sun plate 2 -> G_plate_sun_2
	#--------------------------------------------------
	if (!all(is.na(meta1[, c("Sun_plate_2_1_1", "Sun_plate_2_1_2")]))
			& meta1$Sun_plate_valid_from_2 < last(df$TIMESTAMP)) {
		
		v1 <- meta1$Sun_plate_2_1_1
		v2 <- meta1$Sun_plate_2_1_2
		t0 <- meta1$Sun_plate_valid_from_2
		
		vars_use <- c(v1, v2)
		vars_use <- vars_use[!is.na(vars_use)]
		vars_use <- vars_use[vars_use %in% names(df)]
		
		if (length(vars_use) > 0) {
			df$G_plate_sun_2 <- rowMeans(df[, vars_use, drop = FALSE], na.rm = TRUE)
			df$G_plate_sun_2[rowSums(!is.na(df[, vars_use, drop = FALSE])) == 0] <- NA
			
			if (!is.na(t0)) {
				df$G_plate_sun_2[df$TIMESTAMP < t0] <- NA
			}
		}
	
	}
	
	#--------------------------------------------------
	# Shade plate 1 -> G_plate_shade
	#--------------------------------------------------
	if (!all(is.na(meta1[, c("Shade_plate_1_1_1", "Shade_plate_1_1_2")]))) {
		
		v1 <- meta1$Shade_plate_1_1_1
		v2 <- meta1$Shade_plate_1_1_2
		t0 <- meta1$Shade_plate_valid_from_1
		
		vars_use <- c(v1, v2)
		vars_use <- vars_use[!is.na(vars_use)]
		vars_use <- vars_use[vars_use %in% names(df)]
		
		if (length(vars_use) > 0) {
			df$G_plate_shade <- rowMeans(df[, vars_use, drop = FALSE], na.rm = TRUE)
			df$G_plate_shade[rowSums(!is.na(df[, vars_use, drop = FALSE])) == 0] <- NA
			
			if (!is.na(t0)) {
				df$G_plate_shade[df$TIMESTAMP < t0] <- NA
			}
		}
		
	}
	
	#--------------------------------------------------
	# Sun SG 1
	#--------------------------------------------------
	if (!all(is.na(meta1[, c("Sun_TS_1_1_1", "Sun_TS_1_1_2","Sun_SWC_1")]))){
		v1 <- meta1$Sun_TS_1_1_1
		v2 <- meta1$Sun_TS_1_1_2
		v3 <- meta1$Sun_SWC_1
		t0 <- meta1$Sun_plate_valid_from_1
		
		TS_1 <- df %>% pull(v1)
		TS_1[df$TIMESTAMP < t0] <- NA
		
		TS_2 <- df %>% pull(v2)
		TS_2[df$TIMESTAMP < t0] <- NA
		
		SWC <- df %>% pull(v3)
		SWC[is.na(SWC)] <- median(SWC,na.rm=T)
		SWC[df$TIMESTAMP < t0] <- NA
		
		cp <- (2.31 * meta1$Organic_fraction  +	
					 	2.5 * (meta1$Bulk_density  / 2.6 - meta1$Organic_fraction) +	
					 	4.18 * SWC/100) * 1e6
		
		dTdt1 <- c(diff(TS_1)[1], diff(TS_1)) / (30 * 60)
		dTdt2 <- c(diff(TS_2)[2], diff(TS_2)) / (30 * 60)
		
		x <- data.frame(timestamp = df$TIMESTAMP, dTdt1 = dTdt1, dTdt2 = dTdt2)
		flag <- despikeLF(x,"dTdt1",iter = 10,	light = NULL,	z = 10,	c = 7, var_thr = c(-0.0022,0.0022))
		flag[is.na(flag)] <- 0
		x$dTdt1 <- ifelse(flag > 0, NA, x$dTdt1)
		
		flag <- despikeLF(x,"dTdt2",iter = 10,	light = NULL,	z = 10,	c = 7, var_thr = c(-0.0022,0.0022))
		flag[is.na(flag)] <- 0
		x$dTdt2 <- ifelse(flag > 0, NA, x$dTdt2)
		
		dTdt <- rowMeans(x[, c("dTdt1","dTdt2"), drop = FALSE], na.rm = TRUE)
		dTdt = zoo::na.approx(dTdt, maxgap = 8, na.rm = FALSE)
		
		df$SG_sun_1 <- cp * meta1$Sun_plate_depth_1 * dTdt
	}
	
	#--------------------------------------------------
	# Sun SG 2
	#--------------------------------------------------
	if (!all(is.na(meta1[, c("Sun_TS_2_1_1", "Sun_TS_2_1_2","Sun_SWC_2")]))
			& meta1$Sun_plate_valid_from_2 < last(df$TIMESTAMP)
			){
		v1 <- meta1$Sun_TS_2_1_1
		v2 <- meta1$Sun_TS_2_1_2
		v3 <- meta1$Sun_SWC_2
		t0 <- meta1$Sun_plate_valid_from_2
		
		TS_1 <- df %>% pull(v1)
		TS_1[df$TIMESTAMP < t0] <- NA
		
		TS_2 <- df %>% pull(v2)
		TS_2[df$TIMESTAMP < t0] <- NA
		
		SWC <- df %>% pull(v3)
		SWC[is.na(SWC)] <- median(SWC,na.rm=T)
		SWC[df$TIMESTAMP < t0] <- NA
		
		cp <- (2.31 * meta1$Organic_fraction  +	
					 	2.5 * (meta1$Bulk_density  / 2.6 - meta1$Organic_fraction) +	
					 	4.18 * SWC/100) * 1e6
		
		dTdt1 <- c(diff(TS_1)[1], diff(TS_1)) / (30 * 60)
		dTdt2 <- c(diff(TS_2)[2], diff(TS_2)) / (30 * 60)
		
		x <- data.frame(timestamp = df$TIMESTAMP, dTdt1 = dTdt1, dTdt2 = dTdt2)
		flag <- despikeLF(x,"dTdt1",iter = 3,	light = NULL,	z = 10,	c = 7, var_thr = c(-0.0022,0.0022))
		flag[is.na(flag)] <- 0
		x$dTdt1 <- ifelse(flag > 0, NA, x$dTdt1)
		
		flag <- despikeLF(x,"dTdt2",iter = 3,	light = NULL,	z = 10,	c = 7, var_thr = c(-0.0022,0.0022))
		flag[is.na(flag)] <- 0
		x$dTdt2 <- ifelse(flag > 0, NA, x$dTdt2)
		
		dTdt <- rowMeans(x[, c("dTdt1","dTdt2"), drop = FALSE], na.rm = TRUE)
		dTdt = zoo::na.approx(dTdt, maxgap = 8, na.rm = FALSE)
		
		df$SG_sun_2 <- cp * meta1$Sun_plate_depth_2 * dTdt
	}
	
	#--------------------------------------------------
	# Shade SG 
	#--------------------------------------------------
	if (!all(is.na(meta1[, c("Shade_TS_1_1_1", "Shade_TS_1_1_2","Shade_SWC_1")]))){
		v1 <- meta1$Shade_TS_1_1_1
		v2 <- meta1$Shade_TS_1_1_2
		v3 <- meta1$Shade_SWC_1
		t0 <- meta1$Shade_plate_valid_from_1
		
		TS_1 <- df %>% pull(v1)
		TS_1[df$TIMESTAMP < t0] <- NA
		
		TS_2 <- df %>% pull(v2)
		TS_2[df$TIMESTAMP < t0] <- NA
		
		SWC <- df %>% pull(v3)
		SWC[is.na(SWC)] <- median(SWC,na.rm=T)
		SWC[df$TIMESTAMP < t0] <- NA
		
		cp <- (2.31 * meta1$Organic_fraction  +	
					 	2.5 * (meta1$Bulk_density  / 2.6 - meta1$Organic_fraction) +	
					 	4.18 * SWC/100) * 1e6
		
		dTdt1 <- c(diff(TS_1)[1], diff(TS_1)) / (30 * 60)
		dTdt2 <- c(diff(TS_2)[2], diff(TS_2)) / (30 * 60)
		
		x <- data.frame(timestamp = df$TIMESTAMP, dTdt1 = dTdt1, dTdt2 = dTdt2)
		flag <- despikeLF(x,"dTdt1",iter = 3,	light = NULL,	z = 10,	c = 7, var_thr = c(-0.0022,0.0022))
		flag[is.na(flag)] <- 0
		x$dTdt1 <- ifelse(flag > 0, NA, x$dTdt1)
		
		flag <- despikeLF(x,"dTdt2",iter = 3,	light = NULL,	z = 10,	c = 7, var_thr = c(-0.0022,0.0022))
		flag[is.na(flag)] <- 0
		x$dTdt2 <- ifelse(flag > 0, NA, x$dTdt2)
		
		dTdt <- rowMeans(x[, c("dTdt1","dTdt2"), drop = FALSE], na.rm = TRUE)
		dTdt = zoo::na.approx(dTdt, maxgap = 8, na.rm = FALSE)
		
		df$SG_shade <- cp * meta1$Shade_plate_depth_1 * dTdt
	}
	
	# Sun plate: mean of available sun components
	sun_vars <- c("G_plate_sun_1", "G_plate_sun_2")
	sun_vars <- sun_vars[sun_vars %in% names(df)]
	
	if (length(sun_vars) > 0) {
		df$G_plate_sun <- rowMeans(df[, sun_vars, drop = FALSE], na.rm = TRUE)
		df$G_plate_sun[rowSums(!is.na(df[, sun_vars, drop = FALSE])) == 0] <- NA
	}
	
	# Sun SG: mean of available sun components
	sun_vars <- c("SG_sun_1", "SG_sun_2")
	sun_vars <- sun_vars[sun_vars %in% names(df)]
	
	if (length(sun_vars) > 0) {
		df$SG_sun <- rowMeans(df[, sun_vars, drop = FALSE], na.rm = TRUE)
		df$SG_sun[rowSums(!is.na(df[, sun_vars, drop = FALSE])) == 0] <- NA
	}
	
	
	df$G_sun_calorimetric <- df$SG_sun + df$G_plate_sun
	df$G_shade_calorimetric <- df$SG_shade + df$G_plate_shade

	# # range cut
	# flag <- despikeLF(as.data.frame(df %>% mutate(timestamp = TIMESTAMP)),"G_sun",
	# 									iter = 1,	light = NULL,	z = 10,	c = Inf, var_thr = c(-150,250))
	# flag[is.na(flag)] <- 0
	# df$G_sun <- ifelse(flag > 0, NA, df$G_sun)
	# df$G_sun = zoo::na.approx(df$G_sun, maxgap = 8, na.rm = FALSE)
	# 
	# flag <- despikeLF(as.data.frame(df %>% mutate(timestamp = TIMESTAMP)),"G_shade",
	# 									iter = 1,	light = NULL,	z = 10,	c = Inf, var_thr = c(-150,250))
	# flag[is.na(flag)] <- 0
	# df$G_shade <- ifelse(flag > 0, NA, df$G_shade)
	# df$G_shade = zoo::na.approx(df$G_shade, maxgap = 8, na.rm = FALSE)
	
	df$G_calorimetric <- meta1$G_sun_weight * df$G_sun_calorimetric + (1 - meta1$G_sun_weight) * df$G_shade_calorimetric
	df$SG <- meta1$G_sun_weight * df$SG_sun + (1 - meta1$G_sun_weight) * df$SG_shade
	df$G_plate <- meta1$G_sun_weight * df$G_plate_sun + (1 - meta1$G_sun_weight) * df$G_plate_shade
	
	# save
	site_list[[i]] <- df
	
	# To do: soil heat flux figure save

}

################################################################################
#### Bulk Precipitataion gap filling 
################################################################################
# to do: need to develop


################################################################################
#### Soil surface soil heat flux (without storage term) - harmonic
################################################################################

# #--------------------------------------------------
# # harmonic reconstruction helper (experiment...)
# #--------------------------------------------------
# predict_harmonic_surface <- function(x, hour, zd, depth = 0.08, harmonics = c(1, 2, 4)) {
# 
# 	if (all(is.na(x)) || all(is.na(hour)) || all(is.na(zd))) {
# 		return(rep(NA_real_, length(x)))
# 	}
# 
# 	zd_use <- unique(na.omit(zd))
# 	if (length(zd_use) == 0) return(rep(NA_real_, length(x)))
# 	zd_use <- zd_use[1]
# 
# 	dat <- tibble(
# 		x = x,
# 		hour = hour
# 	) %>%
# 		mutate(
# 			theta1 = 2 * pi * hour / 24,
# 			theta2 = 2 * pi * 2 * hour / 24,
# 			theta4 = 2 * pi * 4 * hour / 24
# 		)
# 
# 	# fit only available rows
# 	fit_dat <- dat %>% filter(!is.na(x), !is.na(hour))
# 	if (nrow(fit_dat) < 10) return(rep(NA_real_, length(x)))
# 
# 	fit <- lm(
# 		x ~ cos(theta1) + sin(theta1) +
# 			cos(theta2) + sin(theta2) +
# 			cos(theta4) + sin(theta4),
# 		data = fit_dat
# 	)
# 
# 	cf <- coef(fit)
# 
# 	get_coef <- function(name) {
# 		ifelse(name %in% names(cf), cf[[name]], 0)
# 	}
# 
# 	c0 <- get_coef("(Intercept)")
# 
# 	a1 <- get_coef("cos(theta1)")
# 	b1 <- get_coef("sin(theta1)")
# 	a2 <- get_coef("cos(theta2)")
# 	b2 <- get_coef("sin(theta2)")
# 	a4 <- get_coef("cos(theta4)")
# 	b4 <- get_coef("sin(theta4)")
# 
# 	# harmonic-specific correction
# 	# if x_plate = A cos(theta - phi), then
# 	# x_surface = A*exp(depth/zd*sqrt(n)) * cos(theta - phi + depth/zd*sqrt(n))
# 	# coefficient form:
# 	# corrected cos/sin coefficients can be obtained by rotating the vector
# 
# 	correct_pair <- function(a, b, n, zd_use, depth) {
# 		amp_factor <- exp(depth / zd_use * sqrt(n))
# 		delta <- depth / zd_use * sqrt(n)
# 
# 		# x_surface(theta) = amp_factor * x_plate(theta + delta)
# 		a_new <- amp_factor * (a * cos(delta) + b * sin(delta))
# 		b_new <- amp_factor * (-a * sin(delta) + b * cos(delta))
# 
# 		list(a = a_new, b = b_new)
# 	}
# 
# 	c1 <- correct_pair(a1, b1, n = 1, zd_use = zd_use, depth = depth)
# 	c2 <- correct_pair(a2, b2, n = 2, zd_use = zd_use, depth = depth)
# 	c4 <- correct_pair(a4, b4, n = 4, zd_use = zd_use, depth = depth)
# 
# 	pred <- dat %>%
# 		mutate(
# 			G0 =
# 				c0 +
# 				c1$a * cos(theta1) + c1$b * sin(theta1) +
# 				c2$a * cos(theta2) + c2$b * sin(theta2) +
# 				c4$a * cos(theta4) + c4$b * sin(theta4)
# 		) %>%
# 		pull(G0)
# 
# 	pred
# }


#--------------------------------------------------
# peak hour estimation function
#--------------------------------------------------
get_peak_hour <- function(x, hour, frac = 0.9, hour_min = 10, hour_max = 20) {
	ok <- !is.na(x) & !is.na(hour) & hour >= hour_min & hour <= hour_max
	
	if (!any(ok)) return(NA_real_)
	
	x2 <- x[ok]
	h2 <- hour[ok]
	
	xmax <- max(x2, na.rm = TRUE)
	if (!is.finite(xmax)) return(NA_real_)
	
	idx <- which(x2 >= frac * xmax)
	if (length(idx) == 0) return(NA_real_)
	
	mean(h2[idx])
}


# Fill damping-depth edges explicitly when only one finite estimate exists.
# No estimate means no reconstruction; do not invent a depth or a zero lag.
level2_fill_damping <- function(x) {
  x[!is.finite(x)] <- NA_real_
  valid <- which(!is.na(x))
  if (!length(valid)) return(rep(NA_real_, length(x)))
  if (length(valid) == 1L) return(rep(x[valid], length(x)))
  zoo::na.approx(x, na.rm = FALSE, rule = 2)
}

level2_smooth_damping <- function(x) {
  x <- level2_fill_damping(x)
  # Incremental windows shorter than 30 days cannot support the original
  # 30-day mean. Keep their interpolated daily estimates without smoothing.
  if (length(x) < 30L || all(is.na(x))) return(x)
  level2_fill_damping(zoo::rollmean(x, k = 30, fill = NA, align = "center"))
}

level2_lead_or_na <- function(x, n) {
  # A missing daily peak/depth leaves the optional Leuning estimate unavailable.
  if (length(n) != 1L || !is.finite(n) || n < 0 || n != floor(n) || n >= length(x))
    return(rep(NA_real_, length(x)))
  dplyr::lead(x, n = as.integer(n))
}

#--------------------------------------------------
# Harmonic based soil heat flux
#--------------------------------------------------

for(i in all_sites) {
	df <- site_list[[i]]
	
	dirs_use <- dirs %>% filter(site == i)
	meta1 <- meta_G %>% filter(site == i)
	
	df <- df %>%
		mutate(
			DATE = as.Date(TIMESTAMP),
			HOUR = hour(TIMESTAMP) + minute(TIMESTAMP)/60,
			MON = month(TIMESTAMP),
			YEAR = year(TIMESTAMP)
		)
	
	### dynamic damping depth model
	dd <- df %>%
		group_by(DATE) %>%
		summarise(
			#peak time of G_plate
			peak_hour_sun = get_peak_hour(G_plate_sun, HOUR),
			peak_hour_shade = get_peak_hour(G_plate_shade, HOUR),
			
			#phase difference
			dt_sun = peak_hour_sun - 12,
			dt_shade = peak_hour_shade - 12,

			dt_sun = ifelse(dt_sun <= 0.5, 0.5, ifelse(dt_sun >= 6, 6,dt_sun)),
			dt_shade = ifelse(dt_shade <= 0.5, 0.5, ifelse(dt_shade >= 6, 6,dt_shade)),

			# damping depth
			zd_sun = 24/(2*pi)/dt_sun*meta1$Sun_plate_depth_1,
			zd_shade = 24/(2*pi)/dt_shade*meta1$Shade_plate_depth_1,
			
			# constrain typical daily damping range
			zd_sun = ifelse(zd_sun >= 0.3, 0.3, ifelse(zd_sun <= 0.1, 0.1, zd_sun)),
			zd_shade = ifelse(zd_shade >= 0.3, 0.3, ifelse(zd_shade <= 0.1, 0.1, zd_shade)),
			
			.groups = "drop"
		) %>%
		mutate(
			# -----------------------------------
			# gap fill by linear interpolation
			# -----------------------------------
			# Preserve the 30-day smoothing where supported; handle short
			# incremental windows and missing daily peaks without losing all rows.
			zd_sun   = level2_smooth_damping(zd_sun),
			zd_shade = level2_smooth_damping(zd_shade),
			
			# recalculate dT
			dt_sun   = 24 / (2 * pi) * meta1$Sun_plate_depth_1   / zd_sun,
			dt_shade = 24 / (2 * pi) * meta1$Shade_plate_depth_1 / zd_shade,
			
			# phase difference hour/30min
			dn_sun   = round(dt_sun / 0.5),
			dn_shade = round(dt_shade / 0.5)
		)
	
	df <- df %>% left_join(dd)
	
	# simple surface soil heat flux model
	df <- df %>%
		group_by(DATE) %>%
		mutate(
			# Method 1: first order (Leuning model)
			G_sun_leuning = level2_lead_or_na(G_plate_sun,first(dn_sun)) * exp(meta1$Sun_plate_depth_1 / zd_sun),
			G_shade_leuning = level2_lead_or_na(G_plate_shade,first(dn_shade)) * exp(meta1$Shade_plate_depth_1 / zd_shade),
			
			# # Method 2: harmonic decomposition
			# G_sun_Kim = predict_harmonic_surface(G_plate_sun, HOUR, zd = zd_sun, depth = meta1$Sun_plate_depth_1),
			# G_shade_Kim = predict_harmonic_surface(G_plate_shade, HOUR, zd = zd_shade, depth = meta1$Sun_plate_depth_1),
			# 
			# edge_flag = HOUR < 1 | HOUR >= 23,
			# G_sun_Kim = ifelse(edge_flag, NA, G_sun_Kim),
			# G_shade_Kim = ifelse(edge_flag, NA, G_shade_Kim),
			
		) %>%
		ungroup() %>%
		mutate(
			G_sun_leuning = na.approx(G_sun_leuning, maxgap = 24, na.rm = FALSE),
			G_shade_leuning = na.approx(G_shade_leuning, maxgap = 24, na.rm = FALSE),
			# G_sun_Kim = na.approx(G_sun_Kim, maxgap = 24, na.rm = FALSE),
			# G_shade_Kim = na.approx(G_shade_Kim, maxgap = 24, na.rm = FALSE),
			
		)
	
	df$G_leuning <- meta1$G_sun_weight * df$G_sun_leuning + (1 - meta1$G_sun_weight) * df$G_shade_leuning
	# df$G_Kim <- meta1$G_sun_weight * df$G_sun_Kim + (1 - meta1$G_sun_weight) * df$G_shade_Kim
	
	# df$G_PI <- rowMeans(df[, c("G_calorimetric","G_leuning"), drop = FALSE], na.rm = TRUE)
	df$G_PI <- df$G_calorimetric
	
	# save
	site_list[[i]] <- df
	
	# To do: soil heat flux figure save
	a <- pipeline_figure_data(df) %>%
		mutate(HOUR = hour(TIMESTAMP) + minute(TIMESTAMP)/60) %>%
		group_by(HOUR) %>%
		summarise(G_calorimetric = mean(G_calorimetric,na.rm=T),
							G_leuning = mean(G_leuning,na.rm=T),
							G_plate = mean(G_plate,na.rm=T),
							G_PI = mean(G_PI,na.rm=T)
							) %>%
		select(HOUR,G_calorimetric,G_leuning,G_plate,G_PI) %>%
		gather(key = "var",value = "G",-HOUR) %>%
		ggplot(aes(HOUR,G, color = var,shape = var)) +
		geom_hline(yintercept = 0) +
		geom_point() +
		geom_line() + theme_bw() +
		labs(y = expression(G~(W~m^-2))) +
		geom_vline(xintercept = 12, lty = 2, alpha = 0.2) +
		ggtitle(i)
	# The soil comparison object is not exported; named diagnostic PDFs above are retained.
	
}




################################################################################
#### save
################################################################################

write_custom_csv <- pipeline_write

# PA_PI_F retains the logger PA scale (kPa); filling does not change its unit.
df_unit <- data.frame(TIMESTAMP = c("TS",NA), TA_PI_F = c('deg C',"derived"), 
											TA_PI_F_method = c("#","derived"), RH_PI_F = c('%',"derived"), 
											RH_PI_F_method = c("#","derived"), PA_PI_F = c('kPa',"derived"), 
											PA_PI_F_method = c("#","derived"), P_PI_F = c('mm',"derived"), 
											P_PI_F_method = c("#","derived"), VPD_PI_F = c("hPa","derived"),
											SW_IN_PI_F = c("W m-2","derived"), SW_IN_PI_F_method = c("#","derived"), 
											LW_IN_PI_F = c("W m-2","derived"), LW_IN_PI_F_method = c("#","derived"), 
											G_plate_sun = c("W m-2","derived"), SG_sun = c("W m-2","derived"), 
											G_plate_shade = c("W m-2","derived"), SG_shade = c("W m-2","derived"),
											G_plate = c("W m-2","derived"), SG = c("W m-2","derived"), 
											G_sun_calorimetric = c("W m-2","derived"), G_shade_calorimetric = c("W m-2","derived"), 
											G_sun_leuning = c("W m-2","derived"), G_shade_leuning = c("W m-2","derived"),
											G_calorimetric = c("W m-2","derived"), G_leuning = c("W m-2","derived"),
											G_PI = c("W m-2","derived")
											)

for(i in all_sites) {
	df <- site_list[[i]]
	dirs_use <- dirs %>% filter(site == i)
	
	df_use <- df %>%
		select(TIMESTAMP, TA_PI_F, TA_PI_F_method, RH_PI_F, RH_PI_F_method, 
					 PA_PI_F, PA_PI_F_method, P_PI_F, P_PI_F_method, VPD_PI_F, SW_IN_PI_F, SW_IN_PI_F_method, 
					 LW_IN_PI_F, LW_IN_PI_F_method, G_plate_sun, SG_sun, G_plate_shade, 
					 SG_shade, G_plate, SG, G_sun_calorimetric, G_shade_calorimetric, 
					 G_sun_leuning, G_shade_leuning, G_calorimetric, G_leuning, G_PI
					 )
	
	print(which(colnames(df_use) != colnames(df_unit)))
	
	write_custom_csv(
		df = df_use,
		units_df = df_unit,
		file = paste0(pipeline_path(base_dir, dirs_use$dir_output),"/",i,"_Level_2_PI_vars_logger_data_",date(Sys.time()),".csv"),
		description = "Gap-filled variables and soil heat flux.Gap filling method represents: 0 = Original, 1 = Short interpolation, 2 = Same-site regression, 3 = Other-site regression, 4 = Final fill / MDC"
	)
	
}


