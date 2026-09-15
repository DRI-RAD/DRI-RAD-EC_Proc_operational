
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
library(REddyProc)

# load direction meta
dirs <- pipeline_options$dirs

# NAS direction for data
base_dir <- pipeline_options$base_dir

################################################################################
#### Define interesting site 
################################################################################
# site id for runing this code
# currently supported site: ECDP EDVG EDVP ERVA ERVP

# single site
# site_id <- "EDVG"
# all sites
for(site_id in pipeline_options$sites){

dirs_use <- dirs %>% filter(site == site_id)

################################################################################
#### Load logger data 
################################################################################
#--------------------------------------------------
# level 1
#--------------------------------------------------
files <- list.files(
  path = pipeline_path(base_dir, dirs_use$dir_output),
  pattern = "level_1_QAQCed_logger_data",
  full.names = TRUE,
  ignore.case = TRUE
)


latest_file <- pipeline_latest(files)

met_unit <- read_csv(latest_file, skip = 1, n_max = 2)
met <- read_csv(latest_file, skip = 4, col_names = FALSE)
colnames(met) <- colnames(met_unit)
# Restrict input to the requested calculation window before processing.
met <- pipeline_slice(met)

# PotRad
met <- met %>% select(TIMESTAMP, PotRad)
met_unit <- met_unit %>% select(TIMESTAMP, PotRad)


################################################################################
#### Load LI710 
################################################################################

if(!is.na(dirs_use$dir_LI710_old)){ 
  ### Load LI710 old data
  LI710_unit_old <- read_delim(pipeline_path(base_dir, dirs_use$dir_LI710_old),skip = 1, n_max = 2)
  LI710_old <- read_delim(pipeline_path(base_dir, dirs_use$dir_LI710_old),skip = 4, col_names = F)
  colnames(LI710_old) <- colnames(LI710_unit_old)
  # Restrict input to the requested calculation window before processing.
  LI710_old <- pipeline_slice(LI710_old)
  
  ### Load LI710 data
  LI710_unit <- read_delim(pipeline_path(base_dir, dirs_use$dir_LI710),skip = 1, n_max = 2)
  LI710 <- read_delim(pipeline_path(base_dir, dirs_use$dir_LI710),skip = 4, col_names = F)
  colnames(LI710) <- colnames(LI710_unit)
  # Restrict input to the requested calculation window before processing.
  LI710 <- pipeline_slice(LI710)
  
  LI710 <- full_join(LI710_old, LI710)
  LI710_unit <- full_join(LI710_unit_old, LI710_unit)
  
} else{
  ### Load LI710 data
  LI710_unit <- read_delim(pipeline_path(base_dir, dirs_use$dir_LI710),skip = 1, n_max = 2)
  LI710 <- read_delim(pipeline_path(base_dir, dirs_use$dir_LI710),skip = 4, col_names = F)
  colnames(LI710) <- colnames(LI710_unit)  
  # Restrict input to the requested calculation window before processing.
  LI710 <- pipeline_slice(LI710)
  
}

LI710[LI710==9999999] <- NA
LI710[LI710==-9999] <- NA

if(!"data_qc" %in% colnames(LI710)){
  LI710$data_qc <- 0
  LI710_unit$data_qc <- NA
}

LI710 <- LI710 %>%
  select(c('TIMESTAMP','LE_710','H_710','diag','flow','tilt','data_qc')) %>%
  distinct(TIMESTAMP, .keep_all = TRUE) %>% 
  mutate(TIMESTAMP_30 = floor_date(TIMESTAMP, "30 minutes")) %>%
  group_by(TIMESTAMP_30) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  rename(TIMESTAMP = TIMESTAMP_30)

LI710_unit <- LI710_unit %>% select(c('TIMESTAMP','LE_710','H_710','diag','flow','tilt','data_qc'))

colnames(LI710) <- c('TIMESTAMP','LE_710','H_710','diag_710','flow_710','tilt_710','data_qc_710')
colnames(LI710_unit) <- c('TIMESTAMP','LE_710','H_710','diag_710','flow_710','tilt_710','data_qc_710')

# li710 observations until 2024 May is unconvincing 
if(site_id =="EDVG") {
	LI710 <- LI710 %>%
		filter(TIMESTAMP > as.POSIXct("2024-10-26 00:00"))
	
}

# left join to HH data
LI710 <- left_join(met,LI710)
LI710_unit <- left_join(met_unit,LI710_unit)


################################################################################
#### QC LI710 column
################################################################################
# flow: 0 good, 1 outside range
# diag: 0 good, 1 bad
# desipke: 0 good, 2 spike, 1 physical range out

LI710 <- LI710 %>%
	mutate(LE_710_QC_diag = NA,
				 LE_710_QC_despike = NA, # despike + physical range
				 LE_710_QC_longrun = NA,
				 LE_710_QC_flow = NA,
				 
				 H_710_QC_diag = NA, # diag based filtering
				 H_710_QC_despike = NA, # despike + physical range
				 H_710_QC_longrun = NA,
				 H_710_QC_flow = NA
				 )



################################################################################
#### Despike and physical range, fetch filtering, and filter long-run
################################################################################
# physical range 
if(site_id %in% c('EDVG', 'ERVA')){
	# physical range for agriculture
	LE_range = c(-30, 1100)
	H_range = c(-300, 1100)
	co2_flux_range = c(-60,30)
	
} else{
	# physical range for phreatophytes
	LE_range = c(-30, 600)
	H_range = c(-100, 1100)
	co2_flux_range = c(-30,10)
}

# despike + physical range filtering
LI710$LE_710_QC_despike <- despikeLF(as.data.frame(LI710 %>% mutate(timestamp = TIMESTAMP)),
															"LE_710",iter = 10,light = "PotRad", z = 10, var_thr = LE_range)
LI710$H_710_QC_despike <- despikeLF(as.data.frame(LI710 %>% mutate(timestamp = TIMESTAMP)),
														 "H_710",iter = 10,light = "PotRad", z = 10, var_thr = H_range)


# filterLongRuns (replace runs, i.e sequences of numerically equal values, by NA)
LI710_longrun <- LI710 %>%
	select(LE_710,H_710) %>%
	filterLongRuns(c("LE_710","H_710"), minNRunLength = 8)

LI710$LE_710_QC_longrun <- ifelse(is.na(LI710_longrun$LE_710) & !is.na(LI710$LE_710),1,0)
LI710$H_710_QC_longrun <- ifelse(is.na(LI710_longrun$H_710) & !is.na(LI710$H_710),1,0)

# flow rate (Average flow for a 30-minute period is less than 125 sccm or greater than 330 sccm) 
LI710$LE_710_QC_flow <- ifelse((LI710$flow_710 < 125 | LI710$flow_710 > 330), 1, 0)
LI710$H_710_QC_flow <- 0

# Diagnostics
# Code	Description
# 1	Average flow for a 30-minute period is less than 125 sccm or greater than 330 sccm
# 8	Voltage is less than or equal to 1.6 volts for more than 50 % of time for a 30-minute period
# 16	Temperature is greater than 65 °C or less than −50 °C for more than 5 % of time for a 30-minute period
# 128	Poor sonic signals persist for more than 10% of time for a 30-minute period
# 512	High humidity shutdown
# 1024	Cold temperature shutdown

has_bad_diag <- function(diag, bad_bits = c(1, 8, 16, 128, 512, 1024)) {
  rowSums(sapply(bad_bits, function(x) bitwAnd(diag, x) > 0)) > 0
}

LI710$data_qc_710[is.na(LI710$data_qc_710)] <- 0

LI710$LE_710_QC_diag <- ifelse(has_bad_diag(LI710$diag_710) & LI710$data_qc_710 > 15, 1, 0)
LI710$H_710_QC_diag <- ifelse(has_bad_diag(LI710$diag_710) & LI710$data_qc_710 > 15, 1, 0)



################################################################################
#### Visualization
################################################################################

plot_flux_qc_pdf <- function(df,
														 var = "LE_710",
														 site_id = "",
														 file_out = "",
														 n_panel_per_page = 6,
														 months_per_panel = 2) {
	
	stopifnot("TIMESTAMP" %in% names(df))
	stopifnot(var %in% names(df))
	
	#--------------------------------------------------
	# QC column names
	#--------------------------------------------------
	qc_diag    <- paste0(var, "_QC_diag")
	qc_despike <- paste0(var, "_QC_despike")
	qc_longrun <- paste0(var, "_QC_longrun")
	qc_flow <- paste0(var, "_QC_flow")
	
	qc_cols <- c(qc_diag, qc_despike, qc_longrun, qc_flow)
	missing_qc <- qc_cols[!qc_cols %in% names(df)]
	
	if (length(missing_qc) > 0) {
		stop("Missing QC columns: ", paste(missing_qc, collapse = ", "))
	}
	
	#--------------------------------------------------
	# prepare
	#--------------------------------------------------
	df <- df %>%
		dplyr::filter(!is.na(TIMESTAMP)) %>%
		dplyr::arrange(TIMESTAMP)
	
	df_good <- df %>%
		dplyr::filter(
			.data[[qc_diag]] == 0,
			.data[[qc_despike]] == 0,
			.data[[qc_longrun]] == 0,
			.data[[qc_flow]] == 0
		)
	
	if (nrow(df_good) == 0 || all(is.na(df_good[[var]]))) {
		stop("No valid good-quality data found for ", var)
	}
	
	y_rng <- range(df_good[[var]], na.rm = TRUE)
	
	#--------------------------------------------------
	# time breaks
	#--------------------------------------------------
	t_start <- floor_date(min(df$TIMESTAMP), "month")
	t_end   <- ceiling_date(max(df$TIMESTAMP), "month")
	breaks  <- seq(t_start, t_end + months(months_per_panel), 
								 by = paste(months_per_panel, "months"))
	
	#--------------------------------------------------
	# y label
	#--------------------------------------------------
	y_lab <- switch(
		var,
		"LE_710" = "Latent heat flux (LE_710)",
		"H_710" = "Sensible heat flux (H_710)",
		var
	)
	
	#--------------------------------------------------
	# generate plots
	#--------------------------------------------------
	all_plots <- list()
	k <- 1
	
	for (i in seq_len(length(breaks) - 1)) {
		
		s <- breaks[i]
		e <- breaks[i + 1]
		
		df_sub <- df %>%
			dplyr::filter(TIMESTAMP >= s, TIMESTAMP < e)
		
		if (nrow(df_sub) == 0) next
		
		df_sub_good <- df_sub %>%
			dplyr::filter(
				.data[[qc_diag]] == 0,
				.data[[qc_despike]] == 0,
				.data[[qc_longrun]] == 0,
				.data[[qc_flow]] == 0
			)
		
		p <- ggplot2::ggplot(df_sub_good, ggplot2::aes(x = TIMESTAMP, y = .data[[var]])) +
			ggplot2::theme_bw() +
			ggplot2::geom_line(na.rm = TRUE) +
			ggplot2::geom_point(size = 1, na.rm = TRUE) +
			ggplot2::geom_point(
					ggplot2::aes(y = .data[[var]], color = "diag"),
					data = df_sub %>% dplyr::filter(.data[[qc_diag]] >= 1),
					size = 1,
					na.rm = TRUE
				) +	
			ggplot2::geom_point(
				ggplot2::aes(y = .data[[var]], color = "despike & physical range"),
				data = df_sub %>% dplyr::filter(.data[[qc_despike]] >= 1),
				size = 1,
				na.rm = TRUE
			) +
			ggplot2::geom_point(
				ggplot2::aes(y = .data[[var]], color = "longrun"),
				data = df_sub %>% dplyr::filter(.data[[qc_longrun]] >= 1),
				size = 1,
				na.rm = TRUE
			) +
			ggplot2::geom_point(
				ggplot2::aes(y = .data[[var]], color = "flow"),
				data = df_sub %>% dplyr::filter(.data[[qc_flow]] >= 1),
				size = 1,
				na.rm = TRUE
			) +
			ggplot2::coord_cartesian(ylim = y_rng) +
			ggplot2::scale_color_manual(
				values = c(
					"diag" = "purple3",
					"despike & physical range" = "orange2",
					"longrun" = "gray",
					"flow" = "dodgerblue3"
				),
				drop = FALSE
			) +
			ggplot2::labs(
				title = paste0(
					site_id, " | ", var, " | ",
					format(s, "%Y-%m"), " to ",
					format(e - lubridate::seconds(1), "%Y-%m")
				),
				x = NULL,
				y = y_lab,
				caption = paste(
					"Black points indicate data passing all QC filters.",
					"Colored points represent data flagged by individual QC tests:",
					"flow, flag, despiking/physical range, and longrun."
				)
			) +
			ggplot2::theme(
				plot.title = ggplot2::element_text(face = "bold", size = 10),
				plot.caption = ggplot2::element_text(size = 7, hjust = 0),
				axis.text.x = ggplot2::element_text(size = 7),
				axis.text.y = ggplot2::element_text(size = 7),
				axis.title.y = ggplot2::element_text(size = 8),
				panel.grid.minor = ggplot2::element_blank(),
				plot.margin = ggplot2::margin(3, 6, 3, 3)
			)
		
		all_plots[[k]] <- p
		k <- k + 1
	}
	
	if (length(all_plots) == 0) {
		stop("No plots were generated.")
	}
	
	#--------------------------------------------------
	# split pages
	#--------------------------------------------------
	page_id <- ceiling(seq_along(all_plots) / n_panel_per_page)
	pages <- split(all_plots, page_id)
	
	grDevices::pdf(file_out, width = 16, height = 20, onefile = TRUE, useDingbats = FALSE)
	
	for (i in seq_along(pages)) {
		
		plotlist_i <- pages[[i]]
		
		if (length(plotlist_i) < n_panel_per_page) {
			plotlist_i <- c(
				plotlist_i,
				rep(list(NULL), n_panel_per_page - length(plotlist_i))
			)
		}
		
		page_plot <- cowplot::plot_grid(
			plotlist = plotlist_i,
			nrow = n_panel_per_page,
			align = "v"
		)
		
		print(page_plot)
	}
	
	grDevices::dev.off()
	
	message("Saved: ", file_out)
}


plot_flux_qc_pdf(LI710, var = "LE_710", 
								 file_out = paste0(pipeline_path(base_dir, dirs_use$dir_output),"/figures/",site_id,"_Level_3_LE_710_flag.pdf")
)
plot_flux_qc_pdf(LI710, var = "H_710", 
								 file_out = paste0(pipeline_path(base_dir, dirs_use$dir_output),"/figures/",site_id,"_Level_3_H_710_flag.pdf")
)


################################################################################
#### save
################################################################################

write_custom_csv <- pipeline_write


LI710_unit_final <- cbind(LI710_unit,
                       data.frame(LE_710_QC_diag = c("#","derived"), LE_710_QC_despike = c("#","derived"),LE_710_QC_longrun = c("#","derived"),LE_710_QC_flow = c("#","derived"),
                                  H_710_QC_diag = c("#","derived"), H_710_QC_despike = c("#","derived"),H_710_QC_longrun = c("#","derived"),H_710_QC_flow = c("#","derived")
                       ))

print(which(colnames(LI710) != colnames(LI710_unit_final)))

write_custom_csv(
	df = LI710 %>% select(-PotRad),
	units_df = LI710_unit_final %>% select(-PotRad),
	file = paste0(pipeline_path(base_dir, dirs_use$dir_output),"/",site_id,"_Level_3_li710_with_qaqc_",date(Sys.time()),".csv"),
	description = "LI710 flux data and QAQC flag"
)

}
