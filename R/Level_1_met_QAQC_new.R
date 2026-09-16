# Run this stage through run_pipeline.R so period and output state are explicit.
if (is.null(getOption("ec.pipeline"))) stop("Use run_pipeline() from run_pipeline.R.")
pipeline_options <- getOption("ec.pipeline")

################################################################################
#### Load library and meta data
################################################################################
library(tidyverse)
library(lubridate)
library(zoo)
library(openeddy)
library(cowplot)
library(REddyProc)

# load direction meta
dirs <- pipeline_options$dirs

# load manual flag data 
manual_flag <- read_csv(file.path(.pipeline_root, "meta_files", "manual_flag.csv"), skip = 2) %>% filter(flag_type == "remove")
manual_flag$start <- as.POSIXct(manual_flag$start, format = "%m/%d/%Y %H:%M", tz = "UTC")
manual_flag$end <- as.POSIXct(manual_flag$end, format = "%m/%d/%Y %H:%M", tz = "UTC")

# NAS direction for data
base_dir <- pipeline_options$base_dir

################################################################################
#### Define interesting site and load meta
################################################################################
# site id for runing this code
# currently supported site: ECDP EDVG EDVP ERVA ERVP ECSM

# single site
# The runner supplies the selected site.
# all sites
for(site_id in pipeline_options$sites){
	

dirs_use <- dirs %>% filter(site == site_id)

################################################################################
#### Load CSFormat data
################################################################################
CSFormat_unit <- read_delim(pipeline_path(base_dir, dirs_use$dir_met),skip = 1, n_max = 2)
CSFormat <- read_delim(pipeline_path(base_dir, dirs_use$dir_met),skip = 4, col_names = F)
colnames(CSFormat) <- colnames(CSFormat_unit)
# Restrict input to the requested calculation window before processing.
CSFormat <- pipeline_slice(CSFormat)

# Full time sequence
CSFormat <- CSFormat %>%
	distinct(TIMESTAMP, .keep_all = TRUE)
TIMESTAMP <- seq(from = min(CSFormat$TIMESTAMP), to = max(CSFormat$TIMESTAMP), 60*30)
CSFormat <- left_join(data.frame(TIMESTAMP), CSFormat)

# remove Easyflux output etc.
need_to_remove <- 
	c('RECORD', 'FC_mass', 'FC_QC', 'FC_samples', 'LE', 'LE_QC', 'LE_samples', 'H', 
		'H_QC', 'H_samples', 'G', 'SG', 'energy_closure', 'poor_enrg_clsur', 'Bowen_ratio', 
		'TAU', 'TAU_QC', 'USTAR', 'TSTAR', 'TKE', 'e_amb', 'e_sat_amb', 'e', 'e_sat', 
		'e_probe', 'e_sat_probe', 'H2O_density_probe', 'VPD', 'Ux', 'Ux_SIGMA', 'Uy', 
		'Uy_SIGMA','Uz', 'Uz_SIGMA', 'T_SONIC', 'T_SONIC_SIGMA', 'sonic_azimuth', 'WS', 
		'WS_RSLT','WD_SONIC', 'WD_SIGMA', 'WD', 'WS_MAX', 'CO2_density', 'CO2_density_SIGMA',
		'H2O_density', 'H2O_density_SIGMA', 'G_1_1_1', 'G_1_1_2', 'G_1_1_3', 'G_1_1_4', 
		'G_1_1_5', 'G_1_1_6', 'SG_1_1_1', 'SG_1_1_2', 'SG_1_1_3', 'FETCH_MAX', 'FETCH_90', 
		'FETCH_55', 'FETCH_40', 'UPWND_DIST_INTRST', 'FP_DIST_INTRST', 'FP_EQUATION', 
		'slowsequence_count', 'slow_main_count')

need_to_remove <- need_to_remove[which(need_to_remove %in% colnames(CSFormat))]

CSFormat <- CSFormat %>%
	select(!need_to_remove)
CSFormat_unit <- CSFormat_unit %>%
	select(!need_to_remove)
head(CSFormat)
head(CSFormat_unit)

df <- CSFormat
df_unit <- CSFormat_unit

if(site_id =="EDVG"){
	df <- df %>%
		filter(TIMESTAMP >= as.POSIXct("2023-09-26 11:30:00", tz = "UTC"))
}

################################################################################
#### Load LI710 (if available) - merge met variable only 
################################################################################
if(!is.na(dirs_use$dir_LI710)){
  
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
  
  # Full time sequence for 30 min
	LI710 <- LI710 %>%
		select(c('TIMESTAMP','Atm_press','AirTC','RH_710')) %>%
		distinct(TIMESTAMP, .keep_all = TRUE) %>% 
		mutate(TIMESTAMP_30 = floor_date(TIMESTAMP, "30 minutes")) %>%
		group_by(TIMESTAMP_30) %>%
		summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
							.groups = "drop") %>%
		rename(TIMESTAMP = TIMESTAMP_30)
	
	LI710_unit <- LI710_unit %>% select(c('TIMESTAMP','Atm_press','AirTC','RH_710'))
	
	colnames(LI710) <- c('TIMESTAMP','PA_710','TA_710','RH_710')
	colnames(LI710_unit) <- c('TIMESTAMP','PA_710','TA_710','RH_710')
	
	# left join to HH data
	df <- left_join(df, LI710)
	df_unit <- left_join(df_unit, LI710_unit)
}

################################################################################
#### Load Teros and merge (if available)
################################################################################
if(!is.na(dirs_use$dir_Teros)){
	### Load Teros data
	Teros_unit <- read_delim(pipeline_path(base_dir, dirs_use$dir_Teros),skip = 1, n_max = 2)
	Teros <- read_delim(pipeline_path(base_dir, dirs_use$dir_Teros),skip = 4, col_names = F)
	colnames(Teros) <- colnames(Teros_unit)
	# Restrict input to the requested calculation window before processing.
	Teros <- pipeline_slice(Teros)
	
	# Full time sequence for 30 min
	Teros_30min <- Teros %>%
		select(!RECORD) %>%
		distinct(TIMESTAMP, .keep_all = TRUE) %>% 
		mutate(TIMESTAMP_30 = floor_date(TIMESTAMP, "30 minutes")) %>%
		group_by(TIMESTAMP_30) %>%
		summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
							.groups = "drop") %>%
		rename(TIMESTAMP = TIMESTAMP_30)
	
	Teros_unit <- Teros_unit %>% select(!RECORD)
	
	# left join to HH data
	df <- left_join(df, Teros_30min)
	df_unit <- left_join(df_unit, Teros_unit)
}

################################################################################
#### Load bulk precipitation data, process 
#### and merge (if available)
################################################################################
# bulk precipitation (Pluvio or SG400)
if(!is.na(dirs_use$dir_BulkP)){
	### Load precip data
	BulkP_unit <- read_delim(pipeline_path(base_dir, dirs_use$dir_BulkP),skip = 1, n_max = 2)
	BulkP <- read_delim(pipeline_path(base_dir, dirs_use$dir_BulkP),skip = 4, col_names = F)
	colnames(BulkP) <- colnames(BulkP_unit)
	# Restrict input to the requested calculation window before processing.
	BulkP <- pipeline_slice(BulkP, predecessor = TRUE)
	
	# Full time sequence for 30 min
	if(site_id=='ERVP'){
		BulkP <- BulkP %>%
			mutate(DATE = as.Date(TIMESTAMP)) %>%
			filter(DATE > as.Date("2024-12-17")) %>%
			select(!DATE)
	} 
	# Full time sequence for 30 min
	BulkP <- BulkP %>%
		distinct(TIMESTAMP, .keep_all = TRUE)
	TIMESTAMP <- seq(from = min(BulkP$TIMESTAMP), to = max(BulkP$TIMESTAMP), 60*30)
	BulkP <- left_join(data.frame(TIMESTAMP), BulkP)
	
	BulkP <- BulkP %>%
		mutate(
			DATE = as.Date(TIMESTAMP),
			YEAR = year(TIMESTAMP))
	
	# calculate raw precipitation from cumlative
	if(dirs_use$BulkP_type == "Pluvio"){
		# precipitation calculation
		BulkP <- BulkP %>%
			mutate(
				cum_signal = Accu_Total_NRT,
				P_raw  = cum_signal - lag(cum_signal),
				P_raw  = ifelse(P_raw < -4, 0, P_raw),  # remove weight drop
				P_raw = replace_na(P_raw, 0),
				P_cumsum_raw = cumsum(P_raw)
			)
	} 
	
	if(dirs_use$BulkP_type == "SG400"){
		# precipitation calculation
		BulkP <- BulkP %>%
			mutate(
				cum_signal = SG400_Pc,
				P_raw  = cum_signal - lag(cum_signal),
				P_raw  = ifelse(P_raw < -10, 0, P_raw),  # remove weight drop
				P_raw = replace_na(P_raw, 0),
				P_cumsum_raw = cumsum(P_raw)
			)
	} 
	
	# function for daily P threshold
	find_best_threshold_by_year <- function(daily_tbl, annual_target) {
		years <- sort(unique(daily_tbl$YEAR))
		
		res <- lapply(years, function(y) {
			x <- daily_tbl %>% filter(YEAR == y)
			target <- annual_target %>% filter(YEAR == y) %>% pull(raw_total)
			
			candidates <- sort(unique(c(0, x$day_net)))
			
			eval_tbl <- tibble(
				threshold = candidates,
				processed_total = sapply(candidates, function(thr) {
					sum(x$day_net[x$day_net >= thr], na.rm = TRUE)
				})
			) %>%
				mutate(
					target = target,
					abs_diff = abs(processed_total - target)
				) %>%
				arrange(abs_diff, threshold)
			
			eval_tbl[1, ] %>% mutate(YEAR = y)
		})
		
		bind_rows(res) %>%
			select(YEAR, threshold, processed_total, target, abs_diff)
	}
	##
	daily_tbl <- BulkP %>%
		group_by(YEAR, DATE) %>%
		summarise(
			day_net = pmax(last(P_cumsum_raw[!is.na(P_cumsum_raw)]) - first(P_cumsum_raw[!is.na(P_cumsum_raw)]), 0),
			.groups = "drop"
		)
	
	annual_target <- BulkP %>%
		group_by(YEAR) %>%
		summarise(
			raw_total = sum(P_raw, na.rm = TRUE),
			.groups = "drop"
		)
	
	best_thr_year <- find_best_threshold_by_year(daily_tbl, annual_target)
	print(best_thr_year)
	n_peak <- 3

	BulkP <- BulkP %>%
		left_join(best_thr_year %>% select(YEAR,threshold)) %>%
		group_by(DATE) %>%
		group_modify(~{
			df <- .x
			
			day_net <- pmax(
				last(df$P_cumsum_raw[!is.na(df$P_cumsum_raw)]) -
					first(df$P_cumsum_raw[!is.na(df$P_cumsum_raw)]),
				0
			)
			
			thr <- df$threshold[1]
			p <- ifelse(df$P_raw > 0, df$P_raw, 0)
			out <- rep(0, length(p))
			
			if (day_net < thr || sum(p, na.rm = TRUE) <= 0) {
				df$day_net <- day_net
				df$pos_sum <- sum(p, na.rm = TRUE)
				df$P_Bulk <- 0
				return(df)
			}
			
			# top n_peak index
			ord <- order(p, decreasing = TRUE)
			peak_idx <- ord[seq_len(min(n_peak, sum(p > 0)))]
			
			peak_sum <- sum(p[peak_idx], na.rm = TRUE)
			
			if (peak_sum >= day_net) {
				out[peak_idx] <- p[peak_idx] * day_net / peak_sum
				
			} else {
				out[peak_idx] <- p[peak_idx]
				
				remain_net <- day_net - peak_sum
				remain_idx <- setdiff(which(p > 0), peak_idx)
				remain_sum <- sum(p[remain_idx], na.rm = TRUE)
				
				if (remain_net > 0 && remain_sum > 0) {
					out[remain_idx] <- p[remain_idx] * remain_net / remain_sum
				}
			}
			
			df$day_net <- day_net
			df$pos_sum <- sum(p, na.rm = TRUE)
			df$P_Bulk <- out
			df
		}) %>%
		ungroup() %>%
		mutate(
			P_cumsum = cumsum(P_Bulk)
		)
	
	BulkP <- BulkP %>%
		mutate(P_raw = ifelse(is.na(cum_signal),NA, P_raw),
					 P_Bulk = ifelse(is.na(cum_signal),NA, P_Bulk),
					 P_cumsum_raw = ifelse(is.na(cum_signal),NA, P_cumsum_raw),
					 P_cumsum = ifelse(is.na(cum_signal),NA, P_cumsum),
					 )
	
	# figure
	a <- pipeline_figure_data(BulkP) %>%
		ggplot(aes(TIMESTAMP)) +
		theme_bw() +
		geom_point(aes(y = P_raw, color = "raw"), size = 2) +
		geom_point(aes(y = P_Bulk, color = "corrected"), size = 1) +
		labs(title = "Bulk Precipitation",
				 y = "Bulk P (mm/30min)"
				 ) +
		scale_color_manual(values = c("black","red"))
	
	b <- pipeline_figure_data(BulkP) %>%
		ggplot(aes(TIMESTAMP)) +
		theme_bw() +
		geom_line(aes(y = cum_signal, color = "sensor output"), linewidth = 1, lty = 2) +
		geom_line(aes(y = P_cumsum_raw, color = "raw"), linewidth = 2) +
		geom_line(aes(y = P_cumsum, color = "corrected"), linewidth = 1) +
				labs(title = "P Accumulation",
				 y = "P Accumulation (mm)"
		)+
		scale_color_manual(values = c("black","red","orange"))
	
	p <- plot_grid(a,b,ncol = 1)
	
	if(dirs_use$BulkP_type == "SG400"){
	ggsave(plot = p,
				 filename = paste0(pipeline_path(base_dir, dirs_use$dir_output),"/figures/",site_id,"_Level_1_Figure_Bulk_P.png"),
				 width = 12, height = 8, dpi = 150
				 )
	}

	
	# Pluvio do not need correction.
	if(dirs_use$BulkP_type == "Pluvio"){
	  BulkP <- BulkP %>% select(TIMESTAMP, P_raw)
	  colnames(BulkP) <- c("TIMESTAMP","P_Bulk")
	  BulkP_unit <- BulkP_unit %>% select(TIMESTAMP)
	  BulkP_unit$P_Bulk <- c("mm",NA)
	} 
	
	# SG400 need correction
	if(dirs_use$BulkP_type == "SG400"){
	  BulkP <- BulkP %>% select(TIMESTAMP, P_Bulk)
	  colnames(BulkP) <- c("TIMESTAMP","P_Bulk")
	  BulkP_unit <- BulkP_unit %>% select(TIMESTAMP)
	  BulkP_unit$P_Bulk <- c("mm",NA)
	} 
	
	# left join to merged data
	df <- left_join(df, BulkP)
	df_unit <- left_join(df_unit, BulkP_unit)
}

#################################################
#### gap statistics
#################################################
df_orig <- df

# gap stat
# To do: need to produce each year?
df_stat <- data.frame(vars = colnames(df)[2:ncol(df)], 
											total = nrow(df),
											gap = as.numeric(colSums(is.na(df[, 2:ncol(df)]))),
											gap_after_QC = NA) 

#################################################
#### QC 0: filterLongRuns
#################################################
# replace runs, i.e sequences of numerically equal values, by NA
# To do: under development
# need to select most variables except for precipitation. will update in next round

vars_for_longrun_filter <- 
  c("SW_IN","SW_OUT","PPFD_IN","LW_IN","LW_OUT","T_nr",
    "TA_1_1_1", "TA_1_1_2","TA_1_1_3","RH_1_1_1","RH_1_1_2","RH_1_1_3",
    "PA","G_plate_1_1_1","G_plate_1_1_2","G_plate_1_1_3","G_plate_1_1_4",
    "G_plate_1_1_5","G_plate_1_1_6","TS_1_1_1","TS_1_1_2","TS_1_1_3",
    "TS_CS65X_1_1_1","TS_CS65X_1_1_2","TS_CS65X_1_1_3","SWC_1_1_1","SWC_1_1_2","SWC_1_1_3",
    "T12_1_Temp_Avg", "T12_1_VWC_Avg", "T12_1_EC_Avg", "T12_2_Temp_Avg", "T12_2_VWC_Avg", 
    "T12_2_EC_Avg", "T12_3_Temp_Avg", "T12_3_VWC_Avg", "T12_3_EC_Avg", "T12_4_Temp_Avg", 
    "T12_4_VWC_Avg",  "T12_4_EC_Avg", "T54_1_Temp_Avg", "T54_1_VWC_Avg", "T54_2_Temp_Avg", 
    "T54_2_VWC_Avg", "T54_3_Temp_Avg", "T54_3_VWC_Avg", "T54_4_Temp_Avg", "T54_4_VWC_Avg" 
  )

vars_for_longrun_filter <- vars_for_longrun_filter[which(vars_for_longrun_filter %in% colnames(df))]

df <- df %>%
  filterLongRuns(colNames =vars_for_longrun_filter, minNRunLength = 48*2)

#################################################
#### QC 1: Radiation
#################################################
# To do: implement clear-sky radiation (Rso) to constrain the upper bound and assess potential bias.
# df <- df %>%
# 	mutate(SW_IN = ifelse(SW_IN > Rso + 50, NA, SW_IN))

# Occasional spikes were observed in PPFD_IN while SW_IN remained stable.
# Such cases are considered non-physical, and PPFD_IN is removed when it exceeds
# a threshold relative to SW_IN based on an empirical conversion factor.

# df <- df %>%
# 	mutate(PPFD_IN = ifelse(PPFD_IN > (SW_IN + 100) / 0.5, NA, PPFD_IN)) #conversion factor should be 0.48~0.5x

# Shortwave radiation: remove non-physical values (negative, SW_OUT > SW_IN).
# Longwave radiation: apply T_nr-based correction when T_nr is inconsistent with air temperature.
# Net radiation components are subsequently adjusted accordingly.

df <- df %>%
	mutate(PPFD_IN = ifelse(PPFD_IN < 0, 0, PPFD_IN),
				 SW_IN = ifelse(SW_IN < 0, 0, SW_IN),
				 SW_OUT = ifelse(SW_IN == 0, 0, 
				 								ifelse(SW_OUT < 0, 0, 
				 											 ifelse(SW_OUT > SW_IN, SW_IN, SW_OUT))),
				 T_nr_error = abs(T_nr - 273.15 - TA_1_1_3) > 15,
				 T_nr = ifelse(T_nr_error, TA_1_1_3 + 273.15, T_nr),
				 LW_IN = ifelse(T_nr_error, R_LW_in_meas + 5.6718e-8 * T_nr^4, LW_IN),
				 LW_OUT = ifelse(T_nr_error, R_LW_out_meas + 5.6718e-8 * T_nr^4, LW_OUT),
				 NETRAD = SW_IN - SW_OUT + LW_IN - LW_OUT,
				 ALB = ifelse(SW_IN == 0 , 0, SW_OUT/SW_IN * 100)
	) %>% 
	select(!T_nr_error)


#################################################
#### QC 2: Manual flag (details in the manual_flag.csv file)
#################################################
# select manual flag for the site
manual_flag_use <- manual_flag %>% filter(site == site_id)

# Manual flag QA/QC
apply_manual_flag <- function(df, flag_tbl, time_col = "TIMESTAMP") {
	
	for (i in seq_len(nrow(flag_tbl))) {
		
		# 1. flag datetime
		row <- flag_tbl[i, ]
		start_time <- row$start
		end_time   <- row$end
		
		# 2. vars parsing
		vars <- strsplit(row$vars, ",\\s*")[[1]]
		vars <- vars[vars %in% names(df)]
		
		if (length(vars) == 0) next
		
		# 3. NA
		idx <- df[[time_col]] >= start_time & df[[time_col]] <= end_time
		df[idx, vars] <- NA
	}
	
	return(df)
}

df <- apply_manual_flag(df, manual_flag_use)

################################################################################
#### QC 3: Despike (MAD approach) + physical range
################################################################################
# Apply a conservative MAD-based despiking (despikeLF).
# This method is sensitive to abrupt changes and non-stationary behavior,
# therefore it is NOT applied to variables with strong physical variability
# or discontinuities (e.g., precipitation, wind direction, shortwave radiation).


calc_box_range <- function(x,k = 3) {
	x <- x[!is.na(x)]
	q1 <- quantile(x, 0.25)
	q3 <- quantile(x, 0.75)
	iqr <- q3 - q1
	
	lower <- q1 - k * iqr
	upper <- q3 + k * iqr
	
	return(as.numeric(c(lower, upper)))
}

despike_vars <- c("TA_1_1_1", "T_DP_1_1_1",
									"TA_1_1_2", "T_DP_1_1_2",
									"TA_1_1_3", "T_DP_1_1_3",	"PA",
									"TS_1_1_1","TS_1_1_2","TS_1_1_3")

for (v in despike_vars) {
	if (!v %in% names(df)) next
	
	var_thr_use <- calc_box_range(df[[v]], k = 3)
	
	flag <- despikeLF(
		as.data.frame(df %>% mutate(timestamp = TIMESTAMP)),
		v,
		iter = 10,
		light = NULL,
		z = 10,
		c = 7,
		var_thr = var_thr_use
	)
	
	flag[is.na(flag)] <- 0
	df[[v]] <- ifelse(flag > 0, NA, df[[v]])
}

range_only_vars <- c("SW_IN","SW_OUT","PPFD_IN","LW_IN","LW_OUT")

for (v in range_only_vars) {
	if (!v %in% names(df)) next
	
	var_thr_use <- calc_box_range(df[[v]], k = 4)
	
	flag <- despikeLF(
		as.data.frame(df %>% mutate(timestamp = TIMESTAMP)),
		v,
		iter = 1,
		light = NULL,
		z = 10,
		c = Inf,
		var_thr = var_thr_use
	)
	
	flag[is.na(flag)] <- 0
	df[[v]] <- ifelse(flag > 0, NA, df[[v]])
}

range_only_vars <- c("G_plate_1_1_1","G_plate_1_1_2",
										 "G_plate_1_1_3","G_plate_1_1_4",
										 "G_plate_1_1_5","G_plate_1_1_6",
										 "TS_CS65X_1_1_1","TS_CS65X_1_1_2","TS_CS65X_1_1_3",
										 "WS_05103_mean")

for (v in range_only_vars) {
	if (!v %in% names(df)) next
	
	var_thr_use <- calc_box_range(df[[v]], k = 5)
	
	flag <- despikeLF(
		as.data.frame(df %>% mutate(timestamp = TIMESTAMP)),
		v,
		iter = 1,
		light = NULL,
		z = 10,
		c = Inf,
		var_thr = var_thr_use
	)
	
	flag[is.na(flag)] <- 0
	df[[v]] <- ifelse(flag > 0, NA, df[[v]])
}


### Additional physical range
df <- df %>%
	mutate(PPFD_IN = ifelse(PPFD_IN > 2200, NA, PPFD_IN), 
				 across(any_of(c("RH_1_1_1","RH_1_1_2","RH_1_1_3")),
				 			 ~ ifelse(.x > 110 | .x < 0, NA, .x)),
				 across(any_of(c("SWC_1_1_1","SWC_1_1_2","SWC_1_1_3")),
				 			 ~ ifelse(.x > 100 | .x < 0, NA, .x))
				 ) 



df <- df %>%
  mutate(NETRAD = SW_IN - SW_OUT + LW_IN - LW_OUT,
         ALB = ifelse(SW_IN == 0 , 0, SW_OUT/SW_IN * 100)
  ) 

#################################################
#### QC 6: Multivariate comparison
#################################################
# To do: under development


################################################################################
#### potential radiation
################################################################################

calc_asce_radiation <- function(
    timestamp,
    lat,
    lon,
    elev,
    utc_offset,
    interval_mins = 30,
    Gsc = 4.92
) {
  
  # TIMESTAMP is stored as POSIXct with tz = "UTC",
  # but the displayed clock time represents local standard time.
  # Therefore, no timezone conversion should be applied.
  ts_mid <- timestamp
  
  # Extract the day of year and local standard clock time
  J <- yday(ts_mid)
  
  t_clock <- hour(ts_mid) +
    minute(ts_mid) / 60 +
    second(ts_mid) / 3600
  
  # Convert latitude to radians
  phi <- lat * pi / 180
  
  # Convert longitude to positive degrees west of Greenwich.
  # Input longitude is assumed to be negative west of Greenwich.
  Lm <- -lon
  
  # Longitude of the center of the local standard time zone,
  # expressed as positive degrees west of Greenwich
  Lz <- -15 * utc_offset
  
  # Inverse relative distance between Earth and the Sun
  # ASCE Eq. 50
  dr <- 1 + 0.033 * cos(
    2 * pi * J / 365
  )
  
  # Solar declination
  # ASCE Eq. 51
  delta <- 0.409 * sin(
    2 * pi * J / 365 - 1.39
  )
  
  # Seasonal correction for solar time
  # ASCE Eq. 57 and Eq. 58
  b <- 2 * pi * (J - 81) / 364
  
  Sc <- 0.1645 * sin(2 * b) -
    0.1255 * cos(b) -
    0.025 * sin(b)
  
  # Solar time angle at the midpoint of the interval
  # ASCE Eq. 55
  omega <- pi / 12 * (
    t_clock +
      0.06667 * (Lz - Lm) +
      Sc -
      12
  )
  
  # Half-width of the averaging interval in radians
  half_interval_angle <- pi * (interval_mins / 60) / 24
  
  # Solar time angles at the beginning and end of the interval
  # ASCE Eq. 53 and Eq. 54
  omega1 <- omega - half_interval_angle
  omega2 <- omega + half_interval_angle
  
  # Sunset hour angle
  # ASCE Eq. 59
  sunset_argument <- -tan(phi) * tan(delta)
  
  # Restrict the argument to the valid range of acos()
  sunset_argument <- pmax(
    -1,
    pmin(1, sunset_argument)
  )
  
  omega_s <- acos(sunset_argument)
  
  # Restrict the integration period to daylight hours
  # ASCE Eq. 56
  omega1 <- pmax(
    -omega_s,
    pmin(omega1, omega_s)
  )
  
  omega2 <- pmax(
    -omega_s,
    pmin(omega2, omega_s)
  )
  
  # Set the integration period to zero when the interval
  # falls completely outside daylight hours
  omega1 <- ifelse(
    omega1 > omega2,
    omega2,
    omega1
  )
  
  # Extraterrestrial radiation integrated over the interval
  # ASCE Eq. 48
  # Units: MJ m-2 interval-1
  Ra_MJ <- (12 / pi) * Gsc * dr * (
    (omega2 - omega1) * sin(phi) * sin(delta) +
      cos(phi) * cos(delta) *
      (sin(omega2) - sin(omega1))
  )
  
  Ra_MJ <- pmax(Ra_MJ, 0)
  
  # Convert interval-integrated extraterrestrial radiation
  # to average flux density in W m-2
  PotRad <- Ra_MJ * 1e6 / (interval_mins * 60)
  
  # Calculate clear-sky solar radiation
  # ASCE Eq. 47
  ClearSkyRad <- (
    0.75 + 2e-5 * elev
  ) * PotRad
  
  tibble(
    TIMESTAMP = timestamp,
    PotRad = PotRad,
    ClearSkyRad = ClearSkyRad
  )
}

lat  <- dirs_use$latitude
lon  <- dirs_use$longitude
elev <- dirs_use$elevation

rad_asce <- calc_asce_radiation(
  timestamp = df$TIMESTAMP,
  lat = lat,
  lon = lon,
  elev = elev,
  utc_offset = -8,
  interval_mins = 30
)

df <- df %>%
  left_join(rad_asce)

df_unit <- df_unit %>%
  mutate(
    PotRad = c("W m-2", "Avg"),
    ClearSkyRad = c("W m-2", "Avg")
  )


#################################################
#### Stat update and figures
#################################################
# gap stat update
df_stat$gap_after_QC <- as.numeric(colSums(is.na(df[, 2:(ncol(df)-2)])))
write_csv(df_stat,
						file = paste0(pipeline_path(base_dir, dirs_use$dir_output),"/",site_id,"_Level_1_gap_info.csv")
						)

# figures

#--------------------------------------------------
# function
#--------------------------------------------------
plot_qc_compare_pdf <- function(df_orig, df, vars,
																file_out = "qc_compare.pdf",
																chunk_size = 48 * 183,
																n_panel_per_page = 3,
																max_n_per_panel = 3000,
																clear_sky_comparison=FALSE) {
	
	stopifnot("TIMESTAMP" %in% names(df_orig), "TIMESTAMP" %in% names(df))
	
	vars_use <- vars[vars %in% names(df_orig) & vars %in% names(df)]
	
	if (length(vars_use) == 0) {
		stop("No variables in vars were found in both df_orig and df.")
	}
	
	# align by timestamp
	dat_all <- df_orig %>%
		select(TIMESTAMP, all_of(vars_use)) %>%
		rename_with(~ paste0(.x, "_orig"), -TIMESTAMP) %>%
		left_join(
			df %>%
				select(TIMESTAMP, all_of(vars_use)) %>%
				rename_with(~ paste0(.x, "_qc"), -TIMESTAMP),
			by = "TIMESTAMP"
		) %>%
		arrange(TIMESTAMP)
	
	if (clear_sky_comparison) {
	  dat_all <- dat_all %>%
	    left_join(df %>% select(TIMESTAMP,ClearSkyRad))
	}
	
	
	n <- nrow(dat_all)
	starts <- seq(1, n, by = chunk_size)
	
	# page unit plot list
	all_plots <- list()
	k <- 1
	
	for (v in vars_use) {
		v_orig <- paste0(v, "_orig")
		v_qc   <- paste0(v, "_qc")
		
		for (s in starts) {
			e <- min(s + chunk_size - 1, n)
			dat_sub <- dat_all[s:e, ]
			
			#--------------------------------------------
			# lighten: downsample
			#--------------------------------------------
			if (nrow(dat_sub) > max_n_per_panel) {
				step_n <- ceiling(nrow(dat_sub) / max_n_per_panel)
				dat_sub <- dat_sub[seq(1, nrow(dat_sub), by = step_n), ]
			}
			
			p <- ggplot(dat_sub, aes(x = TIMESTAMP)) +
				geom_line(aes(y = .data[[v_orig]]), color = "red", na.rm = TRUE, linewidth = 0.25) +
				geom_line(aes(y = .data[[v_qc]]),   color = "black", na.rm = TRUE, linewidth = 0.3)
			
			if (clear_sky_comparison) {
				p <- p +
					geom_line(aes( y = ClearSkyRad), color = "orange",na.rm = TRUE, linewidth = 0.25)
			}
			
			p <- p +
				labs(
					title = paste0(v, "  |  rows ", s, "-", e),
					subtitle = paste0(
						format(min(dat_sub$TIMESTAMP, na.rm = TRUE), "%Y-%m-%d %H:%M"),
						" to ",
						format(max(dat_sub$TIMESTAMP, na.rm = TRUE), "%Y-%m-%d %H:%M")
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
					plot.margin = margin(3, 6, 3, 3)
				)
			
			all_plots[[k]] <- p
			k <- k + 1
		}
	}
	
	# split plots into pages
	page_id <- ceiling(seq_along(all_plots) / n_panel_per_page)
	pages <- split(all_plots, page_id)
	
	pdf(file_out, width = 12, height = 14, onefile = TRUE, useDingbats = FALSE)
	
	for (i in seq_along(pages)) {
		page_plot <- plot_grid(
			plotlist = pages[[i]],
			ncol = 1,
			align = "v",
			axis = "lr"
		)
		
		grid::grid.draw(page_plot)
	}
	
	dev.off()
	
	message("Saved: ", file_out)
}

#--------------------------------------------------
# run
#--------------------------------------------------
key_vars <- c("SW_IN","SW_OUT","PPFD_IN",
              "LW_IN","LW_OUT","T_nr",
              "TA_1_1_3","RH_1_1_3","PA",
              "G_plate_1_1_1","G_plate_1_1_2",
              "G_plate_1_1_3","G_plate_1_1_4",
              "G_plate_1_1_5","G_plate_1_1_6",
              "TS_1_1_1","TS_1_1_2","TS_1_1_3",
              "TS_CS65X_1_1_1","TS_CS65X_1_1_2","TS_CS65X_1_1_3",
              "SWC_1_1_1","SWC_1_1_2","SWC_1_1_3", "WS_05103_mean")

plot_qc_compare_pdf(
	df_orig = pipeline_figure_data(df_orig),
	df = pipeline_figure_data(df),
	vars = key_vars,
	file_out = paste0(pipeline_path(base_dir, dirs_use$dir_output), "/figures/", site_id, "_Level_1_Figure_met_QAQC.pdf"),
	chunk_size = 48 * 365,
	n_panel_per_page = 3,
	max_n_per_panel = 2000
)



plot_qc_compare_pdf(
  df_orig = pipeline_figure_data(df_orig),
  df = pipeline_figure_data(df),
  vars = "SW_IN",
  file_out = paste0(pipeline_path(base_dir, dirs_use$dir_output), "/figures/", site_id, "_Level_1_Figure_met_QAQC_SW_IN_details.pdf"),
  chunk_size = 48 * 14,
  n_panel_per_page = 4,
  max_n_per_panel = 2000,
  clear_sky_comparison = T
)

#################################################
#### save
#################################################
head(df)
head(df_unit)

write_custom_csv <- pipeline_write

print(which(colnames(df) != colnames(df_unit)))

write_custom_csv(
	df = df,
	units_df = df_unit,
	file = paste0(pipeline_path(base_dir, dirs_use$dir_output),"/",site_id,"_Level_1_QAQCed_logger_data_",date(Sys.time()),".csv"),
	description = "Data from CSFormat (with easy flux output removed) have been merged with soil and bulk precipitation data. This represents QA/QCed met + soil data without EC output."
)


# loop all sites end
}
