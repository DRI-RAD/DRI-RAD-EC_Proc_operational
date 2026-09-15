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
library(ggpmisc)
library(gridExtra)
library(Metrics)

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
# The runner supplies the selected site.
# all sites
for(site_id in pipeline_options$sites){

dirs_use <- dirs %>% filter(site == site_id)

################################################################################
#### Load processed data 
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

met_unit1 <- read_csv(latest_file, skip = 1, n_max = 2)
met1 <- read_csv(latest_file, skip = 4, col_names = FALSE)
colnames(met1) <- colnames(met_unit1)
# Restrict input to the requested calculation window before processing.
met1 <- pipeline_slice(met1)

met1 <- met1 %>% select(!c(T_DP_1_1_1,T_DP_1_1_2,T_DP_1_1_3,T_nr,R_LW_in_meas,
													 R_LW_out_meas, sun_azimuth, sun_elevation,hour_angle,
													 sun_declination, air_mass_coeff, daytime, WS_05103_rslt,
													 WS_05103_5s_TMx, WS_05103_5s_Min,WS_05103_5s_TMn
													 ))
met_unit1 <- met_unit1 %>% select(!c(T_DP_1_1_1,T_DP_1_1_2,T_DP_1_1_3,T_nr,R_LW_in_meas,
													 R_LW_out_meas, sun_azimuth, sun_elevation,hour_angle,
													 sun_declination, air_mass_coeff, daytime, WS_05103_rslt,
													 WS_05103_5s_TMx, WS_05103_5s_Min,WS_05103_5s_TMn
))

if("P_bulk_raw" %in% colnames(met1)){
	met1 <- met1 %>% select(!c(P_bulk_raw))
	met_unit1 <- met_unit1 %>% select(!c(P_bulk_raw))
}


#--------------------------------------------------
# level 2
#--------------------------------------------------
files <- list.files(
	path = pipeline_path(base_dir, dirs_use$dir_output),
	pattern = "level_2_PI_vars_logger_data",
	full.names = TRUE,
	ignore.case = TRUE
)


latest_file <- pipeline_latest(files)

met_unit2 <- read_csv(latest_file, skip = 1, n_max = 2)
met2 <- read_csv(latest_file, skip = 4, col_names = FALSE)
colnames(met2) <- colnames(met_unit2)
# Restrict input to the requested calculation window before processing.
met2 <- pipeline_slice(met2)

df <- left_join(met1,met2, by = 'TIMESTAMP')
df_unit <- left_join(met_unit1, met_unit2, by = 'TIMESTAMP')



################################################################################
#### Load EC and QC flag
################################################################################

#--------------------------------------------------
# level 3 EC QC flag
#--------------------------------------------------
files <- list.files(
	path = pipeline_path(base_dir, dirs_use$dir_output),
	pattern = "Level_3_EC_with_qaqc",
	full.names = TRUE,
	ignore.case = TRUE
)


latest_file <- pipeline_latest(files)

ec_unit <- read_csv(latest_file, skip = 1, n_max = 2)
ec <- read_csv(latest_file, skip = 4, col_names = FALSE)
colnames(ec) <- colnames(ec_unit)
# Restrict input to the requested calculation window before processing.
ec <- pipeline_slice(ec)

EC_period <- c(ec$TIMESTAMP[min(which(!is.na(ec$LE)))] - 30*60, ec$TIMESTAMP[nrow(ec)] + 30*60)

df <- left_join(df,ec, by = 'TIMESTAMP')
df_unit <- left_join(df_unit, ec_unit, by = 'TIMESTAMP')

################################################################################
#### Load LI710 
################################################################################

if(!is.na(dirs_use$dir_LI710)){

	#--------------------------------------------------
	# level 3 LI710 QC flag
	#--------------------------------------------------
	files <- list.files(
		path = pipeline_path(base_dir, dirs_use$dir_output),
		pattern = "Level_3_li710_with_qaqc",
		full.names = TRUE,
		ignore.case = TRUE
	)
	
	
	latest_file <- pipeline_latest(files)
	
	li710_unit <- read_csv(latest_file, skip = 1, n_max = 2)
	li710 <- read_csv(latest_file, skip = 4, col_names = FALSE)
	colnames(li710) <- colnames(li710_unit)
	# Restrict input to the requested calculation window before processing.
	li710 <- pipeline_slice(li710)
	
	LI710_period <- c(li710$TIMESTAMP[min(which(!is.na(li710$LE_710)))] - 30*60, li710$TIMESTAMP[nrow(li710)] + 30*60)
	
	df <- left_join(df,li710, by = 'TIMESTAMP')
	df_unit <- left_join(df_unit, li710_unit, by = 'TIMESTAMP')
	
}

################################################################################
#### QAQC flag applications
################################################################################
df$LE_QC_despike[is.na(df$LE_QC_despike)] <- 0
df$H_QC_despike[is.na(df$H_QC_despike)] <- 0
df$FC_QC_despike[is.na(df$FC_QC_despike)] <- 0

df <- df %>%
	mutate(# statistical and technical filtering
		     LE_filtered1 = ifelse(LE_QC_despike >= 1 | LE_QC_sig_str >= 1 | LE_QC_fetch >= 1 | LE_QC_longrun >= 1,
																		NA, LE),
				 H_filtered1 = ifelse(H_QC_despike >= 1 | H_QC_sig_str >= 1 | H_QC_fetch >= 1 | H_QC_longrun >= 1,
				 													 NA, H),
				 FC_filtered1 = ifelse(FC_QC_despike >= 1 | FC_QC_sig_str >= 1 | FC_QC_fetch >= 1 | FC_QC_longrun >= 1,
				 													 NA, FC),
				 
				 # micrometeorological filtering (apply ustar threshold only NEE)
				 # LE_filtered2 = ifelse(LE_QC_SSITC >= 2 | (LE_QC_SSITC >= 1 & LE_QC_ustar >= 1),
				 # 													 NA, LE_filtered1),
				 # H_filtered2 = ifelse(H_QC_SSITC >= 2 | (H_QC_SSITC >= 1 & H_QC_ustar >= 1),
				 # 													NA, H_filtered1),
				 LE_filtered2 = ifelse(LE_QC_SSITC >= 2, NA, LE_filtered1),
				 H_filtered2 = ifelse(H_QC_SSITC >= 2,	NA, H_filtered1),
				 FC_filtered2 = ifelse(FC_QC_SSITC >= 2 | FC_QC_ustar >= 1, NA, FC_filtered1),
				 
				 # LI710 filter
				 LE_710_filtered = ifelse(LE_710_QC_diag >= 1 | LE_710_QC_despike >= 1 | LE_710_QC_longrun >= 1 | LE_710_QC_flow >= 1,
				 												 NA, LE_710),
				 H_710_filtered = ifelse(H_710_QC_diag >= 1 | H_710_QC_despike >= 1 | H_710_QC_longrun >= 1 | H_710_QC_flow >= 1,
				 												 NA, H_710)
				 )


################################################################################
#### Storage term apply
################################################################################
if(sum(!is.na(df$SLE)) > 0){
  flag <- despikeLF(as.data.frame(df %>% mutate(timestamp = TIMESTAMP)),
                    "SLE",iter = 10,light = "PotRad", z = 10)
  flag[is.na(flag)] <- 0
  df$SLE <- ifelse(flag >= 1, NA, df$SLE)
  flag <- despikeLF(as.data.frame(df %>% mutate(timestamp = TIMESTAMP)),
                    "SH",iter = 10,light = "PotRad", z = 10)
  flag[is.na(flag)] <- 0
  df$SH <- ifelse(flag >= 1, NA, df$SH)
  flag <- despikeLF(as.data.frame(df %>% mutate(timestamp = TIMESTAMP)),
                    "SC",iter = 10,light = "PotRad", z = 10)
  flag[is.na(flag)] <- 0
  df$SC <- ifelse(flag >= 1, NA, df$SC)
}

df <- df %>%
	mutate(
		LE_strg_added = ifelse(is.na(SLE), LE_filtered2, LE_filtered2 + SLE),
		H_strg_added = ifelse(is.na(SH), H_filtered2, H_filtered2 + SH),
		FC_strg_added = ifelse(is.na(SC), FC_filtered2, FC_filtered2 + SC),
		)



################################################################################
#### MDS gap filling available energy
################################################################################
EddyData <- df %>%
	mutate(Year = year(TIMESTAMP),
				 DoY = yday(TIMESTAMP),
				 Hour = hour(TIMESTAMP) + minute(TIMESTAMP)/60,
				 NEE = FC_strg_added,
				 LE = LE_strg_added,
				 H = H_strg_added, 
				 Rg = SW_IN_PI_F,
				 Tair = TA_PI_F,
				 rH = RH_PI_F,
				 VPD = VPD_PI_F,
				 Ustar = USTAR
	) %>%
	select(TIMESTAMP, Year,DoY,Hour,NEE,LE,H,Rg,Tair,rH,VPD,Ustar, NETRAD, G_PI) %>%
	as.data.frame()

#+++ Add time stamp in POSIX time format
EddyDataWithPosix <- EddyData %>% 
	fConvertTimeToPosix('YDH', Year = 'Year', Day = 'DoY', Hour = 'Hour')
#+++ Initalize R5 reference class sEddyProc for processing of eddy data
#+++ with all variables needed for processing later
EProc <- sEddyProc$new(
	site_id, EddyDataWithPosix, c('NEE','Rg','Tair','VPD', 'Ustar', 'LE', 'H','NETRAD',"G_PI"))

#++ Fill NEE gaps with MDS gap filling algorithm (without prior ustar filtering)
EProc$sMDSGapFill('NETRAD', FillAll = F)
EProc$sMDSGapFill('G_PI', FillAll = F)

MDSout <- data.frame(TIMESTAMP = EddyData$TIMESTAMP, 
										 EProc$sExportResults() %>% 
										 	select(NETRAD_f, G_PI_f) %>%
										 	rename(NETRAD_PI_F = NETRAD_f, G_PI_F = G_PI_f)
)

df <- left_join(df, MDSout)

################################################################################
#### MDS gap filling EC
################################################################################
EddyData <- df %>%
	filter(TIMESTAMP >= EC_period[1] & TIMESTAMP <= EC_period[2]) %>%
	mutate(Year = year(TIMESTAMP),
				 DoY = yday(TIMESTAMP),
				 Hour = hour(TIMESTAMP) + minute(TIMESTAMP)/60,
				 NEE = FC_strg_added,
				 LE = LE_strg_added,
				 H = H_strg_added, 
				 Rg = SW_IN_PI_F,
				 Tair = TA_PI_F,
				 rH = RH_PI_F,
				 VPD = VPD_PI_F,
				 Ustar = USTAR
	) %>%
	select(TIMESTAMP, Year,DoY,Hour,NEE,LE,H,Rg,Tair,rH,VPD,Ustar, NETRAD, G_PI) %>%
	as.data.frame()

#+++ Add time stamp in POSIX time format
EddyDataWithPosix <- EddyData %>% 
	fConvertTimeToPosix('YDH', Year = 'Year', Day = 'DoY', Hour = 'Hour')
#+++ Initalize R5 reference class sEddyProc for processing of eddy data
#+++ with all variables needed for processing later
EProc <- sEddyProc$new(
	site_id, EddyDataWithPosix, c('NEE','Rg','Tair','VPD', 'Ustar', 'LE', 'H','NETRAD',"G_PI"))

#++ Fill NEE gaps with MDS gap filling algorithm (without prior ustar filtering)
EProc$sMDSGapFill('LE', FillAll = T)
EProc$sMDSGapFill('H', FillAll = F)
EProc$sMDSGapFill('NEE', FillAll = F)

MDSout <- data.frame(TIMESTAMP = EddyData$TIMESTAMP, 
										 EProc$sExportResults() %>% 
										   mutate(LE_fall = ifelse(LE_fall_qc == 1, LE_fall, NA)) %>%
										 	select(LE_f, LE_fall,H_f, NEE_f) %>%
										 	rename( LE_PI_F = LE_f, LE_MDS= LE_fall, H_PI_F = H_f, FC_PI_F = NEE_f)
										 )

df <- left_join(df, MDSout)


################################################################################
#### MDS gap filling LI710
################################################################################
EddyData <- df %>%
	filter(TIMESTAMP >= LI710_period[1] & TIMESTAMP <= LI710_period[2]) %>%
	mutate(Year = year(TIMESTAMP),
				 DoY = yday(TIMESTAMP),
				 Hour = hour(TIMESTAMP) + minute(TIMESTAMP)/60,
				 LE = LE_710_filtered,
				 H = H_710_filtered, 
				 Rg = SW_IN_PI_F,
				 Tair = TA_PI_F,
				 rH = RH_PI_F,
				 VPD = VPD_PI_F,
				 Ustar = USTAR
	) %>%
	select(TIMESTAMP, Year,DoY,Hour,LE,H,Rg,Tair,rH,VPD,Ustar) %>%
	as.data.frame()

#+++ Add time stamp in POSIX time format
EddyDataWithPosix <- EddyData %>% 
	fConvertTimeToPosix('YDH', Year = 'Year', Day = 'DoY', Hour = 'Hour')
#+++ Initalize R5 reference class sEddyProc for processing of eddy data
#+++ with all variables needed for processing later
EProc <- sEddyProc$new(
	site_id, EddyDataWithPosix, c('Rg','Tair','VPD', 'Ustar', 'LE', 'H'))

#++ Fill NEE gaps with MDS gap filling algorithm (without prior ustar filtering)
EProc$sMDSGapFill('LE', FillAll = T)
EProc$sMDSGapFill('H', FillAll = F)
MDSout <- data.frame(TIMESTAMP = EddyData$TIMESTAMP, 
										 EProc$sExportResults() %>% 
										   mutate(LE_fall = ifelse(LE_fall_qc == 1, LE_fall, NA)) %>%
										 	select(LE_f, LE_fall,H_f) %>%
										 	rename(LE_710_PI_F = LE_f, LE_710_MDS= LE_fall, H_710_PI_F = H_f)
)

df <- left_join(df, MDSout)



################################################################################
#### EBR correction (modified Mauder et al., 2013) 
################################################################################
# preserving daytime Bowen ratio, while correcting turbulent flux only in daytime
# EBR_d*daytime_sum(LE+H) + nighttime_sum(LE+H) = daytime_sum(Rn - G) + nighttime_sum(Rn - G)
# therefore, EBR_d = (daytime_sum(Rn - G) + nighttime_sum(Rn - G) - nighttime_sum(LE+H))/daytime_sum(LE+H)

use <- df %>%
  mutate(DATE = date(TIMESTAMP),
         is_day = PotRad > 50
  ) %>%
  filter(!is.na(NETRAD_PI_F) & !is.na(G_PI_F) &!is.na(LE_PI_F) &!is.na(H_PI_F) & !is.na(is_day)) %>%
  group_by(DATE) %>%
  summarise(
    day_AE = sum((NETRAD_PI_F - G_PI_F)[is_day]),
    night_AE = sum((NETRAD_PI_F - G_PI_F)[!is_day]),
    day_turb = sum((LE_PI_F + H_PI_F)[is_day]),
    night_turb = sum((LE_PI_F + H_PI_F)[!is_day]),
    S_n = night_AE - night_turb,
    S_d = -S_n,
    EBR_d = day_turb/(day_AE - S_d),
    EBR_d = ifelse(EBR_d < 0, NA, ifelse(EBR_d < 0.3, 0.3, ifelse(EBR_d > 1.7, 1.7, EBR_d)))
  ) %>%
  select(DATE,EBR_d)

df <- df %>%
  mutate(DATE = date(TIMESTAMP),
         is_day = PotRad > 50) %>%
  left_join(use) %>%
  mutate(
    EBR_d = ifelse(is_day,EBR_d,NA),
    LE_PI_CORR = ifelse(!is.na(EBR_d), LE_PI_F / EBR_d, LE_PI_F),
    H_PI_CORR  = ifelse(!is.na(EBR_d), H_PI_F  / EBR_d, H_PI_F)
  ) %>% select(!c('DATE','is_day'))

use <- df %>%
  mutate(DATE = date(TIMESTAMP),
         is_day = PotRad > 50
  ) %>%
  filter(!is.na(NETRAD_PI_F) & !is.na(G_PI_F) &!is.na(LE_710_PI_F) &!is.na(H_710_PI_F) & !is.na(is_day)) %>%
  group_by(DATE) %>%
  summarise(
    day_AE = sum((NETRAD_PI_F - G_PI_F)[is_day], na.rm = TRUE),
    night_AE = sum((NETRAD_PI_F - G_PI_F)[!is_day], na.rm = TRUE),
    day_turb = sum((LE_710_PI_F + H_710_PI_F)[is_day], na.rm = TRUE),
    night_turb = sum((LE_710_PI_F + H_710_PI_F)[!is_day], na.rm = TRUE),
    EBR_d_710 = day_turb/(day_AE + night_AE - night_turb),
    EBR_d_710 = ifelse(EBR_d_710 < 0, NA, ifelse(EBR_d_710 < 0.3, 0.3, ifelse(EBR_d_710 > 1.7, 1.7, EBR_d_710)))
  ) %>%
  select(DATE,EBR_d_710)

df <- df %>%
  mutate(DATE = date(TIMESTAMP),
         is_day = PotRad > 50) %>%
  left_join(use) %>%
  mutate(
    EBR_d_710 = ifelse(is_day,EBR_d_710,NA),
    LE_710_PI_CORR = ifelse(!is.na(EBR_d_710),LE_710_PI_F / EBR_d_710, LE_710_PI_F),
    H_710_PI_CORR = ifelse(!is.na(EBR_d_710), H_710_PI_F / EBR_d_710, H_710_PI_F),
  ) %>% select(!c('DATE','is_day'))


################################################################################
#### energy balance approach: VBR and residual energy
################################################################################
# Variance Bowen ratio 
# T_SONIC_SIGMA  : sonic temperature sd [K]
# H2O_SIGMA : water vapor sd [mmolH2O mol-1]


df <- df %>%
	mutate(lv = latent.heat.vaporization(TA_PI_F),
				 q_SIGMA =  0.622 * (H2O_SIGMA / 1000) / (1 - (1 - 0.622) * H2O / 1000),
				 VBR = 1004.834/lv * (T_SONIC_SIGMA/q_SIGMA - 0.51 * (TA_PI_F + 273.15)),
				 VBR_Ts = 1004.834/lv * (T_SONIC_SIGMA/q_SIGMA),
				 LE_VBR = (NETRAD - G_PI)/(1+VBR),
				 LE_VBR_Ts = (NETRAD - G_PI)/(1+VBR_Ts),
				 
				 #residual energy LE
				 # H_VBR = un_H * H_scf / (1+0.51*1004.834*(TA_PI_F + 273.15)/(lv*VBR)),
				 LE_RE = NETRAD - G_PI - H_strg_added
				 )


### LE_VBR LE_RE gap filling
EddyData <- df %>%
	filter(TIMESTAMP >= EC_period[1] & TIMESTAMP <= EC_period[2]) %>%
	mutate(Year = year(TIMESTAMP),
				 DoY = yday(TIMESTAMP),
				 Hour = hour(TIMESTAMP) + minute(TIMESTAMP)/60,
				 LE_VBR = LE_VBR,
				 LE_VBR_Ts = LE_VBR_Ts,
				 LE_RE= LE_RE,
				 Rg = SW_IN_PI_F,
				 Tair = TA_PI_F,
				 rH = RH_PI_F,
				 VPD = VPD_PI_F,
				 Ustar = USTAR
	) %>%
	select(TIMESTAMP, Year,DoY,Hour,LE_VBR,LE_VBR_Ts,LE_RE,Rg,Tair,rH,VPD,Ustar) %>%
	as.data.frame()

#+++ Add time stamp in POSIX time format
EddyDataWithPosix <- EddyData %>% 
	fConvertTimeToPosix('YDH', Year = 'Year', Day = 'DoY', Hour = 'Hour')
#+++ Initalize R5 reference class sEddyProc for processing of eddy data
#+++ with all variables needed for processing later
EProc <- sEddyProc$new(
	site_id, EddyDataWithPosix, c('Rg','Tair','VPD', 'Ustar', 'LE_RE','LE_VBR','LE_VBR_Ts'))

#++ Fill NEE gaps with MDS gap filling algorithm
EProc$sMDSGapFill('LE_VBR', FillAll = T)
EProc$sMDSGapFill('LE_VBR_Ts', FillAll = T)
EProc$sMDSGapFill('LE_RE', FillAll = T)
MDSout <- data.frame(TIMESTAMP = EddyData$TIMESTAMP, 
										 EProc$sExportResults() %>% 
										   mutate(LE_RE_fall = ifelse(LE_RE_fall_qc == 1, LE_RE_fall, NA),
										          LE_VBR_fall = ifelse(LE_VBR_fall_qc == 1, LE_VBR_fall, NA),
										          LE_VBR_fall = ifelse(LE_VBR_Ts_fall_qc == 1, LE_VBR_Ts_fall, NA),
										          ) %>%
										 	select(LE_RE_f,LE_RE_fall, LE_VBR_f, LE_VBR_fall, LE_VBR_Ts_f, LE_VBR_Ts_fall) %>%
										 	rename(LE_RE_F = LE_RE_f,LE_RE_MDS = LE_RE_fall, 
										 	       LE_VBR_F = LE_VBR_f,  LE_VBR_MDS = LE_VBR_fall,
										 	       LE_VBR_Ts_F = LE_VBR_Ts_f,  LE_VBR_Ts_MDS = LE_VBR_Ts_fall
										 	       )
)

df <- left_join(df, MDSout)


################################################################################
#### ET calculation (mm/ 30min)
################################################################################
df <- df %>%
	mutate(ET_EC = LE_strg_added/lv * 60*30,
				 ET_EC_F = LE_PI_F/lv * 60*30,
				 ET_EC_EBC = LE_PI_CORR/lv * 60*30,
				 ET_710 = LE_710_filtered/lv * 60*30,
				 ET_710_F = LE_710_PI_F/lv * 60*30,
				 ET_710_EBC = LE_710_PI_CORR/lv * 60*30,
				 ET_VBR = LE_VBR/lv * 60*30,
				 ET_VBR_F = LE_VBR_F/lv * 60*30,
				 ET_VBR_Ts = LE_VBR_Ts/lv * 60*30,
				 ET_VBR_Ts_F = LE_VBR_Ts_F/lv * 60*30,
				 ET_RE = LE_RE/lv * 60*30,
				 ET_RE_F = LE_RE_F/lv * 60*30,
				 ) %>%
	select(!c(lv,q_SIGMA))


# ################################################################################
# #### Diagnostic figure 1a: Energy balance 
# ################################################################################
# # ec
# ebr_val <- df %>%
# 	mutate(
# 		RnG = NETRAD - G_PI,
# 		LEH = LE_strg_added + H_strg_added
# 	) %>%
# 	filter(!is.na(RnG) & !is.na(LEH)) %>%
# 	summarise(EBR_EC = sum(LEH) / sum(RnG)) %>%
# 	pull(EBR_EC) %>% round(3)
# 
# a <- df %>%
# 	ggplot(aes(NETRAD - G_PI,LE_strg_added + H_strg_added)) +
# 	theme_bw() +
# 	geom_point(alpha = 0.2) +
# 	geom_abline(lty = 2) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "LE + H (W m-2)",
# 			 x = "Rn - G (W m-2)",
# 			 title = paste0(site_id,": Eddy covariance (half hourly)")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
# 	annotate("text",x = -Inf, y = Inf,label = paste0("EBR = ", ebr_val),	hjust = -0.2,vjust = 4, size = 3.5)
# 
# b <- df %>%
# 	mutate(DATE = date(TIMESTAMP)) %>%
# 	group_by(DATE) %>%
# 	summarise_all(mean) %>%
# 	ggplot(aes(NETRAD_PI_F - G_PI_F,LE_PI_F + H_PI_F)) +
# 	theme_bw() +
# 	geom_point() +
# 	geom_abline(lty = 2) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "LE + H (W m-2)",
# 			 x = "Rn - G (W m-2)",
# 			 title = paste0(site_id,": Eddy covariance (daily)")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0)
# 
# c <- df %>%
# 	mutate(
# 		HOUR = hour(TIMESTAMP) + minute(TIMESTAMP) / 60,
# 		`Rn - G` = NETRAD - G_PI,
# 		`LE + H` = LE_strg_added + H_strg_added
# 	) %>%
# 	filter(!is.na(`Rn - G`) & !is.na(`LE + H`)) %>%
# 	group_by(HOUR) %>%
# 	summarise(
# 		`Rn - G` = mean(`Rn - G`, na.rm = TRUE),
# 		`LE + H` = mean(`LE + H`, na.rm = TRUE),
# 		LE = mean(LE_strg_added, na.rm = TRUE),
# 		H = mean(H_strg_added, na.rm = TRUE),
# 		.groups = "drop"
# 	) %>%
# 	pivot_longer(
# 		cols = c(`Rn - G`, `LE + H`,LE,H),
# 		names_to = "var",
# 		values_to = "value"
# 	) %>%
# 	ggplot(aes(HOUR, value, color = var, shape = var)) +
# 	theme_bw() +
# 	geom_hline(yintercept = 0) +
# 	geom_line() +
# 	geom_point(size = 2) +
# 	scale_color_manual(
# 		name = NULL,
# 		breaks = c("Rn - G", "LE + H","LE", "H"),
# 		values = c(
# 			"Rn - G" = "#000000",        # black
# 			"LE + H" = "#999999",        # gray
# 			"LE" = "#0072B2",        # blue
# 			"H" = "#D55E00" # vermillion
# 		)
# 	) +
# 	scale_shape_manual(
# 		name = NULL,
# 		breaks = c("Rn - G", "LE + H","LE", "H"),
# 		values = c(
# 			"Rn - G" = 1,   # open circle
# 			"LE + H" = 16,  # filled circle
# 			"LE" = 18, 
# 			"H" = 17  # triangle
# 		)
# 	) +
# 	labs(y = "Energy flux (W m-2)",
# 			 x = "",
# 			 title = paste0(site_id, ": EC mean diurnal cycle")
# 	) +
# 	theme(legend.position = "bottom",
# 				axis.title.x = element_blank())
# 
# 
# d <- df %>%
# 	filter(TIMESTAMP >= EC_period[1] & TIMESTAMP <= EC_period[2]) %>%
# 	mutate(DATE = date(TIMESTAMP)) %>%
# 	group_by(DATE) %>%
# 	summarise_all(mean,na.rm = F) %>% 
# 	ggplot(aes(DATE)) +
# 	theme_bw() +
# 	geom_hline(yintercept = 0, lty = 2) +
# 	geom_line(aes(y = NETRAD_PI_F, color = "Rn")) +
# 	geom_line(aes(y = G_PI_F, color = "G")) +
# 	geom_line(aes(y = LE_PI_F, color = "LE")) +
# 	geom_line(aes(y = H_PI_F, color = "H")) +
# 	scale_color_manual(
# 		name = NULL,
# 		breaks = c("Rn", "G","LE", "H"),
# 		values = c(
# 			"Rn" = "#000000",        # black
# 			"G" = "#999999",        # gray
# 			"LE" = "#0072B2",        # blue
# 			"H" = "#D55E00" # vermillion
# 		)
# 	) +
# 	labs(y = "Energy flux (W m-2)",
# 			 x = "",
# 			 title = paste0(site_id, ": EC daily mean energy flux (EBR uncorrected)")
# 	)
# 
# p1 <- plot_grid(plot_grid(a,b,c,ncol = 3),d,ncol = 1)
# 
# 
# # 710
# 
# ebr_val <- df %>%
# 	mutate(
# 		RnG = NETRAD - G_PI,
# 		LEH = LE_710_filtered + H_710_filtered
# 	) %>%
# 	filter(!is.na(RnG) & !is.na(LEH)) %>%
# 	summarise(EBR_EC = sum(LEH) / sum(RnG)) %>%
# 	pull(EBR_EC) %>% round(3)
# 
# a <- df %>%
# 	ggplot(aes(NETRAD - G_PI,LE_710_filtered + H_710_filtered)) +
# 	theme_bw() +
# 	geom_point(alpha = 0.2) +
# 	geom_abline(lty = 2) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "LE + H (W m-2)",
# 			 x = "Rn - G (W m-2)",
# 			 title = paste0(site_id,": LI710 (half hourly)")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
# 	annotate("text",x = -Inf, y = Inf,label = paste0("EBR = ", ebr_val),	hjust = -0.2,vjust = 4, size = 3.5)
# 
# 
# b <- df %>%
# 	mutate(DATE = date(TIMESTAMP)) %>%
# 	group_by(DATE) %>%
# 	summarise_all(mean) %>%
# 	ggplot(aes(NETRAD_PI_F - G_PI_F,LE_710_PI_F + H_710_PI_F)) +
# 	theme_bw() +
# 	geom_point() +
# 	geom_abline(lty = 2) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "LE + H (W m-2)",
# 			 x = "Rn - G (W m-2)",
# 			 title = paste0(site_id,": LI710 (daily)")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0)
# 
# 
# c <- df %>%
# 	mutate(
# 		HOUR = hour(TIMESTAMP) + minute(TIMESTAMP) / 60,
# 		`Rn - G` = NETRAD - G_PI,
# 		`LE + H` = LE_710_filtered + H_710_filtered
# 	) %>%
# 	filter(!is.na(`Rn - G`) & !is.na(`LE + H`)) %>%
# 	group_by(HOUR) %>%
# 	summarise(
# 		`Rn - G` = mean(`Rn - G`, na.rm = TRUE),
# 		`LE + H` = mean(`LE + H`, na.rm = TRUE),
# 		LE = mean(LE_710_filtered, na.rm = TRUE),
# 		H = mean(H_710_filtered, na.rm = TRUE),
# 		.groups = "drop"
# 	) %>%
# 	pivot_longer(
# 		cols = c(`Rn - G`, `LE + H`,LE,H),
# 		names_to = "var",
# 		values_to = "value"
# 	) %>%
# 	ggplot(aes(HOUR, value, color = var, shape = var)) +
# 	theme_bw() +
# 	geom_hline(yintercept = 0) +
# 	geom_line() +
# 	geom_point(size = 2) +
# 	scale_color_manual(
# 		name = NULL,
# 		breaks = c("Rn - G", "LE + H","LE", "H"),
# 		values = c(
# 			"Rn - G" = "#000000",        # black
# 			"LE + H" = "#999999",        # gray
# 			"LE" = "#0072B2",        # blue
# 			"H" = "#D55E00" # vermillion
# 		)
# 	) +
# 	scale_shape_manual(
# 		name = NULL,
# 		breaks = c("Rn - G", "LE + H","LE", "H"),
# 		values = c(
# 			"Rn - G" = 1,   # open circle
# 			"LE + H" = 16,  # filled circle
# 			"LE" = 18, 
# 			"H" = 17  # triangle
# 		)
# 	) +
# 	labs(y = "Energy flux (W m-2)",
# 			 x = "",
# 			 title = paste0(site_id, ": 710 mean diurnal cycle")
# 	) +
# 	theme(legend.position = "bottom",
# 				axis.title.x = element_blank())
# 
# 
# d <- df %>%
# 	filter(TIMESTAMP >= EC_period[1] & TIMESTAMP <= EC_period[2]) %>%
# 	mutate(DATE = date(TIMESTAMP)) %>%
# 	group_by(DATE) %>%
# 	summarise_all(mean,na.rm = F) %>% 
# 	ggplot(aes(DATE)) +
# 	theme_bw() +
# 	geom_hline(yintercept = 0, lty = 2) +
# 	geom_line(aes(y = NETRAD_PI_F, color = "Rn")) +
# 	geom_line(aes(y = G_PI_F, color = "G")) +
# 	geom_line(aes(y = LE_710_PI_F, color = "LE")) +
# 	geom_line(aes(y = H_710_PI_F, color = "H")) +
# 	scale_color_manual(
# 		name = NULL,
# 		breaks = c("Rn", "G","LE", "H"),
# 		values = c(
# 			"Rn" = "#000000",        # black
# 			"G" = "#999999",        # gray
# 			"LE" = "#0072B2",        # blue
# 			"H" = "#D55E00" # vermillion
# 		)
# 	) +
# 	labs(y = "Energy flux (W m-2)",
# 			 x = "",
# 			 title = paste0(site_id, ": 710 daily mean energy flux (EBR uncorrected)")
# 	)
# 
# p2 <- plot_grid(plot_grid(a,b,c,ncol = 3),d,ncol = 1)
# 
# ggsave(plot = plot_grid(p1,p2,ncol = 1),
# 			 filename = paste0(pipeline_path(base_dir, dirs_use$dir_output),"/figures/",site_id,"_Level_4_Figure_1a_energy_balance.jpeg"),
# 			 width = 11, height = 13, dpi = 150
# )
# 
# 
# 
# ################################################################################
# #### Diagnostic figure 1b: Energy balance (corrected)
# ################################################################################
# # ec
# ebr_val <- df %>%
# 	mutate(
# 		RnG = NETRAD - G_PI,
# 		LEH = LE_PI_CORR + H_PI_CORR
# 	) %>%
# 	filter(!is.na(RnG) & !is.na(LEH)) %>%
# 	summarise(EBR_EC = sum(LEH) / sum(RnG)) %>%
# 	pull(EBR_EC) %>% round(3)
# 
# a <- df %>%
# 	ggplot(aes(NETRAD - G_PI,LE_PI_CORR + H_PI_CORR)) +
# 	theme_bw() +
# 	geom_point(alpha = 0.2) +
# 	geom_abline(lty = 2) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "LE + H (W m-2)",
# 			 x = "Rn - G (W m-2)",
# 			 title = paste0(site_id,": corrected EC (half hourly)")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
# 	annotate("text",x = -Inf, y = Inf,label = paste0("EBR = ", ebr_val),	hjust = -0.2,vjust = 4, size = 3.5)
# 
# b <- df %>%
# 	mutate(DATE = date(TIMESTAMP)) %>%
# 	group_by(DATE) %>%
# 	summarise_all(mean) %>%
# 	ggplot(aes(NETRAD_PI_F - G_PI_F,LE_PI_CORR + H_PI_CORR)) +
# 	theme_bw() +
# 	geom_point() +
# 	geom_abline(lty = 2) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "LE + H (W m-2)",
# 			 x = "Rn - G (W m-2)",
# 			 title = paste0(site_id,": corrected EC (daily)")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0)
# 
# c <- df %>%
# 	mutate(
# 		HOUR = hour(TIMESTAMP) + minute(TIMESTAMP) / 60,
# 		`Rn - G` = NETRAD - G_PI,
# 		`LE + H` = LE_PI_CORR + H_PI_CORR
# 	) %>%
# 	filter(!is.na(`Rn - G`) & !is.na(`LE + H`)) %>%
# 	group_by(HOUR) %>%
# 	summarise(
# 		`Rn - G` = mean(`Rn - G`, na.rm = TRUE),
# 		`LE + H` = mean(`LE + H`, na.rm = TRUE),
# 		LE = mean(LE_PI_CORR, na.rm = TRUE),
# 		H = mean(H_PI_CORR, na.rm = TRUE),
# 		.groups = "drop"
# 	) %>%
# 	pivot_longer(
# 		cols = c(`Rn - G`, `LE + H`,LE,H),
# 		names_to = "var",
# 		values_to = "value"
# 	) %>%
# 	ggplot(aes(HOUR, value, color = var, shape = var)) +
# 	theme_bw() +
# 	geom_hline(yintercept = 0) +
# 	geom_line() +
# 	geom_point(size = 2) +
# 	scale_color_manual(
# 		name = NULL,
# 		breaks = c("Rn - G", "LE + H","LE", "H"),
# 		values = c(
# 			"Rn - G" = "#000000",        # black
# 			"LE + H" = "#999999",        # gray
# 			"LE" = "#0072B2",        # blue
# 			"H" = "#D55E00" # vermillion
# 		)
# 	) +
# 	scale_shape_manual(
# 		name = NULL,
# 		breaks = c("Rn - G", "LE + H","LE", "H"),
# 		values = c(
# 			"Rn - G" = 1,   # open circle
# 			"LE + H" = 16,  # filled circle
# 			"LE" = 18, 
# 			"H" = 17  # triangle
# 		)
# 	) +
# 	labs(y = "Energy flux (W m-2)",
# 			 x = "",
# 			 title = paste0(site_id, ": corrected EC mdc")
# 	) +
# 	theme(legend.position = "bottom",
# 				axis.title.x = element_blank())
# 
# 
# d <- df %>%
# 	filter(TIMESTAMP >= EC_period[1] & TIMESTAMP <= EC_period[2]) %>%
# 	mutate(DATE = date(TIMESTAMP)) %>%
# 	group_by(DATE) %>%
# 	summarise_all(mean,na.rm = F) %>% 
# 	ggplot(aes(DATE)) +
# 	theme_bw() +
# 	geom_hline(yintercept = 0, lty = 2) +
# 	geom_line(aes(y = NETRAD_PI_F, color = "Rn")) +
# 	geom_line(aes(y = G_PI_F, color = "G")) +
# 	geom_line(aes(y = LE_PI_CORR, color = "LE")) +
# 	geom_line(aes(y = H_PI_CORR, color = "H")) +
# 	scale_color_manual(
# 		name = NULL,
# 		breaks = c("Rn", "G","LE", "H"),
# 		values = c(
# 			"Rn" = "#000000",        # black
# 			"G" = "#999999",        # gray
# 			"LE" = "#0072B2",        # blue
# 			"H" = "#D55E00" # vermillion
# 		)
# 	) +
# 	labs(y = "Energy flux (W m-2)",
# 			 x = "",
# 			 title = paste0(site_id, ": EC daily mean energy flux (EBR corrected)")
# 	)
# 
# p1 <- plot_grid(plot_grid(a,b,c,ncol = 3),d,ncol = 1)
# 
# 
# # 710
# 
# ebr_val <- df %>%
# 	mutate(
# 		RnG = NETRAD - G_PI,
# 		LEH = LE_710_PI_CORR + H_710_PI_CORR
# 	) %>%
# 	filter(!is.na(RnG) & !is.na(LEH)) %>%
# 	summarise(EBR_EC = sum(LEH) / sum(RnG)) %>%
# 	pull(EBR_EC) %>% round(3)
# 
# a <- df %>%
# 	ggplot(aes(NETRAD - G_PI,LE_710_PI_CORR + H_710_PI_CORR)) +
# 	theme_bw() +
# 	geom_point(alpha = 0.2) +
# 	geom_abline(lty = 2) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "LE + H (W m-2)",
# 			 x = "Rn - G (W m-2)",
# 			 title = paste0(site_id,": corrected LI710 (half hourly)")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
# 	annotate("text",x = -Inf, y = Inf,label = paste0("EBR = ", ebr_val),	hjust = -0.2,vjust = 4, size = 3.5)
# 
# 
# b <- df %>%
# 	mutate(DATE = date(TIMESTAMP)) %>%
# 	group_by(DATE) %>%
# 	summarise_all(mean) %>%
# 	ggplot(aes(NETRAD_PI_F - G_PI_F,LE_710_PI_CORR + H_710_PI_CORR)) +
# 	theme_bw() +
# 	geom_point() +
# 	geom_abline(lty = 2) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "LE + H (W m-2)",
# 			 x = "Rn - G (W m-2)",
# 			 title = paste0(site_id,": corrected LI710 (daily)")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0)
# 
# 
# c <- df %>%
# 	mutate(
# 		HOUR = hour(TIMESTAMP) + minute(TIMESTAMP) / 60,
# 		`Rn - G` = NETRAD - G_PI,
# 		`LE + H` = LE_710_PI_CORR + H_710_PI_CORR
# 	) %>%
# 	filter(!is.na(`Rn - G`) & !is.na(`LE + H`)) %>%
# 	group_by(HOUR) %>%
# 	summarise(
# 		`Rn - G` = mean(`Rn - G`, na.rm = TRUE),
# 		`LE + H` = mean(`LE + H`, na.rm = TRUE),
# 		LE = mean(LE_710_PI_CORR, na.rm = TRUE),
# 		H = mean(H_710_PI_CORR, na.rm = TRUE),
# 		.groups = "drop"
# 	) %>%
# 	pivot_longer(
# 		cols = c(`Rn - G`, `LE + H`,LE,H),
# 		names_to = "var",
# 		values_to = "value"
# 	) %>%
# 	ggplot(aes(HOUR, value, color = var, shape = var)) +
# 	theme_bw() +
# 	geom_hline(yintercept = 0) +
# 	geom_line() +
# 	geom_point(size = 2) +
# 	scale_color_manual(
# 		name = NULL,
# 		breaks = c("Rn - G", "LE + H","LE", "H"),
# 		values = c(
# 			"Rn - G" = "#000000",        # black
# 			"LE + H" = "#999999",        # gray
# 			"LE" = "#0072B2",        # blue
# 			"H" = "#D55E00" # vermillion
# 		)
# 	) +
# 	scale_shape_manual(
# 		name = NULL,
# 		breaks = c("Rn - G", "LE + H","LE", "H"),
# 		values = c(
# 			"Rn - G" = 1,   # open circle
# 			"LE + H" = 16,  # filled circle
# 			"LE" = 18, 
# 			"H" = 17  # triangle
# 		)
# 	) +
# 	labs(y = "Energy flux (W m-2)",
# 			 x = "",
# 			 title = paste0(site_id, ": corrected 710 mdc")
# 	) +
# 	theme(legend.position = "bottom",
# 				axis.title.x = element_blank())
# 
# 
# d <- df %>%
# 	filter(TIMESTAMP >= EC_period[1] & TIMESTAMP <= EC_period[2]) %>%
# 	mutate(DATE = date(TIMESTAMP)) %>%
# 	group_by(DATE) %>%
# 	summarise_all(mean,na.rm = F) %>% 
# 	ggplot(aes(DATE)) +
# 	theme_bw() +
# 	geom_hline(yintercept = 0, lty = 2) +
# 	geom_line(aes(y = NETRAD_PI_F, color = "Rn")) +
# 	geom_line(aes(y = G_PI_F, color = "G")) +
# 	geom_line(aes(y = LE_710_PI_CORR, color = "LE")) +
# 	geom_line(aes(y = H_710_PI_CORR, color = "H")) +
# 	scale_color_manual(
# 		name = NULL,
# 		breaks = c("Rn", "G","LE", "H"),
# 		values = c(
# 			"Rn" = "#000000",        # black
# 			"G" = "#999999",        # gray
# 			"LE" = "#0072B2",        # blue
# 			"H" = "#D55E00" # vermillion
# 		)
# 	) +
# 	labs(y = "Energy flux (W m-2)",
# 			 x = "",
# 			 title = paste0(site_id, ": 710 daily mean energy flux (EBR corrected)")
# 	)
# 
# p2 <- plot_grid(plot_grid(a,b,c,ncol = 3),d,ncol = 1)
# 
# ggsave(plot = plot_grid(p1,p2,ncol = 1),
# 			 filename = paste0(pipeline_path(base_dir, dirs_use$dir_output),"/figures/",site_id,"_Level_4_Figure_1b_energy_balance_corrected.jpeg"),
# 			 width = 11, height = 13, dpi = 150
# )
# 
# ################################################################################
# #### Diagnostic figure 2: water balance 
# ################################################################################
# offset <- 5   
# scale_factor <- 0.1
# if(site_id == "ERVP"){EC_period[1] <- as.POSIXct("2025-01-10 10:00:00 UTC")}
# 
# a <- df %>%
# 	filter(TIMESTAMP >= EC_period[1] & TIMESTAMP <= EC_period[2]) %>%
# 	mutate(date = date(TIMESTAMP)) %>%
# 	group_by(date) %>%
# 	summarise_if(is.numeric,sum) %>%
# 	rowwise() %>%
# 	mutate(ET_max = max(ET_EC_F,ET_EC_EBC,na.rm=T),
# 				 ET_min = min(ET_EC_F,ET_EC_EBC,na.rm=T)
# 				 ) %>%
# 	ggplot(aes(x = date)) +
# 	geom_line(aes(y = offset - P_PI_F * scale_factor),
# 						color = "blue") +
# 	geom_ribbon(aes(ymax = ET_max, ymin = ET_min)) +
# 	scale_y_continuous(
# 		name = expression("EC ET (mm "*d^-1*")"),
# 		sec.axis = sec_axis(~ (offset - .)/scale_factor,
# 												name = "Precipitation (mm)")
# 	) +
# 	theme_bw() +
# 	labs(title = paste0(site_id,": Daily EC ET and Precipitation"))
# 
# b <- df %>%
# 	filter(TIMESTAMP >= EC_period[1] & TIMESTAMP <= EC_period[2]) %>%
# 	mutate(date = date(TIMESTAMP)) %>%
# 	group_by(date) %>%
# 	summarise(P_PI_F = sum(P_PI_F),
# 						ET_EC_F = sum(ET_EC_F),
# 						ET_EC_EBC = sum(ET_EC_EBC),
# 						SWC_1_1_1 = mean(SWC_1_1_1,na.rm = TRUE),
# 						SWC_1_1_2 = mean(SWC_1_1_2,na.rm = TRUE),
# 						.groups = "drop"
# 	) %>%
# 	mutate(SWC = rowMeans(across(c(SWC_1_1_1, SWC_1_1_2)), na.rm = TRUE),
# 				 SSM = (SWC - first(SWC))*0.08 * 1000 / 100
# 				 ) %>%
# 	mutate(ET_cum = cumsum(ET_EC_F),
# 				 ET_cum_corr = cumsum(ET_EC_EBC),
# 				 Pr_cum = cumsum(P_PI_F)
# 				 ) %>%
# 	ggplot(aes(date)) +
# 	geom_ribbon(aes(ymax = ET_cum_corr, ymin = ET_cum, color = "EC ET"),alpha = 0.5) +
# 	geom_line(aes(y = Pr_cum,	color = "P")) +
# 	geom_line(aes(y = SSM,color = "SSM")) +
# 	theme_bw() +
#   geom_text(
#     data = ~ dplyr::slice_tail(.x, n = 1),
#     aes(y = Pr_cum, label = round(Pr_cum,1)),
#     hjust = -0.1, color = "#0072B2"
#   ) +
#   geom_text(
#     data = ~ dplyr::slice_tail(.x, n = 1),
#     aes(y = ET_cum_corr, label = round(ET_cum_corr,1)),
#     hjust = -0.1, color = "#000000"
#   ) +
#   geom_text(
#     data = ~ dplyr::slice_tail(.x, n = 1),
#     aes(y = ET_cum, label = round(ET_cum,1)),
#     hjust = -0.1, color = "#000000"
#   ) +
#   geom_text(
#     data = ~ dplyr::slice_tail(.x, n = 1),
#     aes(y = SSM, label = round(SSM,1)),
#     hjust = -0.1, color = "#D55E00"
#   ) +
# 	scale_color_manual(
# 		breaks = c("EC ET", "P","SSM"),
# 		values = c(
# 			"EC ET" = "#000000",        # black
# 			"P" = "#0072B2",        # blue
# 			"SSM" = "#D55E00" # vermillion
# 		)
# 	) +
# 	labs(title = paste0(site_id,": Cumulative EC ET vs precipitation"),
# 			 y = expression("Cumulative flux (mm)"))
# 
# c <- df %>%
# 	filter(TIMESTAMP >= LI710_period[1] & TIMESTAMP <= LI710_period[2]) %>%
# 	filter(TIMESTAMP >= EC_period[1] & TIMESTAMP <= EC_period[2]) %>%
# 	mutate(date = date(TIMESTAMP)) %>%
# 	group_by(date) %>%
# 	summarise_if(is.numeric,sum) %>%
# 	rowwise() %>%
# 	mutate(ET_max = max(ET_710_F,ET_710_EBC,na.rm=T),
# 				 ET_min = min(ET_710_F,ET_710_EBC,na.rm=T)
# 	) %>%
# 	ggplot(aes(x = date)) +
# 	geom_line(aes(y = offset - P_PI_F * scale_factor),
# 						color = "blue") +
# 	geom_ribbon(aes(ymax = ET_max, ymin = ET_min)) +
# 	scale_y_continuous(
# 		name = expression("710 ET (mm "*d^-1*")"),
# 		sec.axis = sec_axis(~ (offset - .)/scale_factor,
# 												name = "Precipitation (mm)")
# 	) +
# 	theme_bw() +
# 	labs(title = paste0(site_id,": Daily 710 ET and Precipitation"))
# 
# d <- df %>%
# 	filter(TIMESTAMP >= LI710_period[1] & TIMESTAMP <= LI710_period[2]) %>%
# 	filter(TIMESTAMP >= EC_period[1] & TIMESTAMP <= EC_period[2]) %>%
# 	mutate(date = date(TIMESTAMP)) %>%
# 	group_by(date) %>%
# 	summarise(P_PI_F = sum(P_PI_F),
# 						ET_EC_F = sum(ET_710_F),
# 						ET_EC_EBC = sum(ET_710_EBC),
# 						SWC_1_1_1 = mean(SWC_1_1_1,na.rm = TRUE),
# 						SWC_1_1_2 = mean(SWC_1_1_2,na.rm = TRUE),
# 						.groups = "drop"
# 	) %>%
# 	mutate(SWC = rowMeans(across(c(SWC_1_1_1, SWC_1_1_2)), na.rm = TRUE),
# 				 SSM = (SWC - first(SWC))*0.08 * 1000 / 100
# 	) %>%
# 	mutate(ET_cum = cumsum(ET_EC_F),
# 				 ET_cum_corr = cumsum(ET_EC_EBC),
# 				 Pr_cum = cumsum(P_PI_F)
# 	) %>%
# 	ggplot(aes(date)) +
# 	geom_ribbon(aes(ymax = ET_cum_corr, ymin = ET_cum, color = "710 ET"),alpha = 0.5) +
# 	geom_line(aes(y = Pr_cum,	color = "P")) +
# 	geom_line(aes(y = SSM,color = "SSM")) +
# 	theme_bw() +
#   geom_text(
#     data = ~ dplyr::slice_tail(.x, n = 1),
#     aes(y = Pr_cum, label = round(Pr_cum,1)),
#     hjust = -0.1, color = "#0072B2"
#   ) +
#   geom_text(
#     data = ~ dplyr::slice_tail(.x, n = 1),
#     aes(y = ET_cum_corr, label = round(ET_cum_corr,1)),
#     hjust = -0.1, color = "#000000"
#   ) +
#   geom_text(
#     data = ~ dplyr::slice_tail(.x, n = 1),
#     aes(y = ET_cum, label = round(ET_cum,1)),
#     hjust = -0.1, color = "#000000"
#   ) +
#   geom_text(
#     data = ~ dplyr::slice_tail(.x, n = 1),
#     aes(y = SSM, label = round(SSM,1)),
#     hjust = -0.1, color = "#D55E00"
#   ) +
# 	scale_color_manual(
# 		breaks = c("710 ET", "P","SSM"),
# 		values = c(
# 			"710 ET" = "#000000",        # black
# 			"P" = "#0072B2",        # blue
# 			"SSM" = "#D55E00" # vermillion
# 		)
# 	) +
# 	labs(title = paste0(site_id,": Cumulative 710 ET vs precipitation"),
# 			 y = expression("Cumulative flux (mm)"))
# 
# ggsave(plot = plot_grid(a,b,c,d,ncol = 1),
# 			 filename = paste0(pipeline_path(base_dir, dirs_use$dir_output),"/figures/",site_id,"_Level_4_Figure_2_water_balance.jpeg"),
# 			 width = 11, height = 13, dpi = 150
# )
# 
# 
# ################################################################################
# #### Diagnostic figure 3: ET methods intercomparison  
# ################################################################################
# a <- df %>%
# 	ggplot(aes(ET_EC,ET_710)) +
# 	theme_bw() +
# 	geom_abline(lty = 2) +
# 	geom_point(alpha = 0.2) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "710 ET (mm/30min)",
# 			 x = "EC ET (mm/30min)",
# 			 title = paste0(site_id,": EC vs 710")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
# 	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)
# 
# b <- df %>%
# 	filter(!is.na(ET_710)) %>%
# 	ggplot(aes(ET_EC,ET_710_EBC)) +
# 	theme_bw() +
# 	geom_abline(lty = 2) +
# 	geom_point(alpha = 0.2) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "corrected 710 ET (mm/30min)",
# 			 x = "EC ET (mm/30min)",
# 			 title = paste0(site_id,": EC vs EBC 710")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
# 	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)
# 
# c <- df %>%
# 	ggplot(aes(ET_EC,ET_VBR)) +
# 	theme_bw() +
# 	geom_abline(lty = 2) +
# 	geom_point(alpha = 0.2) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "VBR ET (mm/30min)",
# 			 x = "EC ET (mm/30min)",
# 			 title = paste0(site_id,": EC vs VBR")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
# 	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)
# 
# d <- df %>%
# 	ggplot(aes(ET_EC,ET_RE)) +
# 	theme_bw() +
# 	geom_abline(lty = 2) +
# 	geom_point(alpha = 0.2) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "RE ET (mm/30min)",
# 			 x = "EC ET (mm/30min)",
# 			 title = paste0(site_id,": EC vs RE")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
# 	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)
# 
# e <- df %>%
# 	filter(!is.na(ET_EC)) %>%
# 	ggplot(aes(ET_EC_EBC,ET_710)) +
# 	theme_bw() +
# 	geom_abline(lty = 2) +
# 	geom_point(alpha = 0.2) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "710 ET (mm/30min)",
# 			 x = "corrected EC ET (mm/30min)",
# 			 title = paste0(site_id,": EBC EC vs 710")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
# 	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)
# 
# f <- df %>%
# 	filter(!is.na(ET_EC) & !is.na(ET_710)) %>%
# 	ggplot(aes(ET_EC_EBC,ET_710_EBC)) +
# 	theme_bw() +
# 	geom_abline(lty = 2) +
# 	geom_point(alpha = 0.2) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "corrected 710 ET (mm/30min)",
# 			 x = "corrected EC ET (mm/30min)",
# 			 title = paste0(site_id,": EBC EC vs EBC 710")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
# 	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)
# 
# g <- df %>%
# 	filter(!is.na(ET_EC)) %>%
# 	ggplot(aes(ET_EC_EBC,ET_VBR)) +
# 	theme_bw() +
# 	geom_abline(lty = 2) +
# 	geom_point(alpha = 0.2) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "VBR ET (mm/30min)",
# 			 x = "corrected EC ET (mm/30min)",
# 			 title = paste0(site_id,": EBC EC vs VBR")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
# 	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)
# 
# h <- df %>%
# 	filter(!is.na(ET_EC)) %>%
# 	ggplot(aes(ET_EC_EBC,ET_RE)) +
# 	theme_bw() +
# 	geom_abline(lty = 2) +
# 	geom_point(alpha = 0.2) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "RE ET (mm/30min)",
# 			 x = "corrected EC ET (mm/30min)",
# 			 title = paste0(site_id,": EBC EC vs RE")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
# 	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)
# 
# p1 <- plot_grid(a,b,c,d,e,f,g,h,ncol = 4)
# 
# 
# a <- df %>%
# 	mutate(date = date(TIMESTAMP)) %>%
# 	group_by(date) %>%
# 	summarise_if(is.numeric,sum) %>%
# 	ggplot(aes(ET_EC_F,ET_710_F)) +
# 	theme_bw() +
# 	geom_abline(lty = 2) +
# 	geom_point(alpha = 1) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "710 ET (mm/day)",
# 			 x = "EC ET (mm/day)",
# 			 title = paste0(site_id,": EC vs 710")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
# 	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)
# 
# b <- df %>%
# 	mutate(date = date(TIMESTAMP)) %>%
# 	group_by(date) %>%
# 	summarise_if(is.numeric,sum) %>%
# 	ggplot(aes(ET_EC_F,ET_710_EBC)) +
# 	theme_bw() +
# 	geom_abline(lty = 2) +
# 	geom_point(alpha = 1) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "corrected 710 ET (mm/day)",
# 			 x = "EC ET (mm/day)",
# 			 title = paste0(site_id,": EC vs EBC 710")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
# 	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)
# 
# c <- df %>%
# 	mutate(date = date(TIMESTAMP)) %>%
# 	group_by(date) %>%
# 	summarise_if(is.numeric,sum) %>%
# 	ggplot(aes(ET_EC_F,ET_VBR_F)) +
# 	theme_bw() +
# 	geom_abline(lty = 2) +
# 	geom_point(alpha = 1) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "VBR ET (mm/day)",
# 			 x = "EC ET (mm/day)",
# 			 title = paste0(site_id,": EC vs VBR")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
# 	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)
# 
# d <- df %>%
# 	mutate(date = date(TIMESTAMP)) %>%
# 	group_by(date) %>%
# 	summarise_if(is.numeric,sum) %>%
# 	ggplot(aes(ET_EC_F,ET_RE_F)) +
# 	theme_bw() +
# 	geom_abline(lty = 2) +
# 	geom_point(alpha = 1) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "RE ET (mm/day)",
# 			 x = "EC ET (mm/day)",
# 			 title = paste0(site_id,": EC vs RE")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
# 	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)
# 
# e <- df %>%
# 	mutate(date = date(TIMESTAMP)) %>%
# 	group_by(date) %>%
# 	summarise_if(is.numeric,sum) %>%
# 	ggplot(aes(ET_EC_EBC,ET_710_F)) +
# 	theme_bw() +
# 	geom_abline(lty = 2) +
# 	geom_point(alpha = 1) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "710 ET (mm/day)",
# 			 x = "corrected EC ET (mm/day)",
# 			 title = paste0(site_id,": EBC EC vs 710")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
# 	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)
# 
# f <- df %>%
# 	mutate(date = date(TIMESTAMP)) %>%
# 	group_by(date) %>%
# 	summarise_if(is.numeric,sum) %>%
# 	ggplot(aes(ET_EC_EBC,ET_710_EBC)) +
# 	theme_bw() +
# 	geom_abline(lty = 2) +
# 	geom_point(alpha = 1) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "corrected 710 ET (mm/day)",
# 			 x = "corrected EC ET (mm/day)",
# 			 title = paste0(site_id,"EBC EC vs EBC 710")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
# 	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)
# 
# g <- df %>%
# 	mutate(date = date(TIMESTAMP)) %>%
# 	group_by(date) %>%
# 	summarise_if(is.numeric,sum) %>%
# 	ggplot(aes(ET_EC_EBC,ET_VBR_F)) +
# 	theme_bw() +
# 	geom_abline(lty = 2) +
# 	geom_point(alpha = 1) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "VBR ET (mm/day)",
# 			 x = "corrected EC ET (mm/day)",
# 			 title = paste0(site_id,": EBC EC vs VBR")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
# 	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)
# 
# h <- df %>%
# 	mutate(date = date(TIMESTAMP)) %>%
# 	group_by(date) %>%
# 	summarise_if(is.numeric,sum) %>%
# 	ggplot(aes(ET_EC_EBC,ET_RE_F)) +
# 	theme_bw() +
# 	geom_abline(lty = 2) +
# 	geom_point(alpha = 1) +
# 	geom_vline(xintercept = 0,lty = 2) + geom_hline(yintercept = 0,lty = 2) +
# 	geom_smooth(method = "lm", formula = "y ~ x + 0") +
# 	labs(y = "RE ET (mm/day)",
# 			 x = "corrected EC ET (mm/day)",
# 			 title = paste0(site_id,": EBC EC vs RE")) +
# 	stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
# 	stat_poly_eq(use_label("adj.R2"), formula = y ~ x,vjust = 1.7)
# 
# p2 <- plot_grid(a,b,c,d,e,f,g,h,ncol = 4)
# 
# ggsave(plot = p1,
# 			 filename = paste0(pipeline_path(base_dir, dirs_use$dir_output),"/figures/",site_id,"_Level_4_Figure_3_ET_intercomparison.jpeg"),
# 			 width = 11, height = 6.5, dpi = 150
# )
# 
# ggsave(plot = p2,
#        filename = paste0(pipeline_path(base_dir, dirs_use$dir_output),"/figures/",site_id,"_Level_4_Figure_4_ET_intercomparison_daily.jpeg"),
#        width = 11, height = 6.5, dpi = 150
# )

################################################################################
#### save
################################################################################

write_custom_csv <- pipeline_write

# units for new variables
df_unit_addition <- 
	data.frame(LE_filtered1 = c("[W+1m-2]","derived"),H_filtered1 = c("[W+1m-2]","derived"),
						 FC_filtered1 = c("[µmol+1s-1m-2]","derived"),LE_filtered2 = c("[W+1m-2]","derived"),
						 H_filtered2 = c("[W+1m-2]","derived"), FC_filtered2 = c("[µmol+1s-1m-2]","derived"),
						 LE_710_filtered = c("[W+1m-2]","derived"),H_710_filtered = c("[W+1m-2]","derived"),
						 LE_strg_added = c("[W+1m-2]","derived"),H_strg_added = c("[W+1m-2]","derived"),        
						 FC_strg_added = c("[µmol+1s-1m-2]","derived"), 
						 NETRAD_PI_F = c("[W+1m-2]","derived"), G_PI_F = c("[W+1m-2]","derived"),
						 LE_PI_F = c("[W+1m-2]","derived"), LE_MDS = c("[W+1m-2]","derived"),
						 H_PI_F = c("[W+1m-2]","derived"),FC_PI_F = c("[µmol+1s-1m-2]","derived"),
						 LE_710_PI_F = c("[W+1m-2]","derived"), LE_710_MDS = c("[W+1m-2]","derived"), H_710_PI_F = c("[W+1m-2]","derived"),
						 EBR_d = c("ratio","derived"), LE_PI_CORR = c("[W+1m-2]","derived"),H_PI_CORR = c("[W+1m-2]","derived"),
						 EBR_d_710 = c("ratio","derived"), LE_710_PI_CORR = c("[W+1m-2]","derived"),H_710_PI_CORR = c("[W+1m-2]","derived"),      
						 VBR = c("-","derived"),VBR_Ts = c("-","derived"),
						 LE_VBR = c("[W+1m-2]","derived"),LE_VBR_Ts = c("[W+1m-2]","derived"),
						 LE_RE = c("[W+1m-2]","derived"),
						 LE_RE_F = c("[W+1m-2]","derived"), LE_RE_MDS = c("[W+1m-2]","derived"),
						 LE_VBR_F = c("[W+1m-2]","derived"),  LE_VBR_MDS = c("[W+1m-2]","derived"),        
						 LE_VBR_Ts_F = c("[W+1m-2]","derived"),  LE_VBR_Ts_MDS = c("[W+1m-2]","derived"),
						 ET_EC = c("mm","derived"),ET_EC_F = c("mm","derived"),ET_EC_EBC = c("mm","derived"),
						 ET_710 = c("mm","derived"),ET_710_F = c("mm","derived"), ET_710_EBC = c("mm","derived"),
						 ET_VBR = c("mm","derived"),ET_VBR_F = c("mm","derived"),
						 ET_VBR_Ts = c("mm","derived"),ET_VBR_Ts_F = c("mm","derived"),
						 ET_RE = c("mm","derived"),ET_RE_F = c("mm","derived")
	)
df_unit_final <- cbind(df_unit,df_unit_addition)

print(which(colnames(df) != colnames(df_unit_final)))

write_custom_csv(
	df = df,
	units_df = df_unit_final,
	file = paste0(pipeline_path(base_dir, dirs_use$dir_output),"/",site_id,"_Level_4_post_processed_data_",date(Sys.time()),".csv"),
	description = "Fully processed half hourly data"
)

}
