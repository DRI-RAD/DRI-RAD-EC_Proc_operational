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

df <- left_join(df,ec, by = 'TIMESTAMP')
df_unit <- left_join(df_unit, ec_unit, by = 'TIMESTAMP')

################################################################################
#### Load LI710 
################################################################################

{ # LI710 is optional; the runner supplies an NA placeholder when absent.

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
pipeline_validate_mds(EddyData$TIMESTAMP)
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
pipeline_validate_mds(EddyData$TIMESTAMP)
EProc <- sEddyProc$new(
	site_id, EddyDataWithPosix, c('NEE','Rg','Tair','VPD', 'Ustar', 'LE', 'H','NETRAD',"G_PI"))

#++ Fill NEE gaps with MDS gap filling algorithm (without prior ustar filtering)
EProc$sMDSGapFill('LE', FillAll = T)
EProc$sMDSGapFill('H', FillAll = F)
EProc$sMDSGapFill('NEE', FillAll = F)
# Record actual gap-fill quality without changing scientific outputs.
pipeline_capture_mds(EddyData, EProc$sExportResults(), c(LE = "LE", H = "H", NEE = "FC"))

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
	mutate(Year = year(TIMESTAMP),
				 DoY = yday(TIMESTAMP),
				 Hour = hour(TIMESTAMP) + minute(TIMESTAMP)/60,
				 LE = as.numeric(LE_710_filtered),
				 H = as.numeric(H_710_filtered),
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
pipeline_validate_mds(EddyData$TIMESTAMP)
#++ Fill NEE gaps with MDS gap filling algorithm (without prior ustar filtering)
# Preserve missing LI710 channels rather than attempting MDS with no observations.
mds_li <- as.data.frame(matrix(NA_real_, nrow = nrow(EddyData), ncol = 5))
names(mds_li) <- c("LE_f", "LE_fall", "H_f", "LE_fqc", "H_fqc")
if (any(is.finite(EddyData$LE)) || any(is.finite(EddyData$H))) {
EProc <- sEddyProc$new(
	site_id, EddyDataWithPosix, c('Rg','Tair','VPD', 'Ustar', 'LE', 'H'))
if (any(is.finite(EddyData$LE))) EProc$sMDSGapFill('LE', FillAll = TRUE)
if (any(is.finite(EddyData$H))) EProc$sMDSGapFill('H', FillAll = FALSE)
export_li <- EProc$sExportResults()
for (name in intersect(names(mds_li), names(export_li))) mds_li[[name]] <- export_li[[name]]
if ("LE_fall_qc" %in% names(export_li))
  mds_li$LE_fall <- ifelse(export_li$LE_fall_qc == 1, mds_li$LE_fall, NA_real_)
}
pipeline_capture_mds(EddyData, mds_li, c(LE = "LE_710", H = "H_710"))
MDSout <- data.frame(TIMESTAMP = EddyData$TIMESTAMP, LE_710_PI_F = mds_li$LE_f,
                     LE_710_MDS = mds_li$LE_fall, H_710_PI_F = mds_li$H_f)
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
# VBR without the _Ts suffix is TEST ONLY; not used for the study analysis.
# VBR_Ts was used for analysis. Keep each MDS estimate and its QC mask separate.
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
pipeline_validate_mds(EddyData$TIMESTAMP)
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
										          LE_VBR_Ts_fall = ifelse(LE_VBR_Ts_fall_qc == 1, LE_VBR_Ts_fall, NA),
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


################################################################################
#### Diagnostic figure 1a: Energy balance
################################################################################
# Plot only newly saved rows by default; calculations retain all MDS context.
# Optional soil-water curves can be absent at sites without those sensors.
# Period and cumulative figures are rendered during transactional publication.
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
