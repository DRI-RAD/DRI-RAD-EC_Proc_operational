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

# signal strength
use <- met %>% select(TIMESTAMP, CO2_sig_strgth_Min, H2O_sig_strgth_Min,PotRad)
use_unit <- met_unit %>% select(TIMESTAMP, CO2_sig_strgth_Min, H2O_sig_strgth_Min,PotRad)

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

met_unit <- read_csv(latest_file, skip = 1, n_max = 2)
met <- read_csv(latest_file, skip = 4, col_names = FALSE)
colnames(met) <- colnames(met_unit)
# Restrict input to the requested calculation window before processing.
met <- pipeline_slice(met)

# Gap-filled vars for day night partitioning and ustar threshold
met <- met %>% select(TIMESTAMP, SW_IN_PI_F,TA_PI_F,RH_PI_F,VPD_PI_F)
met_unit <- met_unit %>% select(TIMESTAMP, SW_IN_PI_F,TA_PI_F,RH_PI_F,VPD_PI_F)

met <- left_join(met,use, by = 'TIMESTAMP')
met_unit <- left_join(met_unit, use_unit, by = 'TIMESTAMP')


################################################################################
#### Load Easyflux output
################################################################################
Easyflux_unit <- read_delim(pipeline_path(base_dir, dirs_use$dir_Easyflux),skip = 1, n_max = 2)
Easyflux <- read_delim(pipeline_path(base_dir, dirs_use$dir_Easyflux),skip = 4, col_names = F)
colnames(Easyflux) <- colnames(Easyflux_unit)
# Restrict input to the requested calculation window before processing.
Easyflux <- pipeline_slice(Easyflux)

# Full time sequence
Easyflux <- Easyflux %>%
  distinct(TIMESTAMP, .keep_all = TRUE)
TIMESTAMP <- seq(from = min(Easyflux$TIMESTAMP), to = max(Easyflux$TIMESTAMP), 60*30)
Easyflux <- left_join(data.frame(TIMESTAMP), Easyflux)

# select Easyflux output etc.
need_to_select <- 
  c('TIMESTAMP', 'LE', 'LE_SSITC_TEST', 'H', 'H_SSITC_TEST', 'FC', 'FC_SSITC_TEST', 'TAU', 'TAU_SSITC_TEST',
    'CO2','CO2_SIGMA','H2O','H2O_SIGMA',"T_SONIC","T_SONIC_SIGMA",
    'WS',"WS_MAX",'WD','USTAR',"ZL",'MO_LENGTH', "U_SIGMA","V_SIGMA","W_SIGMA",
    "FETCH_MAX","FETCH_90","FETCH_55","FETCH_40")


Easyflux <- Easyflux %>%
  select(need_to_select)
Easyflux_unit <- Easyflux_unit %>%
  select(need_to_select)


# EasyFlux QC: Foken et al. 2012 1-9 grade
# Convert to simplified 0-2 scale comparable to Mauder & Foken 2004 / AmeriFlux QC
if(max(Easyflux$LE_SSITC_TEST,na.rm=T) > 2){
  # original method
  # Easyflux <- Easyflux %>%
  #   mutate(H_SSITC_TEST = case_when(H_SSITC_TEST %in% 1:3 ~ 0, H_SSITC_TEST %in% 4:6 ~ 1, H_SSITC_TEST %in% 7:9 ~ 2, TRUE ~ NA_real_),
  #          LE_SSITC_TEST = case_when(LE_SSITC_TEST %in% 1:3 ~ 0, LE_SSITC_TEST %in% 4:6 ~ 1, LE_SSITC_TEST %in% 7:9 ~ 2, TRUE ~ NA_real_),
  #          FC_SSITC_TEST = case_when(FC_SSITC_TEST %in% 1:3 ~ 0, FC_SSITC_TEST %in% 4:6 ~ 1, FC_SSITC_TEST %in% 7:9 ~ 2, TRUE ~ NA_real_),
  #          TAU_SSITC_TEST = case_when(TAU_SSITC_TEST %in% 1:3 ~ 0, TAU_SSITC_TEST %in% 4:6 ~ 1, TAU_SSITC_TEST %in% 7:9 ~ 2, TRUE ~ NA_real_)
  #   )
  # 
  # more conservative method (as we use despike SSITC can be use conservatively)
  Easyflux <- Easyflux %>%
    mutate(H_SSITC_TEST = case_when(H_SSITC_TEST %in% 1:3 ~ 0, H_SSITC_TEST %in% 4:8 ~ 1, H_SSITC_TEST %in% 9 ~ 2, TRUE ~ NA_real_),
           LE_SSITC_TEST = case_when(LE_SSITC_TEST %in% 1:3 ~ 0, LE_SSITC_TEST %in% 4:8 ~ 1, LE_SSITC_TEST %in% 9 ~ 2, TRUE ~ NA_real_),
           FC_SSITC_TEST = case_when(FC_SSITC_TEST %in% 1:3 ~ 0, FC_SSITC_TEST %in% 4:8 ~ 1, FC_SSITC_TEST %in% 9 ~ 2, TRUE ~ NA_real_),
           TAU_SSITC_TEST = case_when(TAU_SSITC_TEST %in% 1:3 ~ 0, TAU_SSITC_TEST %in% 4:8 ~ 1, TAU_SSITC_TEST %in% 9 ~ 2, TRUE ~ NA_real_)
    )
}


add_suffix <- 
  c('TIMESTAMP', 'LE_EasyFlux', 'LE_SSITC_TEST_EasyFlux', 'H_EasyFlux', 'H_SSITC_TEST_EasyFlux', 
    'FC_EasyFlux', 'FC_SSITC_TEST_EasyFlux', 'TAU_EasyFlux', 'TAU_SSITC_TEST_EasyFlux',
    'CO2_EasyFlux','CO2_SIGMA_EasyFlux','H2O_EasyFlux','H2O_SIGMA_EasyFlux',"T_SONIC_EasyFlux","T_SONIC_SIGMA_EasyFlux",
    'WS_EasyFlux',"WS_MAX_EasyFlux",'WD_EasyFlux','USTAR_EasyFlux',"ZL_EasyFlux",'MO_LENGTH_EasyFlux', 
    "U_SIGMA_EasyFlux","V_SIGMA_EasyFlux","W_SIGMA_EasyFlux",
    "FETCH_MAX_EasyFlux","FETCH_90_EasyFlux","FETCH_55_EasyFlux","FETCH_40_EasyFlux")

colnames(Easyflux) <- add_suffix
colnames(Easyflux_unit) <- add_suffix

head(Easyflux)
head(Easyflux_unit)

df <- met %>%
  left_join(Easyflux, by = "TIMESTAMP")
df_unit <- met_unit %>%
  left_join(Easyflux_unit, by = "TIMESTAMP")

################################################################################
#### Load EddyPro output 
################################################################################
# Read site/month_year/output/eddypro_*_full_output_*.csv by field names.
# The runner validates and loads the files once, then shares the selected window.
if(TRUE){
  EP <- pipeline_slice(pipeline_options$eddypro$data)
  EP_unit <- pipeline_options$eddypro$units

  # sigma calculation and convert sigma unit to mole fraction instead of molar density.
  EP <- EP %>%
    mutate(co2_sigma = sqrt(co2_var) * co2_mole_fraction / co2_molar_density, 
           h2o_sigma = sqrt(h2o_var) * h2o_mole_fraction / h2o_molar_density,
           ts_sigma = sqrt(ts_var),
           u_sigma = sqrt(u_var),
           v_sigma = sqrt(v_var),
           w_sigma = sqrt(w_var)
    )
  
  EP_unit <- EP_unit %>%
    mutate(co2_sigma = co2_mole_fraction, 
           h2o_sigma = h2o_mole_fraction,
           ts_sigma = sonic_temperature,
           u_sigma = wind_speed,
           v_sigma = wind_speed,
           w_sigma = wind_speed
    )
  
  need_to_select <- 
    c('TIMESTAMP', 'LE', 'qc_LE', 'H', 'qc_H', 'co2_flux', 'qc_co2_flux', 'Tau', 'qc_Tau',
      'co2_mole_fraction','co2_sigma','h2o_mole_fraction','h2o_sigma',"sonic_temperature","ts_sigma",
      'wind_speed',"max_wind_speed",'wind_dir','u*',"(z-d)/L",'L', "u_sigma","v_sigma","w_sigma",
      "x_90%","x_70%","x_50%","x_30%","x_10%",
      "rand_err_LE","rand_err_H","rand_err_co2_flux","rand_err_Tau",
      "LE_strg","H_strg","co2_strg"
      )
  
  EP <- EP %>% 
    select(need_to_select)
  EP_unit <- EP_unit %>% 
    select(need_to_select)
  EP_unit[2,] <- NA
  
  add_suffix <- 
    c('TIMESTAMP', 'LE_EddyPro', 'LE_SSITC_TEST_EddyPro', 'H_EddyPro', 'H_SSITC_TEST_EddyPro', 
      'FC_EddyPro', 'FC_SSITC_TEST_EddyPro', 'TAU_EddyPro', 'TAU_SSITC_TEST_EddyPro',
      'CO2_EddyPro','CO2_SIGMA_EddyPro','H2O_EddyPro','H2O_SIGMA_EddyPro',"T_SONIC_EddyPro","T_SONIC_SIGMA_EddyPro",
      'WS_EddyPro',"WS_MAX_EddyPro",'WD_EddyPro','USTAR_EddyPro',"ZL_EddyPro",'MO_LENGTH_EddyPro', 
      "U_SIGMA_EddyPro","V_SIGMA_EddyPro","W_SIGMA_EddyPro",
      "FETCH_90_EddyPro","FETCH_70_EddyPro","FETCH_50_EddyPro", "FETCH_30_EddyPro", "FETCH_10_EddyPro",
      'LE_rand_err_EddyPro','H_rand_err_EddyPro','FC_rand_err_EddyPro','TAU_rand_err_EddyPro',
      'SLE_EddyPro','SH_EddyPro','SC_EddyPro'  
      )
  
  colnames(EP) <- add_suffix
  colnames(EP_unit) <- add_suffix
  
  df <- df %>%
    left_join(EP, by = "TIMESTAMP")
  df_unit <- df_unit %>%
    left_join(EP_unit, by = "TIMESTAMP")
  
  
  # comparison figures
  # physical range 
  if(site_id %in% c('EDVG', 'ERVA')){
    # physical range for agriculture
    LE_range = c(-30, 1100)
    H_range = c(-300, 1100)
    FC_range = c(-60,30)
    
  } else{
    # physical range for phreatophytes
    LE_range = c(-30, 600)
    H_range = c(-100, 1100)
    FC_range = c(-30,10)
  }
 
  a <- pipeline_figure_data(df) %>%
    filter(LE_EasyFlux < LE_range[2] &  LE_EasyFlux >= LE_range[1] & 
             LE_EddyPro < LE_range[2] &  LE_EddyPro >= LE_range[1]) %>%
    ggplot(aes(LE_EasyFlux,LE_EddyPro)) +
    geom_point() +
    geom_abline(lty = 2) +
    theme_bw() +
    geom_smooth(method = "lm",formula = 'y~x+0') +
    stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
    labs(x = "EasyFlux LE (W m-2)",
         y = "EddyPro LE (W m-2)",
         title = paste0(site_id, ": LE")
    )
  
  b <- pipeline_figure_data(df) %>%
    filter(H_EasyFlux < H_range[2] &  H_EasyFlux >= H_range[1] & 
             H_EddyPro < H_range[2] &  H_EddyPro >= H_range[1]) %>%
    ggplot(aes(H_EasyFlux,H_EddyPro)) +
    geom_point() +
    geom_abline(lty = 2) +
    theme_bw() +
    geom_smooth(method = "lm",formula = 'y~x+0') +
    stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
    labs(x = "EasyFlux H (W m-2)",
         y = "EddyPro H (W m-2)",
         title = paste0(site_id, ": H")
    )
    
  
  c <- pipeline_figure_data(df) %>%
    filter(FC_EasyFlux < FC_range[2] &  FC_EasyFlux >= FC_range[1] & 
             FC_EddyPro < FC_range[2] &  FC_EddyPro >= FC_range[1]) %>%
    ggplot(aes(FC_EasyFlux,FC_EddyPro)) +
    geom_point() +
    geom_abline(lty = 2) +
    theme_bw() +
    geom_smooth(method = "lm",formula = 'y~x+0') +
    stat_poly_eq(use_label("eq"), formula = y ~ x + 0) +
    labs(x = "EasyFlux FC (umolCO2 m-2 s-1)",
         y = "EddyPro FC (umolCO2 m-2 s-1)",
         title = paste0(site_id, ": FC")
    )
  
  
  ggsave(plot = plot_grid(a,b,c,ncol = 3),
         filename = paste0(pipeline_path(base_dir, dirs_use$dir_output),"/figures/",site_id,"_Level_3_Eddypro_vs_EasyFlux.jpeg"),
         width = 11, height = 4, dpi = 150
  )
  
  
}


#################################################
#### combined EC file
#################################################
head(df)
head(df_unit)

write_custom_csv <- pipeline_write


write_custom_csv(
  df = df %>% select( -c(SW_IN_PI_F,TA_PI_F,RH_PI_F,VPD_PI_F,CO2_sig_strgth_Min,H2O_sig_strgth_Min,PotRad)),
  units_df = df_unit %>% select( -c(SW_IN_PI_F,TA_PI_F,RH_PI_F,VPD_PI_F,CO2_sig_strgth_Min,H2O_sig_strgth_Min,PotRad)),
  file = paste0(pipeline_path(base_dir, dirs_use$dir_output),"/",site_id,"_Level_3_EC_data_",date(Sys.time()),".csv"),
  description = "Data from EasyFlux output have been merged with EddyPro output."
)

################################################################################
#### Use Eddypro if available, otherwise use Easyflux
################################################################################
eddy_cols <- names(df)[grepl("_EddyPro$", names(df))]

EC <- df %>%
  mutate(
    processing = if (length(eddy_cols) == 0) {
      "EasyFlux"
    } else {
      if_else(
        rowSums(!is.na(across(all_of(eddy_cols)))) > 0,
        "EddyPro",
        "EasyFlux"
      )
    }
  )

EC_unit <- df_unit %>% mutate(processing = NA)

# final variables
EC <- EC %>%
  mutate(LE = ifelse(processing == "EddyPro", LE_EddyPro,LE_EasyFlux),
         LE_SSITC_TEST = ifelse(processing == "EddyPro", LE_SSITC_TEST_EddyPro,LE_SSITC_TEST_EasyFlux),
         H = ifelse(processing == "EddyPro", H_EddyPro,H_EasyFlux),
         H_SSITC_TEST = ifelse(processing == "EddyPro", H_SSITC_TEST_EddyPro,H_SSITC_TEST_EasyFlux),
         FC = ifelse(processing == "EddyPro", FC_EddyPro,FC_EasyFlux),
         FC_SSITC_TEST = ifelse(processing == "EddyPro", FC_SSITC_TEST_EddyPro,FC_SSITC_TEST_EasyFlux),
         TAU = ifelse(processing == "EddyPro", TAU_EddyPro,TAU_EasyFlux),
         TAU_SSITC_TEST = ifelse(processing == "EddyPro", TAU_SSITC_TEST_EddyPro,TAU_SSITC_TEST_EasyFlux),
         CO2 = ifelse(processing == "EddyPro", CO2_EddyPro,CO2_EasyFlux),
         CO2_SIGMA = ifelse(processing == "EddyPro", CO2_SIGMA_EddyPro,CO2_SIGMA_EasyFlux),
         H2O = ifelse(processing == "EddyPro", H2O_EddyPro,H2O_EasyFlux),
         H2O_SIGMA = ifelse(processing == "EddyPro", H2O_SIGMA_EddyPro,H2O_SIGMA_EasyFlux),
         T_SONIC = ifelse(processing == "EddyPro", T_SONIC_EddyPro,T_SONIC_EasyFlux),
         T_SONIC_SIGMA = ifelse(processing == "EddyPro", T_SONIC_SIGMA_EddyPro,T_SONIC_SIGMA_EasyFlux),
         WS = ifelse(processing == "EddyPro", WS_EddyPro,WS_EasyFlux),
         WS_MAX = ifelse(processing == "EddyPro", WS_MAX_EddyPro,WS_MAX_EasyFlux),
         WD = ifelse(processing == "EddyPro", WD_EddyPro,WD_EasyFlux),
         USTAR = ifelse(processing == "EddyPro", USTAR_EddyPro,USTAR_EasyFlux),
         ZL = ifelse(processing == "EddyPro", ZL_EddyPro,ZL_EasyFlux),
         MO_LENGTH = ifelse(processing == "EddyPro", MO_LENGTH_EddyPro,MO_LENGTH_EasyFlux),
         U_SIGMA = ifelse(processing == "EddyPro", U_SIGMA_EddyPro,U_SIGMA_EasyFlux),
         V_SIGMA = ifelse(processing == "EddyPro", V_SIGMA_EddyPro,V_SIGMA_EasyFlux),
         W_SIGMA = ifelse(processing == "EddyPro", W_SIGMA_EddyPro,W_SIGMA_EasyFlux),
         FETCH_90 = ifelse(processing == "EddyPro", FETCH_90_EddyPro,FETCH_90_EasyFlux),
         FETCH_70 = ifelse(processing == "EddyPro", FETCH_70_EddyPro,NA),
         SLE = ifelse(processing == "EddyPro", SLE_EddyPro,NA),
         SH = ifelse(processing == "EddyPro", SH_EddyPro,NA),
         SC = ifelse(processing == "EddyPro", SC_EddyPro,NA),
         ) %>%
  select(c('TIMESTAMP','SW_IN_PI_F','TA_PI_F','RH_PI_F','VPD_PI_F','CO2_sig_strgth_Min','H2O_sig_strgth_Min','PotRad',
           'processing','LE', 'LE_SSITC_TEST', 'H', 'H_SSITC_TEST','FC', 'FC_SSITC_TEST', 'TAU', 'TAU_SSITC_TEST',
           'CO2','CO2_SIGMA','H2O','H2O_SIGMA',"T_SONIC","T_SONIC_SIGMA",
           'WS',"WS_MAX",'WD','USTAR',"ZL",'MO_LENGTH', "U_SIGMA","V_SIGMA","W_SIGMA",
           "FETCH_90", "FETCH_70",'SLE','SH','SC'))



EC_unit <- EC_unit %>%
  select(c('TIMESTAMP','SW_IN_PI_F','TA_PI_F','RH_PI_F','VPD_PI_F','CO2_sig_strgth_Min','H2O_sig_strgth_Min','PotRad',
           'processing','LE_EasyFlux', 'LE_SSITC_TEST_EasyFlux', 'H_EasyFlux', 'H_SSITC_TEST_EasyFlux',
           'FC_EasyFlux', 'FC_SSITC_TEST_EasyFlux', 'TAU_EasyFlux', 'TAU_SSITC_TEST_EasyFlux',
           'CO2_EasyFlux','CO2_SIGMA_EasyFlux','H2O_EasyFlux','H2O_SIGMA_EasyFlux',"T_SONIC_EasyFlux","T_SONIC_SIGMA_EasyFlux",
           'WS_EasyFlux',"WS_MAX_EasyFlux",'WD_EasyFlux','USTAR_EasyFlux',"ZL_EasyFlux",'MO_LENGTH_EasyFlux', 
           "U_SIGMA_EasyFlux","V_SIGMA_EasyFlux","W_SIGMA_EasyFlux",
           "FETCH_90_EasyFlux")) 

colnames(EC_unit) <- c('TIMESTAMP','SW_IN_PI_F','TA_PI_F','RH_PI_F','VPD_PI_F','CO2_sig_strgth_Min','H2O_sig_strgth_Min','PotRad',
                       'processing','LE', 'LE_SSITC_TEST', 'H', 'H_SSITC_TEST','FC', 'FC_SSITC_TEST', 'TAU', 'TAU_SSITC_TEST',
                       'CO2','CO2_SIGMA','H2O','H2O_SIGMA',"T_SONIC","T_SONIC_SIGMA",
                       'WS',"WS_MAX",'WD','USTAR',"ZL",'MO_LENGTH', "U_SIGMA","V_SIGMA","W_SIGMA",
                       "FETCH_90")

EC_unit <- EC_unit %>%
  mutate(FETCH_70 = FETCH_90,
         SLE = LE,
         SH = H,
         SC = FC
         )

################################################################################
#### QC column
################################################################################
# SSITC: 0 good, 1 acceptable, 2 bad
# signal strength: 0 good, 1 bad (need to update threshold)
# despike: 0 good, 2 spike, 1 physical range out
# fetch: 70% fetch outside of region of interest
# To do: WD filter - need to know orientation for all sites. sh 

EC <- EC %>%
	mutate(LE_QC_SSITC = LE_SSITC_TEST, # Steady State and Integral Turbulence Characteristics test
				 LE_QC_sig_str = ifelse(H2O_sig_strgth_Min < 0.6, 1, 0), # signal strength of IRGA
				 LE_QC_despike = NA, # despike + physical range
				 LE_QC_fetch = NA, # fetch filtering 
				 LE_QC_ustar = NA,
				 LE_QC_longrun = NA,
				 # LE_QC_WD = NA, # data from backside of sonic
				 
				 H_QC_SSITC = H_SSITC_TEST, # Steady State and Integral Turbulence Characteristics test
				 H_QC_sig_str = 0, # signal strength of IRGA
				 H_QC_despike = NA, # despike + physical range
				 H_QC_fetch = NA, # fetch filtering
				 H_QC_ustar = NA,
				 H_QC_longrun = NA,
				 # H_QC_WD = NA, # data from backside of sonic
				 
				 FC_QC_SSITC = FC_SSITC_TEST, # Steady State and Integral Turbulence Characteristics test
				 FC_QC_sig_str = ifelse(CO2_sig_strgth_Min < 0.6, 1, 0), # signal strength of IRGA
				 FC_QC_despike = NA, # despike + physical range
				 FC_QC_fetch = NA, # fetch filtering
				 FC_QC_ustar = NA,
				 FC_QC_longrun = NA,
				 # FC_QC_WD = NA # data from backside of sonic
				 )


################################################################################
#### Despike and physical range, fetch filtering, and filter long-run
################################################################################
# physical range 
if(site_id %in% c('EDVG', 'ERVA')){
	# physical range for agriculture
	LE_range = c(-30, 1100)
	H_range = c(-300, 1100)
	FC_range = c(-60,30)
		
} else{
	# physical range for phreatophytes
	LE_range = c(-30, 600)
	H_range = c(-100, 1100)
	FC_range = c(-30,10)
}

# despike + physical range filtering
EC$LE_QC_despike <- pipeline_despike(as.data.frame(EC %>% mutate(timestamp = TIMESTAMP)),
					"LE",iter = 10,light = "PotRad", z = 10, var_thr = LE_range)
EC$H_QC_despike <- pipeline_despike(as.data.frame(EC %>% mutate(timestamp = TIMESTAMP)),
															"H",iter = 10,light = "PotRad", z = 10, var_thr = H_range)
EC$FC_QC_despike <- pipeline_despike(as.data.frame(EC %>% mutate(timestamp = TIMESTAMP)),
														 "FC",iter = 10,light = "PotRad", z = 10, var_thr = FC_range)

# filter out if 70% fetch is outside of defined region of interest
if(site_id %in% c('EDVG', 'ERVA')){
  # fetch filter applying only for agriculture
  roi <- read_csv(pipeline_path(base_dir, dirs_use$dir_fetch))
  EC$LE_QC_fetch <- fetch_filter(as.data.frame(EC %>% mutate(timestamp = TIMESTAMP)),
                                 "FETCH_70", "WD", roi$Fetch)
  EC$LE_QC_fetch[is.na(EC$LE_QC_fetch)] <- 0
  
} else{
  EC$LE_QC_fetch <- 0
}
EC$H_QC_fetch <- EC$LE_QC_fetch
EC$FC_QC_fetch <- EC$LE_QC_fetch

# filterLongRuns (replace runs, i.e sequences of numerically equal values, by NA)
EC_longrun <- EC %>%
	select(LE,H,FC) %>%
	filterLongRuns(c("LE","H","FC"), minNRunLength = 8)

EC$LE_QC_longrun <- ifelse(is.na(EC_longrun$LE) & !is.na(EC$LE),1,0)
EC$H_QC_longrun <- ifelse(is.na(EC_longrun$H) & !is.na(EC$H),1,0)
EC$FC_QC_longrun <- ifelse(is.na(EC_longrun$FC) & !is.na(EC$FC),1,0)

################################################################################
#### ReddyProc ustar threshold 
################################################################################
# Will use this function when working on FC data...
# EddyData <- EC %>%
# 	mutate(Year = year(TIMESTAMP),
# 				 DoY = yday(TIMESTAMP),
# 				 Hour = hour(TIMESTAMP) + minute(TIMESTAMP)/60,
# 				 NEE = ifelse(FC_QC_despike >= 1 | FC_QC_fetch >= 1,NA,FC),
# 				 LE = LE,
# 				 H = H, 
# 				 Rg = SW_IN_PI_F,
# 				 Tair = TA_PI_F,
# 				 rH = RH_PI_F,
# 				 VPD = VPD_PI_F,
# 				 Ustar = USTAR
# 				 ) %>%
# 	select(Year,DoY,Hour,NEE,LE,H,Rg,Tair,rH,VPD,Ustar) %>%
# 	as.data.frame()
# 
# #+++ Add time stamp in POSIX time format
# EddyDataWithPosix <- EddyData %>% 
# 	filterLongRuns("NEE") %>% 
# 	fConvertTimeToPosix('YDH', Year = 'Year', Day = 'DoY', Hour = 'Hour')
# #+++ Initalize R5 reference class sEddyProc for processing of eddy data
# #+++ with all variables needed for processing later
# EProc <- sEddyProc$new(
# 	site_id, EddyDataWithPosix, c('NEE','Rg','Tair','VPD', 'Ustar'))
# seasonFactor <- factor(rep("all", nrow(EddyDataWithPosix)))
# yearOfSeasonFactor <- c(all = unique(EddyData$Year)[1])
# 
# 
# # to do: need to investigate variable ustar and different threshold...
# set.seed(0815)
# uStarRes <- EProc$sEstUstarThresholdDistribution(seasonFactor = seasonFactor)
# print(uStarRes)
# threshold <- uStarRes$`50%`[1]
# # threshold <- min(uStarRes$uStarTh$uStar)
# 
# EC <- EC %>%
# 	mutate(LE_QC_ustar = ifelse(USTAR < threshold & SW_IN_PI_F < 30, 1, 0), 
# 				 H_QC_ustar = LE_QC_ustar, 
# 				 FC_QC_ustar = LE_QC_ustar
# 				 )

EC <- EC %>%
  mutate(LE_QC_ustar = 0,
         H_QC_ustar = 0,
         FC_QC_ustar = 0
  )

################################################################################
#### Visualization
################################################################################

plot_flux_qc_pdf <- function(df,
														 var = "LE",
														 site_id = "",
														 file_out = "",
														 n_panel_per_page = 6,
														 months_per_panel = 2) {
	
	stopifnot("TIMESTAMP" %in% names(df))
	stopifnot(var %in% names(df))
	
	#--------------------------------------------------
	# QC column names
	#--------------------------------------------------
	qc_ssitc   <- paste0(var, "_QC_SSITC")
	qc_sig_str <- paste0(var, "_QC_sig_str")
	qc_despike <- paste0(var, "_QC_despike")
	qc_fetch   <- paste0(var, "_QC_fetch")
	qc_ustar   <- paste0(var, "_QC_ustar")
	qc_longrun   <- paste0(var, "_QC_longrun")
	
	
	qc_cols <- c(qc_ssitc, qc_sig_str, qc_despike, qc_fetch,qc_ustar, qc_longrun)
	missing_qc <- qc_cols[!qc_cols %in% names(df)]
	
	if (length(missing_qc) > 0) {
		stop("Missing QC columns: ", paste(missing_qc, collapse = ", "))
	}
	
	#--------------------------------------------------
	# prepare
	#--------------------------------------------------
	df <- df %>%
		filter(!is.na(TIMESTAMP)) %>%
		arrange(TIMESTAMP)
	
	df_good <- df %>%
		filter(
			.data[[qc_ssitc]] <= 1,
			.data[[qc_sig_str]] == 0,
			.data[[qc_despike]] == 0,
			.data[[qc_fetch]] == 0,
			.data[[qc_ustar]] == 0,
			.data[[qc_longrun]] == 0,
		)
	
	if (nrow(df_good) == 0 || all(is.na(df_good[[var]]))) {
		# Missing or rejected observations are a valid diagnostic outcome.
		# Save an explicit status page instead of aborting the data publication.
		grDevices::pdf(file_out, width = 12, height = 7, useDingbats = FALSE)
		on.exit(grDevices::dev.off(), add = TRUE)
		graphics::plot.new()
		graphics::title(main = paste(site_id, var, "QAQC status"))
		graphics::text(0.5, 0.65, "No observations passed all QC tests in this period.")
		graphics::text(0.5, 0.5, paste("Rows:", nrow(df),
			"| Available flux:", sum(is.finite(df[[var]])),
			"| Missing despike flags:", sum(is.na(df[[qc_despike]]))))
		if (nrow(df)) graphics::text(0.5, 0.35,
			paste(format(min(df$TIMESTAMP)), "to", format(max(df$TIMESTAMP))))
		return(invisible(file_out))
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
		"LE" = "Latent heat flux (LE)",
		"H" = "Sensible heat flux (H)",
		"FC" = expression(CO[2]~flux),
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
			filter(TIMESTAMP >= s, TIMESTAMP < e)
		
		if (nrow(df_sub) == 0) next
		
		df_sub_good <- df_sub %>%
			filter(
				.data[[qc_ssitc]] <= 1,
				.data[[qc_sig_str]] == 0,
				.data[[qc_despike]] == 0,
				.data[[qc_fetch]] == 0,
				.data[[qc_ustar]] == 0,
				.data[[qc_longrun]] == 0,
			)
		
		p <- ggplot(df_sub_good, aes(x = TIMESTAMP, y = .data[[var]])) +
			theme_bw() +
			geom_line(na.rm = TRUE) +
			geom_point(size = 1, na.rm = TRUE) +
			geom_point(
				aes(y = .data[[var]], color = "SSITC"),
				data = df_sub %>% filter(.data[[qc_ssitc]] >= 2),
				size = 1,
				na.rm = TRUE
			) +
			geom_point(
				aes(y = .data[[var]], color = "signal_str"),
				data = df_sub %>% filter(.data[[qc_sig_str]] >= 1),
				size = 1,
				na.rm = TRUE
			) +
			geom_point(
				aes(y = .data[[var]], color = "despike & physical range"),
				data = df_sub %>% filter(.data[[qc_despike]] >= 1),
				size = 1,
				na.rm = TRUE
			) +
			geom_point(
				aes(y = .data[[var]], color = "fetch70%"),
				data = df_sub %>% filter(.data[[qc_fetch]] >= 1),
				size = 1,
				na.rm = TRUE
			) +
			# geom_point(
			# 	aes(y = .data[[var]], color = "ustar"),
			# 	data = df_sub %>% filter(.data[[qc_ustar]] >= 1),
			# 	size = 1,
			# 	na.rm = TRUE
			# ) +
			geom_point(
				aes(y = .data[[var]], color = "longrun"),
				data = df_sub %>% filter(.data[[qc_longrun]] >= 1),
				size = 1,
				na.rm = TRUE
			) +
			coord_cartesian(ylim = y_rng) +
			scale_color_manual(
				values = c(
					"fetch70%" = "dodgerblue3",
					"despike & physical range" = "orange2",
					"signal_str" = "purple3",
					"SSITC" = "red3",
					# "ustar" = "green",
					"longrun" = "gray"
				),
				drop = FALSE
			) +
			labs(
				title = paste0(
					site_id, " | ", var, " | ",
					format(s, "%Y-%m"), " to ",
					format(e - seconds(1), "%Y-%m")
				),
				x = NULL,
				y = y_lab,
				caption = paste(
					"Black points indicate data passing all QC filters.",
					"Colored points represent data flagged by individual QC tests:",
					"fetch (<70%), despiking/physical range, signal strength,",
					"and SSITC quality control."
				)
			) +
			theme(
				plot.title = element_text(face = "bold", size = 10),
				plot.caption = element_text(size = 7, hjust = 0),
				axis.text.x = element_text(size = 7),
				axis.text.y = element_text(size = 7),
				axis.title.y = element_text(size = 8),
				panel.grid.minor = element_blank(),
				plot.margin = margin(3, 6, 3, 3)
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
	
	pdf(file_out, width = 16, height = 20, onefile = TRUE, useDingbats = FALSE)
	
	for (i in seq_along(pages)) {
		
		plotlist_i <- pages[[i]]
		
		if (length(plotlist_i) < n_panel_per_page) {
			plotlist_i <- c(
				plotlist_i,
				rep(list(NULL), n_panel_per_page - length(plotlist_i))
			)
		}
		
		page_plot <- plot_grid(
			plotlist = plotlist_i,
			nrow = n_panel_per_page,
			align = "v"
		)
		
		print(page_plot)
	}
	
	dev.off()
	
	message("Saved: ", file_out)
}


plot_flux_qc_pdf(pipeline_figure_data(EC), var = "LE",
								 file_out = paste0(pipeline_path(base_dir, dirs_use$dir_output),"/figures/",site_id,"_Level_3_LE_flag.pdf")
								 	)
plot_flux_qc_pdf(pipeline_figure_data(EC), var = "H",
								 file_out = paste0(pipeline_path(base_dir, dirs_use$dir_output),"/figures/",site_id,"_Level_3_H_flag.pdf")
								 )
plot_flux_qc_pdf(pipeline_figure_data(EC), var = "FC",
								 file_out = paste0(pipeline_path(base_dir, dirs_use$dir_output),"/figures/",site_id,"_Level_3_FC_flag.pdf")
								 )


################################################################################
#### save
################################################################################

write_custom_csv <- pipeline_write

EC_unit_final <- cbind(EC_unit,
                 data.frame(LE_QC_SSITC = c("#","derived"),LE_QC_sig_str = c("#","derived"),
											LE_QC_despike = c("#","derived"),LE_QC_fetch = c("#","derived"), LE_QC_ustar = c("#","derived"), LE_QC_longrun = c("#","derived"),
											H_QC_SSITC = c("#","derived"),H_QC_sig_str = c("#","derived"), H_QC_despike = c("#","derived"),
											H_QC_fetch = c("#","derived"),H_QC_ustar = c("#","derived"), H_QC_longrun = c("#","derived"), FC_QC_SSITC = c("#","derived"),
											FC_QC_sig_str = c("#","derived"),FC_QC_despike = c("#","derived"),
											FC_QC_fetch = c("#","derived"), FC_QC_ustar = c("#","derived"), FC_QC_longrun = c("#","derived")
											))

print(which(colnames(EC) != colnames(EC_unit_final)))

write_custom_csv(
	df = EC %>% select( -c(SW_IN_PI_F,TA_PI_F,RH_PI_F,VPD_PI_F,CO2_sig_strgth_Min,H2O_sig_strgth_Min,PotRad)),
	units_df = EC_unit_final %>% select( -c(SW_IN_PI_F,TA_PI_F,RH_PI_F,VPD_PI_F,CO2_sig_strgth_Min,H2O_sig_strgth_Min,PotRad)),
	file = paste0(pipeline_path(base_dir, dirs_use$dir_output),"/",site_id,"_Level_3_EC_with_qaqc_",date(Sys.time()),".csv"),
	description = "flux data and QAQC flag"
)


# loop all sites end
}
