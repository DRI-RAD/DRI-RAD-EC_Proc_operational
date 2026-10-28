# Synthetic tests: export semantics, unit conversions, timestamp and file format.
source("run_level5_ameriflux.R", encoding = "UTF-8")
fail <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
t <- pipeline_time(c("2024-02-28 23:30:00", "2024-02-29 00:00:00", "2024-02-29 01:00:00"))
d <- data.frame(TIMESTAMP = t, processing = c("EddyPro", "EasyFlux", "EddyPro"),
  LE = c(100,200,300), H = c(10,20,30), FC = c(-1,-2,-3),
  PA = c(90000,91000,92000), PA_PI_F = c(90,91,92),
  TA_1_1_3 = c(280,281,282), RH_1_1_3 = c(0.2,0.3,0.4),
  SWC_1_1_1 = c(0.1,0.2,0.3), T_SONIC = c(300,25,301), T_SONIC_SIGMA = c(2,3,4),
  G_plate = 10, SG = 3, G_PI = 13, G_PI_F = 14, NETRAD = 100, NETRAD_PI_F = 101,
  LE_PI_F = c(110,220,330), H_PI_F = 40, FC_PI_F = 5,
  SLE = 5, SH = 6, SC = 0.2, VPD_PI_F = c(1,2,3), ALB = c(20,30,40),
  CO2 = 400, H2O = 10, USTAR = c(0.001,0.002,0.003), WS = 3,
  TS_1_1_1 = 5, PPFD_IN = 20, P = c(1,NA,3),
  P_Bulk = c(2,NA,4), P_PI_F = c(2,0,4),
  LE_710_filtered = c(90, NA, 270), H_710_filtered = c(9, NA, 27),
  LE_710_PI_F = c(90, 180, 270), H_710_PI_F = c(9, 18, 27),
  SW_IN = c(100, NA, 300), SW_IN_PI_F = c(100, 200, 300))
for (v in c("LE", "H", "FC")) {
  d[[paste0(v,"_QC_despike")]] <- c(0,1,NA_real_)
  d[[paste0(v,"_QC_fetch")]] <- c(0,0,1)
  # Deliberately bad SSITC, signal strength, long run and USTAR MUST NOT mask row1.
  for (s in c("_QC_SSITC", "_SSITC_TEST", "_QC_sig_str", "_QC_longrun", "_QC_ustar"))
    d[[paste0(v,s)]] <- c(2,0,0)
  d[[paste0(v,"_filtered1")]] <- NA_real_
}
units <- setNames(rep("#", ncol(d)), names(d))
units[c("LE","H","SLE","SH","G_plate","SG","G_PI","G_PI_F","NETRAD","NETRAD_PI_F","LE_PI_F","H_PI_F")] <- "[W+1m-2]"
units[c("FC","SC","FC_PI_F")] <- "[\u00b5mol+1s-1m-2]"
units[c("PA")] <- "Pa"; units['PA_PI_F'] <- 'kPa'
units[c("TA_1_1_3","T_SONIC_SIGMA")] <- "K"
units['T_SONIC'] <- "deg C"; units['RH_1_1_3'] <- 'fraction'
units['SWC_1_1_1'] <- 'm3/m3'; units['ALB'] <- 'fraction' # legacy header, L1 computes percent
units['VPD_PI_F'] <- 'kPa'; units['TS_1_1_1'] <- 'unknown'
units['CO2'] <- 'ppm'; units['H2O'] <- 'mmol/mol'; units[c('WS','USTAR')] <- 'm/s'
units['PPFD_IN'] <- 'umol/m2/s';units['P'] <- 'mm'
units[c('P_Bulk','P_PI_F')] <- 'mm'
units[c('LE_710_filtered','H_710_filtered','LE_710_PI_F','H_710_PI_F','SW_IN','SW_IN_PI_F')] <- 'W m-2'
u <- as.data.frame(as.list(units), stringsAsFactors = FALSE)
item <- list(data=d, units=u)
p <- amf_prepare(item)
x <- p$data
stopifnot(nrow(x)==4, x$TIMESTAMP_START[2]=='202402282330', x$TIMESTAMP_END[3]=='202402290030',
  x$LE_1_1_1[1]==100, x$H_1_1_1[1]==10, x$FC[1]==-1, x$LE_SSITC_TEST_1_1_1[1]==2,
  all(is.na(x$LE_1_1_1[2:4])), x$SLE_1_1_1[1]==5,
  x$G[1]==10,x$SG[1]==3,
  x$PA[1]==90,abs(x$TA_1_1_3[1]-6.85)<1e-9,
  x$RH_1_1_3[1]==20,x$SWC_1_1_1[1]==10,x$ALB[1]==20,
  abs(x$T_SONIC[1]-26.85)<1e-9,x$T_SONIC[2]==25,x$T_SONIC_SIGMA[1]==2,
  !'TS_1_1_1'%in%names(x),
  setequal(grep('_PI', names(x), value=TRUE), c('LE_PI_F','H_PI_F','FC_PI_F')),
  x$LE_PI_F[2]==220, x$H_PI_F[2]==40, x$FC_PI_F[2]==5,
  is.na(x$LE_PI_F[3]), !any(c('LE','H') %in% names(x)),
  x$LE_1_1_2[1]==90, x$H_1_1_2[4]==27,
  x$LE_F_1_1_2[2]==180, x$H_F_1_1_2[2]==18,
  x$P_1_1_1[1]==2, x$P_1_1_2[1]==1, is.na(x$P_1_1_1[2]),
  is.na(x$P_1_1_2[2]), x$P_F_1_1_1[2]==0,
  !'P' %in% names(x), !'FC_F_1_1_1' %in% names(x), !'LE_F_1_1_1' %in% names(x),
  is.na(x$LE_1_1_2[2]), is.na(x$SW_IN[2]),
  any(p$audit$source=='TS_1_1_1' & p$audit$status=='omitted'))
review <- amf_prepare(item, include_pi=TRUE)$data
stopifnot(review$LE_PI_F[2]==220, review$G_PI[1]==13, review$VPD_PI_F[1]==10)
stopifnot(!anyDuplicated(names(review)),
  identical(amf_prepare(item, include_storage=FALSE)$data$LE_PI_F, x$LE_PI_F))
no_ec_filled <- item
no_ec_filled$data[c('LE_PI_F','H_PI_F','FC_PI_F')] <- NULL
stopifnot(!any(c('LE_PI_F','H_PI_F','FC_PI_F') %in% names(amf_prepare(no_ec_filled)$data)))
absent <- item
absent$data$LE_710_filtered <- absent$data$H_710_filtered <- NULL
stopifnot(!any(c('LE_1_1_2','H_1_1_2') %in% names(amf_prepare(absent)$data)))
absent$data$LE_710_PI_F <- absent$data$H_710_PI_F <- absent$data$P_Bulk <- NULL
stopifnot(!any(c('LE_F_1_1_2','H_F_1_1_2','P_1_1_1') %in% names(amf_prepare(absent)$data)),
          'P_1_1_2' %in% names(amf_prepare(absent)$data))
fail(amf_prepare(item, include_pi=NA))
p2 <- amf_prepare(item, unit_overrides=c(TS_1_1_1='deg C'))
stopifnot(p2$data$TS_1_1_1[1]==5)
fail(amf_prepare(item,local_standard_time=FALSE))
bad<-item;bad$data$TIMESTAMP[2]<-bad$data$TIMESTAMP[1];fail(amf_prepare(bad))
bad<-item;bad$data$TIMESTAMP[1]<-bad$data$TIMESTAMP[1]+60;fail(amf_prepare(bad))
bad<-item;bad$data$LE_QC_fetch<-NULL;fail(amf_prepare(bad))
bad<-item;bad$units$FC<-'wrong';fail(amf_prepare(bad))
fail(amf_prepare(item,start='2025-01-01 00:00:00'))
stopifnot(is.na(amf_convert(Inf,'Pa','kPa','PA')),is.na(amf_convert(-9999,'Pa','kPa','PA')))
dir <- tempfile('amf-test-');dir.create(dir)
file <- amf_write(p,'ECSM',dir)
lines<-readLines(file)
stopifnot(length(lines)==5,!any(grepl('"',lines)),any(grepl('-9999',lines)),
  basename(file)=='ECSM_Level_5_ameriflux_ready.csv',file.exists(file.path(dirname(file),'metadata','variable_mapping.csv')))
# End-to-end runner with real pipeline CSV parsing and input preservation.
site_dir<-file.path(dir,'site');dir.create(site_dir)
input<-file.path(site_dir,'ECSM_Level_4_post_processed_data_test.csv')
pipeline_write_table(d,rbind(u,u),input,'test')
hash<-tools::md5sum(input)
config<-file.path(dir,'config.csv');write.csv(data.frame(site='ECSM',dir_output=site_dir),config,row.names=FALSE)
preview<-run_level5_ameriflux('ECSM',base_dir=dir,config_file=config,dry_run=TRUE)
stopifnot(!dir.exists(file.path(site_dir,'AmeriFlux')),nrow(preview$ECSM$data)==4)
written<-run_level5_ameriflux('ECSM',base_dir=dir,config_file=config)
stopifnot(file.exists(written$ECSM),identical(hash,tools::md5sum(input)))
# Automatic export still runs when L4 has no pending observations. It exports
# the complete cumulative period, not just the bounds of this pipeline call.
source("run_pipeline.R")
for (stage in c("L1", "L2", "L3_EC")) {
  pipeline_write_table(d, rbind(u,u),
    file.path(site_dir, paste0("ECSM_", pipeline_patterns[[stage]], "_test.csv")), "test")
}
write.csv(data.frame(site='ECSM',dir_output=site_dir,dir_eddypro=NA_character_), config, row.names=FALSE)
args <- list(sites='ECSM', stages='L4', base_dir=dir, config_file=config,
             start=t[2], end=t[3])
count_exports <- function() length(list.files(file.path(site_dir, 'AmeriFlux'),
  pattern='ameriflux_ready.csv$', recursive=TRUE))
before <- count_exports()
preview <- do.call(run_pipeline, c(args, list(dry_run=TRUE)))
stopifnot(count_exports()==before, attr(preview,'ameriflux')$ECSM=='planned')
disabled <- do.call(run_pipeline, c(args, list(export_ameriflux=FALSE)))
stopifnot(count_exports()==before, length(attr(disabled,'ameriflux'))==0)
automatic <- do.call(run_pipeline, args)
export_file <- attr(automatic, 'ameriflux')$ECSM
exported <- read.csv(export_file, colClasses='character')
stopifnot(automatic[['ECSM:L4']]=='skipped', count_exports()==before+1,
  nrow(exported)==4, exported$TIMESTAMP_END[1]=='202402282330',
  identical(hash, tools::md5sum(input)))
again <- do.call(run_pipeline, args)
stopifnot(count_exports()==before+2, file.exists(attr(again,'ameriflux')$ECSM))
fail(do.call(run_pipeline, c(args,list(export_ameriflux=NA))))
unlink(dir,recursive=TRUE)
cat('PASS: requested flux masks, separate storage/PI products, unit provenance, missing intervals/leap day, errors, CSV writing, runner and immutable inputs.\n')
