# Regression for auxiliary audit columns added after the initial gap table.
# Evaluate the production assignment, including reordered and single columns.
script <- if (file.exists('R/Level_1_met_QAQC_new.R')) 'R/Level_1_met_QAQC_new.R' else 'Level_1_met_QAQC_new.R'
code <- readLines(script)
statement <- grep('^df_stat\\$gap_after_QC <-', code, value=TRUE)
stopifnot(length(statement)==1L)
update <- parse(text=statement)
df <- as.data.frame(setNames(rep(list(c(1,NA,2)),79),paste0('v',1:79)))
df_stat <- data.frame(vars=names(df),gap=1L)
for(i in 1:8) df[[paste0('audit',i)]] <- NA_real_
df$TIMESTAMP <- 1:3;df$PotRad<-0;df$ClearSkyRad<-0
df$v2 <- c(NA,NA,2)
df <- df[,rev(names(df))]
eval(update)
stopifnot(nrow(df_stat)==79,df_stat$gap_after_QC[2]==2,
 all(df_stat$gap_after_QC[-2]==1))
df_stat<-data.frame(vars='v2');eval(update)
stopifnot(df_stat$gap_after_QC==2)
df_stat<-data.frame(vars='missing')
stopifnot(inherits(tryCatch(eval(update),error=identity),'error'))
invisible(parse(script))
cat('PASS: 79 original variables plus 8 audit columns, reordered columns, single variable, missing-variable failure and stage syntax.\n')
