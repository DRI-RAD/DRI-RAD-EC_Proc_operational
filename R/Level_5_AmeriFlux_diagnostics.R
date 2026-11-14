# Diagnostics are built from the written submission CSV, not a different L4 view.
# Base R graphics + local HTML/PNG assets require no browser/network packages.
amf_html <- function(x) {
  x <- as.character(x); x[is.na(x)] <- ""
  for (pair in list(c("&","&amp;"),c("<","&lt;"),c(">","&gt;"),c('"',"&quot;"),c("'","&#39;")))
    x <- gsub(pair[1],pair[2],x,fixed=TRUE)
  x
}

amf_draw_variable <- function(t, y, variable, unit, path) {
  # Cairo avoids slow Windows GDI drawing for dense multi-year point clouds.
  device_type <- if (capabilities('cairo')) 'cairo' else getOption('bitmapType','windows')
  grDevices::png(path, width=1800, height=1050, res=150, type=device_type)
  on.exit(grDevices::dev.off(),add=TRUE)
  graphics::layout(matrix(c(1,1,2,3),2,2,byrow=TRUE),heights=c(1.35,1))
  graphics::par(mar=c(4,4.5,3,1),mgp=c(2.6,.8,0))
  valid <- is.finite(y); n <- sum(valid)
  xr <- range(t); if(diff(as.numeric(xr))==0)xr <- xr+c(-900,900)
  yr <- if(n)range(y[valid]) else c(0,1)
  if(diff(yr)==0)yr <- yr+c(-1,1)*max(abs(yr[1])*.03,1)
  graphics::plot(t,rep(NA_real_,length(t)),ylim=yr,xlim=xr,
    xlab="Local standard time (interval end)",ylab=unit,
    main=paste(variable, "| submitted half-hour values"))
  if(n)graphics::points(t[valid],y[valid],pch=20,cex=.24,col=grDevices::adjustcolor('#176D8C',alpha.f=.4))
  else graphics::text(mean(as.numeric(xr)),mean(yr),"All values missing",col='#8C3C29')
  graphics::mtext(sprintf("Valid: %s / %s | Missing: %.2f%% | No thinning, clipping or gap filling",format(n,big.mark=','),format(length(y),big.mark=','),100*mean(!valid)),side=3,line=.2,cex=.7)
  graphics::grid(col='#E5E9ED')
  discrete <- grepl('QC|SSITC|METHOD|DIAG',variable) && length(unique(y[valid]))<=32
  if(!n){ graphics::plot.new();graphics::title('Distribution');graphics::text(.5,.5,'No finite values') }
  else if(discrete){
    tab <- table(y[valid]);graphics::barplot(tab,col='#176D8C',border=NA,
      main='Code frequencies (no assumed good/bad meaning)',xlab='Exported code',ylab='Count',cex.names=.75)
  } else {
    graphics::hist(y[valid],breaks=40,col='#176D8C',border='white',main='Distribution',xlab=unit)
  }
  # Attribute interval-end midnight to the month containing its interval start.
  month <- format(t-1800,'%Y-%m',tz='UTC')
  months <- unique(month)
  completeness <- vapply(months,function(m)100*mean(valid[month==m]),numeric(1))
  bp <- graphics::barplot(completeness,ylim=c(0,100),col='#578D70',border=NA,
    main='Monthly data availability',ylab='Valid exported intervals (%)',names.arg=rep('',length(months)))
  shown <- unique(round(seq(1,length(months),length.out=min(8,length(months)))))
  graphics::axis(1,at=bp[shown],labels=months[shown],cex.axis=.7,las=2)
  graphics::mtext('Denominator: intervals present in this export, including partial months',side=3,line=.1,cex=.65)
}

plot_level5_ameriflux <- function(file, output_dir=NULL) {
  file <- normalizePath(file,winslash='/',mustWork=TRUE)
  if(is.null(output_dir))output_dir <- file.path(dirname(file),'diagnostics')
  output_dir <- normalizePath(output_dir,winslash='/',mustWork=FALSE)
  before <- tools::md5sum(file)
  d <- readr::read_csv(file,show_col_types=FALSE,progress=FALSE,na=c('-9999','NA',''),
    col_types=readr::cols(.default=readr::col_double(),
      TIMESTAMP_START=readr::col_character(),TIMESTAMP_END=readr::col_character()))
  if(nrow(readr::problems(d)))stop('Nonnumeric/invalid input in submission CSV; inspect readr parsing problems.')
  if(!all(c('TIMESTAMP_START','TIMESTAMP_END')%in%names(d)) || !nrow(d))stop('Expected a nonempty Level 5 CSV with interval timestamps.')
  t <- as.POSIXct(d$TIMESTAMP_END,format='%Y%m%d%H%M',tz='UTC')
  ts <- as.POSIXct(d$TIMESTAMP_START,format='%Y%m%d%H%M',tz='UTC')
  if(anyNA(t)||anyNA(ts)||anyDuplicated(t)||is.unsorted(t)||any(as.numeric(t-ts,units='secs')!=1800))
    stop('Invalid half-hour interval timestamps in submission CSV.')
  vars <- setdiff(names(d),c('TIMESTAMP_START','TIMESTAMP_END'))
  if(!length(vars))stop('No submitted data variables to plot.')
  mapping_file <- file.path(dirname(file),'metadata','variable_mapping.csv')
  mapping <- if(file.exists(mapping_file))read.csv(mapping_file,stringsAsFactors=FALSE,check.names=FALSE) else data.frame()
  if(nrow(mapping))mapping <- mapping[mapping$status=='exported',,drop=FALSE]
  if(nrow(mapping)&&anyDuplicated(mapping$target))stop('Duplicate exported mapping targets.')
  dir.create(file.path(output_dir,'figures'),recursive=TRUE,showWarnings=FALSE)
  results <- vector('list',length(vars));cards <- character(length(vars))
  fmt <- function(v)ifelse(is.finite(v),format(signif(v,5),trim=TRUE),'NA')
  for(i in seq_along(vars)){
    v <- vars[i];y <- d[[v]]; valid <- is.finite(y)
    idx <- if(nrow(mapping))match(v,mapping$target) else NA_integer_
    meta <- function(col,default='')if(!is.na(idx)&&col%in%names(mapping))mapping[[col]][idx] else default
    unit <- meta('output_unit','Unit not supplied');source <- meta('source');note <- meta('note');stage <- meta('source_stage')
    relative <- sprintf('figures/%03d_%s.png',i,gsub('[^A-Za-z0-9_-]','_',v))
    problem <- tryCatch({amf_draw_variable(t,y,v,unit,file.path(output_dir,relative));''},error=function(e)conditionMessage(e))
    row <- data.frame(variable=v,unit=unit,source=source,n=nrow(d),valid=sum(valid),missing=sum(!valid),
      missing_percent=100*mean(!valid),min=if(any(valid))min(y[valid]) else NA_real_,
      max=if(any(valid))max(y[valid]) else NA_real_,mean=if(any(valid))mean(y[valid]) else NA_real_,
      median=if(any(valid))median(y[valid]) else NA_real_,
      first_valid=if(any(valid))format(min(t[valid]),'%Y-%m-%d %H:%M',tz='UTC') else '',
      last_valid=if(any(valid))format(max(t[valid]),'%Y-%m-%d %H:%M',tz='UTC') else '',
      plot=if(nzchar(problem))'' else relative,status=if(nzchar(problem))'error' else 'ok',error=problem)
    results[[i]]<-row
    cards[i]<-paste0('<details class="variable" data-search="',amf_html(tolower(paste(v,source,unit,note))),
      '" data-missing="',if(all(!valid))'all' else if(any(!valid))'some' else 'none','"><summary><strong>',amf_html(v),
      '</strong> <span>',amf_html(unit),' | missing ',sprintf('%.2f',row$missing_percent),'% | range ',fmt(row$min),' to ',fmt(row$max),
      if(nzchar(problem))' | PLOT ERROR' else '', '</span></summary><p>Source: ',amf_html(source),' | ',amf_html(stage),
      '</p><p>',amf_html(note),'</p>',if(nzchar(problem))paste0('<p class="error">',amf_html(problem),'</p>') else
      paste0('<a href="',relative,'" target="_blank"><img loading="lazy" src="',relative,'" alt="',amf_html(v),' diagnostic plots"></a>'),'</details>')
    if(i%%25==0||i==length(vars))message('AmeriFlux diagnostics: ',i,'/',length(vars),' variables')
  }
  summary<-do.call(rbind,results)
  write.csv(summary,file.path(output_dir,'variable_summary.csv'),row.names=FALSE,na='')
  table_rows<-vapply(seq_len(nrow(summary)),function(i){r<-summary[i,];paste0('<tr><td>',amf_html(r$variable),'</td><td>',amf_html(r$unit),'</td><td>',r$valid,'</td><td>',sprintf('%.2f',r$missing_percent),'</td><td>',fmt(r$min),'</td><td>',fmt(r$max),'</td><td>',fmt(r$mean),'</td></tr>')},character(1))
  html <- c('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>AmeriFlux submission diagnostics</title>',
    '<style>body{font:16px system-ui,sans-serif;color:#203447;background:#f6f8fa;margin:0}main{max-width:1200px;margin:auto;padding:28px}h1{font-size:28px}header,p{line-height:1.5}.controls{position:sticky;top:0;background:#e6edf3;padding:14px;z-index:1;display:flex;gap:10px;flex-wrap:wrap}input,select,button{font:inherit;padding:8px;border:1px solid #a9bac7;border-radius:5px}input{flex:1;min-width:220px}details{background:white;border:1px solid #d3dfe7;margin:9px 0;border-radius:6px;padding:12px}summary{cursor:pointer}summary span{font-size:14px;color:#516475}img{width:100%;height:auto}table{width:100%;border-collapse:collapse;font-size:13px}th,td{text-align:left;padding:7px;border-bottom:1px solid #ddd}.scroll{overflow:auto;max-height:400px}.error{color:#9b3420}[hidden]{display:none!important}a{color:#126b91}</style></head><body><main>',
    paste0('<header><h1>AmeriFlux submission diagnostics</h1><p>',amf_html(basename(file)),'<br>',format(min(ts),'%Y-%m-%d %H:%M',tz='UTC'),' through ',format(max(t),'%Y-%m-%d %H:%M',tz='UTC'),' | local standard time</p>'),
    paste0('<p>',nrow(d),' half-hour intervals | ',length(vars),' submitted variables | ',sum(summary$valid==0),' all-missing variables | ',sum(summary$status=='error'),' plot errors.</p>'),
    '<p>Plots use the actual exported values; -9999 is missing. No thinning, outlier removal, filling or QC is applied here. Timestamp columns define the horizontal axis. Click a variable to expand, then click its image for full resolution.</p>',
    paste0('<p>Source CSV MD5: <code>',unname(before),'</code> | <a href="variable_summary.csv">Download complete summary CSV</a></p></header>'),
    '<details><summary>All-variable numerical summary</summary><div class="scroll"><table><thead><tr><th>Variable</th><th>Unit</th><th>Valid</th><th>Missing %</th><th>Min</th><th>Max</th><th>Mean</th></tr></thead><tbody>',table_rows,'</tbody></table></div></details>',
    '<div class="controls"><input id="search" type="search" placeholder="Search variable, source, unit or processing"><select id="missing"><option value="any">All variables</option><option value="all">All missing</option><option value="some">Partially missing</option><option value="none">Complete</option></select><button id="expand">Expand visible</button><button id="collapse">Collapse all</button><span id="count" aria-live="polite"></span></div>',cards,
    '<script>const cards=[...document.querySelectorAll("details.variable")];const search=document.getElementById("search"),missing=document.getElementById("missing");function filter(){const q=search.value.trim().toLowerCase();let n=0;cards.forEach(c=>{c.hidden=!(c.dataset.search.includes(q)&&(missing.value==="any"||c.dataset.missing===missing.value));if(!c.hidden)n++;});document.getElementById("count").textContent=n+" / "+cards.length;}search.addEventListener("input",filter);missing.addEventListener("change",filter);document.getElementById("expand").onclick=()=>cards.forEach(c=>{if(!c.hidden)c.open=true});document.getElementById("collapse").onclick=()=>cards.forEach(c=>c.open=false);filter();</script></main></body></html>')
  connection<-file(file.path(output_dir,'index.html'),open='w',encoding='UTF-8')
  tryCatch(writeLines(html,connection,useBytes=FALSE),finally=close(connection))
  stopifnot(identical(before,tools::md5sum(file)))
  if(any(summary$status=='error'))warning('Some variable plots failed; see diagnostics/variable_summary.csv. Submission CSV is unchanged.',call.=FALSE)
  invisible(file.path(output_dir,'index.html'))
}
