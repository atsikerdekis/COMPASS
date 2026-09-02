############
### INIT ###
############
source(paste0(path_code,"04.preprocess.R"))

plot_title <- "Aerosol"
plot_level <- 1000

figure_box <- FALSE
field_show_box <- FALSE
lonmin <- 0
lonmax <- 360
latmin <- -90
latmax <- 90
projection <- "+proj=robin"
gridlines <- 10
coastlineWorldFine_lwd <- 1

seqDate <- format(seq.Date(from=as.Date(sDate,format="%Y%m%d"),to=as.Date(eDate,format="%Y%m%d"),by="day"),"%Y%m%d")

########################
### VARIABLE FAMILIES ###
########################
plot_types <- c("mmr","ddp","sdm","wdl","wdc","mss","mss_from_mr","ngt")

plot_type_title <- c(
  mmr = "Mass mixing ratio",
  ddp = "Dry deposition",
  sdm = "Sedimentation",
  wdl = "Large-scale wet dep.",
  wdc = "Convective wet dep.",
  mss = "Column mass burden",
  mss_from_mr = "Column burden from mixing ratio",
  ngt = "Negative fixer"
)

dep_fluxes <- c("ddp","sdm","wdl","wdc","ngt")
dep_flux_colors <- c(ddp="#846040",sdm="#D78C6A",wdl="#74ACE8",wdc="#3E7DD1",ngt="#3FC13A")
requested_individual_variables <- variables_requested[!startsWith(variables_requested,"dep_")]

########################
### HELPER FUNCTIONS ###
########################
axis_ticks <- function(x,n=5,small_thresh=0.01,large_thresh=100) {
  rng <- range(x,na.rm=TRUE,finite=TRUE)
  ticks <- pretty(rng,n=n)
  max_abs <- max(abs(ticks))
  min_abs <- if (all(ticks == 0)) 0 else min(abs(ticks[ticks != 0]))
  if (max_abs >= large_thresh || (min_abs > 0 && min_abs < small_thresh)) labels <- sprintf("%.2e",ticks) else {
    step <- min(diff(ticks))
    ndigits <- if (step > 0) max(2,-floor(log10(step))) else 2
    labels <- formatC(ticks,format="f",digits=ndigits)
  }
  list(breaks=ticks,labels=labels)
}

get_type_variables <- function(variable_table,type) {
  if (type == "mss") {
    x <- unique(variable_table$logical_name[
      startsWith(variable_table$logical_name,"mss_") &
      !startsWith(variable_table$logical_name,"mss_from_mr_")
    ])
  } else {
    x <- unique(variable_table$logical_name[
      startsWith(variable_table$logical_name,paste0(type,"_"))
    ])
  }
  x[x %in% requested_individual_variables]
}

nice_positive_max <- function(x,n=6) {
  xmax <- max(x,na.rm=TRUE)
  if (!is.finite(xmax) || xmax <= 0) return(1)
  max(pretty(c(0,xmax),n=n))
}

positive_breaks <- function(x,ncolors=200) {
  seq(0,nice_positive_max(x),length.out=ncolors+1)
}

difference_breaks <- function(x,ncolors=200) {
  xmax <- max(abs(x),na.rm=TRUE)
  if (!is.finite(xmax) || xmax == 0) xmax <- 1
  xmax <- max(abs(pretty(c(-xmax,xmax),n=6)))

  breaks <- seq(-xmax,xmax,length.out=ncolors+1)
  breaks[which.min(abs(breaks))] <- 0

  breaks
}

positive_axis_ticks <- function(x,n=10) {
  xmax <- nice_positive_max(x,n)
  ticks <- pretty(c(0,xmax),n=n)
  ticks <- ticks[ticks >= 0 & ticks <= xmax]
  ticks <- unique(c(0,ticks,xmax))
  max_abs <- max(abs(ticks))
  min_abs <- if (all(ticks == 0)) 0 else min(abs(ticks[ticks != 0]))
  if (max_abs >= 100 || (min_abs > 0 && min_abs < 0.01)) labels <- sprintf("%.2e",ticks) else {
    step <- min(diff(ticks))
    ndigits <- if (step > 0) max(2,-floor(log10(step))) else 2
    labels <- formatC(ticks,format="f",digits=ndigits)
  }
  list(breaks=ticks,labels=labels)
}

read_lon_lat <- function(file) {
  nc <- nc_open(file)
  lon <- ncvar_get(nc,"longitude")
  lat <- ncvar_get(nc,"latitude")
  nc_close(nc)
  list(lon=lon,lat=lat)
}

get_plot_units <- function(type) {
  if (type == "mmr") return("kg kg^-1")
  if (type %in% c("mss","mss_from_mr")) return("kg m^-2")
  if (type %in% c("ddp","sdm","wdl","wdc","ngt")) return("kg m^-2 s^-1")
  " "
}

get_plot_category <- function(logical_name) {

  suffix <- get_variable_suffix(logical_name)

  ### HAM vs AER: only common species comparisons are allowed
  if (exptype1 != exptype2) return("per_species")

  ### AER vs AER
  if (exptype1 == "AER" && exptype2 == "AER") {
    if (suffix %in% common_HAM_AER_species) return("per_species")
    if (suffix %in% aer_names) return("per_tracer")
    stop("Could not determine AER plot category for ",logical_name)
  }

  ### HAM vs HAM
  if (suffix %in% c("soluble","insoluble")) return("per_solubility")
  if (suffix %in% ham_modes) return("per_mode")
  if (grepl(paste0("_(",paste(ham_modes,collapse="|"),")$"),suffix)) return("per_tracer")
  return("per_species")
}

read_daily_variable <- function(expname,exptype,logical_name,variable_table,type,date) {
  file <- variable_file(logical_name,variable_table,expname,date)
  if (!file.exists(file)) stop("Input file not found: ",file)
  if (type == "mmr") return(read_variable_level(file,exptype,logical_name,variable_table,plot_level))
  if (type == "mss_from_mr") { return(read_column_burden(file,exptype,logical_name,variable_table)) }
  read_variable(file,exptype,logical_name,variable_table)
}

read_plot_data <- function(expname,exptype,logical_name,variable_table,type) {
  field_day <- list()
  for (d in seq_along(seqDate)) field_day[[d]] <- read_daily_variable(expname,exptype,logical_name,variable_table,type,seqDate[d])
  nx <- dim(field_day[[1]])[1]
  ny <- dim(field_day[[1]])[2]
  nt <- sum(sapply(field_day,function(x) dim(x)[3]))
  data <- array(unlist(field_day),dim=c(nx,ny,nt))
  list(nx=nx,ny=ny,nt=nt,data=data)
}

get_display_type <- function(logical_name) {
  toupper(get_variable_suffix(logical_name))
}

############################
### LOOP DIAGNOSTIC TYPES ###
############################
for (type in plot_types) {

  vars1 <- get_type_variables(variables_exp1,type)
  vars2 <- get_type_variables(variables_exp2,type)

  if (length(vars1) == 0 && length(vars2) == 0) next

  if (length(vars1) != length(vars2)) {
    stop("Different number of ",type," variables between experiments. ",expname1,": ",paste(vars1,collapse=", ")," | ",expname2,": ",paste(vars2,collapse=", "))
  }

  #############################
  ### LOOP VARIABLE PAIRS ###
  #############################
  for (v in seq_along(vars1)) {

    variable1 <- vars1[v]
    variable2 <- vars2[v]

    message("---> Plotting ",variable1," vs ",variable2)

    #################
    ### READ DATA ###
    #################
    field_day1 <- list()
    field_day2 <- list()

    for (d in seq_along(seqDate)) {
      field_day1[[d]] <- read_daily_variable(expname1,exptype1,variable1,variables_exp1,type,seqDate[d])
      field_day2[[d]] <- read_daily_variable(expname2,exptype2,variable2,variables_exp2,type,seqDate[d])
    }

    ########################
    ### COMBINE TIME AXIS ###
    ########################
    nx1 <- dim(field_day1[[1]])[1]
    ny1 <- dim(field_day1[[1]])[2]
    nx2 <- dim(field_day2[[1]])[1]
    ny2 <- dim(field_day2[[1]])[2]

    if (nx1 != nx2 || ny1 != ny2) stop("Spatial dimensions differ between experiments for ",variable1," and ",variable2)

    nt1 <- sum(sapply(field_day1,function(x) dim(x)[3]))
    nt2 <- sum(sapply(field_day2,function(x) dim(x)[3]))

    if (nt1 != nt2) stop("Time dimensions differ between experiments for ",variable1," and ",variable2)

    data1 <- array(unlist(field_day1),dim=c(nx1,ny1,nt1))
    data2 <- array(unlist(field_day2),dim=c(nx2,ny2,nt2))

    ################
    ### LON / LAT ###
    ################
    file1 <- variable_file(variable1,variables_exp1,expname1,seqDate[1])
    ll <- read_lon_lat(file1)
    field_lon <- ll$lon
    field_lat <- ll$lat

    units <- get_plot_units(type)

    ####################
    ### TEMPORAL MEAN ###
    ####################
    field_var1 <- apply(data1,c(1,2),mean,na.rm=TRUE)
    field_var2 <- apply(data2,c(1,2),mean,na.rm=TRUE)

    ###################
    ### GLOBAL MEAN ###
    ###################
    if (type %in% c("mss","mss_from_mr")) {
      tmean_var1 <- global_mass_tg(data1,field_lon,field_lat)
      tmean_var2 <- global_mass_tg(data2,field_lon,field_lat)
    } else {
      tmean_var1 <- apply(data1,3,mean,na.rm=TRUE)
     tmean_var2 <- apply(data2,3,mean,na.rm=TRUE)
    }


    nt <- dim(data1)[3]

    tmean_tim <- seq.POSIXt(
      from=as.POSIXct(paste0(substr(sDate,1,4),"-",substr(sDate,5,6),"-",substr(sDate,7,8)," 00:00:00"),tz="UTC"),
      by="3 hours",
      length.out=nt
    )

    massdiag1 <- NULL
    massdiag2 <- NULL

    if (exists("massdiag_compare") && massdiag_compare && type %in% c("mss","mss_from_mr")) {
      massdiag1 <- read_massdiag_series(expname1,variable1,variables_exp1,seqDate)
      massdiag2 <- read_massdiag_series(expname2,variable2,variables_exp2,seqDate)
    }

    ###################
    ### DAILY CYCLE ###
    ###################
    hour <- as.integer(format(tmean_tim,"%H"))
    hours <- c(0,3,6,9,12,15,18,21)

    dhourmean_var1 <- sapply(hours,function(h) mean(tmean_var1[hour == h],na.rm=TRUE))
    dhourmean_var2 <- sapply(hours,function(h) mean(tmean_var2[hour == h],na.rm=TRUE))

    ###################
    ### MAP BREAKS ###
    ###################
    field_breaks <- positive_breaks(c(field_var1,field_var2),ncolors=200)
    field_breaks_diff <- difference_breaks(field_var2-field_var1,ncolors=200)

    ###################
    ### OUTPUT FILE ###
    ###################
    plot_category <- get_plot_category(variable1)
    plot_dir <- paste0(path_plot,plot_category,"/")
    dir.create(plot_dir,recursive=TRUE,showWarnings=FALSE)

    file_out <- paste0(plot_dir,gsub(" ","",plot_title),"_",variable1,"_vs_",variable2,"_",expname1,"-",expname2,"_",sDate,"-",eDate,".png")

    dpi <- 300
    png(file_out,width=(0.2+3*3.9+0.8+0.8)*dpi,height=(0.23+0.15+2+2.5)*dpi)

    layout(
      mat=matrix(c(1,1,1,1,1,1,2:13,14,14,14,14,15,15),4,6,byrow=TRUE),
      widths=c(0.2,3.9,3.9,0.8,3.9,0.8),
      heights=c(0.23,0.15,2,2.5)
    )

    ################
    ### HEADINGS ###
    ################
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0("Experiments: ",expname1," VS ",expname2,"   |   Type: ",get_display_type(variable1),"   |   Period: ",sDate,"-",eDate),col="grey50",cex=6,family="Century Gothic"); abline(h=c(0,1),col="grey50",lwd=3)
    par(mai=c(0,0,0,0)); plot.new()
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0(expname1," (",exptype1,")"),col="grey20",cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0(expname2," (",exptype2,")"),col="grey20",cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new()
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0(expname2," - ",expname1),col="grey20",cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new()

    title_text <- plot_type_title[[type]]
    if (type == "mmr") title_text <- paste0(title_text," @ ",plot_level," hPa")

    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,title_text,col="grey20",cex=5,family="Century Gothic",srt=90)

    #################
    ### MAP EXP 1 ###
    #################
    MapNC(filename_topo="",figure_box=figure_box,field_show_box=field_show_box,coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,field_value=field_var1,field_lon=field_lon,field_lat=field_lat,field_pallete_name="TROPOMI_NEW",field_breaks=field_breaks,field_units=units,field_pallete_starting_alpha=100,field_show_legend=FALSE)

    #################
    ### MAP EXP 2 ###
    #################
    MapNC(filename_topo="",figure_box=figure_box,field_show_box=field_show_box,coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,field_value=field_var2,field_lon=field_lon,field_lat=field_lat,field_pallete_name="TROPOMI_NEW",field_breaks=field_breaks,field_units=units,field_pallete_starting_alpha=100,field_show_legend=TRUE,field_legend_mai_right=1.8,field_legend_nlabels=7)

    ######################
    ### DIFFERENCE MAP ###
    ######################
    MapNC(filename_topo="",figure_box=figure_box,field_show_box=field_show_box,coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,field_value=field_var2-field_var1,field_lon=field_lon,field_lat=field_lat,field_pallete_name="MNMB",field_breaks=field_breaks_diff,field_units=units,field_pallete_starting_alpha=100,field_show_legend=TRUE,field_legend_mai_right=1.8,field_legend_nlabels=7)

    ##################
    ### TIMESERIES ###
    ##################
    par(mai=c(2,2,0,0.4),family="Century Gothic")

    x <- seq_along(tmean_tim)
    ts_values <- c(tmean_var1,tmean_var2)
    if (!is.null(massdiag1)) ts_values <- c(ts_values,massdiag1$value,massdiag2$value)
    yseq <- positive_axis_ticks(ts_values,n=10)
    plot(x,type="n",axes=FALSE,ann=FALSE,ylim=c(0,max(yseq$breaks)),yaxs="i")
    mtext("Time",side=1,line=12,cex=3.5)
    if (type %in% c("mss","mss_from_mr")) mtext("Global mass (Tg)",side=2,line=11,cex=3.5)

    IDx_labels <- which(format(tmean_tim,"%H") == "00" & format(tmean_tim,"%d") %in% c("01","05","10","15","20","25"))
    IDx_labels <- unique(c(1,IDx_labels,length(tmean_tim)))

    axis(1,at=x[IDx_labels],labels=format(tmean_tim[IDx_labels],"%Y-%m-%d"),cex.axis=4,line=4,lty=0)
    axis(1,at=x[IDx_labels],labels=FALSE,tck=0.01)
    axis(1,at=x[IDx_labels],labels=FALSE,tck=-0.01)
    axis(2,at=yseq$breaks,labels=yseq$labels,las=1,cex.axis=3)

    box(lwd=2)
    abline(h=yseq$breaks,lwd=1,col="grey")
    abline(v=x[IDx_labels],lwd=1,col="grey")

    lines(x,tmean_var1,lwd=5,col="blue")
    lines(x,tmean_var2,lwd=5,col="red")

    if (!is.null(massdiag1)) {
      massdiag_x1 <- as.numeric(difftime(massdiag1$time,tmean_tim[1],units="hours"))/3 + 1
      massdiag_x2 <- as.numeric(difftime(massdiag2$time,tmean_tim[1],units="hours"))/3 + 1
      valid1 <- massdiag_x1 >= 1 & massdiag_x1 <= length(tmean_tim) & is.finite(massdiag1$value)
      valid2 <- massdiag_x2 >= 1 & massdiag_x2 <= length(tmean_tim) & is.finite(massdiag2$value)
      points(massdiag_x1[valid1],massdiag1$value[valid1],pch=4,cex=4,lwd=4,col="blue")
      points(massdiag_x2[valid2],massdiag2$value[valid2],pch=4,cex=4,lwd=4,col="red")
    }

    if (is.null(massdiag1)) {
      legend("top",legend=c(paste0(expname1," (",variable1,")"),paste0(expname2," (",variable2,")")),lwd=7,pch=19,col=c("blue","red"),cex=3)
    } else {
      legend("top",legend=c(paste0(expname1," OUTPUT"),paste0(expname1," MASSDIA"),paste0(expname2," OUTPUT"),paste0(expname2," MASSDIA")),lwd=c(7,NA,7,NA),pch=c(19,4,19,4),pt.lwd=c(1,4,1,4),col=c("blue","blue","red","red"),cex=2.7,ncol=2)
    }

    ###################
    ### DAILY CYCLE ###
    ###################
    par(mai=c(2,2,0,0.4),family="Century Gothic")
    yseq <- positive_axis_ticks(c(dhourmean_var1,dhourmean_var2),n=10)
    plot(1:8,type="n",axes=FALSE,ann=FALSE,ylim=c(0,max(yseq$breaks)),yaxs="i")
    mtext("Time (3 hourly UTC)",side=1,line=12,cex=3.5)

    axis(1,at=1:8,labels=c("00","03","06","09","12","15","18","21"),cex.axis=4,line=4,lty=0)
    axis(1,at=1:8,labels=FALSE,tck=0.01)
    axis(1,at=1:8,labels=FALSE,tck=-0.01)
    axis(2,at=yseq$breaks,labels=yseq$labels,las=1,cex.axis=3)

    box(lwd=2)
    abline(h=yseq$breaks,lwd=1,col="grey")
    abline(v=1:8,lwd=1,col="grey")

    lines(1:8,dhourmean_var1,lwd=8); points(1:8,dhourmean_var1,pch=19,cex=2.8)
    lines(1:8,dhourmean_var2,lwd=8); points(1:8,dhourmean_var2,pch=19,cex=2.8)
    lines(1:8,dhourmean_var1,lwd=5,col="blue"); points(1:8,dhourmean_var1,pch=19,cex=2.2,col="blue")
    lines(1:8,dhourmean_var2,lwd=5,col="red"); points(1:8,dhourmean_var2,pch=19,cex=2.2,col="red")

    legend("top",legend=c(paste0(expname1," (",variable1,")"),paste0(expname2," (",variable2,")")),lwd=7,pch=19,col=c("blue","red"),cex=3)

    dev.off()

    ################
    ### COMPRESS ###
    ################
    file_tmp <- paste0(file_out,".tmp.png")
    compress(file_in=file_out,file_out=file_tmp)
    file.rename(file_tmp,file_out)

    message("---> Figure: ",file_out)
  }
}

############################
### COMPOSITE DEP_* PLOTS ###
############################
if (length(dep_variables) > 0) {

  for (dep_name in dep_variables) {

    dep_suffix <- sub("^dep_","",dep_name)
    flux_variables <- paste0(dep_fluxes,"_",dep_suffix)

    available_fluxes <- dep_fluxes[
      sapply(flux_variables,function(v) v %in% variables_exp1$logical_name && v %in% variables_exp2$logical_name)
    ]

    if (length(available_fluxes) == 0) {
      message("---> Skipping ",dep_name,": no common deposition fluxes available.")
      next
    }

    message("---> Plotting composite ",dep_name," using ",paste(available_fluxes,collapse=", "))

    dep_data <- list()

    for (flux in available_fluxes) {
      logical_name <- paste0(flux,"_",dep_suffix)

      d1 <- read_plot_data(expname1,exptype1,logical_name,variables_exp1,flux)
      d2 <- read_plot_data(expname2,exptype2,logical_name,variables_exp2,flux)

      if (d1$nx != d2$nx || d1$ny != d2$ny) stop("Spatial dimensions differ between experiments for ",logical_name)
      if (d1$nt != d2$nt) stop("Time dimensions differ between experiments for ",logical_name)

      file1 <- variable_file(logical_name,variables_exp1,expname1,seqDate[1])
      ll <- read_lon_lat(file1)

      field_var1 <- apply(d1$data,c(1,2),mean,na.rm=TRUE)
      field_var2 <- apply(d2$data,c(1,2),mean,na.rm=TRUE)
      tmean_var1 <- apply(d1$data,3,mean,na.rm=TRUE)
      tmean_var2 <- apply(d2$data,3,mean,na.rm=TRUE)

      tmean_tim <- seq.POSIXt(
        from=as.POSIXct(paste0(substr(sDate,1,4),"-",substr(sDate,5,6),"-",substr(sDate,7,8)," 00:00:00"),tz="UTC"),
        by="3 hours",
        length.out=d1$nt
      )

      hour <- as.integer(format(tmean_tim,"%H"))
      hours <- c(0,3,6,9,12,15,18,21)

      dep_data[[flux]] <- list(
        field_var1=field_var1,
        field_var2=field_var2,
        field_lon=ll$lon,
        field_lat=ll$lat,
        units=get_plot_units(flux),
        field_breaks=positive_breaks(c(field_var1,field_var2),ncolors=200),
        field_breaks_diff=difference_breaks(field_var2-field_var1,ncolors=200),
        tmean_var1=tmean_var1,
        tmean_var2=tmean_var2,
        tmean_tim=tmean_tim,
        dhourmean_var1=sapply(hours,function(h) mean(tmean_var1[hour == h],na.rm=TRUE)),
        dhourmean_var2=sapply(hours,function(h) mean(tmean_var2[hour == h],na.rm=TRUE))
      )
    }

    ###################
    ### OUTPUT FILE ###
    ###################
    plot_category <- get_plot_category(dep_name)
    plot_dir <- paste0(path_plot,plot_category,"/")
    dir.create(plot_dir,recursive=TRUE,showWarnings=FALSE)

    file_out <- paste0(plot_dir,gsub(" ","",plot_title),"_",dep_name,"_",expname1,"-",expname2,"_",sDate,"-",eDate,".png")

    ################
    ### LAYOUT ###
    ################
    nflux <- length(available_fluxes)

    # Exact extension of the individual six-column layout:
    # header row, experiment-heading row, one six-cell map row per flux, bottom TS/DC row.
    mat <- matrix(0,nrow=nflux+3,ncol=6)
    mat[1,] <- 1
    mat[2,] <- 2:7

    next_id <- 8
    for (i in seq_len(nflux)) {
      mat[2+i,] <- next_id:(next_id+5)
      next_id <- next_id+6
    }

    ts_id <- next_id
    dc_id <- next_id+1
    mat[nflux+3,] <- c(ts_id,ts_id,ts_id,ts_id,dc_id,dc_id)

    dpi <- 300
    png(
      file_out,
      width=(0.2+3*3.9+0.8+0.8)*dpi,
      height=(0.23+0.15+2*nflux+2.5)*dpi
    )

    layout(
      mat=mat,
      widths=c(0.2,3.9,3.9,0.8,3.9,0.8),
      heights=c(0.23,0.15,rep(2,nflux),2.5)
    )

    ################
    ### HEADINGS ###
    ################
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0("Experiments: ",expname1," VS ",expname2,"   |   Type: ",get_display_type(dep_name),"   |   Period: ",sDate,"-",eDate),col="grey50",cex=6,family="Century Gothic"); abline(h=c(0,1),col="grey50",lwd=3)
    par(mai=c(0,0,0,0)); plot.new()
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0(expname1," (",exptype1,")"),col="grey20",cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0(expname2," (",exptype2,")"),col="grey20",cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new()
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0(expname2," - ",expname1),col="grey20",cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new()

    #################
    ### MAP ROWS ###
    #################
    for (flux in available_fluxes) {
      z <- dep_data[[flux]]

      par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,plot_type_title[[flux]],col="grey20",cex=5,family="Century Gothic",srt=90)

      MapNC(filename_topo="",figure_box=figure_box,field_show_box=field_show_box,coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,field_value=z$field_var1,field_lon=z$field_lon,field_lat=z$field_lat,field_pallete_name="TROPOMI_NEW",field_breaks=z$field_breaks,field_units=z$units,field_pallete_starting_alpha=100,field_show_legend=FALSE)

      MapNC(filename_topo="",figure_box=figure_box,field_show_box=field_show_box,coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,field_value=z$field_var2,field_lon=z$field_lon,field_lat=z$field_lat,field_pallete_name="TROPOMI_NEW",field_breaks=z$field_breaks,field_units=z$units,field_pallete_starting_alpha=100,field_show_legend=TRUE,field_legend_mai_right=1.8,field_legend_nlabels=7)

      MapNC(filename_topo="",figure_box=figure_box,field_show_box=field_show_box,coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,field_value=z$field_var2-z$field_var1,field_lon=z$field_lon,field_lat=z$field_lat,field_pallete_name="MNMB",field_breaks=z$field_breaks_diff,field_units=z$units,field_pallete_starting_alpha=100,field_show_legend=TRUE,field_legend_mai_right=1.8,field_legend_nlabels=7)
    }

    ##################
    ### TIMESERIES ###
    ##################
    tmean_tim <- dep_data[[available_fluxes[1]]]$tmean_tim
    x <- seq_along(tmean_tim)
    all_ts <- unlist(lapply(available_fluxes,function(flux) c(dep_data[[flux]]$tmean_var1,dep_data[[flux]]$tmean_var2)))

    par(mai=c(2,2,0,0.4),family="Century Gothic")

    yseq <- axis_ticks(all_ts,n=10)
    pad <- diff(range(yseq$breaks))*0.05
    if (!is.finite(pad) || pad == 0) pad <- max(abs(yseq$breaks),na.rm=TRUE)*0.05
    if (!is.finite(pad) || pad == 0) pad <- 1

    plot(x,type="n",axes=FALSE,ann=FALSE,ylim=c(min(yseq$breaks)-pad,max(yseq$breaks)+pad),yaxs="i")
    mtext("Time",side=1,line=12,cex=3.5)
    if (type %in% c("mss","mss_from_mr")) mtext("Global mass (Tg)",side=2,line=11,cex=3.5)

    IDx_labels <- which(format(tmean_tim,"%H") == "00" & format(tmean_tim,"%d") %in% c("01","05","10","15","20","25"))
    IDx_labels <- unique(c(1,IDx_labels,length(tmean_tim)))

    axis(1,at=x[IDx_labels],labels=format(tmean_tim[IDx_labels],"%Y-%m-%d"),cex.axis=4,line=4,lty=0)
    axis(1,at=x[IDx_labels],labels=FALSE,tck=0.01)
    axis(1,at=x[IDx_labels],labels=FALSE,tck=-0.01)
    axis(2,at=yseq$breaks,labels=yseq$labels,las=1,cex.axis=3)

    box(lwd=2)
    abline(h=yseq$breaks,lwd=1,col="grey")
    abline(v=x[IDx_labels],lwd=1,col="grey")

    for (flux in available_fluxes) {
      lines(x,dep_data[[flux]]$tmean_var1,lwd=5,col=dep_flux_colors[flux],lty=1)
      lines(x,dep_data[[flux]]$tmean_var2,lwd=5,col=dep_flux_colors[flux],lty=2)
    }

    legend("topleft",legend=toupper(available_fluxes),lwd=6,col=dep_flux_colors[available_fluxes],lty=1,cex=2.3,bty="n")
    legend("topright",legend=c(expname1,expname2),lwd=6,col="grey20",lty=c(1,2),cex=2.3,bty="n")

    ###################
    ### DAILY CYCLE ###
    ###################
    all_dc <- unlist(lapply(available_fluxes,function(flux) c(dep_data[[flux]]$dhourmean_var1,dep_data[[flux]]$dhourmean_var2)))

    par(mai=c(2,2,0,0.4),family="Century Gothic")

    yseq <- axis_ticks(all_dc,n=10)
    pad <- diff(range(yseq$breaks))*0.05
    if (!is.finite(pad) || pad == 0) pad <- max(abs(yseq$breaks),na.rm=TRUE)*0.05
    if (!is.finite(pad) || pad == 0) pad <- 1

    plot(1:8,type="n",axes=FALSE,ann=FALSE,ylim=c(min(yseq$breaks)-pad,max(yseq$breaks)+pad),yaxs="i")
    mtext("Time (3 hourly UTC)",side=1,line=12,cex=3.5)

    axis(1,at=1:8,labels=c("00","03","06","09","12","15","18","21"),cex.axis=4,line=4,lty=0)
    axis(1,at=1:8,labels=FALSE,tck=0.01)
    axis(1,at=1:8,labels=FALSE,tck=-0.01)
    axis(2,at=yseq$breaks,labels=yseq$labels,las=1,cex.axis=3)

    box(lwd=2)
    abline(h=yseq$breaks,lwd=1,col="grey")
    abline(v=1:8,lwd=1,col="grey")

    for (flux in available_fluxes) {
      lines(1:8,dep_data[[flux]]$dhourmean_var1,lwd=5,col=dep_flux_colors[flux],lty=1)
      lines(1:8,dep_data[[flux]]$dhourmean_var2,lwd=5,col=dep_flux_colors[flux],lty=2)
    }

    legend("topleft",legend=toupper(available_fluxes),lwd=6,col=dep_flux_colors[available_fluxes],lty=1,cex=2.3,bty="n")
    legend("topright",legend=c(expname1,expname2),lwd=6,col="grey20",lty=c(1,2),cex=2.3,bty="n")

    dev.off()

    ################
    ### COMPRESS ###
    ################
    file_tmp <- paste0(file_out,".tmp.png")
    compress(file_in=file_out,file_out=file_tmp)
    file.rename(file_tmp,file_out)

    message("---> Composite figure: ",file_out)
  }
}

