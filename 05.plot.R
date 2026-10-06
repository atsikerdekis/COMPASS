############
### INIT ###
############
source(paste0(path_code,"04.preprocess.R"))

plot_title <- "Aerosol"
plot_level <- 1000

figure_box <- FALSE
field_show_box <- FALSE
gridlines <- 10
coastlineWorldFine_lwd <- 1

### Region settings
if (!exists("region")) region <- "global"
if (!exists("regions")) regions <- list(global=c(-180,180,-90,90))
if (!region %in% names(regions)) stop("Unknown region '",region,"'. Available regions: ",paste(names(regions),collapse=", "))

region_box <- regions[[region]]
regional_mode <- tolower(region) != "global"

if (regional_mode) {
  lonmin <- region_box[1]
  lonmax <- region_box[2]
  latmin <- region_box[3]
  latmax <- region_box[4]
  loncenter <- mean(c(lonmin,lonmax))
  latcenter <- mean(c(latmin,latmax))
  projection <- paste0("+proj=ortho +lon_0=",loncenter," +lat_0=",latcenter)
  gridlines <- 5
} else {
  lonmin <- 0
  lonmax <- 360
  latmin <- -90
  latmax <- 90
  projection <- "+proj=robin"
}

seqDate <- format(seq.Date(from=as.Date(sDate,format="%Y%m%d"),to=as.Date(eDate,format="%Y%m%d"),by="day"),"%Y%m%d")
massdiag_hours <- c(6,12,18)

########################
### VARIABLE FAMILIES ###
########################
plot_types <- c("mmr","ddp","sdm","wdl","wdc","mss","mss_from_mr","od","ngt")

plot_type_title <- c(
  mmr = "Mass mixing ratio",
  ddp = "Dry deposition",
  sdm = "Sedimentation",
  wdl = "Large-scale wet dep.",
  wdc = "Convective wet dep.",
  mss = "Column mass burden",
  mss_from_mr = "Column burden from mixing ratio",
  od = "Optical depth",
  ngt = "Negative fixer"
)

dep_fluxes <- c("ddp","sdm","wdl","wdc","ngt")
dep_flux_colors <- c(ddp="#846040",sdm="#D78C6A",wdl="#74ACE8",wdc="#3E7DD1",ngt="#3FC13A",wdep="#8B6FB5")
requested_individual_variables <- variables_requested[!startsWith(variables_requested,"dep_")]
optical_plot_variables <- if (exists("optical_variables")) optical_variables else character(0)
column_plot_variables <- if (exists("column_variables")) column_variables else character(0)
wat_plot_variables <- if (exists("wat_variables")) wat_variables else character(0)
panel_species_plot_variables <- if (exists("panel_species_variables")) panel_species_variables else character(0)

ts_panel_mai <- c(1.95,2.45,0.05,0.18)
ts_xlabel_line <- 4.8
ts_ylabel_line <- 6.2
ts_xlabel_cex <- 2.1
ts_ylabel_cex <- 2.1
ts_xaxis_cex <- 4.0
ts_yaxis_cex <- 3.0
map_ylabel_cex <- 2.2
map_title_cex <- 1.25

panel_header_cex <- 4.8
panel_map_title_cex <- 3.6
panel_stats_cex <- 2.8
panel_legend_cex <- 2.35
panel_map_width <- 3.9
panel_legend_width <- 0.55
panel_map_height <- 2.1
panel_header_height <- 0.50

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
    x <- unique(variable_table$logical_name[startsWith(variable_table$logical_name,"mss_") & !startsWith(variable_table$logical_name,"mss_from_mr_")])
  } else {
    x <- unique(variable_table$logical_name[startsWith(variable_table$logical_name,paste0(type,"_"))])
  }
  x[x %in% requested_individual_variables]
}

nice_positive_max <- function(x,n=6) {
  xmax <- max(x,na.rm=TRUE)
  if (!is.finite(xmax) || xmax <= 0) return(1)
  max(pretty(c(0,xmax),n=n))
}

positive_breaks <- function(x,ncolors=200) seq(0,nice_positive_max(x),length.out=ncolors+1)

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
  if (type == "od") return(" ")
  " "
}

get_plot_category <- function(logical_name) {
  suffix <- get_variable_suffix(logical_name)
  if (exptype1 != exptype2) return("per_species")
  if (exptype1 == "AER" && exptype2 == "AER") {
    if (suffix %in% common_HAM_AER_species) return("per_species")
    if (suffix %in% aer_names) return("per_tracer")
    stop("Could not determine AER plot category for ",logical_name)
  }
  if (suffix %in% c("soluble","insoluble")) return("per_solubility")
  if (suffix %in% ham_modes) return("per_mode")
  if (grepl(paste0("_(",paste(ham_modes,collapse="|"),")$"),suffix)) return("per_tracer")
  return("per_species")
}

read_daily_variable <- function(expname,exptype,logical_name,variable_table,type,date) {
  file <- variable_file(logical_name,variable_table,expname,date)
  if (!file.exists(file)) stop("Input file not found: ",file)
  if (type == "mmr") return(read_variable_level(file,exptype,logical_name,variable_table,plot_level))
  if (type == "mss_from_mr") return(read_column_burden(file,exptype,logical_name,variable_table))
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

get_display_type <- function(logical_name) toupper(get_variable_suffix(logical_name))
get_massdiag_pch <- function(exptype) if (exptype == "HAM") 4 else 1
get_midnight_idx <- function(tt) which(format(tt,"%H") == "00")
region_title <- if (regional_mode) paste0("   |   Region: ",region) else ""
region_file_tag <- if (regional_mode) paste0("_",gsub("[^A-Za-z0-9_-]","_",region)) else ""

region_dir <- if (regional_mode) gsub("[^A-Za-z0-9_-]","_",region) else "global"
make_plot_dir <- function(category) {
  out <- paste0(path_plot,category,"/",region_dir,"/")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  out
}

timeseries_axis_indices <- function(tt,nmax=7) {
  n <- length(tt)
  if (n <= nmax) return(seq_len(n))
  unique(round(seq(1,n,length.out=nmax)))
}

add_mean_square <- function(values,col,cex=3.8) {
  m <- mean(values,na.rm=TRUE)
  if (!is.finite(m)) return(invisible(NULL))
  usr <- par("usr")
  xmean <- usr[1] + 0.018*(usr[2]-usr[1])
  points(xmean,m,pch=15,cex=cex*1.18,col="black",xpd=FALSE)
  points(xmean,m,pch=15,cex=cex,col=col,xpd=FALSE)
}

read_lifetime_components <- function(expname,exptype,suffix,variable_table) {

  mss_name <- paste0("mss_",suffix)
  mss_data <- read_plot_data(expname,exptype,mss_name,variable_table,"mss")

  dep_total <- NULL
  available_fluxes <- character(0)

  for (flux in dep_fluxes) {
    dep_name <- paste0(flux,"_",suffix)
    if (!dep_name %in% variable_table$logical_name) next

    d <- read_plot_data(expname,exptype,dep_name,variable_table,flux)

    if (d$nx != mss_data$nx || d$ny != mss_data$ny || d$nt != mss_data$nt)
      stop("MSS/deposition dimensions differ for lifetime_",suffix)

    if (is.null(dep_total)) dep_total <- d$data else dep_total <- dep_total+d$data
    available_fluxes <- c(available_fluxes,flux)
  }

  if (is.null(dep_total))
    stop("No deposition diagnostics available for lifetime_",suffix)

  lifetime_field <- mss_data$data/(dep_total*86400)
  lifetime_field[!is.finite(lifetime_field) | lifetime_field < 0] <- NA_real_

  list(
    mss=mss_data$data,
    dep=dep_total,
    lifetime=lifetime_field,
    nx=mss_data$nx,
    ny=mss_data$ny,
    nt=mss_data$nt,
    fluxes=available_fluxes
  )
}


get_panel_species_label <- function(x) {
  switch(x,
         total = "TOTAL",
         du = "DU",
         ss = "SS",
         pom = "POM",
         bc = "BC",
         so4 = "SO4",
         ni = "NI",
         am = "AM",
         wat = "WAT",
         blank = "",
         toupper(x))
}

panel_tropomi_colors <- function(n) {
  colorRampPalette(c(
    "#EBF7FD","#C9EAFB","#ADE1F6","#9BD0EE","#89C3E6","#7EBAE2",
    "#71AFDC","#60A3D5","#62B49B","#8ACE64","#D0DF6F","#FAE771",
    "#FACD64","#F7B95B","#F8A750","#FB9548","#F6813F","#EA5B3E",
    "#CA1112","#A30008","#820005","#640203","#4E0002","#360202",
    "#240302","#100202"
  ))(n)
}

panel_mnmb_colors <- function(n) {
  colorRampPalette(c("navy","blue","lightskyblue","white","palevioletred1","red","red4"))(n)
}

panel_fge_colors <- function(n) {
  colorRampPalette(rev(c("darkred","red","orange","yellow","greenyellow","green3","green4")))(n)
}

panel_palette_colors <- function(name,n) {
  if (name == "MNMB") return(panel_mnmb_colors(n))
  if (name == "FGE") return(panel_fge_colors(n))
  panel_tropomi_colors(n)
}

panel_safe_positive_breaks <- function(fields,ncolors=200) {
  vals <- unlist(fields,use.names=FALSE)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(seq(0,1,length.out=ncolors+1))
  positive_breaks(vals,ncolors=ncolors)
}

panel_format_stat <- function(x,percent=FALSE) {
  if (!is.finite(x)) return("NA")
  if (percent) return(paste0(formatC(x,format="f",digits=1),"%"))
  format(x,scientific=TRUE,digits=2)
}

panel_draw_map <- function(field_value,lon,lat,breaks,palette_name,label,
                           percent_stats=FALSE) {
  if (regional_mode)
    field_value <- mask_region_field(field_value,lon,lat,region_box)

  MapNC(
    filename_topo="",
    figure_box=figure_box,
    field_show_box=field_show_box,
    coastlineWorldFine_lwd=coastlineWorldFine_lwd,
    gridlines=gridlines,
    projection=projection,
    lonmax=lonmax,lonmin=lonmin,
    latmax=latmax,latmin=latmin,
    drawMapBox=regional_mode,
    field_value=field_value,
    field_lon=lon,
    field_lat=lat,
    field_pallete_name=palette_name,
    field_breaks=breaks,
    field_units="",
    field_pallete_starting_alpha=100,
    field_show_legend=FALSE,
    field_value_mean_global=FALSE,
    title_main=""
  )

  if (nzchar(label)) {
    usr <- par("usr")
    dx <- usr[2]-usr[1]
    dy <- usr[4]-usr[3]

    ### Lower, smaller title than the standard MapNC title.
    tx <- usr[1]+0.020*dx
    ty <- usr[4]-0.075*dy
    text(tx,ty,label,adj=c(0,1),cex=panel_map_title_cex,
         font=2,family="Century Gothic",col="white")
    text(tx+0.002*dx,ty-0.002*dy,label,adj=c(0,1),cex=panel_map_title_cex,
         font=2,family="Century Gothic",col="grey15")

    mn <- mean(field_value,na.rm=TRUE)
    sdv <- sd(field_value,na.rm=TRUE)
    if (is.nan(mn)) mn <- NA_real_
    if (is.nan(sdv)) sdv <- NA_real_

    sy <- usr[3]+0.055*dy
    text(usr[1]+0.035*dx,sy,
         paste0("MN\n",panel_format_stat(mn,percent_stats)),
         adj=c(0,0),cex=panel_stats_cex,font=2,family="Century Gothic")
    text(usr[2]-0.035*dx,sy,
         paste0("SD\n",panel_format_stat(sdv,percent_stats)),
         adj=c(1,0),cex=panel_stats_cex,font=2,family="Century Gothic")
  }
}

panel_draw_legend <- function(breaks,palette_name,units="",percent=FALSE) {
  cols <- panel_palette_colors(palette_name,length(breaks)-1)

  ### Use the same palette renderer as MapNC so the 4x3 panels have the
  ### standard COMPASS narrow legend with triangular end caps.
  par(bg="#FFFFFFFF",family="Century Gothic")
  par(mai=c(0.35,0.05,0.35,0.70))

  try({
    drawPalette(
      at=seq_along(breaks),
      labels=breaks,
      fullpage=TRUE,
      col=cols,
      las=1,
      mai=c(0,0.05,0.10,0.55),
      drawTriangles=TRUE,
      cex=0
    )
  },silent=TRUE)

  mylegend_at <- unique(round(seq(1,length(breaks),length.out=5)))
  mylegend_values <- breaks[mylegend_at]

  if (min(breaks,na.rm=TRUE) < 0 && max(breaks,na.rm=TRUE) > 0) {
    izero <- which.min(abs(mylegend_values))
    mylegend_values[izero] <- 0
  }

  if (percent) {
    mylegend_labels <- paste0(formatC(mylegend_values,format="f",digits=0),"%")
  } else {
    mylegend_labels <- format(mylegend_values,scientific=TRUE,digits=2)
  }

  axis(
    4,
    at=mylegend_at,
    labels=mylegend_labels,
    las=1,
    cex.axis=panel_legend_cex,
    tick=TRUE,
    family="Century Gothic"
  )

  if (nzchar(units) && !percent)
    mtext(units,side=3,line=1.2,adj=1,
          cex=panel_legend_cex/1.8,family="Century Gothic")
}

panel_validation_plot_settings <- function(variable,kind) {
  if (variable == "AOD550") {
    if (kind == "obs") return(list(
      breaks=c(0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.50,0.60,0.80,1.00),
      palette="TROPOMI_NEW",units=""))
    if (kind == "me") return(list(
      breaks=seq(-0.8,0.8,length.out=201),palette="MNMB",units=""))
    return(list(breaks=seq(0,1.0,length.out=201),palette="FGE",units=""))
  }

  if (variable == "AE550to860") {
    if (kind == "obs") return(list(
      breaks=seq(0,2.0,length.out=201),palette="TROPOMI_NEW",units=""))
    if (kind == "me") return(list(
      breaks=seq(-1,1,length.out=201),palette="MNMB",units=""))
    return(list(breaks=seq(0,1.0,length.out=201),palette="FGE",units=""))
  }

  stop("No validation plot settings for ",variable)
}

plot_species_map_grid <- function(field_list,validation,expname,exptype,panel_name,units=" ") {
  title_text <- get_panel_species_definition(panel_name)$title
  plot_dir <- make_plot_dir("optics")
  file_out <- paste0(
    plot_dir,"Panel_",panel_name,"_",expname,region_file_tag,"_",sDate,"-",eDate,".png"
  )

  if (length(field_list) != 9)
    stop("Species panel plot requires exactly 9 fields; found ",length(field_list),
         " for ",panel_name," / ",expname)

  ratio_panel <- panel_name %in% c("aodratio_per_species","massratio_per_species")

  ### Convert all ratio panels, including TOTAL, to percent before computing
  ### statistics and plotting. This gives one 0--100% scale everywhere.
  if (ratio_panel)
    field_list <- lapply(field_list,function(x) x*100)

  populated <- field_list[names(field_list) != "blank"]
  species_breaks <- if (ratio_panel) seq(0,100,length.out=201) else
    panel_safe_positive_breaks(populated,200)

  species_units <- if (ratio_panel) "%" else units

  ### Title + 4 map rows. Every map has its own narrow legend cell.
  panel_layout <- matrix(
    c(
      1,1,1,1,1,1,
      2,3,4,5,6,7,
      8,9,10,11,12,13,
      14,15,16,17,18,19,
      20,21,22,23,24,25
    ),
    nrow=5,byrow=TRUE
  )

  dpi <- 300
  panel_width <- 3*(panel_map_width+panel_legend_width)
  panel_height <- panel_header_height+4*panel_map_height
  png(file_out,width=panel_width*dpi,height=panel_height*dpi)
  layout(
    panel_layout,
    widths=rep(c(panel_map_width,panel_legend_width),3),
    heights=c(panel_header_height,rep(panel_map_height,4))
  )

  par(mai=c(0,0,0,0),family="Century Gothic")
  plot.new()
  text(
    0.5,0.5,
    paste0(
      "Experiment: ",expname," (",exptype,")   |   Type: ",title_text,
      "   |   Validation: ",panel_validation_variable,
      " (",paste(panel_satellites,collapse="+"),")",
      region_title,"   |   Period: ",sDate,"-",eDate
    ),
    col="grey40",cex=panel_header_cex,family="Century Gothic"
  )
  abline(h=c(0,1),col="grey50",lwd=2)

  ### Row 1: satellite ensemble, model mean error, model MAE.
  validation_fields <- list(
    "SAT ENSEMBLE"=validation$obs,
    "ME"=validation$me,
    "MAE"=validation$mae
  )
  validation_kinds <- c("obs","me","mae")
  for (i in seq_along(validation_fields)) {
    settings <- panel_validation_plot_settings(panel_validation_variable,validation_kinds[i])
    panel_draw_map(
      validation_fields[[i]],validation$lon,validation$lat,
      settings$breaks,settings$palette,names(validation_fields)[i],FALSE
    )
    panel_draw_legend(settings$breaks,settings$palette,settings$units,FALSE)
  }

  ### Rows 2--4: species fields, one identical scale for all populated panels.
  for (nm in names(field_list)) {
    label <- get_panel_species_label(nm)
    if (nm == "blank") {
      par(mai=c(0,0,0,0)); plot.new()
      par(mai=c(0,0,0,0)); plot.new()
      next
    }

    panel_draw_map(
      field_list[[nm]],field_lon,field_lat,
      species_breaks,"TROPOMI_NEW",label,ratio_panel
    )
    panel_draw_legend(
      species_breaks,"TROPOMI_NEW",species_units,ratio_panel
    )
  }

  dev.off()
  tmp <- paste0(file_out,".tmp.png")
  compress(file_in=file_out,file_out=tmp)
  file.rename(tmp,file_out)
  message("---> Panel figure: ",file_out)
}

############################
### LOOP DIAGNOSTIC TYPES ###
############################
for (type in plot_types) {
  vars1 <- get_type_variables(variables_exp1,type)
  vars2 <- get_type_variables(variables_exp2,type)
  if (length(vars1) == 0 && length(vars2) == 0) next
  if (length(vars1) != length(vars2)) stop("Different number of ",type," variables between experiments. ",expname1,": ",paste(vars1,collapse=", ")," | ",expname2,": ",paste(vars2,collapse=", "))

  for (v in seq_along(vars1)) {
    variable1 <- vars1[v]
    variable2 <- vars2[v]
    message("---> Plotting ",variable1," vs ",variable2,if (regional_mode) paste0(" for ",region) else "")

    field_day1 <- list()
    field_day2 <- list()
    for (d in seq_along(seqDate)) {
      field_day1[[d]] <- read_daily_variable(expname1,exptype1,variable1,variables_exp1,type,seqDate[d])
      field_day2[[d]] <- read_daily_variable(expname2,exptype2,variable2,variables_exp2,type,seqDate[d])
    }

    nx1 <- dim(field_day1[[1]])[1]; ny1 <- dim(field_day1[[1]])[2]
    nx2 <- dim(field_day2[[1]])[1]; ny2 <- dim(field_day2[[1]])[2]
    if (nx1 != nx2 || ny1 != ny2) stop("Spatial dimensions differ between experiments for ",variable1," and ",variable2)
    nt1 <- sum(sapply(field_day1,function(x) dim(x)[3])); nt2 <- sum(sapply(field_day2,function(x) dim(x)[3]))
    if (nt1 != nt2) stop("Time dimensions differ between experiments for ",variable1," and ",variable2)
    data1 <- array(unlist(field_day1),dim=c(nx1,ny1,nt1))
    data2 <- array(unlist(field_day2),dim=c(nx2,ny2,nt2))

    file1 <- variable_file(variable1,variables_exp1,expname1,seqDate[1])
    ll <- read_lon_lat(file1)
    field_lon <- ll$lon
    field_lat <- ll$lat
    units <- get_plot_units(type)

    field_var1 <- apply(data1,c(1,2),mean,na.rm=TRUE)
    field_var2 <- apply(data2,c(1,2),mean,na.rm=TRUE)

    if (regional_mode) {
      field_plot1 <- mask_region_field(field_var1,field_lon,field_lat,region_box)
      field_plot2 <- mask_region_field(field_var2,field_lon,field_lat,region_box)
    } else {
      field_plot1 <- field_var1
      field_plot2 <- field_var2
    }

    if (type %in% c("mss","mss_from_mr")) {
      if (regional_mode) {
        tmean_var1 <- regional_mass_tg(data1,field_lon,field_lat,region_box)
        tmean_var2 <- regional_mass_tg(data2,field_lon,field_lat,region_box)
      } else {
        tmean_var1 <- global_mass_tg(data1,field_lon,field_lat)
        tmean_var2 <- global_mass_tg(data2,field_lon,field_lat)
      }
    } else {
      if (regional_mode) {
        tmean_var1 <- regional_mean(data1,field_lon,field_lat,region_box)
        tmean_var2 <- regional_mean(data2,field_lon,field_lat,region_box)
      } else {
        tmean_var1 <- apply(data1,3,mean,na.rm=TRUE)
        tmean_var2 <- apply(data2,3,mean,na.rm=TRUE)
      }
    }

    nt <- dim(data1)[3]
    tmean_tim <- seq.POSIXt(from=as.POSIXct(paste0(substr(sDate,1,4),"-",substr(sDate,5,6),"-",substr(sDate,7,8)," 00:00:00"),tz="UTC"),by="3 hours",length.out=nt)


    ### MASSDIA is a global diagnostic and is intentionally disabled for regional plots
    massdiag1 <- NULL
    massdiag2 <- NULL
    if (!regional_mode && exists("massdiag_compare") && massdiag_compare && type %in% c("mss","mss_from_mr")) {
      massdiag1 <- read_massdiag_series_hours(expname1,variable1,variables_exp1,seqDate,hours=massdiag_hours,column="TOT_MASS")
      massdiag2 <- read_massdiag_series_hours(expname2,variable2,variables_exp2,seqDate,hours=massdiag_hours,column="TOT_MASS")
    }

    hour <- as.integer(format(tmean_tim,"%H"))
    hours <- c(0,3,6,9,12,15,18,21)
    dhourmean_var1 <- sapply(hours,function(h) mean(tmean_var1[hour == h],na.rm=TRUE))
    dhourmean_var2 <- sapply(hours,function(h) mean(tmean_var2[hour == h],na.rm=TRUE))

    ts_plot <- aggregate_pair_for_plot(tmean_var1,tmean_var2,tmean_tim,variable1)
    plot_tmean_tim <- ts_plot$time
    plot_tmean_var1 <- ts_plot$value1
    plot_tmean_var2 <- ts_plot$value2

    ### MASSDIA stays at its native diagnostic timestamps; hide it when the
    ### requested output time series is daily/monthly to avoid mixing cadences.
    if (ts_plot$aggregation != "3hourly") {
      massdiag1 <- NULL
      massdiag2 <- NULL
    }

    field_breaks <- positive_breaks(c(field_plot1,field_plot2),ncolors=200)
    field_breaks_diff <- difference_breaks(field_plot2-field_plot1,ncolors=200)

    plot_category <- get_plot_category(variable1)
    plot_dir <- make_plot_dir(plot_category)
    file_out <- paste0(plot_dir,gsub(" ","",plot_title),"_",variable1,"_vs_",variable2,"_",expname1,"-",expname2,region_file_tag,"_",sDate,"-",eDate,".png")

    dpi <- 300
    png(file_out,width=(0.2+3*3.9+0.8+0.8)*dpi,height=(0.23+0.15+2+2.5)*dpi)
    layout(mat=matrix(c(1,1,1,1,1,1,2:13,14,14,14,14,15,15),4,6,byrow=TRUE),widths=c(0.2,3.9,3.9,0.8,3.9,0.8),heights=c(0.23,0.15,2,2.5))

    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0("Experiments: ",expname1," VS ",expname2,"   |   Type: ",get_display_type(variable1),region_title,"   |   Period: ",sDate,"-",eDate),col="grey50",cex=6,family="Century Gothic"); abline(h=c(0,1),col="grey50",lwd=3)
    par(mai=c(0,0,0,0)); plot.new()
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0(expname1," (",exptype1,")"),col="grey20",cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0(expname2," (",exptype2,")"),col="grey20",cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new()
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0(expname2," - ",expname1),col="grey20",cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new()
    title_text <- plot_type_title[[type]]
    if (type == "mmr") title_text <- paste0(title_text," @ ",plot_level," hPa")
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,title_text,col="grey20",cex=map_ylabel_cex,family="Century Gothic",srt=90)

    MapNC(filename_topo="",figure_box=figure_box,field_show_box=field_show_box,coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,drawMapBox=regional_mode,field_value=field_plot1,field_lon=field_lon,field_lat=field_lat,field_pallete_name="TROPOMI_NEW",field_breaks=field_breaks,field_units=units,field_pallete_starting_alpha=100,field_show_legend=FALSE)
    MapNC(filename_topo="",figure_box=figure_box,field_show_box=field_show_box,coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,drawMapBox=regional_mode,field_value=field_plot2,field_lon=field_lon,field_lat=field_lat,field_pallete_name="TROPOMI_NEW",field_breaks=field_breaks,field_units=units,field_pallete_starting_alpha=100,field_show_legend=TRUE,field_legend_mai_right=1.8,field_legend_nlabels=7)
    MapNC(filename_topo="",figure_box=figure_box,field_show_box=field_show_box,coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,drawMapBox=regional_mode,field_value=field_plot2-field_plot1,field_lon=field_lon,field_lat=field_lat,field_pallete_name="MNMB",field_breaks=field_breaks_diff,field_units=units,field_pallete_starting_alpha=100,field_show_legend=TRUE,field_legend_mai_right=1.8,field_legend_nlabels=7)

    par(mai=ts_panel_mai,family="Century Gothic")
    x <- seq_along(plot_tmean_tim)
    ts_values <- c(plot_tmean_var1,plot_tmean_var2)
    if (!is.null(massdiag1)) ts_values <- c(ts_values,massdiag1$value,massdiag2$value)
    yseq <- positive_axis_ticks(ts_values,n=10)

    plot(x,type="n",axes=FALSE,ann=FALSE,ylim=c(0,max(yseq$breaks)),yaxs="i")
    mtext("Time",side=1,line=ts_xlabel_line,cex=ts_xlabel_cex)

    if (type %in% c("mss","mss_from_mr"))
      mtext(if (regional_mode) paste0(region," mass (Tg)") else "Global mass (Tg)",
            side=2,line=ts_ylabel_line,cex=ts_ylabel_cex)

    IDx_labels <- timeseries_axis_indices(plot_tmean_tim)
    axis(1,at=x[IDx_labels],
         labels=format(plot_tmean_tim[IDx_labels],"%Y-%m-%d"),
         cex.axis=ts_xaxis_cex,line=0,lty=0)
    axis(1,at=x[IDx_labels],labels=FALSE,tck=0.01)
    axis(1,at=x[IDx_labels],labels=FALSE,tck=-0.01)
    axis(2,at=yseq$breaks,labels=yseq$labels,las=1,cex.axis=ts_yaxis_cex)

    box(lwd=2)
    abline(h=yseq$breaks,lwd=1,col="grey")
    abline(v=x[IDx_labels],lwd=1,col="grey")

    lines(x,plot_tmean_var1,lwd=5,col="blue")
    points(x,plot_tmean_var1,pch=19,cex=1.8,col="blue")
    lines(x,plot_tmean_var2,lwd=5,col="red")
    points(x,plot_tmean_var2,pch=19,cex=1.8,col="red")

    ### Period means: pch=15 at the far-left side of the time-series panel.
    add_mean_square(plot_tmean_var1,"blue")
    add_mean_square(plot_tmean_var2,"red")

    if (!is.null(massdiag1)) {
      massdiag_x1 <- as.numeric(difftime(massdiag1$time,plot_tmean_tim[1],units="hours"))/3 + 1
      massdiag_x2 <- as.numeric(difftime(massdiag2$time,plot_tmean_tim[1],units="hours"))/3 + 1
      valid1 <- massdiag_x1 >= 1 & massdiag_x1 <= length(plot_tmean_tim) & is.finite(massdiag1$value)
      valid2 <- massdiag_x2 >= 1 & massdiag_x2 <= length(plot_tmean_tim) & is.finite(massdiag2$value)
      pch1 <- get_massdiag_pch(exptype1)
      pch2 <- get_massdiag_pch(exptype2)
      points(massdiag_x1[valid1],massdiag1$value[valid1],
             pch=pch1,cex=2.8,lwd=2.5,col="blue")
      points(massdiag_x2[valid2],massdiag2$value[valid2],
             pch=pch2,cex=2.8,lwd=2.5,col="red")
    }

    if (is.null(massdiag1)) {
      legend("top",legend=c(expname1,expname2), lwd=4,col=c("blue","red"),cex=1.7, bty="n")
    } else {
      legend("top",
             legend=c(paste0(expname1," OUTPUT"),
                      paste0(expname1," MASSDIA"),
                      paste0(expname2," OUTPUT"),
                      paste0(expname2," MASSDIA")),
             lwd=c(5,NA,5,NA),
             pch=c(19,get_massdiag_pch(exptype1),19,get_massdiag_pch(exptype2)),
             pt.lwd=c(1,2.5,1,2.5),
             col=c("blue","blue","red","red"),cex=2.5,ncol=2)
    }

    par(mai=ts_panel_mai,family="Century Gothic")
    yseq <- positive_axis_ticks(c(dhourmean_var1,dhourmean_var2),n=10)
    plot(1:8,type="n",axes=FALSE,ann=FALSE,ylim=c(0,max(yseq$breaks)),yaxs="i")
    mtext("Time (3 hourly UTC)",side=1,line=ts_xlabel_line,cex=ts_xlabel_cex)
    axis(1,at=1:8,labels=c("00","03","06","09","12","15","18","21"),cex.axis=ts_xaxis_cex,line=0,lty=0)
    axis(1,at=1:8,labels=FALSE,tck=0.01); axis(1,at=1:8,labels=FALSE,tck=-0.01)
    axis(2,at=yseq$breaks,labels=yseq$labels,las=1,cex.axis=ts_yaxis_cex)
    box(lwd=2); abline(h=yseq$breaks,lwd=1,col="grey"); abline(v=1:8,lwd=1,col="grey")
    lines(1:8,dhourmean_var1,lwd=5,col="blue")
    points(1:8,dhourmean_var1,pch=19,cex=1.8,col="blue")
    lines(1:8,dhourmean_var2,lwd=5,col="red")
    points(1:8,dhourmean_var2,pch=19,cex=1.8,col="red")
    legend("top",legend=c(expname1,expname2),lwd=5,col=c("blue","red"),cex=3)

    dev.off()
    file_tmp <- paste0(file_out,".tmp.png")
    compress(file_in=file_out,file_out=file_tmp)
    file.rename(file_tmp,file_out)
    message("---> Figure: ",file_out)
  }
}


############################
### RELATIVE HUMIDITY PLOTS ###
############################
if (length(rh_variables) > 0) {

  rh_height_label <- function(logical_name) {
    logical_name <- strip_time_aggregation(logical_name)
    if (logical_name == "rh") return("~surface")
    x <- sub("^rh_","",logical_name)
    paste0("~",sub("m$","",x)," m")
  }

  for (logical_name in rh_variables) {

    def <- get_rh_definition(logical_name)
    message("---> Plotting ",logical_name," (ML",def$model_level,", ",rh_height_label(logical_name),")",
            if (regional_mode) paste0(" for ",region) else "")

    field_day1 <- list()
    field_day2 <- list()

    for (d in seq_along(seqDate)) {
      field_day1[[d]] <- read_relative_humidity(expname1,seqDate[d],logical_name)
      field_day2[[d]] <- read_relative_humidity(expname2,seqDate[d],logical_name)
    }

    nx1 <- dim(field_day1[[1]])[1]; ny1 <- dim(field_day1[[1]])[2]
    nx2 <- dim(field_day2[[1]])[1]; ny2 <- dim(field_day2[[1]])[2]

    if (nx1 != nx2 || ny1 != ny2)
      stop("Spatial dimensions differ between experiments for ",logical_name)

    nt1 <- sum(sapply(field_day1,function(x) dim(x)[3]))
    nt2 <- sum(sapply(field_day2,function(x) dim(x)[3]))

    if (nt1 != nt2)
      stop("Time dimensions differ between experiments for ",logical_name)

    data1 <- array(unlist(field_day1),dim=c(nx1,ny1,nt1))
    data2 <- array(unlist(field_day2),dim=c(nx2,ny2,nt2))

    file1 <- rh_ml_file(expname1,seqDate[1],logical_name)
    ll <- read_lon_lat(file1)
    field_lon <- ll$lon
    field_lat <- ll$lat

    field_var1 <- apply(data1,c(1,2),mean,na.rm=TRUE)
    field_var2 <- apply(data2,c(1,2),mean,na.rm=TRUE)

    if (regional_mode) {
      field_plot1 <- mask_region_field(field_var1,field_lon,field_lat,region_box)
      field_plot2 <- mask_region_field(field_var2,field_lon,field_lat,region_box)
      tmean_var1 <- regional_mean(data1,field_lon,field_lat,region_box)
      tmean_var2 <- regional_mean(data2,field_lon,field_lat,region_box)
    } else {
      field_plot1 <- field_var1
      field_plot2 <- field_var2
      tmean_var1 <- apply(data1,3,mean,na.rm=TRUE)
      tmean_var2 <- apply(data2,3,mean,na.rm=TRUE)
    }

    nt <- dim(data1)[3]
    tmean_tim <- seq.POSIXt(
      from=as.POSIXct(paste0(substr(sDate,1,4),"-",substr(sDate,5,6),"-",substr(sDate,7,8)," 00:00:00"),tz="UTC"),
      by="3 hours",length.out=nt
    )

    hour <- as.integer(format(tmean_tim,"%H"))
    hours <- c(0,3,6,9,12,15,18,21)
    dhourmean_var1 <- sapply(hours,function(h) mean(tmean_var1[hour == h],na.rm=TRUE))
    dhourmean_var2 <- sapply(hours,function(h) mean(tmean_var2[hour == h],na.rm=TRUE))

    ts_plot <- aggregate_pair_for_plot(tmean_var1,tmean_var2,tmean_tim,logical_name)
    plot_tmean_tim <- ts_plot$time
    plot_tmean_var1 <- ts_plot$value1
    plot_tmean_var2 <- ts_plot$value2

    ### RH uses a fixed 0-100% colour scale so experiments/heights are directly comparable.
    field_breaks <- seq(0,100,length.out=201)
    field_breaks_diff <- difference_breaks(field_plot2-field_plot1,ncolors=200)

    plot_dir <- make_plot_dir("meteorology")

    file_out <- paste0(
      plot_dir,"Meteorology_",logical_name,"_vs_",logical_name,"_",
      expname1,"-",expname2,region_file_tag,"_",sDate,"-",eDate,".png"
    )

    dpi <- 300
    png(file_out,width=(0.2+3*3.9+0.8+0.8)*dpi,height=(0.23+0.15+2+2.5)*dpi)
    layout(
      mat=matrix(c(1,1,1,1,1,1,2:13,14,14,14,14,15,15),4,6,byrow=TRUE),
      widths=c(0.2,3.9,3.9,0.8,3.9,0.8),
      heights=c(0.23,0.15,2,2.5)
    )

    par(mai=c(0,0,0,0))
    plot.new()
    text(
      0.5,0.5,
      paste0("Experiments: ",expname1," VS ",expname2,
             "   |   Type: RH (",rh_height_label(logical_name),")",
             region_title,"   |   Period: ",sDate,"-",eDate),
      col="grey50",cex=6,family="Century Gothic"
    )
    abline(h=c(0,1),col="grey50",lwd=3)

    par(mai=c(0,0,0,0)); plot.new()
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0(expname1," (",exptype1,")"),col="grey20",cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0(expname2," (",exptype2,")"),col="grey20",cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new()
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0(expname2," - ",expname1),col="grey20",cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new()

    par(mai=c(0,0,0,0))
    plot.new()
    text(0.5,0.5,paste0("Relative Humidity (",rh_height_label(logical_name),")"),
         col="grey20",cex=map_ylabel_cex,family="Century Gothic",srt=90)

    MapNC(
      filename_topo="",figure_box=figure_box,field_show_box=field_show_box,
      coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,
      projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,
      drawMapBox=regional_mode,
      field_value=field_plot1,field_lon=field_lon,field_lat=field_lat,
      field_pallete_name="TROPOMI_NEW",field_breaks=field_breaks,field_units="%",
      field_pallete_starting_alpha=100,field_show_legend=FALSE
    )

    MapNC(
      filename_topo="",figure_box=figure_box,field_show_box=field_show_box,
      coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,
      projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,
      drawMapBox=regional_mode,
      field_value=field_plot2,field_lon=field_lon,field_lat=field_lat,
      field_pallete_name="TROPOMI_NEW",field_breaks=field_breaks,field_units="%",
      field_pallete_starting_alpha=100,field_show_legend=TRUE,
      field_legend_mai_right=1.8,field_legend_nlabels=6
    )

    MapNC(
      filename_topo="",figure_box=figure_box,field_show_box=field_show_box,
      coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,
      projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,
      drawMapBox=regional_mode,
      field_value=field_plot2-field_plot1,field_lon=field_lon,field_lat=field_lat,
      field_pallete_name="MNMB",field_breaks=field_breaks_diff,field_units="%",
      field_pallete_starting_alpha=100,field_show_legend=TRUE,
      field_legend_mai_right=1.8,field_legend_nlabels=7
    )

    ### Time series
    par(mai=ts_panel_mai,family="Century Gothic")
    x <- seq_along(plot_tmean_tim)
    plot(x,type="n",axes=FALSE,ann=FALSE,ylim=c(0,100),yaxs="i")
    mtext("Time",side=1,line=ts_xlabel_line,cex=ts_xlabel_cex)
    mtext("Relative humidity (%)",side=2,line=ts_ylabel_line,cex=ts_ylabel_cex)

    IDx_labels <- timeseries_axis_indices(plot_tmean_tim)
    axis(1,at=x[IDx_labels],
         labels=format(plot_tmean_tim[IDx_labels],"%Y-%m-%d"),
         cex.axis=ts_xaxis_cex,line=0,lty=0)
    axis(1,at=x[IDx_labels],labels=FALSE,tck=0.01)
    axis(1,at=x[IDx_labels],labels=FALSE,tck=-0.01)
    axis(2,at=seq(0,100,10),labels=seq(0,100,10),las=1,cex.axis=ts_yaxis_cex)

    box(lwd=2)
    abline(h=seq(0,100,10),lwd=1,col="grey")
    abline(v=x[IDx_labels],lwd=1,col="grey")

    lines(x,plot_tmean_var1,lwd=5,col="blue")
    points(x,plot_tmean_var1,pch=19,cex=1.9,col="blue")
    lines(x,plot_tmean_var2,lwd=5,col="red")
    points(x,plot_tmean_var2,pch=19,cex=1.9,col="red")
    add_mean_square(plot_tmean_var1,"blue")
    add_mean_square(plot_tmean_var2,"red")
    legend("top",legend=c(expname1,expname2),lwd=5,col=c("blue","red"),cex=3)

    ### Diurnal cycle
    par(mai=ts_panel_mai,family="Century Gothic")
    plot(1:8,type="n",axes=FALSE,ann=FALSE,ylim=c(0,100),yaxs="i")
    mtext("Time (3 hourly UTC)",side=1,line=ts_xlabel_line,cex=ts_xlabel_cex)
    mtext("Relative humidity (%)",side=2,line=ts_ylabel_line,cex=ts_ylabel_cex)

    axis(1,at=1:8,labels=c("00","03","06","09","12","15","18","21"),cex.axis=ts_xaxis_cex,line=0,lty=0)
    axis(1,at=1:8,labels=FALSE,tck=0.01)
    axis(1,at=1:8,labels=FALSE,tck=-0.01)
    axis(2,at=seq(0,100,10),labels=seq(0,100,10),las=1,cex.axis=ts_yaxis_cex)

    box(lwd=2)
    abline(h=seq(0,100,10),lwd=1,col="grey")
    abline(v=1:8,lwd=1,col="grey")

    lines(1:8,dhourmean_var1,lwd=5,col="blue")
    points(1:8,dhourmean_var1,pch=19,cex=1.9,col="blue")
    lines(1:8,dhourmean_var2,lwd=5,col="red")
    points(1:8,dhourmean_var2,pch=19,cex=1.9,col="red")
    legend("top",legend=c(expname1,expname2),lwd=5,col=c("blue","red"),cex=3)

    dev.off()

    file_tmp <- paste0(file_out,".tmp.png")
    compress(file_in=file_out,file_out=file_tmp)
    file.rename(file_tmp,file_out)
    message("---> RH figure: ",file_out)
  }
}



############################
### OPTICAL PROPERTY PLOTS ###
############################
if (length(optical_plot_variables) > 0) {
  for (logical_name in optical_plot_variables) {
    message("---> Plotting ",logical_name,if (regional_mode) paste0(" for ",region) else "")
    d1 <- list(); d2 <- list()
    for (d in seq_along(seqDate)) {
      d1[[d]] <- read_optical_property(expname1,seqDate[d],logical_name)
      d2[[d]] <- read_optical_property(expname2,seqDate[d],logical_name)
    }
    nx <- dim(d1[[1]])[1]; ny <- dim(d1[[1]])[2]
    nt <- sum(sapply(d1,function(x) dim(x)[3]))
    data1 <- array(unlist(d1),dim=c(nx,ny,nt)); data2 <- array(unlist(d2),dim=c(nx,ny,nt))
    ll <- read_lon_lat(optics_file(expname1,seqDate[1])); field_lon <- ll$lon; field_lat <- ll$lat
    optical_base <- strip_time_aggregation(logical_name)
    species_ts1 <- species_ts2 <- species_dc1 <- species_dc2 <- NULL

    if (optical_base == "aod550_species") {
      species_names <- get_od_species_suffixes(include_soa=TRUE)
      species_ts1 <- list(); species_ts2 <- list(); species_dc1 <- list(); species_dc2 <- list()
      comp1 <- list(); comp2 <- list()
      for (suffix in species_names) {
        tmp1 <- list(); tmp2 <- list()
        for (d in seq_along(seqDate)) {
          tmp1[[d]] <- read_species_od_components(expname1,seqDate[d],exptype1)[[suffix]]
          tmp2[[d]] <- read_species_od_components(expname2,seqDate[d],exptype2)[[suffix]]
        }
        data_comp1 <- array(unlist(tmp1),dim=c(nx,ny,nt))
        data_comp2 <- array(unlist(tmp2),dim=c(nx,ny,nt))
        if (regional_mode) {
          comp1[[suffix]] <- regional_mean(data_comp1,field_lon,field_lat,region_box)
          comp2[[suffix]] <- regional_mean(data_comp2,field_lon,field_lat,region_box)
        } else {
          comp1[[suffix]] <- apply(data_comp1,3,mean,na.rm=TRUE)
          comp2[[suffix]] <- apply(data_comp2,3,mean,na.rm=TRUE)
        }
      }
      species_ts1 <- comp1
      species_ts2 <- comp2
    }

    if (optical_base %in% c("ae550to865","ssa550","mec550")) {
      ### Representative period map: derive from period-mean physical
      ### components instead of averaging the instantaneous ratio.
      f1 <- read_period_optical_map(expname1,seqDate,logical_name)
      f2 <- read_period_optical_map(expname2,seqDate,logical_name)
    } else {
      f1 <- apply(data1,c(1,2),mean,na.rm=TRUE)
      f2 <- apply(data2,c(1,2),mean,na.rm=TRUE)
    }

    if (regional_mode) {
      fp1 <- mask_region_field(f1,field_lon,field_lat,region_box); fp2 <- mask_region_field(f2,field_lon,field_lat,region_box)
      ts1 <- regional_mean(data1,field_lon,field_lat,region_box); ts2 <- regional_mean(data2,field_lon,field_lat,region_box)
    } else { fp1 <- f1; fp2 <- f2; ts1 <- apply(data1,3,mean,na.rm=TRUE); ts2 <- apply(data2,3,mean,na.rm=TRUE) }
    tt <- seq.POSIXt(from=as.POSIXct(paste0(substr(sDate,1,4),"-",substr(sDate,5,6),"-",substr(sDate,7,8)," 00:00:00"),tz="UTC"),by="3 hours",length.out=nt)
    hr <- as.integer(format(tt,"%H")); hrs <- c(0,3,6,9,12,15,18,21)
    dc1 <- sapply(hrs,function(h) mean(ts1[hr==h],na.rm=TRUE)); dc2 <- sapply(hrs,function(h) mean(ts2[hr==h],na.rm=TRUE))
    ts_plot <- aggregate_pair_for_plot(ts1,ts2,tt,logical_name)
    tt <- ts_plot$time
    ts1 <- ts_plot$value1
    ts2 <- ts_plot$value2
    br <- positive_breaks(c(fp1,fp2),200); brd <- difference_breaks(fp2-fp1,200)
    units <- if (strip_time_aggregation(logical_name)=="mec550") "m2 g^-1" else " "
    plot_dir <- make_plot_dir("optics")
    file_out <- paste0(plot_dir,"Optics_",logical_name,"_",expname1,"-",expname2,region_file_tag,"_",sDate,"-",eDate,".png")
    dpi <- 300; png(file_out,width=(0.2+3*3.9+0.8+0.8)*dpi,height=(0.23+0.15+2+2.5)*dpi)
    layout(mat=matrix(c(1,1,1,1,1,1,2:13,14,14,14,14,15,15),4,6,byrow=TRUE),widths=c(0.2,3.9,3.9,0.8,3.9,0.8),heights=c(0.23,0.15,2,2.5))
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0("Experiments: ",expname1," VS ",expname2,"   |   Type: ",get_optical_definition(logical_name)$title,region_title,"   |   Period: ",sDate,"-",eDate),col="grey50",cex=6,family="Century Gothic"); abline(h=c(0,1),col="grey50",lwd=3)
    par(mai=c(0,0,0,0)); plot.new(); par(mai=c(0,0,0,0)); plot.new(); text(.5,.5,paste0(expname1," (",exptype1,")"),cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new(); text(.5,.5,paste0(expname2," (",exptype2,")"),cex=4.5,family="Century Gothic"); par(mai=c(0,0,0,0)); plot.new(); par(mai=c(0,0,0,0)); plot.new(); text(.5,.5,paste0(expname2," - ",expname1),cex=4.5,family="Century Gothic"); par(mai=c(0,0,0,0)); plot.new(); par(mai=c(0,0,0,0)); plot.new(); text(.5,.5,get_optical_definition(logical_name)$title,cex=map_ylabel_cex,family="Century Gothic",srt=90)
    MapNC(filename_topo="",figure_box=figure_box,field_show_box=field_show_box,coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,drawMapBox=regional_mode,field_value=fp1,field_lon=field_lon,field_lat=field_lat,field_pallete_name="TROPOMI_NEW",field_breaks=br,field_units=units,field_pallete_starting_alpha=100,field_show_legend=FALSE)
    MapNC(filename_topo="",figure_box=figure_box,field_show_box=field_show_box,coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,drawMapBox=regional_mode,field_value=fp2,field_lon=field_lon,field_lat=field_lat,field_pallete_name="TROPOMI_NEW",field_breaks=br,field_units=units,field_pallete_starting_alpha=100,field_show_legend=TRUE,field_legend_mai_right=1.8,field_legend_nlabels=7)
    MapNC(filename_topo="",figure_box=figure_box,field_show_box=field_show_box,coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,drawMapBox=regional_mode,field_value=fp2-fp1,field_lon=field_lon,field_lat=field_lat,field_pallete_name="MNMB",field_breaks=brd,field_units=units,field_pallete_starting_alpha=100,field_show_legend=TRUE,field_legend_mai_right=1.8,field_legend_nlabels=7)
    par(mai=ts_panel_mai,family="Century Gothic")
    x <- seq_along(tt)
    species_cols <- c(total="black",ss="#1f77b4",du="#c49c0b",pom="#2ca02c",bc="#9467bd",so4="#d62728",ni="#ff7f0e",am="#17becf",soa="#e377c2")
    if (optical_base == "aod550_species") {
      species_names <- names(species_ts1)
      agg1 <- lapply(species_names,function(s) aggregate_time_series(species_ts1[[s]], seq.POSIXt(from=as.POSIXct(paste0(substr(sDate,1,4),"-",substr(sDate,5,6),"-",substr(sDate,7,8)," 00:00:00"),tz="UTC"),by="3 hours",length.out=length(species_ts1[[s]])), get_time_aggregation(logical_name)[1])$value)
      agg2 <- lapply(species_names,function(s) aggregate_time_series(species_ts2[[s]], seq.POSIXt(from=as.POSIXct(paste0(substr(sDate,1,4),"-",substr(sDate,5,6),"-",substr(sDate,7,8)," 00:00:00"),tz="UTC"),by="3 hours",length.out=length(species_ts2[[s]])), get_time_aggregation(logical_name)[1])$value)
      names(agg1) <- species_names; names(agg2) <- species_names
      yt <- axis_ticks(c(ts1,ts2,unlist(agg1),unlist(agg2)),10)
      plot(x,type="n",axes=FALSE,ann=FALSE,ylim=range(yt$breaks),yaxs="i")
      mtext("Time",1,ts_xlabel_line,cex=ts_xlabel_cex)
      mtext("AOD550 (total + species)",2,ts_ylabel_line,cex=ts_ylabel_cex)
      idx <- timeseries_axis_indices(tt)
      axis(1,at=x[idx],labels=format(tt[idx],"%Y-%m-%d"),cex.axis=ts_xaxis_cex,line=0,lty=0)
      axis(2,at=yt$breaks,labels=yt$labels,las=1,cex.axis=ts_yaxis_cex)
      box(); abline(h=yt$breaks,col="grey")
      lines(x,ts1,lwd=5,col=species_cols["total"],lty=1); points(x,ts1,pch=19,cex=1.8,col=species_cols["total"])
      lines(x,ts2,lwd=5,col=species_cols["total"],lty=2); points(x,ts2,pch=19,cex=1.8,col=species_cols["total"])
      for (s in species_names) {
        lines(x,agg1[[s]],lwd=3,col=species_cols[s],lty=1)
        lines(x,agg2[[s]],lwd=3,col=species_cols[s],lty=2)
      }
      add_mean_square(ts1,species_cols["total"])
      add_mean_square(ts2,"grey40")
      legend("topright",legend=c(expname1,expname2),lwd=4,lty=c(1,2),col="black",cex=2.3,bty="n")
      legend("topleft",legend=c("total",species_names),lwd=c(5,rep(3,length(species_names))),col=species_cols[c("total",species_names)],cex=1.9,bty="n",ncol=3)
      for (s in species_names) {
        species_dc1[[s]] <- sapply(hrs,function(h) mean(species_ts1[[s]][hr==h],na.rm=TRUE))
        species_dc2[[s]] <- sapply(hrs,function(h) mean(species_ts2[[s]][hr==h],na.rm=TRUE))
      }
      yd <- axis_ticks(c(dc1,dc2,unlist(species_dc1),unlist(species_dc2)),10)
      par(mai=ts_panel_mai,family="Century Gothic")
      plot(1:8,type="n",axes=FALSE,ann=FALSE,ylim=range(yd$breaks),yaxs="i")
      mtext("Time (3 hourly UTC)",1,ts_xlabel_line,cex=ts_xlabel_cex)
      mtext("AOD550 (total + species)",2,ts_ylabel_line,cex=ts_ylabel_cex)
      axis(1,at=1:8,labels=sprintf("%02d",hrs),cex.axis=ts_xaxis_cex,line=0,lty=0)
      axis(2,at=yd$breaks,labels=yd$labels,las=1,cex.axis=ts_yaxis_cex)
      box(); abline(h=yd$breaks,col="grey")
      lines(1:8,dc1,lwd=5,col=species_cols["total"],lty=1); points(1:8,dc1,pch=19,cex=1.8,col=species_cols["total"])
      lines(1:8,dc2,lwd=5,col=species_cols["total"],lty=2); points(1:8,dc2,pch=19,cex=1.8,col=species_cols["total"])
      for (s in species_names) {
        lines(1:8,species_dc1[[s]],lwd=3,col=species_cols[s],lty=1)
        lines(1:8,species_dc2[[s]],lwd=3,col=species_cols[s],lty=2)
      }
      legend("topright",legend=c(expname1,expname2),lwd=4,lty=c(1,2),col="black",cex=2.3,bty="n")
    } else {
      yt <- axis_ticks(c(ts1,ts2),10)
      plot(x,type="n",axes=FALSE,ann=FALSE,ylim=range(yt$breaks),yaxs="i")
      mtext("Time",1,ts_xlabel_line,cex=ts_xlabel_cex)
      mtext(get_optical_definition(logical_name)$title,2,ts_ylabel_line,cex=ts_ylabel_cex)
      idx <- timeseries_axis_indices(tt)
      axis(1,at=x[idx],labels=format(tt[idx],"%Y-%m-%d"),cex.axis=ts_xaxis_cex,line=0,lty=0)
      axis(2,at=yt$breaks,labels=yt$labels,las=1,cex.axis=ts_yaxis_cex)
      box(); abline(h=yt$breaks,col="grey")
      lines(x,ts1,lwd=5,col="blue")
      points(x,ts1,pch=19,cex=1.8,col="blue")
      lines(x,ts2,lwd=5,col="red")
      points(x,ts2,pch=19,cex=1.8,col="red")
      add_mean_square(ts1,"blue")
      add_mean_square(ts2,"red")
      legend("top",c(expname1,expname2),lwd=4,col=c("blue","red"),cex=1.7,bty="n")
      par(mai=ts_panel_mai,family="Century Gothic"); yd <- axis_ticks(c(dc1,dc2),10); plot(1:8,type="n",axes=FALSE,ann=FALSE,ylim=range(yd$breaks),yaxs="i"); mtext("Time (3 hourly UTC)",1,ts_xlabel_line,cex=ts_xlabel_cex); mtext(get_optical_definition(logical_name)$title,2,ts_ylabel_line,cex=ts_ylabel_cex); axis(1,at=1:8,labels=sprintf("%02d",hrs),cex.axis=ts_xaxis_cex,line=0,lty=0); axis(2,at=yd$breaks,labels=yd$labels,las=1,cex.axis=ts_yaxis_cex); box(); abline(h=yd$breaks,col="grey"); lines(1:8,dc1,lwd=5,col="blue"); points(1:8,dc1,pch=19,cex=1.8,col="blue"); lines(1:8,dc2,lwd=5,col="red"); points(1:8,dc2,pch=19,cex=1.8,col="red"); legend("top",c(expname1,expname2),lwd=4,col=c("blue","red"),cex=1.7,bty="n")
    }
    dev.off(); tmp <- paste0(file_out,".tmp.png"); compress(file_in=file_out,file_out=tmp); file.rename(tmp,file_out)
  }
}



############################
### MULTI-PANEL SPECIES MAPS ###
############################
if (length(panel_species_plot_variables) > 0) {
  for (logical_name in panel_species_plot_variables) {
    panel_base <- tolower(strip_time_aggregation(logical_name))
    message("---> Plotting ",logical_name,if (regional_mode) paste0(" for ",region) else "")

    ### The optical validation file is guaranteed by 02.start.R for every
    ### multi-panel figure, so it is a stable source of the model grid.
    ref_file <- optics_file(expname1,seqDate[1])
    if (!file.exists(ref_file))
      stop("Optics reference file for multi-panel species plot not found: ",ref_file)
    ll <- read_lon_lat(ref_file)
    field_lon <- ll$lon
    field_lat <- ll$lat

    maps1 <- build_panel_species_maps(expname1,exptype1,seqDate,panel_base)
    maps2 <- build_panel_species_maps(expname2,exptype2,seqDate,panel_base)

    ### Validation is experiment-specific because ME/MAE use that experiment,
    ### while the satellite ensemble itself is identical in both figures.
    validation1 <- build_panel_satellite_validation(expname1,seqDate,panel_validation_variable)
    validation2 <- build_panel_satellite_validation(expname2,seqDate,panel_validation_variable)

    units <- ""
    if (panel_base == "mass_per_species") units <- "kg m^-2"
    if (panel_base == "mec_per_species") units <- "m2 g^-1"
    if (panel_base %in% c("aodratio_per_species","massratio_per_species")) units <- "%"

    plot_species_map_grid(maps1,validation1,expname1,exptype1,panel_base,units=units)
    plot_species_map_grid(maps2,validation2,expname2,exptype2,panel_base,units=units)
  }
}

############################
### WATER-AOD PLOTS ###
############################
if (length(wat_plot_variables) > 0) {
  if (exptype1 != "HAM" || exptype2 != "HAM") stop("Water-AOD variables are currently supported only for HAM/HAM comparisons.")
  for (logical_name in wat_plot_variables) {
    message("---> Plotting ",logical_name,if (regional_mode) paste0(" for ",region) else "")
    d1 <- list(); d2 <- list()
    for (d in seq_along(seqDate)) {
      d1[[d]] <- read_water_aod(expname1,seqDate[d],logical_name)
      d2[[d]] <- read_water_aod(expname2,seqDate[d],logical_name)
    }
    nx <- dim(d1[[1]])[1]; ny <- dim(d1[[1]])[2]
    nt <- sum(sapply(d1,function(x) dim(x)[3]))
    data1 <- array(unlist(d1),dim=c(nx,ny,nt)); data2 <- array(unlist(d2),dim=c(nx,ny,nt))
    ll <- read_lon_lat(wat_ml_file(expname1,seqDate[1])); field_lon <- ll$lon; field_lat <- ll$lat
    f1 <- apply(data1,c(1,2),mean,na.rm=TRUE); f2 <- apply(data2,c(1,2),mean,na.rm=TRUE)
    if (regional_mode) {
      fp1 <- mask_region_field(f1,field_lon,field_lat,region_box); fp2 <- mask_region_field(f2,field_lon,field_lat,region_box)
      ts1 <- regional_mean(data1,field_lon,field_lat,region_box); ts2 <- regional_mean(data2,field_lon,field_lat,region_box)
    } else { fp1 <- f1; fp2 <- f2; ts1 <- apply(data1,3,mean,na.rm=TRUE); ts2 <- apply(data2,3,mean,na.rm=TRUE) }
    tt_native <- seq.POSIXt(from=as.POSIXct(paste0(substr(sDate,1,4),"-",substr(sDate,5,6),"-",substr(sDate,7,8)," 00:00:00"),tz="UTC"),by="3 hours",length.out=nt)
    hr <- as.integer(format(tt_native,"%H")); hrs <- c(0,3,6,9,12,15,18,21)
    dc1 <- sapply(hrs,function(h) mean(ts1[hr==h],na.rm=TRUE)); dc2 <- sapply(hrs,function(h) mean(ts2[hr==h],na.rm=TRUE))
    ts_plot <- aggregate_pair_for_plot(ts1,ts2,tt_native,logical_name)
    tt <- ts_plot$time; ts1 <- ts_plot$value1; ts2 <- ts_plot$value2
    br <- positive_breaks(c(fp1,fp2),200); brd <- difference_breaks(fp2-fp1,200)
    plot_dir <- make_plot_dir("optics")
    file_out <- paste0(plot_dir,"WAT_",logical_name,"_",expname1,"-",expname2,region_file_tag,"_",sDate,"-",eDate,".png")
    dpi <- 300; png(file_out,width=(0.2+3*3.9+0.8+0.8)*dpi,height=(0.23+0.15+2+2.5)*dpi)
    layout(mat=matrix(c(1,1,1,1,1,1,2:13,14,14,14,14,15,15),4,6,byrow=TRUE),widths=c(0.2,3.9,3.9,0.8,3.9,0.8),heights=c(0.23,0.15,2,2.5))
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0("Experiments: ",expname1," VS ",expname2,"   |   Type: ",get_wat_definition(logical_name)$title,region_title,"   |   Period: ",sDate,"-",eDate),col="grey50",cex=6,family="Century Gothic"); abline(h=c(0,1),col="grey50",lwd=3)
    par(mai=c(0,0,0,0)); plot.new(); par(mai=c(0,0,0,0)); plot.new(); text(.5,.5,paste0(expname1," (",exptype1,")"),cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new(); text(.5,.5,paste0(expname2," (",exptype2,")"),cex=4.5,family="Century Gothic"); par(mai=c(0,0,0,0)); plot.new(); par(mai=c(0,0,0,0)); plot.new(); text(.5,.5,paste0(expname2," - ",expname1),cex=4.5,family="Century Gothic"); par(mai=c(0,0,0,0)); plot.new(); par(mai=c(0,0,0,0)); plot.new(); text(.5,.5,get_wat_definition(logical_name)$title,cex=map_ylabel_cex,family="Century Gothic",srt=90)
    MapNC(filename_topo="",figure_box=figure_box,field_show_box=field_show_box,coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,drawMapBox=regional_mode,field_value=fp1,field_lon=field_lon,field_lat=field_lat,field_pallete_name="TROPOMI_NEW",field_breaks=br,field_units=" ",field_pallete_starting_alpha=100,field_show_legend=FALSE)
    MapNC(filename_topo="",figure_box=figure_box,field_show_box=field_show_box,coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,drawMapBox=regional_mode,field_value=fp2,field_lon=field_lon,field_lat=field_lat,field_pallete_name="TROPOMI_NEW",field_breaks=br,field_units=" ",field_pallete_starting_alpha=100,field_show_legend=TRUE,field_legend_mai_right=1.8,field_legend_nlabels=7)
    MapNC(filename_topo="",figure_box=figure_box,field_show_box=field_show_box,coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,drawMapBox=regional_mode,field_value=fp2-fp1,field_lon=field_lon,field_lat=field_lat,field_pallete_name="MNMB",field_breaks=brd,field_units=" ",field_pallete_starting_alpha=100,field_show_legend=TRUE,field_legend_mai_right=1.8,field_legend_nlabels=7)
    par(mai=ts_panel_mai,family="Century Gothic")
    x <- seq_along(tt); yt <- axis_ticks(c(ts1,ts2),10)
    plot(x,type="n",axes=FALSE,ann=FALSE,ylim=range(yt$breaks),yaxs="i"); mtext("Time",1,ts_xlabel_line,cex=ts_xlabel_cex); mtext(get_wat_definition(logical_name)$title,2,ts_ylabel_line,cex=ts_ylabel_cex)
    idx <- timeseries_axis_indices(tt); axis(1,at=x[idx],labels=format(tt[idx],"%Y-%m-%d"),cex.axis=ts_xaxis_cex,line=0,lty=0); axis(2,at=yt$breaks,labels=yt$labels,las=1,cex.axis=ts_yaxis_cex); box(); abline(h=yt$breaks,col="grey")
    lines(x,ts1,lwd=5,col="blue"); points(x,ts1,pch=19,cex=1.8,col="blue"); lines(x,ts2,lwd=5,col="red"); points(x,ts2,pch=19,cex=1.8,col="red"); add_mean_square(ts1,"blue"); add_mean_square(ts2,"red"); legend("top",c(expname1,expname2),lwd=4,col=c("blue","red"),cex=1.7,bty="n")
    par(mai=ts_panel_mai,family="Century Gothic"); yd <- axis_ticks(c(dc1,dc2),10); plot(1:8,type="n",axes=FALSE,ann=FALSE,ylim=range(yd$breaks),yaxs="i"); mtext("Time (3 hourly UTC)",1,ts_xlabel_line,cex=ts_xlabel_cex); mtext(get_wat_definition(logical_name)$title,2,ts_ylabel_line,cex=ts_ylabel_cex); axis(1,at=1:8,labels=sprintf("%02d",hrs),cex.axis=ts_xaxis_cex,line=0,lty=0); axis(2,at=yd$breaks,labels=yd$labels,las=1,cex.axis=ts_yaxis_cex); box(); abline(h=yd$breaks,col="grey"); lines(1:8,dc1,lwd=5,col="blue"); points(1:8,dc1,pch=19,cex=1.8,col="blue"); lines(1:8,dc2,lwd=5,col="red"); points(1:8,dc2,pch=19,cex=1.8,col="red"); legend("top",c(expname1,expname2),lwd=4,col=c("blue","red"),cex=1.7,bty="n")
    dev.off(); tmp <- paste0(file_out,".tmp.png"); compress(file_in=file_out,file_out=tmp); file.rename(tmp,file_out)
  }
}

############################
### TOTAL-COLUMN GAS PLOTS ###
############################
if (length(column_plot_variables) > 0) {

  for (logical_name in column_plot_variables) {

    def <- get_column_definition(logical_name)

    message("---> Plotting ",logical_name,
            if (regional_mode) paste0(" for ",region) else "")

    d1 <- list()
    d2 <- list()

    for (d in seq_along(seqDate)) {
      d1[[d]] <- read_total_column(expname1,seqDate[d],logical_name)
      d2[[d]] <- read_total_column(expname2,seqDate[d],logical_name)
    }

    nx1 <- dim(d1[[1]])[1]; ny1 <- dim(d1[[1]])[2]
    nx2 <- dim(d2[[1]])[1]; ny2 <- dim(d2[[1]])[2]

    if (nx1 != nx2 || ny1 != ny2)
      stop("Spatial dimensions differ between experiments for ",logical_name)

    nt1 <- sum(sapply(d1,function(x) dim(x)[3]))
    nt2 <- sum(sapply(d2,function(x) dim(x)[3]))

    if (nt1 != nt2)
      stop("Time dimensions differ between experiments for ",logical_name)

    data1 <- array(unlist(d1),dim=c(nx1,ny1,nt1))
    data2 <- array(unlist(d2),dim=c(nx2,ny2,nt2))

    ll <- read_lon_lat(column_file(expname1,seqDate[1]))
    field_lon <- ll$lon
    field_lat <- ll$lat

    f1 <- apply(data1,c(1,2),mean,na.rm=TRUE)
    f2 <- apply(data2,c(1,2),mean,na.rm=TRUE)

    if (regional_mode) {
      fp1 <- mask_region_field(f1,field_lon,field_lat,region_box)
      fp2 <- mask_region_field(f2,field_lon,field_lat,region_box)
      ts1 <- regional_mean(data1,field_lon,field_lat,region_box)
      ts2 <- regional_mean(data2,field_lon,field_lat,region_box)
    } else {
      fp1 <- f1
      fp2 <- f2
      ts1 <- apply(data1,3,mean,na.rm=TRUE)
      ts2 <- apply(data2,3,mean,na.rm=TRUE)
    }

    nt <- dim(data1)[3]
    tt <- seq.POSIXt(
      from=as.POSIXct(paste0(substr(sDate,1,4),"-",substr(sDate,5,6),"-",
                             substr(sDate,7,8)," 00:00:00"),tz="UTC"),
      by="3 hours",length.out=nt
    )

    hr <- as.integer(format(tt,"%H"))
    hrs <- c(0,3,6,9,12,15,18,21)
    dc1 <- sapply(hrs,function(h) mean(ts1[hr == h],na.rm=TRUE))
    dc2 <- sapply(hrs,function(h) mean(ts2[hr == h],na.rm=TRUE))

    ts_plot <- aggregate_pair_for_plot(ts1,ts2,tt,logical_name)
    tt <- ts_plot$time
    ts1 <- ts_plot$value1
    ts2 <- ts_plot$value2

    br <- positive_breaks(c(fp1,fp2),200)
    brd <- difference_breaks(fp2-fp1,200)

    plot_dir <- make_plot_dir("per_species")
    file_out <- paste0(
      plot_dir,"GasColumn_",logical_name,"_",expname1,"-",expname2,
      region_file_tag,"_",sDate,"-",eDate,".png"
    )

    dpi <- 300
    png(file_out,width=(0.2+3*3.9+0.8+0.8)*dpi,
        height=(0.23+0.15+2+2.5)*dpi)

    layout(
      mat=matrix(c(1,1,1,1,1,1,2:13,14,14,14,14,15,15),4,6,byrow=TRUE),
      widths=c(0.2,3.9,3.9,0.8,3.9,0.8),
      heights=c(0.23,0.15,2,2.5)
    )

    par(mai=c(0,0,0,0))
    plot.new()
    text(
      0.5,0.5,
      paste0("Experiments: ",expname1," VS ",expname2,
             "   |   Type: ",def$title,
             region_title,"   |   Period: ",sDate,"-",eDate),
      col="grey50",cex=6,family="Century Gothic"
    )
    abline(h=c(0,1),col="grey50",lwd=3)

    par(mai=c(0,0,0,0)); plot.new()
    par(mai=c(0,0,0,0)); plot.new()
    text(.5,.5,paste0(expname1," (",exptype1,")"),cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new()
    text(.5,.5,paste0(expname2," (",exptype2,")"),cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new()
    par(mai=c(0,0,0,0)); plot.new()
    text(.5,.5,paste0(expname2," - ",expname1),cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new()
    par(mai=c(0,0,0,0)); plot.new()
    text(.5,.5,def$title,cex=map_ylabel_cex,family="Century Gothic",srt=90)

    MapNC(
      filename_topo="",figure_box=figure_box,field_show_box=field_show_box,
      coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,
      projection=projection,lonmax=lonmax,lonmin=lonmin,
      latmax=latmax,latmin=latmin,drawMapBox=regional_mode,
      field_value=fp1,field_lon=field_lon,field_lat=field_lat,
      field_pallete_name="TROPOMI_NEW",field_breaks=br,
      field_units=def$units,field_pallete_starting_alpha=100,
      field_show_legend=FALSE
    )

    MapNC(
      filename_topo="",figure_box=figure_box,field_show_box=field_show_box,
      coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,
      projection=projection,lonmax=lonmax,lonmin=lonmin,
      latmax=latmax,latmin=latmin,drawMapBox=regional_mode,
      field_value=fp2,field_lon=field_lon,field_lat=field_lat,
      field_pallete_name="TROPOMI_NEW",field_breaks=br,
      field_units=def$units,field_pallete_starting_alpha=100,
      field_show_legend=TRUE,field_legend_mai_right=1.8,field_legend_nlabels=7
    )

    MapNC(
      filename_topo="",figure_box=figure_box,field_show_box=field_show_box,
      coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,
      projection=projection,lonmax=lonmax,lonmin=lonmin,
      latmax=latmax,latmin=latmin,drawMapBox=regional_mode,
      field_value=fp2-fp1,field_lon=field_lon,field_lat=field_lat,
      field_pallete_name="MNMB",field_breaks=brd,
      field_units=def$units,field_pallete_starting_alpha=100,
      field_show_legend=TRUE,field_legend_mai_right=1.8,field_legend_nlabels=7
    )

    ### Time series
    par(mai=ts_panel_mai,family="Century Gothic")
    x <- seq_along(tt)
    yt <- axis_ticks(c(ts1,ts2),10)

    plot(x,type="n",axes=FALSE,ann=FALSE,
         ylim=range(yt$breaks),yaxs="i")

    mtext("Time",1,ts_xlabel_line,cex=ts_xlabel_cex)
    mtext(paste0(def$title," (",def$units,")"),2,ts_ylabel_line,cex=ts_ylabel_cex)

    idx <- timeseries_axis_indices(tt)

    axis(1,at=x[idx],labels=format(tt[idx],"%Y-%m-%d"),
         cex.axis=ts_xaxis_cex,line=0,lty=0)
    axis(2,at=yt$breaks,labels=yt$labels,las=1,cex.axis=ts_yaxis_cex)

    box()
    abline(h=yt$breaks,col="grey")
    abline(v=x[idx],col="grey")

    lines(x,ts1,lwd=5,col="blue")
    points(x,ts1,pch=19,cex=1.8,col="blue")
    lines(x,ts2,lwd=5,col="red")
    points(x,ts2,pch=19,cex=1.8,col="red")
    add_mean_square(ts1,"blue")
    add_mean_square(ts2,"red")

    legend("top",c(expname1,expname2),
           lwd=5,col=c("blue","red"),cex=3)

    ### Diurnal cycle
    par(mai=ts_panel_mai,family="Century Gothic")
    yd <- axis_ticks(c(dc1,dc2),10)

    plot(1:8,type="n",axes=FALSE,ann=FALSE,
         ylim=range(yd$breaks),yaxs="i")

    mtext("Time (3 hourly UTC)",1,ts_xlabel_line,cex=ts_xlabel_cex)
    mtext(paste0(def$title," (",def$units,")"),2,ts_ylabel_line,cex=ts_ylabel_cex)

    axis(1,at=1:8,labels=sprintf("%02d",hrs),
         cex.axis=ts_xaxis_cex,line=0,lty=0)
    axis(2,at=yd$breaks,labels=yd$labels,las=1,cex.axis=ts_yaxis_cex)

    box()
    abline(h=yd$breaks,col="grey")
    abline(v=1:8,col="grey")

    lines(1:8,dc1,lwd=5,col="blue")
    points(1:8,dc1,pch=19,cex=1.8,col="blue")
    lines(1:8,dc2,lwd=5,col="red")
    points(1:8,dc2,pch=19,cex=1.8,col="red")

    legend("top",c(expname1,expname2),
           lwd=5,col=c("blue","red"),cex=3)

    dev.off()

    tmp <- paste0(file_out,".tmp.png")
    compress(file_in=file_out,file_out=tmp)
    file.rename(tmp,file_out)

    message("---> Total-column figure: ",file_out)
  }
}

#########################
### PRECIPITATION PLOTS ###
#########################
if (length(precip_variables) > 0) {

  for (logical_name in precip_variables) {

    def <- get_precip_definition(logical_name)

    message("---> Plotting ",logical_name,
            if (regional_mode) paste0(" for ",region) else "")

    field_day1 <- list()
    field_day2 <- list()

    for (d in seq_along(seqDate)) {
      field_day1[[d]] <- read_precipitation(expname1,seqDate[d],logical_name)
      field_day2[[d]] <- read_precipitation(expname2,seqDate[d],logical_name)
    }

    nx1 <- dim(field_day1[[1]])[1]; ny1 <- dim(field_day1[[1]])[2]
    nx2 <- dim(field_day2[[1]])[1]; ny2 <- dim(field_day2[[1]])[2]

    if (nx1 != nx2 || ny1 != ny2)
      stop("Spatial dimensions differ between experiments for ",logical_name)

    nt1 <- sum(sapply(field_day1,function(x) dim(x)[3]))
    nt2 <- sum(sapply(field_day2,function(x) dim(x)[3]))

    if (nt1 != nt2)
      stop("Time dimensions differ between experiments for ",logical_name)

    data1 <- array(unlist(field_day1),dim=c(nx1,ny1,nt1))
    data2 <- array(unlist(field_day2),dim=c(nx2,ny2,nt2))

    file1 <- precip_file(expname1,seqDate[1])
    ll <- read_lon_lat(file1)
    field_lon <- ll$lon
    field_lat <- ll$lat

    field_var1 <- apply(data1,c(1,2),mean,na.rm=TRUE)
    field_var2 <- apply(data2,c(1,2),mean,na.rm=TRUE)

    if (regional_mode) {
      field_plot1 <- mask_region_field(field_var1,field_lon,field_lat,region_box)
      field_plot2 <- mask_region_field(field_var2,field_lon,field_lat,region_box)
      tmean_var1 <- regional_mean(data1,field_lon,field_lat,region_box)
      tmean_var2 <- regional_mean(data2,field_lon,field_lat,region_box)
    } else {
      field_plot1 <- field_var1
      field_plot2 <- field_var2
      tmean_var1 <- apply(data1,3,mean,na.rm=TRUE)
      tmean_var2 <- apply(data2,3,mean,na.rm=TRUE)
    }

    ### The eight values per forecast correspond to the intervals:
    ### 00-03, 03-06, ..., 21-24 UTC. Associate each amount with the
    ### beginning of its 3-hour interval to retain COMPASS's 00...21 axis.
    nt <- dim(data1)[3]
    tmean_tim <- seq.POSIXt(
      from=as.POSIXct(paste0(substr(sDate,1,4),"-",substr(sDate,5,6),"-",substr(sDate,7,8)," 00:00:00"),tz="UTC"),
      by="3 hours",length.out=nt
    )

    hour <- as.integer(format(tmean_tim,"%H"))
    hours <- c(0,3,6,9,12,15,18,21)

    dhourmean_var1 <- sapply(hours,function(h) mean(tmean_var1[hour == h],na.rm=TRUE))
    dhourmean_var2 <- sapply(hours,function(h) mean(tmean_var2[hour == h],na.rm=TRUE))

    ts_plot <- aggregate_pair_for_plot(tmean_var1,tmean_var2,tmean_tim,logical_name)
    plot_tmean_tim <- ts_plot$time
    plot_tmean_var1 <- ts_plot$value1
    plot_tmean_var2 <- ts_plot$value2

    field_breaks <- positive_breaks(c(field_plot1,field_plot2),ncolors=200)
    field_breaks_diff <- difference_breaks(field_plot2-field_plot1,ncolors=200)

    plot_dir <- make_plot_dir("meteorology")

    file_out <- paste0(
      plot_dir,"Meteorology_",logical_name,"_vs_",logical_name,"_",
      expname1,"-",expname2,region_file_tag,"_",sDate,"-",eDate,".png"
    )

    dpi <- 300
    png(file_out,width=(0.2+3*3.9+0.8+0.8)*dpi,height=(0.23+0.15+2+2.5)*dpi)
    layout(
      mat=matrix(c(1,1,1,1,1,1,2:13,14,14,14,14,15,15),4,6,byrow=TRUE),
      widths=c(0.2,3.9,3.9,0.8,3.9,0.8),
      heights=c(0.23,0.15,2,2.5)
    )

    par(mai=c(0,0,0,0))
    plot.new()
    text(
      0.5,0.5,
      paste0("Experiments: ",expname1," VS ",expname2,
             "   |   Type: ",def$title,
             region_title,"   |   Period: ",sDate,"-",eDate),
      col="grey50",cex=6,family="Century Gothic"
    )
    abline(h=c(0,1),col="grey50",lwd=3)

    par(mai=c(0,0,0,0)); plot.new()
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0(expname1," (",exptype1,")"),col="grey20",cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0(expname2," (",exptype2,")"),col="grey20",cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new()
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0(expname2," - ",expname1),col="grey20",cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new()

    par(mai=c(0,0,0,0))
    plot.new()
    text(0.5,0.5,paste0(def$title," (3 h)"),
         col="grey20",cex=map_ylabel_cex,family="Century Gothic",srt=90)

    MapNC(
      filename_topo="",figure_box=figure_box,field_show_box=field_show_box,
      coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,
      projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,
      drawMapBox=regional_mode,
      field_value=field_plot1,field_lon=field_lon,field_lat=field_lat,
      field_pallete_name="TROPOMI_NEW",field_breaks=field_breaks,field_units="mm / 3 h",
      field_pallete_starting_alpha=100,field_show_legend=FALSE
    )

    MapNC(
      filename_topo="",figure_box=figure_box,field_show_box=field_show_box,
      coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,
      projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,
      drawMapBox=regional_mode,
      field_value=field_plot2,field_lon=field_lon,field_lat=field_lat,
      field_pallete_name="TROPOMI_NEW",field_breaks=field_breaks,field_units="mm / 3 h",
      field_pallete_starting_alpha=100,field_show_legend=TRUE,
      field_legend_mai_right=1.8,field_legend_nlabels=7
    )

    MapNC(
      filename_topo="",figure_box=figure_box,field_show_box=field_show_box,
      coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,
      projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,
      drawMapBox=regional_mode,
      field_value=field_plot2-field_plot1,field_lon=field_lon,field_lat=field_lat,
      field_pallete_name="MNMB",field_breaks=field_breaks_diff,field_units="mm / 3 h",
      field_pallete_starting_alpha=100,field_show_legend=TRUE,
      field_legend_mai_right=1.8,field_legend_nlabels=7
    )

    ### Time series
    par(mai=ts_panel_mai,family="Century Gothic")
    x <- seq_along(plot_tmean_tim)
    yseq <- positive_axis_ticks(c(plot_tmean_var1,plot_tmean_var2),n=10)

    plot(x,type="n",axes=FALSE,ann=FALSE,
         ylim=c(0,max(yseq$breaks)),yaxs="i")

    mtext("Time",side=1,line=ts_xlabel_line,cex=ts_xlabel_cex)
    mtext("Precipitation (mm / 3 h)",side=2,line=ts_ylabel_line,cex=ts_ylabel_cex)

    IDx_labels <- timeseries_axis_indices(plot_tmean_tim)

    axis(1,at=x[IDx_labels],
         labels=format(plot_tmean_tim[IDx_labels],"%Y-%m-%d"),
         cex.axis=ts_xaxis_cex,line=0,lty=0)
    axis(1,at=x[IDx_labels],labels=FALSE,tck=0.01)
    axis(1,at=x[IDx_labels],labels=FALSE,tck=-0.01)
    axis(2,at=yseq$breaks,labels=yseq$labels,las=1,cex.axis=ts_yaxis_cex)

    box(lwd=2)
    abline(h=yseq$breaks,lwd=1,col="grey")
    abline(v=x[IDx_labels],lwd=1,col="grey")

    lines(x,plot_tmean_var1,lwd=5,col="blue")
    points(x,plot_tmean_var1,pch=19,cex=1.9,col="blue")
    lines(x,plot_tmean_var2,lwd=5,col="red")
    points(x,plot_tmean_var2,pch=19,cex=1.9,col="red")
    add_mean_square(plot_tmean_var1,"blue")
    add_mean_square(plot_tmean_var2,"red")
    legend("top",legend=c(expname1,expname2),lwd=5,col=c("blue","red"),cex=3)

    ### Diurnal cycle
    par(mai=ts_panel_mai,family="Century Gothic")
    yseq <- positive_axis_ticks(c(dhourmean_var1,dhourmean_var2),n=10)

    plot(1:8,type="n",axes=FALSE,ann=FALSE,
         ylim=c(0,max(yseq$breaks)),yaxs="i")

    mtext("3-hour interval starting UTC",side=1,line=12,cex=3.5)
    mtext("Precipitation (mm / 3 h)",side=2,line=ts_ylabel_line,cex=ts_ylabel_cex)

    axis(1,at=1:8,labels=c("00","03","06","09","12","15","18","21"),cex.axis=ts_xaxis_cex,line=0,lty=0)
    axis(1,at=1:8,labels=FALSE,tck=0.01)
    axis(1,at=1:8,labels=FALSE,tck=-0.01)
    axis(2,at=yseq$breaks,labels=yseq$labels,las=1,cex.axis=ts_yaxis_cex)

    box(lwd=2)
    abline(h=yseq$breaks,lwd=1,col="grey")
    abline(v=1:8,lwd=1,col="grey")

    lines(1:8,dhourmean_var1,lwd=5,col="blue")
    points(1:8,dhourmean_var1,pch=19,cex=1.9,col="blue")
    lines(1:8,dhourmean_var2,lwd=5,col="red")
    points(1:8,dhourmean_var2,pch=19,cex=1.9,col="red")
    legend("top",legend=c(expname1,expname2),lwd=5,col=c("blue","red"),cex=3)

    dev.off()

    file_tmp <- paste0(file_out,".tmp.png")
    compress(file_in=file_out,file_out=file_tmp)
    file.rename(file_tmp,file_out)

    message("---> Precipitation figure: ",file_out)
  }
}


###################################
### DEPOSITION LIFETIME PLOTS ###
###################################
if (length(lifetime_variables) > 0) {

  for (lifetime_name in lifetime_variables) {

    suffix <- sub("^lifetime_","",strip_time_aggregation(lifetime_name))
    mss_name <- paste0("mss_",suffix)

    message("---> Plotting ",lifetime_name,
            if (regional_mode) paste0(" for ",region) else "")

    d1 <- read_lifetime_components(
      expname1,exptype1,suffix,variables_exp1
    )
    d2 <- read_lifetime_components(
      expname2,exptype2,suffix,variables_exp2
    )

    if (d1$nx != d2$nx || d1$ny != d2$ny)
      stop("Spatial dimensions differ between experiments for ",lifetime_name)

    if (d1$nt != d2$nt)
      stop("Time dimensions differ between experiments for ",lifetime_name)

    file1 <- variable_file(mss_name,variables_exp1,expname1,seqDate[1])
    ll <- read_lon_lat(file1)
    field_lon <- ll$lon
    field_lat <- ll$lat

    ### Spatial lifetime uses the period-mean burden divided by the
    ### period-mean total deposition loss at each grid cell.
    field_var1 <- apply(d1$mss,c(1,2),mean,na.rm=TRUE) /
                  (apply(d1$dep,c(1,2),mean,na.rm=TRUE)*86400)
    field_var2 <- apply(d2$mss,c(1,2),mean,na.rm=TRUE) /
                  (apply(d2$dep,c(1,2),mean,na.rm=TRUE)*86400)

    field_var1[!is.finite(field_var1) | field_var1 < 0] <- NA_real_
    field_var2[!is.finite(field_var2) | field_var2 < 0] <- NA_real_

    if (regional_mode) {
      field_plot1 <- mask_region_field(field_var1,field_lon,field_lat,region_box)
      field_plot2 <- mask_region_field(field_var2,field_lon,field_lat,region_box)

      mass1 <- regional_mass_tg(d1$mss,field_lon,field_lat,region_box)
      mass2 <- regional_mass_tg(d2$mss,field_lon,field_lat,region_box)
      dep1 <- regional_flux_tg_day(d1$dep,field_lon,field_lat,region_box)
      dep2 <- regional_flux_tg_day(d2$dep,field_lon,field_lat,region_box)
    } else {
      field_plot1 <- field_var1
      field_plot2 <- field_var2

      mass1 <- global_mass_tg(d1$mss,field_lon,field_lat)
      mass2 <- global_mass_tg(d2$mss,field_lon,field_lat)
      dep1 <- global_flux_tg_day(d1$dep,field_lon,field_lat)
      dep2 <- global_flux_tg_day(d2$dep,field_lon,field_lat)
    }

    ### Regional/global lifetime is ratio of integrated burden to integrated
    ### total deposition loss, giving days because dep1/dep2 are Tg day^-1.
    tmean_var1 <- mass1/dep1
    tmean_var2 <- mass2/dep2
    tmean_var1[!is.finite(tmean_var1) | tmean_var1 < 0] <- NA_real_
    tmean_var2[!is.finite(tmean_var2) | tmean_var2 < 0] <- NA_real_

    tmean_tim <- seq.POSIXt(
      from=as.POSIXct(
        paste0(substr(sDate,1,4),"-",substr(sDate,5,6),"-",
               substr(sDate,7,8)," 00:00:00"),
        tz="UTC"
      ),
      by="3 hours",length.out=d1$nt
    )

    hour <- as.integer(format(tmean_tim,"%H"))
    hours <- c(0,3,6,9,12,15,18,21)

    dhourmean_var1 <- sapply(
      hours,function(h) mean(tmean_var1[hour == h],na.rm=TRUE)
    )
    dhourmean_var2 <- sapply(
      hours,function(h) mean(tmean_var2[hour == h],na.rm=TRUE)
    )

    ts_plot <- aggregate_pair_for_plot(
      tmean_var1,tmean_var2,tmean_tim,lifetime_name
    )
    plot_tmean_tim <- ts_plot$time
    plot_tmean_var1 <- ts_plot$value1
    plot_tmean_var2 <- ts_plot$value2

    field_breaks <- positive_breaks(
      c(field_plot1,field_plot2),ncolors=200
    )
    field_breaks_diff <- difference_breaks(
      field_plot2-field_plot1,ncolors=200
    )

    plot_category <- get_plot_category(mss_name)
    plot_dir <- make_plot_dir(plot_category)

    file_out <- paste0(
      plot_dir,"DepositionLifetime_",lifetime_name,"_",
      expname1,"-",expname2,region_file_tag,"_",sDate,"-",eDate,".png"
    )

    dpi <- 300
    png(
      file_out,
      width=(0.2+3*3.9+0.8+0.8)*dpi,
      height=(0.23+0.15+2+2.5)*dpi
    )

    layout(
      mat=matrix(
        c(1,1,1,1,1,1,2:13,14,14,14,14,15,15),
        4,6,byrow=TRUE
      ),
      widths=c(0.2,3.9,3.9,0.8,3.9,0.8),
      heights=c(0.23,0.15,2,2.5)
    )

    par(mai=c(0,0,0,0))
    plot.new()
    text(
      0.5,0.5,
      paste0(
        "Experiments: ",expname1," VS ",expname2,
        "   |   Type: Deposition lifetime ",toupper(suffix),
        region_title,"   |   Period: ",sDate,"-",eDate
      ),
      col="grey50",cex=6,family="Century Gothic"
    )
    abline(h=c(0,1),col="grey50",lwd=3)

    par(mai=c(0,0,0,0)); plot.new()
    par(mai=c(0,0,0,0)); plot.new()
    text(
      0.5,0.5,paste0(expname1," (",exptype1,")"),
      col="grey20",cex=4.5,family="Century Gothic"
    )
    par(mai=c(0,0,0,0)); plot.new()
    text(
      0.5,0.5,paste0(expname2," (",exptype2,")"),
      col="grey20",cex=4.5,family="Century Gothic"
    )
    par(mai=c(0,0,0,0)); plot.new()
    par(mai=c(0,0,0,0)); plot.new()
    text(
      0.5,0.5,paste0(expname2," - ",expname1),
      col="grey20",cex=4.5,family="Century Gothic"
    )
    par(mai=c(0,0,0,0)); plot.new()
    par(mai=c(0,0,0,0)); plot.new()
    text(
      0.5,0.5,"Deposition lifetime",
      col="grey20",cex=map_ylabel_cex,family="Century Gothic",srt=90
    )

    MapNC(
      filename_topo="",figure_box=figure_box,field_show_box=field_show_box,
      coastlineWorldFine_lwd=coastlineWorldFine_lwd,
      gridlines=gridlines,projection=projection,
      lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,
      drawMapBox=regional_mode,
      field_value=field_plot1,field_lon=field_lon,field_lat=field_lat,
      field_pallete_name="TROPOMI_NEW",
      field_breaks=field_breaks,field_units="days",
      field_pallete_starting_alpha=100,field_show_legend=FALSE
    )

    MapNC(
      filename_topo="",figure_box=figure_box,field_show_box=field_show_box,
      coastlineWorldFine_lwd=coastlineWorldFine_lwd,
      gridlines=gridlines,projection=projection,
      lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,
      drawMapBox=regional_mode,
      field_value=field_plot2,field_lon=field_lon,field_lat=field_lat,
      field_pallete_name="TROPOMI_NEW",
      field_breaks=field_breaks,field_units="days",
      field_pallete_starting_alpha=100,field_show_legend=TRUE,
      field_legend_mai_right=1.8,field_legend_nlabels=7
    )

    MapNC(
      filename_topo="",figure_box=figure_box,field_show_box=field_show_box,
      coastlineWorldFine_lwd=coastlineWorldFine_lwd,
      gridlines=gridlines,projection=projection,
      lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,
      drawMapBox=regional_mode,
      field_value=field_plot2-field_plot1,
      field_lon=field_lon,field_lat=field_lat,
      field_pallete_name="MNMB",
      field_breaks=field_breaks_diff,field_units="days",
      field_pallete_starting_alpha=100,field_show_legend=TRUE,
      field_legend_mai_right=1.8,field_legend_nlabels=7
    )

    ### Time series
    par(mai=ts_panel_mai,family="Century Gothic")
    x <- seq_along(plot_tmean_tim)
    yseq <- positive_axis_ticks(
      c(plot_tmean_var1,plot_tmean_var2),n=10
    )

    plot(
      x,type="n",axes=FALSE,ann=FALSE,
      ylim=c(0,max(yseq$breaks)),yaxs="i"
    )

    mtext("Time",side=1,line=ts_xlabel_line,cex=ts_xlabel_cex)
    mtext("Deposition lifetime (days)",side=2,line=ts_ylabel_line,cex=ts_ylabel_cex)

    IDx_labels <- timeseries_axis_indices(plot_tmean_tim)

    axis(
      1,at=x[IDx_labels],
      labels=format(plot_tmean_tim[IDx_labels],"%Y-%m-%d"),
      cex.axis=ts_xaxis_cex,line=0,lty=0
    )
    axis(1,at=x[IDx_labels],labels=FALSE,tck=0.01)
    axis(1,at=x[IDx_labels],labels=FALSE,tck=-0.01)
    axis(2,at=yseq$breaks,labels=yseq$labels,las=1,cex.axis=ts_yaxis_cex)

    box(lwd=2)
    abline(h=yseq$breaks,lwd=1,col="grey")
    abline(v=x[IDx_labels],lwd=1,col="grey")

    lines(x,plot_tmean_var1,lwd=5,col="blue")
    points(x,plot_tmean_var1,pch=19,cex=1.8,col="blue")
    lines(x,plot_tmean_var2,lwd=5,col="red")
    points(x,plot_tmean_var2,pch=19,cex=1.8,col="red")

    add_mean_square(plot_tmean_var1,"blue")
    add_mean_square(plot_tmean_var2,"red")

    legend(
      "top",legend=c(expname1,expname2),
      lwd=5,col=c("blue","red"),cex=3
    )

    ### Native 3-hourly diurnal cycle
    par(mai=ts_panel_mai,family="Century Gothic")
    yseq <- positive_axis_ticks(
      c(dhourmean_var1,dhourmean_var2),n=10
    )

    plot(
      1:8,type="n",axes=FALSE,ann=FALSE,
      ylim=c(0,max(yseq$breaks)),yaxs="i"
    )

    mtext("Time (3 hourly UTC)",side=1,line=ts_xlabel_line,cex=ts_xlabel_cex)
    mtext("Deposition lifetime (days)",side=2,line=ts_ylabel_line,cex=ts_ylabel_cex)

    axis(
      1,at=1:8,
      labels=c("00","03","06","09","12","15","18","21"),
      cex.axis=ts_xaxis_cex,line=0,lty=0
    )
    axis(1,at=1:8,labels=FALSE,tck=0.01)
    axis(1,at=1:8,labels=FALSE,tck=-0.01)
    axis(2,at=yseq$breaks,labels=yseq$labels,las=1,cex.axis=ts_yaxis_cex)

    box(lwd=2)
    abline(h=yseq$breaks,lwd=1,col="grey")
    abline(v=1:8,lwd=1,col="grey")

    lines(1:8,dhourmean_var1,lwd=5,col="blue")
    points(1:8,dhourmean_var1,pch=19,cex=1.8,col="blue")
    lines(1:8,dhourmean_var2,lwd=5,col="red")
    points(1:8,dhourmean_var2,pch=19,cex=1.8,col="red")

    legend(
      "top",legend=c(expname1,expname2),
      lwd=5,col=c("blue","red"),cex=3
    )

    dev.off()

    file_tmp <- paste0(file_out,".tmp.png")
    compress(file_in=file_out,file_out=file_tmp)
    file.rename(file_tmp,file_out)

    message("---> Deposition lifetime figure: ",file_out)
  }
}

############################
### COMPOSITE DEP_* PLOTS ###
############################
if (length(dep_variables) > 0) {
  for (dep_name in dep_variables) {
    dep_suffix <- get_variable_suffix(dep_name)
    flux_variables <- paste0(dep_fluxes,"_",dep_suffix)
    available_fluxes <- dep_fluxes[sapply(flux_variables,function(v) v %in% variables_exp1$logical_name && v %in% variables_exp2$logical_name)]

    if (length(available_fluxes) == 0) {
      message("---> Skipping ",dep_name,": no common deposition fluxes available.")
      next
    }

    message("---> Plotting composite ",dep_name," using ",paste(available_fluxes,collapse=", "),if (regional_mode) paste0(" for ",region) else "")
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

      if (regional_mode) {
        field_var1 <- mask_region_field(field_var1,ll$lon,ll$lat,region_box)
        field_var2 <- mask_region_field(field_var2,ll$lon,ll$lat,region_box)
        tmean_var1 <- regional_flux_tg_day(d1$data,ll$lon,ll$lat,region_box)
        tmean_var2 <- regional_flux_tg_day(d2$data,ll$lon,ll$lat,region_box)
      } else {
        tmean_var1 <- global_flux_tg_day(d1$data,ll$lon,ll$lat)
        tmean_var2 <- global_flux_tg_day(d2$data,ll$lon,ll$lat)
      }

      tmean_tim <- seq.POSIXt(from=as.POSIXct(paste0(substr(sDate,1,4),"-",substr(sDate,5,6),"-",substr(sDate,7,8)," 00:00:00"),tz="UTC"),by="3 hours",length.out=d1$nt)
      hour <- as.integer(format(tmean_tim,"%H"))
      hours <- c(0,3,6,9,12,15,18,21)

      massdiag1 <- NULL
      massdiag2 <- NULL
      massdiag_column <- c(ddp="DDEP_FLX",sdm="SEDM_FLX",ngt="NEGA_FIX")[flux]

      if (!regional_mode && !is.na(massdiag_column) && exists("massdiag_compare") && massdiag_compare) {
        massdiag1 <- read_massdiag_series_hours(expname1,logical_name,variables_exp1,seqDate,hours=massdiag_hours,column=massdiag_column)
        massdiag2 <- read_massdiag_series_hours(expname2,logical_name,variables_exp2,seqDate,hours=massdiag_hours,column=massdiag_column)
        massdiag1$value <- -massdiag1$value
        massdiag2$value <- -massdiag2$value
      }

      dhourmean_var1 <- sapply(hours,function(h) mean(tmean_var1[hour == h],na.rm=TRUE))
      dhourmean_var2 <- sapply(hours,function(h) mean(tmean_var2[hour == h],na.rm=TRUE))

      ts_plot <- aggregate_pair_for_plot(tmean_var1,tmean_var2,tmean_tim,dep_name)

      if (ts_plot$aggregation != "3hourly") {
        massdiag1 <- NULL
        massdiag2 <- NULL
      }

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
        plot_var1=ts_plot$value1,
        plot_var2=ts_plot$value2,
        plot_tim=ts_plot$time,
        aggregation=ts_plot$aggregation,
        massdiag1=massdiag1,
        massdiag2=massdiag2,
        dhourmean_var1=dhourmean_var1,
        dhourmean_var2=dhourmean_var2
      )
    }

    wet_massdiag1 <- NULL
    wet_massdiag2 <- NULL

    if (!regional_mode && dep_data[[available_fluxes[1]]]$aggregation == "3hourly" &&
        all(c("wdl","wdc") %in% available_fluxes) &&
        exists("massdiag_compare") && massdiag_compare) {
      wet_logical_name <- paste0("wdl_",dep_suffix)
      wet_massdiag1 <- read_massdiag_series_hours(expname1,wet_logical_name,variables_exp1,seqDate,hours=massdiag_hours,column="WDEP_FLX")
      wet_massdiag2 <- read_massdiag_series_hours(expname2,wet_logical_name,variables_exp2,seqDate,hours=massdiag_hours,column="WDEP_FLX")
      wet_massdiag1$value <- -wet_massdiag1$value
      wet_massdiag2$value <- -wet_massdiag2$value
    }

    plot_category <- get_plot_category(dep_name)
    plot_dir <- make_plot_dir(plot_category)
    file_out <- paste0(plot_dir,gsub(" ","",plot_title),"_",dep_name,"_",expname1,"-",expname2,region_file_tag,"_",sDate,"-",eDate,".png")

    nflux <- length(available_fluxes)
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
    png(file_out,width=(0.2+3*3.9+0.8+0.8)*dpi,height=(0.23+0.15+2*nflux+2.5)*dpi)
    layout(mat=mat,widths=c(0.2,3.9,3.9,0.8,3.9,0.8),heights=c(0.23,0.15,rep(2,nflux),2.5))

    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0("Experiments: ",expname1," VS ",expname2,"   |   Type: ",get_display_type(dep_name),region_title,"   |   Period: ",sDate,"-",eDate),col="grey50",cex=6,family="Century Gothic"); abline(h=c(0,1),col="grey50",lwd=3)
    par(mai=c(0,0,0,0)); plot.new()
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0(expname1," (",exptype1,")"),col="grey20",cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0(expname2," (",exptype2,")"),col="grey20",cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new()
    par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,paste0(expname2," - ",expname1),col="grey20",cex=4.5,family="Century Gothic")
    par(mai=c(0,0,0,0)); plot.new()

    for (flux in available_fluxes) {
      z <- dep_data[[flux]]
      par(mai=c(0,0,0,0)); plot.new(); text(0.5,0.5,plot_type_title[[flux]],col="grey20",cex=map_ylabel_cex,family="Century Gothic",srt=90)

      MapNC(filename_topo="",figure_box=figure_box,field_show_box=field_show_box,coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,drawMapBox=regional_mode,field_value=z$field_var1,field_lon=z$field_lon,field_lat=z$field_lat,field_pallete_name="TROPOMI_NEW",field_breaks=z$field_breaks,field_units=z$units,field_pallete_starting_alpha=100,field_show_legend=FALSE)
      MapNC(filename_topo="",figure_box=figure_box,field_show_box=field_show_box,coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,drawMapBox=regional_mode,field_value=z$field_var2,field_lon=z$field_lon,field_lat=z$field_lat,field_pallete_name="TROPOMI_NEW",field_breaks=z$field_breaks,field_units=z$units,field_pallete_starting_alpha=100,field_show_legend=TRUE,field_legend_mai_right=1.8,field_legend_nlabels=7)
      MapNC(filename_topo="",figure_box=figure_box,field_show_box=field_show_box,coastlineWorldFine_lwd=coastlineWorldFine_lwd,gridlines=gridlines,projection=projection,lonmax=lonmax,lonmin=lonmin,latmax=latmax,latmin=latmin,drawMapBox=regional_mode,field_value=z$field_var2-z$field_var1,field_lon=z$field_lon,field_lat=z$field_lat,field_pallete_name="MNMB",field_breaks=z$field_breaks_diff,field_units=z$units,field_pallete_starting_alpha=100,field_show_legend=TRUE,field_legend_mai_right=1.8,field_legend_nlabels=7)
    }

    tmean_tim <- dep_data[[available_fluxes[1]]]$plot_tim
    x <- seq_along(tmean_tim)
    all_ts <- unlist(lapply(
      available_fluxes,
      function(flux) c(dep_data[[flux]]$plot_var1,dep_data[[flux]]$plot_var2)
    ))

    for (flux in available_fluxes) {
      if (!is.null(dep_data[[flux]]$massdiag1))
        all_ts <- c(all_ts,dep_data[[flux]]$massdiag1$value,dep_data[[flux]]$massdiag2$value)
    }

    wet_output1 <- NULL
    wet_output2 <- NULL
    if (all(c("wdl","wdc") %in% available_fluxes)) {
      wet_output1 <- dep_data[["wdl"]]$plot_var1 + dep_data[["wdc"]]$plot_var1
      wet_output2 <- dep_data[["wdl"]]$plot_var2 + dep_data[["wdc"]]$plot_var2
      all_ts <- c(all_ts,wet_output1,wet_output2)
      if (!is.null(wet_massdiag1))
        all_ts <- c(all_ts,wet_massdiag1$value,wet_massdiag2$value)
    }

    par(mai=ts_panel_mai,family="Century Gothic")
    yseq <- axis_ticks(all_ts,n=10)
    pad <- diff(range(yseq$breaks))*0.05
    if (!is.finite(pad) || pad == 0) pad <- max(abs(yseq$breaks),na.rm=TRUE)*0.05
    if (!is.finite(pad) || pad == 0) pad <- 1

    plot(x,type="n",axes=FALSE,ann=FALSE,ylim=c(min(yseq$breaks)-pad,max(yseq$breaks)+pad),yaxs="i")
    mtext("Time",side=1,line=ts_xlabel_line,cex=ts_xlabel_cex)
    mtext(if (regional_mode) paste0(region," flux (Tg/day)") else "Global flux (Tg/day)",side=2,line=ts_ylabel_line,cex=ts_ylabel_cex)

    IDx_labels <- timeseries_axis_indices(tmean_tim)
    axis(1,at=x[IDx_labels],labels=format(tmean_tim[IDx_labels],"%Y-%m-%d"),cex.axis=ts_xaxis_cex,line=0,lty=0)
    axis(1,at=x[IDx_labels],labels=FALSE,tck=0.01)
    axis(1,at=x[IDx_labels],labels=FALSE,tck=-0.01)
    axis(2,at=yseq$breaks,labels=yseq$labels,las=1,cex.axis=ts_yaxis_cex)
    box(lwd=2)
    abline(h=yseq$breaks,lwd=1,col="grey")
    abline(v=x[IDx_labels],lwd=1,col="grey")

    for (flux in available_fluxes) {
      lines(x,dep_data[[flux]]$plot_var1,lwd=5,col=dep_flux_colors[flux],lty=1)
      points(x,dep_data[[flux]]$plot_var1,pch=19,cex=1.5,col=dep_flux_colors[flux])
      lines(x,dep_data[[flux]]$plot_var2,lwd=5,col=dep_flux_colors[flux],lty=2)
      points(x,dep_data[[flux]]$plot_var2,pch=19,cex=1.5,col=dep_flux_colors[flux])
      add_mean_square(dep_data[[flux]]$plot_var1,dep_flux_colors[flux],cex=2.3)
      add_mean_square(dep_data[[flux]]$plot_var2,dep_flux_colors[flux],cex=2.3)
    }

    pch1 <- get_massdiag_pch(exptype1)
    pch2 <- get_massdiag_pch(exptype2)

    for (flux in available_fluxes) {
      z <- dep_data[[flux]]
      if (!is.null(z$massdiag1)) {
        mx1 <- as.numeric(difftime(z$massdiag1$time,tmean_tim[1],units="hours"))/3 + 1
        mx2 <- as.numeric(difftime(z$massdiag2$time,tmean_tim[1],units="hours"))/3 + 1
        valid1 <- mx1 >= 1 & mx1 <= length(tmean_tim) & is.finite(z$massdiag1$value)
        valid2 <- mx2 >= 1 & mx2 <= length(tmean_tim) & is.finite(z$massdiag2$value)
        points(mx1[valid1],z$massdiag1$value[valid1],pch=pch1,cex=2.6,lwd=2.5,col=dep_flux_colors[flux])
        points(mx2[valid2],z$massdiag2$value[valid2],pch=pch2,cex=2.6,lwd=2.5,col=dep_flux_colors[flux])
      }
    }

    if (!is.null(wet_output1)) {
      lines(x,wet_output1,lwd=5,col=dep_flux_colors["wdep"],lty=1)
      points(x,wet_output1,pch=19,cex=1.5,col=dep_flux_colors["wdep"])
      lines(x,wet_output2,lwd=5,col=dep_flux_colors["wdep"],lty=2)
      points(x,wet_output2,pch=19,cex=1.5,col=dep_flux_colors["wdep"])
      add_mean_square(wet_output1,dep_flux_colors["wdep"],cex=2.3)
      add_mean_square(wet_output2,dep_flux_colors["wdep"],cex=2.3)

      if (!is.null(wet_massdiag1)) {
        mx1 <- as.numeric(difftime(wet_massdiag1$time,tmean_tim[1],units="hours"))/3 + 1
        mx2 <- as.numeric(difftime(wet_massdiag2$time,tmean_tim[1],units="hours"))/3 + 1
        valid1 <- mx1 >= 1 & mx1 <= length(tmean_tim) & is.finite(wet_massdiag1$value)
        valid2 <- mx2 >= 1 & mx2 <= length(tmean_tim) & is.finite(wet_massdiag2$value)
        points(mx1[valid1],wet_massdiag1$value[valid1],pch=pch1,cex=2.6,lwd=2.5,col=dep_flux_colors["wdep"])
        points(mx2[valid2],wet_massdiag2$value[valid2],pch=pch2,cex=2.6,lwd=2.5,col=dep_flux_colors["wdep"])
      }
    }

    legend_fluxes <- toupper(available_fluxes)
    legend_colors <- dep_flux_colors[available_fluxes]
    if (!is.null(wet_output1)) {
      legend_fluxes <- c(legend_fluxes,"WDEP")
      legend_colors <- c(legend_colors,dep_flux_colors["wdep"])
    }

    legend("topleft",legend=legend_fluxes,lwd=6,col=legend_colors,lty=1,cex=2.3,bty="n")
    if (regional_mode) {
      legend("topright",legend=c(expname1,expname2),lwd=6,pch=19,col="grey20",lty=c(1,2),cex=2.3,bty="n")
    } else {
      legend("topright",legend=c(paste0(expname1," OUTPUT"),paste0(expname2," OUTPUT"),paste0(expname1," MASSDIA"),paste0(expname2," MASSDIA")),lwd=c(6,6,NA,NA),pch=c(19,19,get_massdiag_pch(exptype1),get_massdiag_pch(exptype2)),col=c("grey20","grey20","grey20","grey20"),lty=c(1,2,NA,NA),pt.lwd=c(1,1,2.5,2.5),cex=2.0,bty="n")
    }

    all_dc <- unlist(lapply(available_fluxes,function(flux) c(dep_data[[flux]]$dhourmean_var1,dep_data[[flux]]$dhourmean_var2)))
    par(mai=ts_panel_mai,family="Century Gothic")

    yseq <- axis_ticks(all_dc,n=10)
    pad <- diff(range(yseq$breaks))*0.05
    if (!is.finite(pad) || pad == 0) pad <- max(abs(yseq$breaks),na.rm=TRUE)*0.05
    if (!is.finite(pad) || pad == 0) pad <- 1

    plot(1:8,type="n",axes=FALSE,ann=FALSE,ylim=c(min(yseq$breaks)-pad,max(yseq$breaks)+pad),yaxs="i")
    mtext("Time (3 hourly UTC)",side=1,line=ts_xlabel_line,cex=ts_xlabel_cex)
    mtext(if (regional_mode) paste0(region," flux (Tg/day)") else "Global flux (Tg/day)",side=2,line=ts_ylabel_line,cex=ts_ylabel_cex)

    axis(1,at=1:8,labels=c("00","03","06","09","12","15","18","21"),cex.axis=ts_xaxis_cex,line=0,lty=0)
    axis(1,at=1:8,labels=FALSE,tck=0.01)
    axis(1,at=1:8,labels=FALSE,tck=-0.01)
    axis(2,at=yseq$breaks,labels=yseq$labels,las=1,cex.axis=ts_yaxis_cex)

    box(lwd=2)
    abline(h=yseq$breaks,lwd=1,col="grey")
    abline(v=1:8,lwd=1,col="grey")

    for (flux in available_fluxes) {
      lines(1:8,dep_data[[flux]]$dhourmean_var1,lwd=5,col=dep_flux_colors[flux],lty=1)
      points(1:8,dep_data[[flux]]$dhourmean_var1,pch=19,cex=1.5,col=dep_flux_colors[flux])
      lines(1:8,dep_data[[flux]]$dhourmean_var2,lwd=5,col=dep_flux_colors[flux],lty=2)
      points(1:8,dep_data[[flux]]$dhourmean_var2,pch=19,cex=1.5,col=dep_flux_colors[flux])
    }

    legend("topleft",legend=toupper(available_fluxes),lwd=6,col=dep_flux_colors[available_fluxes],lty=1,cex=2.3,bty="n")
    legend("topright",legend=c(expname1,expname2),lwd=6,col="grey20",lty=c(1,2),cex=2.3,bty="n")

    dev.off()
    file_tmp <- paste0(file_out,".tmp.png")
    compress(file_in=file_out,file_out=file_tmp)
    file.rename(file_tmp,file_out)

    message("---> Composite figure: ",file_out)
  }
}
