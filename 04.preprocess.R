############################
### PREPROCESS FUNCTIONS ###
############################

### Return the NetCDF variable name produced by MARS for a GRIB code
grib_to_ncname <- function(grib) {
  x <- strsplit(grib,"\\.")[[1]]
  paste0("p",x[2],sprintf("%03d",as.integer(x[1])))
}

### Return whether a logical variable is stored in PL or SFC files
variable_stream <- function(logical_name,variable_table) {

  rows <- variable_table[variable_table$logical_name == logical_name,,drop=FALSE]
  if (nrow(rows) == 0) return(NA)

  columns <- unique(rows$grib_column)

  if (all(columns == "grib")) return("pl")
  if (all(columns %in% c("gribddp","gribsdm","gribwdl","gribwdc","gribmss","gribngt"))) return("sfc")

  stop("Logical variable '",logical_name,"' uses incompatible GRIB streams: ",paste(columns,collapse=", "))
}

### Get source file for a logical variable
variable_file <- function(logical_name,variable_table,expname,date) {

  stream <- variable_stream(logical_name,variable_table)

  if (stream == "pl") return(paste0(path_data,expname,"/CAMS_",expname,"_forecast00to21by03_0.7x0.7_pl_",date,".nc"))
  if (stream == "sfc") return(paste0(path_data,expname,"/CAMS_",expname,"_forecast00to21by03_0.7x0.7_sfc_",date,".nc"))

  stop("Could not determine stream for variable '",logical_name,"'")
}

### Read one logical variable from a NetCDF file
read_variable <- function(file,exptype,logical_name,variable_table) {

  rows <- variable_table[variable_table$logical_name == logical_name,,drop=FALSE]
  if (nrow(rows) == 0) stop("Variable '",logical_name,"' is not available for ",exptype)

  nc <- nc_open(file)
  on.exit(nc_close(nc))

  field <- NULL

  for (i in seq_len(nrow(rows))) {

    ncname <- grib_to_ncname(rows$grib[i])

    if (!ncname %in% names(nc$var)) stop("NetCDF variable '",ncname,"' for ",logical_name," not found in ",file)

    temp <- ncvar_get(nc,ncname)

    if (is.null(field)) field <- temp else field <- field + temp
  }

  field
}

### Read one pressure-level logical variable at a selected pressure level
read_variable_level <- function(file,exptype,logical_name,variable_table,level) {

  rows <- variable_table[variable_table$logical_name == logical_name,,drop=FALSE]
  if (nrow(rows) == 0) stop("Variable '",logical_name,"' is not available for ",exptype)

  nc <- nc_open(file)
  on.exit(nc_close(nc))

  if (!"level" %in% names(nc$dim) && !"level" %in% names(nc$var)) stop("Pressure-level coordinate not found in ",file)

  levels <- ncvar_get(nc,"level")
  ilev <- which(levels == level)

  if (length(ilev) == 0) stop("Pressure level ",level," hPa not found in ",file)

  field <- NULL

  for (i in seq_len(nrow(rows))) {

    ncname <- grib_to_ncname(rows$grib[i])

    if (!ncname %in% names(nc$var)) stop("NetCDF variable '",ncname,"' for ",logical_name," not found in ",file)

    temp <- ncvar_get(nc,ncname)
    temp <- temp[,,ilev,,drop=FALSE]
    temp <- drop(temp)

    if (is.null(field)) field <- temp else field <- field + temp
  }

  field
}

### Read column burden
read_column_burden <- function(file,exptype,logical_name,variable_table) {

  rows <- variable_table[variable_table$logical_name == logical_name,,drop=FALSE]
  if (nrow(rows) == 0) stop("Variable '",logical_name,"' is not available for ",exptype)

  nc <- nc_open(file)
  on.exit(nc_close(nc))

  levels <- ncvar_get(nc,"level")

  dp <- c(
    50,100,150,200,250,650,1000,1500,2000,2500,
    4000,5000,5000,5000,7500,10000,15000,17500,
    11250,7500,5075
  )

  if (length(levels) != length(dp)) stop("Unexpected number of pressure levels in ",file)

  g <- 9.80665
  burden <- NULL

  for (i in seq_len(nrow(rows))) {

    ncname <- grib_to_ncname(rows$grib[i])

    if (!ncname %in% names(nc$var)) {
      stop("NetCDF variable '",ncname,"' for ",logical_name," not found in ",file)
    }

    q <- ncvar_get(nc,ncname)

    temp <- apply(
      sweep(q,3,dp/g,"*"),
      c(1,2,4),
      sum,
      na.rm=TRUE
    )

    if (is.null(burden)) burden <- temp else burden <- burden + temp
  }

  burden
}

### GRIDCELL AREA
gridcell_area <- function(lon,lat) {
  R <- 6371000
  nlon <- length(lon)
  nlat <- length(lat)
  lat_rad <- lat*pi/180
  lat_edges <- numeric(nlat+1)
  lat_edges[2:nlat] <- (lat_rad[1:(nlat-1)] + lat_rad[2:nlat])/2

  if (lat[1] > lat[nlat]) {
    lat_edges[1] <- pi/2
    lat_edges[nlat+1] <- -pi/2
  } else {
    lat_edges[1] <- -pi/2
    lat_edges[nlat+1] <- pi/2
  }

  dlon <- 2*pi/nlon
  area_lat <- R^2 * dlon * abs(sin(lat_edges[1:nlat]) - sin(lat_edges[2:(nlat+1)]))
  matrix(rep(area_lat,each=nlon),nrow=nlon,ncol=nlat)
}

global_mass_tg <- function(field,lon,lat) {
  area <- gridcell_area(lon,lat)

  if (length(dim(field)) == 2) return(sum(field*area,na.rm=TRUE)/1e9)

  if (length(dim(field)) == 3) {
    result <- numeric(dim(field)[3])
    for (t in seq_len(dim(field)[3])) result[t] <- sum(field[,,t]*area,na.rm=TRUE)/1e9
    return(result)
  }

  stop("global_mass_tg expects a 2-D or 3-D field.")
}

global_flux_tg_day <- function(field,lon,lat) {
  area <- gridcell_area(lon,lat)
  if (length(dim(field)) == 2) return(sum(field*area,na.rm=TRUE)*86400/1e9)

  if (length(dim(field)) == 3) {
    result <- numeric(dim(field)[3])
    for (t in seq_len(dim(field)[3])) result[t] <- sum(field[,,t]*area,na.rm=TRUE)*86400/1e9
    return(result)
  }

  stop("global_flux_tg_day expects a 2-D or 3-D field.")
}

massdiag_file <- function(expname,date) {
  paste0(path_data,expname,"/massdiag/massdia_chem__",expname,"_",date,"00.txt")
}

read_massdiag_value_at_hour <- function(expname,logical_name,variable_table,date,hour,column="TOT_MASS") {
  file <- massdiag_file(expname,date)

  if (!file.exists(file)) {
    message("---> MASSDIA missing: ",file)
    return(NA_real_)
  }

  rows <- variable_table[variable_table$logical_name == logical_name,,drop=FALSE]
  if (nrow(rows) == 0) stop("No resolved tracers for MASSDIA variable ",logical_name)

  tracer_names <- unique(toupper(trimws(rows$massdiag_name)))
  md <- read.table(file,header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  md <- md[as.numeric(md$SIM_HOUR) == hour,,drop=FALSE]

  if (nrow(md) == 0) {
    message("---> MASSDIA missing SIM_HOUR=",hour," in ",file)
    return(NA_real_)
  }

  md_names <- toupper(trimws(md$NAME))
  idx <- match(tracer_names,md_names)

  if (any(is.na(idx))) stop("MASSDIA tracer(s) not found for ",logical_name," in ",file,": ",paste(tracer_names[is.na(idx)],collapse=", "))

  sum(as.numeric(md[[column]][idx]),na.rm=TRUE)
}

read_massdiag_series_hours <- function(expname,logical_name,variable_table,dates,hours=c(6,12,18),column="TOT_MASS") {
  grid <- expand.grid(date=dates,hour=hours,stringsAsFactors=FALSE)
  values <- mapply(function(date,hour) read_massdiag_value_at_hour(expname,logical_name,variable_table,date,hour,column),grid$date,grid$hour)
  times <- as.POSIXct(paste0(substr(grid$date,1,4),"-",substr(grid$date,5,6),"-",substr(grid$date,7,8)," ",sprintf("%02d",grid$hour),":00:00"),tz="UTC")
  list(time=times,value=as.numeric(values))
}


