############################
### PREPROCESS FUNCTIONS ###
############################


### Aggregate a one-dimensional time series for plotting.
### Maps are intentionally unaffected by these modifiers.
aggregate_time_series <- function(values,times,aggregation="3hourly") {

  if (aggregation == "3hourly")
    return(list(time=times,value=values))

  if (aggregation == "daily") {
    key <- format(times,"%Y-%m-%d")
    out <- tapply(values,key,mean,na.rm=TRUE)
    tt <- as.POSIXct(paste0(names(out)," 12:00:00"),tz="UTC")
    return(list(time=tt,value=as.numeric(out)))
  }

  if (aggregation == "monthly") {
    key <- format(times,"%Y-%m")
    out <- tapply(values,key,mean,na.rm=TRUE)
    tt <- as.POSIXct(paste0(names(out),"-15 12:00:00"),tz="UTC")
    return(list(time=tt,value=as.numeric(out)))
  }

  stop("Unsupported time aggregation: ",aggregation)
}

aggregate_pair_for_plot <- function(values1,values2,times,logical_name) {
  aggregation <- get_time_aggregation(logical_name)[1]
  a1 <- aggregate_time_series(values1,times,aggregation)
  a2 <- aggregate_time_series(values2,times,aggregation)
  list(
    time=a1$time,
    value1=a1$value,
    value2=a2$value,
    aggregation=aggregation
  )
}

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
  if (all(columns %in% c("gribddp","gribsdm","gribwdl","gribwdc","gribmss","gribod","gribngt"))) return("sfc")

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


###################################
### RELATIVE HUMIDITY FUNCTIONS ###
###################################

### L137 half-level a/b coefficients needed for the six supported RH levels.
### p_half(n) = a(n) + b(n)*ps
### p_full(k) = 0.5*(p_half(k-1) + p_half(k))
l137_half_coeff <- data.frame(
  n = c(104,105,109,110,113,114,117,118,123,124,136,137),
  a = c(
    12668.257813,11901.339844,
    8880.453125,8163.375000,
    6168.531250,5564.382813,
    3955.960938,3489.234375,
    1659.476563,1387.546875,
    0.000000,0.000000
  ),
  b = c(
    0.549301,0.576692,
    0.680643,0.704669,
    0.770798,0.790717,
    0.843881,0.859432,
    0.922096,0.931881,
    0.997630,1.000000
  )
)

### Each RH level is stored in its own ML file.
rh_ml_file <- function(expname,date,logical_name) {
  fnew <- paste0(path_data,expname,"/CAMS_",expname,
                 "_forecast00to21by03_0.7x0.7_ml_rh_",date,".nc")
  if (file.exists(fnew)) return(fnew)
  model_level <- get_rh_definition(logical_name)$model_level
  fold <- paste0(path_data,expname,"/CAMS_",expname,
                 "_forecast00to21by03_0.7x0.7_ml_",model_level,"_",date,".nc")
  fold
}

wat_ml_file <- function(expname,date) {
  paste0(path_data,expname,"/CAMS_",expname,
         "_forecast00to21by03_0.7x0.7_ml_wat_",date,".nc")
}

rh_lnsp_file <- function(expname,date) {
  paste0(path_data,expname,"/CAMS_",expname,
         "_forecast00to21by03_0.7x0.7_lnsp_",date,".nc")
}

model_level_pressure <- function(ps,model_level) {
  upper <- l137_half_coeff[l137_half_coeff$n == model_level-1,,drop=FALSE]
  lower <- l137_half_coeff[l137_half_coeff$n == model_level,,drop=FALSE]

  if (nrow(upper) != 1 || nrow(lower) != 1)
    stop("Missing L137 a/b coefficients for model level ",model_level)

  p_upper <- upper$a + upper$b*ps
  p_lower <- lower$a + lower$b*ps

  0.5*(p_upper+p_lower)
}

### Saturation vapour pressure [Pa], Buck formulation over liquid water.
saturation_vapour_pressure <- function(T) {
  Tc <- T-273.15
  611.21*exp((18.678-Tc/234.5)*(Tc/(257.14+Tc)))
}

### Relative humidity [%] from q [kg kg-1], T [K], p [Pa].
relative_humidity_from_qtp <- function(q,T,p) {
  epsilon <- 0.621981
  e <- q*p/(epsilon+(1-epsilon)*q)
  es <- saturation_vapour_pressure(T)
  rh <- 100*e/es
  rh[rh < 0] <- NA_real_
  rh
}

### Read/calculate one RH diagnostic for one day.
### The ML files contain exactly one requested model level. MARS may therefore
### omit the singleton vertical coordinate completely; the level is known from
### the logical variable name and filename.
read_relative_humidity <- function(expname,date,logical_name) {

  def <- get_rh_definition(logical_name)
  model_level <- def$model_level

  file_ml <- rh_ml_file(expname,date,logical_name)
  file_lnsp <- rh_lnsp_file(expname,date)

  if (!file.exists(file_ml)) stop("RH model-level file not found: ",file_ml)
  if (!file.exists(file_lnsp)) stop("RH lnsp file not found: ",file_lnsp)

  nc_ml <- nc_open(file_ml)
  on.exit(nc_close(nc_ml),add=TRUE)

  nc_lnsp <- nc_open(file_lnsp)
  on.exit(nc_close(nc_lnsp),add=TRUE)

  T_name <- grib_to_ncname("130.128")
  q_name <- grib_to_ncname("133.128")
  lnsp_name <- grib_to_ncname("152.128")

  if (!T_name %in% names(nc_ml$var)) stop("Temperature variable ",T_name," not found in ",file_ml)
  if (!q_name %in% names(nc_ml$var)) stop("Specific humidity variable ",q_name," not found in ",file_ml)
  if (!lnsp_name %in% names(nc_lnsp$var)) stop("lnsp variable ",lnsp_name," not found in ",file_lnsp)

  T_raw <- ncvar_get(nc_ml,T_name)
  q_raw <- ncvar_get(nc_ml,q_name)
  lnsp <- drop(ncvar_get(nc_lnsp,lnsp_name))

  extract_ml <- function(arr,varname) {
    d <- dim(arr)
    if (length(d) == 3) return(drop(arr))
    if (length(d) != 4) stop("Unexpected dimensions for ",varname," in ",file_ml,": ",paste(d,collapse="x"))
    levs <- if ("level" %in% names(nc_ml$dim) || "level" %in% names(nc_ml$var)) ncvar_get(nc_ml,"level") else NULL
    if (is.null(levs)) stop("Level coordinate missing in ",file_ml)
    ilev <- which(levs == model_level)
    if (length(ilev) != 1) stop("Model level ",model_level," not found in ",file_ml)
    drop(arr[,,ilev,])
  }

  T <- extract_ml(T_raw,T_name)
  q <- extract_ml(q_raw,q_name)

  if (length(dim(T)) != 3 || length(dim(q)) != 3)
    stop("Unexpected T/q dimensions in ",file_ml,
         ". Expected lon x lat x time after selecting the target ML. T=",
         paste(dim(T),collapse="x")," q=",paste(dim(q),collapse="x"))

  if (length(dim(lnsp)) != 3)
    stop("Unexpected lnsp dimensions in ",file_lnsp,
         ". Expected lon x lat x time after dropping singleton dimensions. lnsp=",
         paste(dim(lnsp),collapse="x"))

  ps <- exp(lnsp)
  p <- model_level_pressure(ps,model_level)

  if (!all(dim(T) == dim(q)) || !all(dim(T) == dim(p)))
    stop("T, q and pressure dimensions do not match for ",logical_name,
         ": T=",paste(dim(T),collapse="x"),
         ", q=",paste(dim(q),collapse="x"),
         ", p=",paste(dim(p),collapse="x"))

  relative_humidity_from_qtp(q,T,p)
}



##################################
### OPTICAL PROPERTY FUNCTIONS ###
##################################
optics_file <- function(expname,date) {
  paste0(path_data,expname,"/CAMS_",expname,"_forecast00to21by03_0.7x0.7_sfc_optics_",date,".nc")
}
mec_file <- function(expname,date) {
  paste0(path_data,expname,"/CAMS_",expname,"_forecast00to21by03_0.7x0.7_sfc_mec_",date,".nc")
}
read_direct_optical <- function(expname,date,grib) {
  file <- optics_file(expname,date)
  if (!file.exists(file)) stop("Optics file not found: ",file)
  nc <- nc_open(file); on.exit(nc_close(nc))
  ncname <- grib_to_ncname(grib)
  if (!ncname %in% names(nc$var)) stop("Optical variable ",ncname," not found in ",file)
  drop(ncvar_get(nc,ncname))
}
read_total_aerosol_burden <- function(expname,date) {
  file <- mec_file(expname,date)
  if (!file.exists(file)) stop("MEC burden file not found: ",file)
  nc <- nc_open(file); on.exit(nc_close(nc))
  vars <- names(nc$var)
  vars <- vars[grepl("^p[0-9]+$",vars)]
  if (length(vars) == 0) stop("No aerosol burden variables found in ",file)
  out <- NULL
  for (v in vars) {
    x <- ncvar_get(nc,v)
    if (is.null(out)) out <- x else out <- out+x
  }
  drop(out)
}

resolve_species_od_table <- function(exptype,include_soa=TRUE) {
  vars <- paste0("od_",get_od_species_suffixes(include_soa=include_soa))
  if (exptype == "HAM") return(resolve_variables("HAM",vars))
  if (exptype == "AER") return(resolve_variables("AER",vars))
  stop("Unsupported experiment type for species OD: ",exptype)
}

read_species_od_components <- function(expname,date,exptype) {
  table <- resolve_species_od_table(exptype,include_soa=TRUE)
  out <- list()
  for (suffix in get_od_species_suffixes(include_soa=TRUE)) {
    lname <- paste0("od_",suffix)
    rows <- table[table$logical_name == lname,,drop=FALSE]
    if (nrow(rows) == 0) next
    file <- variable_file(lname,table,expname,date)
    out[[suffix]] <- read_variable(file,exptype,lname,table)
  }
  out
}

read_species_od_sum <- function(expname,date,exptype) {
  comps <- read_species_od_components(expname,date,exptype)
  if (length(comps) == 0) stop("No species OD components found for ",expname)
  out <- NULL
  for (nm in names(comps)) {
    if (is.null(out)) out <- comps[[nm]] else out <- out + comps[[nm]]
  }
  out
}

read_water_aod <- function(expname,date,logical_name) {
  base <- strip_time_aggregation(logical_name)
  file <- wat_ml_file(expname,date)
  if (!file.exists(file)) stop("Water-AOD ML file not found: ",file)
  nc <- nc_open(file)
  on.exit(nc_close(nc))
  varname <- grib_to_ncname("22.210")
  if (!varname %in% names(nc$var)) stop("NetCDF variable ",varname," not found in ",file)
  raw <- ncvar_get(nc,varname)
  d <- dim(raw)
  if (length(d) != 4) stop("Unexpected water-AOD dimensions in ",file,": ",paste(d,collapse="x"))
  levs <- if ("level" %in% names(nc$dim) || "level" %in% names(nc$var)) ncvar_get(nc,"level") else c(2,3,4)
  pick_level <- function(level_value) {
    ilev <- which(levs == level_value)
    if (length(ilev) != 1) stop("Level ",level_value," not found in ",file)
    drop(raw[,,ilev,])
  }
  if (base == "od_wat_ks") return(pick_level(2))
  if (base == "od_wat_as") return(pick_level(3))
  if (base == "od_wat_cs") return(pick_level(4))
  if (base == "od_wat") return(pick_level(2) + pick_level(3) + pick_level(4))
  stop("Unsupported water variable: ",logical_name)
}
read_optical_property <- function(expname,date,logical_name) {
  logical_name <- strip_time_aggregation(logical_name)
  exptype <- if (exists("exptype1") && expname == expname1) exptype1 else if (exists("exptype2") && expname == expname2) exptype2 else NA
  if (logical_name %in% c("aod550","aod550_species")) return(read_direct_optical(expname,date,"207.210"))
  if (logical_name == "aod865") return(read_direct_optical(expname,date,"215.210"))
  if (logical_name == "aaod550") return(read_direct_optical(expname,date,"104.215"))
  if (logical_name == "ssa550") return(read_direct_optical(expname,date,"140.215"))
  if (logical_name == "od_sum_species") return(read_species_od_sum(expname,date,exptype))
  if (logical_name == "ae550to865") {
    a1 <- read_direct_optical(expname,date,"207.210")
    a2 <- read_direct_optical(expname,date,"215.210")
    out <- -log(a1/a2)/log(550/865)
    out[!is.finite(out)] <- NA_real_
    return(out)
  }
  if (logical_name == "mec550") {
    aod <- read_direct_optical(expname,date,"207.210")
    burden <- read_total_aerosol_burden(expname,date)
    out <- aod/(burden*1000)  ### m2 g-1
    out[!is.finite(out)] <- NA_real_
    return(out)
  }
  stop("Unsupported optical property: ",logical_name)
}



### Representative period-mean maps for ratio-derived optical properties.
### Average physical components first, then calculate the ratio/transform.
read_period_optical_map <- function(expname,dates,logical_name) {

  logical_name <- strip_time_aggregation(logical_name)

  stack_fields <- function(field_list) {
    nx <- dim(field_list[[1]])[1]
    ny <- dim(field_list[[1]])[2]
    nt <- sum(sapply(field_list,function(x) dim(x)[3]))
    array(unlist(field_list),dim=c(nx,ny,nt))
  }

  if (logical_name == "ae550to865") {

    aod550 <- stack_fields(
      lapply(dates,function(d) read_direct_optical(expname,d,"207.210"))
    )
    aod865 <- stack_fields(
      lapply(dates,function(d) read_direct_optical(expname,d,"215.210"))
    )

    mean550 <- apply(aod550,c(1,2),mean,na.rm=TRUE)
    mean865 <- apply(aod865,c(1,2),mean,na.rm=TRUE)

    out <- -log(mean550/mean865)/log(550/865)
    out[!is.finite(out)] <- NA_real_
    return(out)
  }

  if (logical_name == "ssa550") {

    aod <- stack_fields(
      lapply(dates,function(d) read_direct_optical(expname,d,"207.210"))
    )
    aaod <- stack_fields(
      lapply(dates,function(d) read_direct_optical(expname,d,"104.215"))
    )

    mean_aod <- apply(aod,c(1,2),mean,na.rm=TRUE)
    mean_aaod <- apply(aaod,c(1,2),mean,na.rm=TRUE)

    out <- 1-mean_aaod/mean_aod
    out[!is.finite(out) | out < 0 | out > 1] <- NA_real_
    return(out)
  }

  if (logical_name == "mec550") {

    aod <- stack_fields(
      lapply(dates,function(d) read_direct_optical(expname,d,"207.210"))
    )
    burden <- stack_fields(
      lapply(dates,function(d) read_total_aerosol_burden(expname,d))
    )

    mean_aod <- apply(aod,c(1,2),mean,na.rm=TRUE)
    mean_burden <- apply(burden,c(1,2),mean,na.rm=TRUE)

    ### burden is kg m^-2; convert to g m^-2 for m2 g^-1.
    out <- mean_aod/(mean_burden*1000)
    out[!is.finite(out) | out < 0] <- NA_real_
    return(out)
  }

  NULL
}

################################
### TOTAL-COLUMN GAS FUNCTIONS ###
################################
column_file <- function(expname,date) {
  paste0(path_data,expname,"/CAMS_",expname,
         "_forecast00to21by03_0.7x0.7_sfc_columns_",date,".nc")
}

read_total_column <- function(expname,date,logical_name) {

  def <- get_column_definition(logical_name)
  file <- column_file(expname,date)

  if (!file.exists(file))
    stop("Total-column file not found: ",file)

  nc <- nc_open(file)
  on.exit(nc_close(nc))

  ncname <- grib_to_ncname(def$grib)

  if (!ncname %in% names(nc$var))
    stop("Total-column variable ",ncname," for ",logical_name,
         " not found in ",file)

  field <- drop(ncvar_get(nc,ncname))

  if (length(dim(field)) != 3)
    stop("Unexpected total-column dimensions in ",file,
         ". Expected lon x lat x time. Found: ",
         paste(dim(field),collapse="x"))

  field
}

#################################
### PRECIPITATION FUNCTIONS ###
#################################

precip_file <- function(expname,date) {
  paste0(path_data,expname,"/CAMS_",expname,
         "_forecast03to24by03_0.7x0.7_sfc_precip_",date,".nc")
}

### Convert accumulated forecast precipitation to 3-hour amounts.
### IFS precipitation fields are accumulated from forecast step 0.
### The downloader retrieves steps 3,6,...,24, so:
###   first field  = accumulation 00-03 UTC
###   later fields = difference between consecutive accumulated steps.
### Returned units are mm per 3-hour interval.
read_precipitation <- function(expname,date,logical_name) {

  def <- get_precip_definition(logical_name)
  file <- precip_file(expname,date)

  if (!file.exists(file))
    stop("Precipitation file not found: ",file)

  nc <- nc_open(file)
  on.exit(nc_close(nc))

  ncname <- grib_to_ncname(def$grib)

  if (!ncname %in% names(nc$var))
    stop("Precipitation variable ",ncname," for ",logical_name,
         " not found in ",file)

  accumulated <- drop(ncvar_get(nc,ncname))

  if (length(dim(accumulated)) != 3)
    stop("Unexpected precipitation dimensions in ",file,
         ". Expected lon x lat x time. Found: ",
         paste(dim(accumulated),collapse="x"))

  nt <- dim(accumulated)[3]

  if (nt != 8)
    stop("Expected 8 precipitation forecast steps (3...24 h) in ",file,
         ", found ",nt)

  amount <- accumulated

  ### Step 3 is already the 00-03 UTC accumulation.
  amount[,,1] <- accumulated[,,1]

  ### Steps 6...24 are deaccumulated against the previous forecast step.
  for (t in 2:nt)
    amount[,,t] <- accumulated[,,t]-accumulated[,,t-1]

  ### Precipitation parameters are archived in metres of water equivalent.
  amount <- amount*1000

  ### Preserve real problems instead of silently clipping them. Only tiny
  ### negative values caused by numerical precision are set to zero.
  tiny_negative <- amount < 0 & amount > -1e-6
  amount[tiny_negative] <- 0

  if (any(amount < -1e-6,na.rm=TRUE))
    warning("Negative 3-hour precipitation amounts found for ",
            logical_name," in ",file)

  amount
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

### Return longitude/latitude indices inside a region box:
### box = c(lonmin,lonmax,latmin,latmax), longitudes interpreted in -180..180
region_indices <- function(lon,lat,box) {
  lon180 <- ((lon + 180) %% 360) - 180
  lonmin <- box[1]; lonmax <- box[2]; latmin <- box[3]; latmax <- box[4]

  if (lonmin <= lonmax) ilon <- which(lon180 >= lonmin & lon180 <= lonmax)
  else ilon <- which(lon180 >= lonmin | lon180 <= lonmax)
  ilat <- which(lat >= latmin & lat <= latmax)

  if (length(ilon) == 0 || length(ilat) == 0)
    stop("No model grid cells found inside region box: ",paste(box,collapse=", "))

  list(lon=ilon,lat=ilat)
}

### Mask values outside a region while keeping the original global grid/dimensions
mask_region_field <- function(field,lon,lat,box) {
  idx <- region_indices(lon,lat,box)
  out <- field

  if (length(dim(field)) == 2) {
    keep <- matrix(FALSE,nrow=length(lon),ncol=length(lat))
    keep[idx$lon,idx$lat] <- TRUE
    out[!keep] <- NA
    return(out)
  }

  if (length(dim(field)) == 3) {
    keep <- matrix(FALSE,nrow=length(lon),ncol=length(lat))
    keep[idx$lon,idx$lat] <- TRUE
    for (t in seq_len(dim(field)[3])) {
      temp <- out[,,t]
      temp[!keep] <- NA
      out[,,t] <- temp
    }
    return(out)
  }

  stop("mask_region_field expects a 2-D or 3-D field.")
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

regional_mass_tg <- function(field,lon,lat,box) {
  area <- gridcell_area(lon,lat)
  idx <- region_indices(lon,lat,box)

  if (length(dim(field)) == 2)
    return(sum(field[idx$lon,idx$lat]*area[idx$lon,idx$lat],na.rm=TRUE)/1e9)

  if (length(dim(field)) == 3) {
    result <- numeric(dim(field)[3])
    for (t in seq_len(dim(field)[3]))
      result[t] <- sum(field[idx$lon,idx$lat,t]*area[idx$lon,idx$lat],na.rm=TRUE)/1e9
    return(result)
  }

  stop("regional_mass_tg expects a 2-D or 3-D field.")
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

regional_flux_tg_day <- function(field,lon,lat,box) {
  area <- gridcell_area(lon,lat)
  idx <- region_indices(lon,lat,box)

  if (length(dim(field)) == 2)
    return(sum(field[idx$lon,idx$lat]*area[idx$lon,idx$lat],na.rm=TRUE)*86400/1e9)

  if (length(dim(field)) == 3) {
    result <- numeric(dim(field)[3])
    for (t in seq_len(dim(field)[3]))
      result[t] <- sum(field[idx$lon,idx$lat,t]*area[idx$lon,idx$lat],na.rm=TRUE)*86400/1e9
    return(result)
  }

  stop("regional_flux_tg_day expects a 2-D or 3-D field.")
}

regional_mean <- function(field,lon,lat,box) {
  idx <- region_indices(lon,lat,box)

  if (length(dim(field)) == 2)
    return(mean(field[idx$lon,idx$lat],na.rm=TRUE))

  if (length(dim(field)) == 3) {
    result <- numeric(dim(field)[3])
    for (t in seq_len(dim(field)[3]))
      result[t] <- mean(field[idx$lon,idx$lat,t],na.rm=TRUE)
    return(result)
  }

  stop("regional_mean expects a 2-D or 3-D field.")
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
