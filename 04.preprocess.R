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

  if (length(stream) != 1 || is.na(stream))
    stop("Could not determine stream for variable '",logical_name,
         "'. Resolve the variable before calling variable_file().")

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

get_od_species_suffixes <- function(include_soa=TRUE) {
  x <- c("ss","du","pom","bc","so4","ni","am")
  if (include_soa) x <- c(x,"soa")
  x
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
  if (base == "wat_ks") return(pick_level(2))
  if (base == "wat_as") return(pick_level(3))
  if (base == "wat_cs") return(pick_level(4))
  if (base == "wat") return(pick_level(2) + pick_level(3) + pick_level(4))
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



#########################################
### MULTI-PANEL SPECIES-MAP FUNCTIONS ###
#########################################

panel_display_species <- function(exptype,panel_base) {
  panel_base <- tolower(strip_time_aggregation(panel_base))

  ### Requested display order. SOA is never shown separately.
  if (panel_base %in% c("aod_per_species","aodratio_per_species") && exptype == "HAM")
    return(c("total","du","ss","pom","bc","so4","ni","am","wat"))

  ### AER and the mass/MEC panels have eight populated slots; the ninth
  ### position is intentionally blank so every figure keeps the same 3x3 grid.
  c("total","du","ss","pom","bc","so4","ni","am","blank")
}

panel_stack_fields <- function(field_list) {
  if (length(field_list) == 0) stop("Cannot stack an empty field list.")
  nx <- dim(field_list[[1]])[1]
  ny <- dim(field_list[[1]])[2]
  nt <- sum(sapply(field_list,function(x) dim(x)[3]))
  array(unlist(field_list),dim=c(nx,ny,nt))
}

panel_mean_stack <- function(field_list) {
  out <- apply(panel_stack_fields(field_list),c(1,2),mean,na.rm=TRUE)
  out[is.nan(out)] <- NA_real_
  out
}

panel_add_fields <- function(a,b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  both_na <- is.na(a) & is.na(b)
  aa <- a; bb <- b
  aa[is.na(aa)] <- 0
  bb[is.na(bb)] <- 0
  out <- aa+bb
  out[both_na] <- NA_real_
  out
}

read_single_species_mass <- function(expname,date,exptype,suffix) {
  lname <- paste0("mss_",suffix)
  table <- resolve_variables(exptype,lname)
  if (nrow(table) == 0)
    stop("Mass burden variable not available: ",lname," for ",exptype)
  file <- variable_file(lname,table,expname,date)
  if (!file.exists(file)) stop("Mass input file not found: ",file)
  read_variable(file,exptype,lname,table)
}

build_panel_species_maps <- function(expname,exptype,dates,panel_base) {
  panel_base <- tolower(strip_time_aggregation(panel_base))
  display_order <- panel_display_species(exptype,panel_base)
  out <- list()

  need_aod <- panel_base %in% c("aod_per_species","aodratio_per_species","mec_per_species")
  need_mass <- panel_base %in% c("mass_per_species","massratio_per_species","mec_per_species")

  total_aod <- NULL
  species_aod <- list()
  species_mass <- list()
  total_mass <- NULL

  if (need_aod) {
    total_aod <- panel_mean_stack(lapply(dates,function(d) read_direct_optical(expname,d,"207.210")))

    comps_per_day <- lapply(dates,function(d) read_species_od_components(expname,d,exptype))
    for (suffix in c("du","ss","pom","bc","so4","ni","am","soa")) {
      vals <- lapply(comps_per_day,function(x) x[[suffix]])
      vals <- vals[!vapply(vals,is.null,logical(1))]
      if (length(vals) > 0) species_aod[[suffix]] <- panel_mean_stack(vals)
    }

    ### For AER, displayed POM is primary OM + SOA.
    if (exptype == "AER")
      species_aod[["pom"]] <- panel_add_fields(species_aod[["pom"]],species_aod[["soa"]])
    species_aod[["soa"]] <- NULL

    if (exptype == "HAM" && panel_base %in% c("aod_per_species","aodratio_per_species"))
      species_aod[["wat"]] <- panel_mean_stack(lapply(dates,function(d) read_water_aod(expname,d,"wat")))
  }

  if (need_mass) {
    ### Read SOA even though it is not displayed separately. It remains part
    ### of the physical total burden; for AER it is folded into displayed POM.
    for (suffix in c("du","ss","pom","bc","so4","ni","am","soa")) {
      vals <- lapply(dates,function(d) read_single_species_mass(expname,d,exptype,suffix))
      species_mass[[suffix]] <- panel_mean_stack(vals)
    }

    total_mass <- NULL
    for (suffix in names(species_mass))
      total_mass <- panel_add_fields(total_mass,species_mass[[suffix]])

    if (exptype == "AER")
      species_mass[["pom"]] <- panel_add_fields(species_mass[["pom"]],species_mass[["soa"]])
    species_mass[["soa"]] <- NULL
  }

  if (panel_base == "aod_per_species") {
    out[["total"]] <- total_aod
    for (nm in setdiff(display_order,c("total","blank"))) out[[nm]] <- species_aod[[nm]]
  }

  if (panel_base == "aodratio_per_species") {
    total_ratio <- total_aod/total_aod
    total_ratio[!is.finite(total_ratio)] <- NA_real_
    out[["total"]] <- total_ratio
    for (nm in setdiff(display_order,c("total","blank"))) {
      x <- species_aod[[nm]]/total_aod
      x[!is.finite(x) | x < 0] <- NA_real_
      out[[nm]] <- x
    }
  }

  if (panel_base == "mass_per_species") {
    out[["total"]] <- total_mass
    for (nm in setdiff(display_order,c("total","blank"))) out[[nm]] <- species_mass[[nm]]
  }

  if (panel_base == "massratio_per_species") {
    total_ratio <- total_mass/total_mass
    total_ratio[!is.finite(total_ratio)] <- NA_real_
    out[["total"]] <- total_ratio
    for (nm in setdiff(display_order,c("total","blank"))) {
      x <- species_mass[[nm]]/total_mass
      x[!is.finite(x) | x < 0] <- NA_real_
      out[[nm]] <- x
    }
  }

  if (panel_base == "mec_per_species") {
    total_mec <- total_aod/(total_mass*1000)
    total_mec[!is.finite(total_mec) | total_mec < 0 | total_aod < 0.05] <- NA_real_
    out[["total"]] <- total_mec

    for (nm in setdiff(display_order,c("total","blank"))) {
      aod <- species_aod[[nm]]
      mass <- species_mass[[nm]]
      x <- aod/(mass*1000)
      x[!is.finite(x) | x < 0 | aod < 0.05] <- NA_real_
      out[[nm]] <- x
    }
  }

  ### Return exactly nine ordered map slots.
  template <- out[["total"]]
  ordered <- list()
  for (nm in display_order) {
    if (nm == "blank" || is.null(out[[nm]])) {
      blank <- template
      blank[] <- NA_real_
      ordered[[nm]] <- blank
    } else {
      ordered[[nm]] <- out[[nm]]
    }
  }
  ordered
}

############################################
### SATELLITE PANEL VALIDATION FUNCTIONS ###
############################################

panel_satellite_file <- function(sensor,date) {
  def <- get_panel_satellite_definition(sensor)
  paste0(panel_satellite_path,def$file_prefix,date,"_3H.nc")
}

read_panel_satellite_field <- function(sensor,date,variable=panel_validation_variable) {
  def <- get_panel_validation_definition(variable)
  file <- panel_satellite_file(sensor,date)
  if (!file.exists(file)) {
    warning("Satellite panel file not found: ",file)
    return(NULL)
  }

  nc <- nc_open(file)
  on.exit(nc_close(nc))

  vname <- paste0(sensor,"_",variable,"_mn")
  if (!vname %in% names(nc$var)) {
    warning("Satellite variable ",vname," not found in ",file)
    return(NULL)
  }

  field <- drop(ncvar_get(nc,vname))
  if (length(dim(field)) != 3)
    stop("Unexpected satellite dimensions for ",vname," in ",file,": ",paste(dim(field),collapse="x"))

  if (variable == "AE550to860") {
    aod_name <- paste0(sensor,"_AOD550_mn")
    if (!aod_name %in% names(nc$var))
      stop("AOD550 field required to filter ",vname," is missing in ",file)
    aod <- drop(ncvar_get(nc,aod_name))
    field[!is.finite(aod) | aod < def$aod_filter] <- NA_real_
  }

  lon <- ncvar_get(nc,"longitude")
  lat <- ncvar_get(nc,"latitude")
  list(field=field,lon=lon,lat=lat)
}

panel_array_ensemble_mean <- function(fields) {
  fields <- fields[!vapply(fields,is.null,logical(1))]
  if (length(fields) == 0) return(NULL)

  refdim <- dim(fields[[1]])
  if (any(!vapply(fields,function(x) identical(dim(x),refdim),logical(1))))
    stop("Satellite ensemble fields do not share the same dimensions.")

  total <- array(0,dim=refdim)
  count <- array(0L,dim=refdim)
  for (x in fields) {
    ok <- is.finite(x)
    total[ok] <- total[ok]+x[ok]
    count[ok] <- count[ok]+1L
  }
  out <- total/count
  out[count == 0] <- NA_real_
  out
}

panel_nearest_lon_index <- function(source_lon,target_lon) {
  src <- source_lon %% 360
  trg <- target_lon %% 360
  vapply(trg,function(x) which.min(abs(((src-x+180) %% 360)-180)),integer(1))
}

panel_nearest_lat_index <- function(source_lat,target_lat) {
  vapply(target_lat,function(x) which.min(abs(source_lat-x)),integer(1))
}

panel_model_to_satellite_grid <- function(field,model_lon,model_lat,sat_lon,sat_lat) {
  ilon <- panel_nearest_lon_index(model_lon,sat_lon)
  ilat <- panel_nearest_lat_index(model_lat,sat_lat)
  field[ilon,ilat,,drop=FALSE]
}

read_panel_validation_model <- function(expname,date,variable=panel_validation_variable) {
  def <- get_panel_validation_definition(variable)
  if (def$model_variable == "aod550") return(read_direct_optical(expname,date,"207.210"))
  if (def$model_variable == "ae550to865") return(read_optical_property(expname,date,"ae550to865"))
  stop("Unsupported panel model validation variable: ",def$model_variable)
}

build_panel_satellite_validation <- function(expname,dates,variable=panel_validation_variable) {
  obs_days <- list()
  mod_days <- list()
  sat_lon <- sat_lat <- NULL

  model_ll <- NULL
  for (d in dates) {
    sensor_data <- lapply(panel_satellites,function(s) read_panel_satellite_field(s,d,variable))
    sensor_data <- sensor_data[!vapply(sensor_data,is.null,logical(1))]
    if (length(sensor_data) == 0) {
      warning("No satellite panel data available for date ",d)
      next
    }

    if (is.null(sat_lon)) {
      sat_lon <- sensor_data[[1]]$lon
      sat_lat <- sensor_data[[1]]$lat
    }

    obs <- panel_array_ensemble_mean(lapply(sensor_data,function(x) x[["field"]]))
    if (is.null(obs)) next

    model <- read_panel_validation_model(expname,d,variable)
    if (is.null(model_ll)) {
      nc <- nc_open(optics_file(expname,d))
      model_ll <- list(lon=ncvar_get(nc,"longitude"),lat=ncvar_get(nc,"latitude"))
      nc_close(nc)
    }

    model_sat <- panel_model_to_satellite_grid(model,model_ll$lon,model_ll$lat,sat_lon,sat_lat)
    nt <- min(dim(obs)[3],dim(model_sat)[3])
    obs <- obs[,,seq_len(nt),drop=FALSE]
    model_sat <- model_sat[,,seq_len(nt),drop=FALSE]
    model_sat[!is.finite(obs)] <- NA_real_

    obs_days[[length(obs_days)+1]] <- obs
    mod_days[[length(mod_days)+1]] <- model_sat
  }

  if (length(obs_days) == 0)
    stop("No satellite data available for panel validation over ",paste(range(dates),collapse="-"))

  obs_stack <- panel_stack_fields(obs_days)
  mod_stack <- panel_stack_fields(mod_days)

  obs_mean <- apply(obs_stack,c(1,2),mean,na.rm=TRUE)
  me <- apply(mod_stack-obs_stack,c(1,2),mean,na.rm=TRUE)
  mae <- apply(abs(mod_stack-obs_stack),c(1,2),mean,na.rm=TRUE)
  obs_mean[is.nan(obs_mean)] <- NA_real_
  me[is.nan(me)] <- NA_real_
  mae[is.nan(mae)] <- NA_real_

  list(obs=obs_mean,me=me,mae=mae,lon=sat_lon,lat=sat_lat)
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
