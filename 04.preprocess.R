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
  if (all(columns %in% c("gribddp","gribsdm","gribwdl","gribwdc","gribmss"))) return("sfc")

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
