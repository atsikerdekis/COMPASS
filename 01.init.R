##################
### INITIALIZE ###
##################
username <- system("whoami",intern=TRUE)

path_function <- paste0("/perm/",username,"/cams2_35/COMPASS/code/function/")
path_plot     <- paste0("/perm/",username,"/cams2_35/COMPASS/plot/",expname1,"_",expname2,"_",sDate,"-",eDate,"/")
path_log      <- paste0("/perm/",username,"/cams2_35/COMPASS/log/")
path_data     <- paste0("/scratch/",username,"/cams2_35/data/")
path_PYTHON   <- paste0("/etc/ecmwf/nfs/dh1_perm_b/",username,"/miniforge3/envs/COMPASS/bin/")
path_R        <- paste0("/etc/ecmwf/nfs/dh1_perm_b/",username,"/miniforge3/envs/COMPASS/bin/")

dir.create(path_plot,recursive=TRUE,showWarnings=FALSE)
dir.create(path_log,recursive=TRUE,showWarnings=FALSE)
dir.create(path_data,recursive=TRUE,showWarnings=FALSE)

#################
### LIBRARIES ###
#################
library(ncdf4)
library(showtext)

font_add("Century Gothic",regular="/home/nktt/fonts/CenturyGothic/centurygothic.ttf")
showtext_auto()

#################
### FUNCTIONS ###
#################
source(paste0(path_function,"MapNC.R"))
source(paste0(path_function,"compress.R"))
source(paste0(path_function,"SubmitJob.R"))

Sys.setenv(PROJ_LIB=paste0("/etc/ecmwf/nfs/dh1_perm_b/",username,"/miniforge3/envs/COMPASS/share/proj"))

########################
### EXPERIMENT TYPES ###
########################
supported_exptypes <- c("AER","HAM")

if (!exptype1 %in% supported_exptypes) stop("Unsupported experiment type: ",exptype1)
if (!exptype2 %in% supported_exptypes) stop("Unsupported experiment type: ",exptype2)

###################
### GRIB TABLES ###
###################
grib_AER <- read.csv(paste0(path_code,grib_table_AER),stringsAsFactors=FALSE,colClasses="character",check.names=FALSE)
grib_HAM <- read.csv(paste0(path_code,grib_table_HAM),stringsAsFactors=FALSE,colClasses="character",check.names=FALSE)

names(grib_AER) <- trimws(names(grib_AER))
names(grib_HAM) <- trimws(names(grib_HAM))

grib_AER[] <- lapply(grib_AER,trimws)
grib_HAM[] <- lapply(grib_HAM,trimws)

#############################
### DIAGNOSTIC GRIB COLUMNS ###
#############################
variable_columns <- c(
  mmr = "grib",
  ddp = "gribddp",
  sdm = "gribsdm",
  wdl = "gribwdl",
  wdc = "gribwdc"
)

required_columns <- c("name",unname(variable_columns))

missing_HAM <- setdiff(required_columns,names(grib_HAM))
missing_AER <- setdiff(required_columns,names(grib_AER))

if (length(missing_HAM) > 0) stop("Missing HAM GRIB-table columns: ",paste(missing_HAM,collapse=", "))
if (length(missing_AER) > 0) stop("Missing AER GRIB-table columns: ",paste(missing_AER,collapse=", "))

##########################
### NAME NORMALIZATION ###
##########################
normalize_csv_name <- function(x) {
  x <- trimws(x)
  x <- gsub("-","_",x)
  x <- gsub(" ","_",x)
  tolower(x)
}

##############################
### HAM VARIABLE RESOLUTION ###
##############################
resolve_HAM_variable <- function(logical_name) {

  prefix <- sub("_.*$","",logical_name)
  if (!prefix %in% names(variable_columns)) stop("Unsupported variable prefix '",prefix,"' in ",logical_name)

  suffix <- sub("^[^_]+_","",logical_name)
  grib_column <- variable_columns[[prefix]]

  idx <- which(normalize_csv_name(grib_HAM$name) == suffix)

  if (length(idx) == 0) stop("HAM variable '",logical_name,"' could not be matched to the HAM GRIB table.")
  if (length(idx) > 1) stop("HAM variable '",logical_name,"' matches multiple rows in the HAM GRIB table.")

  csv_name <- grib_HAM$name[idx]
  grib_code <- grib_HAM[[grib_column]][idx]

  if (is.na(grib_code) || grib_code == "") stop("HAM variable '",logical_name,"' has no value in column '",grib_column,"'.")

  data.frame(logical_name=logical_name,csv_name=csv_name,grib_column=grib_column,grib=grib_code,stringsAsFactors=FALSE)
}

##############################
### AER VARIABLE RESOLUTION ###
##############################
AER_aliases <- list(
  du_acc = c("DUS","DUM","DUL")
)

resolve_AER_variable <- function(logical_name) {

  prefix <- sub("_.*$","",logical_name)
  if (!prefix %in% names(variable_columns)) stop("Unsupported variable prefix '",prefix,"' in ",logical_name)

  suffix <- sub("^[^_]+_","",logical_name)
  grib_column <- variable_columns[[prefix]]

  if (!suffix %in% names(AER_aliases)) stop("AER variable '",logical_name,"' is not yet defined in AER_aliases.")

  csv_names <- AER_aliases[[suffix]]
  result <- data.frame()

  for (csv_name in csv_names) {

    idx <- which(normalize_csv_name(grib_AER$name) == normalize_csv_name(csv_name))

    if (length(idx) == 0) stop("AER CSV name '",csv_name,"' not found for variable ",logical_name)
    if (length(idx) > 1) stop("AER CSV name '",csv_name,"' occurs multiple times.")

    grib_code <- grib_AER[[grib_column]][idx]

    if (is.na(grib_code) || grib_code == "") stop("AER variable '",logical_name,"' / ",csv_name," has no value in column '",grib_column,"'.")

    result <- rbind(result,data.frame(logical_name=logical_name,csv_name=grib_AER$name[idx],grib_column=grib_column,grib=grib_code,stringsAsFactors=FALSE))
  }

  result
}

###################################
### RESOLVE REQUESTED VARIABLES ###
###################################
resolve_variables <- function(exptype,variables) {

  result <- data.frame()

  for (logical_name in variables) {

    if (exptype == "HAM") temp <- resolve_HAM_variable(logical_name)

    if (exptype == "AER") {
      suffix <- sub("^[^_]+_","",logical_name)
      if (!suffix %in% names(AER_aliases)) next
      temp <- resolve_AER_variable(logical_name)
    }

    result <- rbind(result,temp)
  }

  result
}

variables_exp1 <- resolve_variables(exptype1,variables)
variables_exp2 <- resolve_variables(exptype2,variables)

###################
### INFORMATION ###
###################
message("---> ",expname1," (",exptype1,"): ",paste(unique(variables_exp1$logical_name),collapse=", "))
message("---> ",expname2," (",exptype2,"): ",paste(unique(variables_exp2$logical_name),collapse=", "))

if (nrow(variables_exp1) > 0) message("---> ",expname1," GRIBs: ",paste(unique(variables_exp1$grib),collapse="/"))
if (nrow(variables_exp2) > 0) message("---> ",expname2," GRIBs: ",paste(unique(variables_exp2$grib),collapse="/"))
