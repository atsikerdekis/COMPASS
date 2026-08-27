##################
### INITIALIZE ###
##################
username <- system("whoami",intern=TRUE)

path_function <- paste0("/perm/",username,"/cams2_35/COMPASS/code_refactored/function/")
path_plot     <- paste0("/perm/",username,"/cams2_35/COMPASS/plot/",expname1,"_",expname2,"_",sDate,"-",eDate,"/")
path_log      <- paste0("/perm/",username,"/cams2_35/COMPASS/log/")
path_data     <- paste0("/scratch/",username,"/cams2_35/data/")
path_temp     <- paste0("/scratch/",username,"/cams2_35/temp/",expname1,"_",expname2,"_",sDate,"-",eDate,"/")
path_PYTHON   <- paste0("/etc/ecmwf/nfs/dh1_perm_b/",username,"/miniforge3/envs/COMPASS/bin/")
path_R        <- paste0("/etc/ecmwf/nfs/dh1_perm_b/",username,"/miniforge3/envs/COMPASS/bin/")

dir.create(path_plot,recursive=TRUE,showWarnings=FALSE)
dir.create(path_log,recursive=TRUE,showWarnings=FALSE)
dir.create(path_data,recursive=TRUE,showWarnings=FALSE)
dir.create(path_temp,recursive=TRUE,showWarnings=FALSE)

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

################################
### DIAGNOSTIC GRIB COLUMNS ###
################################
variable_columns <- c(
  mmr = "grib",
  ddp = "gribddp",
  sdm = "gribsdm",
  wdl = "gribwdl",
  wdc = "gribwdc",
  mss = "gribmss",
  ngt = "gribngt"
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

################################
### HAM7 MODE DEFINITIONS ###
################################
ham_modes <- c("ns","ks","as","cs","ki","ai","ci")

ham_soluble_modes   <- c("ns","ks","as","cs")
ham_insoluble_modes <- c("ki","ai","ci")

ham_names <- normalize_csv_name(grib_HAM$name)

### Only rows of form SPECIES_MODE are considered aerosol species/mode components.
### This automatically excludes AS_N, KS_N, CDNC, ICNC, etc.
ham_component_rows <- grepl(
  paste0("_(",paste(ham_modes,collapse="|"),")$"),
  ham_names
)

######################################
### EXPAND COMPOSITE DEP VARIABLES ###
######################################
variables_requested <- variables
dep_fluxes <- c("ddp","sdm","wdl","wdc","ngt")
dep_variables <- variables_requested[startsWith(variables_requested,"dep_")]

expand_dep_variable <- function(x) {
  suffix <- sub("^dep_","",x)
  paste0(dep_fluxes,"_",suffix)
}

variables_for_resolution <- variables_requested[!startsWith(variables_requested,"dep_")]

if (length(dep_variables) > 0) {
  for (x in dep_variables) variables_for_resolution <- c(variables_for_resolution,expand_dep_variable(x))
}

variables_for_resolution <- unique(variables_for_resolution)

###################################
### BUILD RESOLVED RESULT TABLE ###
###################################
make_result <- function(logical_name,idx,grib_column,table) {

  if (length(idx) == 0) return(data.frame())

  grib_codes <- table[[grib_column]][idx]

  valid <- !is.na(grib_codes) & grib_codes != ""
  idx <- idx[valid]
  grib_codes <- grib_codes[valid]

  if (length(idx) == 0) {
    stop("Variable '",logical_name,"' has no available values in column '",grib_column,"'.")
  }

  data.frame(
    logical_name = rep(logical_name,length(idx)),
    csv_name     = table$name[idx],
    grib_column  = rep(grib_column,length(idx)),
    grib         = grib_codes,
    stringsAsFactors=FALSE
  )
}

###############
### HELPERS ###
###############
get_variable_prefix <- function(logical_name) {
  if (startsWith(logical_name,"mss_from_mr_")) return("mss_from_mr")
  sub("_.*$","",logical_name)
}

get_variable_suffix <- function(logical_name) {
  if (startsWith(logical_name,"mss_from_mr_")) return(sub("^mss_from_mr_","",logical_name))
  sub("^[^_]+_","",logical_name)
}

###############################
### HAM VARIABLE RESOLUTION ###
###############################
resolve_HAM_variable <- function(logical_name) {
  prefix <- get_variable_prefix(logical_name)
  if (!prefix %in% c(names(variable_columns),"mss_from_mr")) { stop("Unsupported variable prefix '",prefix,"' in ",logical_name) }
  suffix <- get_variable_suffix(logical_name)
  grib_column <- if (prefix == "mss_from_mr") "grib" else variable_columns[[prefix]]

  ############################
  ### 1. EXACT MODE/SPECIES ###
  ############################
  ### Example: ddp_ss_cs -> SS_CS
  if (grepl(paste0("_(",paste(ham_modes,collapse="|"),")$"),suffix)) {

    idx <- which(ham_names == suffix)

    if (length(idx) == 0) stop("HAM variable '",logical_name,"' not found in HAM GRIB table.")
    if (length(idx) > 1) stop("HAM variable '",logical_name,"' matches multiple HAM rows.")

    return(make_result(logical_name,idx,grib_column,grib_HAM))
  }

  ########################
  ### 2. SOLUBILITY SUM ###
  ########################
  ### Example: ddp_soluble
  if (suffix == "soluble") {

    pattern <- paste0("_(",paste(ham_soluble_modes,collapse="|"),")$")
    idx <- which(ham_component_rows & grepl(pattern,ham_names))

    return(make_result(logical_name,idx,grib_column,grib_HAM))
  }

  ### Example: ddp_insoluble
  if (suffix == "insoluble") {

    pattern <- paste0("_(",paste(ham_insoluble_modes,collapse="|"),")$")
    idx <- which(ham_component_rows & grepl(pattern,ham_names))

    return(make_result(logical_name,idx,grib_column,grib_HAM))
  }

  ###################
  ### 3. MODE SUM ###
  ###################
  ### Example: ddp_as -> every species in AS
  if (suffix %in% ham_modes) {

    idx <- which(
      ham_component_rows &
      grepl(paste0("_",suffix,"$"),ham_names)
    )

    return(make_result(logical_name,idx,grib_column,grib_HAM))
  }

  ######################
  ### 4. SPECIES SUM ###
  ######################
  ### Example: ddp_ss -> SS_AS + SS_CS
  species_pattern <- paste0(
    "^",suffix,"_(",
    paste(ham_modes,collapse="|"),
    ")$"
  )

  idx <- which(
    ham_component_rows &
    grepl(species_pattern,ham_names)
  )

  if (length(idx) > 0) {
    return(make_result(logical_name,idx,grib_column,grib_HAM))
  }

  #########################
  ### 5. EXACT FALLBACK ###
  #########################
  ### Allows non-mode fields such as H2OPART if explicitly requested.
  idx <- which(ham_names == suffix)

  if (length(idx) == 1) {
    return(make_result(logical_name,idx,grib_column,grib_HAM))
  }

  stop("HAM variable '",logical_name,"' could not be resolved.")
}

###############################
### AER VARIABLE RESOLUTION ###
###############################

aer_names <- normalize_csv_name(grib_AER$name)

### Common logical species names.
### These are also the species which can be compared HAM vs AER.
AER_species_aliases <- list(
  ss  = c("SSS","SSM","SSL"),
  du  = c("DUS","DUM","DUL"),
  pom = c("OMHPHIL","OMHPHOB"),
  bc  = c("BCHPHIL","BCHPHOB"),
  so4 = c("SU"),
  ni  = c("NIF","NIC"),
  am  = c("AM"),
  soa = c("SOA1","SOA2")
)

common_HAM_AER_species <- names(AER_species_aliases)

resolve_AER_variable <- function(logical_name) {
  prefix <- get_variable_prefix(logical_name)
  if (!prefix %in% c(names(variable_columns),"mss_from_mr")) { stop("Unsupported variable prefix '",prefix,"' in ",logical_name) }
  suffix <- get_variable_suffix(logical_name)
  grib_column <- if (prefix == "mss_from_mr") "grib" else variable_columns[[prefix]]

  ########################
  ### 1. SPECIES TOTAL ###
  ########################
  ### Examples:
  ### ddp_ss  -> SSS + SSM + SSL
  ### ddp_du  -> DUS + DUM + DUL
  ### ddp_pom -> OMHPHIL + OMHPHOB
  if (suffix %in% names(AER_species_aliases)) {

    csv_names <- AER_species_aliases[[suffix]]
    idx <- which(aer_names %in% normalize_csv_name(csv_names))

    if (length(idx) == 0) {
      stop("AER species variable '",logical_name,"' could not be resolved.")
    }

    return(make_result(logical_name,idx,grib_column,grib_AER))
  }

  #######################
  ### 2. EXACT TRACER ###
  #######################
  ### Examples:
  ### ddp_sss
  ### ddp_dum
  ### ddp_omhphil
  ### ddp_bchphob
  ### ddp_nic
  ### ddp_soa1
  ### ddp_vfa
  idx <- which(aer_names == suffix)

  if (length(idx) == 1) {
    return(make_result(logical_name,idx,grib_column,grib_AER))
  }

  if (length(idx) > 1) {
    stop("AER variable '",logical_name,"' matches multiple AER rows.")
  }

  stop("AER variable '",logical_name,"' could not be resolved.")
}

###################################
### HAM VS AER COMPATIBILITY ###
###################################
mixed_aerosol_schemes <- exptype1 != exptype2

if (mixed_aerosol_schemes) {

  for (logical_name in variables_for_resolution) {

    suffix <- get_variable_suffix(logical_name)

    if (!suffix %in% common_HAM_AER_species) {
      stop(
        "HAM vs AER comparisons support per-species variables only. ",
        "Variable '",logical_name,"' is not a common HAM/AER species. ",
        "Available species: ",paste(common_HAM_AER_species,collapse=", ")
      )
    }
  }
}

###################################
### RESOLVE REQUESTED VARIABLES ###
###################################
resolve_variables <- function(exptype,variables) {

  result <- data.frame()

  for (logical_name in variables) {

    if (exptype == "HAM") {
      temp <- resolve_HAM_variable(logical_name)
      result <- rbind(result,temp)
    }

    if (exptype == "AER") {
      temp <- resolve_AER_variable(logical_name)
      result <- rbind(result,temp)
    }
  }

  result
}

################################
### RESOLVE BOTH EXPERIMENTS ###
################################
variables_exp1 <- resolve_variables(exptype1,variables_for_resolution)
variables_exp2 <- resolve_variables(exptype2,variables_for_resolution)

###################
### INFORMATION ###
###################
message("---> ",expname1," (",exptype1,"): ",paste(unique(variables_exp1$logical_name),collapse=", "))
message("---> ",expname2," (",exptype2,"): ",paste(unique(variables_exp2$logical_name),collapse=", "))
if (length(dep_variables) > 0) message("---> Composite deposition plots: ",paste(dep_variables,collapse=", "))

if (nrow(variables_exp1) > 0) message("---> ",expname1," GRIBs: ",paste(unique(variables_exp1$grib),collapse="/"))
if (nrow(variables_exp2) > 0) message("---> ",expname2," GRIBs: ",paste(unique(variables_exp2$grib),collapse="/"))
