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

#####################################
### RELATIVE HUMIDITY DEFINITIONS ###
#####################################
### RH is a derived meteorological diagnostic and is intentionally kept
### outside the aerosol GRIB-table resolver. The selected L137 levels are
### nominal levels nearest the requested approximate heights.
rh_definitions <- data.frame(
  logical_name     = c("rh","rh_500m","rh_1000m","rh_1500m","rh_2000m","rh_3000m"),
  model_level      = c(137,124,118,114,110,105),
  approx_height_m  = c(10,501,987,1460,2081,3089),
  stringsAsFactors = FALSE
)

rh_supported <- rh_definitions$logical_name
get_rh_definition <- function(logical_name) {
  x <- rh_definitions[rh_definitions$logical_name == logical_name,,drop=FALSE]
  if (nrow(x) != 1) stop("Unsupported RH variable: ",logical_name)
  x
}


#################################
### PRECIPITATION DEFINITIONS ###
#################################
### Precipitation is a meteorological surface diagnostic and is intentionally
### kept outside the aerosol GRIB-table resolver.
precip_definitions <- data.frame(
  logical_name = c("precip_total","precip_convective","precip_large_scale"),
  grib         = c("228.128","143.128","142.128"),
  title        = c("Total precipitation","Convective precipitation","Large-scale precipitation"),
  stringsAsFactors = FALSE
)

precip_supported <- precip_definitions$logical_name
get_precip_definition <- function(logical_name) {
  x <- precip_definitions[precip_definitions$logical_name == logical_name,,drop=FALSE]
  if (nrow(x) != 1) stop("Unsupported precipitation variable: ",logical_name)
  x
}


##########################
### OPTICAL DEFINITIONS ###
##########################
optical_definitions <- data.frame(
  logical_name = c("aod550","aod865","ae550to865","aaod550","ssa550","mec550"),
  title        = c("Aerosol optical depth @ 550 nm","Aerosol optical depth @ 865 nm","Angstrom exponent 550-865 nm","Absorption aerosol optical depth @ 550 nm","Single-scattering albedo @ 550 nm","Mass extinction coefficient @ 550 nm"),
  stringsAsFactors = FALSE
)
optical_supported <- optical_definitions$logical_name
get_optical_definition <- function(logical_name) {
  x <- optical_definitions[optical_definitions$logical_name == logical_name,,drop=FALSE]
  if (nrow(x) != 1) stop("Unsupported optical variable: ",logical_name)
  x
}
get_optical_required_gribs <- function(logical_names) {
  gribs <- character(0)
  if (any(logical_names %in% c("aod550","ae550to865","mec550"))) gribs <- c(gribs,"207.210")
  if (any(logical_names %in% c("aod865","ae550to865"))) gribs <- c(gribs,"215.210")
  if (any(logical_names %in% "aaod550")) gribs <- c(gribs,"104.215")
  if (any(logical_names %in% "ssa550")) gribs <- c(gribs,"140.215")
  unique(gribs)
}

######################################
### EXPAND COMPOSITE DEP VARIABLES ###
######################################
variables_requested <- variables

invalid_rh <- variables_requested[startsWith(variables_requested,"rh") & !variables_requested %in% rh_supported]
if (length(invalid_rh) > 0) {
  stop("Unsupported RH variable(s): ",paste(invalid_rh,collapse=", "),
       ". Available RH variables: ",paste(rh_supported,collapse=", "))
}

rh_variables <- variables_requested[variables_requested %in% rh_supported]

invalid_precip <- variables_requested[startsWith(variables_requested,"precip") & !variables_requested %in% precip_supported]
if (length(invalid_precip) > 0) {
  stop("Unsupported precipitation variable(s): ",paste(invalid_precip,collapse=", "),
       ". Available precipitation variables: ",paste(precip_supported,collapse=", "))
}

precip_variables <- variables_requested[variables_requested %in% precip_supported]
optical_variables <- variables_requested[variables_requested %in% optical_supported]

special_variables <- unique(c(rh_variables,precip_variables,optical_variables))
aerosol_variables_requested <- variables_requested[!variables_requested %in% special_variables]

dep_fluxes <- c("ddp","sdm","wdl","wdc","ngt")
dep_variables <- aerosol_variables_requested[startsWith(aerosol_variables_requested,"dep_")]

expand_dep_variable <- function(x) {
  suffix <- sub("^dep_","",x)
  paste0(dep_fluxes,"_",suffix)
}

variables_for_resolution <- aerosol_variables_requested[!startsWith(aerosol_variables_requested,"dep_")]

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
    massdiag_name = if ("long_name" %in% names(table)) table$long_name[idx] else table$name[idx],
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

  result <- data.frame(
    logical_name=character(),
    csv_name=character(),
    massdiag_name=character(),
    grib_column=character(),
    grib=character(),
    stringsAsFactors=FALSE
  )

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
message("---> ",expname1," (",exptype1,"): ",if (nrow(variables_exp1) > 0) paste(unique(variables_exp1$logical_name),collapse=", ") else "no aerosol diagnostics")
message("---> ",expname2," (",exptype2,"): ",if (nrow(variables_exp2) > 0) paste(unique(variables_exp2$logical_name),collapse=", ") else "no aerosol diagnostics")
if (length(dep_variables) > 0) message("---> Composite deposition plots: ",paste(dep_variables,collapse=", "))
if (length(rh_variables) > 0) {
  rh_info <- sapply(rh_variables,function(x) {
    z <- get_rh_definition(x)
    paste0(x," (ML",z$model_level,", ~",z$approx_height_m," m)")
  })
  message("---> Relative humidity: ",paste(rh_info,collapse=", "))
}
if (length(precip_variables) > 0) {
  precip_info <- sapply(precip_variables,function(x) {
    z <- get_precip_definition(x)
    paste0(x," (",z$grib,")")
  })
  message("---> Precipitation: ",paste(precip_info,collapse=", "))
}
if (length(optical_variables) > 0) {
  message("---> Optical properties: ",paste(optical_variables,collapse=", "))
}

if (nrow(variables_exp1) > 0) message("---> ",expname1," GRIBs: ",paste(unique(variables_exp1$grib),collapse="/"))
if (nrow(variables_exp2) > 0) message("---> ",expname2," GRIBs: ",paste(unique(variables_exp2$grib),collapse="/"))

###############
### REGIONS ###
###############
regions <- list(
  global    = c(-180,180,-90,90),
  n_america = c(-158,-50,10,85),
  us        = c(-130,-60,32,48),
  s_america = c(-90,-30,-60,18),
  africa    = c(-20,81,-40,38),
  China     = c(100,135,22,45),
  India     = c(35,95,5,70),
  WUS       = c(-130,-100,32,48),
  EUS       = c(-100,-60,32,48),
  europe    = c(-18,40,30,70),
  desert_aeronet = c(-180,180,-90,90),
  ocean_aeronet  = c(-180,180,-90,90),
  se_asia   = c(65,180,-23,50)
)
