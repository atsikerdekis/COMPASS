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
  od  = "gribod",
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
  logical_name <- strip_time_aggregation(logical_name)
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
  logical_name <- strip_time_aggregation(logical_name)
  x <- precip_definitions[precip_definitions$logical_name == logical_name,,drop=FALSE]
  if (nrow(x) != 1) stop("Unsupported precipitation variable: ",logical_name)
  x
}


##########################
### OPTICAL DEFINITIONS ###
##########################
optical_definitions <- data.frame(
  logical_name = c("aod550","aod865","ae550to865","aaod550","ssa550","mec550","od_sum_species","aod550_species"),
  title        = c("Aerosol optical depth @ 550 nm","Aerosol optical depth @ 865 nm","Angstrom exponent 550-865 nm","Absorption aerosol optical depth @ 550 nm","Single-scattering albedo @ 550 nm","Mass extinction coefficient @ 550 nm","Sum of species optical depth @ 550 nm","AOD550 with species time series"),
  stringsAsFactors = FALSE
)
optical_supported <- optical_definitions$logical_name
get_optical_definition <- function(logical_name) {
  logical_name <- strip_time_aggregation(logical_name)
  x <- optical_definitions[optical_definitions$logical_name == logical_name,,drop=FALSE]
  if (nrow(x) != 1) stop("Unsupported optical variable: ",logical_name)
  x
}

#############################################
### MULTI-PANEL SPECIES MAP DEFINITIONS ###
#############################################
panel_species_definitions <- data.frame(
  logical_name = c("aod_per_species","aodratio_per_species",
                   "mass_per_species","massratio_per_species",
                   "mec_per_species"),
  title = c("AOD per species",
            "Species AOD / total AOD550",
            "Mass burden per species",
            "Species mass / total mass burden",
            "Mass extinction coefficient per species"),
  stringsAsFactors = FALSE
)
panel_species_supported <- panel_species_definitions$logical_name
get_panel_species_definition <- function(logical_name) {
  logical_name <- tolower(strip_time_aggregation(logical_name))
  x <- panel_species_definitions[panel_species_definitions$logical_name == logical_name,,drop=FALSE]
  if (nrow(x) != 1) stop("Unsupported multi-panel species variable: ",logical_name)
  x
}

#############################################
### PANEL SATELLITE VALIDATION DEFINITIONS ###
#############################################
### Optional validation row used by the multi-panel species figures.
### Defaults reproduce the MODIS/VIIRS ensemble used in the
### SatelliteEvaluationTool workflow, while allowing config.R to override them.
if (!exists("panel_validation_variable")) panel_validation_variable <- "AOD550"
if (!exists("panel_satellites")) panel_satellites <- c("MOD","MYD","VIIRS_SNPP")
if (!exists("panel_satellite_path"))
  panel_satellite_path <- "/perm/nktt/SatelliteEvaluationTool/data/output/nc/3H/"

panel_validation_definitions <- data.frame(
  variable       = c("AOD550","AE550to860"),
  model_variable = c("aod550","ae550to865"),
  aod_filter     = c(0.0,0.2),
  title          = c("AOD550","AE550to860"),
  stringsAsFactors = FALSE
)

panel_satellite_definitions <- data.frame(
  sensor = c("MOD","MYD","VIIRS_SNPP"),
  file_prefix = c(
    "MOD04_L2_AOD_0.4x0.4_",
    "MYD04_L2_AOD_0.4x0.4_",
    "VIIRS_SNPP_L2_AOD_0.4x0.4_"
  ),
  stringsAsFactors = FALSE
)

get_panel_validation_definition <- function(variable=panel_validation_variable) {
  x <- panel_validation_definitions[panel_validation_definitions$variable == variable,,drop=FALSE]
  if (nrow(x) != 1)
    stop("Unsupported panel_validation_variable: ",variable,
         ". Available: ",paste(panel_validation_definitions$variable,collapse=", "))
  x
}

get_panel_satellite_definition <- function(sensor) {
  x <- panel_satellite_definitions[panel_satellite_definitions$sensor == sensor,,drop=FALSE]
  if (nrow(x) != 1)
    stop("Unsupported panel satellite: ",sensor,
         ". Available: ",paste(panel_satellite_definitions$sensor,collapse=", "))
  x
}

get_panel_validation_required_gribs <- function(variable=panel_validation_variable) {
  def <- get_panel_validation_definition(variable)
  if (def$model_variable == "aod550") return("207.210")
  if (def$model_variable == "ae550to865") return(c("207.210","215.210"))
  stop("No GRIB dependency mapping for panel validation variable: ",variable)
}

### Fail early if the config contains unsupported validation settings.
get_panel_validation_definition(panel_validation_variable)
invisible(lapply(panel_satellites,get_panel_satellite_definition))

get_optical_required_gribs <- function(logical_names) {
  logical_names <- strip_time_aggregation(logical_names)
  gribs <- character(0)

  ### AOD550 is also needed for representative SSA and MEC maps.
  if (any(logical_names %in% c("aod550","ae550to865","ssa550","mec550","aod550_species")))
    gribs <- c(gribs,"207.210")

  if (any(logical_names %in% c("aod865","ae550to865")))
    gribs <- c(gribs,"215.210")

  ### AAOD550 is needed for representative SSA maps.
  if (any(logical_names %in% c("aaod550","ssa550")))
    gribs <- c(gribs,"104.215")

  ### Keep direct SSA for native time series and diurnal cycle.
  if (any(logical_names %in% "ssa550"))
    gribs <- c(gribs,"140.215")

  unique(gribs)
}


###########################################
### HAM DEDICATED SPECIES AOD @ 550 NM ###
###########################################
### Dedicated HAM7 species-total AOD diagnostics. These are preferred for
### od_<species> over summing tracer/mode gribod fields from the HAM table.
### MARS/COMPASS notation PARAM.TABLE is used here:
### e.g. ECMWF paramId 210208 -> 208.210.
ham_species_od_definitions <- data.frame(
  suffix = c("ss",     "du",     "pom",    "bc",     "so4",    "ni",     "am"),
  grib   = c("208.210","209.210","210.210","211.210","212.210","250.210","251.210"),
  short_name = c("ssaod550","duaod550","omaod550","bcaod550","suaod550","niaod550","amaod550"),
  stringsAsFactors = FALSE
)

get_ham_species_od_definition <- function(suffix) {
  x <- ham_species_od_definitions[ham_species_od_definitions$suffix == suffix,,drop=FALSE]
  if (nrow(x) != 1) return(NULL)
  x
}

### Species used by the derived AOD panel diagnostics.
### This helper must live in 01.init.R because 02.start.R needs it in download mode,
### before 04.preprocess.R is sourced.
get_od_species_suffixes <- function(include_soa=TRUE) {
  x <- c("ss","du","pom","bc","so4","ni","am")
  if (include_soa) x <- c(x,"soa")
  x
}


############################
### HAM WATER DEFINITIONS ###
############################
wat_definitions <- data.frame(
  logical_name = c("wat_ks","wat_as","wat_cs","wat"),
  levelist     = c("2","3","4","2/3/4"),
  title        = c("Water AOD KS","Water AOD AS","Water AOD CS","Water AOD total"),
  stringsAsFactors = FALSE
)
wat_supported <- wat_definitions$logical_name
get_wat_definition <- function(logical_name) {
  logical_name <- strip_time_aggregation(logical_name)
  x <- wat_definitions[wat_definitions$logical_name == logical_name,,drop=FALSE]
  if (nrow(x) != 1) stop("Unsupported water variable: ",logical_name)
  x
}


################################
### TOTAL-COLUMN DEFINITIONS ###
################################
### Gas total-column diagnostics archived on the surface stream.
column_definitions <- data.frame(
  logical_name = c("mss_hno3","mss_nh3"),
  grib         = c("6.218","19.218"),
  title        = c("Total column Nitric Acid (HNO3)",
                   "Total column Ammonia (NH3)"),
  units        = c("kg m^-2","kg m^-2"),
  stringsAsFactors = FALSE
)

column_supported <- column_definitions$logical_name

get_column_definition <- function(logical_name) {
  logical_name <- strip_time_aggregation(logical_name)
  x <- column_definitions[column_definitions$logical_name == logical_name,,drop=FALSE]
  if (nrow(x) != 1) stop("Unsupported total-column variable: ",logical_name)
  x
}


##################################
### TIME-RESOLUTION MODIFIERS ###
##################################
### These suffixes change only the plotted time-series resolution.
### Spatial maps continue to use all native 3-hourly fields.
strip_time_aggregation <- function(x) {
  sub("_(daily|monthly)$","",x)
}

get_time_aggregation <- function(x) {
  out <- rep("3hourly",length(x))
  out[grepl("_daily$",x)] <- "daily"
  out[grepl("_monthly$",x)] <- "monthly"
  out
}

normalize_special_variable_names <- function(x) {
  agg <- ifelse(grepl("_daily$",x), "_daily",
                ifelse(grepl("_monthly$",x), "_monthly", ""))
  base <- strip_time_aggregation(x)
  base_low <- tolower(base)
  if (base_low %in% panel_species_supported) return(paste0(base_low,agg))
  x
}

######################################
### EXPAND COMPOSITE DEP VARIABLES ###
######################################
variables_requested <- vapply(variables,normalize_special_variable_names,character(1))
variable_bases <- strip_time_aggregation(variables_requested)
variable_bases_lower <- tolower(variable_bases)

invalid_rh <- variables_requested[
  startsWith(variable_bases,"rh") & !variable_bases %in% rh_supported
]
if (length(invalid_rh) > 0) {
  stop("Unsupported RH variable(s): ",paste(invalid_rh,collapse=", "),
       ". Available RH variables: ",paste(rh_supported,collapse=", "))
}

rh_variables <- variables_requested[variable_bases %in% rh_supported]

invalid_wat <- variables_requested[
  startsWith(variable_bases,"wat") & !variable_bases %in% wat_supported
]
if (length(invalid_wat) > 0) {
  stop("Unsupported water variable(s): ",paste(invalid_wat,collapse=", "),
       ". Available water variables: ",paste(wat_supported,collapse=", "))
}
wat_variables <- variables_requested[variable_bases %in% wat_supported]

invalid_precip <- variables_requested[
  startsWith(variable_bases,"precip") & !variable_bases %in% precip_supported
]
if (length(invalid_precip) > 0) {
  stop("Unsupported precipitation variable(s): ",paste(invalid_precip,collapse=", "),
       ". Available precipitation variables: ",paste(precip_supported,collapse=", "))
}

precip_variables <- variables_requested[variable_bases %in% precip_supported]
optical_variables <- variables_requested[variable_bases %in% optical_supported]
panel_species_variables <- variables_requested[variable_bases_lower %in% panel_species_supported]
column_variables <- variables_requested[variable_bases %in% column_supported]

### Standalone deposition-lifetime diagnostics.
### Example: lifetime_ni_as, lifetime_ni, lifetime_as.
lifetime_variables <- variables_requested[startsWith(variable_bases,"lifetime_")]

special_variables <- unique(c(
  rh_variables,wat_variables,precip_variables,optical_variables,panel_species_variables,
  column_variables,lifetime_variables
))
aerosol_variables_requested <- variables_requested[!variables_requested %in% special_variables]

dep_fluxes <- c("ddp","sdm","wdl","wdc","ngt")
dep_variables <- aerosol_variables_requested[
  startsWith(strip_time_aggregation(aerosol_variables_requested),"dep_")
]

expand_dep_variable <- function(x) {
  suffix <- sub("^dep_","",strip_time_aggregation(x))
  paste0(dep_fluxes,"_",suffix)
}

expand_lifetime_variable <- function(x) {
  suffix <- sub("^lifetime_","",strip_time_aggregation(x))
  c(paste0("mss_",suffix),paste0(dep_fluxes,"_",suffix))
}

variables_for_resolution <- aerosol_variables_requested[
  !startsWith(strip_time_aggregation(aerosol_variables_requested),"dep_")
]

if (length(dep_variables) > 0) {
  for (x in dep_variables)
    variables_for_resolution <- c(variables_for_resolution,expand_dep_variable(x))
}

if (length(lifetime_variables) > 0) {
  for (x in lifetime_variables)
    variables_for_resolution <- c(variables_for_resolution,expand_lifetime_variable(x))
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
  logical_name <- strip_time_aggregation(logical_name)
  if (startsWith(logical_name,"mss_from_mr_")) return("mss_from_mr")
  sub("_.*$","",logical_name)
}

get_variable_suffix <- function(logical_name) {
  logical_name <- strip_time_aggregation(logical_name)
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

  ########################################
  ### DEDICATED SPECIES-TOTAL HAM AOD ###
  ########################################
  ### For od_ss, od_du, od_pom, od_bc, od_so4, od_ni and od_am, use the
  ### dedicated CAMS/HAM 550-nm species AOD diagnostic instead of summing
  ### modal/tracer optical-depth fields from gribod.
  if (prefix == "od") {
    oddef <- get_ham_species_od_definition(suffix)
    if (!is.null(oddef)) {
      return(data.frame(
        logical_name  = logical_name,
        csv_name      = toupper(suffix),
        massdiag_name = oddef$short_name,
        grib_column   = "gribod",
        grib          = oddef$grib,
        stringsAsFactors = FALSE
      ))
    }
  }

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
if (length(wat_variables) > 0) {
  wat_info <- sapply(wat_variables,function(x) {
    z <- get_wat_definition(x)
    paste0(x," (210022, levels ",z$levelist,")")
  })
  message("---> Water AOD: ",paste(wat_info,collapse=", "))
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
if (length(panel_species_variables) > 0) {
  message("---> Multi-panel species maps: ",paste(panel_species_variables,collapse=", "))
  message("---> Panel validation: ",panel_validation_variable,
          " using ",paste(panel_satellites,collapse="+"),
          " from ",panel_satellite_path)
}
if (length(column_variables) > 0) {
  column_info <- sapply(column_variables,function(x) {
    z <- get_column_definition(x)
    paste0(x," (",z$grib,")")
  })
  message("---> Total-column gases: ",paste(column_info,collapse=", "))
}
if (length(lifetime_variables) > 0) {
  message("---> Deposition lifetime: ",paste(lifetime_variables,collapse=", "))
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
  se_asia   = c(65,180,-23,50),

  ## bias-focused rectangles for 4x3 validation maps
  bias_india              = c(68, 90, 8, 31),
  bias_china              = c(102, 124, 22, 42),
  bias_europe             = c(-8, 25, 43, 58),
  bias_gulf_of_guinea     = c(-18, 12, -2, 12),
  bias_mexico_plume       = c(-118, -96, 14, 30),
  bias_north_south_america= c(-82, -58, -5, 13)
)
