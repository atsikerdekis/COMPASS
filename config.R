### EXPERIMENT SETTINGS
### Experiment 1 (control/reference)
expname1  <- "iyfb"
exptype1  <- "AER"
expclass1 <- "rd"
### Experiment 2
expname2  <- "jayq"
exptype2  <- "HAM"
expclass2 <- "nl"
### Period
sDate   <- "20181201"
eDate   <- "20190131"
### Runtype
runtype <- "download" # plot | download
### Other
massdiag_compare <- FALSE

### DOWNLOAD SETTINGS
NumberOfDownloadJobs <- 1

### GRIB TABLES
grib_table_AER <- "config/bins_aerver8.csv"
grib_table_HAM <- "config/bins_hamm7ver4.0.csv"

### REGION
region <- "China"
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

### VARIABLES
variables <- c(
  "mss_ni"
)

### VARIABLES
#variables <- c(
#  "mss_ss","mss_du","mss_pom","mss_bc","mss_so4","mss_ni","mss_am",
#  "dep_ss","dep_du","dep_pom","dep_bc","dep_so4","dep_ni","dep_am"
#)

#variables <- c(
#  "mss_ss","mss_du","mss_pom","mss_bc","mss_so4","mss_ni","mss_am",
#  "mss_from_mr_ss","mss_from_mr_du","mss_from_mr_pom","mss_from_mr_bc","mss_from_mr_so4","mss_from_mr_ni","mss_from_mr_am"
#)

