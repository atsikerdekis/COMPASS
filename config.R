### EXPERIMENT SETTINGS
### Experiment 1 (control/reference)
expname1  <- "b30r"
exptype1  <- "HAM"
expclass1 <- "nl"
### Experiment 2
expname2  <- "b30z"
exptype2  <- "HAM"
expclass2 <- "nl"
### Period
sDate   <- "20181201"
eDate   <- "20181231"
### Runtype
runtype <- "download" # plot | download

### DOWNLOAD SETTINGS
NumberOfDownloadJobs <- 1

### GRIB TABLES
grib_table_AER <- "config/bins_aerver8.csv"
grib_table_HAM <- "config/bins_hamm7ver4.0.csv"

### VARIABLES
#variables <- c(
#  "mss_ss","mss_du","mss_pom","mss_bc","mss_so4","mss_ni","mss_am",
#  "dep_ss","dep_du","dep_pom","dep_bc","dep_so4","dep_ni","dep_am"
#)

variables <- c(
  "mss_ss","mss_du","mss_pom","mss_bc","mss_so4","mss_ni","mss_am",
  "mss_from_mr_ss",
  "mss_from_mr_du",
  "mss_from_mr_pom",
  "mss_from_mr_bc",
  "mss_from_mr_so4",
  "mss_from_mr_ni",
  "mss_from_mr_am"
)

