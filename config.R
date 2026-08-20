### EXPERIMENT SETTINGS
expname1  <- "b30q"
exptype1  <- "HAM"
expclass1 <- "nl"

expname2  <- "iyfb"
exptype2  <- "AER"
expclass2 <- "rd"

sDate   <- "20181201"
eDate   <- "20181203"
runtype <- "plot" # plot | download

### DOWNLOAD SETTINGS
NumberOfDownloadJobs <- 1

### GRIB TABLES
grib_table_AER <- "config/bins_aerver8.csv"
grib_table_HAM <- "config/bins_hamm7ver4.0.csv"

### VARIABLES
variables <- c(
  "mss_ss",
  "mss_du",
  "mss_pom",
  "mss_bc",
  "mss_so4",
  "mss_ni",
  "mss_am",
  "dep_ss",
  "dep_du",
  "dep_pom",
  "dep_bc",
  "dep_so4",
  "dep_ni",
  "dep_am"
)

