### EXPERIMENT SETTINGS
expname1  <- "b30j"
exptype1  <- "HAM"
expclass1 <- "nl"

expname2  <- "b30k"
exptype2  <- "HAM"
expclass2 <- "nl"

sDate   <- "20181201"
eDate   <- "20181203"
runtype <- "plot"

### DOWNLOAD SETTINGS
NumberOfDownloadJobs <- 1

### GRIB TABLES
grib_table_AER <- "config/bins_aerver8.csv"
grib_table_HAM <- "config/bins_hamm7ver4.0.csv"

### VARIABLES
variables <- c(
  "ddp_ss_cs", "wdl_ss_cs", "wdc_ss_cs", "sdm_ss_cs"
)


