### EXPERIMENT SETTINGS
expname1  <- "b30r"
exptype1  <- "HAM"
expclass1 <- "nl"

expname2  <- "b30r"
exptype2  <- "HAM"
expclass2 <- "nl"

sDate   <- "20181201"
eDate   <- "20181229"
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
  "mss_am_as",
  "mss_ni_as",
  "mss_so4_as",
  "mss_bc_as",
  "mss_pom_as",
  "mss_ss_as",
  "mss_du_as",
  "mss_soa_ns",
  "mss_soa_ks",
  "mss_soa_as",
  "mss_soa_cs",
  "mss_soa_ki",
  "mss_bc_ki",
  "mss_pom_ki",
  "mss_du_ai",
  "mss_so4_ks",
  "mss_bc_ks",
  "mss_pom_ks",
  "mss_du_ci",
  "mss_so4_cs",
  "mss_bc_cs",
  "mss_pom_cs",
  "mss_ss_cs",
  "mss_du_cs",
  "mss_so4_ns",
  "mss_ni_cs",

  "dep_am_as",
  "dep_ni_as",
  "dep_so4_as",
  "dep_bc_as",
  "dep_pom_as",
  "dep_ss_as",
  "dep_du_as",
  "dep_soa_ns",
  "dep_soa_ks",
  "dep_soa_as",
  "dep_soa_cs",
  "dep_soa_ki",
  "dep_bc_ki",
  "dep_pom_ki",
  "dep_du_ai",
  "dep_so4_ks",
  "dep_bc_ks",
  "dep_pom_ks",
  "dep_du_ci",
  "dep_so4_cs",
  "dep_bc_cs",
  "dep_pom_cs",
  "dep_ss_cs",
  "dep_du_cs",
  "dep_so4_ns",
  "dep_ni_cs"
)

