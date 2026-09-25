### EXPERIMENT SETTINGS
### Experiment 1 (control/reference)
expname1  <- "b315"
exptype1  <- "HAM"
expclass1 <- "nl"
### Experiment 2
expname2  <- "jbjc"
exptype2  <- "HAM"
expclass2 <- "rd"
### Period
sDate   <- "20190101"
eDate   <- "20190131"
### Runtype
runtype <- "plot" # plot | download
### Other
massdiag_compare <- FALSE

### DOWNLOAD SETTINGS
NumberOfDownloadJobs <- 1

### GRIB TABLES
grib_table_AER <- "config/bins_aerver8.csv"
grib_table_HAM <- "config/bins_hamm7ver4.0.csv"

### REGION
region <- "global"

### VARIABLES
#variables <- c("precip_total","precip_convective","precip_large_scale")

#variables <- c(
#  "rh"        # ML137, ~surface / ~10 m
#  "rh_500m",   # ML124, ~500 m
#  "rh_1000m",  # ML118, ~1000 m
#  "rh_1500m",  # ML114, ~1500 m
#  "rh_2000m",  # ML110, ~2000 m
#  "rh_3000m"   # ML105, ~3000 m
#)

### VARIABLES
#variables <- c(
#  "mss_ni", "mss_su"
#)

### VARIABLES
#variables <- c(
#  "mss_ss","mss_du","mss_pom","mss_bc","mss_so4","mss_ni","mss_am",
#  "dep_ss","dep_du","dep_pom","dep_bc","dep_so4","dep_ni","dep_am"
#)

#variables <- c(
#  "mss_ss","mss_du","mss_pom","mss_bc","mss_so4","mss_ni","mss_am",
#  "mss_from_mr_ss","mss_from_mr_du","mss_from_mr_pom","mss_from_mr_bc","mss_from_mr_so4","mss_from_mr_ni","mss_from_mr_am"
#)

### mss and dep per HAM7 tracer
variables <- c(
  "mss_so4_as","mss_so4_ks","mss_so4_cs","mss_so4_ns",
  "mss_am_as",
  "mss_ni_as","mss_ni_cs",
#  "mss_bc_as","mss_bc_ki","mss_bc_ks","mss_bc_cs",
#  "mss_pom_as","mss_pom_ki","mss_pom_ks","mss_pom_cs",
#  "mss_soa_ns","mss_soa_ks","mss_soa_as","mss_soa_cs","mss_soa_ki",
#  "mss_ss_as","mss_ss_cs",
#  "mss_du_as","mss_du_ai","mss_du_ci","mss_du_cs",

  "dep_so4_as","dep_so4_ks","dep_so4_cs","dep_so4_ns",
  "dep_am_as",
  "dep_ni_as","dep_ni_cs"
#  "dep_bc_as","dep_bc_ki","dep_bc_ks","dep_bc_cs",
#  "dep_pom_as","dep_pom_ki","dep_pom_ks","dep_pom_cs",
#  "dep_soa_ns","dep_soa_ks","dep_soa_as","dep_soa_cs","dep_soa_ki",
#  "dep_ss_as","dep_ss_cs",
#  "dep_du_as","dep_du_ai","dep_du_ci","dep_du_cs"
)
