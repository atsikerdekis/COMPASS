### EXPERIMENT SETTINGS
expname1  <- "b30q"
exptype1  <- "HAM"
expclass1 <- "nl"

expname2  <- "b30b"
exptype2  <- "HAM"
expclass2 <- "nl"

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

  ##################
  ### PER TRACER ###
  ##################
  "ddp_am_as", "sdm_am_as", "wdl_am_as", "wdc_am_as", "mss_am_as", "ngt_am_as",
  "ddp_ni_as", "sdm_ni_as", "wdl_ni_as", "wdc_ni_as", "mss_ni_as", "ngt_ni_as",
  "ddp_so4_as", "sdm_so4_as", "wdl_so4_as", "wdc_so4_as", "mss_so4_as", "ngt_so4_as",
  "ddp_bc_as", "sdm_bc_as", "wdl_bc_as", "wdc_bc_as", "mss_bc_as", "ngt_bc_as",
  "ddp_pom_as", "sdm_pom_as", "wdl_pom_as", "wdc_pom_as", "mss_pom_as", "ngt_pom_as",
  "ddp_ss_as", "sdm_ss_as", "wdl_ss_as", "wdc_ss_as", "mss_ss_as", "ngt_ss_as",
  "ddp_du_as", "sdm_du_as", "wdl_du_as", "wdc_du_as", "mss_du_as", "ngt_du_as",
  "ddp_soa_as", "sdm_soa_as", "wdl_soa_as", "wdc_soa_as", "mss_soa_as", "ngt_soa_as",
  "ddp_du_ai", "sdm_du_ai", "wdl_du_ai", "wdc_du_ai", "mss_du_ai", "ngt_du_ai",
  "ddp_so4_cs", "sdm_so4_cs", "wdl_so4_cs", "wdc_so4_cs", "mss_so4_cs", "ngt_so4_cs",
  "ddp_bc_cs", "sdm_bc_cs", "wdl_bc_cs", "wdc_bc_cs", "mss_bc_cs", "ngt_bc_cs",
  "ddp_pom_cs", "sdm_pom_cs", "wdl_pom_cs", "wdc_pom_cs", "mss_pom_cs", "ngt_pom_cs",
  "ddp_ss_cs", "sdm_ss_cs", "wdl_ss_cs", "wdc_ss_cs", "mss_ss_cs", "ngt_ss_cs",
  "ddp_du_cs", "sdm_du_cs", "wdl_du_cs", "wdc_du_cs", "mss_du_cs", "ngt_du_cs",
  "ddp_ni_cs", "sdm_ni_cs", "wdl_ni_cs", "wdc_ni_cs", "mss_ni_cs", "ngt_ni_cs",
  "ddp_soa_cs", "sdm_soa_cs", "wdl_soa_cs", "wdc_soa_cs", "mss_soa_cs", "ngt_soa_cs",
  "ddp_du_ci", "sdm_du_ci", "wdl_du_ci", "wdc_du_ci", "mss_du_ci", "ngt_du_ci",
  "ddp_bc_ki", "sdm_bc_ki", "wdl_bc_ki", "wdc_bc_ki", "mss_bc_ki", "ngt_bc_ki",
  "ddp_pom_ki", "sdm_pom_ki", "wdl_pom_ki", "wdc_pom_ki", "mss_pom_ki", "ngt_pom_ki",
  "ddp_soa_ki", "sdm_soa_ki", "wdl_soa_ki", "wdc_soa_ki", "mss_soa_ki", "ngt_soa_ki",
  "ddp_so4_ks", "sdm_so4_ks", "wdl_so4_ks", "wdc_so4_ks", "mss_so4_ks", "ngt_so4_ks",
  "ddp_bc_ks", "sdm_bc_ks", "wdl_bc_ks", "wdc_bc_ks", "mss_bc_ks", "ngt_bc_ks",
  "ddp_pom_ks", "sdm_pom_ks", "wdl_pom_ks", "wdc_pom_ks", "mss_pom_ks", "ngt_pom_ks",
  "ddp_soa_ks", "sdm_soa_ks", "wdl_soa_ks", "wdc_soa_ks", "mss_soa_ks", "ngt_soa_ks",
  "ddp_so4_ns", "sdm_so4_ns", "wdl_so4_ns", "wdc_so4_ns", "mss_so4_ns", "ngt_so4_ns",
  "ddp_soa_ns", "sdm_soa_ns", "wdl_soa_ns", "wdc_soa_ns", "mss_soa_ns", "ngt_soa_ns",

  ###################
  ### PER SPECIES ###
  ###################
  "ddp_am",  "sdm_am",  "wdl_am",  "wdc_am",  "mss_am",  "ngt_am",
  "ddp_ni",  "sdm_ni",  "wdl_ni",  "wdc_ni",  "mss_ni",  "ngt_ni",
  "ddp_so4", "sdm_so4", "wdl_so4", "wdc_so4", "mss_so4", "ngt_so4",
  "ddp_bc",  "sdm_bc",  "wdl_bc",  "wdc_bc",  "mss_bc",  "ngt_bc",
  "ddp_pom", "sdm_pom", "wdl_pom", "wdc_pom", "mss_pom", "ngt_pom",
  "ddp_ss",  "sdm_ss",  "wdl_ss",  "wdc_ss",  "mss_ss",  "ngt_ss",
  "ddp_du",  "sdm_du",  "wdl_du",  "wdc_du",  "mss_du",  "ngt_du",
  "ddp_soa", "sdm_soa", "wdl_soa", "wdc_soa", "mss_soa", "ngt_soa",

  ################
  ### PER MODE ###
  ################
  "ddp_ns", "sdm_ns", "wdl_ns", "wdc_ns", "mss_ns", "ngt_ns",
  "ddp_ks", "sdm_ks", "wdl_ks", "wdc_ks", "mss_ks", "ngt_ks",
  "ddp_as", "sdm_as", "wdl_as", "wdc_as", "mss_as", "ngt_as",
  "ddp_cs", "sdm_cs", "wdl_cs", "wdc_cs", "mss_cs", "ngt_cs",
  "ddp_ki", "sdm_ki", "wdl_ki", "wdc_ki", "mss_ki", "ngt_ki",
  "ddp_ai", "sdm_ai", "wdl_ai", "wdc_ai", "mss_ai", "ngt_ai",
  "ddp_ci", "sdm_ci", "wdl_ci", "wdc_ci", "mss_ci", "ngt_ci",

  ######################
  ### PER SOLUBILITY ###
  ######################
  "ddp_soluble",   "sdm_soluble",   "wdl_soluble",   "wdc_soluble",   "mss_soluble",   "ngt_soluble",
  "ddp_insoluble", "sdm_insoluble", "wdl_insoluble", "wdc_insoluble", "mss_insoluble", "ngt_insoluble"
)

