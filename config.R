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
### VARIABLES
variables <- c(

  # Aitken soluble (AS)
  "ddp_am_as",  "sdm_am_as",  "wdl_am_as",  "wdc_am_as",
  "ddp_ni_as",  "sdm_ni_as",  "wdl_ni_as",  "wdc_ni_as",
  "ddp_so4_as", "sdm_so4_as", "wdl_so4_as", "wdc_so4_as",
  "ddp_bc_as",  "sdm_bc_as",  "wdl_bc_as",  "wdc_bc_as",
  "ddp_pom_as", "sdm_pom_as", "wdl_pom_as", "wdc_pom_as",
  "ddp_ss_as",  "sdm_ss_as",  "wdl_ss_as",  "wdc_ss_as",
  "ddp_du_as",  "sdm_du_as",  "wdl_du_as",  "wdc_du_as",
  "ddp_soa_as", "sdm_soa_as", "wdl_soa_as", "wdc_soa_as",

  # Aitken insoluble (AI)
  "ddp_du_ai",  "sdm_du_ai",  "wdl_du_ai",  "wdc_du_ai",

  # Coarse soluble (CS)
  "ddp_so4_cs", "sdm_so4_cs", "wdl_so4_cs", "wdc_so4_cs",
  "ddp_bc_cs",  "sdm_bc_cs",  "wdl_bc_cs",  "wdc_bc_cs",
  "ddp_pom_cs", "sdm_pom_cs", "wdl_pom_cs", "wdc_pom_cs",
  "ddp_ss_cs",  "sdm_ss_cs",  "wdl_ss_cs",  "wdc_ss_cs",
  "ddp_du_cs",  "sdm_du_cs",  "wdl_du_cs",  "wdc_du_cs",
  "ddp_ni_cs",  "sdm_ni_cs",  "wdl_ni_cs",  "wdc_ni_cs",
  "ddp_soa_cs", "sdm_soa_cs", "wdl_soa_cs", "wdc_soa_cs",

  # Coarse insoluble (CI)
  "ddp_du_ci",  "sdm_du_ci",  "wdl_du_ci",  "wdc_du_ci",

  # Aitken insoluble / KI
  "ddp_bc_ki",  "sdm_bc_ki",  "wdl_bc_ki",  "wdc_bc_ki",
  "ddp_pom_ki", "sdm_pom_ki", "wdl_pom_ki", "wdc_pom_ki",
  "ddp_soa_ki", "sdm_soa_ki", "wdl_soa_ki", "wdc_soa_ki",

  # Accumulation soluble / KS
  "ddp_so4_ks", "sdm_so4_ks", "wdl_so4_ks", "wdc_so4_ks",
  "ddp_bc_ks",  "sdm_bc_ks",  "wdl_bc_ks",  "wdc_bc_ks",
  "ddp_pom_ks", "sdm_pom_ks", "wdl_pom_ks", "wdc_pom_ks",
  "ddp_soa_ks", "sdm_soa_ks", "wdl_soa_ks", "wdc_soa_ks",

  # Nucleation soluble (NS)
  "ddp_so4_ns", "sdm_so4_ns", "wdl_so4_ns", "wdc_so4_ns",
  "ddp_soa_ns", "sdm_soa_ns", "wdl_soa_ns", "wdc_soa_ns"
)

