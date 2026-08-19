### EXPIREMENTS settings
expname1  <- "b30b"      # 4 letter of the experiment name
exptype1  <- "HAM_NI_CS" # Either AER or HAM or HAM_NI_CS
expclass1 <- "nl"        # Class from MARS system
expname2  <- "b30j"      # 4 letter of the experiment name
exptype2  <- "HAM_NI_CS" # Either AER or HAM or HAM_NI_CS
expclass2 <- "nl"        # Class from MARS system
sDate     <- "20181203"  # Starting date in the form of YYYYMMDD
eDate     <- "20181203"  # Ending date in the form of YYYYMMDD
runtype   <- "download"  # Could be either "download" or "plot"
### DOWNLOAD settings
NumberOfDownloadJobs <- 1 # This is per experiment per kind (pl and sfc) so 1 = 4 parallel Jobs! Max is 20 under hpc-login
### PLOT settings
plot_title     <- "Aerosol comparison"
variable_title <- c("AOD 550nm","AE 550-865nm","AAOD 550nm","SSA 550nm","MEC 550nm (AOD/Burden)","Total Column Burden","DU Column Burden","SS Column Burden","POM Column Burden","BC Column Burden","SO4 Column Burden","AM Column Burden","NI Column Burden") # "Mass Extinction Coefficient"
variable_name  <- c("aod550","ae550to865","aodabs550","ssa550","mec550","BU_AERO","BU_DU","BU_SS","BU_POM","BU_BC","BU_SO4","BU_AM","BU_NI")
variable_units <- c("Unitless","Unitless","Unitless","Unitless","kg^-1 m^2","g m^-2","g m^-2","g m^-2","g m^-2","g m^-2","g m^-2","g m^-2","g m^-2")
variable_mul   <- c(1,1,1,1,1,10^3,10^3,10^3,10^3,10^3,10^3,10^3,10^3)
field_pallete_name <- c("AOD_white","AE","AAOD","SSA","AOD_white","TROPOMI_NEW","TROPOMI_NEW","TROPOMI_NEW","TROPOMI_NEW","TROPOMI_NEW","TROPOMI_NEW","TROPOMI_NEW","TROPOMI_NEW")
field_pallete_name_diff <- c("MNMB","MNMB","AAOD_diff","SSA_diff","MNMB","MNMB","MNMB","MNMB","MNMB","MNMB","MNMB","MNMB","MNMB")
field_breaks <- list(
  c(0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.50,0.60,0.80,1.00),
  seq(0.0,2.0,0.2),
  c(0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.50,0.60,0.80,1.00)*0.1,
  seq(0.8,1.0,0.02),
  seq(0,6000,500),
  c(0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.50,0.60,0.80,1.00),
  c(0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.50,0.60,0.80,1.00),
  c(0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.50,0.60,0.80,1.00),
  c(0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.50,0.60,0.80,1.00)*0.1,
  c(0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.50,0.60,0.80,1.00)*0.01,
  c(0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.50,0.60,0.80,1.00)*0.1,
  c(0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.50,0.60,0.80,1.00)*0.01,
  c(0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.50,0.60,0.80,1.00)*0.01
)
field_breaks_diff <- list(
  c( rev(c(0.05,0.10,0.20,0.30,0.40,0.50,0.60)*-1),c(0.05,0.10,0.20,0.30,0.40,0.50,0.60) ),
  seq(-0.9,0.9,0.2),
  c( rev(c(0.05,0.10,0.20,0.30,0.40,0.50,0.60)*-1),c(0.05,0.10,0.20,0.30,0.40,0.50,0.60) )*0.1,
  seq(-0.09,0.09,0.02),
  seq(-3000,3000,500),
  c( rev(c(0.05,0.10,0.20,0.30,0.40,0.50,0.60)*-1),c(0.05,0.10,0.20,0.30,0.40,0.50,0.60) ),
  c( rev(c(0.05,0.10,0.20,0.30,0.40,0.50,0.60)*-1),c(0.05,0.10,0.20,0.30,0.40,0.50,0.60) ),
  c( rev(c(0.05,0.10,0.20,0.30,0.40,0.50,0.60)*-1),c(0.05,0.10,0.20,0.30,0.40,0.50,0.60) ),
  c( rev(c(0.05,0.10,0.20,0.30,0.40,0.50,0.60)*-1),c(0.05,0.10,0.20,0.30,0.40,0.50,0.60) )*0.1,
  c( rev(c(0.05,0.10,0.20,0.30,0.40,0.50,0.60)*-1),c(0.05,0.10,0.20,0.30,0.40,0.50,0.60) )*0.01,
  c( rev(c(0.05,0.10,0.20,0.30,0.40,0.50,0.60)*-1),c(0.05,0.10,0.20,0.30,0.40,0.50,0.60) )*0.1,
  c( rev(c(0.05,0.10,0.20,0.30,0.40,0.50,0.60)*-1),c(0.05,0.10,0.20,0.30,0.40,0.50,0.60) )*0.01,
  c( rev(c(0.05,0.10,0.20,0.30,0.40,0.50,0.60)*-1),c(0.05,0.10,0.20,0.30,0.40,0.50,0.60) )*0.01
)

