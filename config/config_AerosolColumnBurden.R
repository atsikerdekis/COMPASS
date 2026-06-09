### EXP settings
expname1  <- "b2t5"     # 4 letter of the experiment name
exptype1  <- "HAM"      # Either AER or HAM
expclass1 <- "nl"       # Class from MARS system
expname2  <- "it41"     # 4 letter of the experiment name
exptype2  <- "AER"      # Either AER or HAM
expclass2 <- "rd"       # Class from MARS system
sDate     <- "20181201" # Starting date in the form of YYYYMMDD
eDate     <- "20181231" # Ending date in the form of YYYYMMDD
runtype   <- "download" # Could be either "download" or "plot"
### DOWNLOAD settings
NumberOfDownloadJobs <- 4 # This is per experiment!
### PLOT settings
plot_title     <- "Aerosol Column Burden"
variable_title <- c("Total Column Burden","DU Column Burden","SS Column Burden","POM Column Burden","BC Column Burden","SO4 Column Burden","AM Column Burden","NI Column Burden")
variable_name  <- c("BU_AERO","BU_DU","BU_SS","BU_POM","BU_BC","BU_SO4","BU_AM","BU_NI")
variable_units <- c("g m^-2","g m^-2","g m^-2","g m^-2","g m^-2","g m^-2","g m^-2","g m^-2")
variable_mul   <- c(10^3,10^3,10^3,10^3,10^3,10^3,10^3,10^3)
field_pallete_name <- c("TROPOMI_NEW","TROPOMI_NEW","TROPOMI_NEW","TROPOMI_NEW","TROPOMI_NEW","TROPOMI_NEW","TROPOMI_NEW","TROPOMI_NEW")
field_pallete_name_diff <- c("MNMB","MNMB","MNMB","MNMB","MNMB","MNMB","MNMB","MNMB")
field_breaks <- list(
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
  c( rev(c(0.05,0.10,0.20,0.30,0.40,0.50,0.60)*-1),c(0.05,0.10,0.20,0.30,0.40,0.50,0.60) ),
  c( rev(c(0.05,0.10,0.20,0.30,0.40,0.50,0.60)*-1),c(0.05,0.10,0.20,0.30,0.40,0.50,0.60) ),
  c( rev(c(0.05,0.10,0.20,0.30,0.40,0.50,0.60)*-1),c(0.05,0.10,0.20,0.30,0.40,0.50,0.60) )*0.1,
  c( rev(c(0.05,0.10,0.20,0.30,0.40,0.50,0.60)*-1),c(0.05,0.10,0.20,0.30,0.40,0.50,0.60) )*0.01,
  c( rev(c(0.05,0.10,0.20,0.30,0.40,0.50,0.60)*-1),c(0.05,0.10,0.20,0.30,0.40,0.50,0.60) )*0.1,
  c( rev(c(0.05,0.10,0.20,0.30,0.40,0.50,0.60)*-1),c(0.05,0.10,0.20,0.30,0.40,0.50,0.60) )*0.01,
  c( rev(c(0.05,0.10,0.20,0.30,0.40,0.50,0.60)*-1),c(0.05,0.10,0.20,0.30,0.40,0.50,0.60) )*0.01
)

