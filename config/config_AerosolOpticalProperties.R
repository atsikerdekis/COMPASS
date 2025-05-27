plot_title     <- "Aerosol Optical Properties"
variable_title <- c("AOD 550nm","AE 550-865nm","AAOD 550nm","SSA 550nm","MEC 550nm (AOD/Burden)") # "Mass Extinction Coefficient"
variable_name  <- c("aod550","ae550to865","aodabs550","ssa550","mec550")
variable_units <- c("Unitless","Unitless","Unitless","Unitless","kg^-1 m^2")
variable_mul   <- c(1,1,1,1,1)
field_pallete_name <- c("AOD_white","AE","AAOD","SSA","AOD_white")
field_pallete_name_diff <- c("MNMB","MNMB","AAOD_diff","SSA_diff","MNMB")
field_breaks <- list(
  c(0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.50,0.60,0.80,1.00),
  seq(0.0,2.0,0.2),
  c(0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.50,0.60,0.80,1.00)*0.1,
  seq(0.8,1.0,0.02),
  seq(0,6000,500)
)
field_breaks_diff <- list(
  c( rev(c(0.05,0.10,0.20,0.30,0.40,0.50,0.60)*-1),c(0.05,0.10,0.20,0.30,0.40,0.50,0.60) ),
  seq(-0.9,0.9,0.2),
  c( rev(c(0.05,0.10,0.20,0.30,0.40,0.50,0.60)*-1),c(0.05,0.10,0.20,0.30,0.40,0.50,0.60) )*0.1,
  seq(-0.09,0.09,0.02),
  seq(-3000,3000,500)
)

