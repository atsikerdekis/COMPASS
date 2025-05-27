### COMPRESS IMAGE
compress <- function(file_in, file_out) {system(paste0("/usr/bin/convert ",file_in," +dither -colors 256 ",file_out))}
