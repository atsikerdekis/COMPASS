Normalization <- function(x,ymin,ymax) {
  xmax <- max(x,na.rm=T)
  xmin <- min(x,na.rm=T)
  y <- ymin + (ymax-ymin) * ((x-xmin)/(xmax-xmin)) 
  y
}

