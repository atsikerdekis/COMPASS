### Thanos: Adapted from mapImage. Basically added the density and border=NA option to polygon()
mapShade <- function (longitude, latitude, z, zlim, zclip = FALSE, breaks, density=10,
          col, colormap, border = NA, lwd = par("lwd"), lty = par("lty"), 
          missingColor = NA, filledContour = FALSE, gridder = "binMean2D", 
          debug = getOption("oceDebug")) 
{
  if ("none" == .Projection()$type) 
    stop("must create a map first, with mapPlot()\n")
  breaksGiven <- !missing(breaks)
  zlimGiven <- !missing(zlim)
  colGiven <- !missing(col)
  oceDebug(debug, "mapImage(..., ", " missingColor=", missingColor, 
           ", ", " filledContour=", filledContour, ", ", " gridder=", 
           gridder, ", ", ", ...) {\n", sep = "", unindent = 1)
  if ("data" %in% slotNames(longitude)) {
    if (3 == sum(c("longitude", "latitude", "z") %in% names(longitude@data))) {
      z <- longitude@data$z
      latitude <- longitude@data$latitude
      longitude <- longitude@data$longitude
    }
  }
  else {
    names <- names(longitude)
    if ("x" %in% names && "y" %in% names && "z" %in% names) {
      z <- longitude$z
      latitude <- longitude$y
      longitude <- longitude$x
    }
  }
  if (!is.matrix(z)) 
    stop("z must be a matrix")
  breaksGiven <- !missing(breaks)
  if (!missing(colormap)) {
    breaks <- colormap$breaks
    breaksGiven <- TRUE
    col <- colormap$col
    missingColor <- colormap$missingColor
    zclip <- colormap$zclip
    colGiven <- TRUE
  }
  if (!breaksGiven) {
    small <- .Machine$double.eps
    zrange <- range(z, na.rm = TRUE)
    if (missing(zlim)) {
      if (missing(col)) {
        breaks <- pretty(zrange, n = 10)
      }
      else {
        if (is.vector(col)) {
          breaks <- pretty(zrange, n = 1 + length(col))
        }
        else if (is.function(col)) {
          breaks <- pretty(zrange, n = 10)
        }
        else {
          stop("'col' must be a vector or a function")
        }
      }
      breaksOrig <- breaks
    }
    else {
      if (!colGiven) {
        oceDebug(debug, "zlim provided, but not breaks or col\n")
        breaks <- c(zlim[1], pretty(zlim, n = 128), zlim[2])
      }
      else {
        oceDebug(debug, "zlim and col provided, but not breaks\n")
        breaks <- seq(zlim[1], zlim[2], length.out = if (is.function(col)) 
          128
          else 1 + length(col))
      }
      breaksOrig <- breaks
      breaks[1] <- min(zrange[1], breaks[1])
      breaks[length(breaks)] <- max(breaks[length(breaks)], 
                                    zrange[2])
    }
  }
  else {
    breaksOrig <- breaks
    if (1 == length(breaks)) {
      oceDebug(debug, "only 1 break given, so taking that as number of breaks\n")
      breaks <- pretty(z, n = breaks)
    }
  }
  if (missing(col)) {
    col <- oce.colorsPalette(n = length(breaks) - 1)
    oceDebug(debug, "using default col\n")
  }
  if (is.function(col)) {
    col <- col(n = length(breaks) - 1)
    oceDebug(debug, "col is a function\n")
  }
  oceDebug(debug, "zclip:", zclip, "\n")
  oceDebug(debug, vectorShow(breaks))
  oceDebug(debug, vectorShow(col))
  ni <- dim(z)[1]
  nj <- dim(z)[2]
  zmin <- min(z, na.rm = TRUE)
  zmax <- max(z, na.rm = TRUE)
  zrange <- zmax - zmin
  small <- .Machine$double.eps
  if (zclip) {
    oceDebug(debug, "using missingColor for out-of-range values\n")
    if (zlimGiven) {
      z[z < zlim[1]] <- NA
      z[z > zlim[2]] <- NA
    }
  }
  else {
    if (zlimGiven) {
      oceDebug(debug, "using 'zlim' to pin image colours\n")
      pinMIN <- min(zlim, na.rm = TRUE)
      pinMAX <- max(zlim, na.rm = TRUE)
      pinlow <- z <= pinMIN
      z[pinlow] <- pinMIN * (1 + sign(pinMIN) * small)
      pinhigh <- z >= pinMAX
      z[pinhigh] <- pinMAX * (1 - sign(pinMAX) * small)
      oceDebug(debug, "pinned ", sum(pinlow), " low values and ", 
               sum(pinhigh), " high z values, out of a total of ", 
               length(z), " values\n")
    }
    else if (breaksGiven) {
      oceDebug(debug, "using 'breaks' to pin image colours\n")
      pinMIN <- min(breaks, na.rm = TRUE)
      pinMAX <- max(breaks, na.rm = TRUE)
      pinlow <- z <= pinMIN
      z[pinlow] <- pinMIN * (1 + sign(pinMIN) * small)
      pinhigh <- z >= pinMAX
      z[pinhigh] <- pinMAX * (1 - sign(pinMAX) * small)
      oceDebug(debug, "pinned ", sum(pinlow), " low values and ", 
               sum(pinhigh), " high z values, out of a total of ", 
               length(z), " values\n")
    }
    else {
      oceDebug(debug, "not clipping AND NEITHER zlim nor breaks suppled\n")
    }
  }
  poly <- .Call("map_assemble_polygons", longitude, latitude, 
                z, NAOK = TRUE, PACKAGE = "oce")
  xy <- lonlat2map(poly$longitude, poly$latitude)
  xy$x[!is.finite(xy$x)] <- NA
  xy$y[!is.finite(xy$y)] <- NA
  usr12 <- par("usr")[1:2]
  xrange <- range(xy$x, na.rm = TRUE)
  if (xrange[1] > usr12[2]) 
    xy$x <- xy$x - 360
  Z <- as.vector(z)
  r <- .Call("map_check_polygons", xy$x, xy$y, poly$z, diff(par("usr"))[1:2]/5, 
             par("usr"), NAOK = TRUE, PACKAGE = "oce")
  breaksMin <- min(breaks, na.rm = TRUE)
  breaksMax <- max(breaks, na.rm = TRUE)
  if (filledContour) {
    oceDebug(debug, "using filled contours\n")
    zz <- Z
    g <- expand.grid(longitude, latitude)
    longitudeGrid <- g[, 1]
    latitudeGrid <- g[, 2]
    usr <- par("usr")
    N <- sum(usr[1] <= xy$x & xy$x <= usr[2] & usr[3] <= 
               xy$y & xy$y <= usr[4], na.rm = TRUE)
    NN <- sqrt(N/10)
    xg <- seq(usr[1], usr[2], length.out = NN)
    yg <- seq(usr[3], usr[4], length.out = NN)
    xy <- lonlat2map(longitudeGrid, latitudeGrid)
    good <- is.finite(zz) & is.finite(xy$x) & is.finite(xy$y)
    if (!zclip) {
      zz[zz < breaksMin] <- breaksMin
      zz[zz > breaksMax] <- breaksMax
    }
    xx <- xy$x[good]
    yy <- xy$y[good]
    zz <- zz[good]
    xtrim <- par("usr")[1:2]
    ytrim <- par("usr")[3:4]
    inFrame <- xtrim[1] <= xx & xx <= xtrim[2] & ytrim[1] <= 
      yy & yy <= ytrim[2]
    oceDebug(debug, "before trimming, length(xx): ", length(xx), 
             "\n")
    xx <- xx[inFrame]
    yy <- yy[inFrame]
    zz <- zz[inFrame]
    oceDebug(debug, "after trimming, length(xx): ", length(xx), 
             "\n")
    if (gridder %in% c("interp", "akima")) {
      oceDebug(debug, "using interp::interp()\n")
      if (requireNamespace("interp", quietly = TRUE)) {
        i <- try(interp::interp(x = xx, y = yy, z = zz, 
                                xo = xg, yo = yg))
        if (inherits(i, "try-error")) {
          stop("gridder=\"", gridder, "\" failed. Try gridder=\"binMean2D\" instead")
        }
      }
      else {
        stop("must install.packages(\"interp\") for gridder=\"", 
             gridder, "\"")
      }
    }
    else if (gridder == "binMean2D") {
      oceDebug(debug, "using binMean2D()\n")
      binned <- binMean2D(xx, yy, zz, xg, yg, fill = TRUE)
      i <- list(x = binned$xmids, y = binned$ymids, z = binned$result)
    }
    else {
      stop("gridder=\"", gridder, "\" not allowed. Try \"binMean2D\" or \"interp\"")
    }
    if (any(is.finite(i$z))) {
      small <- .Machine$double.eps
      .filled.contour(i$x, i$y, i$z, levels = breaks + 
                        small, col = col)
    }
    else {
      warning("no valid z")
    }
  }
  else {
    oceDebug(debug, "using polygons, as opposed to filled contours\n")
    colFirst <- col[1]
    colLast <- tail(col, 1)
    colorLookup <- function(ij) {
      zval <- Z[ij]
      if (!is.finite(zval)) {
        return(missingColor)
      }
      if (zval < breaksMin) {
        return(if (zclip) missingColor else colFirst)
      }
      if (zval > breaksMax) {
        return(if (zclip) missingColor else colLast)
      }
      w <- which(zval <= breaks)[1]
      if (!is.na(w) && w > 1) {
        return(col[-1 + w])
      }
      else {
        return(missingColor)
      }
    }
    method <- options()$mapPolygonMethod
    if (0 == length(method)) {
      method <- 3
    }
    oceDebug(debug, "method=", method, " (set by options()$mapPolygonMethod or default of 3)\n")
    if (method == 1) {
      colPolygon <- unlist(lapply(1:(ni * nj), colorLookup))
    }
    else if (method == 2) {
      colPolygon <- character(ni * nj)
      for (ij in 1:(ni * nj)) {
        zval <- Z[ij]
        if (!is.finite(zval)) {
          colPolygon[ij] <- missingColor
        }
        else if (zval < breaksMin) {
          colPolygon[ij] <- if (zclip) 
            missingColor
          else colFirst
        }
        else if (zval > breaksMax) {
          colPolygon[ij] <- if (zclip) 
            missingColor
          else colLast
        }
        else {
          w <- which(zval <= breaks)[1]
          colPolygon[ij] <- if (!is.na(w) && w > 1) 
            col[-1 + w]
          else missingColor
        }
      }
    }
    else if (method == 3) {
      colPolygon <- rep(missingColor, ni * nj)
      ii <- findInterval(Z, breaks, left.open = TRUE, all.inside = TRUE)
      colPolygon <- col[ii]
      colPolygon[!is.finite(Z)] <- missingColor
      if (zclip) {
        colPolygon[Z < min(breaks)] <- missingColor
        colPolygon[Z > max(breaks)] <- missingColor
      }
    }
    else {
      stop("unknown options(mapPolygonMethod)")
    }
    polygon(xy$x[r$okPoint & !r$clippedPoint], xy$y[r$okPoint & 
                                                      !r$clippedPoint], col = colPolygon[r$okPolygon & 
                                                                                           !r$clippedPolygon], border = NA, lwd = lwd, lty = lty, density=density, angle=45, fillOddEven = FALSE)
    polygon(xy$x[r$okPoint & !r$clippedPoint], xy$y[r$okPoint & 
                                                      !r$clippedPoint], col = colPolygon[r$okPolygon & 
                                                                                           !r$clippedPolygon], border = NA, lwd = lwd, lty = lty, density=density, angle=135, fillOddEven = FALSE)
  }
  oceDebug(debug, "} # mapImage()\n", unindent = 1)
  invisible(NULL)
}

environment(mapShade) <- asNamespace('oce')
assignInNamespace("mapImage", mapShade, ns = "oce")

