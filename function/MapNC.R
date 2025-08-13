MapNC <- function(
    fieldflag_lon,
    fieldflag_lat,
    fieldflag_value,
    
    shade_lon,
    shade_lat,
    shade_value,
    shade_thresh,
    shade_density=50,
    shade_lwd=0.5,
    contour_value,
    contour_lon,
    contour_lat,
    contour_breaks,
    contour_col="red",
    contour_lwd=2,
    field_value,
    field_lon,
    field_lat,
    field_breaks,
    field_pallete_name="default",
    field_pallete_starting_alpha=1,
    field_show_legend=FALSE,
    field_units="",
    legend_units_potition=1,
    field_legend_cex=3.5,
    field_legend_mai_right=1.5,
    field_legend_mai_top=0.5,
    field_show_box=FALSE,
    
    field_mean=FALSE,
    field_mean_lon,
    field_mean_lat,
    field_mean_value,
    field_mean_pch=15,
    field_mean_cex=6,
    
    field_value_mean_global=T,
    
    points_value,
    points_lon,
    points_lat,
    points_cex=2,
    points_highlight_cex=3,
    points_pch=19,
    points_colhex=NULL,
    points_alphahex="33",
    points_ncol=100,
    points_pallete="",
    points_text,
    points_show_text=FALSE,
    points_highlight=FALSE,
    points_show_legend=TRUE,
    coords_terra=FALSE,
    coords_aqua=FALSE,
    coords_tropomi=FALSE,
    filename_topo="",
    elevation_multiplier=100,
    show_legend=FALSE,
    drawMapBox=FALSE,
    projection="+proj=robin",
    lonmin=-180, 
    lonmax=180, 
    latmin=-90, 
    latmax=90,
    title_main="",
    title_legend="",
    gridlines=15,
    ocean_show=TRUE,
    col_ocean="white",
    col_land="white",
    zmin=-1,
    zmax=1,
    col_title="grey95",
    col_title_shadow="black",
    coastlineWorldFine_lwd=1,
    figure_box=FALSE,
    legend_decimals=T
) {
  
  ###############
  ### LIBRARY ###
  ###############
  library("oce")
  library("ocedata")
  library("sf")
  library("raster")
  library("ncdf4")
  library("colourvalues")
  library("colorspace")
  library("viridis")
  
  #################
  ### FUNCTIONS ###
  #################
  find_color <- function(value, palette) {
    if (value <= palette$zlim[1]) { return(palette$col[1]) } 
    if (value >= palette$zlim[2]) { return(palette$col[length(palette$col)]) }
    if (value > palette$zlim[1] & value < palette$zlim[2]) {ID <- findInterval(value, palette$breaks); return(palette$col[ID])}
  }
  
  #################
  ### COASTLINE ###
  #################
  data("coastlineWorldFine")
  
  ##################
  ### TOPOGRAPHY ###
  ##################
  if (filename_topo!="") {
    ### Read topograhy and create shade
    file.nc <- nc_open(filename_topo)
    mytopo  <- ncvar_get(file.nc, "elevation")*elevation_multiplier
    mylon   <- ncvar_get(file.nc, "longitude")
    mylat   <- ncvar_get(file.nc, "latitude")
    nc_close(file.nc)
    ID <- which(mytopo <= 0, arr.ind=T)
    mytopo[ID]=NA
    # Create raster object for topography
    mytopo_raster <- raster(ncol=dim(mytopo)[1], nrow=dim(mytopo)[2])
    values(mytopo_raster) <- c(mytopo[,dim(mytopo)[2]:1])
    # Create shaded relief
    slope  <- terrain(mytopo_raster, opt='slope')
    aspect <- terrain(mytopo_raster, opt='aspect')
    shade  <- hillShade(slope, aspect, 80, 30)
    # Create a matrix for the shaded relief
    shade <- matrix(data=values(shade), nrow=dim(mytopo)[1],  ncol=dim(mytopo)[2])
    shade <- shade[,dim(shade)[2]:1]
    # Create topography colors
    colormap_shade <- colormap( 
      col=grey(0:100/100), # paste0(grey(0:100/100) ),# , c(rep("99",70), sprintf(fmt="%02d",seq(30,0,-1))) ), 
      breaks=seq(min(shade,na.rm=T),max(shade, na.rm=T), length.out=100),
      missingColor="#FFFFFF00"
    )
  }
  
  ####################
  ### FIELD COLORS ###
  ####################
  ### Breaks for color pallete
  if (missing(field_breaks)) { field_breaks <- c(seq(0,0.700,0.005),1.130,1.560,1.990,2.420,2.850,3.280,3.710,4.140,4.570,5.000) }
  ### Colors for color pallete
  if (field_pallete_name=="default") {
    field_pallete <- colormap( 
      col=colorRampPalette(c("#fffcc6","#ffe579","#ffca4d","#ff9f2e","#ff731e","#fc3712","#e21210","#b7000e","#76000d"))(length(field_breaks)-1),
      breaks=field_breaks,
      missingColor="#FFFFFF00")}
  
  if (field_pallete_name=="distinct") {
    field_pallete <- colormap( 
      #col=colorRampPalette(c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))(length(field_breaks)-1),
      col=colorRampPalette(c("#800000","#9A6324","#808000","#469990","#000075","#e6194B","#f58231","#ffe119","#bfef45","#3cb44b","#42d4f4","#42d4f4","#911eb4","#f032e6"))(length(field_breaks)-1),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    field_pallete$col <- paste0(field_pallete$col, "0D")
  }
  if (field_pallete_name=="diff2_NoAlpha") {
    field_pallete <- colormap( 
      col=colorRampPalette(c("navy","deepskyblue","lightskyblue1","white","khaki1","orange","darkred"))(length(field_breaks)-1),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    field_pallete$col <- paste0(field_pallete$col, "FF")
  }
  if (field_pallete_name=="diff2") {
    field_pallete <- colormap( 
      col=colorRampPalette(c("navy","deepskyblue","lightskyblue1","white","khaki1","orange","darkred")),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)/2) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col[1:(length(field_pallete$col)/2)] <- paste0(field_pallete$col[1:(length(field_pallete$col)/2)], rev(hex_alpha_values))
    field_pallete$col[(1+length(field_pallete$col)/2):length(field_pallete$col)] <- paste0(field_pallete$col[(1+length(field_pallete$col)/2):length(field_pallete$col)], hex_alpha_values)
    #plot(1:24, col=field_pallete$col, pch=19, cex=3)
  }
  if (field_pallete_name=="FGE_diff") {
    field_pallete <- colormap( 
      col=colorRampPalette(rev(c("darkred","orange","white","green3", "green4"))),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)/2) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col[1:(length(field_pallete$col)/2)] <- paste0(field_pallete$col[1:(length(field_pallete$col)/2)], rev(hex_alpha_values))
    field_pallete$col[(1+length(field_pallete$col)/2):length(field_pallete$col)] <- paste0(field_pallete$col[(1+length(field_pallete$col)/2):length(field_pallete$col)], hex_alpha_values)
    #plot(1:24, col=field_pallete$col, pch=19, cex=3)
  }
  if (field_pallete_name=="R_diff") {
    field_pallete <- colormap( 
      col=colorRampPalette(c("darkred","orange","white","green3", "green4")),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)/2) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col[1:(length(field_pallete$col)/2)] <- paste0(field_pallete$col[1:(length(field_pallete$col)/2)], rev(hex_alpha_values))
    field_pallete$col[(1+length(field_pallete$col)/2):length(field_pallete$col)] <- paste0(field_pallete$col[(1+length(field_pallete$col)/2):length(field_pallete$col)], hex_alpha_values)
    #plot(1:24, col=field_pallete$col, pch=19, cex=3)
  }
  if (field_pallete_name=="MNMB") {
    field_pallete <- colormap( 
      col=colorRampPalette(c("navy","blue","lightskyblue", "white","palevioletred1","red","red4")),
      breaks=field_breaks,
      missingColor="grey65")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)/2) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col[1:(length(field_pallete$col)/2)] <- paste0(field_pallete$col[1:(length(field_pallete$col)/2)], rev(hex_alpha_values)[1:(length(field_pallete$col)/2)])
    field_pallete$col[(1+length(field_pallete$col)/2):length(field_pallete$col)] <- paste0(field_pallete$col[(1+length(field_pallete$col)/2):length(field_pallete$col)], hex_alpha_values[1:(length(field_pallete$col)/2)])
  }
  if (field_pallete_name=="FGE") {
    field_pallete <- colormap( 
      col=colorRampPalette(rev(c("darkred","red","orange","yellow","greenyellow", "green3", "green4"))),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)/2) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col[1:(length(field_pallete$col)/2)] <- paste0(field_pallete$col[1:(length(field_pallete$col)/2)], rev(hex_alpha_values))
    field_pallete$col[(1+length(field_pallete$col)/2):length(field_pallete$col)] <- paste0(field_pallete$col[(1+length(field_pallete$col)/2):length(field_pallete$col)], hex_alpha_values)
  }
  if (field_pallete_name=="R") {
    field_pallete <- colormap( 
      col=colorRampPalette(c("darkred","red","orange","yellow","greenyellow", "green3", "green4")),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)/2) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col[1:(length(field_pallete$col)/2)] <- paste0(field_pallete$col[1:(length(field_pallete$col)/2)], rev(hex_alpha_values))
    field_pallete$col[(1+length(field_pallete$col)/2):length(field_pallete$col)] <- paste0(field_pallete$col[(1+length(field_pallete$col)/2):length(field_pallete$col)], hex_alpha_values)
  }
  if (field_pallete_name=="colorblind2") {
    field_pallete <- colormap( 
      col=colorRampPalette(c("white","#E1C6DD","#85C0F9","#EBCB3B","#F5793A","#A95AA1","#382119"))(length(field_breaks)-1),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col <- paste0(field_pallete$col, hex_alpha_values)
  }
  if (field_pallete_name=="AOD") {
    field_pallete <- colormap( 
      col=colorRampPalette(c("#E1C6DD","#85C0F9","#EBCB3B","#F5793A","#A95AA1","#382119"))(length(field_breaks)-1),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col <- paste0(field_pallete$col, hex_alpha_values)
  }
  if (field_pallete_name=="AOD_white") {
    field_pallete <- colormap( 
      col=colorRampPalette(c("white","#E1C6DD","#85C0F9","#EBCB3B","#F5793A","#A95AA1","#382119"))(length(field_breaks)-1),
      breaks=field_breaks,
      missingColor="grey65")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col <- paste0(field_pallete$col, hex_alpha_values)
  }
  if (field_pallete_name=="AE") {
    field_pallete <- colormap( 
      col=colorRampPalette(c("navy","blue","yellow","red","darkred"))(length(field_breaks)-1),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col <- paste0(field_pallete$col, hex_alpha_values)
  }
  if (field_pallete_name=="SSA") {
    field_pallete <- colormap( 
      col=colorRampPalette(magma(length(field_breaks)-1)),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col <- paste0(field_pallete$col, hex_alpha_values)
  }
  if (field_pallete_name=="SSA_diff") {
    field_pallete <- colormap( 
      col=colorRampPalette(c("purple4","purple","plum2","white","khaki1","darkgoldenrod1","darkorange1"))(length(field_breaks)-1),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col <- paste0(field_pallete$col, hex_alpha_values)
  }
  if (field_pallete_name=="AAOD") {
    field_pallete <- colormap( 
      col=colorRampPalette(rev(magma(length(field_breaks)-1))),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col <- paste0(field_pallete$col, hex_alpha_values)
  }
  if (field_pallete_name=="AAOD_diff") {
    field_pallete <- colormap( 
      col=colorRampPalette(rev(c("purple4","purple","plum2","white","khaki1","darkgoldenrod1","darkorange1")))(length(field_breaks)-1),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col <- paste0(field_pallete$col, hex_alpha_values)
  }
  if (field_pallete_name=="Temperature") {
    field_pallete <- colormap( 
      col=colorRampPalette(c("royalblue","cadetblue2","yellow","lightsalmon","orangered", "darkred"))(length(field_breaks)-1),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col <- paste0(field_pallete$col, hex_alpha_values)
  }
  if (field_pallete_name=="Clouds") {
    field_pallete <- colormap( 
      col=colorRampPalette(c("#583dff","lightblue","#C9BDA1","#8B826C","dimgray","grey20"))(length(field_breaks)-1),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col <- paste0(field_pallete$col, hex_alpha_values)
  }
  if (field_pallete_name=="Clouds_diff") {
    field_pallete <- colormap( 
      col=colorRampPalette(c("#583dff","lightblue","white","#a59c9c","grey20"))(length(field_breaks)-1),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col <- paste0(field_pallete$col, hex_alpha_values)
  }
  
  if (field_pallete_name=="Pressure") {
    field_pallete <- colormap( 
      col=colorRampPalette(c("navy","#6200AC","orange","#AF005F","darkred"))(length(field_breaks)-1),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col <- paste0(field_pallete$col, hex_alpha_values)
  }
  if (field_pallete_name=="Pressure_diff") {
    field_pallete <- colormap( 
      col=colorRampPalette(c("#000080","#ca85ff","white","#ff6ca3","#b30000"))(length(field_breaks)-1),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col <- paste0(field_pallete$col, hex_alpha_values)
  }
  if (field_pallete_name=="WindComponent") {
    field_pallete <- colormap( 
      #col=colorRampPalette(c("#F3F0E1","#BCD998","#FCE478","#C87F8B","#A25566"))(length(field_breaks)-1),
      col=colorRampPalette(c("#203e30","#4c9160","white","#b95b6b","#5f3252"))(length(field_breaks)-1),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col <- paste0(field_pallete$col, hex_alpha_values)
  }
  if (field_pallete_name=="WindSpeed") {
    field_pallete <- colormap( 
      col=colorRampPalette(c("white","#ffba9a","#ff7b00","#800080","#3a004d"))(length(field_breaks)-1),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col <- paste0(field_pallete$col, hex_alpha_values)
  }
  if (field_pallete_name=="WindSpeed_diff") {
    field_pallete <- colormap( 
      col=colorRampPalette(c("#b35600","#ffe8a8","white","#bb93f7","#3a004d"))(length(field_breaks)-1),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col <- paste0(field_pallete$col, hex_alpha_values)
  }
  if (field_pallete_name=="Precipitation") {
    field_pallete <- colormap( 
      #col=colorRampPalette(c("#F9EFB7","#ABD4FF","#8E82FC","#3A25F6","#184F5A"))(length(field_breaks)-1),  # Precipitation Option 1
      col=colorRampPalette(c("#e8fce2","#83f03f","#4b7d59","#6fdafc","#014ef4"))(length(field_breaks)-1),  # Precipitation Option 2
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col <- paste0(field_pallete$col, hex_alpha_values)
  }
  if (field_pallete_name=="Precipitation_diff") {
    field_pallete <- colormap( 
      col=colorRampPalette(c("#4b3231","#bb3e18","yellow","white","#83f03f","#01665e","#014ef4"))(length(field_breaks)-1),  # Precipitation Option 2
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col <- paste0(field_pallete$col, hex_alpha_values)
  }
  if (field_pallete_name=="TROPOMI_NEW") {
    field_pallete <- colormap( 
      col=colorRampPalette(c("#EBF7FD","#C9EAFB","#ADE1F6","#9BD0EE","#89C3E6","#7EBAE2","#71AFDC","#60A3D5","#62B49B","#8ACE64","#D0DF6F","#FAE771","#FACD64",
                             "#F7B95B","#F8A750","#FB9548","#F6813F","#EA5B3E","#CA1112","#A30008","#820005","#640203","#4E0002","#360202","#240302","#100202"))(length(field_breaks)-1),
      breaks=field_breaks,
      missingColor="#FFFFFF00")
    # Create a sequence of alpha values from 1 to 0
    alpha_values <- log( seq(field_pallete_starting_alpha, 100, length.out = length(field_pallete$col)) )
    # Normalize the alpha values to be between 0 and 1
    alpha_values <- alpha_values / max(alpha_values)
    # Convert alpha values to hexadecimal
    hex_alpha_values <- sprintf("%02X", round(alpha_values * 255))
    # Apply hex alpha
    field_pallete$col <- paste0(field_pallete$col, hex_alpha_values)
  }
  
  #####################
  ### POINTS COLORS ### -> based on field colorspace 
  #####################
  if (missing(points_value)==FALSE) {
    if (field_pallete_name!="") {
      for (n in 1:length(points_value))
        points_colhex <- c( points_colhex, find_color(points_value[n], field_pallete) )
    }
  }
  
  #############################
  ### SATELLITE COORDINATES ###
  #############################
  if (coords_terra != FALSE) {
    ### TERRA
    # Read Terra coords (Produced using the PyOrbital library in python)
    satellite_coords_terra <- read.table(coords_terra, head=T, sep=",")
    # Calculate differences between consecutive elements
    differences <- diff(c(satellite_coords_terra$Latitude,satellite_coords_terra$Latitude[1]))
    # Check if differences are positive, negative, or zero
    satellite_coords_terra$Status <- ifelse(differences > 0, "Ascending", ifelse(differences < 0, "Descending", "Constant"))
    # Keep only Descending for Terra
    ID <- which(satellite_coords_terra$Status == "Descending")
    satellite_coords_terra <- satellite_coords_terra[ID,]
  }
  if (coords_aqua != FALSE) {
    ### AQUA
    # Read Terra coords (Produced using the PyOrbital library in python)
    satellite_coords_aqua <- read.table(coords_aqua, head=T, sep=",")
    # Calculate differences between consecutive elements
    differences <- diff(c(satellite_coords_aqua$Latitude,satellite_coords_aqua$Latitude[1]))
    # Check if differences are positive, negative, or zero
    satellite_coords_aqua$Status <- ifelse(differences > 0, "Ascending", ifelse(differences < 0, "Descending", "Constant"))
    # Keep only Descending for Terra
    ID <- which(satellite_coords_aqua$Status == "Ascending")
    satellite_coords_aqua <- satellite_coords_aqua[ID,]
  }
  if (coords_tropomi != FALSE) {
    ### TROPOMI
    # Read TROPOMI coords (Produced using the PyOrbital library in python)
    satellite_coords_tropomi <- read.table(coords_tropomi, head=T, sep=",")
    # Calculate differences between consecutive elements
    differences <- diff(c(satellite_coords_tropomi$Latitude,satellite_coords_tropomi$Latitude[1]))
    # Check if differences are positive, negative, or zero
    satellite_coords_tropomi$Status <- ifelse(differences > 0, "Ascending", ifelse(differences < 0, "Descending", "Constant"))
    # Keep only Descending for Terra
    ID <- which(satellite_coords_tropomi$Status == "Ascending")
    satellite_coords_tropomi <- satellite_coords_tropomi[ID,]
  }
  
  ################
  ### PLOT MAP ###
  ################
  if (projection=="+proj=robin") { par(mar=c(0,0,0,0)) }
  if (projection!="+proj=robin") { par(mar=c(3,3,0,0)) }
  ### Add ocean color
  if (ocean_show == TRUE) {
    plot(0, 0, type="n", ann=FALSE, axes=FALSE)
    user <- par("usr")
    rect(user[1], user[3], user[2], user[4], col=col_ocean, border=NA)
    ### Global field value mean
    if (field_value_mean_global==T) {
      points_x <- -0.85
      points_y <- -0.88
      #points(points_x,points_y,col="black", pch=15, cex=30)
      #points(points_x,points_y,col="white", pch=15, cex=29)
      text(points_x,points_y,paste0("MN\n",format(mean(field_value,na.rm=T), scientific=TRUE, digits=2)),col="black", pch=15, cex=4.0, family="Century Gothic")
      points_x <- +0.85
      points_y <- -0.88
      #points(points_x,points_y,col="black", pch=15, cex=30)
      #points(points_x,points_y,col="white", pch=15, cex=29)
      text(points_x,points_y,paste0("SD\n",format(sd(field_value,na.rm=T), scientific=TRUE, digits=2)),col="black", pch=15, cex=4.0, family="Century Gothic")
    }
    par(new=TRUE)
    par(new=TRUE)
  }
  
  ### Map
  mapPlot(coastlineWorldFine,
          grid=FALSE,
          projection=projection, 
          axes=F,
          col=col_land, 
          border="white",
          drawBox=drawMapBox,
          longitudelim=c(lonmin, lonmax), 
          latitudelim =c(latmin, latmax))
  
  ### Orography
  if (filename_topo!="") { mapImage(mylon, mylat, shade, colormap=colormap_shade, border=NA) }
  
  ### Field
  if (missing(field_value)==FALSE) { mapImage(field_lon, field_lat, field_value, colormap=field_pallete ) }
  
  ### Field
  if (missing(fieldflag_value)==FALSE) { 
    flag_pallete <- colormap( 
      col=colorRampPalette(c("#80FA7f","#FFFFFF00"))(length(field_breaks)-1),
      breaks=c(0,1),
      missingColor="#FFFFFF00")
    flag_pallete$col <- paste0(flag_pallete$col, c("FF","00"))
    
    mapImage(fieldflag_lon, fieldflag_lat, fieldflag_value, colormap=flag_pallete ) 
  }
  
  
  
  
  
  ### Density
  #if (missing(shade_value)==FALSE) { mapShade(shade_lon, shade_lat, shade_value, col=c("#FFFFFF00","grey50"), breaks=c(-Inf,shade_thresh,Inf), density=shade_density, lwd=shade_lwd) }
  
  ### Satellite orbit Terra
  if (coords_terra != FALSE) {
    mapLines(longitude=satellite_coords_terra$Longitude, latitude=satellite_coords_terra$Latitude, lty=1, lwd=0.5, col="red")
    N <- length(satellite_coords_terra$Longitude)
    mapPoints(longitude=satellite_coords_terra$Longitude, 
              latitude=satellite_coords_terra$Latitude, pch=19, cex=0.3, col="red")
    ID <- which(substr(satellite_coords_terra$Datetime,16,19)=="0:00" | substr(satellite_coords_terra$Datetime,16,19)=="5:00")
    mapPoints(longitude=satellite_coords_terra$Longitude[ID], 
              latitude=satellite_coords_terra$Latitude[ID], pch=19, col="red")
    mapText(longitude=satellite_coords_terra$Longitude[ID], 
            latitude=satellite_coords_terra$Latitude[ID], 
            labels=paste0(substr(satellite_coords_terra$Datetime[ID],12,16),"            "), cex=2, col="red")
  }
  ### Satellite orbit Aqua
  if (coords_aqua != FALSE) {
    mapLines(longitude=satellite_coords_aqua$Longitude, latitude=satellite_coords_aqua$Latitude, lty=1, lwd=0.5, col="blue")
    N <- length(satellite_coords_aqua$Longitude)
    mapPoints(longitude=satellite_coords_aqua$Longitude, 
              latitude=satellite_coords_aqua$Latitude, pch=19, cex=0.3, col="blue")
    ID <- which(substr(satellite_coords_aqua$Datetime,16,19)=="0:00" | substr(satellite_coords_aqua$Datetime,16,19)=="5:00")
    mapPoints(longitude=satellite_coords_aqua$Longitude[ID], 
              latitude=satellite_coords_aqua$Latitude[ID], pch=19, col="blue")
    mapText(longitude=satellite_coords_aqua$Longitude[ID], 
            latitude=satellite_coords_aqua$Latitude[ID], 
            labels=paste0(substr(satellite_coords_aqua$Datetime[ID],12,16),"            "), cex=2, col="blue")
  }
  ### Satellite orbit TROPOMI
  if (coords_tropomi != FALSE) {
    mapLines(longitude=satellite_coords_tropomi$Longitude, latitude=satellite_coords_tropomi$Latitude, lty=1, lwd=0.5, col="purple")
    N <- length(satellite_coords_tropomi$Longitude)
    mapPoints(longitude=satellite_coords_tropomi$Longitude, 
              latitude=satellite_coords_tropomi$Latitude, pch=19, cex=0.3, col="purple")
    ID <- which(substr(satellite_coords_tropomi$Datetime,16,19)=="0:00" | substr(satellite_coords_tropomi$Datetime,16,19)=="5:00")
    mapPoints(longitude=satellite_coords_tropomi$Longitude[ID], 
              latitude=satellite_coords_tropomi$Latitude[ID], pch=19, col="purple")
    mapText(longitude=satellite_coords_tropomi$Longitude[ID], 
            latitude=satellite_coords_tropomi$Latitude[ID], 
            labels=paste0(substr(satellite_coords_tropomi$Datetime[ID],12,16),"            "), cex=2, col="purple")
  }
  
  ### Coastline
  mapLines(coastlineWorldFine, lwd=coastlineWorldFine_lwd, col="grey20")
  mapGrid(col="grey40", lty=1, lwd=0.5, dlongitude=gridlines, dlatitude=gridlines)
  
  ### Contour
  if (missing(contour_value)==FALSE) { 
    mapContour(contour_lon, contour_lat, contour_value, levels=contour_breaks, col=contour_col, lwd=contour_lwd, debug=T) 
  }
  
  ### Points
  if (missing(points_value)==FALSE) {
    mapPoints(longitude=points_lon, latitude=points_lat, pch=points_pch, col=paste0("#000000"), cex=points_cex+points_cex*0.2) # Black outline
    mapPoints(longitude=points_lon, latitude=points_lat, pch=points_pch, col=paste0("#FFFFFF"), cex=points_cex) # White inside to display alpha properly
    mapPoints(longitude=points_lon, latitude=points_lat, pch=points_pch, col=paste0(points_colhex), cex=points_cex)
  }
  
  ### Points Text
  # if (points_show_text==TRUE) { mapText(longitude=points_lon, latitude=points_lat, labels=points_value, col=paste0("#000000"), cex=points_cex*0.2) }
  if (points_show_text==TRUE) { mapText(longitude=points_lon, latitude=points_lat, labels=points_text, col=paste0("#000000"), cex=points_cex*0.25) }
  
  ### Points Highlight: I have set it up as the last entry in the provided (input) statistics/points
  if (points_highlight==TRUE) {
    N <- length(points_lon)
    mapPoints(longitude=points_lon[N], latitude=points_lat[N], pch=points_pch, col=paste0("#000000"), cex=points_highlight_cex*5+points_highlight_cex*3*0.15) # Black outline
    mapPoints(longitude=points_lon[N], latitude=points_lat[N], pch=points_pch, col=paste0("#FFFFFF"), cex=points_highlight_cex*5) # White inside to display alpha properly
    mapPoints(longitude=points_lon[N], latitude=points_lat[N], pch=points_pch, col=paste0(points_colhex[N]), cex=points_highlight_cex*5)
    mapText(  longitude=points_lon[N], latitude=points_lat[N], labels=points_value[N], col=paste0("#000000"), cex=points_highlight_cex*5*0.2, family="Century Gothic")
  }
  
  ### Field Polygon
  if (field_show_box) {
    field_lon_box <- field_lon[which(field_lon > lonmin & field_lon < lonmax)]
    field_lat_box <- field_lat[which(field_lat > latmin & field_lat < latmax)]
    xn <- length(field_lon_box)
    yn <- length(field_lat_box)
    px <- abs(mean(diff(field_lon_box))/2)
    py <- abs(mean(diff(field_lat_box))/2)
    #print(xn)
    #print(yn)
    #print(px)
    #print(py)
    mapPolygon( lwd=6, border="black",
                longitude=c( rep(field_lon_box[1],yn)-px , field_lon_box-px             , rep(field_lon_box[xn],yn)+px , rev(field_lon_box+px) ),
                latitude =c( field_lat_box+py            , rep(field_lat_box[yn],xn)-py , rev(field_lat_box)-py        , rep(field_lat_box[1],xn)+py )
    )
    mapPolygon( lwd=3, border="green",
                longitude=c( rep(field_lon_box[1],yn)-px , field_lon_box-px             , rep(field_lon_box[xn],yn)+px , rev(field_lon_box+px) ),
                latitude =c( field_lat_box+py            , rep(field_lat_box[yn],xn)-py , rev(field_lat_box)-py        , rep(field_lat_box[1],xn)+py )
    )
  }
  
  ### Field mean
  if (field_mean==TRUE) {
    find_color <- function(value, palette) {
      if (value <= palette$zlim[1]) { return(palette$col[1]) } 
      if (value >= palette$zlim[2]) { return(palette$col[length(palette$col)]) }
      if (value > palette$zlim[1] & value < palette$zlim[2]) {ID <- findInterval(value, palette$breaks); return(palette$col[ID])}
    }
    #print(field_mean_value)
    #print(field_pallete)
    field_mean_col <- find_color(field_mean_value, field_pallete)
    # Get luminance of the background
    lumen <- as(hex2RGB(field_mean_col), "polarLUV")@coords[1] # Library colorspace
    if (lumen >= 50) {text_color <- "#000000"}
    if (lumen <  50) {text_color <- "#FFFFFF"}
    
    mapPoints(longitude=field_mean_lon, latitude=field_mean_lat, pch=field_mean_pch, col=paste0("#000000"), cex=field_mean_cex*3+field_mean_cex*3*0.15) # Black outline
    mapPoints(longitude=field_mean_lon, latitude=field_mean_lat, pch=field_mean_pch, col=paste0("#FFFFFF"), cex=field_mean_cex*3) # White inside to display alpha properly
    mapPoints(longitude=field_mean_lon, latitude=field_mean_lat, pch=field_mean_pch, col=field_mean_col, cex=field_mean_cex*3)
    mapText(  longitude=field_mean_lon, latitude=field_mean_lat, labels=field_mean_value, col=text_color, cex=field_mean_cex*3*0.2)
  }
  
  ### Lon-Lat lines
  par(cex.axis=2)  # Change 2 to the desired size
  mapAxis(side=1, longitude=T, axisStyle=5, cex.axis=4)
  mapAxis(side=2, latitude=T,  axisStyle=5, cex.axis=4)
  
  ### Add title
  #title(main=title_main, cex.main=3, line=-2, adj=0)
  title(main=title_main, cex.main=5, line=-3.0, adj=0.003, col.main=col_title_shadow)
  title(main=title_main, cex.main=5, line=-3.0, adj=0.000, col.main=col_title_shadow)
  title(main=title_main, cex.main=5, line=-3.2, adj=0.003, col.main=col_title_shadow)
  title(main=title_main, cex.main=5, line=-3.2, adj=0.000, col.main=col_title_shadow)
  title(main=title_main, cex.main=5, line=-3.1, adj=0.001, col.main=col_title)
  
  ### Legend Points
  if (missing(points_value)==FALSE & points_show_legend==TRUE) {
    if (points_pallete!="") {
      legend("bottomright", title=title_legend, title.cex=2, cex=2, pch=19, pt.cex=3, legend=rev(mycolorvalues$summary_values), col=rev(mycolorvalues$summary_colours), box.lwd=2)
    }
  }
  
  ### Box around the map
  if (figure_box==TRUE) { box(lwd=1, col="black") }
  
  ### Legend Field
  if (field_show_legend==TRUE) {
    
    par(bg="#FFFFFFFF")
    par(mai=c(0.5,0.5,0.5,0.5))
    try({drawPalette(at=1:length(field_breaks), labels=field_breaks, fullpage=TRUE, col=field_pallete$col, las=1, mai=c(0,0.2,field_legend_mai_top,field_legend_mai_right), drawTriangles=T, cex=0)}, silent =T)
    mylegend_at     <- 1:length(field_breaks)
    mylegend_labels <- field_breaks
    if (legend_decimals==T) { axis(4, at=mylegend_at, labels=format(mylegend_labels, scientific=TRUE, digits=2), cex.axis=field_legend_cex, tick=T, las=1, family="Century Gothic") }
    if (legend_decimals==F) { axis(4, at=mylegend_at, labels=mylegend_labels, cex.axis=field_legend_cex, tick=T, las=1, family="Century Gothic") }
    axis(3,legend_units_potition, field_units, las=1, line=2.5, tick=F, cex.axis=field_legend_cex, family="Century Gothic")
    
    ### Box around the field legend
    box(lwd=2, col="black")
  }
  
  
}


