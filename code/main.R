# -####################################################################################################################-
#' TAP WATER
#' 
#' Runs with R v3.5
# -####################################################################################################################-

library(data.table)
library(ggplot2)
library(scales)
library(lubridate)

# Under Linux, `libudunits2-dev` should be installed first
library(spdep)
library(broom)
library(akima)

# Under Linux, `libgdal-dev ` should be installed first
library(rgdal)      # Lire et reprojeter les cartes
library(plotrix)    # Créer des échelles de couleurs
library(classInt)   # Affecter ces couleurs aux données
library(RColorBrewer)
library(spatstat)
library(png)
library(ggrepel)
library(gridExtra)

ROOT_DATA_DIR <- paste0(getwd(), "/data/")
ROOT_OUTPUT_DIR <- paste0(getwd(), "/output/")

# -####################################################################################################################-
#' Create all the plots
#'
#' @export
#'
#' @examples 
#' tapwater.runAll()
tapwater.runAll <- function(gridSize = 6000, forceUpdate = TRUE) {
  for (species in c("Delta_18_16", "Delta_D_H", "DX")) {
    tapwater.run(species, "measure", gridSize = gridSize, forceUpdate = forceUpdate)
    tapwater.run(species, "model", "precip", gridSize = gridSize, forceUpdate = forceUpdate)
    tapwater.run(species, "model", "ground", gridSize = gridSize, forceUpdate = forceUpdate)
  }
  
  tapwater.runByRange("measure", gridSize = gridSize)
  tapwater.runByRange("model", "precip", gridSize = gridSize)
  tapwater.runByRange("model", "ground", gridSize = gridSize)
}

# -####################################################################################################################-
#' Process data and plot a map related to the input parameters.
#'
#' @param species String. Name of the species to be plotted. Can be either `Delta_18_16`, `Delta_D_H` or `DX`.
#' @param measOrModel String. Data coming either from `measure` or `model`. Default is `measure.`
#' @param modelType String. If `measOrModel == "model"`, specify the model type: either `precip` for modelized 
#'                  precipitations, or `ground` for ground water modelisation.
#' @param gridSize Numeric. Meters between 2 grid cells' center (Lambert93 projection). Default is `6000` (6km). 
#' @param selectedRange List of 2 numerics. Optional. If set, only display the isotopic values which are within the
#'                      given range. The rest of the map will stay in grey.
#' @param forceUpdate Boolean. If `FALSE`, try to get already-processed data from a previous run. Useful to speed the 
#'                    process if the underlying data did not change but the plot parameters were modified.
#' 
#' @export
#'
#' @examples 
#' tapwater.run("Delta_18_16")
#' tapwater.run("Delta_D_H")
#' tapwater.run("DX")
tapwater.run <- function(species, measOrModel = "measure", modelType = "", gridSize = 6000, selectedRange, forceUpdate = FALSE) {
  isotopeDt <- tapwater.getIsotopeDt(species, measOrModel, modelType)
  isotopeDt <- tapwater.toLambert93(isotopeDt)
  interpDt <- tapwater.getInterpData(isotopeDt, gridSize)
  
  # Get Europe's coutries boundaries
  europeDt <- tapwater.getEuropeDt()
  
  # Get France's boundaries (without overseas territories)
  franceDt <- europeDt[id == 22]
  tempFranceDt <- copy(franceDt)
  franceDt <- tempFranceDt[long > 0 & long < 2000000] # Remove overseas France
  
  
  processedFileName <- paste0(ROOT_DATA_DIR, "processed/", species, "_", measOrModel, modelType,  "_", gridSize, ".rds")
  if (!file.exists(processedFileName) | forceUpdate) {
    inLandDt <- tapwater.getInLandData(franceDt, interpDt)
    saveRDS(inLandDt, processedFileName)
  } else {
    inLandDt <- readRDS(processedFileName)
  }
  
  # Get "measure" interpolated inland data as a reference for color scale
  measureDt <- readRDS(paste0(ROOT_DATA_DIR, "processed/", species, "_measure_", gridSize, ".rds"))
  measureDt <- inLandDt

  # Create the plot file.
  plot <- tapwater.getPlot(isotopeDt, inLandDt, measureDt, europeDt, species, measOrModel, modelType, selectedRange)
}

# -####################################################################################################################-
#' Create plots for isotope's range.
#'
#' @examples tapwater.runByRange()
tapwater.runByRange <- function(measOrModel = "measure", modelType = "", gridSize = 6000) {
  boundaries <- seq(-3, -15)
  for (i in seq(2, length(boundaries))) {
    print(paste0("Processing Delta_18_16: [", boundaries[i-1], ';', boundaries[i], "]"))
    tapwater.run("Delta_18_16", measOrModel, modelType, gridSize, selectedRange = c(boundaries[i-1], boundaries[i]))
  }
  
  boundaries <- seq(-105, -21, by=7)
  for (i in seq(2, length(boundaries))) {
    print(paste0("Processing Delta_D_H: [", boundaries[i-1], ';', boundaries[i], "]"))
    tapwater.run("Delta_D_H", measOrModel, modelType, gridSize, selectedRange = c(boundaries[i-1], boundaries[i]))
  }
}

# -####################################################################################################################-
#' Build the map plot and save it as PNG image file.
#'
#' @param isotopeDt 
#' @param interpDt 
#' @param measureDt 
#' @param europeDt 
#' @param species 
#' @param measOrModel 
#' @param modelType 
#' @param selectedRange 
#'
#' @return Creates a PNG file.
tapwater.getPlot <- function(isotopeDt, interpDt, measureDt, europeDt, species, measOrModel, modelType, selectedRange) {
  
  # Define the thresholds for iso-lines, and some label display
  if (species == "Delta_18_16") {
    majorContour <- 2
    minorContour <- 0.5
    speciesString <- "d18O"
    speciesExpression <- expression(paste(delta^18, "O"))
  } else if (species == "Delta_D_H") {
    majorContour <- 20
    minorContour <- 5
    speciesString <- "d2H"
    speciesExpression <- expression(paste(delta^2, "H"))
  } else if (species == "DX") {
    majorContour <- 3
    minorContour <- 1
    speciesString <- "dx"
    speciesExpression <- "dx"
  }
  
  # Filter data for specific range plotting
  if (!missing(selectedRange)) {
    rasterInterpDt <- interpDt[value >= min(selectedRange) & value <= max(selectedRange)]
    selectedPrct <- nrow(rasterInterpDt) / nrow(interpDt) * 100
  } else {
    rasterInterpDt <- interpDt
  }
  
  # Build plot's title
  plotTitle <- paste0(speciesString, " - ",
                      measOrModel, ifelse(modelType != "",
                                          paste0("_", modelType), ""),
                      ifelse(!missing(selectedRange), 
                             paste0(" - [", min(selectedRange), ";", max(selectedRange),
                                    "] highlighted (", round(selectedPrct,1), "%)"),
                             ""))
  
  nbColors = 6
  plot <- ggplot() +
    ggtitle(plotTitle) +
    geom_polygon(data = europeDt,
                 mapping = aes(x = long, y = lat, group = group),
                 fill = "grey80",
                 size = 0.2,
                 color = "black") +
    geom_raster(data = rasterInterpDt,
                mapping = aes(x = lon, y = lat, fill = value),
                interpolate = TRUE,
                na.rm = TRUE) +
    geom_contour(data = interpDt,
                 mapping = aes(x = lon, y = lat, z = value, color = value),
                 binwidth = 0.5,
                 breaks = seq(floor(min(interpDt$value)), ceiling(max(interpDt$value)), minorContour),
                 size = 0.5,
                 color = "black",
                 na.rm = TRUE) +
    geom_contour(data = interpDt,
                 mapping = aes(x = lon, y = lat, z = value, color = value),
                 binwidth = 0.5,
                 breaks = seq(floor(min(interpDt$value)), ceiling(max(interpDt$value)), majorContour),
                 size = 1,
                 color = "black",
                 na.rm = TRUE) +
    geom_point(data = isotopeDt,
               mapping = aes(x = lon, y = lat, fill = value),
               shape = 21,
               color = "black",
               na.rm = TRUE) +
    geom_path(data = europeDt,
                 mapping = aes(x = long, y = lat, group = group),
                 size = 0.2,
                 color = "black") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    guides(fill = guide_colourbar(barwidth = 1, barheight = 40)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 14),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    coord_fixed(xlim = c(50000, 1300000),
                ylim = c(6000000, 7200000))
    
  if (species == "DX") {
    plot <- plot +
      scale_fill_gradientn(colors = rev(rainbow(nbColors)),
                           name = speciesExpression
      )
  } else {
    plot <- plot +
      scale_fill_gradientn(colors = rev(rainbow(nbColors)),
                           name = speciesExpression,
                           limits = tapwater.roundScale(range(measureDt$value), digits = 1),
                           breaks = tapwater.roundScale(as.numeric(quantile(measureDt$value, seq(0, 1, length.out = nbColors))), digits = 1),
                           values = scales::rescale(as.numeric(quantile(measureDt$value, seq(0, 1, length.out = nbColors))))
      )
  }
  
  # Build file name
  filename <- paste0(strftime(now(tz = "GMT"), format = "%Y%m%d", tz = "GMT"), "_",
                     speciesString, "_",
                     measOrModel, modelType,
                     ifelse(!missing(selectedRange), paste0("_[", min(selectedRange), "_", max(selectedRange), "]"), ""))
                     

  png(paste0(ROOT_OUTPUT_DIR, filename, ".png"),
      width = 1500, height = 1200, res = 120)
  print(plot)
  dev.off()
  
  return(plot)
}



# -####################################################################################################################-
#' Get a data.table of Europe's coutries' borders projected with Lambert 93.
#' 
#' Note: Continental France and Corsica are within :
#' * 100000 < X < 1300000
#' * 6000000 < Y < 7200000
#' Plot should be in 1:1 ratio to respect Lambert 93 projection.
#' 
#' Continental France and Corsica are ID = 22
#' 
#' 10m precision.
#'
#' @return
#' @export
#'
#' @examples
tapwater.getEuropeDt <- function() {
  world <- readOGR(dsn=paste0(ROOT_DATA_DIR, "shp/world"), layer = "ne_10m_admin_0_countries")
  europe <- world[world$REGION_UN=="Europe",]
  europe <- spTransform(europe, CRS("+init=epsg:2154")) # Lambert 93
  europeDt <- data.table(tidy(europe))
  europeDt <- europeDt[, piece:= as.numeric(piece)]
  
  # Uncomment lines below to have a quick view of the map.
  # ggplot() +
  #   geom_path(data = europeDt[id == 22],
  #                mapping = aes(x = long, y = lat, group = group),
  #                color = "black") +
  #   coord_fixed() +
  #   xlim(c(100000, 1300000)) +
  #   ylim(c(6000000, 7200000))
  return(europeDt)
}

# -####################################################################################################################-
#' Get interpolated isotopic composition, on the lat/lon grid defined by gridSize
#'
#' @param isotopeDt `data.table` containing 3 columns: `lat`, `lon` and `value`. Localized data points which should be 
#' interpolated over a grid.
#' @param gridSize Numeric. Meters between 2 grid cells' center (Lambert93 projection). Default is `6000` (6km). 
#'
#' @return `data.table`. Gridded data.
tapwater.getInterpData <- function(isotopeDt, gridSize) {
  lon = isotopeDt$lon
  lat = isotopeDt$lat
  value = isotopeDt$value
  
  gridX <- seq(min(lon) - gridSize, max(lon) + gridSize, by = gridSize)
  gridY <- seq(min(lat) - gridSize, max(lat) + gridSize, by = gridSize)
  
  inter_lin = interp(x = lon, 
                     y = lat,
                     z = value,
                     xo = gridX,
                     yo = gridY, 
                     duplicate = "mean",
                     linear = FALSE)
  
  interpDt <- data.table(lon = rep(inter_lin$x, length(gridY)),
                         lat = rep(inter_lin$y, each = length(gridX)),
                         value = as.vector(inter_lin$z))
  
  interpDt <- interpDt[!is.na(value)]
  
  return(interpDt)
}

# -####################################################################################################################-
#' Keep only `gridDt` data which are within `borderDt`.
#'
#' @param borderDt `data.table` containing lat/long points forming surfaces boundaries.
#' @param gridDt `data.table`. Grided data (Columns: lat, lon, value)
#'
#' @return `data.table`. Points of `gridDt` which are within `borderDt` boundaries.
tapwater.getInLandData <- function(borderDt, gridDt) {

  inLandDt <- data.table()
  # Loop over each "piece", one piece being a continuous border, i.e. the mainland or an island.
  for(this.piece in unique(borderDt$piece)) {
    landDt <- borderDt[piece == this.piece]
    bound <- tryCatch(
      {
        owin(poly = data.frame(x = landDt$lat,
                               y = landDt$long))
      },
      error=function(cond) { 
        # owin() expects object coordinates to be anticlockwise, and raise an
        # error if it is not the case.
        owin(poly = data.frame(x = rev(landDt$lat),
                               y = rev(landDt$long)))
      }
    )
    
    isInside <- inside.owin(x = gridDt$lat,
                            y = gridDt$lon,
                            w = bound)
    inLandDt <- rbind(inLandDt, gridDt[isInside])
  }
  
  # Uncomment the code below to get a quick insight of the `inLand` data.
  # ggplot() +
  #   geom_raster(aes(x = lon, y = lat), data = inLandDt) +
  #   geom_path(aes(x = long, y = lat, group = piece), data = borderDt, color = "red")
  
  return(unique(inLandDt))
}

# -####################################################################################################################-
#' Get base isotope data from file.
#'
#' @param species String. Name of the species to be plotted. Can be either `Delta_18_16`, `Delta_D_H` or `DX`.
#' @param measOrModel String. Data coming either from `measure` or `model`.
#' @param modelType String. If `measOrModel == "model"`, specify the model type: either `precip` for modelized 
#' precipitations, or `ground` for ground water modelisation.
#'
#' @return `data.table` containing 3 columns: `lat` and `lon` with the WGS84 coordinates of the data points, and `value`
#' with the isotope's numerical value.
#'
#' @examples tapwater.getIsotopeDt("Delta_18_16", "model", "ground")
tapwater.getIsotopeDt <- function(species, measOrModel, modelType) {
  # Get base data
  filename <- paste0(ROOT_DATA_DIR, "isotope/data_20200505.csv")
  isotopeDt <- fread(filename)
  
  # Filter by type and subtype
  isotopeDt <- isotopeDt[type == measOrModel & subtype == modelType]
  
  # Keep only the requested species
  setnames(isotopeDt, species, "value")
  isotopeDt <- isotopeDt[, .(lon, lat, value)]
  
  # Remove points without data
  isotopeDt <- isotopeDt[!is.na(lon)] # Remove points without GPS coordinates
  isotopeDt <- isotopeDt[!is.na(value)] # Remove points without Value
  
  return(isotopeDt)
}

# -####################################################################################################################-
#' Convenience function to properly round and sort `ggplot2`'s scale's limits and breaks
#'
#' @param values List of numeric. Values to be rounded.
#' @param digits Numeric. Number of digits which should be used to round the `values`.
#'
#' @return List of numeric. Rounded and sorted `values`.
#' @export
#'
#' @examples tapwater.roundScale(c(2.123, 0.001, 5.550, 10))
tapwater.roundScale <- function(values, digits = 1) {
  factor <- 10 ^ digits
  values <- sort(values)
  nbVal <- length(values)
  values[1] <- floor(values[1] * factor) / factor
  values[nbVal] <- ceiling(values[nbVal] * factor) / factor
  if (nbVal > 2) {
    values[2:nbVal-1] <- round(values[2:nbVal-1], digits = digits)
  }
  return(values)
}

# -####################################################################################################################-
#' Transform lat/lon coordinates (WGS84) into Lambert 93 coordinates, which is the official projection for France.
#'
#' @param datatable `data.table` containing at least a `lat` and a `lon` column.
#'
#' @return `data.table`. Same as input, with projection changed to Lambert 93.
#' @export
tapwater.toLambert93 <- function(datatable) {
  coordinates(datatable) <- c("lon", "lat")
  proj4string(datatable) <- CRS("+init=epsg:4326") # WGS 84
  projected <- spTransform(datatable, CRS("+init=epsg:2154")) # Lambert 93
  projectedDt <- as.data.table(projected)
  return(projectedDt)
}

