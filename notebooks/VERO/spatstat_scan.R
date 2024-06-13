library(spatstat)
library(sp)
library(jsonlite)
library(terra)
library(raster)



points = list()

IM = list()
IE = list()
IG = list()

PM = list()
PE = list()
PG = list()

DM = list()
DN = list()
DN2 = list()

BM = list()
BE = list()
BG = list()

cell_names = c("cell1","cell2","cell3","cell4","cell5","cell6")
for (cell in cell_names){

  # Reading the JSON file
  celldata <- fromJSON(paste0(cell,'/preprocess.json'))
  scale <- celldata$scale
  celldata$covariates$binary <- lapply(celldata$covariates$binary, function(listimg) im(listimg, xrange = c(0, ncol(listimg)-1)/scale, yrange = c(0, nrow(listimg)-1)/scale))
  celldata$covariates$proximity <- lapply(celldata$covariates$proximity, function(listimg) im(listimg, xrange = c(0, ncol(listimg)-1)/scale, yrange = c(0, nrow(listimg)-1)/scale))
  celldata$covariates$intensity <- lapply(celldata$covariates$intensity, function(listimg) im(listimg, xrange = c(0, ncol(listimg)-1)/scale, yrange = c(0, nrow(listimg)-1)/scale))
  celldata$covariates$distance <- lapply(celldata$covariates$distance, function(listimg) im(listimg, xrange = c(0, ncol(listimg)-1)/scale, yrange = c(0, nrow(listimg)-1)/scale))

    # CREATE SPATSTAT OBJECTS
  
  nuc_con <- cbind(celldata$contour$nucleus$x,celldata$contour$nucleus$y)
  cell_con <- cbind(celldata$contour$cell$x,celldata$contour$cell$y)
  cell_con <- cell_con[nrow(cell_con):1, ]
  poly_list <- list(cell_con, nuc_con)
  window_with_hole <- owin(poly = poly_list)
  
  points[[cell]] <- ppp(celldata$`point patterns`$peroxisomes$x, celldata$`point patterns`$peroxisomes$y, window=window_with_hole)

  IM[[cell]] <-  celldata$covariates$intensity$mitochondria
  IE[[cell]] <-  celldata$covariates$intensity$er
  IG[[cell]] <-  celldata$covariates$intensity$golgi
  
  PM[[cell]] <-  celldata$covariates$proximity$mitochondria
  PE[[cell]] <-  celldata$covariates$proximity$er
  PG[[cell]] <-  celldata$covariates$proximity$golgi

  BM[[cell]]  <-  celldata$covariates$binary$mitochondria
  BG[[cell]]  <-  celldata$covariates$binary$golgi
  BE[[cell]] <-  celldata$covariates$binary$er
  
  DM[[cell]] <-  celldata$covariates$distance$mitochondria
  DN[[cell]]  <-  celldata$covariates$distance$nucleus
  DN2[[cell]]  <-  celldata$covariates$distance$nucleus2
  
}

fitall <- mppm(Points ~ 1 + BE + BM  + BG + PM + PG + DM + DN + DN2 ,hyperframe(Points = points,IM=IM,IE=IE,IG=IG,BM = BM,BE=BE,BG=BG,PM = PM,PE=PE,PG=PG, DN = DN,DN2=DN2,DM = DM))

stepfit <- step(fitall,direction = c("both"))

