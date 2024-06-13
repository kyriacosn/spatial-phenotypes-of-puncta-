library(sf)
library(INLA)
library(inlabru)
library(fmesher)
library(raster)
library(mgcv)
library(ggplot2)
library(terra)
bru_safe_sp(force = TRUE)


# ----------------------------------------------------------------
#         read the data and choose cells for analysis
# ---------------------------------------------------------------

# Read GeoJSON files into an sf objects

all.peroxisomes <- st_read("geo/peroxisomes.geojson",crs = NA)
all.lipid <- st_read("geo/lipid.geojson",crs = NA)
all.boundaries <- st_read("geo/domain.geojson",crs = NA)
all.meshboxes <- st_read("geo/meshbox.geojson",crs = NA)

# Read covariate tifs as raster objects

rast.mt <- rast("geo/mitochondria.tif")
rast.er <- rast("geo/er.tif")
rast.nd <- rast("geo/nucleus_distance.tif")
rast.fd <- rast("geo/fractional_distance.tif")


# Chose the cells to fit the model on 

fit_cells <- c(1,2,4)

# Filter the domain boundaries, mesh boxes and points inside the domains of fit_cells

boundary <- all.boundaries[fit_cells, ]

meshboxes <- all.meshboxes[fit_cells, ]

inside.indx <- apply(st_intersects(all.peroxisomes, boundary, sparse = FALSE),1, function(x) any(x))
peroxisomes <- all.peroxisomes[inside.indx, ]  

# Plot check
ggplot() + gg(meshboxes)+ gg(boundary) + gg(peroxisomes)


# construct the mesh

mesh <- fm_mesh_2d_inla(boundary = list(boundary,meshboxes), max.edge = c(2, 4),cutoff = 0.2)

# plot check the mesh

ggplot() + gg(mesh) + gg(peroxisomes)

# ----------------------------------------------------------------
#                   model construction
# ----------------------------------------------------------------


# A matern random field with:
#   prior of the standard deviation to be more than 2 with probability 0.01
#   prior and the range less than 1um with probability 0.01

matern <- inla.spde2.pcmatern(mesh, prior.sigma = c(2, 0.01), prior.range = c(1, 0.01))

model<- geometry ~ 
  random.field(geometry, model = matern) +
  mitochondria(rast.mt, model = "linear") +    
  er(rast.er, model = "linear") +
  fractional.distance(rast.fd, model = "linear") +
  Intercept(1)


# ----------------------------------------------------------------
#                     fit the model
# ----------------------------------------------------------------

# fit the model
fit <- lgcp(model, peroxisomes ,samplers = boundary , domain = list(geometry = mesh))


summary(fit)

#inspect the marginal distributions



# ------------------------------------------------------------------
#  predict the mean density of cells based on covariates
# ------------------------------------------------------------------

# choose the cells to predict the mean density
predict_cells <- c(3)

# Filter the domain boundaries, mesh boxes and points inside 
# the domains of predict_cells
predict.boundary <- all.boundaries[predict_cells, ]

predict.meshboxes <- all.meshboxes[predict_cells, ]

predict.inside.indx <- apply(st_intersects(all.peroxisomes, predict.boundary, sparse = FALSE),1, function(x) any(x))

predict.peroxisomes <- all.peroxisomes[predict.inside.indx, ]  


# Plot check
ggplot() +
  gg(predict.meshboxes)+
  gg(predict.boundary) +
  gg(predict.peroxisomes)


prediction.mesh <- fm_mesh_2d_inla(boundary = list(predict.boundary,predict.meshboxes), max.edge = c(2, 4),cutoff = 0.2)

prediction.grid <- fm_pixels(prediction.mesh, mask = predict.boundary)

Sys.time()
predicted.intensity <-predict(fit, prediction.grid, ~ exp(  mitochondria + er + fractional.distance + Intercept), n.samples = 200)
Sys.time()


# plot the predicted mean density 
p <- ggplot() +
  gg(predicted.intensity , geom = "tile") +
  gg(predict.peroxisomes,color = "green") +
  geom_segment(aes(x = 160, y = 90, xend = 170, yend = 90), color = "black", size = 1)+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),  # Removes all axis text (numbers)
    axis.title.x = element_blank(),  # Removes x-axis label/title
    axis.title.y = element_blank()   # Removes y-axis label/title
  )

  
ggsave("cell3_predict.pdf", plot = p, width = 7, height = 3, units = "in")

