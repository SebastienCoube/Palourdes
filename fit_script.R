############################
#Chargement des données...#
##########################

source("Random_field_simulator.R")
source("functions.R")

##... de campagne palourde 2018
library(raster)
data_palourdes <- read.table("fichier global palourde 2003 a maintenant.txt", header = TRUE)
data_palourdes<-data_palourdes[which(data_palourdes$Annee == "2018"),c(1,2,6:8)]
colnames(data_palourdes) = c("X","Y", "g_0.25m2", "point", "Annee")
xy_obs<-cbind(data_palourdes[,1],data_palourdes[,2]) 
data_obs<-as.numeric(data_palourdes[,3])


#coordinates(data_palourdes) <- ~X+Y

##... des strates (en shapefile)
require(rgdal)
bassin <- readOGR("./bassin entier/bassin entier.shp")
bassin <- spTransform(bassin, CRS("+init=epsg:2154")) #lambert 93


# On découpe le shapefile en une grille vide de cellules de 50*50. Comme on est en projection L93, 
# ca représente à peu près des carrés de 50*50m

bassin_grid <- expand.grid(latcoords = seq(from = 366774.4, to = 380174.4, by = 50),
                           lngcoords = seq(from = 6402304.4, to = 6415554.4, by = 50))

grid.pts<-SpatialPointsDataFrame(coords= bassin_grid, data=bassin_grid, proj4string = CRS("+init=epsg:2154"))

plot(grid.pts[bassin,])

coord_bassin <- cbind(grid.pts[bassin,]@coords[,1], grid.pts[bassin,]@coords[,2])


## Pour résumer : coord_bassin = la matrice pour l'interpolation 
## donnesxyz = données ponctuelles à partir desquelles il faut faire le random field. 


fit = fit_mcmc_Vecchia(observed_locs = xy_obs, predicted_locs = coord_bassin, observed_field = data_obs, 
                 n_iterations =  3, m = 12, n_cores = 1, n_chains = 3, field_thinning = .05,  n_delayed_acceptance = 4041, 
                 n_join = 10,
                 pc_prior_range = c((max(xy_obs[,2])-min(xy_obs[,2]))/10, .5),
                 pc_prior_sd = sd(data_obs, .5))
