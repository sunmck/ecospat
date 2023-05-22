library(ecospat)
library(raster)
library(rgbif)
library(maptools)
library(devtools)
library(terra)
library(ade4)
library(biomod2)


# Global occurrence of European rabbit
ocGBIF <- occ_search(scientificName = "Oryctolagus cuniculus", 
                     hasCoordinate = TRUE,
                     basisOfRecord = "Preserved_Specimen",
                     limit = 10000,
                     fields = c("decimalLatitude",
                                "decimalLongitude", 
                                "year",
                                "country", 
                                "countryCode"))


ocOccs <- ocGBIF$data

coordinates(ocOccs) <- c("decimalLongitude",
                         "decimalLatitude")



wclim <- getData("worldclim", var="bio", res=10, path="./")
par(mar = c(0,0, 3, 1))
plot(wclim[["bio1"]], main = "Mean Annual Temperature [°F] (BIO1)")



ocOccs <- cbind(ocOccs, extract(wclim, ocOccs))
ocOccs <- ocOccs[complete.cases(data.frame(ocOccs)), ]

eurExt <- extent(c(-20,35,30,90))
ocEUR <- crop(ocOccs, eurExt)

ausExt <- extent(c(100,160,-45,-10))
ocAUS <- crop(ocOccs, ausExt)

data(wrld_simpl)
par(mar = c(1, 0, 0, 0))
plot(wrld_simpl, border = "gray80")
points(ocEUR, pch = 16, col = 2, cex = 0.3)
points(ocAUS, pch = 16, col = 4, cex = 0.3, add=T)



### Pre-modeling
# crop the environmental data to the native and invasive geographical ranges
eurEnvR <- crop(wclim, eurExt)
ausEnvR <- crop(wclim, ausExt)

eurEnvM <- getValues(eurEnvR)
ausEnvM <- getValues(ausEnvR)

# remove missing values
eurEnvM <- eurEnvM[complete.cases(eurEnvM), ]
ausEnvM <- ausEnvM[complete.cases(ausEnvM), ]

# produce global environmental background data
globalEnvM <- rbind(eurEnvM, ausEnvM)


pca.clim <- dudi.pca(globalEnvM, center = TRUE,
                     scale = TRUE, scannf = FALSE, nf = 2)

global.scores <- pca.clim$li

nativeLS.scores <-
  suprow(pca.clim,
         data.frame(ocEUR)[, colnames(globalEnvM)])$li   
invasiveLS.scores <-
  suprow(pca.clim,
         data.frame(ocAUS)[, colnames(globalEnvM)])$li

nativeEnv.scores <- suprow(pca.clim, eurEnvM)$li
invasiveEnv.scores <- suprow(pca.clim, ausEnvM)$li


data.frame(ocEUR)[, colnames(globalEnvM)]
# calculate the Occurrence Density Grid for both native and invasive species
nativeGrid <- ecospat.grid.clim.dyn(global.scores,
                                    nativeEnv.scores,
                                    nativeLS.scores)

invasiveGrid <- ecospat.grid.clim.dyn(global.scores,
                                      invasiveEnv.scores, 
                                      invasiveLS.scores)

ecospat.plot.niche.dyn(nativeGrid, invasiveGrid, quant = 0.1, interest = 2, title = "Niche Overlap", name.axis1 = "PC1", name.axis2 = "PC2")

# plot variable contributions
ecospat.niche.dyn.index(nativeGrid, invasiveGrid, intersection = 0.1)$dynamic.index.w


ocAUSearly <- subset(ocAUS, year <= 1950)
ocAUSlate <- subset(ocAUS, year > 1950)

geoGrid <- expand.grid(longitude =
                         seq(100, 160, length.out = 250),
                       latitude =
                         seq(-45, -10, length.out = 250))

mask <- subset(wrld_simpl, NAME == "Australia")

earlyGeoGrid <- ecospat.grid.clim.dyn(geoGrid, geoGrid,
                                      coordinates(ocAUSearly),
                                      geomask = mask)

lateGeoGrid <- ecospat.grid.clim.dyn(geoGrid, geoGrid,
                                     coordinates(ocAUSlate),
                                     geomask = mask)

ecospat.plot.niche.dyn(earlyGeoGrid, lateGeoGrid, quant = 0)
plot(wrld_simpl, add = TRUE)


# calculate niche overlap

ecospat.niche.overlap(nativeGrid, invasiveGrid, cor=T)

# perform the Niche Equivalency Test

eq.test <- ecospat.niche.equivalency.test(nativeGrid, invasiveGrid, rep = 100, ncores = 2)

# perform the Niche Similarity Test

sim.test <- ecospat.niche.similarity.test(nativeGrid, invasiveGrid, rep = 100, rand.type = 2, ncores = 2)

# plot Equivalency and Similarity Test

par(mfrow=c(1,2))
ecospat.plot.overlap.test(eq.test, "D", "Equivalency") 
ecospat.plot.overlap.test(sim.test, "D", "Similarity")


# gridding the native niche

grid.clim.t.nat <- ecospat.grid.clim.dyn(glob = globalEnvM[,1],
                                         glob1 = data.frame(eurEnvM[,1]),
                                         data.frame(ocEUR)[,4], R = 1000, th.sp = 0)

# gridding the invasive niche

grid.clim.t.inv <- ecospat.grid.clim.dyn (glob = globalEnvM[,1], 
                                          glob1 = data.frame(ausEnvM[,1]), 
                                          data.frame(ocAUS)[,4], R = 1000, th.sp = 0)

t.dyn <- ecospat.niche.dyn.index (grid.clim.t.nat, grid.clim.t.inv, intersection=0.1)

ecospat.plot.niche.dyn(grid.clim.t.nat, grid.clim.t.inv, quant=0.1, interest=2, title= "Niche Overlap", name.axis1="Average temperature")

# showing the shift of the niche centroid along the temperature gradient (compared to the shift of the available climate in the study area)

ecospat.shift.centroids(data.frame(ocEUR)[,4],
                        data.frame(ocAUS)[,4],
                        data.frame(eurEnvM)[,1],
                        data.frame(ausEnvM)[,1])


### ESM modeling

# format data using the biomod2 package
ocEUR$occ <- 1
ocEUR_coords <- data.frame(coordinates(ocEUR))

nat.biomod <- BIOMOD_FormatingData(resp.var = as.numeric(data.frame(ocEUR)[,23]),
                                   PA.strategy = "random", # pseudo absence selection necessary because only occurrence and no absence                                       data is available
                                   PA.nb.rep = 1, # number of repetitions of pseudo-absence points that are drawn
                                   PA.nb.absences = 1000, # number of pseudo-absence points that will be selected
                                   expl.var = eurEnvR, # environmental data (here mean annual temperature)
                                   resp.xy = ocEUR_coords,
                                   resp.name = "Oryctolagus cuniculus")

plot(nat.biomod)

biomodopt <- bm_DefaultModelingOptions()

# Calibration of simple bivariate models
nat.ESM <- ecospat.ESM.Modeling(data = nat.biomod,
                                models = c('GLM'),
                                NbRunEval = 2,
                                DataSplit = 70,
                                weighting.score = c("AUC"),
                                models.options = biomodopt)

# Ensemble models
nat.ESM.ens <- ecospat.ESM.EnsembleModeling(nat.ESM, 
                                            weighting.score = c("SomersD"), 
                                            threshold = 0)

# Projection of simple bivariate models into new space or time
nat.ESM_proj_current <- ecospat.ESM.Projection(ESM.modeling.output = nat.ESM,
                                               new.env = ausEnvR)

# Projection of calibrated ESMs into new space or time
nat.ESM_proj_current$new.env.raster = TRUE
nat.ESM.ens_proj_current <- ecospat.ESM.EnsembleProjection(ESM.prediction.output = nat.ESM_proj_current,
                                                           ESM.EnsembleModeling.output = nat.ESM.ens,
                                                           chosen.models = 'all')


