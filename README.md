# Niche Modelling with Ecospat in R
#### authors: Sunniva McKeever, Isabella Metz and Janik Hoffmann
#### 2023-05-23


This repository documents the presentation of the ecospat R package within the corresponding seminar of the EAGLE module 'Spatial Modeling and Prediction'. It mainly follows the guidance of Di Cola et al. 2016. 


### The ecospat package

The ecospat package aims to provide novel tools and methods for spatial analyses and the modeling of species niches and distributions in a coherent workflow. It includes functions for species niche quantification and measures of phylogenetic diversity among others. As core modeling features it contains functionalities to use the ESM (Ensemble of Small Models) approach and various implementations of the spatially-explicit modeling of species assemblages (SESAM) framework. Post-modeling analyses include the evaluation of species predictions based on presence-only data via the Boyce Index. Ecospat also offers functions to supplement the 'biomod2' package.


### Goals of the seminar

The goals of this seminar are to give an understanding of the core functionalities of the novel ecospat package by going through to example applications. The first one deals with the quantification of environmental niches and the ESM approach, whereas the second one addresses functions to make spatial analyses of a species community.


## Part 1 Niche modeling of European rabbits


In this part, we obtain species occurrence data of the European rabbit in its native range (Europe) and invasive range (Australia) and perform a niche quantification and spatial prediction of the species distribution in its invasive geographical range by applying the ESM approach.

More information on how European rabbits took over Australia can be found [here](https://education.nationalgeographic.org/resource/how-european-rabbits-took-over-australia/). 

### Background on the niche concept

The ecological niche represents a key concept for determining species distributions and planning biodiversity conservation strategies. It links the geographical occurrence of species with the occupied environmental conditions. However, the niche conservatism casts doubt on whether niches remain constant across time and space. Biological invasions represent a unique case to test this hypothesis.

In regard to SDMs, the assumption is that species niches do not change much across space, in other words, that species occupy the same environmental conditions in new geographical ranges. Apart from that, recent studies observed niche shifts for some terrestrial species and invasive species can, in fact, expand their niches as a consequence of shifts in their realized and/or fundamental niche. Being introduced to new ranges, invasive species often encounter new biotic interactions and dispersal limitations leading to adaption processes to exotic environments (Liu et al. 2020).

SDMs enable to quantify niche dynamics by being calibrated with the species occurrence and environmental data in one range which is then applied to predict the distribution in a new range.


### Data download

Load the required libraries

```{r libraries, warning=FALSE, message=FALSE}
library(ecospat)
library(raster)
library(rgbif)
library(maptools)
library(devtools)
library(terra)
library(ade4)
library(ape)
library(biomod2)
```

Download species data from **GBIF**

```{r gbif, warning=FALSE}
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
```

Extract the records from the GBIF data and assign the coordinates to map observations. This will automatically convert the record matrix to a `SpatialPointsDataFrame` object.

```{r gbif data,  warning=FALSE}
ocOccs <- ocGBIF$data

coordinates(ocOccs) <- c("decimalLongitude",
                         "decimalLatitude")
```

Download the climate data from the worldclim 2.1 database by calling the `getData` function from the raster package. These will later on serve as niche descriptors and model predictors.

```{r worldclim, message=FALSE, warning=FALSE}
wclim <- getData("worldclim", var="bio", res=10, path="./")
par(mar = c(0,0, 3, 1))
plot(wclim[["bio1"]], main = "Mean Annual Temperature [°F] (BIO1)")
```

Extract environmental data for the occurrence records and crop these to the extent of the native and invasive geographical ranges. Then plot the point records on a worldmap using the maptools package.

```{r records, message=FALSE, warning=FALSE}
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
```
![datadownload](https://github.com/sunmck/ecospat/assets/116874799/38bba10c-d671-40e4-a4ba-8655a61dc6e8)

### Pre-modeling

The environmental variables need to be clipped to the respective extents of the native and invasive geographical ranges. So in this case to the extents of Europe and Australia. For the niche quantification, a matrix with the background environmental variables from both ranges, as well as the global environment are needed. After cropping the rasters, `getValues`is applied to convert them to a data frame.

```{r pre modeling I,  warning=FALSE}
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
```

**Statistical pre-analysis**

Before the application of spatial ecology methods, a statistical pre-analysis will be conducted to explore the data used in the further course. 

For a matter of applicability, the environmental variables are renamed according to what they actually mean. The wordclim temperature variables are in '°C*10', to avoid this confusion, we also recalculate the temperature values in a loop.

```{r rename&recalculate,  warning=FALSE}
# rename

worldclim2.1_variables <- c(
                   "ann.mean.temp", 
                   "diurnal.temp.range.mean",
                   "isothermality",
                   "temp.seasonality",
                   "max.temp.warmest_month",
                   "min.temp.coldest_month",
                   "ann.temp.range",
                   "mean.temp.wettest_quarter",
                   "mean.temp.driest_quarter",
                   "mean.temp.warmest_quarter",
                   "mean.temp.coldest_quarter",
                   "ann.precip",
                   "precip.wettest_month",
                   "precip.driest_month",
                   "precip.seasonality",
                   "precip.wettest_quarter",
                   "precip.driest_quarter",
                   "precip.warmest_quarter",
                   "precip.coldest_quarter",
                   "x",
                   "y"
                   )

ocEUR.df <- as.data.frame(ocEUR)
ocAUS.df <- as.data.frame(ocAUS)
ocEUR.df <- ocEUR.df[, 4:24]
ocAUS.df <- ocAUS.df[, 4:24]
colnames(ocEUR.df) <- worldclim2.1_variables
colnames(ocAUS.df) <- worldclim2.1_variables


# recalculate 

for (col in 1:ncol(ocEUR.df)) {
  # Check if column name contains "temp"
  if (grepl("temp", colnames(ocEUR.df)[col])) {
    # Divide column values by 10
    ocEUR.df[, col] <- ocEUR.df[, col] / 10
  }
}

for (col in 1:ncol(ocAUS.df)) {
  # Check if column name contains "temp"
  if (grepl("temp", colnames(ocAUS.df)[col])) {
    # Divide column values by 10
    ocAUS.df[, col] <- ocAUS.df[, col] / 10
  }
}
```

In a next step, we want to calculate basic statistics for each of the climatic variables for both the native and invasive ranges.

```{r statistics,  warning=FALSE}
occEURstats <- apply(ocEUR.df, 2, function(x) c(min = min(x), median = median(x), mean = mean(x), max = max(x), sd = sd(x)))
occAUSstats <- apply(ocAUS.df, 2, function(x) c(min = min(x), median = median(x), mean = mean(x), max = max(x), sd = sd(x)))

# combine both statistics for direct comparison
bioStats <- rbind(occEURstats, occAUSstats)
head(bioStats[, c(1, 7, 12)], 10)
```

Spatial autocorrelation can be measured by applying the `ecospat.mantel.correlogram` function. On the y-axis a spatial autocorrelation measure (ranging from -1 to 1) is plotted against distance intervals on the x-axis (km). If the coefficient is above zero, similar values are more closer together, a negative coefficient implies that dissimilar values tend to be close. If the correlation coefficients are highly positive, it indicates increasing spatial clustering. If the correlation coefficient decreases by distance, it shows spatial segregation, whereas if it it remains positive, it indicates a spatial gradient. The latter is the case for our observations.

```{r spatial autocorr, warning=FALSE}
# native range

ecospat.mantel.correlogram(dfvar=ocEUR.df[c(1:21)],colxy=20:21, n=100, colvar=1:19, 
max=10, nclass=10, nperm=100)

# invasive range

ecospat.mantel.correlogram(dfvar=ocAUS.df[c(1:21)],colxy=20:21, n=100, colvar=1:19, 
max=10, nclass=10, nperm=100)
```

## Core-modeling

**Niche quantification**

For the niche quantification, a Principal Component Analysis (PCA) of the environmental data is carried out.

```{r PCA, results=FALSE, warning=FALSE}
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
```

To explain, `dudi.pca`conducts the PCA on the global data (`wclim_M`). The result is a two-dimensional summary of the total environmental variability. In a next step, the observation data (`nat_occ` and `inv_occ`) are mapped into that two-dimensional space using the `suprow` function. From the output only the ´li´ element is important at the moment. 

Now, we can compare the niches of both species types. For that the PCA scores of the global data, the native and invasive environments and the native and invasive occurrence records are used.

```{r niche quantification, message=FALSE, results=FALSE, warning=FALSE}
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

ecospat.plot.contrib(contrib=pca.clim$co, eigen=pca.clim$eig)
```
![occurrencegrid](https://github.com/sunmck/ecospat/assets/116874799/6be3bb34-4b1b-4058-a803-06e9e37985d2)

The plot shows the environmental conditions in Europe (green line) and Australia (red line). The green area represents the native environmental range, whereas the red area shows the invasive range. The blue area elucidates environments occupied in both ranges, called niche overlap. 

Next, we conduct a niche dynamics analysis (niche categories, climate analogy between ranges) to investigate mechanisms behind the niche differences by applying the `ecospat.niche.dyn.index()` function. The result indicates an expansion factor of around 0.3, a stability factor of 0.68 and no niche unfilling. This elucidates that the species is able to expand into novel climates and that it has already occupied the full native environmental range within the new geographical range.

```{r niche dynamics, warning=FALSE}
ecospat.niche.dyn.index(nativeGrid, invasiveGrid, intersection = 0.1)$dynamic.index.w
```

Instead of environmental conditions, the analysis can be applied to geographical ranges. This does not make sense in the case of native and invasive species, however, it could be worth a try to unveil temporal niche shifts of the invasive European rabbit in Australia. For that the distributions before and after 1950 used.

```{r geographic comparison, warning=FALSE}
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
```

![earlylate](https://github.com/sunmck/ecospat/assets/116874799/610f02c8-be1c-4806-bd81-77c7f123a8e3)


**Further methods and tests**

Measures of niche overlap are used to identify potential overlaps or changes in the niches of, for instance, invasive species. The niche equivalency test and the niche similarity test allow to quantify the statistical significance of detected niche differences against null model niches taken randomly from the background data. 

The function ´ecospat.niche.overlap` takes the differences in occurrence densities to measure Schöner's D index, giving a value between 0 (no overlap) and 1 (full overlap). The niche equivalency test randomly shuffle the occurrences between the two ranges and assesses how equal the niches are. A p-value below 0.05 imply that the niches are significantly more different than expected by random models. The niche similarity tests conducts random shifts within available conditions in the study area to proove if the niches are more or less similar than expected by chance (Di Cola et al. 2017). If the p-value is below 0.05, then the niches are significantly more similar.

The outcome of the three tests below indicate that the native and invasive niches are significantly more equivalent than expected by chance. 

```{r tests, warning=FALSE}
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
```

Now we want to plot the niche dynamics against one gradient (here temperature) with `ecospat.plot.niche.dyn()`

```{r one gradient, warning=FALSE}
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
```
![avgtemp](https://github.com/sunmck/ecospat/assets/116874799/03200aa3-1eec-426e-b084-90d1fbf41e1b)

**Ensemble of Small Models (ESM)**

The implementation of the ESM approach is another key feature of the ecospat package. It is specifically valuable when only a limited number of species records are available (rare species). The SDM algorithm fits multiple small models with at least two variables at each time and sums the predictions weighted by each submodel performance. In doing so, overparametrization and overfitting is prevented. 

```{r ESM, message=FALSE, warning=FALSE, results='hide'}
# format data using the biomod2 package

ocEUR$occ <- 1
ocEUR_coords <- data.frame(coordinates(ocEUR))

eurEnvR <- stack(eurEnvR)

nat.biomod <- BIOMOD_FormatingData(resp.var = as.numeric(data.frame(ocEUR)[,23]),
                                   PA.strategy = "random",
                                   PA.nb.rep = 1,
                                   PA.nb.absences = 1000,
                                   expl.var = eurEnvR,
                                   resp.xy = ocEUR_coords,
                                   resp.name = "Oryctolagus cuniculus")

biomodopt <- bm_DefaultModelingOptions()

# calibration of simple bivariate models

my.ESM <- ecospat.ESM.Modeling(data = nat.biomod,
                                models = c('GLM'),
                                NbRunEval = 2,
                                DataSplit = 70,
                                weighting.score = c("AUC"),
                                models.options = biomodopt)

# Ensemble models

my.ESM_EF <- ecospat.ESM.EnsembleModeling(my.ESM, 
                                          weighting.score = c("SomersD"), 
                                          threshold = 0)

ausEnvR <- stack(ausEnvR)


# projection of simple bivariate models into new space or time

my.ESM_proj_current <- ecospat.ESM.Projection(ESM.modeling.output = my.ESM,
                                              new.env = ausEnvR)

# projection of calibrated ESMs into new space or time

my.ESM_EFproj_current <- ecospat.ESM.EnsembleProjection(
  ESM.prediction.output = my.ESM_proj_current,
  ESM.EnsembleModeling.output = my.ESM_EF)



```
```{r ESM plot, message=FALSE, warning=FALSE}
plot(my.ESM_EFproj_current)
plot(ocAUS, add = T)
```

![ESM](https://github.com/sunmck/ecospat/assets/116874799/71f93881-0608-4256-b256-6fcbe15b5c6a)

## Part 2: Community properties of plant assemblages in the Swiss Alps

In this part, we use the ecospat package to perform a spatial analysis of a plant community in the Swiss Alps by investigating their phylogentic diversity and co-occurrence patterns and by predicting the spatial distribution of assemblages. The former step includes functions that represent unique implementations of the SESAM framework.


### Species assemblage's structure and spatial predictions

**Phylogenetic diversity predictions**

Here, we predict phylogenetic diversity, that is basically the summed evolutionary age of all species in the community measured by the phylogenetic tree. Species richness refers to the numbers of species within an community.

```{r phylogenetic diversity, warning=FALSE}
# load the tree and data set from the ecospat package

fpath <- system.file("extdata", "ecospat.testTree.tre", package = "ecospat")

tree <- read.tree(fpath) 

data <- ecospat.testData[9:52]
```


```{r calculate phylogenetic diversity, warning=FALSE}

# plot the correlation of phylogenetic diversity against the species richness

pd <- ecospat.calculate.pd(tree, data, method = "spanning", type = "species", root = FALSE, average = FALSE, verbose = FALSE )

plot(pd)
```

**SESAM framework: spatial prediction of communities**

Commonly, multiple SDMs of individual species have been stacked to obtain species richness predictions. A novel approach combines this method by a coupling with richness predictions. This is one step in the SESAM framework aiming to reconstruct species assemblages by implementing a couple of filters including the application of a biotic rule which is realized in ecospat. This 'probability ranking rule' (PRR) decides which species should be included in the final richness prediction. In doing so, it ranks the species in decreasing order of predicted probability of occurrence and selects species down the list until expected richness value is met.

We predict community composition with the function `ecospat.SESAM.prr()`. The function returns a data frame including binary predictions of each species per site.

```{r sesam I, results=FALSE, warning=FALSE}
# Two dataframes are required: (1) continuous probabilities from SDMs for all species (columns) in the considered sites (rows) and (2) richness value for each site (first column)

# (1)
proba<-ecospat.testData[,73:92]

#(2)
sr<-as.data.frame(rowSums(proba))

# run the function
ecospat.SESAM.prr(proba, sr)
```

**Evaluating community predictions**

Now we calculate the accuracy of the community predictions with different indices by applying the function `ecospat.Community.Eval()`that requires as input the species records (presence/absence) and their predictions. The result is a list of evaluation metrics calculated for each site

```{r evaluate community predictions, results=FALSE, eval=FALSE, warning=FALSE}
eval <- ecospat.testData[c(53,62,58,70,61,66,65,71,69,43,63,56,68,57,55,60,54,67,59 ,64)]

pred<-ecospat.testData[c(73:92)]

ecospat.CommunityEval (eval, pred, proba=TRUE, ntir=1)
```


**Environmentally-constrained species co-occurrence analyses**

In the end, a co-occurrence analysis is performed to test for non-random patterns of species co-occurrences by applying an environmentally-constrained null model with the function `ecospat.cons_Cscore()`. The results are the C-score index for the observed community, the mean of C-score for the simulated communities, the p-values to evaluate whether the difference of the former two are significant and returns the standardized effect size (SES) for the whole community. 

The C-score is a measure of randomness of the distribution of species pairs across a number of biomes. A high C-score means that there is a higher probability that the distribution of one species has been directly affected by the presence of another species. 

```{r co-occurrence analyses, warning=FALSE}
presence <- ecospat.testData[c(9:24)]
pred <- ecospat.testData[65:82]

nbpermut <- 10000

outpath <- getwd()

ecospat.cons_Cscore(presence, pred, nbpermut, outpath)
```



### References

Di Cola et al. 2016. ecospat: an R package to support spatial analyses and modeling of species niches and distributions. Ecography, 40, 6 (2017). 10.1111/ecog.02671

Liu et al. 2020. Most invasive species largely conserve their climatic niche. PNAS, 117, 38. 10.1073/pnas.2004289117




**Quiz Time**

https://isabellam.shinyapps.io/Quiz/ 








