library("raster")

RR <- extent(-26,48,35,72)

# Loading the rasters
setwd("/Users/stefanomammola/Desktop/PAPERS IN CORSO/Mammola_et al_spider_traits_macroecology_europe/raster")
predictors <- stack(list.files(pattern="*.tif", full.names=TRUE))
projection(predictors) <-"+proj=longlat +datum=WGS84 +ellps=WGS84"
names(predictors) <- c("elev",
                       "ice",
                       "karst",
                       "T_mean",
                       "annual_range",
                       "prec"
                       )
predictors <- crop(predictors,RR)
