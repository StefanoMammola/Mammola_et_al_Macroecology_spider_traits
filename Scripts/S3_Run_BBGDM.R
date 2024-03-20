## ------------------------------------------------------------------------
## ' Functional convergence underground? The scale-dependency of community assembly processes in European cave spiders
## ------------------------------------------------------------------------

#Null model  script
# rsync -r caioroza@puhti.csc.fi:/scratch/project_2005062/Mammola_et_al_Macroecology_spider_traits "/Users/graco-roza/OneDrive - University of Helsinki/Ongoing manuscripts/Mammola_et_al_Macroecology_spider_traits" caioroza@puhti.csc.fi:/scratch/project_2005062/
# cd /scratch/project_2005062/Mammola_et_al_Macroecology_spider_traits/
###### Set up of the script ############################################################################################

.libPaths(c("/projappl/project_2005062/project_rpackages_4.0.5", .libPaths()))
# Packages .............................................................................................................
if(!require("pacman")) {install.packages("pacman")}
pacman::p_load("tidyverse","hypervolume","doSNOW", "future","gdm", "surveillance","plotrix", "glue", "magrittr")
#....................................................................................................................... 

#Program cores for the batch ...........................................................................................
options(future.availableCores.methods = "Slurm")

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
jj <- as.numeric(slurm_arrayid)

#....................................................................................................................... 

#directories to save and load files ....................................................................................
input_dir<- "Data/" #Folder where the data are located
object_dir<-"objects/" #Sub-folder where data are located
output_dir<- "objects/" #Folder where to save the data
raster_dir <- "rasters/" # Folder where the raster data are stored
results_dir <- "Results/" #Folder where to save figures, tables, etc 
table_dir<-"tables/" # Folder where the excel files are saved


#....................................................................................................................... 

# Start the script for each dataset ====================================================================================


#### Load Functions
source("Scripts/Functions/BBGDM_functions.R")

load(paste0(input_dir,object_dir,"BBGDM_input.R"))

####### GET FOCAL DISSIMILARITY METRIC & ADD COLUMN OF SITE NAMES 

comm <- comm_parsed[rowSums(comm_parsed) > 2,] 

ntasks<- nrow(comm)

if (jj == 1) {
  trait <- trait_axis
} else {
  trait <- trait_axis[sample(nrow(trait_axis)), ]
}
rownames(trait) <-  rownames(trait_axis)

cores<-6
cl <- makeSOCKcluster(cores)
registerDoSNOW(cl)

alpha.FD <- foreach(
  i = 1:ntasks,
  .combine = hypervolume::hypervolume_join,
  .multicombine = TRUE,
  .errorhandling = 'remove'
) %dopar% {
  .libPaths(c("/projappl/project_2005062/project_rpackages", .libPaths()))
  if(!require("hypervolume")) {install.packages("hypervolume")}
  abun <- comm[i, which(comm[i, ] > 0)]

  hv <- hypervolume::hypervolume_box(data= trait[names(comm[i, which(comm[i, ] > 0)]),1:3],
                                     verbose = FALSE,
                                     name = rownames(comm)[i]
                                     #kde.bandwidth = estimate_bandwidth(trait[names(comm[i, which(comm[i, ] > 0)]),1:3])
                                     )
}
stopCluster(cl)

#get number of tasks for beta diversity 
name_sites<-sapply(seq_along(alpha.FD@HVList), function(i) alpha.FD@HVList[[i]]@Name)


ntasks<- length(name_sites)

#set parallel clusters again
cl <- makeSOCKcluster(cores)
registerDoSNOW(cl)
{
pairwise_beta <-  foreach(i=1:ntasks,
                          .combine=appendList) %:%
  foreach(j=i:ntasks, .combine=appendList) %dopar% {
    if(!require("hypervolume")) {install.packages("hypervolume")}
    hyperSet <- hypervolume::hypervolume_set(
      alpha.FD@HVList[[i]], alpha.FD@HVList[[j]],      
      check.memory = FALSE, verbose = FALSE, num.points.max = 1000)
    union <- hyperSet[[4]]@Volume
    unique1 <- hyperSet[[5]]@Volume
    unique2 <- hyperSet[[6]]@Volume
    
    union <- 2 * union - unique1 - unique2
    Btotal <- (unique1 + unique2)/union
    Brepl <- 2 * min(unique1, unique2)/union
    Brich <- abs(unique1 - unique2)/union
    output<-list(Btotal=Btotal,Brepl=Brepl,Brich=Brich)
    return(output)
  }    
stopCluster(cl)
beepr::beep()  
}  
output_beta<-list()
output_beta$Btotal <- matrix(NA,ntasks,ntasks)
output_beta$Brepl <- matrix(NA,ntasks,ntasks)
output_beta$Brich <- matrix(NA,ntasks,ntasks)

output_beta$Btotal[lower.tri(output_beta$Btotal,diag=TRUE)] <-round(pairwise_beta$Btotal,3)
output_beta$Brepl [lower.tri(output_beta$Btotal,diag=TRUE)] <- round(pairwise_beta$Brich,3)
output_beta$Brich [lower.tri(output_beta$Btotal,diag=TRUE)] <- round(pairwise_beta$Brepl,3)
Fbeta <- lapply(output_beta, as.dist)
names(Fbeta$Btotal) <- names(Fbeta$Brepl) <- names(Fbeta$Brich) <- name_sites

expl.var<- predictors %>%  
  rename(site = ID,
         x= decimalLongitude,
         y = decimalLatitude) 


model <- lapply(Fbeta, function(x) {
bbgdm(data= make_tab(x, pred = expl.var), geo=TRUE, bootstraps=999, ncores=6, knots=NULL, splines= NULL)
})

###### EXTRACT OUTPUTS

###### Get Wald Test results and signs of predictors
var_names<- c("Geographic","entranceSize","development",
              "negativeDrop","elev","ice","karst","T_mean","annual_range","prec")

coefficients <- lapply(model, function(x) tryCatch(x %>% 
                                                     pluck("median.coefs") %>% 
                                                     matrix(ncol=length(.)/3) %>% 
                                                     t() %>% 
                                                     data.frame %>% 
                                                     set_rownames(var_names) %>% 
                                                     rowSums() %>% 
                                                     data.frame() %>% 
                                                     set_colnames("effect_size"),
                                                   error=function(e) NA)) 

##### Get spline values

splines<- lapply(model, function(mod) tryCatch({
index<- which(!sapply(mod$gdms,is.null, simplify=TRUE))[1] #retrieve index for baseline model 
gdms <- mod$gdms
PSAMPLE <- 200
preddata <- rep(0,times=PSAMPLE)
preds <- length(gdms[[index]]$predictors)
pred.names <- gdms[[index]]$predictors
predmax <- 0
splineindex <- 1
numsplines <- 3

curves <- list()
for(i in 1:preds){  
  
  predplotdat <- lapply(gdms,function(x).C( "GetPredictorPlotData", 
                                            pdata = as.double(preddata),
                                            as.integer(PSAMPLE),
                                            as.double(x$coefficients[splineindex:(splineindex+numsplines-1)]),
                                            as.double(x$knots[splineindex:(splineindex+numsplines-1)]),
                                            as.integer(numsplines),
                                            PACKAGE = "gdm"))
  
  ## generating posterior stats from runs.
  quant.preds <- apply(plyr::ldply(predplotdat, function(x) c(x$pdata)),2,
                       function(x)quantile(x,c(.05,.5,.95),na.rm=T))
  varNam <- gdms[[index]]$predictors[i]
  env_grad <- matrix(seq(from=gdms[[index]]$knots[splineindex], to=gdms[[index]]$knots[(splineindex+numsplines-1)], length=PSAMPLE),ncol=1)
  colnames(env_grad) <- pred.names[i]
  #I added some modifications to this piece of code so we can make a table ready for ggplot
  curves[[i]] <- data.frame(env_grad,t(quant.preds))
  curves[[i]]$variable <- rep(colnames(curves[[i]])[1],PSAMPLE)
  colnames(curves[[i]])[1]<- "value"
  splineindex <- splineindex + numsplines
  rownames(curves[[i]])<-NULL 
}
res<- curves %>%  bind_rows()},
error = function(e) res<-NA)
)

r2_all_models<- lapply(model, function(x) sapply(x$gdms, function(i) i$explained))
r2_median <-lapply(r2_all_models, function(x) median(unlist(x)))

results <- list(curves = splines, #save models 
                coefficient = coefficients, #waldtest results 
                r2 = r2_median # r2 results
                ) 

results %>%  transpose() %>% pluck("Btotal","curves") %>% 
  ggplot(aes(x=value,y=X50.))+
  geom_line() + 
  geom_ribbon(aes(ymin=X5. ,ymax= X95.), alpha=.3,colour="gray")+
  facet_wrap(~variable, scales="free_x")
  
saveRDS(results, file=glue::glue("Models/Null_BBGDM_{jj}.rds"))

site.pos<-matrix(NA,nrow=length(alpha.FD@HVList), ncol=3)

for (i in 1:length(alpha.FD@HVList)) {
  
site.pos[i,] <- get_centroid(alpha.FD@HVList[[i]])   

}

rownames(site.pos) <- sapply(alpha.FD@HVList, function(x) x@Name)

gamma_space<- hypervolume_gaussian(trait[colSums(comm) > 0,],
                              weight= occupancy)



occupancy<- colSums(comm)[colSums(comm) > 0 ]

gamma_space@RandomPoints %>% 
  data.frame() %>% 
  ggplot(aes(Pco1,Pco2)) + 
  geom_density2d_filled() + 
  geom_point(data=data.frame(trait), aes(Pco1,Pco2), col="red")+
geom_point(data=data.frame(site.pos), aes(X1,X2), col="white")
