
#Null model  script
# rsync -r caioroza@puhti.csc.fi:/scratch/project_2005062/Mammola_et_al_Macroecology_spider_traits "/Users/graco-roza/OneDrive - University of Helsinki/Ongoing manuscripts/Mammola_et_al_Macroecology_spider_traits" caioroza@puhti.csc.fi:/scratch/project_2005062/
# cd /scratch/project_2005062/Mammola_et_al_Macroecology_spider_traits/
###### Set up of the script ############################################################################################

.libPaths(c("/projappl/project_2005062/project_rpackages_4.0.5", .libPaths()))
# Packages .............................................................................................................
if(!require("pacman")) {install.packages("pacman")}
pacman::p_load("tidyverse","hypervolume","doSNOW", "future","gdm", "surveillance","plotrix")
#....................................................................................................................... 

#Program cores for the batch ...........................................................................................
options(future.availableCores.methods = "Slurm")
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
comm_parsed<- comm_parsed[1:10,]
####### GET FOCAL DISSIMILARITY METRIC & ADD COLUMN OF SITE NAMES 

cores<- 40
cl <- makeSOCKcluster(cores)
registerDoSNOW(cl)
results$Trial_1$Btotal$model$gdms[[1]]$
results<- foreach(i = 1, .combine = combine_custom) %dopar% {
  
  .libPaths(c("/projappl/project_2005062/project_rpackages_4.0.5", .libPaths()))
  if(!require("pacman")) {install.packages("pacman")}
  pacman::p_load("tidyverse","hypervolume","doSNOW", "future","gdm")
  
  if (i == 1) {
    trait_shuffled <- trait_axis
  } else {
    trait_shuffled <- trait_axis[sample(nrow(trait_axis)), ]
  }
  rownames(trait_shuffled) <-  rownames(trait_axis)
  
  fd <- beta_fd(comm = comm_parsed, trait = trait_shuffled)
  
  output<- lapply(fd, function(x) tryCatch(run_bbgdm(x), error=function(e) NA))
  return(output)
}
stopCluster(cl)
names(results) <- paste0("Trial_",1:length(results))

saveRDS(results, file="Models/Null_BBGDM.rds")
