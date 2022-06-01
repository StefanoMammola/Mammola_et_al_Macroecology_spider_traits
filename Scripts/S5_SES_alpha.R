
#Null model  script
# rsync -r caioroza@puhti.csc.fi:/scratch/project_2005062/Mammola_et_al_Macroecology_spider_traits "/Users/graco-roza/OneDrive - University of Helsinki/Ongoing manuscripts/Mammola_et_al_Macroecology_spider_traits" caioroza@puhti.csc.fi:/scratch/project_2005062/
# cd /scratch/project_2005062/Mammola_et_al_Macroecology_spider_traits/
###### Set up of the script ############################################################################################

.libPaths(c("/projappl/project_2005062/project_rpackages_4.0.5", .libPaths()))
# Packages .............................................................................................................
if(!require("pacman")) {install.packages("pacman")}
pacman::p_load("tidyverse","hypervolume","doSNOW", "future","BAT","magrittr")
#....................................................................................................................... 

#Program cores for the batch ...........................................................................................
options(future.availableCores.methods = "Slurm")

#....................................................................................................................... 

#directories to save and load files ....................................................................................
input_dir<- "Alpha/" #Folder where the data are located
output_dir<- "F_metrics/" #Folder where to save the data

#....................................................................................................................... 

names <-  input_dir %>%
  list.files %>%
  tools::file_path_sans_ext()

names.file<-  input_dir %>%
  list.files(full.names = TRUE) 

res<-list()
for ( i in 1:length(names.file)){
print(i)
  res[[i]] <- names.file[i] %>% 
  readRDS %>%
         get_volume
}

names(res) <- names

saveRDS(res, paste0(output_dir,"SES_FR.rds"))