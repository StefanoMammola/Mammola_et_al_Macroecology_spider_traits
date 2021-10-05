
# Set working directory to the main folder where the data are -----------------------------------------------------
dir<- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(dir)


# Set global variables --------------------------------------------------------------------------------------------

input_dir<- "Data/" #Folder where the data are located
object_dir<-"objects/" #Sub-folder where data are located
output_dir<- "objects/" #Folder where to save the data
raster_dir <- "rasters/" # Folder where the raster data are stored
results_dir <- "Results/" #Folder where to save figures, tables, etc 
table_dir<-"tables/" # Folder where the excel files are saved

run_fdist <- FALSE  #should we calculate functional dissimilarity? This is computationally expensive

# Load libraries and data------------------------------------------------------------------------------------------

#packages
library("Amelia") #Checking missing data
library("arakno") #Correcting taxonomies
library("BAT") #Biodiversity assessment tools
library("GGally") #Plot variable pairs in ggplot
library("raster") # Used for raster operations 
library("tidyverse") #Data wrangling

load(paste0(input_dir,object_dir,"Parsed_data.R"))

#spatial data
site   <- read.csv(paste0(input_dir,table_dir, "Cave_description.csv"), sep="\t", dec=",", as.is = FALSE) 

# Functional similarity between species ---------------------------------------------------------------------------

#Select traits to be included 

trait_m <- trait_parsed %>%
  select(
    Pigment,
    Eyeless,
    Eyes_regression,
    ends_with(c("_web", "_hunter")),
    Food_specialist,
    AME_type,
    AME,
    ALE,
    PLE, 
    PME,
    Femur_elongation,
    Sexual_Size_Dimorphism,
    Body_length_avg,
    Femur_elongation,
    Prosoma_shape) %>%
  mutate_at(c("AME_type","Eyeless","Eyes_regression",
              "Capture_web", "Sensing_web", "No_web","Tube_web" ,"Sheet_web",            
              "Space_web", "Orb_web", "Ambush_hunter","Active_hunter","Food_specialist"), as.factor) %>%
  mutate_at("Pigment", ordered, levels = rev(c(
    "Depigmented", "Partly", "Variable", "Fully"
  ))) 

# General exploration and preparation ----------------------------------------

# Checking continuous variables
continuous <- c(15:22) #Id continuous variables

par(mfrow= c(3,3))
for (i in continuous)
  dotchart(trait_m[,i], main= colnames(trait_m)[i]) ; rm(i)

# Missing data 
trait_m %>% Amelia::missmap()

# Standardize continuous traits
trait_m <- BAT::standard(trait = trait_m, method = "standard", convert = continous)

# Collinearity
theme_set(theme_classic())
ggpairs(trait_m, columns = continuous, ggplot2::aes(colour=trait_parsed$Ecological_classification))  

# Estimating Gower distance  --------------------------------------

# Creating groups for traits to obtain optimization
groups_traits <-
  c(rep(1, 3), #Pigment, Eyeless, Eyes_regression
    rep(2, 10), # web and hunting (hunting strategy)
    rep(1,6), # eyes, Femur elongation
    rep(3,3) # body
  )

length(groups_traits) == ncol(trait_m) #check if groups have same length(). Should be TRUE

if (run_fdist == TRUE){
fdist <- gawdis::gawdis(data.frame(trait_m) %>%
                          mutate_at(c("Eyeless","Eyes_regression",
                                      "Capture_web", "Sensing_web", "No_web","Tube_web" ,"Sheet_web",            
                                      "Space_web", "Orb_web", "Ambush_hunter","Active_hunter","Food_specialist"), as.numeric),  #binary traits need to be numeric
                        groups = groups_traits, # grouping variable defined above
                        w.type = "optimized",
                        opti.maxiter = 300,
                        groups.weight = TRUE)

saveRDS(f_dist,paste0(input_dir,output_dir,"functional_distance.rds")) #storing the data
}

# Principal Coordinate Analysis -----------------------------------------------------------------------------------

f_dist <- readRDS(paste0(input_dir,output_dir,"functional_distance.rds")) #(re)loading the data

## Extracting trait axis from Principal component Analysis

trait_axis<- f_dist %>% 
  ape::pcoa() %>% #run principal component analysis using package ape
  pluck("vectors") %>% #select list object using purrr package
  as_tibble()  %>%  #convert to data.frame
  select(1:4) %>% #select columns 
  rename_with(~paste0("Pco",1:4)) %>%  # rename all columns
  mutate(species_names = trait_parsed$Genus_species) %>% #create column with species rownames
  column_to_rownames("species_names") %>%  #select the column species_names to be the rownames of the table
 as_tibble()

# Extracting centroid per family

centroid <- 
  trait_axis %>% 
  add_column(family = trait_parsed$Family) %>%
  group_by(family) %>%
 # summarise(Centroid_x = mean(PCo1), Centroid_y = mean(PCo2))
 summarise(across(starts_with("Pco"), ~ mean(.x, na.rm = TRUE)))

#Plot species coordinates with family centroids 
myCol <- viridis::viridis(n = 7, option = "B")

(plot_density <-  trait_axis %>%  
    ggplot(aes(Pco1, Pco2)) +
    stat_density_2d(
      aes(fill = ..level..),
      contour = T,
      adjust = 1.5,
      geom = "polygon",
      colour = NA,
      alpha = .5,
      show.legend = FALSE
    ) +
    geom_density_2d(colour = "black", adjust = 1.5,alpha = .1) +
    scale_fill_gradientn(colours = rev(myCol)) +
    ggnewscale::new_scale_fill() +
    ylim(-.5, .5) +
    xlim(-.6, .4) +
    theme(legend.position = "none") +
    geom_point(
      data = centroid,
      fill = "white",
      size = 3,
      shape = 21,
      colour = "black",
      stroke = 2
    ) +
    geom_point(
      data = centroid,
      colour = "black",
      size = 1
    ) +
    theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),
      text = element_blank(),
      axis.ticks = element_blank()
    ) +
    ggrepel::geom_text_repel(
      data = subset(centroid),
      aes(label = family),
      segment.linetype = 2,
      box.padding = 0.5,
      max.overlaps = Inf,
      segment.size = .2
    ) +
    theme_bw()
)
ggsave(plot_density, filename = paste0(results_dir,"functional_pcoa.pdf"), device=cairo_pdf)


# Extracting the broad scale variables  ---------------------------------------------------------------------------

 
predictors <- raster::stack(list.files(path = paste0(input_dir, raster_dir), pattern="*.tif", full.names = TRUE))

raster::projection(predictors) <-"+proj=longlat +datum=WGS84 +ellps=WGS84"
names(predictors) <- c("elev",
                       "ice",
                       "karst",
                       "T_mean",
                       "annual_range",
                       "prec")

predictors_regional<- extract(predictors, site[,4:5]) 

predictors <-
  site %>%
  dplyr::select(
    decimalLongitude,
    decimalLatitude,
    entranceNumber,
    entranceSize,
    development,
    ID,
    positiveDrop,
    negativeDrop
  ) %>% # select columns 
  add_column(data.frame(predictors_regional)) %>% # add regional predictors
  mutate_at(vars(entranceSize, development), log1p) %>% #log transform variables 
  column_to_rownames("ID")

ggpairs(predictors) # check collinearity

save(comm_parsed, trait_axis, predictors, paste0(input_dir, object_dir,"BBGDM_input.R"))
