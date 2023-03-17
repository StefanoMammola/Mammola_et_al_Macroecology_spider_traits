## ------------------------------------------------------------------------
## '
## ------------------------------------------------------------------------

# Functional trait dimension of spider community

## Change scale progressively

## BBGDM - compare troglophiles and troglobionts

## For each scale: null model where we randomly shuffle species and calculated the observed FD vs expected FD... 
## Pool of species for the null modelling is the cell in the upper level in which each cave occur.
## Model SES vs scale.

# Competition should act at lower scales, environmental filtering at upper scales.

# Beta Taxonomy 

# Beta Dispersion

# Beta Evenness

## ------------------------------------------------------------------------
# 'R script to reproduce the analyses'
## ------------------------------------------------------------------------

# Analysis performed with R (v. R 4.1.0) and R studio (v. 1.4.1103)

# Setting working directory -----------------------------------------------

setwd("/Users/stefanomammola/Desktop/PAPERS IN CORSO/Mammola_et al_spider_traits_macroecology_europe") #change with your directory

# Loading R packages ------------------------------------------------------

library("Amelia") 
library("arakno")
library("BAT")
library("ggplot2")
#library("ggpointdensity")
#library("ggnewscale")
#library("gridExtra")
library("hypervolume")
#library("StatMatch")
library("psych")
library("summarytools")
library("tidyverse")

# Custom function for kernel.beta to have parallel --------------------

#################################################################
# Loading the data ----------------------------------------------
#################################################################

### loading the community data ### 
comm   <- read.csv("Community_composition.csv", sep="\t", dec=",", as.is = FALSE) 
site   <- read.csv("Cave_description.csv", sep="\t", dec=",", as.is = FALSE) 

### loading the trait data ###
trait  <- read.csv("Database_Mammola_et_al_2021_Figshare.csv", sep="\t", dec=",", as.is = FALSE) 

### get the World Spider Catalog ###
arakno::wsc() 

#################################################################
# Prepare data and check taxonomy -------------------------------
#################################################################

### check taxonomy in the comm matrix ###

# re-define Ecological_classification (adaptation) in the comm database as in the trait database
colnames(comm)[7] <- "Ecological_classification"
comm$Ecological_classification <- ifelse(comm$Ecological_classification == 1, "Troglobiont", "Troglophile")

# remove "sp."
comm$Species <- ifelse(comm$Species == "sp.", "", as.character(comm$Species))

# select taxa
tax_comm  <- c(paste(comm$Genus,comm$Species,sep = " "))

# changing problematic taxa (glitch in arakno)
tax_comm[16]  <- "Tegenaria bosnica"     #arakno::checknames("Tegenaria animata") gives an error, manually replacing the synonym
tax_comm[183] <- "Porrhomma rosenhaueri" #arakno::checknames("Porrhomma myops") gives an error, manually replacing the synonym

#Updating taxomy
tax_comm <- data.frame(arakno::checknames(tax_comm, full = TRUE))

tax_comm_checked <- c()

for(i in 1:nrow(tax_comm)){  
  
  if(tax_comm$Note[i] == "Nomenclature change" | tax_comm$Note[i] == "Ok") { #if name change or OK use the right name
    
    tax_comm_checked <- c(tax_comm_checked, tax_comm$Best.match[i])
    
  } else { #for species of unknown identity (sp.), use traits for a congeneric species with similar level of adaptation
    
    sp_adaptation <- comm$Ecological_classification[i]
    sample_sp <- trait[trait$Genus == as.character(comm$Genus[i]) & trait$Ecological_classification == sp_adaptation,]
    sampled_sp <- ifelse(nrow(sample_sp) > 0, as.character(sample(sample_sp$Genus_species)[1]), tax_comm$Best.match[i]) #if match, randomly sample a species to assign its traits; else no match
    
    tax_comm_checked <- c(tax_comm_checked, sampled_sp)
    
  }
  
}

#double check
(tax_comm_checked <- data.frame(tax_comm, checked = tax_comm_checked)) 

## save the new community data matrix and transpose it (site x species matrix)

#selecting only sites
comm2 <- comm[9:ncol(comm)] 

#removing row & column names
rownames(comm2) <- colnames(comm2) <- NULL

#transpose & rename rows and colums
comm2 <- t(comm2) ; dim(comm2) #475 caves, 325 species

rownames(comm2) <- colnames(comm)[9:ncol(comm)]
colnames(comm2) <- tax_comm_checked$checked

#Summing up rows with duplicated name after taxponomic correction colums
comm2 <- t(rowsum(t(comm2), group = colnames(comm2), na.rm = T))
comm2 <- data.frame(comm2) ; colnames(comm2) <- gsub('\\.', ' ', colnames(comm2))

### Check taxonomy in the trait matrix ###
trait$Genus_species <- as.character(trait$Genus_species)

#renaming "sp."
trait$Genus_species[383] <- "Troglohyphantes vignai"    #paper in prep.
trait$Genus_species[384] <- "Troglohyphantes vignai"    #paper in prep.
trait$Genus_species[385] <- "Troglohyphantes angelicus" #paper in prep.

# changing problematic taxa (glitch in arakno)
trait$Genus_species[25]  <- "Tegenaria bosnica" #arakno::checknames("Tegenaria animata") gives an error, manually replacing the synonym
trait$Genus_species[53]  <- "Dysdera pavani"

#Check & update taxonomy
tax_trait <- data.frame(arakno::checknames(trait$Genus_species, full = TRUE))
trait$Genus_species <- tax_trait$Best.match

# Wrap up -----------------------------------------------------------------

# sort alphabetically
comm2  <- comm2[, order(colnames(comm2))]

# sort alphabetically
trait2 <- trait[-c(25,383:385),] #dropping duplicated levels
rownames(trait2) <- trait2$Genus_species

trait2  <- trait2[order(rownames(trait2)),]

# clean the work space
rm(trait, tax_comm, sample_sp, sp_adaptation, i, tax_trait, tax_comm_checked, wscdata, sampled_sp)

#################################################################
# Prepare the traits matrix -------------------------------------
#################################################################

# Only selecting species traits for species in comm -----------------------

comm2   <- comm2[,colnames(comm2) %in% rownames(trait2)]
trait2  <- trait2[rownames(trait2) %in% colnames(comm2),]

all(rownames(trait2) == colnames(comm2)) # should be TRUE


# Select traits that will be used -----------------------------------------

trait_m <- trait2 %>%
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
    Prosoma_shape) 

trait_m <- trait_m %>%
  mutate_at(c("AME_type","Eyeless","Eyes_regression",
              "Capture_web", "Sensing_web", "No_web","Tube_web" ,"Sheet_web",            
              "Space_web", "Orb_web", "Ambush_hunter","Active_hunter","Food_specialist"), as.factor) %>%
  mutate_at("Pigment", ordered, levels = rev(c(
    "Depigmented", "Partly", "Variable", "Fully"
  ))) 

# Data exploration and preparation ----------------------------------------

# General summary
summarytools::dfSummary(trait_m, style='grid')

# Checking continous variables
continous <- c(15:22) #Id continous variables

par(mfrow= c(3,3))
for (i in continous)
  dotchart(trait_m[,i], main= colnames(trait_m)[i]) ; rm(i)
 
# Missing data 
Amelia::missmap(trait_m)

# Standardize continuous traits
trait_m <- BAT::standard(trait = trait_m, method = "standard", convert = continous)

# Collinearity
psych::pairs.panels(trait_m[,colnames(trait_m[continous])]) 

# Estimating Gower distance  --------------------------------------

# Creating groups for traits to obtain optimization
groups_traits <-
  c(rep(1, 3), #Pigment, Eyeless, Eyes_regression
    rep(2, 10), # web and hunting (hunting strategy)
    rep(1,6), # eyes, Femur elongation
    rep(3,3) # body
  )

length(groups_traits) == ncol(trait_m) #check if groups have same length(). Should be TRUE

# # DO NOT RUN (takes several minutes)
 fdist <- gawdis::gawdis(data.frame(trait_m) %>%
                           mutate_at(c("Eyeless","Eyes_regression",
                                       "Capture_web", "Sensing_web", "No_web","Tube_web" ,"Sheet_web",            
                                       "Space_web", "Orb_web", "Ambush_hunter","Active_hunter","Food_specialist"), as.numeric),  #binary traits need to be numeric
                  groups = groups_traits, # grouping variable defined above
                  w.type = "optimized",
                  opti.maxiter = 300,
                  ord="podani",
                  groups.weight = TRUE)
# saveRDS(fdist,"functional_distance.rds") #storing the data

fdist <- readRDS("functional_distance.rds") ; rm(groups_traits)

## Converting the trait space as a Principal Coordinate Analysis
scores <- data.frame(stats::cmdscale(fdist, k = 3)) ; colnames(scores) <- c("PCo1", "PCo2", "PCo3")
trait_m2 <- scores ; rownames(trait_m2) <- trait2$Genus_species

## Plot
centroid <- scores[,1:2] %>%
  add_column(family = trait2$Family) %>%
  group_by(family) %>%
  summarise(cen.1 = mean(PCo1), cen.2 = mean(PCo2))

#with points
theme_set(theme_bw())
(plot_center <- ggplot(scores[,1:2], aes(PCo1, PCo2)) +
    scale_fill_gradient(low = "blue", high = "orange") +
    ggpointdensity::geom_pointdensity(size = 3, alpha = .6) +
    geom_point(data = centroid,
               aes(cen.1, cen.2),
               colour = "black",
               shape = 19) +
    theme(legend.position = "none") +
    ggrepel::geom_label_repel(data = centroid, aes(x = cen.1, y = cen.2, label =
                                                     family)) +
    scale_colour_viridis_b())

#with density
myCol <- viridis::viridis(n = 7, option = "B")

(plot_density <-  ggplot(scores[,1:2], aes(PCo1, PCo2)) +
    stat_density_2d(
      aes(fill = ..level..),
      contour = T,
      adjust = 1.5,
      geom = "polygon",
      colour = NA,
      alpha = .5
    ) +
    geom_density_2d(colour = "black", adjust = 1.5,alpha = .1) +
    scale_fill_gradientn(colours = rev(myCol)) +
    ggnewscale::new_scale_fill() +
    ylim(-.5, .5) +
    xlim(-.6, .4) +
    theme(legend.position = "none") +
    
    geom_point(
      data = centroid,
      aes(x = cen.1, y = cen.2),
      fill = "white",
      size = 3,
      shape = 21,
      colour = "black",
      stroke = 2
    ) +
    geom_point(
      data = centroid,
      aes(x = cen.1, y = cen.2),
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
      aes(cen.1, cen.2, label = family),
      segment.linetype = 2,
      box.padding = 0.5,
      max.overlaps = Inf,
      segment.size = .2
    )
)

#################################################################
# Calculating functional diversity ------------------------------
#################################################################

# Selecting community with > 2 species
comm3 <- comm2[rowSums(comm2) > 2,]

# Matching wth Site matrix
site2 <- site[site$ID %in%  rownames(comm3),]

# Generating full hypervolume --------------------------------------------
hv_tot <- BAT::kernel.build(comm = rep(1, nrow(trait_m2)), trait = trait_m2,
                          method = "gaussian", axes = 0, abund = FALSE, cores = 1)

# Visualize
plot(hv_tot, 
     show.data=TRUE,
     show.random=TRUE,
     show.centroid=FALSE,
     col="black", 
     num.points.max.random = 3000,
      # plot.function.additional=function(i,j){
      #  # points( x= data[data$Adapt=="TB", i], y=data[data$Adapt =="TB", j],cex=1,pch=17)
      #  # points( x= data[data$Adapt=="TF", i], y=data[data$Adapt =="TF", j],cex=1,pch=19)
      #  
      #  points( x= data[, i], y=data[, j],cex=1.2,pch=16,col="black")
      #  
      #  
      #  points( x= data[data$Family=="Agelenidae", i], y=data[data$Family =="Agelenidae", j],cex=1,pch=16,col="brown")
      #  points( x= data[data$Family=="Linyphiidae", i], y=data[data$Family =="Linyphiidae", j],cex=1,pch=16,col="red")
      #  #points( x= data[data$Family=="Dysderidae", i], y=data[data$Family =="Dysderidae", j],cex=1,pch=16,col="orange")
      #  
      #  points( x= data[data$Family=="Pholcidae", i], y=data[data$Family =="Pholcidae", j],cex=1,pch=16,col="orange")
      #  
      #  points( x= data[data$Family=="Leptonetidae", i], y=data[data$Family =="Leptonetidae", j],cex=1,pch=16,col="blue")
      #  points( x= data[data$Family=="Nesticidae", i], y=data[data$Family =="Nesticidae", j],cex=1,pch=16,col="turquoise")
      #  
      #  
      #  points( x= data[data$Family=="Tetragnathidae", i], y=data[data$Family =="Tetragnathidae", j],cex=1,pch=16,col="grey")
      #  
      #  points( x= data[data$Family=="Pimoidae", i], y=data[data$Family =="Pimoidae", j],cex=1,pch=16,col="yellow")
      #  
      #  points( x= data[data$Family=="Symphytognathidae", i], y=data[data$Family =="Symphytognathidae", j],cex=1,pch=16,col="purple")
      #  
      # }
)

# Generating hypervolumes -------------------------------------------------

# DO NOT RUN (takes several minutes)
# hv <- BAT::kernel.build(comm = comm3, trait = trait_m2, 
#                         method = "gaussian", axes = 0, abund = FALSE, cores = 10)
# 
# saveRDS(hv,"hypervolumes.rds") #storing the data
hv <- readRDS("hypervolumes.rds")

#Calculating the three dimension of FD (richness, divergence, regularity)
FDrich <- BAT::kernel.alpha(hv)
FDdisp <- BAT::kernel.dispersion(hv, func = "dissimilarity", frac = 0.1)
FDeve  <- BAT::kernel.evenness(hv)

results_alpha <- data.frame(ID = site2$ID, FDrich, FDdisp, FDeve) ; rm(FDrich, FDdisp, FDeve)

par(mfrow = c(1,3), mar = rep(2,4))
hist(results_alpha$FDrich, main = "FD richness")
hist(results_alpha$FDdisp, main = "FD divergence")
hist(results_alpha$FDeve,  main = "FD regularity")

#What is the correlation among the three metrics?
psych::pairs.panels(results_alpha[,2:4]) 

#######################################
# Calculating functional beta diversity
#######################################

#FDbeta <- kernel.beta(comm = hv, func = "jaccard", comp = FALSE)

#Custum way to have parallel

library("doSNOW")
library("future")
appendList <- function (x, val) {
  stopifnot(is.list(x), is.list(val))
  xnames <- names(x)
  for (v in names(val)) {
    x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]]))
      appendList(x[[v]], val[[v]])
    else c(x[[v]], val[[v]])
  }
  x
}

print(as.character(Sys.time()), stdout())
ntasks <- length(hv@HVList) # Count hypervolumes that succeed in the estimation

cl <- makeSOCKcluster(10) # create the parallel structure for 10 cores
registerDoSNOW(cl)  # initialize the parallel

print(as.character(Sys.time()), stdout())

output_beta <-  foreach(i = 1:ntasks,
                        .combine = appendList) %:%
  foreach(j = i:ntasks,
          .combine = appendList) %dopar% {
      
            .libPaths(c("/projappl/project_2004675/project_rpackages", .libPaths())) #set again the library path
            if(!require("hypervolume")) {install.packages("hypervolume")} #load package in the parallel environment
            
            hyperSet <- hypervolume::hypervolume_set(hv@HVList[[i]], hv@HVList[[j]],
                                                     check.memory = FALSE, verbose = FALSE, num.points.max = 1000)
            union   <- hyperSet[[4]]@Volume
            unique1 <- hyperSet[[5]]@Volume
            unique2 <- hyperSet[[6]]@Volume
            
            union <- 2 * union - unique1 - unique2
            Btotal <- (unique1 + unique2)/union
            Brepl <- 2 * min(unique1, unique2)/union
            Brich <- abs(unique1 - unique2)/union
            output<- list(Btotal = Btotal, Brepl = Brepl, Brich = Brich)
          }
stopCluster(cl)

TRpa<-list()
TRpa$Btotal <- matrix(NA,ntasks,ntasks)
TRpa$Brepl <- matrix(NA,ntasks,ntasks)
TRpa$Brich <- matrix(NA,ntasks,ntasks)

TRpa$Btotal[lower.tri(TRpa$Btotal,diag=TRUE)] <-round(output_beta$Btotal,3)
TRpa$Brepl[lower.tri(TRpa$Btotal,diag=TRUE)] <- round(output_beta$Brich,3)
TRpa$Brich[lower.tri(TRpa$Btotal,diag=TRUE)] <- round(output_beta$Brepl,3)
TRpa<-lapply(TRpa, as.dist)



hist(TRpa$Btotal)
hist(TRpa$Brepl)
hist(TRpa$Brich)

library("fossil")

geo_dist <- fossil::earth.dist(data.frame(site2$decimalLongitude,site2$decimalLongitude), dist = TRUE)

plot(as.vector(geo_dist),as.vector(TRpa$Btotal))
plot(as.vector(geo_dist),as.vector(TRpa$Brepl))
plot(as.vector(geo_dist),as.vector(TRpa$Brich))

vegan::mantel(geo_dist, TRpa$Btota, method="spearman", parallel = 10)
vegan::mantel(geo_dist, TRpa$Brepl, method="spearman", parallel = 10)
vegan::mantel(geo_dist, TRpa$Brich, method="spearman", parallel = 10)

