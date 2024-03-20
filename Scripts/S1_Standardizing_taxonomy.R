## ------------------------------------------------------------------------
## ' Functional convergence underground? The scale-dependency of community assembly processes in European cave spiders
## ------------------------------------------------------------------------

# Set working directory to the main folder where the data are -----------------------------------------------------
dir<- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(dir)

# Load libraries and data------------------------------------------------------------------------------------------

#packages
library("Amelia") #Checking missing data
library("arakno") #Correcting taxonomies
library("BAT") #Biodiversity assessment tools
library("GGally") #Plot variable pairs in ggplot
library("tidyverse") #Data wrangling

input_dir<- "Data/"
table_dir<-"tables/"
output_dir<- "objects/"

#community data 
comm   <- read.csv(paste0(input_dir,table_dir,"Community_composition.csv"), sep="\t", dec=",", as.is = FALSE) 

#trait data
trait  <- read.csv(paste0(input_dir,table_dir,"Database_Mammola_et_al_2021_Figshare.csv"), sep="\t", dec=",", as.is = FALSE) 

#World Spider Catalog 
arakno::wsc() 


# Check taxonomy  -------------------------------------------------------------------------------------------------

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

#selecting only sites
comm_subset <- comm[9:ncol(comm)] 

#removing row & column names
rownames(comm_subset) <- colnames(comm_subset) <- NULL

#transpose & rename rows and columns
comm_subset <- t(comm_subset)
dim(comm_subset) #475 caves, 325 species

rownames(comm_subset) <- colnames(comm)[9:ncol(comm)]
colnames(comm_subset) <- tax_comm_checked$checked

#Summing up rows with duplicated name after taxonomic correction columns
comm_parsed <- t(rowsum(t(comm_subset), group = colnames(comm_subset), na.rm = T))
comm_parsed <- data.frame(comm_parsed) ; colnames(comm_parsed) <- gsub('\\.', ' ', colnames(comm_parsed))
comm_parsed  <- comm_parsed[, order(colnames(comm_parsed))]

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
trait_parsed <- trait[-c(25,383:385),] #dropping duplicated levels
rownames(trait_parsed) <- trait_parsed$Genus_species #assign row names
trait_parsed  <- trait_parsed[order(rownames(trait_parsed)),] #order row names 


comm_parsed   <- comm_parsed[,colnames(comm_parsed) %in% rownames(trait_parsed)]
trait_parsed  <- trait_parsed[rownames(trait_parsed) %in% colnames(comm_parsed),]

all(rownames(trait_parsed) == colnames(comm_parsed)) # should be TRUE

save(comm_parsed, trait_parsed, file=paste0(input_dir,output_dir,"Parsed_data.R"))