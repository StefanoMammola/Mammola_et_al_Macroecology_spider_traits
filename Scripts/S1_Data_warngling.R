


# Set working directory to the main folder where the data are -----------------------------------------------------
dir<- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(dir)



# Load libraries and data------------------------------------------------------------------------------------------

#packages
library("Amelia") #Checking missing data
library("arakno") #Correcting taxonomies
library("BAT") #Biodiversity assessment tools
library("GGally") #Plot variable pairs in ggplots
library("tidyverse") #Data wrangling

#community data 
comm   <- read.csv("Community_composition.csv", sep="\t", dec=",", as.is = FALSE) 

#trait data
trait  <- read.csv("Database_Mammola_et_al_2021_Figshare.csv", sep="\t", dec=",", as.is = FALSE) 

#spatial data
site   <- read.csv("Cave_description.csv", sep="\t", dec=",", as.is = FALSE) 

#World Spider Catalog 
arakno::wsc() 