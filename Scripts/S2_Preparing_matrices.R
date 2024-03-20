## ------------------------------------------------------------------------
## ' Functional convergence underground? The scale-dependency of community assembly processes in European cave spiders
## ------------------------------------------------------------------------

# Set working directory to the main folder where the data are -----------------------------------------------------
dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(dir)

#remotes::install_github("jthomasmock/gtExtras")
# Set global variables --------------------------------------------------------------------------------------------

input_dir <- "Data/" #Folder where the data are located
object_dir <- "objects/" #Sub-folder where data are located
output_dir <- "objects/" #Folder where to save the data
raster_dir <- "rasters/" # Folder where the raster data are stored
results_dir <-
  "Results/" #Folder where to save figures, tables, etc
table_dir <- "tables/" # Folder where the excel files are saved

run_fdist <-
  FALSE  #should we calculate functional dissimilarity? This is computationally expensive

# Load libraries and data------------------------------------------------------------------------------------------

#packages

pacman::p_load(
  "Amelia" #Checking missing data
  ,
  "arakno" #Correcting taxonomies
  ,
  "BAT"    #Biodiversity assessment tools
  ,
  "GGally" #Plot variable pairs in ggplot
  ,
  "raster" # Used for raster operations
  ,
  "tidyverse" #Data wrangling
  ,
  "gt" #make tables
  ,
  "gtExtras" #format gt tables
  ,
  "colorspace" #colors manipulation
  ,
  "gtsummary"
)



# Load data -------------------------------------------------------------------------------------------------------

load(paste0(input_dir, object_dir, "Parsed_data.R"))

#spatial data
site   <-
  read.csv(
    paste0(input_dir, table_dir, "Cave_description.csv"),
    sep = "\t",
    dec = ",",
    as.is = FALSE
  )

# Functional similarity between species ---------------------------------------------------------------------------

#Select traits to be included

trait_m <- trait_parsed %>%
  dplyr::select(
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
    Prosoma_shape
  ) %>%
  mutate_at(
    c(
      "AME_type",
      "Eyeless",
      "Eyes_regression",
      "Capture_web",
      "Sensing_web",
      "No_web",
      "Tube_web" ,
      "Sheet_web",
      "Space_web",
      "Orb_web",
      "Ambush_hunter",
      "Active_hunter",
      "Food_specialist"
    ),
    as.factor
  ) %>%
  mutate_at("Pigment", ordered, levels = rev(c(
    "Depigmented", "Partly", "Variable", "Fully"
  )))

# General exploration and preparation ----------------------------------------

#Table with trait distribution of species

#We decide to make 2 tables, one for continuous traits and one for categorical
#The reason is that continuous trait are represented by mean,sd, min, max and histograms.
#Categorical traits are represented by frequency of levels, and barplots


# Continuous traits -----------------------------------------------------------------------------------------------


# function for plotting histogram lines
plot_hist <- function(name, df) {
  plot_object <-
    ggplot(data = df,
           aes(x = value)) +
    geom_histogram(colour = NA, fill = "mediumorchid4") +
    theme_void()
  return(plot_object)
}

# function for plotting percentage charts
plot_barchart <-
  function(name, df) {
    plot_object <- df %>%
      ggplot() +
      geom_bar(aes(x = 1, y = name),
               stat = "identity",
               fill = lighten("mediumorchid4", .9)) +
      geom_bar(aes(x = Freq, y = name), stat = "identity", fill = "mediumorchid4") +
      theme_void()
    return(plot_object)
  }

data_continuous <-  trait_m %>%
  select(AME:Prosoma_shape) %>%
  rownames_to_column("species") %>%
  pivot_longer(!species, names_to = "Trait", values_to = "value") %>%
  select(-species) %>%
  group_by(Trait) %>%
  drop_na() %>%
  summarise(
    mean = mean(value),
    sd = sd(value),
    min = min(value),
    max = max(value),
    n = NA,
    Freq = NA,
    data = list(data.frame(value = value))
  ) %>%
  dplyr::mutate(inset = map2(Trait, data, plot_hist))  %>%
  mutate(
    mean_sd = glue::glue(
      "{formatC( round( mean, 2 ), format='f', digits=2 )} ± {formatC( round( sd, 2 ), format='f', digits=2 )}"
    ),
    range = glue::glue(
      "{formatC( round( min, 2 ), format='f', digits=2 )} – {formatC( round( max, 2 ), format='f', digits=2 )}"
    ),
    N_Freq = NA,
    .keep = "unused"
  )

data_categorical <-
  trait_m %>%
  select(Pigment:AME_type) %>%
  mutate(Pigment = factor(
    as.character(Pigment),
    levels = c("Fully", "Variable", "Partly", "Depigmented")
  )) %>%
  rownames_to_column("species") %>%
  #mutate(across(c(Pigment,AME_type),~ factor(as.numeric(.x)))) %>%
  as_tibble() %>%
  pivot_longer(!species, names_to = "Trait", values_to = "value") %>%
  mutate(value = case_when(value == 0 ~ "FALSE", value == 1 ~ "TRUE", TRUE ~ as.character(value))) %>%
  select(-species) %>%
  mutate(
    Trait2 = fct_collapse(
      Trait,
      "Pigmentation" = "Pigment",
      "Absence of eyes" = "Eyeless",
      "Eyes regression" = "Eyes_regression",
      "Hunting strategy" = paste0(c("Capture_", "Sensing_", "No_"), "web"),
      "Guild" = c(paste0(
        c("Tube_", "Sheet_", "Space_", "Orb_"), "web"
      ), "Ambush_hunter", "Active_hunter"),
      "Food specialization" = "Food_specialist",
      "Eye ontology" = "AME_type"
    )
  ) %>%
  janitor::tabyl(Trait2, Trait, value) %>%
  bind_rows(.id = "Type") %>%
  filter(Type != 'NA_') %>%
  pivot_longer(!c(Type, Trait2)) %>%
  filter(value != 0) %>%
  arrange(Trait2, name, Type, value) %>%
  group_by(name) %>%
  mutate(Freq = value / sum(value)) %>%
  filter(Type != "FALSE") %>%
  add_column(
    mean = NA,
    sd = NA,
    min = 0,
    max = 1
  ) %>%
  relocate(Trait2, name, Type, value, Freq) %>%
  mutate(name = case_when(name == "AME_type" ~ Type, name == "Pigment" ~ Type, TRUE ~ name)) %>%
  group_by(name) %>%
  mutate(data = list(data.frame(name = name, Freq = Freq))) %>%
  mutate(inset = map2(name, data, plot_barchart)) %>%
  select(-Type) %>%
  rename(n = value,
         Trait = Trait2,
         Category = name) %>%
  mutate(
    mean_sd = NA,
    range = "Yes/No" ,
    N_Freq = glue::glue(
      "{formatC( as.integer(n), format='f', digits=0)} ({formatC( round(Freq*100,2), format='f', digits=2)}%)"
    ),
    .keep = "unused"
  )



table_data <- 
  list(Continuous = data_continuous,
                   Categorical = data_categorical) %>% 
  bind_rows() %>% 
  relocate(Trait, Category) %>%
  mutate(Category2 = ifelse(is.na(Category), Trait, Category),
         Trait = ifelse(is.na(Category), "", Trait)) %>%
  add_column(ggplot1 = NA) %>%
  mutate(id = row_number()) %>%
  mutate(Category = Category2, .keep = "unused") %>%
  select(-c(mean,sd,min,max)) %>% 
  mutate(Category = fct_recode(as.factor(Category),"Anterior Lateral Eyes (ALE)" = "ALE"
                               ,"Anterior Median Eyes (AME)" = "AME"
                               ,"Average Body length" = "Body_length_avg"
                               ,"Femur elongation" = "Femur_elongation"
                               ,"Posterior Lateral Eyes (PLE)" = "PLE"
                               ,"Posterior Median Eyes (PME)" = "PME"
                               ,"Prosoma shape" = "Prosoma_shape"
                               ,"Sexual size dimorphism" = "Sexual_Size_Dimorphism"
                               ,"Active hunter" = "Active_hunter"
                               ,"Ambush hunter" = "Ambush_hunter"
                               ,"Orb web" = "Orb_web"
                               ,"Sheet web" = "Sheet_web"
                               ,"Space web" = "Space_web"
                               ,"Tube web" = "Tube_web"
                               ,"Absent due to adaptation" = "Absent_Adaptation"
                               ,"Absent due to ontology" = "Absent_Ontology"
                               ,"Capture web" = "Capture_web"
                               ,"No web" = "No_web"
                               ,"Sensing web" = "Sensing_web"
                               ,"Absence of eyes" = "Eyeless"
                               ,"Eye regression" = "Eyes_regression"
                               ,"Specialized diet (stenophagous)" = "Food_specialist"
                               ,"Depigmented" = "Depigmented"
                               ,"Pigmented" = "Fully"
                               ,"Partly pigmented" = "Partly"
                               ,"Variable" = "Variable") 
         
         ) %>% 
    mutate(Category = factor(Category, levels=Category[c(3,7,8,1,2,5,6,4,21,22,17,15,16,25,26,27,24,9,10,11,12,13,14,18,20,19,23)])) %>% 
  arrange(Category) %>% 
  add_column("Gower grouping" =  c( 
    rep("Morphology",3)
   ,rep("Adaptation",14)
   ,rep("Ecology", 10)
  )) %>%   
    #add gower groups to the table 
  mutate(id = row_number()) 

table_data$Trait[c(1:10,nrow(table_data))] <- c("", "  ", "   ", "    ","     ","      ","       ","        ", "         ","          ","           ")

table_1 <-
  table_data%>%
  gt::gt(id= "spiders",
        # rowname_col = "Category",
         groupname_col =  "Trait") %>%
    cols_align(align = "left", column = Category) %>% 
    cols_label(Category = "Functional trait"  
               ,mean_sd = "Mean ± SD"
               ,N_Freq = "N (%)"
               ,ggplot1 = "Distribution"
               ,"Gower grouping" = "Group"
                 ) %>% 
  tab_style(
    style = cell_borders(
      sides = c("top"),
      color = c("white"),
      weight = px(2)
    ),
    locations = cells_body(rows = c(12,13,15:17,19:23,25:26))
  ) %>% 
      gt::text_transform(
        locations = cells_body(columns = c(ggplot1)),
        fn = function(x) {
          map(table_data$inset,
              ggplot_image,
              height = px(20),
              aspect_ratio = 6)
        }
      ) %>%
      gt::cols_hide(c(data, inset, Trait, n, Freq, id)) %>% 
      fmt_missing(columns = 1:9,
                  rows = everything(),
                  missing_text = "") %>%
      opt_align_table_header("left") %>%
    tab_style(
      style = cell_borders(
        sides = c("top"),
        color = c("white"),
        weight = px(2)
      ),
      locations = cells_row_groups(groups = table_data$Trait[c(1:10,27)])
    ) %>% 
    tab_style(
      style = cell_borders(
        sides = c("top"),
        color = c("#d3d3d3"),
        weight = px(2)
      ),
      locations = cells_row_groups(groups = table_data$Trait[c(11:26)])
    ) %>% 
    tab_style(
      style = cell_borders(
        sides = c("bottom"),
        color = c("#D3D3D3"),
        weight = px(2)
      ),
      locations = cells_row_groups(groups = everything())
    ) %>%
    tab_style(
      style = cell_borders(
        sides = c("bottom"),
        color = c("white"),
        weight = px(2)
      ),
      locations = cells_row_groups(groups = unique(table_data$Trait)[11:14])
    ) %>% 
       tab_style(
         style = list("font-variant: small-caps;"),
         locations = cells_body(columns = "Category", rows= 11:26)) %>%
       opt_all_caps(all_caps = TRUE) %>%
       opt_table_font(font = list(google_font("Chivo"),
                                  default_fonts()),
                      weight = 300) %>%
       tab_style(locations = cells_body(columns = "Category", rows= c(1:10,27)),
                 style = list(
                   cell_text(
                     font = google_font(name = 'Chivo'),
                     weight = 400,
                     transform = "uppercase",
                     size = "80%"
                   )
                 )) %>%
       tab_style(
         style = cell_borders(
           sides = "top",
           color = "black",
           weight = px(0)
         ),
         locations = gt::cells_column_labels(columns = gt::everything())
       ) %>%
    tab_style(
      style = list(
        cell_text(align = "left")
      ),
      locations = cells_stub(rows = TRUE)
    ) %>% 
      tab_options(
        table.width = px(900),
        heading.title.font.weight = "bold",
        table.align = "center",
        column_labels.background.color = "white",
        heading.border.bottom.style = "none",
        table.border.top.width = px(3),
        table.border.top.style = "none",
        table.border.bottom.style = "none",
        column_labels.font.weight = "normal",
        column_labels.border.top.style = "none",
        column_labels.border.bottom.width = px(2),
        column_labels.border.bottom.color = "black",
        row_group.border.top.style = "none",
        row_group.border.top.color = "black",
        row_group.border.bottom.style = "none",
        row_group.border.bottom.width = px(0),
        row_group.border.bottom.color = "white",
        stub.border.color = "#D3D3D3",
        stub.border.width = px(2),
        source_notes.font.size = 12,
        source_notes.border.lr.style = "none",
        table.font.size = 16,
        heading.align = "left"
      ) %>% 
  gt::opt_vertical_padding(scale=0.1) %>% 
    cols_align("left",1) %>% 
    opt_css(css = 
    "tbody tr:last-child {\n    border-bottom: 2px solid #ffffff00;\n      }\n\n", 
add = TRUE)

table_1 %>%  gtsave("table_traits.rtf")

modified_table <- readLines("table_traits.html")
groups<-c('Eye ontology','Pigmentation','Guild','Hunting strategy')
for (i in 1:length(groups)){
focal_group<- c('Eye ontology','Pigmentation','Guild','Hunting strategy')[i]
rownumber<- which(str_detect(test, as.character(focal_group)))
modified_table[rownumber] <- glue::glue('<td colspan=\"6\" class=\"gt_group_heading\" style=\"padding-top: 10px; padding-bottom: 10px; border-top-width: 2px; border-top-style: solid; border-top-color: #d3d3d3; border-bottom-width: 2px; border-bottom-style: solid; border-bottom-color: white;\">{focal_group}</td>')
}

writeLines(modified_table, "table_traits.html")

#End of table making ------------------------------------------------------------------------

    # Checking continuous variables
    continuous <- c(15:22) #Id continuous variables
    
    par(mfrow = c(3, 3))
    for (i in continuous) {
      dotchart(trait_m[, i], main = colnames(trait_m)[i])
      rm(i)
    }
    
    # Missing data
    trait_m %>% Amelia::missmap()
    
    # Standardize continuous traits
    trait_m <-
      BAT::standard(trait = trait_m,
                    method = "standard",
                    convert = continuous)
    
    # Collinearity
    theme_set(theme_classic())
    ggpairs(trait_m,
            columns = continuous,
            ggplot2::aes(colour = trait_parsed$Ecological_classification))
    names(trait_m)
    # Estimating Gower distance  --------------------------------------
    
    # Creating groups for traits to obtain optimization
    groups_traits <-
      c(rep(3, 3), #Pigment, Eyeless, Eyes_regression (adaptation)
        rep(2, 10), # web and hunting (hunting strategy)
        rep(3, 6), # eyes, Femur elongation (adaptation)
        rep(1, 3)) # body (morphology and morphometry))
        
        length(groups_traits) == ncol(trait_m) #check if groups have same length(). Should be TRUE
        
        if (run_fdist == TRUE) {
          set.seed(123)
          fdist <- gawdis::gawdis(
            data.frame(trait_m) %>%
              mutate_at(
                c(
                  "Eyeless",
                  "Eyes_regression",
                  "Capture_web",
                  "Sensing_web",
                  "No_web",
                  "Tube_web" ,
                  "Sheet_web",
                  "Space_web",
                  "Orb_web",
                  "Ambush_hunter",
                  "Active_hunter",
                  "Food_specialist"
                ),
                as.numeric
              ),
            #binary traits need to be numeric
            groups = groups_traits,
            # grouping variable defined above
            w.type = "optimized",
            opti.maxiter = 300,
            groups.weight = TRUE
          )
          
          f_dist <- gawdis::gawdis(
            data.frame(trait_m) %>%
              mutate_at(
                c(
                  "Eyeless",
                  "Eyes_regression",
                  "Capture_web",
                  "Sensing_web",
                  "No_web",
                  "Tube_web" ,
                  "Sheet_web",
                  "Space_web",
                  "Orb_web",
                  "Ambush_hunter",
                  "Active_hunter",
                  "Food_specialist"
                ),
                as.numeric
              ),
            #binary traits need to be numeric
            groups = groups_traits,
            # grouping variable defined above
            w.type = "equal",
            opti.maxiter = 300,
            groups.weight = FALSE)
          
          
          saveRDS(f_dist,
                  paste0(input_dir, output_dir, "functional_distance.rds")) #storing the data
        }
        
        # Principal Coordinate Analysis -----------------------------------------------------------------------------------
        
        f_dist <-
          readRDS(paste0(input_dir, output_dir, "functional_distance.rds")) #(re)loading the data
        ## Extracting trait axis from Principal component Analysis
        trait_axis <- f_dist %>%
          ape::pcoa() %>% #run principal component analysis using package ape
          pluck("vectors") %>% #select list object using purrr package
          as_tibble()  %>%  #convert to data.frame
          dplyr::select(1:4) %>% #select columns
          rename_with( ~ paste0("Pco", 1:4)) %>%  # rename all columns
          mutate(species_names = trait_parsed$Genus_species) %>% #create column with species rownames
          column_to_rownames("species_names") #select the column species_names to be the rownames of the table
        
        library(vegan)
        source("Scripts/Functions/get_position.R")

        ord <- cmdscale(f_dist,k=3)
        coordinates <- data.frame(ord)
        colnames(coordinates) <- c("PCo1", "PCo2","PCo3")
        centroid <- coordinates %>% as_tibble() #%>%
        #   add_column(family = data$Family) %>%
        #   group_by(family) %>%
        #   summarise(cen.1 = mean(PC1), cen.2 = mean(PC2))
        (fit <- vegan::envfit(ord, trait_m, na.rm = TRUE,choices=c(1,2)))
        
        {
          plot(ord)
          trait_position <- get_position(fit, add = TRUE)
        }
        foo<-ape::pcoa(f_dist,correction="lingoes")
        foo$trace.cor
        procrustes(foo, scale = FALSE)
        trait_pos <- data.frame(trait_position) %>%
          rownames_to_column("Trait") %>%
          filter(!grepl("0",Trait)) |> 
          mutate(Trait =  gsub("1","", Trait)) |> 
          mutate(
            Trait = 
            forcats::fct_recode(Trait,
            "Eyes regression" = "Eyes_regression",
        #    "Eye reduction ratio" = "Eye_reduction_ratio",
            "Tube web" = "Tube_web",
            "Sheet web" = "Sheet_web",
            "Space web" = "Space_web",
            "Orb web" = "Orb_web",
            "Capture web" = "Capture_web",
            "No web" = "No_web",
          #  "Ambush hunter" = "Ambush_hunter",
            "Active hunter" = "Active_hunter",
            "Food specialist" = "Food_specialist",
            "Size dimorphism" = "Sexual_Size_Dimorphism",
            "Body size" = "Body_length_avg",
         #   "Leg elongation" = "Leg_elongation",
            "Femur elongation" = "Femur_elongation",
            "Prosoma shape" = "Prosoma_shape",
            "Pigm.(F)" = "PigmentFully",
            "Pigm.(V)" = "PigmentVariable",
            "Pigm.(P)" = "PigmentPartly",
            "Depigm." = "PigmentDepigmented"
         #   "Soil" = "VerticalitySoil",
          #  "Wall" =  "VerticalityWall",
         #   "Soil+Wall" = "VerticalitySoil + Wall"
          )) |> 
          filter(!Trait %in% c("AME_typeAbsent_Adaptation","AME_typeAbsent_Ontology", "AME_typePresent")) 

          ggplot(coordinates, aes(PCo1, PCo2)) +
            stat_density_2d(
              aes(fill = ..level..),
              geom = "polygon",
              colour = NA,
              alpha = .3,
              h = .25
            ) +
            geom_point(size=.1)+
            geom_hline(aes(yintercept = 0), linetype = 3, colour = "gray70") +
            geom_vline(aes(xintercept = 0), linetype = 3, colour = "gray70") +
           # scale_fill_gradientn(colours = rev(myCol)) +
            ggrepel::geom_text_repel(data = trait_pos |> filter(Dim1 > -0.05), aes(x = Dim1, y = Dim2, label = Trait), xlim=c(0.25,NA), direction="y") +
            ggrepel::geom_text_repel(data = trait_pos |> filter(Dim1 < -0.05), aes(x = Dim1, y = Dim2, label = Trait), xlim=c(NA,-.3), direction="y") +
            geom_point(
              data = trait_pos,
              aes(x = Dim1, y = Dim2),
              shape = 21,
              fill = "white",
              size = 2.5) +
            geom_point(
              data = trait_pos,
              aes(x = Dim1, y = Dim2),
              shape = 19,
              colour = "black",
              size = 1
            ) +
            theme(
              panel.background = element_rect(
                fill = NA,
                colour = "black",
                size = 1,
                linetype = "solid"
              ),
              panel.grid = element_blank(),
              legend.position = "none"
            )  +
            labs(x = "PCoA 1 (43%)", y = "PCoA 2 (35%)") +
            ylim(-.46, .36) + xlim(-.46, .36) +
            #coord_fixed()+
            scale_fill_viridis_c(option="C")
          
          plot_bs <-  trait  %>%
            ggplot(aes(x = expm1(Body_length_avg))) +
            geom_density(fill = myCol[3],
                         alpha = 1,
                         colour = "black") +
            labs(title = "Body Size (mm)", x = NULL, y = "Density") +
            theme(
              plot.background = element_blank(),
              panel.grid = element_blank(),
              axis.line.y = element_line(),
              axis.line.x = element_line(),
              rect = element_blank(),
              axis.title.y = element_text(size = 10, colour = "black"),
              axis.text.x = element_text(size = 10, colour = "black")
            )       
        # Extracting centroid per family
        
        centroid <-
          trait_axis %>%
          add_column(family = trait_parsed$Family) %>%
          group_by(family) %>%
          # summarise(Centroid_x = mean(PCo1), Centroid_y = mean(PCo2))
          summarise(across(starts_with("Pco"), ~ mean(.x, na.rm = TRUE)))
        
        #Plot species coordinates with family centroids
        myCol <- viridis::viridis(n = 7, option = "B")
        
        (
          plot_density <-  trait_axis %>%
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
            geom_density_2d(
              colour = "black",
              adjust = 1.5,
              alpha = .1
            ) +
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
        ggsave(
          plot_density,
          filename = paste0(results_dir, "functional_pcoa.pdf"),
          device = cairo_pdf
        )
        
        
        # Extracting the broad scale variables  ---------------------------------------------------------------------------
        
        
        predictors <-
          raster::stack(list.files(
            path = paste0(input_dir, raster_dir),
            pattern = "*.tif",
            full.names = TRUE
          ))
        
        raster::projection(predictors) <-
          "+proj=longlat +datum=WGS84 +ellps=WGS84"
        names(predictors) <- c("elev",
                               "ice",
                               "karst",
                               "T_mean",
                               "annual_range",
                               "prec")
        
        predictors_regional <- raster::extract(predictors, site[, 4:5])
        
        predictors <-
          site %>%
          dplyr::select(ID,
                        decimalLongitude,
                        decimalLatitude,
                        #entranceNumber,
                        entranceSize,
                        development,
                        #positiveDrop,
                        negativeDrop) %>% # select columns
          add_column(data.frame(predictors_regional)) %>% # add regional predictors
          mutate(across(
            c(development, negativeDrop, , ice, karst, entranceSize,  elev, prec),
            log1p
          )) %>%  #log transform variables
          mutate(T_mean = T_mean / 10) #%>%
        #mutate(karst = log10(karst+1))
        
        predictors %>%
          dplyr::select(where(is.numeric)) %>%
          pivot_longer(everything(), names_to = "var", values_to = "value") %>%
          ggplot(aes(x = value)) +
          geom_histogram() +
          facet_wrap( ~ var, scale = "free") +
          theme_minimal()
        
        ggpairs(predictors, columns = 2:ncol(predictors)) # check collinearity
        
        save(comm_parsed,
             trait_axis,
             predictors,
             file = paste0(input_dir, object_dir, "BBGDM_input.R"))
        