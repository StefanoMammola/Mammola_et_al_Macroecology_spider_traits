pacman::p_load(
  "tidyverse",
  "hypervolume",
  "doSNOW",
  "future",
  "gdm",
  "surveillance",
  "plotrix",
  "glue",
  "magrittr",
  "BAT",
  "ggpp",
  "scales",
  "patchwork",
  "gghighlight",
 "nlme",
 "ggridges",
 "janitor",
 "gt",
 "gtExtras"
  )


#directories to save and load files ....................................................................................
input_dir <- "Data/" #Folder where the data are located
object_dir <- "objects/" #Sub-folder where data are located
output_dir <- "null_model/" #Folder where to save the data
raster_dir <- "rasters/" # Folder where the raster data are stored
results_dir <-
  "Results/" #Folder where to save figures, tables, etc
table_dir <- "tables/" # Folder where the excel files are saved
plot_colour <- "mediumorchid4"

# Funtions  -------------------------------------------------------------------------------------------------------

#SES curves
ses_tidy <- function(obj, var,param = FALSE) {
  res <- data.frame(
    value =
      BAT::ses(
        obs = obj %>%  pull(all_of(var)) %>% extract(1) ,
        est = obj %>%  pull(all_of(var)) %>% extract(-1),
        param = param
      )
  ) %>%  set_rownames(paste0(c("SES_", "p_"), var))
  return(res)
}



#' -----------------------------------------------------------------------------------------------------------------
# Alpha Diversity ----
#' -----------------------------------------------------------------------------------------------------------------

alpha.null <- readRDS(paste0(input_dir,object_dir,"SES_FR.rds"))
load(paste0(input_dir,object_dir,"BBGDM_input.R"))

alpha_ses <- alpha.null %>%
  bind_rows() %>%
  pivot_longer(everything(), names_to = "caves", values_to = "alpha") %>%
  split(.$caves)  %>%
  map( ~ .x %>%
         ses_tidy("alpha", param=FALSE) %>%
         rownames_to_column("metric")) %>%
  bind_rows(.id= "ID") %>% 
  filter(metric == "SES_alpha")

# Frequency of over/underdispersion
table_disp<- alpha_ses %>% 
  mutate(d = value > 0) %>% 
  pull("d") %>%
  janitor::tabyl() %>%
  mutate(percent = percent*100) 

figure_2b <-
  alpha.null %>%  
  bind_rows() %>%
  pivot_longer(everything(), names_to = "caves",values_to = "alpha") %>% 
  group_split(caves) %>% 
  map(~ .x %>% 
        ses_tidy("alpha",param=FALSE) %>%  rownames_to_column("metric")) %>% 
  bind_rows() %>% 
  mutate(ID = sort(rep(1:367, 2))) %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  mutate(signif = !between(p_alpha, 0.025, 0.975),
         dispersion = ifelse(SES_alpha > 0, TRUE, FALSE)) %>% 
  ggplot()+
  geom_density_ridges_gradient(aes(x=SES_alpha,  y = 1, fill=..x..), colour= NA, show.legend=FALSE,
                               jittered_points = TRUE,
                               position = position_points_jitter(width = 0.05, height = 0),
                               point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7, point_colour="gray10")+
   geom_vline(xintercept = 0, linetype = 3 , size=.5) + 
  geom_segment(aes(xend = 0.2, x = 2.5 , y = 1.75, yend = 1.75), arrow = arrow(length=unit(0.15,"cm"), ends="first", type = "open"), colour = "gray40") + 
  geom_segment(aes(xend = -0.7, x = -3 , y = 1.75, yend = 1.75), arrow = arrow(length=unit(0.15,"cm"), ends="first", type = "open"), colour = "gray40") + 
  annotate("text", y = 1.8, x = -2, label = glue::glue("Underdispersion ({round(table_disp,2)$percent[1]}%)"), hjust=.5)+
  annotate("text", y = 1.8, x = 1.3, label = glue::glue("Overdispersion ({round(table_disp,2)$percent[2]}%)"), hjust=.5)+
     scale_fill_gradient2( low = "#0d0887FF", 
                         mid = "white",
                         high = "orange")+
  theme_classic() + 
  theme(axis.title.y = element_blank()) + 
  labs(x= "Standardized effect size (SES)")


model_data <- alpha_ses %>%
  left_join( 
  predictors, by ="ID"
  ) %>%  
  filter(ID != "CAVE_190") %>% 
  mutate(across(entranceSize:prec, scale))

pred_vars <- model_data %>%  select(-ID, -metric, - value, -decimalLatitude, decimalLongitude)

formula1 <- as.formula(paste("value ~",paste(glue::glue("s({colnames(pred_vars)})",), collapse="+")))
formula2 <- as.formula(paste("value ~",paste(glue::glue("{colnames(pred_vars)}",), collapse="+")))

#'check duplicated coordinates  ------------------------------------------
names_keep<- model_data %>%  
  distinct(decimalLatitude,decimalLatitude, .keep_all = TRUE) %>% 
  pull(ID)

setdiff(model_data$ID, names_keep)
#' ----------------------------------------------------------------------

mod2<- lm(formula2, data = model_data)
mod3<- gls(value ~ entranceSize + development + negativeDrop + elev + ice + 
             karst + T_mean + annual_range + prec, data= model_data, na.action = na.omit,
               correlation = corExp(form = ~decimalLatitude  + decimalLongitude, nugget = TRUE))

summary(mod2)
summary(mod3)

performance::check_autocorrelation(mod3)
performance::check_collinearity(mod3)

table_gls <- mod3 %>%  
  broom.mixed::tidy() %>% 
  add_column(VIF = c(NA,car::vif(mod3))) %>% 
  mutate(term = recode_factor( factor(term),  
                                   "(Intercept)" = "Intercept_ ",
                                   development = "Development_[m]",
                                   elev = "Elevation_[m]",
                                   entranceSize="Entrance size_[m²]",
                                   Geographic="Geographic distance_[km]",
                                   ice="LGM ice distance_[km]",
                                   karst="Karst area_[km²]",
                                   negativeDrop="Negative drop_[m]",
                                   annual_range="Annual range_[ºC]",
                                   prec="Precipitation_[mm]",
                                   T_mean="Temperature_[ºC]"
  )) %>% 
  separate(term, into = c("term","unit"), sep="_") %>% 
  gt::gt() %>% 
  tab_header(
    title = md("**Generalised Least Squares**"),
    subtitle = md("*Response* = Standard effect sizes of Functional richness")
  ) %>%
  fmt_number(
    columns = c(3,4,5,7),
    decimals = 2
  ) %>%
  fmt_number(
    columns = 6,
    decimals = 6
  ) %>% 
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = p.value,
      rows = p.value <= 0.05
    )
  )  %>% 
  cols_label(
    estimate = "Estimate",
    std.error = "Standard Error",
    statistic = "t-statistic",
    p.value = "p-value",
    VIF="VIF"
  ) %>% 
  cols_align("left", term) %>% 
  tab_footnote(
    footnote = "Variance Inflation Factor.",
    locations = cells_column_labels(columns = VIF)
  ) %>% 
  sub_missing(columns=everything(), rows = everything(), missing_text = "") %>% 
  gt_theme_538() %>% 
  gt_merge_stack(col1 = term, col2 = unit, palette = c("gray40","gray60"))
  
table_gls %>%  gtsave("table_gls.rtf")

  

# Map of Europe with SES ------------------------------------------------------------------------------------------

library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)

local_rich <- comm_parsed %>% 
  rownames_to_column ("ID") %>% 
  column_to_rownames("ID") %>% 
  vegan::specnumber() %>% 
  data.frame(richness=.) %>% 
  rownames_to_column ("ID") 
  
world <- ne_countries(scale = "medium", returnclass = "sf")


figure_2a<- ggplot(world) +
  geom_sf(size=.5, colour="gray80", fill="white") +
  geom_point(data = model_data %>% 
               left_join(local_rich, by="ID"), aes(x = decimalLongitude, y= decimalLatitude, fill=value, size= richness), shape=21, alpha=.8) + 
  scale_fill_gradient2( low = "#0d0887FF", 
                        mid = "white",
                        high = "orange", limits = c(-2.5,2.5))+
  coord_sf(xlim = c(-10,35), ylim = c(32,70), expand = TRUE) + 
  theme_classic() + 
  labs(x="Longitude (degrees)", y ="Latitude (degrees)", size = "Species richness", fill = "Standardized effect size (SES) ") + 
  theme(legend.position = "bottom") +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, barwidth=10),
         size = guide_legend(title.position="top", title.hjust = 0.5))



figure_2c<- mod3 %>% 
summary %>% 
  extract2("tTable") %>% 
  data.frame %>% rownames_to_column("Variable") %>% 
  slice(-1) %>% 
  dplyr::rename(SE = 3, t_value = 4, p_value = 5) %>% 
  mutate(Variable = recode_factor( Variable,  
                                   development = "Development [m]",
                                   elev = "Elevation [m]",
                                   entranceSize="Entrance size [m²]",
                                   Geographic="Geographic distance [km]",
                                   ice="LGM ice distance [km]",
                                   karst="Karst area [km²]",
                                   negativeDrop="Negative drop [m]",
                                   annual_range="Annual range [ºC]",
                                   prec="Precipitation [mm]",
                                   T_mean="Temperature [ºC]"
  )) %>% 
  mutate(Variable = fct_reorder(Variable, Value)) %>% 
  mutate(signif = case_when(p_value < 0.05 ~ "TRUE", TRUE ~ "FALSE")) %>% 
  
 ggplot(aes(y=Variable, x=Value,colour = signif, fill = signif)) +
    geom_vline(lty = 3, size = 0.7, col = "grey50", xintercept = 0) +
    geom_errorbarh(aes(xmin = Value-SE, xmax = Value+SE), height=0.1, show.legend= FALSE) +
    geom_text(aes(y=Variable, x=Value, label = round(Value,2)),
              show.legend = FALSE,
              vjust = -1, size = 3) +
    geom_point(size = 2, pch = 21, show.legend= FALSE) +
  scale_colour_manual(values= rev(c("gray10","gray70")))+
   scale_fill_manual(values= rev(c("gray10","gray70"))) + 
  labs(x= "Estimate (± SE)",
       subtitle =glue::glue("N = {summary(mod3) %>% 
                            extract2('dims') %>% 
                            pluck('N')}; R² = {round(
                            as.numeric(
                            performance::r2(mod3)
                            ),2)}"))+
  theme_classic() +
    theme(axis.title.y = element_blank())


Figure_2 <- {figure_2a  | 
  {figure_2b   / figure_2c}} +
  plot_annotation(tag_levels = "a", theme = theme(plot.title= element_text(face="bold")))  +
 plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        plot.tag = element_text(face= "bold"))


ggsave("Results/Figure2.pdf", Figure_2, width = 14, height = 10, device = cairo_pdf)


# Beta total ------------------------------------------------------------------------------------------------------

comp <- "Btotal"

# Load objects  ---------------------------------------------------------------------------------------------------
names <-  paste0(input_dir,output_dir) %>%
  list.files %>%
  tools::file_path_sans_ext()

models <-
  paste0(input_dir,output_dir) %>%  list.files(full.names = TRUE) %>% map(readRDS) %>%  set_names(names)


#Rank variables by effect size ---------------
order <- 
  models %>%
  #access the first iteration and extract the curves for the beta total component
  pluck("Null_BBGDM_1", "coefficient", comp) %>% 
  arrange(desc(effect_size)) %>% 
  rownames_to_column("var") %>%
  select("var") %>% 
  mutate(var = factor(var, levels=var)) %>% 
  mutate(var = recode_factor( var,  
                              development = "Development\n[m]",
                              elev = "Elevation\n[m]",
                              entranceSize="Entrance\nsize [m²]",
                              Geographic="Geographic\ndistance [km]",
                              ice="LGM ice\ndistance [km]",
                              karst="Karst\narea [km²]",
                              negativeDrop="Negative\ndrop [m]",
                              annual_range="Annual range\n[ºC]",
                              prec="Precipitation\n[mm]",
                              T_mean="Temperature\n[ºC]")) %>%  pull(var)

# Get curves for original variables (Panel 1) ---------------------------------------------------------------------

panel_1 <- models %>%
  #access the first iteration and extract the curves for the beta total component
  pluck("Null_BBGDM_1", "curves", comp) %>%
  mutate(variable = recode_factor( variable,  
                                   development = "Development\n[m]",
                                   elev = "Elevation\n[m]",
                                   entranceSize="Entrance\nsize [m²]",
                                   Geographic="Geographic\ndistance [km]",
                                   ice="LGM ice\ndistance [km]",
                                   karst="Karst\narea [km²]",
                                   negativeDrop="Negative\ndrop [m]",
                                   annual_range="Annual range\n[ºC]",
                                   prec="Precipitation\n[mm]",
                                   T_mean="Temperature\n[ºC]"
  )) %>% 
  mutate(variable = factor(variable, levels=order)
  )  %>% 
  #plot the mean value of the splines
  ggplot(aes(value, X50.)) +
  geom_line() +
  #plot the error
  geom_ribbon(aes(ymin=X5., ymax=X95.), alpha=.2) +
  #split plot by variable but allow only a single column
  facet_wrap( ~ variable, scale = "free", ncol = 1)  +
  theme_classic() +
  theme(strip.text = element_blank(),
        axis.text = element_text(size= 12),
        axis.title = element_text(size=14), 
        plot.background = element_rect(fill=NA, colour= NA),
        panel.background = element_rect(fill = NA, colour = NA ))+
  scale_y_continuous(breaks=breaks_extended(4))


# Object to plot histogram ( panel 2 ) ----------------------------------------------------------------------------

Null_coef <- models %>%
  #map over each iteration (first level of the list)
  map(~ .x %>% 
        #within each iteration, extract the coefficient data ( second level of the list)
        pluck("coefficient") 
        ) %>%
  #transpose the list such that the beta diversity component becomes the first level of the list
  purrr::transpose() %>%
  #map over each beta diversity component (current first level of the list) 
  map(
    ~ .x %>%
  #additional mapping to go through each iteration again    
      map( ~ .x %>% 
             #convert coefficient information into data frame
             data.frame %>%
             #convert row names into columns (this is actually the variables name)
             rownames_to_column("var")   
           ) %>%
    #remove empy objects (these were generated by failed runs in the super computer)
      compact %>%
    #combine second-order lists and import names to identify each iteration
      bind_rows(.id = "iteration")
  ) %>%
  #combine first-order lists and import names to identify each beta diversity component
  bind_rows(.id = "Component") %>%
  #remove variables generated by failed runs 
  filter(var != 1) %>%
  #keep only beta diversity component, variables and effect size (i.e., sum of spline coefficients)
  #select(Component, var, effect_size) %>% 
  #rename column "effect_size" to "coefficient" to not confound with the effect size estimated afterwards
  rename(coefficient  = "effect_size") %>% 
  group_by(Component, var) %>% 
  add_count(coefficient) %>% 
  ungroup()  %>% 
  mutate(var = recode_factor( var,  
                              development = "Development\n[m]",
                              elev = "Elevation\n[m]",
                              entranceSize="Entrance\nsize [m²]",
                              Geographic="Geographic\ndistance [km]",
                              ice="LGM ice\ndistance [km]",
                              karst="Karst\narea [km²]",
                              negativeDrop="Negative\ndrop [m]",
                              annual_range="Annual range\n[ºC]",
                              prec="Precipitation\n[mm]",
                              T_mean="Temperature\n[ºC]"
  )) %>% 
  mutate(var = factor(var, levels=order))



ES_p_value <- #create object for p_values of the null distribution (non-parametric)
  Null_coef %>%
  #divide object into lists based on the beta diversity component 
  split(f = as.factor(.$Component)) %>%
  #go through each beta diversity component (first level list)
  map( ~ .x %>%
         #divide the objects within each beta diversity component into lists based on variable name
         split(f = as.factor(.$var)) %>%
         #go trhough each variable object (second level list) to estimate effect size (SE) and p value of the null distribution
         map(
           ~ data.frame(
             value = BAT::ses(
               obs = .x$coefficient[1], #here we select only the value for actual trait composition
               est = .x$coefficient[-1], #here we select the values for null trait composition (999 runs)
               param = FALSE #set this to TRUE to gave a parametric SES estimation
             )
           ) %>%
             #row names are the SES and p value, convert it into a column called "Metric"
             rownames_to_column("Metric")
         )) %>%
  # go trough each variable object (second level list) and combine into a data frame 
  # import name to identify the variables (.id  = var)
  map(~ bind_rows(.x, .id = "var")) %>%
  #combine the beta diversity components (first level list) into a single data frame
  #import name to identify the components (.id = "Component")
  bind_rows(.id = "Component") %>%
  #spread the column "metric" into two new columns (SES, p_value)
  pivot_wider(names_from = Metric, values_from = value) %>%
  #rename column "p-value"
  rename(p = "p-value") %>% 
  full_join(
    #add observed values
  Null_coef %>%
    filter(iteration %in% "Null_BBGDM_1") %>% 
    select(-.) 
  ) %>%   mutate(var = recode_factor( var,  
                                      development = "Development\n[m]",
                                      elev = "Elevation\n[m]",
                                      entranceSize="Entrance\nsize [m²]",
                                      Geographic="Geographic\ndistance [km]",
                                      ice="LGM ice\ndistance [km]",
                                      karst="Karst\narea [km²]",
                                      negativeDrop="Negative\ndrop [m]",
                                      annual_range="Annual range\n[ºC]",
                                      prec="Precipitation\n[mm]",
                                      T_mean="Temperature\n[ºC]"
  )) %>% 
  mutate(var = factor(var, levels=order))

panel_2<- 
Null_coef %>%
  filter(Component %in% comp) %>%
  ggplot(aes(x = coefficient)) +
  geom_histogram(fill="gray10") +
  geom_segment(
    data = ES_p_value %>% filter(Component %in% comp),
    aes(x = coefficient, xend = coefficient, y = 0, yend=250, colour=p > 0.975),
    linetype = 1
  ) +
  geom_point(
    data = ES_p_value %>% filter(Component %in% comp),
    aes(x = coefficient,y=250, colour = p>0.975)
  ) +
  geom_text_npc(data = ES_p_value %>% filter(Component %in% comp),
                aes(
                  npcx = .5,
                  npcy = .95,
                  label = glue::glue("SES = {round(ses,3)}\n p-value = {round(p,3)}")
                ), size=4, hjust=0) +
  facet_wrap( ~ var, scale = "free_x", ncol = 1, shrink=TRUE, strip.position="right")  +
  theme_classic()+
  theme(strip.text = element_blank(),
        legend.position="none",
        axis.text = element_text(size= 12),
        axis.title = element_text(size=14),
        plot.background = element_rect(fill=NA, colour= NA),
        panel.background = element_rect(fill = NA, colour = NA)) +
  scale_y_continuous(breaks=breaks_extended(3), limits=c(0,500))+
  scale_colour_manual(values= rev(c(plot_colour, colorspace::lighten(plot_colour,.6))))

models %>%
  #access the first iteration and extract the curves for the beta total component
  pluck("Null_BBGDM_1", "curves", comp) %>% 
  distinct(variable)

Null_curves <-
  models %>%
  map( ~ .x %>%  pluck("curves")) %>%
  purrr::transpose() %>%
  map(
    ~ .x %>%
      map(~ .x %>% data.frame) %>%
      compact %>%
      bind_rows(.id = "iteration") %>%
      group_by(variable, iteration) %>%
      mutate(code = row_number(value)) %>%
      ungroup()
  ) %>%
  map(
    ~ .x %>%
      split(f = as.factor(.$variable)) %>%
      #variable list
      map(
        ~ .x %>%
          split(f = as.factor(.$code)) %>%
          #code list
          map(
            ~ rbind(ses_tidy(.x, "X5."),
                    ses_tidy(.x, "X50."),
                    ses_tidy(.x, "X95.")) %>%
              rownames_to_column("metric")
          ) %>%
          bind_rows(.id = "code")
      ) %>%
      #variable list
      bind_rows(.id = "variable")
  ) %>%
  #Component list
  bind_rows(.id = "component")  %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  mutate(code = as.numeric(code)) %>%
  full_join(
    models %>%
      pluck(1, "curves", 1) %>%
      group_by(variable) %>%
      mutate(code = row_number(value)) %>%
      ungroup() %>%
      select(value, code, variable)) %>% 
  mutate(variable = recode_factor( variable,  
                              development = "Development\n[m]",
                              elev = "Elevation\n[m]",
                              entranceSize="Entrance\nsize [m²]",
                              Geographic="Geographic\ndistance [km]",
                              ice="LGM ice\ndistance [km]",
                              karst="Karst\narea [km²]",
                              negativeDrop="Negative\ndrop [m]",
                              annual_range="Annual range\n[ºC]",
                              prec="Precipitation\n[mm]",
                              T_mean="Temperature\n[ºC]"
  )) %>% 
  mutate(variable = factor(variable, levels=order))

Null_curves |> filter(variable == "Annual range\n[ºC]")

panel_3 <- Null_curves %>%
  filter(component == comp, code > 1) %>%
  mutate(signif =case_when(p_X50. < 0.025 ~ TRUE,
                             p_X50. > .975 ~ TRUE , TRUE ~ FALSE)) %>% 
  ggplot(aes(
    x = value,
    y = SES_X50.,
  )) +
  geom_hline(yintercept = 0, linetype=3)+
  geom_line(size=2,colour = colorspace::lighten(plot_colour,.8)) +
  geom_line(data=~.x |> filter(signif == TRUE), colour=plot_colour,size=2)+
#  gghighlight(if_any(p_X50.) > .975,
#              unhighlighted_params = list(colour = colorspace::lighten(plot_colour,.8), size=2)) +
  facet_wrap(
    ~ variable,
    scale = "free_x",
    ncol = 1,
    shrink = TRUE,
    strip.position = "right"
  )  +
  theme_classic() +
  scale_linetype_manual("Significant",
                        values = c(1, 2),
                        labels = c("True", "False")) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 12),
        axis.text = element_text(size= 12),
        axis.title= element_text(size=14)) +
  scale_y_continuous(breaks = breaks_extended(4)) 


figure_3<- (panel_1 + labs(y="Change in Beta diversity",x="Gradient") + 
                 panel_2 + labs(y="Counts", x = "Total change in beta diversity")  +
                 panel_3 + labs(y="Standardized effect size (SES)",x="Gradient")+
                 plot_annotation(tag_level = "a") ) & 
  theme(plot.tag = element_text(face="bold"))



ggsave(figure_3, filename="Results/Figure3.pdf", height=17, width=11, device=cairo_pdf)
 