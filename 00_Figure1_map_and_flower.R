# ======================================================================
# Script: 00_Figure1_map_and_flower.R
# Purpose: Generate Figure 1: geographic sampling map of populations by genetic group alongside a floral schema panel.
#
# This script is part of the reproducible analysis accompanying the manuscript.
# It is intended to be run from the project root directory so that relative paths
# (e.g. 'data/' and 'outputs/') resolve correctly.
#
# Commenting convention:
# - Section headers are delimited by '=' or '-' rulers.
# - Short, action-oriented comments precede the code blocks they describe.
# - Existing code lines are left unmodified; only comment lines are added.
# ======================================================================
# ----------------------------------------------------------------------
# Package requirements
# ----------------------------------------------------------------------
# Load all R packages required for data import, manipulation, modelling,
# and figure/table generation.
library(tidyverse)
library(sf)
library(rnaturalearth)
library(magick)


# ----------------------------------------------------------------------
# LOAD DATA
# ----------------------------------------------------------------------
#### LOAD GENETIC DATA  ####
genetic_groups <- openxlsx::read.xlsx("data/genetic groups.xlsx", 
                                      sheet = 2) %>% 
  mutate(genetic_group_k5 = str_to_sentence(genetic_group_k5))

#### LOAD POPULATIONS DATA  ####
datos_loc <- openxlsx::read.xlsx("data/Populations.xlsx", 
                                 sheet = 2) %>% 
  right_join(genetic_groups)

aux <- datos_loc %>% 
  sf::st_as_sf(coords = c("Longitude", "Latitude"), 
               crs = 4326)

countries <- rnaturalearth::ne_countries(scale = "large",
                                         returnclass = "sf")
countries_filtered <- countries %>% 
  filter(name_en %in% c("Spain", 
                        "France", 
                        "Portugal", 
                        "Morocco",
                        "Algeria",
                        "Andorra"))

countries_filtered <- st_transform(countries_filtered, 
                                  crs = 32630)

aux <- st_transform(aux, 
                    crs = 32630) %>% 
  mutate(genetic_group_k5 = factor(genetic_group_k5,
                                levels = c("Algarve", 
                                           "Cádiz",
                                           "Doñana",
                                           "Hercynian",
                                           "Tetraploid")))


ocean <- st_polygon(list(cbind(c(seq(-160000, 1060000, len = 100), 
                                 rep(1060000, 100), 
                                 seq(1060000, -160000, len = 100), 
                                 rep(-160000, 100)),
                               c(rep(3800000, 100),
                                 seq(3800000, 5000000, len = 100),
                                 rep(5000000, 100),
                                 seq(5000000, 3800000, len = 100))))) %>% 
  st_sfc(crs = 32630) %>% 
  st_as_sf() 


clade_colors <- c("#009E73", "#F0E442", "#CC79A7", "#E69F00", "#6388b4")


# Construct the plot object.
final_map <- ggplot(data = countries_filtered) +
  geom_sf(data = ocean, fill = "#8080ff80") +
  geom_sf(fill = "lightgrey", 
          color = "black") + 
  theme_minimal() +
  geom_sf(data = aux,
             aes(fill = genetic_group_k5, 
                 shape = genetic_group_k5),
             alpha = 0.5,
             color = "black",
             size = 2.5) +
  coord_sf(xlim = c(-110000, 1050000),
           ylim = c(3900000, 4900000),
           expand = FALSE) +
  scale_shape_manual(values = c("Tetraploid" = 21, 
                                "Algarve" = 22, 
                                "Doñana" = 22,  
                                "Cádiz" = 22, 
                                "Hercynian" = 22),
                     name = "Group") +

  scale_fill_manual(values = clade_colors,
                    name = "Group") +
  theme_minimal()

final_map


img_flower <- cowplot::ggdraw() +
  cowplot::draw_image("data/flower_schema.svg")

final_plot <- cowplot::plot_grid(final_map,
          img_flower,
          ncol = 2,
          rel_widths = c(2, 1.25),
          rel_heights = c(0.9, 1),
          labels = c("a", "b"))


final_plot <- ggpubr::ggarrange(final_map,
                                img_flower,
                                ncol = 2,
                                labels = c("a", "b"),
                                widths = c(2, 1.3))

# Export figure files in publication-ready formats.
ggsave("outputs/figures/Figure_1.pdf",
       final_plot,
       width = 7,
       height = 3)
ggsave("outputs/figures/Figure_1.png",
       final_plot,
       width = 7,
       height = 3)
