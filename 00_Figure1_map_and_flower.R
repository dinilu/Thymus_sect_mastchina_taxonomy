# Figure 1: study area map + flower schematic
# This script reads sampling locations and genetic group assignments, builds a publication-quality map
# of the Western Mediterranean study region, and combines it with a flower schematic into a two-panel figure.
# Outputs are saved as PDF and PNG in outputs/figures/.
# Reproducibility notes:
# - Requires the input Excel files in data/ and an SVG schematic.
# - Uses UTM zone 30N (EPSG:32630) for mapping and consistent distances.
# - Colours are chosen to be colour-blind friendly.
#
# Load required packages
# - tidyverse: data wrangling + ggplot2
# - sf: spatial objects and coordinate transforms
# - rnaturalearth: country boundaries
# - magick: image handling (used indirectly; kept for compatibility)
#
library(tidyverse)
library(sf)
library(rnaturalearth)
library(magick)

# Input data
# 1) Genetic group assignments (five clades/groups)
# 2) Population table with coordinates
# These are joined to associate each sampled population/individual with its genetic group.
#
# Read genetic group assignments
# - Sheet 2 is used.
# - Group labels are converted to sentence case for consistent plotting/legend labels.
#
genetic_groups <- openxlsx::read.xlsx("data/genetic groups.xlsx", 
                                      sheet = 2) %>% 
  mutate(genetic_group_k5 = str_to_sentence(genetic_group_k5))

# Read sampling locations (populations)
# - Sheet 2 is used.
# - Joined with genetic group information to retain group labels for each location.
#
datos_loc <- openxlsx::read.xlsx("data/Populations.xlsx", 
                                 sheet = 2) %>% 
  right_join(genetic_groups)

# Convert tabular coordinates to an sf point layer
# - Coordinates are provided as Longitude/Latitude in WGS84 (EPSG:4326).
#
aux <- datos_loc %>% 
  sf::st_as_sf(coords = c("Longitude", "Latitude"), 
               crs = 4326)

# Download country boundaries (Natural Earth)
# - Scale 'large' provides higher-resolution coastlines suitable for figures.
#
countries <- rnaturalearth::ne_countries(scale = "large",
                                         returnclass = "sf")

# Filter to the set of countries shown in the study area map
# - Adjust this list if you need a different geographic extent.
#
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


# Create a simple 'ocean' background polygon
# - This is a rectangular polygon covering the plotting window.
# - It is used to provide a subtle blue background for sea areas.
#
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


# Define a colour-blind friendly palette for genetic groups
# - Ensure the order matches the factor levels defined above.
#
daltonic_selection <- c("#009E73", "#F0E442", "#CC79A7", "#E69F00", "#6388b4")


# Build the map panel (Figure 1a)
# - Base layer: ocean background
# - Land layer: country polygons
# - Points: sampling locations coloured/shaped by genetic group
# - coord_sf limits define the figure extent; adjust to zoom in/out.
#
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
  scale_shape_manual(values = c("Tetraploid" = 21,  # Círculo lleno
                                "Algarve" = 22,  # Cuadrado
                                "Doñana" = 22,  # Cuadrado
                                "Cádiz" = 22,  # Cuadrado
                                "Hercynian" = 22),
                     name = "Group") +  # Cuadrado
  scale_fill_manual(values = daltonic_selection,
                    name = "Group") +
  theme_minimal()

# Render the map in the R graphics device (interactive check).
#
final_map


# Load the flower schematic (Figure 1b)
# - The schematic should be an SVG located at data/flower_schema.svg.
# - Using cowplot keeps the graphic as a grob for arrangement with the map.
#
img_flower <- cowplot::ggdraw() +
  cowplot::draw_image("data/flower_schema.svg")

# Combine panels into a single two-column figure
# Two alternative approaches are included below:
# - cowplot::plot_grid (commented as an example)
# - ggpubr::ggarrange (used as the final assignment to final_plot)
# The final object saved is the one assigned last to 'final_plot'.
#
# Example with PNG (for fun, the OP's avatar - I love the raccoon)
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

final_plot

# Export figure to disk
# - Save both PDF (vector) and PNG (raster) versions.
# - Ensure outputs/figures/ exists before running (or create it).
#
ggsave("outputs/figures/Figure_1.pdf",
       final_plot,
       width = 7,
       height = 3)
# Export figure to disk
# - Save both PDF (vector) and PNG (raster) versions.
# - Ensure outputs/figures/ exists before running (or create it).
#
ggsave("outputs/figures/Figure_1.png",
       final_plot,
       width = 7,
       height = 3)
