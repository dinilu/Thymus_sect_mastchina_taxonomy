
library(tidyverse)
library(vegan)


# Load genetic groups data ----------------------------------------------

genetic_groups <- "data/genetic groups.xlsx" %>% 
  openxlsx::read.xlsx(2) %>% 
  rename(clade_2 = genetic_group_k2) %>% 
  mutate(clade_2 = case_match(clade_2,
                              "diploid" ~ "2x",
                              "tetraploid" ~ "4x")) %>% 
  mutate(clade_2 = factor(clade_2, 
                          levels = c("2x", 
                                     "4x"))) %>% 
  rename(clade_5 = genetic_group_k5) %>% 
  mutate(clade_5 = str_to_sentence(clade_5)) %>% 
  mutate(clade_5 = factor(clade_5,
                          levels = c("Algarve",
                                     "Cádiz",
                                     "Doñana",
                                     "Hercynian",
                                     "Tetraploid")))


# Load releve data, manipulate them, and check --------------------------

# Read file names
files <- list.files("data/releve/", 
                    pattern = "inventarios tomillos", 
                    full.names = TRUE)

# Read data for each file in a list format
releve <- lapply(files, 
                 FUN = openxlsx::read.xlsx, 
                 sheet = "PLOTS", 
                 skipEmptyRows = FALSE)

# Set names of elements in the list to match population names (e.g. TA1)
names(releve) <- files |>
  stringr::str_remove("data/releve/inventarios tomillos ") |>
  stringr::str_remove(".xlsx") |>
  tibble::tibble() |>
  tidyr::separate("stringr::str_remove(...)", 
                  into = c("Population_ID", "Population_name"),
                  sep = "_") |> 
  dplyr::pull(Population_ID)

# Split data from each population into three independent relevés
split_releve <- function(df) {
  n <- rowSums(is.na(df[, 2:101])) == ncol(df[, 2:101])
  cs <- cumsum(n) + 1
  s <- split(df[!n, ], cs[!n])
  s[[2]] <- s[[2]][-1,]
  s[[3]] <- s[[3]][-1,]
  names(s) <- c("P1", "P2", "P3")
  return(s)
}

releve_plots <- lapply(releve, FUN = split_releve)

# Calculate cover from each species in each relevé
calculate_cover <- function(df) {
  sp_names <- df[,1]
  df <- df[,2:101] %>% 
    dplyr::mutate(across(everything(), 
                         as.numeric)) %>% 
    tibble::tibble()
  data.frame(species = sp_names, 
             cover = rowSums(df, na.rm = TRUE))
}

releve_cover <- lapply(releve_plots, 
                       FUN = function(x){
                         lapply(x, FUN = calculate_cover)})


# Merge with rbind all dataframes in the nested list format
# It returns a long format dataframe
releve_df <- releve_cover %>% 
  lapply(FUN = purrr::list_rbind, 
         names_to = "releve") %>% 
  purrr::list_rbind(names_to = "Population_ID")

# Standardize some names (e.g. Hojarasca, Thymus, etcetera).
releve_df <- releve_df %>% 
  dplyr::mutate(species = dplyr::if_else(
    stringr::str_detect(species, "FRUTO LARGO"),
    "[No determinado]",
    species)) %>% 
  dplyr::mutate(species = dplyr::if_else( # Combine "Hojarasca *" under the same class "Hojarasca".
    stringr::str_detect(species, "Hojarasca"),
    stringr::str_extract(species, "Hojarasca"),
    species)) %>% 
  dplyr::mutate(species = dplyr::if_else(
    stringr::str_detect(species, "Thymus"),
    stringr::str_extract(species, "Thymus"),
    species)) %>% 
  dplyr::group_by(Population_ID, releve, species) %>% 
  dplyr::summarize(cover = sum(cover)) %>% 
  dplyr::mutate(cover = dplyr::if_else( # Limit cover to 100 in Hojarasca
    cover > 100,
    100, 
    cover
  ))

# Check taxonomy
# releve_df %>% pull(species) %>% unique() %>% sort()
# 
# releve_df %>% filter(species == "Hojarasca Pino")

# Fill data for species not found in relevés 
releve_df <- releve_df %>% 
  tidyr::pivot_wider(id_cols = c("Population_ID", "releve"),
                     names_from = species,
                     values_from = cover,
                     values_fill = 0) %>% 
  tidyr::pivot_longer(cols = -c(Population_ID, releve),
                      names_to = "species", 
                      values_to = "cover")  %>% 
  left_join(genetic_groups) 

# Calculate mean coverage at population level
releve_populations <- releve_df %>% 
  dplyr::group_by(Population_ID, species) %>% 
  dplyr::summarise(cover = mean(cover)) %>% 
  mutate(cover = round(cover, digits = 0))

# Transform into community matrix format
population_cm <- releve_populations %>%
  tidyr::pivot_wider(id_cols = "Population_ID",
                     names_from = species,
                     values_from = cover) %>% 
  ungroup()

# Select substrate variables
population_cm_substrate <- population_cm %>% 
  select(c("Population_ID", 
           "Bloque", "Piedra", "Grava", "Arena", "Arcilla",
           "Madera", "Madera pino", "Humus", "Hojarasca", 
           "Excrementos", "Cenizas")) %>% 
  # Remove species presence in less than 2 populations
  # mutate(across(everything(), ~ ifelse(. > 0, 1, 0))) %>% 
  # select(where(~ sum(.x) >= 5)) %>% 
  select(where(~ sum(ifelse(. > 0, 1, 0)) >= 5))

# Select species variables
population_cm_species <- population_cm %>% 
  select(-c("Population_ID", "[No determinado]",
            "Bloque", "Piedra", "Grava", "Arena", "Arcilla",
            "Madera", "Madera pino", "Humus", "Hojarasca", 
            "Excrementos", "Cenizas")) %>%
  # Remove species presence in less than 2 populations
  # mutate(across(everything(), ~ ifelse(. > 0, 1, 0))) %>% 
  # select(where(~ sum(.x) >= 5)) %>% 
  select(where(~ sum(ifelse(. > 0, 1, 0)) >= 5)) %>%
  # Apply quadratic transformation to reduce weight of more frequent species
  # mutate(across(everything(), ~sqrt(.))) %>% 
  # Apply the hellinger transformation instead to see if something improve
  decostand(method = "hellinger") %>% 
  # Get back population names
  mutate(Population_ID = population_cm$Population_ID) %>% 
  select(Population_ID, everything()) 


# Run substrate analysis --------------------------------------------------

# Calculate popluation distances
population_substrate_dist <- population_cm_substrate %>% 
  column_to_rownames(var = "Population_ID") %>% 
  vegan::vegdist(method = "euclidean")

# PCoA 
run_pcoa <- function(dist, plot_gg, plot_pl, width = 7, height = 7) {
  # dist <- population_substrate_dist
  # plot_gg <- "outputs/figures/ordination_substrate_genetic_groups.pdf"
  # plot_pl <- "outputs/figures/ordination_substrate_ploidy.pdf"
  # width <- 7
  # height <- 7
  # pcoa <- cmdscale(dist, k = 2, eig = TRUE)
  pcoa <- ape::pcoa(dist)
  
  # Create data.frame with the first two axis of PCoA
  # df_pcoa <- data.frame(axis1 = pcoa$points[,1],
  #                       axis2 = pcoa$points[,2], 
  #                       Population_ID = rownames(pcoa$points)) %>% 
  df_pcoa <- data.frame(axis1 = pcoa$vectors[,1],
                        axis2 = pcoa$vectors[,2], 
                        Population_ID = rownames(pcoa$vectors)) %>% 
    left_join(genetic_groups, by = "Population_ID") 
  
  getCHarea_gg <- function(x){
    area <- sp::Polygon(coords = x[,1:2], hole = F)@area
    data.frame(clade_5 = unique(x$clade_5), charea = area)
  }
  
  getCHarea_pl <- function(x){
    area <- sp::Polygon(coords = x[,1:2], hole = F)@area
    data.frame(clade_2 = unique(x$clade_2), charea = area)
  }
  
  gg_hull <- df_pcoa %>% 
    group_by(clade_5) %>% 
    slice(chull(axis1, axis2))
  
  gg_ch_area <- gg_hull %>% 
    group_map(~ getCHarea_gg(.), .keep = TRUE) %>% 
    dplyr::bind_rows()
  
  
  # Draw biplot for genetic groups
  
  plot1 <- ggplot(df_pcoa, 
                  aes(x = axis1, 
                      y = axis2, 
                      color = clade_5,
                      fill = clade_5)) +
    geom_point(alpha = 0.3,
               aes(shape = clade_5)) +  # Puntos de cada parcela
    geom_polygon(data = gg_hull, 
                 alpha = 0.05) + 
    geom_point(data = df_pcoa %>% group_by(clade_5) %>% 
               summarize(axis1 = mean(axis1),
                         axis2 = mean(axis2)),
               aes(shape = clade_5),
               size = 5) +
    scale_colour_manual(values =  c(Algarve = "#009E73", 
                                    Cádiz = "#F0E442", 
                                    Doñana = "#CC79A7",
                                    Hercynian = "#E69F00", 
                                    Tetraploid = "#6388b4"
                        )) + 
    scale_fill_manual(values =  c(Algarve = "#009E73", 
                                  Cádiz = "#F0E442", 
                                  Doñana = "#CC79A7",
                                  Hercynian = "#E69F00", 
                                  Tetraploid = "#6388b4"
                      )) +
    scale_shape_manual(values = c(Tetraploid = 22,  # cuadrado
                                  Cádiz = 21,    # círculo
                                  Doñana = 21,
                                  Algarve = 21,
                                  Hercynian = 21
    ))  +
    labs(x = paste0("PCoA axis 1 (",
                    sprintf("%.1f",
                            100 * pcoa$values[1, 2]), 
                    "%)"),
         y = paste0("PCoA axis 2 (", 
                    sprintf("%.1f",
                            100 * pcoa$values[2, 2]),
                    "%)"),
         color = "Contribution") +
    theme_bw() +
    # theme(legend.position = "none") +
    guides(color = guide_legend(title="Group"),
           shape = guide_legend(title="Group"),
           fill = guide_legend(title="Group"))
  
  ggsave(plot_gg, plot1, width = width, height = height)
  
  
  pl_hull <- df_pcoa %>% 
    group_by(clade_2) %>% 
    slice(chull(axis1, axis2))
  
  pl_ch_area <- pl_hull %>% 
    group_map(~ getCHarea_pl(.), .keep = TRUE) %>% 
    dplyr::bind_rows()
  
  # Draw biplot for ploidy levels
  plot2 <- ggplot(df_pcoa, 
                  aes(x = axis1, 
                      y = axis2, 
                      color = clade_2, 
                      fill = clade_2)) +
    geom_point(alpha = 0.3,
               aes(shape = clade_2)) +  # Puntos de cada parcela
    geom_polygon(data = pl_hull, 
                 alpha = 0.05) + 
    geom_point(data = df_pcoa %>% group_by(clade_2) %>% 
                 summarize(axis1 = mean(axis1),
                           axis2 = mean(axis2)),
               aes(shape = clade_2),
               size = 5) +
    scale_colour_manual(values =  c("2x" = "#ffae34",
                                    "4x" = "#6388b4"
                        )) + 
    scale_fill_manual(values =  c("2x" = "#ffae34",
                                  "4x" = "#6388b4"
                      )) +
    scale_shape_manual(values = c("2x" = 21,
                                  "4x" = 22)) +
    theme_bw() +
    labs(x = paste0("PCoA axis 1 (",
                    sprintf("%.1f",
                            100 * pcoa$values[1, 2]), 
                    "%)"),
         y = paste0("PCoA axis 2 (", 
                    sprintf("%.1f",
                            100 * pcoa$values[2, 2]),
                    "%)"),
         color = "Contribution") +
    guides(color = guide_legend(title="Ploidy"),
           shape = guide_legend(title="Ploidy"),
           fill = guide_legend(title="Ploidy"))
    # theme(legend.position = "bottom") + 
    
  ggsave(plot_pl, plot2, width = width, height = height)
  return(list(pcoa = pcoa, 
              df_pcoa = df_pcoa, 
              gg_ch_area = gg_ch_area, 
              pl_ch_area = pl_ch_area,
              gg_plot = plot1,
              pl_plot = plot2))
}

pcoa_substrate <- run_pcoa(population_substrate_dist,
                           "outputs/figures/ordination_substrate_genetic_groups.pdf",
                           "outputs/figures/ordination_substrate_ploidy.pdf")


# Run Releve Analysis -----------------------------------------------------   

# Calculate popluation distances
population_dist <- population_cm_species %>% 
  column_to_rownames(var = "Population_ID") %>% 
  vegan::vegdist(method = "bray")

# PCoA 
population_pcoa <- run_pcoa(population_dist,
                            plot_gg = "outputs/figures/ordination_genetic_groups.pdf",
                            plot_pl = "outputs/figures/ordination_ploidy.pdf")


# Analysis with companion species -----------------------------------------

#### WITH COMPANION SPECIES

# Read data for each file in a list format
companion <- lapply(files, 
                    FUN = openxlsx::read.xlsx, 
                    sheet = "ACOMPAÑANTES", 
                    skipEmptyRows = FALSE)

names(companion) <- files |>
  stringr::str_remove("data/releve/inventarios tomillos ") |>
  stringr::str_remove(".xlsx") |>
  tibble::tibble() |>
  tidyr::separate("stringr::str_remove(...)", 
                  into = c("Population_ID", "Population_name"),
                  sep = "_") |>
  dplyr::pull(Population_ID)


companion_df <- companion %>% 
  purrr::list_rbind(names_to = "Population_ID") %>% 
  rename(species = SPP.ACOMPAÑANTE) %>% 
  mutate(cover = 1) %>% 
  select(Population_ID, species, cover) %>% 
  group_by(Population_ID, species) %>% 
  summarize(cover = max(cover)) %>% 
  ungroup() %>% 
  left_join(genetic_groups) %>% 
  select(Population_ID, species, cover)

relcom_populations <- releve_populations %>% 
  bind_rows(companion_df) %>% 
  group_by(Population_ID, species) %>% 
  summarize(cover = max(cover)) %>% 
  ungroup() 

# Transform into long format
relcom_cm <- relcom_populations %>% 
  tidyr::pivot_wider(id_cols = c("Population_ID"),
                     names_from = species,
                     values_from = cover,
                     values_fill = 0) %>% 
  tidyr::pivot_longer(cols = -c(Population_ID),
                      names_to = "species", 
                      values_to = "cover") %>%
  tidyr::pivot_wider(id_cols = "Population_ID",
                     names_from = species,
                     values_from = cover) %>% 
  ungroup()

relcom_cm_substrate <- relcom_cm %>% 
  select(c("Population_ID", 
           "Bloque", "Piedra", "Grava", "Arena", "Arcilla",
           "Madera", "Madera pino", "Humus", "Hojarasca", 
           "Excrementos", "Cenizas"))

relcom_cm_species <- relcom_cm %>% 
  select(-c("Population_ID", "[No determinado]",
            "Bloque", "Piedra", "Grava", "Arena", "Arcilla",
            "Madera", "Madera pino", "Humus", "Hojarasca", 
            "Excrementos", "Cenizas")) %>%
  # Remove species presence in less than 2 populations
  select(where(~ sum(ifelse(. > 0, 1, 0)) >= 5)) %>%
  # mutate(across(everything(), ~ ifelse(. > 0, 1, 0))) %>% 
  # select(where(~ sum(.x) >= 5)) %>% 
  # Aplica transformación cuadrática para reducir efecto de especies más frecuentes
  # mutate(across(everything(), ~sqrt(.))) %>% 
  # Apply the hellinger transformation instead to see if something improve
  decostand(method = "hellinger") %>% 
  # Get back population names
  mutate(Population_ID = relcom_cm$Population_ID) %>% 
  select(Population_ID, everything()) 


relcom_dist <- relcom_cm_species %>% 
  column_to_rownames(var = "Population_ID") %>% 
  vegan::vegdist(method = "bray")

# PCoA 
population_companion_pcoa <- run_pcoa(relcom_dist,
                                      plot_gg = "outputs/figures/Figure_6a.pdf",
                                      plot_pl = "outputs/figures/Figure_6b.pdf")

fig_6 <- ggpubr::ggarrange(population_companion_pcoa[["pl_plot"]],
                  population_companion_pcoa[["gg_plot"]])

fig_6

ggsave("outputs/figures/Figure_6.pdf", 
       fig_6,
       width = 9,
       height = 3)
ggsave("outputs/figures/Figure_6.png", 
       fig_6,
       width = 9,
       height = 3)





# Summarize communities at the genetic groups level -----------------------

relcom_gg <- relcom_cm_species %>% 
  left_join(genetic_groups) %>% 
  select(-c(Population_name, clade_2)) %>% 
  group_by(clade_5) %>% 
  summarize(across(2:45, ~ mean(.x)))

ifelse(relcom_gg[,-1] > 0, 1, 0)


ifelse(relcom_gg[,-1] > 0, 1, 0) %>% 
  apply(MARGIN = 1, FUN = sum) %>% 
  setNames(relcom_gg$clade_5)




# Species indicator analysis ----------------------------------------------

## Load library
library(indicspecies)

# Calculate IndVal to find indicative species
indval <- population_cm_species %>% 
  column_to_rownames(var = "Population_ID") %>% 
  multipatt(population_pcoa$df_pcoa$clade_5, 
            func = "IndVal.g", 
            control = how(nperm = 999))

# Ver especies con significancia estadística (p < 0.05)
summary(indval)


# Calcular IndVal para encontrar especies indicadoras
indval <- population_cm_species %>% 
  column_to_rownames(var = "Population_ID") %>% 
  multipatt(population_pcoa$df_pcoa$clade_2, 
            func = "IndVal.g", 
            control = how(nperm = 999))

# Ver especies con significancia estadística (p < 0.05)
summary(indval)






# Calcular IndVal para encontrar especies indicadoras
indval <- relcom_cm_species %>% 
  column_to_rownames(var = "Population_ID") %>% 
  multipatt(population_companion_pcoa$df_pcoa$clade_5, 
            func = "IndVal.g", 
            control = how(nperm = 999))

# Ver especies con significancia estadística (p < 0.05)
summary(indval)


# Calcular IndVal para encontrar especies indicadoras
indval <- relcom_cm_species %>% 
  column_to_rownames(var = "Population_ID") %>% 
  multipatt(population_companion_pcoa$df_pcoa$clade_2, 
            func = "IndVal.g", 
            control = how(nperm = 999))

# Ver especies con significancia estadística (p < 0.05)
summary(indval)

