# ======================================================================
# Script: 02_Main_morphology_population_table.R
# Purpose: Summarise morphological, pollen and stomatal traits at the population level; export population-level summary tables.
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
# devtools::install_github("psyteachr/introdataviz")
# ----------------------------------------------------------------------
# Package requirements
# ----------------------------------------------------------------------
# Load all R packages required for data import, manipulation, modelling,
# and figure/table generation.
library(tidyverse)
library(reshape2)
library(FSA)
library(factoextra)
library(ggpubr)


# ----------------------------------------------------------------------
# LOAD DATA
# ----------------------------------------------------------------------

#### LOAD GENETIC GROUPS 
genetic_groups <- "data/genetic groups.xlsx" %>% 
  openxlsx::read.xlsx(2) %>% 
  mutate(genetic_group_k2 = case_match(genetic_group_k2,
                                       "tetraploid" ~ "4x",
                                       "diploid" ~ "2x")) %>% 
  mutate(genetic_group_k5 = str_to_sentence(genetic_group_k5))


#### LOAD MORPHOMETRIC DATA  
morph_measures <- "data/Morphometric_measures.xlsx" %>% 
  openxlsx::read.xlsx(2)

# Inspect object structure and summary statistics.
head(morph_measures)  
summary(morph_measures)


# Combine morphometric and genetic group data
morph_measures <- morph_measures %>% 
  left_join(genetic_groups) %>% 
  rename(clade_5 = genetic_group_k5) %>% 
  rename(clade_2 = genetic_group_k2)

# Create column with ploidy level
morph_measures <- morph_measures %>% 
  mutate(across(c(Population_ID, 
                  Sex, 
                  clade_2,
                  clade_5), 
                as.factor))

# Check populations and genetic groups
table(morph_measures$clade_5, morph_measures$Population_ID)
table(morph_measures$clade_2, morph_measures$Population_ID)
table(morph_measures$clade_2, morph_measures$clade_5)


# LOAD POLLEN AND STOMA DATA

# Load pollen data
pollen_measures <- openxlsx::read.xlsx("data/pollen.xlsx", 2)

# Inspect object structure and summary statistics.
head(pollen_measures)  
summary(pollen_measures)

# Combine pollen data and genetic groups
pollen_measures <- pollen_measures %>% 
  left_join(genetic_groups) %>% 
  rename(clade_5 = genetic_group_k5) %>% 
  rename(clade_2 = genetic_group_k2)

# Load stoma data
stoma_measures <- openxlsx::read.xlsx("data/stoma.xlsx", 2)

# Inspect object structure and summary statistics.
head(stoma_measures)  
summary(stoma_measures)

# Combine stoma data and genetic groups
stoma_measures <- stoma_measures %>% 
  left_join(genetic_groups) %>% 
  rename(clade_5 = genetic_group_k5) %>% 
  rename(clade_2 = genetic_group_k2)


# ----------------------------------------------------------------------
# ANALYSE VARIABLES
# ----------------------------------------------------------------------

#### morph

morph_vars <- c("co",
                "ca", 
                "ltooth",
                "stooth", 
                "hair",
                "fl_infl",
                "d_infl",
                "prop_ltooth")

stoma_vars <- c("stoma_w", "stoma_l")

pollen_vars <- c("d_max", "d_min")

# Create morph summary by population
morph_pop_summary <- lapply(morph_vars, 
                            FUN = \(x)Rmisc::summarySE(morph_measures, 
                                                       measurevar = x, 
                                                       groupvars = c("Population_ID"), 
                                                       na.rm = TRUE)) %>% 
  lapply(FUN = \(x){colnames(x)[3] <- "mean"; x})

# Create morph minmax by genetic group
morph_pop_minmax <- lapply(morph_vars, 
                           FUN = \(x)summarise(group_by(morph_measures,
                                                        Population_ID),
                                               min = min(get(x),
                                                         na.rm = TRUE),
                                               max = max(get(x),
                                                         na.rm = TRUE)))

# Set names in both summaries
names(morph_pop_summary) <- 
  names(morph_pop_minmax) <- 
  morph_vars


# Pivot both summaries to wider format
morph_pop_summary <- melt(morph_pop_summary, 
                            measure.vars = c("N", "mean", "sd", "se", "ci")) %>% 
  pivot_wider(names_from = "variable",
              values_from = "value") %>% 
  rename(variable = L1)

morph_pop_minmax <- melt(morph_pop_minmax, 
                           measure.vars = c("min", "max")) %>% 
  pivot_wider(names_from = "variable",
              values_from = "value") %>% 
  rename(variable = L1)


# Create pollen summary by genetic group
pollen_pop_summary <- lapply(pollen_vars, 
                               FUN = \(x)Rmisc::summarySE(pollen_measures, 
                                                          measurevar = x, 
                                                          groupvars = c("Population_ID"), 
                                                          na.rm = TRUE)) %>% 
  lapply(FUN = \(x){colnames(x)[3] <- "mean"; x})


# Create pollen minmax by genetic group
pollen_pop_minmax <- lapply(pollen_vars,
                            FUN = \(x)summarise(group_by(pollen_measures, 
                                                         Population_ID),
                                                min = min(get(x),
                                                          na.rm = TRUE),
                                                max = max(get(x), 
                                                          na.rm = TRUE)))


# Set names in both summaries
names(pollen_pop_summary) <- 
  names(pollen_pop_minmax) <- 
  pollen_vars

# Pivot both summaries to wider format
pollen_pop_summary <- melt(pollen_pop_summary, 
                             measure.vars = c("N", "mean", "sd", "se", "ci")) %>% 
  pivot_wider(names_from = "variable",
              values_from = "value") %>% 
  rename(variable = L1)

pollen_pop_minmax <- melt(pollen_pop_minmax, 
                            measure.vars = c("min", "max")) %>% 
  pivot_wider(names_from = "variable",
              values_from = "value") %>% 
  rename(variable = L1)


## stoma ##

# Create stoma summary by genetic group
stoma_pop_summary <- lapply(stoma_vars, 
                              FUN = \(x)Rmisc::summarySE(stoma_measures, 
                                                         measurevar = x, 
                                                         groupvars = c("Population_ID"), 
                                                         na.rm = TRUE)) %>% 
  lapply(FUN = \(x){colnames(x)[3] <- "mean"; x})

# Create stoma minmax by ploidy level
stoma_pop_minmax <- lapply(stoma_vars, 
                           FUN = \(x)summarise(group_by(stoma_measures, 
                                                        Population_ID),
                                               min = min(get(x),
                                                         na.rm = TRUE),
                                               max = max(get(x), 
                                                         na.rm = TRUE)))

# Set names in both summaries
names(stoma_pop_summary) <- 
  names(stoma_pop_minmax) <- 
  stoma_vars

# Pivot both summaries to wider format
stoma_pop_summary <- melt(stoma_pop_summary, 
                            measure.vars = c("N", "mean", "sd", "se", "ci")) %>% 
  pivot_wider(names_from = "variable",
              values_from = "value") %>% 
  rename(variable = L1)

stoma_pop_minmax <- melt(stoma_pop_minmax, 
                           measure.vars = c("min", "max")) %>% 
  pivot_wider(names_from = "variable",
              values_from = "value") %>% 
  rename(variable = L1)


# Combine all summaries in a single object
all_summary <- bind_rows(morph_pop_summary,
                         stoma_pop_summary,
                         pollen_pop_summary) %>% 
  select(variable, Population_ID, N, mean, sd, se, ci) %>% 
  arrange(factor(variable, levels = c(morph_vars,
                                      pollen_vars,
                                      stoma_vars)))

# Combine all minmax in a single object
all_minmax <- bind_rows(morph_pop_minmax,
                        stoma_pop_minmax,
                        pollen_pop_minmax)  %>% 
  arrange(factor(variable, levels = c(morph_vars,
                                      pollen_vars,
                                      stoma_vars)))

# Combine global summaries and minmax in a single object
all_summary <- all_minmax %>% 
  left_join(all_summary) %>% 
  left_join(genetic_groups) %>% 
  select(variable, Population_ID, Population_name, genetic_group_k5, N, min, max, mean, sd, se, ci) %>% 
  rename(pop_id = Population_ID, 
         pop_name = Population_name,
         genetic_group = genetic_group_k5) %>% 
  mutate(pop_name = str_to_title(pop_name)) %>% 
  arrange(factor(variable, levels = c(morph_vars,
                                      pollen_vars,
                                      stoma_vars)),
          genetic_group,
          pop_id)


# ----------------------------------------------------------------------
# EXPORT RESULTS
# ----------------------------------------------------------------------

all_summary %>% 
  mutate(across(min:ci, \(x)round(x, digits = 2))) %>% 
  openxlsx::write.xlsx("outputs/tables/Table_S3_populations.xlsx")


all_summary %>%
  pivot_longer(cols = c(min, max)) %>% 
  mutate(name = factor(name, levels = c("min", "max"))) %>% 
  ggplot(aes(x = genetic_group,
             y = value)) +
    geom_boxplot() +
  facet_grid(rows = vars(variable),
             cols = vars(name),
             scale = "free")
  
# Export figure files in publication-ready formats.
ggsave("outputs/figures/temp_morph_populations.pdf",
       width = 5,
       height = 10)
