# ======================================================================
# Script: 01_Main_morphology.R
# Purpose: Summarise morphological, pollen and stomatal traits by genetic group and ploidy level; export summary tables for the manuscript.
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
require(tidyverse)
require(reshape2)
require(FSA)
require(factoextra)
library(ggpubr)


# ----------------------------------------------------------------------
# LOAD DATA
# ----------------------------------------------------------------------
#### LOAD GENETIC GROUPS DATA
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


#### JOIN GENETIC GROUPS AND MORPHOMETRIC DATA
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


#### LOAD POLLEN DATA
pollen_measures <- openxlsx::read.xlsx("data/pollen.xlsx", 2)

# Inspect object structure and summary statistics.
head(pollen_measures)  
summary(pollen_measures)


#### LOAD STOMA DATA
stoma_measures <- openxlsx::read.xlsx("data/stoma.xlsx", 2)

# Inspect object structure and summary statistics.
head(stoma_measures)  
summary(stoma_measures)


#### JOIN POLLEN AND STOMA DATA WITH GENETIC GROUP DATA
pollen_measures <- pollen_measures %>% 
  left_join(genetic_groups) %>% 
  rename(clade_5 = genetic_group_k5) %>% 
  rename(clade_2 = genetic_group_k2)

stoma_measures <- stoma_measures %>% 
  left_join(genetic_groups) %>% 
  rename(clade_5 = genetic_group_k5) %>% 
  rename(clade_2 = genetic_group_k2)


# ----------------------------------------------------------------------
# ANALYSE VARIABLES
# ----------------------------------------------------------------------
#### ANALYSE VARIABLES ####

## morphology ##

# set morphological variables' names
morph_vars <- c("co",
                "ca", 
                "ltooth",
                "stooth", 
                "hair",
                "fl_infl",
                "d_infl",
                "prop_ltooth"
)


# Create morph summary by genetic group
morph_clade_summary <- lapply(morph_vars,
                              FUN = \(x)Rmisc::summarySE(morph_measures, 
                                                         measurevar = x, 
                                                         groupvars = c("clade_5"), 
                                                         na.rm = TRUE)) %>% 
  lapply(FUN = \(x){colnames(x)[3] <- "mean"; x})


# Create morph summary by ploidy level
morph_ploidy_summary <- lapply(morph_vars, 
                               FUN = \(x)Rmisc::summarySE(morph_measures, 
                                                          measurevar = x, 
                                                          groupvars = c("clade_2"), 
                                                          na.rm = TRUE)) %>% 
  lapply(FUN = \(x){colnames(x)[3] <- "mean"; x})


# Create morph minmax by genetic group
morph_clade_minmax <- lapply(morph_vars,
                             FUN = \(x)summarise(group_by(morph_measures, 
                                                          clade_5),
                                                 min = min(get(x),
                                                           na.rm = TRUE),
                                                 max = max(get(x),
                                                           na.rm = TRUE)))


# Create morph minmax by ploidy level
morph_ploidy_minmax <- lapply(morph_vars, 
                              FUN = \(x)summarise(group_by(morph_measures, 
                                                           clade_2),
                                                  min = min(get(x),
                                                            na.rm = TRUE),
                                                  max = max(get(x),
                                                            na.rm = TRUE)))


# Set names in both summaries
names(morph_clade_summary) <-
  names(morph_ploidy_summary) <- 
  names(morph_clade_minmax) <- 
  names(morph_ploidy_minmax) <- 
  morph_vars


# Pivot all summaries to wider format
morph_clade_summary <- melt(morph_clade_summary, 
                     measure.vars = c("N", "mean", "sd", "se", "ci")) %>% 
  pivot_wider(names_from = "variable",
              values_from = "value") %>% 
  rename(clade = clade_5) %>% 
  rename(variable = L1)

morph_ploidy_summary <- melt(morph_ploidy_summary, 
                            measure.vars = c("N", "mean", "sd", "se", "ci")) %>% 
  pivot_wider(names_from = "variable",
              values_from = "value") %>% 
  rename(clade = clade_2) %>% 
  rename(variable = L1)

morph_clade_minmax <- melt(morph_clade_minmax, 
                            measure.vars = c("min", "max")) %>% 
  pivot_wider(names_from = "variable",
              values_from = "value") %>% 
  rename(clade = clade_5) %>% 
  rename(variable = L1)

morph_ploidy_minmax <- melt(morph_ploidy_minmax, 
                           measure.vars = c("min", "max")) %>% 
  pivot_wider(names_from = "variable",
              values_from = "value") %>% 
  rename(clade = clade_2) %>% 
  rename(variable = L1)



## pollen ##

# Set names of pollen variables
pollen_vars <- c("d_max", "d_min")


# Create pollen summary by genetic group
pollen_clade_summary <- lapply(pollen_vars, 
                               FUN = \(x)Rmisc::summarySE(pollen_measures, 
                                                          measurevar = x, 
                                                          groupvars = c("clade_5"), 
                                                          na.rm = TRUE)) %>% 
  lapply(FUN = \(x){colnames(x)[3] <- "mean"; x})


# Create pollen summary by ploidy

pollen_ploidy_summary <- lapply(pollen_vars, 
                                FUN = \(x)Rmisc::summarySE(pollen_measures, 
                                                           measurevar = x, 
                                                           groupvars = c("clade_2"), 
                                                           na.rm = TRUE)) %>% 
  lapply(FUN = \(x){colnames(x)[3] <- "mean"; x})


# Create morph minmax by genetic group
pollen_clade_minmax <- lapply(pollen_vars, 
                              FUN = \(x)summarise(group_by(pollen_measures,
                                                           clade_5),
                                                  min = min(get(x), 
                                                            na.rm = TRUE),
                                                  max = max(get(x),
                                                            na.rm = TRUE)))


# Create morph minmax by ploidy level
pollen_ploidy_minmax <- lapply(pollen_vars,
                               FUN = \(x)summarise(group_by(pollen_measures,
                                                            clade_2),
                                                   min = min(get(x), 
                                                             na.rm = TRUE),
                                                   max = max(get(x),
                                                             na.rm = TRUE)))


# Set names in all summaries
names(pollen_clade_summary) <- 
  names(pollen_ploidy_summary) <- 
  names(pollen_clade_minmax) <- 
  names(pollen_ploidy_minmax) <- 
  pollen_vars


# Pivot all summaries to wider format
pollen_clade_summary <- melt(pollen_clade_summary, 
                      measure.vars = c("N", "mean", "sd", "se", "ci")) %>% 
  pivot_wider(names_from = "variable",
              values_from = "value") %>% 
  rename(clade = clade_5) %>% 
  rename(variable = L1)

pollen_ploidy_summary <- melt(pollen_ploidy_summary, 
                             measure.vars = c("N", "mean", "sd", "se", "ci")) %>% 
  pivot_wider(names_from = "variable",
              values_from = "value") %>% 
  rename(clade = clade_2) %>% 
  rename(variable = L1)

pollen_clade_minmax <- melt(pollen_clade_minmax, 
                           measure.vars = c("min", "max")) %>% 
  pivot_wider(names_from = "variable",
              values_from = "value") %>% 
  rename(clade = clade_5) %>% 
  rename(variable = L1)

pollen_ploidy_minmax <- melt(pollen_ploidy_minmax, 
                            measure.vars = c("min", "max")) %>% 
  pivot_wider(names_from = "variable",
              values_from = "value") %>% 
  rename(clade = clade_2) %>% 
  rename(variable = L1)



## stoma ##

# Set names for stoma variables
stoma_vars <- c("stoma_w", "stoma_l")


# Create stoma summary by genetic group
stoma_clade_summary <- lapply(stoma_vars, 
                              FUN = \(x)Rmisc::summarySE(stoma_measures, 
                                                         measurevar = x, 
                                                         groupvars = c("clade_5"), 
                                                         na.rm = TRUE)) %>% 
  lapply(FUN = \(x){colnames(x)[3] <- "mean"; x})

# Create stoma summary by ploidy
stoma_ploidy_summary <- lapply(stoma_vars,
                               FUN = \(x)Rmisc::summarySE(stoma_measures, 
                                                          measurevar = x, 
                                                          groupvars = c("clade_2"), 
                                                          na.rm = TRUE)) %>% 
  lapply(FUN = \(x){colnames(x)[3] <- "mean"; x})

# Create stoma minmax by genetic group
stoma_clade_minmax <- lapply(stoma_vars, 
                             FUN = \(x)summarise(group_by(stoma_measures, 
                                                          clade_5),
                                                 min = min(get(x), 
                                                           na.rm = TRUE),
                                                 max = max(get(x),
                                                           na.rm = TRUE)))

# Create stoma minmax by ploidy level
stoma_ploidy_minmax <- lapply(stoma_vars, FUN = \(x)summarise(group_by(stoma_measures, clade_2),
                                                                min = min(get(x), na.rm = TRUE),
                                                                max = max(get(x), na.rm = TRUE)))

# Set names in all summaries
names(stoma_clade_summary) <- 
  names(stoma_ploidy_summary) <-
  names(stoma_clade_minmax) <- 
  names(stoma_ploidy_minmax) <-
  stoma_vars


# Pivot all summaries to wider format
stoma_clade_summary <- melt(stoma_clade_summary, 
                             measure.vars = c("N", "mean", "sd", "se", "ci")) %>% 
  pivot_wider(names_from = "variable",
              values_from = "value") %>% 
  rename(clade = clade_5) %>% 
  rename(variable = L1)

stoma_ploidy_summary <- melt(stoma_ploidy_summary, 
                            measure.vars = c("N", "mean", "sd", "se", "ci")) %>% 
  pivot_wider(names_from = "variable",
              values_from = "value") %>% 
  rename(clade = clade_2) %>% 
  rename(variable = L1)

stoma_clade_minmax <- melt(stoma_clade_minmax, 
                            measure.vars = c("min", "max")) %>% 
  pivot_wider(names_from = "variable",
              values_from = "value") %>% 
  rename(clade = clade_5) %>% 
  rename(variable = L1)

stoma_ploidy_minmax <- melt(stoma_ploidy_minmax, 
                             measure.vars = c("min", "max")) %>% 
  pivot_wider(names_from = "variable",
              values_from = "value") %>% 
  rename(clade = clade_2) %>% 
  rename(variable = L1)



# Combine all summaries in a single object
all_summary <- bind_rows(morph_ploidy_summary, 
                         stoma_ploidy_summary, 
                         pollen_ploidy_summary,
                         morph_clade_summary,
                         stoma_clade_summary,
                         pollen_clade_summary) %>% 
  select(variable, clade, N, mean, sd, se, ci) %>% 
  arrange(factor(variable, levels = c(morph_vars,
                                pollen_vars,
                                stoma_vars)),
                 factor(clade, levels = c("4x", 
                                          "2x",
                                          "Tetraploid",
                                          "Hercynian",
                                          "Algarve",
                                          "Doñana",
                                          "Cádiz")))

# Combine all minmax in a single object
all_minmax <- bind_rows(morph_ploidy_minmax, 
                        stoma_ploidy_minmax, 
                        pollen_ploidy_minmax,
                        morph_clade_minmax,
                        stoma_clade_minmax,
                        pollen_clade_minmax)  %>% 
  arrange(factor(variable, levels = c(morph_vars,
                                      pollen_vars,
                                      stoma_vars)),
          factor(clade, levels = c("4x", 
                                   "2x",
                                   "Tetraploid",
                                   "Hercynian",
                                   "Algarve",
                                   "Doñana",
                                   "Cádiz")))

# Left join global minmax and summary objects in a single object 
all_summary <- all_minmax %>% left_join(all_summary) %>% 
  select(variable, clade, N, min, max, mean, sd, se, ci)



# ----------------------------------------------------------------------
# EXPORT RESULTS
# ----------------------------------------------------------------------
# Write table to excel 
all_summary %>% 
  mutate(across(min:ci, \(x)round(x, digits = 2))) %>% 
  openxlsx::write.xlsx("outputs/tables/Table_S2.xlsx")

# Remove temporary objects
rm(morph_ploidy_summary, 
   stoma_ploidy_summary, 
   pollen_ploidy_summary,
   morph_clade_summary,
   stoma_clade_summary,
   pollen_clade_summary,
   all_summary)

rm(morph_ploidy_minmax, 
   stoma_ploidy_minmax, 
   pollen_ploidy_minmax,
   morph_clade_minmax,
   stoma_clade_minmax,
   pollen_clade_minmax,
   all_minmax)


# ----------------------------------------------------------------------
# LM MODELS
# ----------------------------------------------------------------------
#### LM MODELS ####
# Test for differences with corolla length and test normality and homocedasticity
hist(morph_measures$co) # Not normal
shapiro.test(morph_measures$co) # Not normal

hist(log(morph_measures$co)) # Log transformation improve normality, but not entirely
shapiro.test(log(morph_measures$co)) # Not normal. We continue just to double check on the residuals.

aov_co <- aov(log(co) ~ clade_5 * Sex, data = morph_measures)
anova_co <- drop1(aov_co, ~., test = "F") ##anova with typeIII error
anova_co # There are differences between clades and sexes. Hermaphrodites are always bigger within the same clade
tukey <- TukeyHSD(aov(log(morph_measures$co) ~ morph_measures$clade_5 * morph_measures$Sex))
tukey

print(tukey.cld_co)

hist(residuals(aov_co))
shapiro.test(residuals(aov_co)) # Residuals are not normal
lmtest::bptest(aov_co) # Residuals are not homoskedastic

rm(aov_co,
   tukey)


# ----------------------------------------------------------------------
# ANALYSE DATA WITH NO PARAMETRIC APPROACHES
# ----------------------------------------------------------------------

#Define function to do the comparisons with Dunn test.
extract_differences_from_dunnTest <- function(var, group_var, data){
  form <- as.formula(paste0(var, " ~ ", group_var))
  post_hoc_clade <- dunnTest(form,
                                data = data,
                                method = "bonferroni") 
  
  post_hoc_clade_res <- post_hoc_clade$res
  
  if(nrow(post_hoc_clade_res) > 1){
    clade_cld <- rcompanion::cldList(comparison = post_hoc_clade_res$Comparison,
                                     p.value = post_hoc_clade_res$P.adj,
                                     threshold = 0.05)[1:2]
    clade_cld
  } else {
    post_hoc_clade_res
  }
}

# Run the comparisons for morphology across genetic groups and ploidy levels
morph_clade_diff_letters <- lapply(morph_vars, FUN = extract_differences_from_dunnTest, "clade_5", morph_measures)
morph_ploidy_diff_letters <- lapply(morph_vars, FUN = extract_differences_from_dunnTest, "clade_2", morph_measures)
names(morph_clade_diff_letters) <- names(morph_ploidy_diff_letters) <- morph_vars

# Run the comparisons for pollen across genetic groups and ploidy levels
pollen_clade_diff_letters <- lapply(pollen_vars, FUN = extract_differences_from_dunnTest, "clade_5", pollen_measures)
pollen_ploidy_diff_letters <- lapply(pollen_vars, FUN = extract_differences_from_dunnTest, "clade_2", pollen_measures)
names(pollen_clade_diff_letters) <- names(pollen_ploidy_diff_letters) <- pollen_vars

# Run the comparisons for stoma across genetic groups and ploidy levels
stoma_clade_diff_letters <- lapply(stoma_vars, FUN = extract_differences_from_dunnTest, "clade_5", stoma_measures)
stoma_ploidy_diff_letters <- lapply(stoma_vars, FUN = extract_differences_from_dunnTest, "clade_2", stoma_measures)
names(stoma_clade_diff_letters) <- names(stoma_ploidy_diff_letters) <- stoma_vars

# Combine results for genetic groups
clade_diff_letters <- list(morph_clade_diff_letters, 
     pollen_clade_diff_letters, 
     stoma_clade_diff_letters) %>% 
  melt(id.vars =c("Group", "Letter")) %>% 
  select(Group, Letter, L2) %>% 
  rename(name = L2,
         clade_5 = Group,
         label = Letter)

# Combine results for ploidy levels
ploidy_diff_letters <- list(morph_ploidy_diff_letters,
                     pollen_ploidy_diff_letters,
                     stoma_ploidy_diff_letters) %>% 
  melt(id.vars =c("Comparison", "Z", "P.unadj", "P.adj")) %>% 
  select(Comparison, Z, P.unadj, P.adj, L2) %>% 
  rename(name = L2)

# ----------------------------------------------------------------------
# PLOTS
# ----------------------------------------------------------------------

# Define colors for genetic groups in all plots
clade_colors <- c("#009E73", "#F0E442", "#CC79A7", "#E69F00", "#6388b4")

# Transform data in longer format
morph_lf <- morph_measures %>%
  select(Population_ID, Individual_ID, clade_2, clade_5, all_of(morph_vars)) %>% 
  pivot_longer(cols = all_of(morph_vars))

pollen_lf <- pollen_measures %>%
  select(Population_ID, Individual_ID, clade_2, clade_5, all_of(pollen_vars)) %>% 
  pivot_longer(cols = all_of(pollen_vars))

stoma_lf <- stoma_measures %>%
  select(Population_ID, Individual_ID, clade_2, clade_5, all_of(stoma_vars)) %>% 
  pivot_longer(cols = all_of(stoma_vars))

# Combine all objects in a single object
data_lf <- bind_rows(morph_lf,
                     pollen_lf,
                     stoma_lf) %>% 
  mutate(name = factor(name, levels = c(morph_vars,
                                        stoma_vars,
                                        pollen_vars)))
  
clade_diff_letters <- data_lf %>% 
  group_by(name, clade_5) %>% 
  summarise(max = max(value, na.rm = TRUE) + 
              (0.5 * max(value, na.rm = TRUE))) %>% 
  right_join(clade_diff_letters) %>% 
  mutate(name = factor(name, levels = c(morph_vars,
                                        stoma_vars,
                                        pollen_vars)))

ploidy_diff_letters <- data_lf %>%
  group_by(name) %>% 
  summarise(max = max(value, na.rm = TRUE) + 
              (0.25 * max(value, na.rm = TRUE))) %>% 
  left_join(ploidy_diff_letters) %>% 
  mutate(name = factor(name, 
                       levels = c(morph_vars,
                                  stoma_vars,
                                  pollen_vars))) %>% 
  mutate(label = case_when(
    P.adj > 0.05 ~ "ns",
    P.adj <= 0.05 & P.adj > 0.01 ~ "*",
    P.adj <= 0.01 & P.adj > 0.001 ~ "**",
    P.adj <= 0.001 ~ "***")) %>% 
  mutate(clade_2 = "Diploid")



#### CLADE PLOT ####

# Define variables for main figure...
main_vars <- c("co", "ca", "ltooth", "hair", "fl_infl", "d_infl", "prop_ltooth")

# ... and for the supplementary figures
suppl_vars <- c("stoma_w", "stoma_l", "d_max", "d_min")

# Construct the plot object.
measures_ploidy_plot <- ggplot(data_lf %>% 
                                 filter(name %in% main_vars), 
                               aes(x = name, 
                                   y = value, 
                                   fill = clade_2)) +
  introdataviz::geom_split_violin(alpha = .4, 
                                  trim = FALSE) + 
  geom_boxplot(width = 0.15, 
               alpha = 0.6, 
               show.legend = FALSE) + 
  stat_summary(fun = mean, 
               geom = "point", 
               color= "black", 
               position = position_dodge(.15)) + 
  scale_y_continuous(name = "Length",
                     expand = expansion(mult = c(0.05, 
                                                 0.15))) +  
  scale_fill_manual(values = c("2x" = "#ffae34",
                               "4x" = "#6388b4")) +  
  labs(fill = "Ploidy") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  geom_text(data = ploidy_diff_letters %>% 
              filter(name %in% main_vars),
            aes(label = label),
            y = Inf,
            vjust = 1.4,
            size = 4) +
  facet_wrap(facets = vars(name), 
             ncol = 4,
             scales = "free") 

# Move legend to empty cell in the grid
measures_ploidy_plot <- measures_ploidy_plot %>%
  lemon::reposition_legend(x = 0.5,
                           y = 0.7,
                           just = c("center", "center"),
                           panel = "panel-4-2")

measures_ploidy_plot

# Export figure files in publication-ready formats.
ggsave("outputs/figures/Figure_2.pdf",
       measures_ploidy_plot,
       width = 8,
       height = 4)
ggsave("outputs/figures/Figure_2.png",
       measures_ploidy_plot,
       width = 8,
       height = 4)


# Construct the plot object.
measures_ploidy_plot_suppl <- ggplot(data_lf %>% 
                                 filter(name %in% suppl_vars), 
                               aes(x = name, 
                                   y = value, 
                                   fill = clade_2)) +
  introdataviz::geom_split_violin(alpha = .4, 
                                  trim = FALSE) +  
  geom_boxplot(width = 0.15, 
               alpha = 0.6, 
               show.legend = FALSE) +
  stat_summary(fun = mean, 
               geom = "point", 
               color= "black", 
               position = position_dodge(.15)) +  
  scale_y_continuous(name = "Length",
                     expand = expansion(mult = c(0.05, 
                                                 0.15))) +
  scale_fill_manual(values = c("2x" = "#ffae34",
                               "4x" = "#6388b4")) +  
  labs(fill = "Ploidy") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  geom_text(data = ploidy_diff_letters %>% 
              filter(name %in% suppl_vars),
            aes(label = label),
            y = Inf,
            vjust = 1.4,
            size = 4) +
  facet_wrap(facets = vars(name), 
             ncol = 4,
             scales = "free") 


measures_ploidy_plot_suppl


#### CLADE PLOT ####

# Construct the plot object.
measures_clade_plot <- ggplot(data_lf %>%
                                filter(name %in% main_vars),
                              aes(x = clade_5,
                                  y = value, 
                                  fill = clade_5)) +
  geom_violin(alpha = .4, 
              trim = FALSE) +
  geom_boxplot(width = 0.15, 
               alpha = 0.6, 
               show.legend = FALSE) +
  scale_fill_manual(values = clade_colors) + 
  stat_summary(fun.data = "mean_sdl", 
               fun.args = list(mult = 1),
               geom = "point", 
               size = 2, 
               position = position_dodge(0.9)) +
  labs(fill = "Group") +
  theme_bw() +
  scale_y_continuous(name = "Length",
                     expand = expansion(mult = c(0.05, 
                                                 0.15))) +  
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 1, 
                                   hjust=1),
        axis.title.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  geom_text(data = clade_diff_letters  %>%
              filter(name %in% main_vars) %>% 
              filter(clade_5 %in% c("Algarve",
                                  "Cádiz",
                                  "Doñana",
                                  "Hercynian",
                                  "Tetraploid")),
            aes(label = label),
            y = Inf,
            vjust = 1.2) +
  facet_wrap(facets = vars(name), 
             ncol = 4,
             scales = "free_y")

# Move legend to empty cell in the grid
measures_clade_plot <- measures_clade_plot %>%
  lemon::reposition_legend(x = 0.5,
                           y = 0.25,
                           just = c("center", "center"),
                           panel = "panel-4-2")

measures_clade_plot

# Export figure files in publication-ready formats.
ggsave("outputs/figures/Figure_3.pdf",
       measures_clade_plot,
       width = 8,
       height = 4)
ggsave("outputs/figures/Figure_3.png",
       measures_clade_plot,
       width = 8,
       height = 4)



# Construct the plot object.
measures_clade_plot_suppl_b <- ggplot(data_lf %>%
                                filter(name %in% suppl_vars),
                              aes(x = clade_5,
                                  y = value, 
                                  fill = clade_5)) +
  geom_violin(alpha = .4, 
              trim = FALSE) +
  geom_boxplot(width = 0.15, 
               alpha = 0.6, 
               show.legend = FALSE) +
  scale_fill_manual(values = clade_colors) + 
  stat_summary(fun.data = "mean_sdl", 
               fun.args = list(mult = 1),
               geom = "point", 
               size = 2, 
               position = position_dodge(0.9)) +
  labs(fill = "Group") +
  theme_bw() +
  scale_y_continuous(name = "Length",
                     expand = expansion(mult = c(0.05, 
                                                 0.15))) +  # Etiqueta eje Y
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 1, 
                                   hjust=1),
        axis.title.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  geom_text(data = clade_diff_letters  %>%
              filter(name %in% suppl_vars) %>% 
              filter(clade_5 %in% c("Algarve",
                                    "Cádiz",
                                    "Doñana",
                                    "Hercynian",
                                    "Tetraploid")),
            aes(label = label),
            y = Inf,
            vjust = 1.2) +
  facet_wrap(facets = vars(name), 
             ncol = 4,
             scales = "free_y")

measures_clade_plot_suppl_b

# Combine the two supplementary figures in a single plot
measures_plot_suppl <- ggpubr::ggarrange(measures_ploidy_plot_suppl, 
                                         measures_clade_plot_suppl_b,
                                         ncol = 1, 
                                         align = "v",
                                         heights = c(1, 1.3))

measures_plot_suppl

# Export figure files in publication-ready formats.
ggsave("outputs/figures/Figure_S1.pdf",
       measures_plot_suppl,
       width = 9,
       height = 5)
ggsave("outputs/figures/Figure_S1.png",
       measures_plot_suppl,
       width = 9,
       height = 5)


# ----------------------------------------------------------------------
# SEX EFFECT
# ----------------------------------------------------------------------
#### Sex effect on corolla comparisons ####

sex_diff_letters <- morph_measures %>% 
  dunnTest(co ~ Sex, data = ., method = "bonferroni") %>% 
  .[["res"]] %>% 
  mutate(label = case_when(
    P.adj > 0.05 ~ "ns",
    P.adj <= 0.05 & P.adj > 0.01 ~ "*",
    P.adj <= 0.01 & P.adj > 0.001 ~ "**",
    P.adj <= 0.001 ~ "***")) %>% 
  mutate(Sex = "")

# Construct the plot object.
plot_sex_1 <- morph_measures %>% ggplot(aes(x = Sex,
                         y = co, 
                         fill = Sex)) +
    geom_violin(alpha = .4, trim = FALSE) +
    geom_boxplot(width = 0.15, alpha = 0.6, show.legend = FALSE) +
    scale_fill_manual(values = c("coral", "turquoise")) + 
    stat_summary(fun.data = "mean_sdl", 
                 fun.args = list(mult = 1),
                 geom = "point", 
                 size = 2, 
                 position = position_dodge(0.9)) +
# Add layers and styling to the plot.
  scale_y_continuous(name = "Corolla",
                     expand = expansion(mult = c(0.05, 
                                                 0.15))) +  # Etiqueta eje Y
# Add layers and styling to the plot.
  theme_bw() +
    theme(axis.title.x = element_blank()) +
    annotate("text",
             x = 1.5,
             y = Inf,
             vjust = 3,
             label = sex_diff_letters$label)

plot_sex_1

# Wrangle and prepare data for downstream analysis.
morph_clade_F_diff_letters <- lapply("co", FUN = extract_differences_from_dunnTest, "clade_5", morph_measures %>% filter(Sex == "F"))
morph_ploidy_F_diff_letters <- lapply("co", FUN = extract_differences_from_dunnTest, "clade_2", morph_measures %>% filter(Sex == "F"))
morph_clade_H_diff_letters <- lapply("co", FUN = extract_differences_from_dunnTest, "clade_5", morph_measures %>% filter(Sex == "H"))
morph_ploidy_H_diff_letters <- lapply("co", FUN = extract_differences_from_dunnTest, "clade_2", morph_measures %>% filter(Sex == "H"))

names(morph_clade_F_diff_letters) <- names(morph_ploidy_F_diff_letters) <- "F"
names(morph_clade_H_diff_letters) <- names(morph_ploidy_H_diff_letters) <- "H"

clade_sex_diff_letters <- list(morph_clade_F_diff_letters, 
                           morph_clade_H_diff_letters) %>% 
  melt(id.vars =c("Group", "Letter")) %>% 
  select(Group, Letter, L2) %>% 
  rename(Sex = L2,
         clade_5 = Group,
         label = Letter)

ploidy_sex_diff_letters <- list(morph_ploidy_F_diff_letters,
                                morph_ploidy_H_diff_letters) %>% 
  melt(id.vars =c("Comparison", "Z", "P.unadj", "P.adj")) %>% 
  select(Comparison, Z, P.unadj, P.adj, L2) %>% 
  rename(Sex = L2) %>% 
  mutate(clade_2 = "2x") %>% 
  mutate(label = case_when(
    P.adj > 0.05 ~ "ns",
    P.adj <= 0.05 & P.adj > 0.01 ~ "*",
    P.adj <= 0.01 & P.adj > 0.001 ~ "**",
    P.adj <= 0.001 ~ "***"))


# Construct the plot object.
plot_sex_2 <- ggplot(morph_measures, aes(x = "co", 
                    y = co, 
                    fill = clade_2)) +
  introdataviz::geom_split_violin(alpha = .4, 
                                  trim = FALSE) +
  geom_boxplot(width = 0.15, 
               alpha = 0.6, 
               show.legend = FALSE) +
  stat_summary(fun = mean, 
               geom = "point", 
               color= "black", 
               position = position_dodge(.15)) +  
  scale_y_continuous(name = "Corolla",
                     expand = expansion(mult = c(0.05, 
                                                 0.15))) +  
  scale_fill_manual(values = c("2x" = "#ffae34",
                               "4x" = "#6388b4")) +  
  labs(fill = "Ploidy") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  geom_text(data = ploidy_sex_diff_letters,
            aes(label = label),
            y = Inf,
            vjust = 1.2) +
  facet_wrap(facets = vars(Sex), 
             ncol = 4)

plot_sex_2

# Construct the plot object.
plot_sex_3 <- ggplot(morph_measures, aes(x = clade_5, 
                           y = co, 
                           fill = clade_5)) +
  geom_violin(alpha = .4, trim = FALSE) +
  geom_boxplot(width = 0.15, 
               alpha = 0.6, 
               show.legend = FALSE) +
  stat_summary(fun.data = "mean_sdl", 
               fun.args = list(mult = 1),
               geom = "point", 
               size = 2, 
               position = position_dodge(0.9)) +
  scale_y_continuous(name = "Corolla",
                     expand = expansion(mult = c(0.05, 
                                                 0.15))) +
  scale_fill_manual(values = clade_colors) + 
  labs(fill = "Group") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 1, 
                                   hjust=1),
        axis.title.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  geom_text(data = clade_sex_diff_letters %>% 
              filter(clade_5 %in% c("Algarve",
                                    "Cádiz",
                                    "Doñana",
                                    "Hercynian",
                                    "Tetraploid")),
            aes(label = label),
            y = Inf,
            vjust = 1.1) +
  facet_wrap(facets = vars(Sex), 
             ncol = 4)

plot_sex_3

plot_sex <- ggarrange(ggarrange(plot_sex_1, 
                    plot_sex_2, 
                    ncol = 2, 
                    labels = c("A", "B")),
          plot_sex_3, 
          nrow = 2,
          labels = c("","C"))

plot_sex

# Export figure files in publication-ready formats.
ggsave("outputs/figures/Figure_S2.pdf",
       plot_sex,
       width = 8,
       height = 6.5)
ggsave("outputs/figures/Figure_S2.png",
       plot_sex,
       width = 8,
       height = 6.5)


# ----------------------------------------------------------------------
# PCA ANALYSIS
# ----------------------------------------------------------------------

#### genetic groups ####

# Run PCA
pca <- prcomp(~ co + ca + ltooth + hair + d_infl + prop_ltooth,
              data = morph_measures, 
              scale. = TRUE, 
              na.action = na.omit)

df.pca <- pca$x

data_pca <- morph_measures %>% 
  select(co, ca, ltooth, hair, d_infl, prop_ltooth, clade_2, clade_5) %>% 
  na.omit()

rownames(pca$rotation) <- c("Corolla", 
                            "Calyx",
                            "Longest calyx tooth", 
                            "Calyx hair length", 
                            "Inflorescence diameter", 
                            "Proportion of longest calyx tooth")

pca_ploidy_plot <- fviz_pca_biplot(pca, 
                          label = "var", 
                          habillage = data_pca$clade_2,
                          col.var = "black",
                          repel = TRUE, 
                          geom.ind = "point", 
                          col.ind = data_pca$clade_2,     
                          fill.ind = data_pca$clade_2,   
                          shape.ind = data_pca$clade_2,  
                          alpha.ind = 0.3,
                          addEllipses = TRUE,
                          ellipse.level=0.95, 
                          ellipse.alpha = 0.05,
                          ellipse.type = "norm",
                          ellipse.fill = data_pca$clade_2,
                          title = "",
                          mean.point.size = 5) +
  scale_shape_manual(values = c("2x" = 22, "4x" = 21)) +
  scale_fill_manual(values = c("2x" = "#ffae34", "4x" = "#6388b4")) +
  scale_color_manual(values = c("2x" = "#ffae34", "4x" = "#6388b4")) +  # para el borde y la elipse
  theme_bw() + 
  labs(x = "Principal Component 1 (69.2%)", 
       y = "Principal Component 2 (12.4%)") +
  guides(fill = guide_legend(title="Ploidy"),
         color = guide_legend(title="Ploidy"),
         shape = guide_legend(title="Ploidy"))

pca_ploidy_plot

# Export figure files in publication-ready formats.
ggsave("outputs/figures/Figure_S3.pdf", 
       pca_ploidy_plot,
       width = 6,
       height = 6)
ggsave("outputs/figures/Figure_S3.png", 
       pca_ploidy_plot,
       width = 6,
       height = 6)  



pca_clades_plot <- fviz_pca_biplot(pca, 
                                   label = "var", 
                                   habillage = data_pca$clade_5,
                                   col.var = "black",
                                   repel = TRUE, 
                                   geom.ind = "point", 
                                   col.ind = data_pca$clade_5,    
                                   fill.ind = data_pca$clade_5,   
                                   shape.ind = data_pca$clade_5,  
                                   alpha.ind = 0.3,
                                   addEllipses = TRUE,
                                   ellipse.level=0.95, 
                                   ellipse.alpha = 0.05,
                                   ellipse.type = "norm",
                                   ellipse.fill = data_pca$clade_5,
                                   title = "",
                                   mean.point.size = 5) +
  scale_shape_manual(values = c("Tetraploid" = 21,
                                "Hercynian" = 22,
                                "Algarve" = 22,
                                "Doñana" = 22,
                                "Cádiz" = 22)) + 
  scale_fill_manual(values = c("Tetraploid" = "#6388b4", 
                               "Hercynian" = "#E69F00", 
                               "Algarve" = "#009E73",
                               "Doñana" = "#CC79A7", 
                               "Cádiz" = "#F0E442")) + 
  scale_color_manual(values = c("Tetraploid" = "#6388b4",
                                "Hercynian" = "#E69F00", 
                                "Algarve" = "#009E73",
                                "Doñana" = "#CC79A7", 
                                "Cádiz" = "#F0E442")) + 
  labs(x = "Principal Component 1 (69.2%)", 
       y = "Principal Component 2 (12.4%)") +
  guides(fill = guide_legend(title="Group"),
         color = guide_legend(title="Group"),
         shape = guide_legend(title="Group"))

pca_clades_plot


# Export figure files in publication-ready formats.
ggsave("outputs/figures/Figure_S4.pdf", 
       pca_clades_plot,
       width = 7,
       height = 6)
ggsave("outputs/figures/Figure_S4.png", 
       pca_clades_plot,
       width = 7,
       height = 6) 
