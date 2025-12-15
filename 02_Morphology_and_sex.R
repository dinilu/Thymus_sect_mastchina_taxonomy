# ======================================================================
# Script: 02_Morphology_and_sex.R
# Purpose: Assess sex-related differences in morphological traits within clades and ploidy levels; produce summary tables and non-parametric post hoc tests.
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
#### LOAD GENETIC GROUP DATA
genetic_groups <- "data/genetic groups.xlsx" %>% 
  openxlsx::read.xlsx(2) %>% 
  mutate(genetic_group_k2 = case_match(genetic_group_k2,
             "tetraploid" ~ "4x",
             "diploid" ~ "2x")) %>% 
  mutate(genetic_group_k5 = str_to_sentence(genetic_group_k5))



#### LOAD MORPHOMETRIC DATA  ####
morph_measures <- "data/Morphometric_measures.xlsx" %>% 
  openxlsx::read.xlsx(2)

# Inspect object structure and summary statistics.
head(morph_measures)  
# Inspect object structure and summary statistics.
summary(morph_measures)

#### COMBINE GENETIC GROUP AND MORPHOMETRIC DATA  ####
morph_measures <- morph_measures %>% 
  left_join(genetic_groups) %>% 
  rename(clade_5 = genetic_group_k5) %>% 
  rename(clade_2 = genetic_group_k2) %>% 
  mutate(across(c(Population_ID, # Change to factors.
                  Sex, 
                  clade_2,
                  clade_5), 
                as.factor))



# ----------------------------------------------------------------------
# ANALYSE VARIABLES
# ----------------------------------------------------------------------
#### ANALYSE VARIABLES ####

## morph ##

morph_vars <- c("co",
                "ca", 
                "ltooth",
                "stooth", 
                "hair",
                "fl_infl",
                "d_infl",
                "prop_ltooth"
)


# Create morph summary by sex and genetic group
morph_clade_summary <- lapply(morph_vars, 
                              FUN = \(x)Rmisc::summarySE(morph_measures, 
                                                         measurevar = x, 
                                                         groupvars = c("clade_5", 
                                                                       "Sex"), 
                                                         na.rm = TRUE)) %>% 
  lapply(FUN = \(x){colnames(x)[4] <- "mean"; x})


# Create morph summary by sex and ploidy level
morph_ploidy_summary <- lapply(morph_vars, 
                               FUN = \(x)Rmisc::summarySE(morph_measures, 
                                                          measurevar = x, 
                                                          groupvars = c("clade_2",
                                                                        "Sex"), 
                                                          na.rm = TRUE)) %>% 
  lapply(FUN = \(x){colnames(x)[4] <- "mean"; x})

# Create morph minmax by sex and genetic group
morph_clade_minmax <- lapply(morph_vars, 
                             FUN = \(x)summarise(group_by(morph_measures, 
                                                          clade_5, Sex),
                                                             min = min(get(x),
                                                                       na.rm = TRUE),
                                                             max = max(get(x),
                                                                       na.rm = TRUE)))


# Create morph minmax by sex and ploidy level
morph_ploidy_minmax <- lapply(morph_vars,
                              FUN = \(x)summarise(group_by(morph_measures, 
                                                           clade_2, Sex),
                                                  min = min(get(x),
                                                            na.rm = TRUE),
                                                  max = max(get(x),
                                                            na.rm = TRUE)))


# Set names in all summaries and minmax
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



# Combine all summaries in a single object
all_summary <- bind_rows(morph_ploidy_summary, 
                         morph_clade_summary) %>% 
  select(variable, clade, Sex, N, mean, sd, se, ci) %>% 
  arrange(factor(variable, levels = c(morph_vars)),
                 factor(clade, levels = c("4x", 
                                          "2x",
                                          "Tetraploid",
                                          "Hercynian",
                                          "Algarve",
                                          "Doñana",
                                          "Cádiz")))

# Combine all minmax in a single object
all_minmax <- bind_rows(morph_ploidy_minmax, 
                        morph_clade_minmax)  %>% 
  arrange(factor(variable, levels = c(morph_vars)),
          factor(clade, levels = c("4x", 
                                   "2x",
                                   "Tetraploid",
                                   "Hercynian",
                                   "Algarve",
                                   "Doñana",
                                   "Cádiz")))

# Combine global summary and minmax in a single object
all_summary <- all_minmax %>%
  left_join(all_summary) %>% 
  select(variable, clade, Sex, N, min, max, mean, sd, se, ci)


# ----------------------------------------------------------------------
# EXPORT RESULTS
# ----------------------------------------------------------------------

# Write table to excel 
all_summary %>% 
  mutate(across(min:ci, \(x)round(x, digits = 2))) %>% 
  openxlsx::write.xlsx("outputs/tables/Table_morph_sex.xlsx")

# Remove temporary objects
rm(morph_ploidy_summary, 
   morph_clade_summary,
   all_summary)

rm(morph_ploidy_minmax, 
   morph_clade_minmax,
   all_minmax)



# ----------------------------------------------------------------------
# RUN ANALYSIS WITH NO PARAMETRIC APPROACHES
# ----------------------------------------------------------------------

# Define function to run the comparisons with a Dunn test
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

# Wrangle and prepare data for downstream analysis.
morph_Algarve_clade_diff_letters <- lapply(morph_vars, FUN = extract_differences_from_dunnTest, "Sex", morph_measures %>% filter(clade_5 == "Algarve"))
morph_Cadiz_clade_diff_letters <- lapply(morph_vars, FUN = extract_differences_from_dunnTest, "Sex", morph_measures %>% filter(clade_5 == "Cádiz"))
morph_Doñana_clade_diff_letters <- lapply(morph_vars, FUN = extract_differences_from_dunnTest, "Sex", morph_measures %>% filter(clade_5 == "Doñana"))
morph_Hercynian_clade_diff_letters <- lapply(morph_vars, FUN = extract_differences_from_dunnTest, "Sex", morph_measures %>% filter(clade_5 == "Hercynian"))
morph_Tetraploid_clade_diff_letters <- lapply(morph_vars, FUN = extract_differences_from_dunnTest, "Sex", morph_measures %>% filter(clade_5 == "Tetraploid"))

# Wrangle and prepare data for downstream analysis.
morph_2x_ploidy_diff_letters <- lapply(morph_vars, FUN = extract_differences_from_dunnTest, "Sex", morph_measures %>% filter(clade_2 == "2x"))
morph_4x_ploidy_diff_letters <- lapply(morph_vars, FUN = extract_differences_from_dunnTest, "Sex", morph_measures %>% filter(clade_2 == "4x"))

# Set names for the temporary objects
names(morph_Algarve_clade_diff_letters) <- 
  names(morph_Cadiz_clade_diff_letters) <- 
  names(morph_Doñana_clade_diff_letters) <- 
  names(morph_Hercynian_clade_diff_letters) <-
  names(morph_Tetraploid_clade_diff_letters) <- 
  names(morph_2x_ploidy_diff_letters)  <- 
  names(morph_4x_ploidy_diff_letters) <- morph_vars

# Combine comparisons in a single object
clade_diff_letters <- list("Algarve" = morph_Algarve_clade_diff_letters,
                           "Cádiz" =  morph_Cadiz_clade_diff_letters,
                           "Doñana" = morph_Doñana_clade_diff_letters,
                           "Hercynian" =  morph_Hercynian_clade_diff_letters,
                           "Tetraploid" = morph_Tetraploid_clade_diff_letters) %>% 
  melt(id.vars =c("Comparison", "Z", "P.unadj", "P.adj")) %>% 
  rename(name = L2, 
         clade_5 = L1)

ploidy_diff_letters <- list("2x" = morph_2x_ploidy_diff_letters,
                            "4x" = morph_4x_ploidy_diff_letters) %>% 
  melt(id.vars =c("Comparison", "Z", "P.unadj", "P.adj")) %>% 
  select(Comparison, Z, P.unadj, P.adj, L2, L1) %>% 
  rename(name = L2, 
         clade_2 = L1)

# ----------------------------------------------------------------------
# DRAW FIGURES
# ----------------------------------------------------------------------

# Define colores for todos los gráficos
clade_colors <- c("#009E73", "#F0E442", "#CC79A7", "#E69F00", "#6388b4")

# Pivot data object in longer format
morph_lf <- morph_measures %>%
  select(Population_ID, Individual_ID, Sex, clade_2, clade_5, all_of(morph_vars)) %>% 
  pivot_longer(cols = all_of(morph_vars)) %>% 
  mutate(name = factor(name, levels = morph_vars))

# Prepare letters for plotting
clade_diff_letters <- clade_diff_letters %>% 
  mutate(name = factor(name, 
                       levels = morph_vars)) %>% 
  mutate(label = case_when(
    P.adj > 0.05 ~ "ns",
    P.adj <= 0.05 & P.adj > 0.01 ~ "*",
    P.adj <= 0.01 & P.adj > 0.001 ~ "**",
    P.adj <= 0.001 ~ "***")) %>% 
  mutate(Sex = "F")

ploidy_diff_letters <- ploidy_diff_letters %>% 
  mutate(name = factor(name, 
                       levels = morph_vars)) %>% 
  mutate(label = case_when(
    P.adj > 0.05 ~ "ns",
    P.adj <= 0.05 & P.adj > 0.01 ~ "*",
    P.adj <= 0.01 & P.adj > 0.001 ~ "**",
    P.adj <= 0.001 ~ "***")) %>% 
  mutate(Sex = "F")


#### PLOIDY PLOTS ####

# Set variables for main figure or supplementary information
main_vars <- c("co", "ca", "ltooth", "hair", "fl_infl", "d_infl", "prop_ltooth")

suppl_vars <- c("stoma_w", "stoma_l", "d_max", "d_min")

# Construct the plot object.
measures_ploidy_plot <- ggplot(morph_lf %>% 
         filter(name %in% main_vars), 
       aes(x = clade_2, 
           y = value, 
           fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, 
              trim = FALSE) + 
  geom_boxplot(width = 0.5, 
               alpha = 0.6) +
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "point",
               size = 2,
               position = position_dodge(0.5)) +
  scale_y_continuous(name = "Length",
                     expand = expansion(mult = c(0.05, 
                                                 0.15))) +
  scale_fill_manual(values = c("F" = "gold3",
                               "H" = "seagreen3")) +  
  labs(fill = "Sex") +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  geom_text(data = ploidy_diff_letters %>% 
              filter(name %in% main_vars),
            aes(label = label),
            y = Inf,
            vjust = 1.4,
            size = 4) +
  facet_wrap(vars(name), 
             ncol = 4,
             scales = "free_y") 


measures_ploidy_plot <- measures_ploidy_plot %>%
  lemon::reposition_legend(x = 0.5,
                           y = 0.7,
                           just = c("center", "center"),
                           panel = "panel-4-2")

measures_ploidy_plot



#### CLADE PLOT ####

# Create plot object
measures_clade_plot <- ggplot(morph_lf %>% 
                                filter(name %in% main_vars), 
                              aes(x = clade_5, 
                                  y = value, 
                                  fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, 
                                  trim = FALSE) + 
  geom_boxplot(width = 0.5, 
               alpha = 0.6) +
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "point",
               size = 2,
               position = position_dodge(0.5)) +
  scale_y_continuous(name = "Length",
                     expand = expansion(mult = c(0.05, 
                                                 0.15))) +
  scale_fill_manual(values = c("F" = "gold3",
                               "H" = "seagreen3")) +  
  labs(fill = "Sex") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 1, 
                                   hjust=1),
        axis.title.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  geom_text(data = clade_diff_letters %>% 
              filter(name %in% main_vars),
            aes(label = label),
            y = Inf,
            vjust = 1.2,
            size = 4) +
  facet_wrap(vars(name), 
             ncol = 2,
             scales = "free_y") 


measures_clade_plot <- measures_clade_plot %>%
  lemon::reposition_legend(x = 0.5,
                           y = 0.4,
                           just = c("center", "center"),
                           panel = "panel-2-4")

measures_clade_plot


# ----------------------------------------------------------------------
# Sex effect on corolla comparisons
# ----------------------------------------------------------------------

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
  scale_y_continuous(name = "Corolla",
                     expand = expansion(mult = c(0.05, 
                                                 0.15))) +
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

plot_sex <- ggarrange(measures_ploidy_plot,
                      ggarrange(plot_sex_2,
                                plot_sex_3,
                                widths = c(1, 1.5),
                                ncol = 2,
                                labels = c("B", "C")),
                      nrow = 2,
                      heights = c(2, 1),
                      labels = c("A",""))

plot_sex

# Export figure files in publication-ready formats.
ggsave("outputs/figures/Figure_S2.pdf",
       plot_sex,
       width = 8,
       height = 7.5)
ggsave("outputs/figures/Figure_S2.png",
       plot_sex,
       width = 8,
       height = 7.5)


