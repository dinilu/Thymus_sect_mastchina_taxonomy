# ======================================================================
# Script: 03_Abiotic.R
# Purpose: Analyse abiotic (soil) variables across genetic groups and ploidy levels; generate Figure 4 and ordination outputs.
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
# Load libraries
# ----------------------------------------------------------------------
# Package requirements
# ----------------------------------------------------------------------
# Load all R packages required for data import, manipulation, modelling,
# and figure/table generation.
library(tidyverse)
library(ggpubr)
library(factoextra)
library(FSA)
library(reshape2)
library(introdataviz)


# ----------------------------------------------------------------------
# Data import and preprocessing
# ----------------------------------------------------------------------
# Loading data ------------------------------------------------------------

# Load genetic group data

pops_data <- "data/genetic groups.xlsx" %>% 
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


# Load soil data
soil_vars <- c("pH", "EC", "OM", "WRP", "WRC")

soil_data <- "data/Soil_samples.xlsx" %>%  
  openxlsx::read.xlsx(2) %>%  
  filter(Sample_location == "b") %>% 
  mutate(Population_ID = as.factor(Population_ID)) %>% 
  left_join(pops_data) %>% 
  mutate(across(pH:WRC, log)) 


# Check for the total number of data in each population
test <- soil_data %>%  
  group_by(Population_ID) %>%  
  summarize(count = n())

test

# Soil analysis -----------------------------------------------------------

# Compare significance between groups 
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

soil_clade_diff_letters <- lapply(soil_vars, 
                                  FUN = extract_differences_from_dunnTest, 
                                  "clade_5", 
                                  soil_data)
soil_ploidy_diff_letters <- lapply(soil_vars, 
                                  FUN = extract_differences_from_dunnTest, 
                                  "clade_2", 
                                  soil_data)

names(soil_clade_diff_letters) <- 
  names(soil_ploidy_diff_letters) <- 
  soil_vars


soil_clade_diff_letters <- soil_clade_diff_letters %>% 
  melt(id.vars =c("Group", "Letter")) %>% 
  rename(name = L1,
         clade_5 = Group,
         label = Letter) %>% 
  mutate(name = factor(name, levels = soil_vars)) %>% 
  mutate(taxa_treatment = "Group") %>%
  mutate(taxa_treatment = factor(taxa_treatment, 
                                 levels = c("Ploidy",
                                            "Group"))) %>% 
  rename(taxa_name = clade_5) %>% 
  mutate(taxa_name = factor(taxa_name, 
                            levels = c("2x", 
                                       "4x",
                                       "Algarve",
                                       "Cádiz",
                                       "Doñana",
                                       "Hercynian",
                                       "Tetraploid")))

soil_ploidy_diff_letters <- soil_ploidy_diff_letters %>% 
  melt(id.vars =c("Comparison", "Z", "P.unadj", "P.adj")) %>% 
  select(Comparison, Z, P.unadj, P.adj, L1) %>% 
  rename(name = L1) %>% 
  mutate(clade_2 = "2x")%>% 
  mutate(name = factor(name, 
                       levels = soil_vars)) %>% 
  mutate(label = case_when(
    P.adj > 0.05 ~ "ns",
    P.adj <= 0.05 & P.adj > 0.01 ~ "*",
    P.adj <= 0.01 & P.adj > 0.001 ~ "**",
    P.adj <= 0.001 ~ "***")) %>% 
  mutate(taxa_treatment = "Ploidy") %>%
  mutate(taxa_treatment = factor(taxa_treatment, 
                                 levels = c("Ploidy",
                                            "Group"))) %>% 
  rename(taxa_name = clade_2) %>%
  mutate(taxa_name = factor(taxa_name, 
                            levels = c("2x",
                                       "4x",
                                       "Algarve",
                                       "Cádiz",
                                       "Doñana",
                                       "Hercynian",
                                       "Tetraploid"))) %>% 
  select(taxa_name, label, name, taxa_treatment)



# ----------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------
#### Figures ####

# Define colors for all plots

clade_colors <- c(Algarve = "#009E73", # Algarve
                   Cádiz = "#F0E442", # Cadiz
                   Doñana = "#CC79A7", # Doñana
                   Hercynian = "#E69F00", # Hercynian
                   Tetraploid = "#6388b4") # Tetraploid

ploidy_colors <- c("2x" = "#E69F00", # 2x 
                   "4x" = "#6388b4") # 4x


# Prepare final data for plotting
soil_data_lf <- soil_data %>%
  select(Population_ID, Sample_ID, clade_2, clade_5, all_of(soil_vars)) %>% 
  pivot_longer(cols = all_of(soil_vars)) %>% 
  mutate(name = factor(name, levels = soil_vars))

soil_ploidy_lf <- soil_data_lf %>% 
  mutate(taxa_treatment = "Ploidy") %>% 
  rename(taxa_name = clade_2)

soil_clade_lf <- soil_data_lf %>% 
  mutate(taxa_treatment = "Group") %>% 
  rename(taxa_name = clade_5)

soil_longer_data <- soil_ploidy_lf %>% 
  bind_rows(soil_clade_lf) %>% 
  mutate(taxa_treatment = factor(taxa_treatment, 
                                 levels = c("Ploidy", 
                                            "Group"))) %>% 
  mutate(taxa_name = factor(taxa_name, 
                            levels = c("2x",
                                       "4x",
                                       "Algarve",
                                       "Cádiz",
                                       "Doñana",
                                       "Hercynian",
                                       "Tetraploid")))


# Create combined fiture 
fig_4 <- ggplot(soil_longer_data) +
  introdataviz::geom_split_violin(data = soil_longer_data %>% 
                                    filter(taxa_treatment == "Ploidy"),
                                  aes(x = "2x vs 4x", 
                                      y = value, 
                                      fill = taxa_name),
                                  alpha = 0.4, 
                                  trim = FALSE) +
  geom_boxplot(data = soil_longer_data %>%
                 filter(taxa_treatment == "Ploidy"),
               aes(x = "2x vs 4x",
                   y = value,
                   fill = taxa_name),
               width = 0.15, 
               alpha = 0.6,
               position = position_dodge(0.15)) +
  stat_summary(
    data = soil_longer_data %>%
      filter(taxa_treatment == "Ploidy"),
    aes(x = "2x vs 4x",
        y = value,
        group = taxa_name),
    fun = mean,
    geom = "point",
    color = "black",
    position = position_dodge(0.15),
    size = 1.2) +
  geom_violin(data = soil_longer_data %>% 
                filter(taxa_treatment == "Group"), 
              aes(x = taxa_name, 
                  y = value,
                  fill = taxa_name), 
              alpha = 0.4, 
              trim = FALSE) +
  geom_boxplot(data = soil_longer_data %>%
                 filter(taxa_treatment == "Group"),
               aes(x = taxa_name,
                   y = value,
                   fill = taxa_name),
               width = 0.08, 
               alpha = 0.6) +
  stat_summary(
    data = soil_longer_data %>%
      filter(taxa_treatment == "Group"),
    aes(x = taxa_name, y = value),
    fun = mean,
    geom = "point",
    color = "black",
    position = position_dodge(0.9),
    size = 1.2) +
  geom_text(data = soil_clade_diff_letters,
            aes(x = taxa_name,
                label = label),
            y = Inf,
            vjust = 1) +
  geom_text(data = soil_ploidy_diff_letters,
            aes(x = 1,
                label = label),
            y = Inf,
            vjust = 1.2,
            size = 5) +
  facet_grid(name ~ taxa_treatment,
             scales = "free",
             labeller = labeller(
               taxa_treatment = c("Ploidy" = "Ploidy levels", 
                                  "Group" = "Genetic groups"))) +
  theme_bw() +
  labs(x = "", 
       y = "Values (log transf.)",
       fill = "") +
  scale_fill_manual(values = c(ploidy_colors,
                               clade_colors)) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 
                                                 0.25))) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1))

fig_4


# Export figure files in publication-ready formats.
ggsave("outputs/figures/Figure_4.pdf",
       fig_4,
       width = 7,
       height = 7.5)
ggsave("outputs/figures/Figure_4.png",
       fig_4,
       width = 7,
       height = 7.5)

# Plot multivariate ordination based on soil variables
soil_pca_data <- soil_data %>% 
  select(-Comments) %>% 
  na.omit() %>% 
  mutate(shape_group = ifelse(clade_5 == "Tetraploid", 22, 21))

soil_pca <- soil_pca_data %>% 
  select(all_of(soil_vars)) %>% 
  prcomp(scale. = TRUE) 

soil_pca_summary <- soil_pca %>% 
  summary() %>% 
  .$importance

# Screeplot 
fig_S5_a <- fviz_screeplot(soil_pca,
               addlabels = TRUE,
               title = "") +
  labs(y = "Explained variance (%)") +
  ylim(0, 65) +
  theme_bw()

fig_S5_a

# PCA contribution among variables
fig_S5_b <- fviz_pca_var(soil_pca, 
             col.var="contrib",
             gradient.cols = c("#0072B2",
                               "#F0E442",
                               "#E69F00"),
             repel = TRUE,
             title = "") +
  theme_bw() + 
  labs(x = paste0("Principal Component 1 (",
                  sprintf("%.1f",
                          100 * soil_pca_summary[2, 1]), 
                  "%)"),
       y = paste0("Principal Component 2 (", 
                  sprintf("%.1f",
                          100 * soil_pca_summary[2, 2]),
                  "%)"),
       color = "Contribution")

fig_S5_b


# PCA among ploidy levels
fig_S5_c <- fviz_pca_ind(soil_pca, 
                label = "var", 
                habillage = soil_pca_data$clade_2, 
                addEllipses = TRUE,
                ellipse.type = "convex",
                ellipse.alpha = 0.05, 
                title = "",
                shape.ind = soil_pca_data$clade_2,
                col.ind = soil_pca_data$clade_2,  
                fill.ind = soil_pca_data$clade_2, 
                mean.point.size = 5,
                alpha.ind = 0.3) +
  theme_bw() + 
  guides(color = guide_legend(title="Taxa"),
         shape = guide_legend(title="Taxa"),
         fill = guide_legend(title="Taxa")) +
  scale_color_manual(values = ploidy_colors) +
  scale_fill_manual(values = ploidy_colors) +
  scale_shape_manual(values = c("4x" = 22, 
                                "2x" = 21)) +   
  theme(legend.position = "none") +
  labs(x = paste0("Principal Component 1 (",
                  sprintf("%.1f",
                          100 * soil_pca_summary[2, 1]), 
                  "%)"),
       y = paste0("Principal Component 2 (", 
                  sprintf("%.1f",
                          100 * soil_pca_summary[2, 2]),
                  "%)"),
       color = "Contribution")

fig_S5_c


# PCA among genetic groups
fig_S5_d <- fviz_pca_ind(soil_pca, 
                         label = "var", 
                         geom = "point",
                         habillage = soil_pca_data$clade_5,  
                         addEllipses = TRUE,
                         ellipse.type = "convex", 
                         ellipse.alpha = 0.05,
                         title = "",
                         shape.ind = soil_pca_data$clade_5,
                         col.ind = soil_pca_data$clade_5,  
                         fill.ind = soil_pca_data$clade_5, 
                         mean.point.size = 5,
                         alpha.ind = 0.3) +
  scale_shape_manual(values = c("Tetraploid" = 22,
                                "Cádiz" = 21,    
                                "Doñana" = 21,
                                "Algarve" = 21,
                                "Hercynian" = 21
  ))  +
  scale_color_manual(values = clade_colors)   +
  scale_fill_manual(values = clade_colors)   +
  theme_bw() + 
  guides(color = guide_legend(title="Taxa"),
         fill = guide_legend(title = "Taxa")
  ) +
  theme(legend.position = "none") +
  labs(x = paste0("Principal Component 1 (",
                  sprintf("%.1f",
                          100 * soil_pca_summary[2, 1]),
                  "%)"),
       y = paste0("Principal Component 2 (",
                  sprintf("%.1f", 
                          100 * soil_pca_summary[2, 2]), 
                  "%)"))

fig_S5_d


# Combine all plots in a single figure
fig_S5 <- ggpubr::ggarrange(fig_S5_a,
                           fig_S5_b,
                           fig_S5_c,
                           fig_S5_d,
                           ncol = 2,
                           nrow = 2,
                           labels = c("A", "B", "C", "D"))

fig_S5

# Export figure files in publication-ready formats.
ggsave("outputs/figures/Figure_S5.pdf", 
       fig_S5,
       width = 7, 
       height = 7)
ggsave("outputs/figures/Figure_S5.png", 
       fig_S5,
       width = 7, 
       height = 7)



# Bioclim analysis --------------------------------------------------------

# Load bioclimatic Chelsa data

# Load presaved version of the bioclimatic data
bioclim_vars <- c(c("ai"),
                  paste0("bio", 1:19), 
                  c("cmi_mean"),
                  c("gdd5"),
                  c("sfcWind_mean"))

bioclim_data <- "data/Abiotic_population_means.xlsx" %>% 
  openxlsx::read.xlsx("climate") %>% 
  select(c(Population_ID, all_of(bioclim_vars))) %>% 
  mutate(Population_ID = as.factor(Population_ID)) %>% 
  left_join(pops_data)


# Compare significance between groups 

# Plot differences for bioclim variables
bioclim_pca_data <- bioclim_data %>%
  select(all_of(bioclim_vars))

bioclim_pca <- bioclim_pca_data %>% 
  prcomp(scale. = TRUE) 

bioclim_pca_vars <- paste0("PC", 1:5)

bioclim_data <- bioclim_data %>% 
  cbind(., bioclim_pca$x %>% 
          as.data.frame() %>% 
          select(paste0("PC", 1:5)))

bioclim_longer_data <- bioclim_data %>% 
  pivot_longer(all_of(bioclim_pca_vars),
               names_to = "variable",
               values_to = "value")


bioclim_clade_diff_letters <- lapply(bioclim_pca_vars, 
                                  FUN = extract_differences_from_dunnTest, 
                                  "clade_5", 
                                  bioclim_data)
bioclim_ploidy_diff_letters <- lapply(bioclim_pca_vars, 
                                   FUN = extract_differences_from_dunnTest, 
                                   "clade_2", 
                                   bioclim_data)

names(bioclim_clade_diff_letters) <- names(bioclim_ploidy_diff_letters) <- bioclim_pca_vars


bioclim_clade_diff_letters <- bioclim_clade_diff_letters %>% 
  melt(id.vars =c("Group", "Letter")) %>% 
  rename(name = L1,
         clade_5 = Group,
         label = Letter) %>% 
  mutate(name = factor(name, levels = bioclim_pca_vars)) %>% 
  mutate(taxa_treatment = "Group") %>%
  mutate(taxa_treatment = factor(taxa_treatment, 
                                 levels = c("Ploidy",
                                            "Group"))) %>% 
  rename(taxa_name = clade_5) %>% 
  mutate(taxa_name = factor(taxa_name, 
                            levels = c("2x", 
                                       "4x",
                                       "Algarve",
                                       "Cádiz",
                                       "Doñana",
                                       "Hercynian",
                                       "Tetraploid")))


bioclim_ploidy_diff_letters <- bioclim_ploidy_diff_letters %>% 
  melt(id.vars =c("Comparison", "Z", "P.unadj", "P.adj")) %>% 
  select(Comparison, Z, P.unadj, P.adj, L1) %>% 
  rename(name = L1) %>% 
  mutate(clade_2 = "2x") %>% 
  mutate(name = factor(name, 
                       levels = bioclim_pca_vars)) %>% 
  mutate(label = case_when(
    P.adj > 0.05 ~ "ns",
    P.adj <= 0.05 & P.adj > 0.01 ~ "*",
    P.adj <= 0.01 & P.adj > 0.001 ~ "**",
    P.adj <= 0.001 ~ "***")) %>% 
  mutate(taxa_treatment = "Ploidy") %>%
  mutate(taxa_treatment = factor(taxa_treatment, 
                                 levels = c("Ploidy",
                                            "Group"))) %>% 
  rename(taxa_name = clade_2) %>%
  mutate(taxa_name = factor(taxa_name, 
                            levels = c("2x",
                                       "4x",
                                       "Algarve",
                                       "Cádiz",
                                       "Doñana",
                                       "Hercynian",
                                       "Tetraploid"))) %>% 
  select(taxa_name, label, name, taxa_treatment)



# Prepare final data for plotting
bioclim_data_lf <- bioclim_data %>%
  select(Population_ID, clade_2, clade_5, all_of(bioclim_pca_vars)) %>% 
  pivot_longer(cols = all_of(bioclim_pca_vars)) %>% 
  mutate(name = factor(name, levels = bioclim_pca_vars))

bioclim_ploidy_lf <- bioclim_data_lf %>% 
  mutate(taxa_treatment = "Ploidy") %>% 
  rename(taxa_name = clade_2)

bioclim_clade_lf <- bioclim_data_lf %>% 
  mutate(taxa_treatment = "Group") %>% 
  rename(taxa_name = clade_5)

bioclim_longer_data <- bioclim_ploidy_lf %>% 
  bind_rows(bioclim_clade_lf) %>% 
  mutate(taxa_treatment = factor(taxa_treatment, 
                                 levels = c("Ploidy", 
                                            "Group"))) %>% 
  mutate(taxa_name = factor(taxa_name, 
                            levels = c("2x",
                                       "4x",
                                       "Algarve",
                                       "Cádiz",
                                       "Doñana",
                                       "Hercynian",
                                       "Tetraploid")))


# Create plot
fig_5 <- ggplot(bioclim_longer_data) +
  introdataviz::geom_split_violin(data = bioclim_longer_data %>% 
                                    filter(taxa_treatment == "Ploidy"),
                                  aes(x = "2x vs 4x", 
                                      y = value, 
                                      fill = taxa_name),
                                  alpha = 0.4, 
                                  trim = FALSE) +
  geom_boxplot(data = bioclim_longer_data %>%
                 filter(taxa_treatment == "Ploidy"),
               aes(x = "2x vs 4x",
                   y = value,
                   fill = taxa_name),
               width = 0.15, 
               alpha = 0.6,
               position = position_dodge(0.15)) +
  stat_summary(
    data = bioclim_longer_data %>%
      filter(taxa_treatment == "Ploidy"),
    aes(x = "2x vs 4x",
        y = value,
        group = taxa_name),
    fun = mean,
    geom = "point",
    color = "black",
    position = position_dodge(0.15),
    size = 1.2) +
  geom_violin(data = bioclim_longer_data %>% 
                filter(taxa_treatment == "Group"), 
              aes(x = taxa_name, 
                  y = value,
                  fill = taxa_name), 
              alpha = 0.4, 
              trim = FALSE) +
  geom_boxplot(data = bioclim_longer_data %>%
                 filter(taxa_treatment == "Group"),
               aes(x = taxa_name,
                   y = value,
                   fill = taxa_name),
               width = 0.08, 
               alpha = 0.6) +
  stat_summary(
    data = bioclim_longer_data %>%
      filter(taxa_treatment == "Group"),
    aes(x = taxa_name, y = value),
    fun = mean,
    geom = "point",
    color = "black",
    position = position_dodge(0.9),
    size = 1.2) +
  geom_text(data = bioclim_clade_diff_letters,
            aes(x = taxa_name,
                label = label),
            y = Inf,
            vjust = 1) +
  geom_text(data = bioclim_ploidy_diff_letters,
            aes(x = 1,
                label = label),
            y = Inf,
            vjust = 1.1,
            size = 4) +
  facet_grid(name ~ taxa_treatment,
             scales = "free",
             labeller = labeller(
               taxa_treatment = c("Ploidy" = "Ploidy levels", 
                                  "Group" = "Genetic groups"))) +
  theme_bw() +
  labs(x = "", 
       y = "Values (log transf.)",
       fill = "") +
  scale_fill_manual(values = c(ploidy_colors,
                               clade_colors)) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 
                                                 0.25))) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1))

fig_5

# Export figure files in publication-ready formats.
ggsave("outputs/figures/Figure_5.pdf",
       fig_5,
       width = 7,
       height = 7.5)
ggsave("outputs/figures/Figure_5.png",
       fig_5,
       width = 7,
       height = 7.5)


# PCA analysis

# Screeplot
fig_s6_a <- fviz_screeplot(bioclim_pca,
                          addlabels = TRUE,
                          title = "") +
  theme_bw() +
  labs(y = "Explained variance (%)") +
  ylim(0, 55)

fig_s6_a

bioclim_pca_var <- get_pca_var(bioclim_pca)

corrplot::corrplot(bioclim_pca_var$cos2,
         is.corr=FALSE)


# PCA contribution
bioclim_pca_summary <- bioclim_pca %>% 
  summary() %>% 
  .$importance

fig_s6_b <- fviz_pca_var(bioclim_pca, 
                        col.var="contrib",
                        gradient.cols = c("#0072B2",
                                          "#F0E442",
                                          "#E69F00"),
                        title = "",
                        repel = TRUE) +
  theme_bw() +
  labs(x = paste0("Principal Component 1 (",
                  sprintf("%.1f",
                          100 * bioclim_pca_summary[2, 1]), "%)"),
       y = paste0("Principal Component 2 (", 
                  sprintf("%.1f",
                          100 * bioclim_pca_summary[2, 2]), "%)")) +
  labs(color = "Contribution")  


fig_s6_b

# PCA among ploidy levels
fig_s6_c <- fviz_pca_ind(bioclim_pca, 
                         label = "var", 
                         habillage = bioclim_data$clade_2, 
                         addEllipses = TRUE,
                         ellipse.type = "convex",
                         ellipse.alpha = 0.05, 
                         title = "",
                         shape.ind = bioclim_data$clade_2,
                         col.ind = bioclim_data$clade_2,  
                         fill.ind = bioclim_data$clade_2, 
                         mean.point.size = 5,
                         alpha.ind = 0.3) +
  theme_bw() + 
  guides(color = guide_legend(title="Taxa"),
         shape = guide_legend(title="Taxa"),
         fill = guide_legend(title="Taxa")) +
  scale_color_manual(values = ploidy_colors) +
  scale_fill_manual(values = ploidy_colors) +
  scale_shape_manual(values = c("4x" = 22,
                                "2x" = 21)) +   
  theme(legend.position = "none") +
  labs(x = paste0("Principal Component 1 (",
                  sprintf("%.1f",
                          100 * bioclim_pca_summary[2, 1]), 
                  "%)"),
       y = paste0("Principal Component 2 (", 
                  sprintf("%.1f",
                          100 * bioclim_pca_summary[2, 2]),
                  "%)"),
       color = "Contribution")

fig_s6_c


# PCA among genetic groups
fig_s6_d <- fviz_pca_ind(bioclim_pca, 
                         label = "var",  
                         geom = "point",
                         habillage = bioclim_data$clade_5,   
                         addEllipses = TRUE,
                         ellipse.type = "convex",  
                         ellipse.alpha = 0.05,
                         title = "",
                         shape.ind = bioclim_data$clade_5,
                         col.ind = bioclim_data$clade_5,  
                         fill.ind = bioclim_data$clade_5, 
                         mean.point.size = 5,
                         alpha.ind = 0.3) +
  scale_shape_manual(values = c("Tetraploid" = 22,
                                "Cádiz" = 21,    
                                "Doñana" = 21,
                                "Algarve" = 21,
                                "Hercynian" = 21 ))  +
  scale_color_manual(values = clade_colors)   +
  scale_fill_manual(values = clade_colors)   +
  theme_bw() + 
  guides(color = guide_legend(title="Taxa"),
         fill = guide_legend(title = "Taxa")) +
  theme(legend.position = "none") +
  labs(x = paste0("Principal Component 1 (",
                  sprintf("%.1f",
                          100 * soil_pca_summary[2, 1]),
                  "%)"),
       y = paste0("Principal Component 2 (",
                  sprintf("%.1f", 
                          100 * soil_pca_summary[2, 2]), 
                  "%)"))

fig_s6_d


# Combine all plots in a single figure
fig_s6 <- ggpubr::ggarrange(fig_s6_a,
                            fig_s6_b,
                            fig_s6_c,
                            fig_s6_d,
                            ncol = 2,
                            nrow = 2,
                            labels = c("A", "B", "C", "D"))

fig_s6

# Export figure files in publication-ready formats.
ggsave("outputs/figures/Figure_S6.pdf", 
       fig_s6,
       width = 7, 
       height = 7)
ggsave("outputs/figures/Figure_S6.png", 
       fig_s6,
       width = 7, 
       height = 7)

