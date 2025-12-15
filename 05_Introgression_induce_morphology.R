# Load libraries
library(tidyverse)
library(openxlsx)
library(glmmTMB)
library(sjPlot)
# library(ggeffects)
library(bbmle)
# library(DHARMa)


# Loading data ------------------------------------------------------------

# Load genetic groups

genetic_groups <- openxlsx::read.xlsx("data/genetic groups.xlsx", 
                                      2) %>% 
  mutate(genetic_group_k2 = case_match(genetic_group_k2,
                                       "tetraploid" ~ "4x",
                                       "diploid" ~ "2x")) %>% 
  mutate(genetic_group_k5 = str_to_sentence(genetic_group_k5)) %>% 
  rename(clade_2 = genetic_group_k2,
         clade_5 = genetic_group_k5)


# Load morphological data

morpho_data <- "data/Morphometric_measures.xlsx" %>% 
  read.xlsx(2) %>% 
  group_by(Population_ID, Individual_ID) %>% 
  summarise(across(Sex, first), across(co:prop_ltooth, mean)) %>% 
  ungroup()


# Load genetic ancestry data

ancestry <- "data/Ancestry.xlsx" %>% 
  read.xlsx(1) %>% 
  mutate(across(group1:group5, ~.x*100 %>%
                  as_tibble() %>%
                  pull(1)))


# Combine bioclim, soil, introgresion and morphology data

data_df <- morpho_data %>%
  inner_join(ancestry) %>% 
  inner_join(genetic_groups) %>%
  dplyr::select(-stooth) %>% 
  # na.omit() %>% 
  mutate(Sex = as.factor(Sex))


# Define groups of variables

morpho_vars <- c("co", "ca", "ltooth", "hair", "d_infl", "prop_ltooth")

intro_vars <- paste0("group", c(3, 5, 4, 2, 1))

clades <- c("Algarve", "Cádiz", "Doñana", "Hercynian", "Tetraploid")


# Exploratory plots

test <- data_df %>%
  pivot_longer(cols = all_of(morpho_vars), 
               names_to = "morpho_var",
               values_to = "morpho_value") %>% 
  pivot_longer(cols = all_of(c(intro_vars)), 
               names_to = "ind_var",
               values_to = "ind_value")

ggplot(test,
       aes(y = morpho_value, 
           x = ind_value, 
           colour = clade_5)) +
  geom_point(alpha = 0.5) +
  facet_grid(cols = vars(ind_var),
             rows = vars(morpho_var),
             scales = "free")
  



fit_gg_models <- function(n, clades, indep_vars, dat, dep_var){
  # n <- 1
  # indep_vars <- intro_vars
  # dat <- data_df
  # dep_var <- morpho_vars[1]
  
  cl <- clades[n]
  
  indep_vars <- indep_vars[-n]
  
  d_cl <- dat %>% 
    filter(clade_5 == cl)
    
  form_0 <- paste0(dep_var, " ~ ", 
                   paste(indep_vars,
                         collapse = " + ")) %>%
    as.formula()
  form_1 <- paste0(dep_var, " ~ ", 
                   paste(indep_vars, 
                         collapse = " + "), 
                   " + (1 | Sex)") %>% 
    as.formula()
  
  model_0 <- d_cl %>% 
    glmmTMB(form_0,
            data = ., 
            family = gaussian(),
            na.action = "na.exclude")
  
  model_1 <- d_cl %>% 
    glmmTMB(form_1,
            data = .,
            family = gaussian(),
            na.action = "na.exclude")
  
  mod_comp <- AICtab(model_0, model_1)
  
  model_lm <- d_cl %>% 
    lm(form_0,
            data = .)
  
  list(model_0, model_1, mod_comp, model_lm)
}

ca_gg_mods <- lapply(1:5, FUN = fit_gg_models, clades, intro_vars, data_df, "ca")
co_gg_mods <- lapply(1:5, FUN = fit_gg_models, clades, intro_vars, data_df, "co")
dinfl_gg_mods <- lapply(1:5, FUN = fit_gg_models, clades, intro_vars, data_df, "d_infl")
finfl_gg_mods <- lapply(1:5, FUN = fit_gg_models, clades, intro_vars, data_df, "fl_infl")
hair_gg_mods <- lapply(1:5, FUN = fit_gg_models, clades, intro_vars, data_df, "hair")
ltooth_gg_mods <- lapply(1:5, FUN = fit_gg_models, clades, intro_vars, data_df, "ltooth")
propltooth_gg_mods <- lapply(1:5, FUN = fit_gg_models, clades, intro_vars, data_df, "prop_ltooth")

names(ca_gg_mods) <- 
  names(co_gg_mods) <- 
  names(dinfl_gg_mods) <- 
  names(finfl_gg_mods) <-
  names(hair_gg_mods) <-
  names(ltooth_gg_mods) <- 
  names(propltooth_gg_mods) <- 
  clades

lapply(ca_gg_mods, FUN = \(x)x[[3]])
lapply(co_gg_mods, FUN = \(x)x[[3]])
lapply(dinfl_gg_mods, FUN = \(x)x[[3]])
lapply(finfl_gg_mods, FUN = \(x)x[[3]])
lapply(hair_gg_mods, FUN = \(x)x[[3]])
lapply(ltooth_gg_mods, FUN = \(x)x[[3]])
lapply(propltooth_gg_mods, FUN = \(x)x[[3]])


tab_model(lapply(ca_gg_mods, FUN = \(x)x[[4]]))
tab_model(lapply(co_gg_mods, FUN = \(x)x[[4]]))
tab_model(lapply(co_gg_mods, FUN = \(x)x[[2]]))
tab_model(lapply(dinfl_gg_mods, FUN = \(x)x[[4]]))
tab_model(lapply(finfl_gg_mods, FUN = \(x)x[[4]]))
tab_model(lapply(hair_gg_mods, FUN = \(x)x[[4]]))
tab_model(lapply(ltooth_gg_mods, FUN = \(x)x[[4]]))
tab_model(lapply(propltooth_gg_mods, FUN = \(x)x[[4]]))


tab_model(lapply(ca_gg_mods, FUN = \(x)x[[4]]), file = "outputs/ancestry_models/ca_mods.html")
tab_model(lapply(co_gg_mods, FUN = \(x)x[[4]]), file = "outputs/ancestry_models/co_mods_lm.html")
tab_model(lapply(co_gg_mods, FUN = \(x)x[[2]]), file = "outputs/ancestry_models/co_mods_glmm.html")
tab_model(lapply(dinfl_gg_mods, FUN = \(x)x[[4]]), file = "outputs/ancestry_models/dinfl_mods.html")
tab_model(lapply(finfl_gg_mods, FUN = \(x)x[[4]]), file = "outputs/ancestry_models/finfl_mods.html")
tab_model(lapply(hair_gg_mods, FUN = \(x)x[[4]]), file = "outputs/ancestry_models/hair_mods.html")
tab_model(lapply(ltooth_gg_mods, FUN = \(x)x[[4]]), file = "outputs/ancestry_models/ltooth_mods.html")
tab_model(lapply(propltooth_gg_mods, FUN = \(x)x[[4]]), file = "outputs/ancestry_models/propltooth_mods.html")


qq_plot_model <- function(model, name, out_dir = "outputs/diagnostics/QQplots",
                          is_glmm = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  res <- if (is_glmm) residuals(model, type = "response") else residuals(model)
  
  file <- file.path(out_dir, paste0("QQ_", name, ".png"))
  png(file, width = 2000, height = 2000, res = 300)
  qqnorm(res, main = paste("Q–Q plot:", name))
  qqline(res, col = "red")
  dev.off()
}

gg_mods <- list(
  ca          = ca_gg_mods,
  co          = co_gg_mods,
  d_infl      = dinfl_gg_mods,
  fl_infl     = finfl_gg_mods,
  hair        = hair_gg_mods,
  ltooth      = ltooth_gg_mods,
  prop_ltooth = propltooth_gg_mods
)

invisible(
  lapply(names(gg_mods), function(trait) {
    lapply(names(gg_mods[[trait]]), function(clade) {
      mods <- gg_mods[[trait]][[clade]]
      
       ## glmmTMB with Sex random effect
      qq_plot_model(mods[[2]],
                    name   = paste(trait, clade, "glmm_SexRE", sep = "_"),
                    is_glmm = TRUE)
      
      ## linear model
      qq_plot_model(mods[[4]],
                    name   = paste(trait, clade, "lm_noSex", sep = "_"),
                    is_glmm = FALSE)
    })
  })
)




#### UNIVARIATE MODELS ####

fit_gg_single_models <- function(n, indep_vars, dat, dep_var){
  # n <- 1
  # indep_vars <- intro_vars
  # dat <- data_df
  # dep_var <- morpho_vars
  
  clades <- c("Algarve", "Cádiz", "Doñana", "Hercynian", "Tetraploid")
  
  cl <- clades[n]
  # indep_vars <- indep_vars[-n]
  
  d_cl <- dat %>% 
    filter(clade_5 == cl)
  
  form_0_g1 <- paste0(dep_var,
                      " ~ ", 
                      indep_vars[1]) %>%
    as.formula()
  form_0_g2 <- paste0(dep_var,
                    " ~ ", 
                    indep_vars[2]) %>%
    as.formula()
  form_0_g3 <- paste0(dep_var,
                      " ~ ", 
                      indep_vars[3]) %>%
    as.formula()
  form_0_g4 <- paste0(dep_var,
                      " ~ ", 
                      indep_vars[4]) %>%
    as.formula()
  form_0_g5 <- paste0(dep_var,
                      " ~ ", 
                      indep_vars[5]) %>%
    as.formula()
  form_1_g1 <- paste0(dep_var,
                      " ~ ", 
                      indep_vars[1],
                      " + (1 | Sex)") %>%
    as.formula()
  form_1_g2 <- paste0(dep_var,
                      " ~ ", 
                      indep_vars[2],
                      " + (1 | Sex)") %>%
    as.formula()
  form_1_g3 <- paste0(dep_var,
                      " ~ ", 
                      indep_vars[3],
                      " + (1 | Sex)") %>%
    as.formula()
  form_1_g4 <- paste0(dep_var,
                      " ~ ", 
                      indep_vars[4],
                      " + (1 | Sex)") %>%
    as.formula()
  form_1_g5 <- paste0(dep_var,
                      " ~ ", 
                      indep_vars[5],
                      " + (1 | Sex)") %>%
    as.formula()

  model_0_g1 <- d_cl %>% 
    glmmTMB(form_0_g1,
            data = ., 
            family = gaussian(),
            na.action = "na.exclude")
  model_0_g2 <- d_cl %>% 
    glmmTMB(form_0_g2,
            data = ., 
            family = gaussian(),
            na.action = "na.exclude")
  model_0_g3 <- d_cl %>% 
    glmmTMB(form_0_g3,
            data = ., 
            family = gaussian(),
            na.action = "na.exclude")
  model_0_g4 <- d_cl %>% 
    glmmTMB(form_0_g4,
            data = ., 
            family = gaussian(),
            na.action = "na.exclude")
  model_0_g5 <- d_cl %>% 
    glmmTMB(form_0_g5,
            data = ., 
            family = gaussian(),
            na.action = "na.exclude")
  
  model_1_g1 <- d_cl %>% 
    glmmTMB(form_1_g1,
            data = ., 
            family = gaussian(),
            na.action = "na.exclude")
  model_1_g2 <- d_cl %>% 
    glmmTMB(form_1_g2,
            data = ., 
            family = gaussian(),
            na.action = "na.exclude")
  model_1_g3 <- d_cl %>% 
    glmmTMB(form_1_g3,
            data = ., 
            family = gaussian(),
            na.action = "na.exclude")
  model_1_g4 <- d_cl %>% 
    glmmTMB(form_1_g4,
            data = ., 
            family = gaussian(),
            na.action = "na.exclude")
  model_1_g5 <- d_cl %>% 
    glmmTMB(form_1_g5,
            data = ., 
            family = gaussian(),
            na.action = "na.exclude")
  
  mod_comp_g1 <- AICtab(model_0_g1, model_1_g1)
  mod_comp_g2 <- AICtab(model_0_g2, model_1_g2)
  mod_comp_g3 <- AICtab(model_0_g3, model_1_g3)
  mod_comp_g4 <- AICtab(model_0_g4, model_1_g4)
  mod_comp_g5 <- AICtab(model_0_g5, model_1_g5)
  
  model_lm_g1 <- d_cl %>% 
    lm(form_0_g1,
       data = .)
  model_lm_g2 <- d_cl %>% 
    lm(form_0_g2,
       data = .)
  model_lm_g3 <- d_cl %>% 
    lm(form_0_g3,
       data = .)
  model_lm_g4 <- d_cl %>% 
    lm(form_0_g4,
       data = .)
  model_lm_g5 <- d_cl %>% 
    lm(form_0_g5,
       data = .)
  
  list(list(model_0_g1, model_1_g1, mod_comp_g1, model_lm_g1), 
       list(model_0_g2, model_1_g2, mod_comp_g2, model_lm_g2), 
       list(model_0_g3, model_1_g3, mod_comp_g3, model_lm_g3), 
       list(model_0_g4, model_1_g4, mod_comp_g4, model_lm_g4), 
       list(model_0_g5, model_1_g5, mod_comp_g5, model_lm_g5))
}

ca_gg_single_mods <- lapply(1:5, FUN = fit_gg_single_models, intro_vars, data_df, "ca")
co_gg_single_mods <- lapply(1:5, FUN = fit_gg_single_models, intro_vars, data_df, "co")
dinfl_gg_single_mods <- lapply(1:5, FUN = fit_gg_single_models, intro_vars, data_df, "d_infl")
finfl_gg_single_mods <- lapply(1:5, FUN = fit_gg_single_models, intro_vars, data_df, "fl_infl")
hair_gg_single_mods <- lapply(1:5, FUN = fit_gg_single_models, intro_vars, data_df, "hair")
ltooth_gg_single_mods <- lapply(1:5, FUN = fit_gg_single_models, intro_vars, data_df, "ltooth")
propltooth_gg_single_mods <- lapply(1:5, FUN = fit_gg_single_models, intro_vars, data_df, "prop_ltooth")

names(ca_gg_single_mods) <- 
  names(co_gg_single_mods) <- 
  names(dinfl_gg_single_mods) <- 
  names(finfl_gg_single_mods) <-
  names(hair_gg_single_mods) <-
  names(ltooth_gg_single_mods) <- 
  names(propltooth_gg_single_mods) <- 
  clades

lapply(ca_gg_single_mods, FUN = \(x)lapply(x, FUN = \(y)y[[3]]))
lapply(co_gg_single_mods, FUN = \(x)lapply(x, FUN = \(y)y[[3]]))
lapply(dinfl_gg_single_mods, FUN = \(x)lapply(x, FUN = \(y)y[[3]]))
lapply(finfl_gg_single_mods, FUN = \(x)lapply(x, FUN = \(y)y[[3]]))
lapply(ltooth_gg_single_mods, FUN = \(x)lapply(x, FUN = \(y)y[[3]]))
lapply(propltooth_gg_single_mods, FUN = \(x)lapply(x, FUN = \(y)y[[3]]))


tab_model(lapply(ca_gg_single_mods[[1]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/ca_algarve.html")
tab_model(lapply(ca_gg_single_mods[[2]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/ca_cadiz.html")
tab_model(lapply(ca_gg_single_mods[[3]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/ca_doñana.html")
tab_model(lapply(ca_gg_single_mods[[4]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/ca_hercynian.html")
tab_model(lapply(ca_gg_single_mods[[5]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/ca_tetraploid.html")

tab_model(lapply(co_gg_single_mods[[1]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/co_algarve.html")
tab_model(lapply(co_gg_single_mods[[2]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/co_cadiz.html")
tab_model(lapply(co_gg_single_mods[[3]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/co_doñana.html")
tab_model(lapply(co_gg_single_mods[[4]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/co_hercynian.html")
tab_model(lapply(co_gg_single_mods[[5]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/co_tetraploid.html")

tab_model(lapply(dinfl_gg_single_mods[[1]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/dinfl_algarve.html")
tab_model(lapply(dinfl_gg_single_mods[[2]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/dinfl_cadiz.html")
tab_model(lapply(dinfl_gg_single_mods[[3]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/dinfl_doñana.html")
tab_model(lapply(dinfl_gg_single_mods[[4]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/dinfl_hercynian.html")
tab_model(lapply(dinfl_gg_single_mods[[5]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/dinfl_tetraploid.html")

tab_model(lapply(finfl_gg_single_mods[[1]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/finfl_algarve.html")
tab_model(lapply(finfl_gg_single_mods[[2]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/finfl_cadiz.html")
tab_model(lapply(finfl_gg_single_mods[[3]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/finfl_doñana.html")
tab_model(lapply(finfl_gg_single_mods[[4]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/finfl_hercynian.html")
tab_model(lapply(finfl_gg_single_mods[[5]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/finfl_tetraploid.html")

tab_model(lapply(hair_gg_single_mods[[1]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/hair_algarve.html")
tab_model(lapply(hair_gg_single_mods[[2]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/hair_cadiz.html")
tab_model(lapply(hair_gg_single_mods[[3]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/hair_doñana.html")
tab_model(lapply(hair_gg_single_mods[[4]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/hair_hercynian.html")
tab_model(lapply(hair_gg_single_mods[[5]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/hair_tetraploid.html")

tab_model(lapply(ltooth_gg_single_mods[[1]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/ltooth_algarve.html")
tab_model(lapply(ltooth_gg_single_mods[[2]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/ltooth_cadiz.html")
tab_model(lapply(ltooth_gg_single_mods[[3]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/ltooth_doñana.html")
tab_model(lapply(ltooth_gg_single_mods[[4]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/ltooth_hercynian.html")
tab_model(lapply(ltooth_gg_single_mods[[5]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/ltooth_tetraploid.html")

tab_model(lapply(propltooth_gg_single_mods[[1]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/propltooth_algarve.html")
tab_model(lapply(propltooth_gg_single_mods[[2]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/propltooth_cadiz.html")
tab_model(lapply(propltooth_gg_single_mods[[3]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/propltooth_doñana.html")
tab_model(lapply(propltooth_gg_single_mods[[4]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/propltooth_hercynian.html")
tab_model(lapply(propltooth_gg_single_mods[[5]], FUN = \(x)x[[4]]), 
          file = "outputs/ancestry_models/single_models/propltooth_tetraploid.html")




## 1. Helper to make and save Q–Q plots -----------------------------

qq_plot_model <- function(model, name,
                          out_dir = "outputs/diagnostics/QQplots_single",
                          is_glmm = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  res <- if (is_glmm) residuals(model, type = "pearson") else residuals(model)
  
  file <- file.path(out_dir, paste0("QQ_", name, ".png"))
  png(file, width = 2000, height = 2000, res = 300)
  qqnorm(res, main = paste("Q–Q plot:", name))
  qqline(res, col = "red")
  dev.off()
}

## 2. Wrap all univariate clade-level models -----------------------

gg_single_mods <- list(
  ca          = ca_gg_single_mods,
  co          = co_gg_single_mods,
  d_infl      = dinfl_gg_single_mods,
  fl_infl     = finfl_gg_single_mods,
  hair        = hair_gg_single_mods,
  ltooth      = ltooth_gg_single_mods,
  prop_ltooth = propltooth_gg_single_mods
)

# Genetic group labels in the same order used to build the models
intro_vars <- paste0("group", c(3, 5, 4, 2, 1))

## 3. Loop over trait × clade × genetic group ----------------------

invisible(
  lapply(names(gg_single_mods), function(trait) {
    trait_list <- gg_single_mods[[trait]]        # e.g. ca_gg_single_mods
    
    lapply(names(trait_list), function(clade) {  # e.g. "Algarve", "Cádiz", ...
      clade_list <- trait_list[[clade]]          # length 5: one per genetic group
      
      lapply(seq_along(clade_list), function(g_idx) {
        group_name <- intro_vars[g_idx]
        mod_set    <- clade_list[[g_idx]]
        
        # mod_set structure:
        # [[1]] = glmmTMB univariate, no Sex random effect
        # [[2]] = glmmTMB univariate, + (1 | Sex)
        # [[3]] = AIC comparison
        # [[4]] = lm univariate
        
        # glmmTMB with Sex random effect
        qq_plot_model(
          model   = mod_set[[2]],
          name    = paste(trait, clade, group_name, "glmm_SexRE", sep = "_"),
          is_glmm = TRUE
        )
        
        # linear model
        qq_plot_model(
          model   = mod_set[[4]],
          name    = paste(trait, clade, group_name, "lm_noSex", sep = "_"),
          is_glmm = FALSE
        )
      })
    })
  })
)


