# ======================================================================
# Script: 06-Variation_between_and_within_groups.R
# Purpose: Partition morphological variance within and among genetic groups using hierarchical models (including Bayesian variance components).
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
# =========================
# 0) Paquetes
# =========================

# ----------------------------------------------------------------------
# Package requirements
# ----------------------------------------------------------------------
# Load all R packages required for data import, manipulation, modelling,
# and figure/table generation.
library(tidyverse)
library(openxlsx)
library(brms)
library(lme4)
library(permuco)
library(posterior)

# Ensure reproducible stochastic procedures (e.g. permutations, MCMC).
set.seed(123)  # reproducibilidad

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

morph_measures <- morph_measures %>% 
  mutate(across(c(Population_ID,
                  Individual_ID,
                  Sex, 
                  clade_2,
                  clade_5), 
                as.factor)) %>% 
  mutate(ca_log = log(ca))


# ----------------------------------------------------------------------
# BRMS approach
# ----------------------------------------------------------------------

priors <- c(
  prior(normal(0, 5), class = "Intercept"),
  prior(exponential(1), class = "sd"),      # all random effects
  prior(exponential(1), class = "sigma")    # residual SD
)


b <- brm(
  ca_log ~ 1 + (1 | clade_5 / Population_ID / Individual_ID),
  data   = morph_measures,
  family = gaussian(),
  prior  = priors,
  chains = 4,
  cores  = 8,
  iter   = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.95)
)


# Inspect object structure and summary statistics.
summary(b)
summary(b)$fixed
summary(b)$random

# Check model
plot(b)

pp_check(b, ndraws = 100)

# Extract parameter estimates and variance components.
vc <- VarCorr(b)
vc


# Extract variance components
draws <- as_draws_df(b)

var_draws <- draws %>%
  transmute(
    var_clade = `sd_clade_5__Intercept`^2,
    var_pop   = `sd_clade_5:Population_ID__Intercept`^2,
    var_ind   = `sd_clade_5:Population_ID:Individual_ID__Intercept`^2,
    var_res   = sigma^2
  ) %>%
  mutate(
    var_total = var_clade + var_pop + var_ind + var_res,
    prop_clade = var_clade / var_total,
    prop_pop   = var_pop   / var_total,
    prop_ind   = var_ind   / var_total,
    prop_res   = var_res   / var_total
  )

summary_props <- var_draws %>%
  summarise(
    clade_mean_var = mean(var_clade),
    clade_lwr_var  = quantile(var_clade, 0.025),
    clade_upr_var  = quantile(var_clade, 0.975),
    
    pop_mean_var = mean(var_pop),
    pop_lwr_var  = quantile(var_pop, 0.025),
    pop_upr_var  = quantile(var_pop, 0.975),
    
    ind_mean_var = mean(var_ind),
    ind_lwr_var  = quantile(var_ind, 0.025),
    ind_upr_var  = quantile(var_ind, 0.975),
    
    res_mean_var = mean(var_res),
    res_lwr_var  = quantile(var_res, 0.025),
    res_upr_var  = quantile(var_res, 0.975),
    
    total_mean_var = mean(var_total),
    total_lwr_var  = quantile(var_total, 0.025),
    total_upr_var  = quantile(var_total, 0.975),

    clade_mean_prop = mean(prop_clade),
    clade_lwr_prop  = quantile(prop_clade, 0.025),
    clade_upr_prop  = quantile(prop_clade, 0.975),
    
    pop_mean_prop = mean(prop_pop),
    pop_lwr_prop  = quantile(prop_pop, 0.025),
    pop_upr_prop  = quantile(prop_pop, 0.975),
    
    ind_mean_prop = mean(prop_ind),
    ind_lwr_prop  = quantile(prop_ind, 0.025),
    ind_upr_prop  = quantile(prop_ind, 0.975),
    
    res_mean_prop = mean(prop_res),
    res_lwr_prop  = quantile(prop_res, 0.025),
    res_upr_prop  = quantile(prop_res, 0.975)
  )

summary_props


var_draws %>%
  select(prop_clade, prop_pop, prop_ind, prop_res) %>%
  pivot_longer(everything(), names_to = "component", values_to = "prop") %>%
  ggplot(aes(x = prop)) +
  geom_density() +
  facet_wrap(~ component, ncol = 2) +
  xlab("Proportion of total variance") +
  ylab("Posterior density")

var_draws <- var_draws %>%
  mutate(
    prop_inter = prop_clade + prop_pop,
    prop_intra = prop_ind + prop_res
  )

var_draws %>%
  summarise(
    inter_mean = mean(prop_inter),
    inter_lwr  = quantile(prop_inter, 0.025),
    inter_upr  = quantile(prop_inter, 0.975),
    intra_mean = mean(prop_intra),
    intra_lwr  = quantile(prop_intra, 0.025),
    intra_upr  = quantile(prop_intra, 0.975)
  )




fit_one_clade <- function(dat_clade) {
  b_cl <- brm(
    ca ~ 1 + (1 | Population_ID / Individual_ID),
    data = dat_clade,
    family = gaussian(),
    prior  = priors,
    chains = 4, 
    cores = 4,
    iter = 4000, 
    warmup = 1000,
    control = ctrl,
    silent = TRUE, refresh = 0
  )
  
  dr <- as_draws_df(b_cl)
  
  sd_names <- grep("^sd_", names(dr), value = TRUE)
  
  pop_sd <- sd_names[grepl("^sd_Population_ID__Intercept$", sd_names)]
  ind_sd <- sd_names[grepl("^sd_Population_ID:Individual_ID__Intercept$", sd_names)]
  
  if (length(pop_sd) != 1 || length(ind_sd) != 1) {
    stop("Could not uniquely identify population/individual SD columns. Inspect sd_names:\n",
         paste(sd_names, collapse = "\n"))
  }
  
  out <- dr %>%
    transmute(
      var_pop = (.data[[pop_sd]])^2,
      var_ind = (.data[[ind_sd]])^2,
      var_res = sigma^2
    ) %>%
    mutate(
      var_total = var_pop + var_ind + var_res,
      prop_pop  = var_pop / var_total,
      prop_ind  = var_ind / var_total,
      prop_res  = var_res / var_total
    )
  
  list(fit = b_cl, draws = out)
}


clades <- levels(morph_measures$clade_5)

# Controls
ctrl <- list(adapt_delta = 0.95)

fits <- map(clades, ~ fit_one_clade(filter(morph_measures, clade_5 == .x)))
names(fits) <- clades

posterior_by_clade <- imap_dfr(
  fits,
  ~ mutate(.x$draws, clade_5 = .y, .before = 1)
)

# Quick summaries (posterior mean and 95% credible interval)
summary_by_clade <- posterior_by_clade %>%
  pivot_longer(cols = starts_with("var_") | starts_with("prop_"),
               names_to = "component", values_to = "value") %>%
  group_by(clade_5, component) %>%
  summarise(
    mean = mean(value),
    lwr  = quantile(value, 0.025),
    upr  = quantile(value, 0.975),
    .groups = "drop"
  )

summary_by_clade


