# =========================
# 0) Paquetes
# =========================

require(tidyverse)
library(openxlsx)
# library(vegan)
# library(effectsize)
# library(performance)
# library(DHARMa)
library(brms)
library(lme4)
library(permuco)
library(posterior)

set.seed(123)  # reproducibilidad

# =========================
# 1) Importar Excel
# =========================
#### LOAD MORPHOMETRIC DATA  ####
genetic_groups <- "data/genetic groups.xlsx" %>% 
  openxlsx::read.xlsx(2) %>% 
  mutate(genetic_group_k2 = case_match(genetic_group_k2,
                                       "tetraploid" ~ "4x",
                                       "diploid" ~ "2x")) %>% 
  mutate(genetic_group_k5 = str_to_sentence(genetic_group_k5))



#### LOAD MORPHOMETRIC DATA  ####
morph_measures <- "data/Morphometric_measures.xlsx" %>% 
  openxlsx::read.xlsx(2)

head(morph_measures)  
summary(morph_measures)


morph_measures <- morph_measures %>% 
  left_join(genetic_groups) %>% 
  rename(clade_5 = genetic_group_k5) %>% 
  rename(clade_2 = genetic_group_k2)

# Create column with ploidy level
morph_measures <- morph_measures %>% 
  mutate(across(c(Population_ID, # Change to factors.
                  Individual_ID,
                  Sex, 
                  clade_2,
                  clade_5), 
                as.factor)) %>% 
  mutate(ca_log = log(ca))


# measures_count <- morph_measures %>% count(clade_5, Population_ID, Individual_ID) %>% arrange(n)
# 
# selected_individuals <- measures_count %>% filter(n > 2) %>% pull(Individual_ID)
# 
# morph_measures <- morph_measures %>% filter(Individual_ID %in% selected_individuals)














#### BRMS approach

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


summary(b)
summary(b)$fixed
summary(b)$random


plot(b)

pp_check(b, ndraws = 100)

vc <- VarCorr(b)
vc


draws <- as_draws_df(b)

# Extract variance components
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
  # dat_clade <- morph_measures %>% filter(clade_5 == "Algarve")
  
  # Fit model: variance among populations + among individuals within populations
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
  
  # Extract SD names programmatically (robust to naming order)
  sd_names <- grep("^sd_", names(dr), value = TRUE)
  
  # Identify population and individual SD columns
  # Typically: sd_Population_ID__Intercept and sd_Individual_ID:Population_ID__Intercept
  pop_sd <- sd_names[grepl("^sd_Population_ID__Intercept$", sd_names)]
  ind_sd <- sd_names[grepl("^sd_Population_ID:Individual_ID__Intercept$", sd_names)]
  
  if (length(pop_sd) != 1 || length(ind_sd) != 1) {
    stop("Could not uniquely identify population/individual SD columns. Inspect sd_names:\n",
         paste(sd_names, collapse = "\n"))
  }
  
  # Posterior draws of variance components + proportions
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





















# =========================
# 2) Preparación de datos
# =========================



# Selección de variables florales numéricas
# Si hay columnas numéricas que NO son florales (p.ej. códigos), exclúyelas manualmente:
X <- morph_measures %>% select(ca)




##################



m <- lmer(ca ~ 1 + (1 | clade_5 / Population_ID / Individual_ID), data = morph_measures)
VarCorr(m)

isSingular(m, tol = 1e-5)


vc <- as.data.frame(VarCorr(m))


vc$prop <- vc$vcov / sum(vc$vcov)
vc

confint(m, method = "profile")

confint(m, method = "boot", nsim = 1000)




check_model(m)         # overall diagnostic plots
check_singularity(m)   # warns about boundary/singular random-effects fits

sim_res <- simulateResiduals(fittedModel = m, n = 1000)
plot(sim_res)                       # overall diagnostic panel
testUniformity(sim_res)             # general misspecification
testDispersion(sim_res)             # over/under-dispersion (more relevant for GLMMs, still informative)
testOutliers(sim_res)               # outlier test


re <- ranef(m, condVar = TRUE)

# QQ plot of random intercepts per grouping level
for (lvl in names(re)) {
  qqnorm(re[[lvl]][[1]], main = paste("QQ random effects:", lvl))
  qqline(re[[lvl]][[1]])
}



infl_pop <- influence(m, group = "clade_5")
plot(infl_pop, which = "cook")

infl_pop <- influence(m, group = "clade_5/Population_ID")
plot(infl_pop, which = "cook")



# =========================
# 3) ANÁLISIS MULTIVARIADO (todas las variables a la vez)
#    PERMANOVA: varianza entre vs dentro (por permutación)
# =========================
# Distancia Euclídea (típica para medidas continuas). Alternativa: "gower" si hay mixtas.
# Recomendación: estandarizar si están en escalas muy diferentes:
# X_scaled <- scale(X)

dist_mat <- dist(X,
                 method = "euclidean")

# permanova <- adonis2(dist_mat ~ clade_5, 
#                      data = morph_measures,
#                      permutations = 9999,
#                      by = "margin")

permanova <- adonis2(dist_mat ~ Population_ID, 
                     data = morph_measures,
                     permutations = 999,
                     by = "margin")
 
print(permanova)

# "R2" en adonis2 es la proporción de variación explicada (entre grupos / total).
# El resto (1 - R2) puede interpretarse como variación residual (dentro de grupos + error).

# =========================
# 4) VARIACIÓN DENTRO DE GRUPOS (dispersión intragrupo)
#    betadisper evalúa si los grupos difieren en su dispersión (variación interna)
# =========================
# bd <- betadisper(dist_mat,
#                  group = morph_measures$clade_5, 
#                  type = "centroid")

bd <- betadisper(dist_mat,
                 group = morph_measures$Population_ID, 
                 type = "centroid")

anova_bd <- anova(bd)                 # test permutacional/anova sobre dispersión

perm_bd  <- permutest(bd,
                      permutations = 9999)

print(anova_bd)
print(perm_bd)

# Resumen por grupo: distancia media al centroide (medida de variación dentro de grupo)
disp_por_grupo <- tapply(bd$distances, 
                         morph_measures$Population_ID, 
                         mean)
print(disp_por_grupo)



TukeyHSD(bd)




