## Dog culling in Tanzania: quantitative analysis
##
## Code for: "Any dog that was roaming around... we deal with it": investigating
## the scale, organisation, implementation and motivations behind dog culling in
## a rabies-endemic setting
##
## The survey data are not public, as they could identify the LGAs and officers
## who took part. They are available on request (see README). Run this script
## from the project folder, with the data in data_shared/.
##
## Files
##   data_shared/occurrence.csv      survey of occurrence (49 LGAs)
##   data_shared/practice.csv        survey of practice (16 LGAs)
##   data_shared/councils_2022.csv   all LGAs, 2022 census (194 units)
##   code/functions/trunc_dist_funs.R  zero-truncated distribution functions
##
## Sections
##   1. Data
##   2. Survey of occurrence: descriptives, logistic regression, power analysis
##   3. Survey of practice: descriptives, scale of culling, Figure 1
##   4. National estimate: events distribution, cull size model, simulation,
##      sensitivity analyses
##   5. Surveyed vs not surveyed LGAs
##   6. Session info
##
## Where each output is used
##   Table 2              2a, 2b
##   Table 3              3a
##   Figure 1             3d
##   Table S2             2d
##   Table S3             5b, 5c
##   Table S4             4b (a), 4c (b)
##   Figure S1            4b (a), 4c (b)
##   Table S5             4d
##   Methods S1           3c, 4c, 4d
##   Results 3.1.1        1, 2a-2c, 5b
##   Results 3.1.2        3b, 4b-4d
##   Supp File 2, model checks   2b (check)
##   Not included here (use named data): Table S1, Table 3 funding source,
##   number of regions with culling


library(brms)
library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(fitdistrplus)
library(truncdist)
library(logistf)
library(ggplot2)
library(patchwork)
library(DHARMa)

# Zero-truncated negative binomial and Poisson functions
source("code/functions/trunc_dist_funs.R")

#### 1. Data ####

## 1a. Survey of occurrence
# 49 LGAs (51 responses: one from a national park authority and one duplicate
# removed; one response covering a district since split into two councils)
occ <- read_csv("data_shared/occurrence.csv", show_col_types = FALSE) %>%
  mutate(
    dpm_text        = tolower(dpm_activity_12m),
    culled_12m      = str_detect(dpm_text, "killing"),
    culled          = as.integer(culled_12m),
    settlement_type = tolower(settlement_type),
    urban           = settlement_type == "urban",
    hdr             = human_pop / dog_pop_2022,
    log2_human_pop  = log2(human_pop),
    log2_dog_pop    = log2(dog_pop_2022),
    rabies_dog      = tolower(rabies_dog_12m) == "yes",
    rabies_human    = tolower(rabies_human_12m) == "yes",
    rabies_problem  = tolower(rabies_perceived_problem) == "yes",
    dog_vaccination = tolower(dog_vaccination_12m) == "yes",
    funding_cull    = tolower(dpm_killing_allocated_funding_12m) == "yes",
    dpm_strategy = case_when(
      str_detect(dpm_text, "killing") &
        str_detect(dpm_text, "spay|tether|education|impound") ~ "Culling + other methods",
      str_detect(dpm_text, "killing") ~ "Only culling",
      dpm_text == "nothing"           ~ "No activity",
      TRUE                            ~ "Non-lethal methods only"
    )
  )
stopifnot(nrow(occ) == 49)
table(occ$culled_12m)
table(occ$settlement_type, occ$culled_12m)

## 1b. Survey of practice
# 16 LGAs; prac_complete = the 15 with cull size and frequency (used in 3b, 4)
prac <- read_csv("data_shared/practice.csv", show_col_types = FALSE) %>%
  mutate(
    Non_specific_Number_culled = tolower(Non_specific_Number_culled),
    annual_total_culled = Number_culled_1Event * Frequency_12m,
    pct_pop_culled      = 100 * annual_total_culled / Dog_population_2022
  )

prac_complete <- prac %>%
  filter(!is.na(Number_culled_1Event),
         !is.na(Frequency_12m),
         !is.na(Non_specific_Number_culled))

## 1c. All LGAs, 2022
# 195 LGAs in the 2022 census, as 194 units (the two councils covered by one
# survey response are combined): 140 rural, 54 urban
# Human population and households: 2022 Population and Housing Census
# Dog population: LGA human:dog ratios applied to 2022 census populations
# Zones: 2022 TDHS-MIS, grouped into North, East, South, West, Central
councils_2022 <- read_csv("data_shared/councils_2022.csv", show_col_types = FALSE)
stopifnot(nrow(councils_2022) == 194)

nat_dog_pop <- sum(councils_2022$dog_pop_2022)   # national dog population
nat_dog_pop

N_rural     <- sum(councils_2022$urban_rural == "rural")
N_urban     <- sum(councils_2022$urban_rural == "urban")
N_districts <- N_rural + N_urban

## 1d. Missing data
# -> Supp File 2, missing data paragraph
made_here <- c("dpm_text", "culled_12m", "culled", "urban", "hdr", "log2_human_pop",
               "log2_dog_pop", "rabies_dog", "rabies_human", "rabies_problem",
               "dog_vaccination", "funding_cull", "dpm_strategy",
               "annual_total_culled", "pct_pop_culled")
occ_q  <- dplyr::select(occ,  -any_of(made_here))
prac_q <- dplyr::select(prac, -any_of(made_here))
missing_data <- bind_rows(
  tibble(survey = "Survey of occurrence", variable = names(occ_q),
         n_missing = colSums(is.na(occ_q)), n = nrow(occ_q)),
  tibble(survey = "Survey of practice", variable = names(prac_q),
         n_missing = colSums(is.na(prac_q)), n = nrow(prac_q))
) %>% filter(n_missing > 0)
missing_data


#### 2. Survey of occurrence ####

## 2a. Descriptives (Table 2)
# -> Table 2, culling / no culling columns
# Populations as median (IQR), rounded, so no single LGA is identifiable
med_iqr2 <- function(x, dp = 0, round_to = NA) {
  q <- quantile(x, c(0.5, 0.25, 0.75), na.rm = TRUE)
  if (!is.na(round_to)) q <- round(q / round_to) * round_to
  f <- function(v) formatC(v, format = "f", digits = dp, big.mark = ",")
  sprintf("%s (%s-%s)", f(q[1]), f(q[2]), f(q[3]))
}
n_pct <- function(x) sprintf("%d (%.0f%%)", sum(x, na.rm = TRUE),
                             100 * mean(x, na.rm = TRUE))

occ %>% group_by(culled_12m) %>%
  summarise(n           = n(),
            human_pop   = med_iqr2(human_pop, round_to = 1000),
            dog_pop     = med_iqr2(dog_pop_2022, round_to = 100),
            hdr         = med_iqr2(hdr, 1),
            urban       = n_pct(urban),
            rabies_dog  = n_pct(rabies_dog),
            rabies_hum  = n_pct(rabies_human),
            rabies_prob = n_pct(rabies_problem),
            vacc        = n_pct(dog_vaccination)) %>%
  print(width = Inf)

## 2b. Univariable logistic regression (Table 2)
# -> Table 2 OR and p columns; Results 3.1.1 (urban OR 5.90, p = 0.017)
predictors <- c(urban           = "Settlement type: urban",
                log2_human_pop  = "Human population (per doubling)",
                log2_dog_pop    = "Dog population (per doubling)",
                hdr             = "Human:dog ratio",
                rabies_dog      = "Rabies reported in dogs",
                rabies_human    = "Rabies reported in humans",
                rabies_problem  = "Rabies perceived as a problem",
                dog_vaccination = "Dog vaccination conducted")

# Firth's penalised regression for complete separation
firth_vars <- c("rabies_problem")

fit_one <- function(v) {
  f <- reformulate(v, response = "culled")
  if (v %in% firth_vars) {
    m <- logistf(f, data = occ)
    term <- names(coef(m))[2]
    data.frame(variable = predictors[[v]], method = "Firth",
               OR = exp(coef(m)[term]),
               lo = exp(m$ci.lower[term]), hi = exp(m$ci.upper[term]),
               p  = m$prob[term])
  } else {
    m <- glm(f, data = occ, family = binomial)
    ci <- exp(confint(m))   # profile likelihood
    data.frame(variable = predictors[[v]], method = "glm",
               OR = exp(coef(m)[2]), lo = ci[2, 1], hi = ci[2, 2],
               p  = summary(m)$coefficients[2, 4])
  }
}

uni <- bind_rows(lapply(names(predictors), fit_one)) %>%
  mutate(term  = names(predictors),
         OR_ci = sprintf("%.2f (%.2f - %.2f)", OR, lo, hi),
         p     = signif(p, 2))
rownames(uni) <- NULL
uni

# 46/49 (Table 2 footnote, Results 3.1.1)
sum(occ$rabies_problem); nrow(occ)

# Other methods used alongside culling
# -> Results 3.1.1 para 2 (tethering 6/23 etc.)
occ %>% filter(culled_12m) %>%
  summarise(n = n(),
            tethering   = sum(str_detect(dpm_text, "tether")),
            spay_neuter = sum(str_detect(dpm_text, "spay")),
            education   = sum(str_detect(dpm_text, "education")),
            impounding  = sum(str_detect(dpm_text, "impound")))

table2_models <- uni %>% dplyr::select(variable, method, OR_ci, p)
table2_models

## 2b (check). Residual checks for Table 2 models (DHARMa)
# -> Supp File 2, model checks paragraph; Methods 2.2.1
# Yes/no outcome, so dispersion and zero-inflation tests don't apply.
# Checks: overall fit (uniformity), outliers, and fit across the range
# of continuous predictors. Firth model (rabies_problem) not supported.
dharma_check <- function(v) {
  m <- glm(reformulate(v, response = "culled"), data = occ, family = binomial)
  s <- simulateResiduals(m, n = 1000, seed = 1)
  tibble(variable   = predictors[[v]],
         uniformity = testUniformity(s, plot = FALSE)$p.value,
         outliers   = testOutliers(s, plot = FALSE)$p.value,
         form       = if (is.numeric(occ[[v]]))
           testQuantiles(s, predictor = occ[[v]], plot = FALSE)$p.value
         else NA_real_)
}
dharma_table <- bind_rows(lapply(setdiff(names(predictors), firth_vars), dharma_check)) %>%
  mutate(across(-variable, ~ signif(.x, 2)))
dharma_table

# Dog population: is the relationship with culling non-linear? (form test p < 0.05)
# -> Supp File 2, model checks (quadratic LRT p = 0.006; 47%, 56%, 38% by thirds)
m_dog  <- glm(culled ~ log2_dog_pop, data = occ, family = binomial)
m_dog2 <- glm(culled ~ log2_dog_pop + I(log2_dog_pop^2), data = occ, family = binomial)
anova(m_dog, m_dog2, test = "LRT")      # p < 0.05 = curved fits better

# Culling by dog population thirds
table(dog_third = cut(occ$dog_pop_2022, quantile(occ$dog_pop_2022, 0:3/3),
                      include.lowest = TRUE, dig.lab = 6),
      culled = occ$culled_12m)

# Plot of residuals against dog population
plotResiduals(simulateResiduals(m_dog, n = 1000, seed = 1), form = occ$log2_dog_pop)

## 2c. Multivariable model (variables with p < 0.05)
# only urban < 0.05 so not fitted (Results 3.1.1)
keep <- uni$term[uni$p < 0.05 & !(uni$term %in% firth_vars)]
keep
if (length(keep) > 1) {
  multi <- glm(reformulate(keep, response = "culled"), data = occ, family = binomial)
  cbind(OR = exp(coef(multi)), exp(confint(multi)))[-1, ]
}

# DPM strategy and culling funding determine the outcome, so descriptive only
# -> Results 3.1.1 (70%/30%, 69%/31%, funding 9/23)
table(occ$dpm_strategy, occ$culled_12m)
table(funding = tolower(occ$dpm_killing_allocated_funding_12m) == "yes",
      culled_12m = occ$culled_12m)


## 2d. Sensitivity power analysis
# -> Table S2. Slow (a few minutes)
# Smallest OR detectable with 80% power (alpha = 0.05, likelihood ratio test), by simulation.
# Uses the observed values of each explanatory variable; only the outcome is simulated
power_for_or <- function(or, x, p_cull = mean(occ$culled_12m), n_sim = 1000, seed = 1) {
  set.seed(seed)
  x_c <- x - mean(x)
  hits <- replicate(n_sim, {
    y <- rbinom(length(x), 1, plogis(qlogis(p_cull) + log(or) * x_c))
    if (length(unique(y)) < 2) return(FALSE)
    m0 <- glm(y ~ 1, family = binomial)
    m1 <- suppressWarnings(glm(y ~ x_c, family = binomial))
    isTRUE(anova(m0, m1, test = "LRT")$`Pr(>Chi)`[2] < 0.05)   # NA (no difference) = not significant
  })
  mean(hits)
}

# Stops at the first OR reaching 80% power (faster; same result)
min_detectable_or <- function(x, ors) {
  for (or in ors) if (power_for_or(or, x) >= 0.8) return(or)
  NA
}

binary_mdor <- sapply(c("urban", "rabies_dog", "rabies_human", "dog_vaccination"),
                      function(v) min_detectable_or(as.numeric(occ[[v]]), seq(1.5, 30, by = 0.5)))

cont_mdor <- c(
  sapply(c("log2_human_pop", "log2_dog_pop"),
         function(v) min_detectable_or(occ[[v]], seq(1.1, 10, by = 0.1))),
  hdr = min_detectable_or(occ$hdr, seq(1.01, 2, by = 0.01)))   # finer grid: per 1-unit change

# NA: 80% power not reached within the range tested
power_table <- tibble(variable = names(c(binary_mdor, cont_mdor)),
                      min_detectable_OR = c(binary_mdor, cont_mdor))
power_table


#### 3. Survey of practice ####

## 3a. Descriptives (Table 3)
# -> Table 3 (funding source grouped from free text, not included)
yn_cols <- names(prac)[sapply(prac, function(x) {
  v <- unique(na.omit(as.character(x)))
  length(v) > 0 && all(v %in% c("Yes", "No"))
})]

table3 <- tibble(
  variable   = yn_cols,
  n_yes      = sapply(yn_cols, function(v) sum(prac[[v]] == "Yes", na.rm = TRUE)),
  n_answered = sapply(yn_cols, function(v) sum(!is.na(prac[[v]]))),
  pct_of_16  = round(100 * n_yes / nrow(prac))
)
print(table3, n = Inf)

table(prac$Rabies_DPM, useNA = "ifany")
table(prac$Method_Shoot_Poison, useNA = "ifany")
table(prac$Who, useNA = "ifany")

## 3b. Scale of culling
# -> Results 3.1.2 para 1 (15/16, 1-23 events, 9-230 dogs, 7 minimum, 0.11-5.24%)
nrow(prac_complete)
range(prac_complete$Frequency_12m)
range(prac_complete$Number_culled_1Event)
sum(prac_complete$Non_specific_Number_culled == "more")
# range of counts reported as exact -> Methods 2.2.1, Methods S1 (9-230)
range(prac_complete$Number_culled_1Event[prac_complete$Non_specific_Number_culled != "more"])
round(range(prac_complete$pct_pop_culled), 2)

## 3c. Association between event frequency and cull size
# -> Methods S1 (association between events and cull size). p_corr not in paper
assoc_nb   <- MASS::glm.nb(Number_culled_1Event ~ log2(Frequency_12m), data = prac_complete)
assoc_null <- MASS::glm.nb(Number_culled_1Event ~ 1, data = prac_complete)
exp(cbind(RR = coef(assoc_nb), confint(assoc_nb)))[2, ]   # rate ratio per doubling of events
anova(assoc_null, assoc_nb)                                # likelihood ratio test p-value

# Same question, different models -> Methods S1 (p = 0.08-0.44; Poisson p < 0.001)
d <- prac_complete %>%
  mutate(events  = Frequency_12m,
         size    = Number_culled_1Event,
         minimum = Non_specific_Number_culled == "more")

glm_row <- function(label, m, m0, term) {
  ci <- suppressMessages(confint(m))[term, ]
  tibble(model = label,
         estimate = exp(coef(m)[[term]]),
         lo = exp(ci[[1]]), hi = exp(ci[[2]]),
         p  = anova(m0, m, test = "Chisq")[2, ncol(anova(m0, m, test = "Chisq"))])
}

glm_compare <- bind_rows(
  glm_row("1 NB: size ~ log2(events)",
          MASS::glm.nb(size ~ log2(events), data = d),
          MASS::glm.nb(size ~ 1, data = d), "log2(events)"),
  glm_row("2 NB: size ~ events",
          MASS::glm.nb(size ~ events, data = d),
          MASS::glm.nb(size ~ 1, data = d), "events"),
  glm_row("3 NB: size ~ log2(events) + minimum",
          MASS::glm.nb(size ~ log2(events) + minimum, data = d),
          MASS::glm.nb(size ~ minimum, data = d), "log2(events)"),
  glm_row("4 NB: events ~ log2(size)",
          MASS::glm.nb(events ~ log2(size), data = d),
          MASS::glm.nb(events ~ 1, data = d), "log2(size)"),
  glm_row("5 Linear: log(size) ~ log2(events)",
          glm(log(size) ~ log2(events), data = d),
          glm(log(size) ~ 1, data = d), "log2(events)"),
  glm_row("6 Poisson: size ~ log2(events)",
          glm(size ~ log2(events), family = poisson, data = d),
          glm(size ~ 1, family = poisson, data = d), "log2(events)")
) %>%
  mutate(across(c(estimate, lo, hi), ~ round(.x, 2)), p = signif(p, 2))
glm_compare

p_corr <- ggplot(prac_complete,
                 aes(Frequency_12m, Number_culled_1Event,
                     shape = Non_specific_Number_culled)) +
  geom_point(size = 2.5) +
  scale_shape_manual(values = c(equal = 16, more = 2),
                     labels = c(equal = "Exact count", more = "Minimum count"),
                     name = NULL) +
  labs(x = "Culling events in previous 12 months", y = "Dogs culled per event") +
  theme_minimal(base_size = 11)
p_corr

## 3d. Figure 1
# -> Figure 1
fig_dat <- prac_complete %>%
  arrange(Number_culled_1Event) %>%
  mutate(plot_id = factor(row_number(), levels = row_number()))

theme_fig1 <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4)
  )

# A: dogs culled per event ("+" = minimum count)
p_top <- ggplot(fig_dat, aes(x = plot_id, y = Number_culled_1Event)) +
  geom_col(fill = "grey80", colour = "black", linewidth = 0.3) +
  geom_text(data = filter(fig_dat, Non_specific_Number_culled == "more"),
            label = "+", vjust = -0.3, size = 4) +
  scale_y_continuous(breaks = seq(0, 300, 40), expand = expansion(mult = c(0, 0.1))) +
  labs(y = "Dogs culled per event") +
  theme_fig1

# B: events per year
p_mid <- ggplot(fig_dat, aes(x = plot_id, y = Frequency_12m)) +
  geom_col(fill = "grey55", colour = "black", linewidth = 0.3) +
  scale_y_continuous(breaks = seq(0, 24, 4), expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Events per year") +
  theme_fig1

# C: percentage of LGA dog population culled per year
p_bottom <- ggplot(fig_dat, aes(x = plot_id, y = pct_pop_culled)) +
  geom_col(fill = "grey15", colour = "black", linewidth = 0.3) +
  scale_y_continuous(limits = c(0, 6), breaks = 0:6, expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Percentage of LGA dog\npopulation culled annually (%)",
       x = "LGA (anonymised)") +
  theme_fig1 +
  theme(axis.text.x = element_text(size = 8),
        axis.ticks.x = element_line(),
        axis.title.x = element_text(size = 10))

fig1 <- p_top / p_mid / p_bottom + plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = "A")
fig1


#### 4. National estimate ####
# Each simulated year: each LGA culls with probability p_rural or p_urban;
# each culling LGA gets a number of events (shifted negative binomial, 4b)
# and one cull size from the Bayesian cull size model (4c), used for all its events.
# Uncertainty: 1,000 bootstrap resamples x 100 simulated years; median and
# 2.5th-97.5th percentiles reported.

## 4a. Inputs
# rural 13/36, urban 10/13 (Methods S1)
cull_size <- prac_complete$Number_culled_1Event
events    <- prac_complete$Frequency_12m
is_min    <- prac_complete$Non_specific_Number_culled == "more"   # minimum ("more than X") counts

p_by_stratum <- occ %>%
  group_by(settlement_type) %>%
  summarise(n = n(), n_culling = sum(culled_12m), .groups = "drop") %>%
  mutate(Pcull = n_culling / n)
p_by_stratum

Pcull_rural      <- p_by_stratum$Pcull[p_by_stratum$settlement_type == "rural"]
Pcull_urban      <- p_by_stratum$Pcull[p_by_stratum$settlement_type == "urban"]
n_surveyed_rural <- p_by_stratum$n[p_by_stratum$settlement_type == "rural"]
n_surveyed_urban <- p_by_stratum$n[p_by_stratum$settlement_type == "urban"]

## 4b. Events per LGA: fitted distributions (Table S4a, Figure S1a)
# Shifted negative binomial (main model) / shifted Poisson: fitted to count - 1
# Zero-truncated negative binomial: compared by AIC and used as a sensitivity analysis
# mean/var -> Methods 2.2.1 (5.2, 29.3); AICs -> Results 3.1.2
fit_events_dist <- function(ev, dist = "trunc_nbinom") {
  switch(dist,
         trunc_nbinom = fitdist(ev, "trunc_nbinom", start = list(size = 1, mu = 2), discrete = TRUE),
         shifted_nb   = fitdist(ev - 1, "nbinom", discrete = TRUE),
         shifted_pois = fitdist(ev - 1, "pois", discrete = TRUE))
}

r_from_fit <- function(fit, dist = "trunc_nbinom") {
  e <- fit$estimate
  switch(dist,
         trunc_nbinom = function(n) rtrunc_nbinom(n, size = e[["size"]], mu = e[["mu"]]),
         shifted_nb   = function(n) rnbinom(n, size = e[["size"]], mu = e[["mu"]]) + 1,
         shifted_pois = function(n) rpois(n, lambda = e[["lambda"]]) + 1)
}

mean(events); var(events)

dists <- c("trunc_nbinom", "shifted_nb", "shifted_pois")
fits_events <- setNames(lapply(dists, function(d) fit_events_dist(events, d)), dists)

table_s4a <- bind_rows(lapply(dists, function(d) {
  f <- fits_events[[d]]
  tibble(dist      = d,
         parameter = names(f$estimate),
         estimate  = unname(f$estimate),
         se        = unname(f$sd),
         loglik    = f$loglik,
         AIC       = f$aic,
         converged = f$convergence == 0)
})) %>%
  mutate(across(c(estimate, se, loglik, AIC), ~ round(.x, 2)))
table_s4a

# Figure S1 (a): observed vs fitted cumulative distributions of events
cdf_fitted <- function(fits, x) {
  e <- lapply(fits, `[[`, "estimate")
  bind_rows(
    tibble(x = x, dist = "Zero-truncated negative binomial",
           F = ptrunc_nbinom(x, size = e$trunc_nbinom[["size"]], mu = e$trunc_nbinom[["mu"]])),
    tibble(x = x, dist = "Shifted negative binomial",
           F = pnbinom(x - 1, size = e$shifted_nb[["size"]], mu = e$shifted_nb[["mu"]])),
    tibble(x = x, dist = "Shifted Poisson",
           F = ppois(x - 1, lambda = e$shifted_pois[["lambda"]])))
}
dist_lines <- scale_linetype_manual(
  values = c("Shifted negative binomial"        = "solid",    # main model
             "Zero-truncated negative binomial" = "dashed",
             "Shifted Poisson"                  = "dotted"),
  breaks = c("Shifted negative binomial", "Zero-truncated negative binomial", "Shifted Poisson"),
  name = NULL)
theme_figS1 <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4))

x_ev   <- 0:max(events)
obs_ev <- tibble(x = x_ev, F = ecdf(events)(x_ev))
p_ev <- ggplot() +
  geom_step(data = obs_ev, aes(x, F, colour = "Observed"), linewidth = 1.2) +
  geom_step(data = cdf_fitted(fits_events, x_ev), aes(x, F, linetype = dist)) +
  scale_colour_manual(values = c("Observed" = "grey60"), name = NULL) +
  dist_lines +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Culling events per LGA in previous 12 months", y = "Cumulative proportion") +
  theme_figS1 +
  theme(legend.position = "top", legend.box = "vertical",
        legend.spacing.y = unit(0, "pt"))

## 4c. Cull size model (Table S4b, Figure S1b, Methods S1)
# Bayesian shifted negative binomial regression (brms): dogs per event - 1,
# with number of events as the predictor. One value per LGA.
# Minimum ("more than X") counts treated as between X and 2X (interval censored).

# Mac only: if Stan fails to compile with the newest macOS SDK, point to an older one
# Sys.setenv(SDKROOT = "/Library/Developer/CommandLineTools/SDKs/MacOSX15.sdk")
# options(buildtools.check = function(action) TRUE)

prac_complete$Number_culled_1Event_shift <- prac_complete$Number_culled_1Event - 1
prac_complete$cens_indicator <- ifelse(is_min, "interval", "none")
prac_complete$upper_shift    <- ifelse(is_min, 2 * cull_size, cull_size) - 1

# Wide normal priors for intercept and slope; brms default for shape
cull_size_prior <- c(set_prior("normal(0,100000)", class = "Intercept"),
                     set_prior("normal(0,100000)", class = "b"))

# Main cull size model
cull_size_model <- brm(
  Number_culled_1Event_shift | cens(cens_indicator, upper_shift) ~ Frequency_12m,
  data = prac_complete, prior = cull_size_prior, family = negbinomial(),
  seed = 123, chains = 4, warmup = 1500, iter = 3500)
summary(cull_size_model)   # -> Table S4b (posterior means, 95% CrI, R-hat)
plot(cull_size_model)

# Probability that more events -> fewer dogs per event
draws <- as_draws_df(cull_size_model)
mean(draws$b_Frequency_12m < 0)

# Why minimum counts are bounded: values implied for them with no upper bound
# -> Methods 2.2.1, Methods S1 (~2-9x those reported)
prac_complete$cens_right <- ifelse(is_min, "right", "none")
cull_size_model_unbounded <- brm(
  Number_culled_1Event_shift | cens(cens_right) ~ Frequency_12m,
  data = prac_complete, prior = cull_size_prior, family = negbinomial(),
  seed = 123, chains = 4, warmup = 1500, iter = 3500, refresh = 0)
set.seed(123)
pp_unb <- posterior_predict(cull_size_model_unbounded, ndraws = 4000) + 1
sapply(which(is_min), function(j) {
  v <- pp_unb[, j]; v <- v[v > cull_size[j]]
  c(reported = cull_size[j], model_median = median(v), model_upper95 = unname(quantile(v, 0.95)))
})

# Negative binomial vs Poisson (leave-one-out cross-validation) -> Table S4b
# Uses the unbounded models: the Poisson model could not be fitted with X-2X bounds
cull_size_model_pois <- brm(
  Number_culled_1Event_shift | cens(cens_right) ~ Frequency_12m,
  data = prac_complete, prior = cull_size_prior, family = poisson(),
  seed = 123, chains = 4, warmup = 1500, iter = 3500)
loo_compare(loo(cull_size_model_unbounded, moment_match = TRUE),
            loo(cull_size_model_pois, moment_match = TRUE))

# Residual checks (DHARMa) -> Methods S1
# DHARMa does not handle censored data, so uses a version with minimum counts as exact
cull_size_model_uncens <- brm(
  Number_culled_1Event_shift ~ Frequency_12m,
  data = prac_complete, prior = cull_size_prior, family = negbinomial(),
  seed = 123, chains = 4, warmup = 1500, iter = 3500)
model_check <- createDHARMa(
  simulatedResponse       = t(posterior_predict(cull_size_model_uncens)),
  observedResponse        = prac_complete$Number_culled_1Event_shift,
  fittedPredictedResponse = apply(t(posterior_epred(cull_size_model_uncens)), 1, median),
  integerResponse = TRUE)
plotQQunif(model_check, pch = 4, col = "navy")
plotResiduals(model_check, quantreg = TRUE, pch = 4, col = "navy")
testUniformity(model_check, plot = FALSE)$p.value
testDispersion(model_check, plot = FALSE)$p.value
testOutliers(model_check, plot = FALSE)$p.value
testQuantiles(model_check, predictor = prac_complete$Frequency_12m, plot = FALSE)$p.value   # NA: too few LGAs

# Figure S1 (b): cull sizes predicted by the model for the 15 LGAs vs observed
# Observed curve: Kaplan-Meier, with minimum counts as right-censored
pp <- posterior_predict(cull_size_model, ndraws = 1000) + 1   # undo the shift
x_dogs <- 0:400
pp_cdf <- apply(pp, 1, function(v) ecdf(v)(x_dogs))           # one column per draw
pred_band <- tibble(x   = x_dogs,
                    med = apply(pp_cdf, 1, median),
                    lo  = apply(pp_cdf, 1, quantile, 0.025),
                    hi  = apply(pp_cdf, 1, quantile, 0.975))
km_lga  <- survival::survfit(survival::Surv(cull_size, !is_min) ~ 1)
obs_lga <- tibble(x = c(0, km_lga$time), F = c(0, 1 - km_lga$surv))

p_dogs <- ggplot() +
  geom_ribbon(data = pred_band, aes(x, ymin = lo, ymax = hi, fill = "Model: 95% interval")) +
  geom_step(data = pred_band, aes(x, med, colour = "Model: median")) +
  geom_step(data = obs_lga, aes(x, F, colour = "Observed"), linewidth = 1.2) +
  scale_fill_manual(values = c("Model: 95% interval" = "grey85"), name = NULL) +
  scale_colour_manual(values = c("Observed" = "grey60", "Model: median" = "black"),
                      breaks = c("Observed", "Model: median"), name = NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Dogs culled per event", y = "Cumulative proportion") +
  theme_figS1 +
  theme(legend.position = "top")

figS1 <- p_ev / p_dogs + plot_annotation(tag_levels = "a")
figS1
ggsave("FigureS1.png", figS1, width = 7, height = 8.5, dpi = 300, bg = "white")


## 4d. Simulation: main estimate and sensitivity analyses (Table S5, Results 3.1.2)
# Each scenario changes one thing from the main model. Scenarios marked (refit)
# refit the cull size model. Each runs 1,000 bootstrap resamples (~15 min in all).
# Results are saved to output/table_s5/ as they finish, so the run can be restarted.

## Helpers ----------------------------------------------------------------

# Posterior draws as intercept, slope and shape (slope = 0 if the model has none)
get_dogs_draws <- function(model) {
  d <- as.data.frame(posterior::as_draws_df(model))
  slope_col <- setdiff(grep("^b_", names(d), value = TRUE), "b_Intercept")
  data.frame(intercept = d$b_Intercept,
             slope     = if (length(slope_col)) d[[slope_col[1]]] else 0,
             shape     = d$shape)
}

# Cull size data for refitting (one row per LGA, shifted by -1)
#   idx        which survey of practice LGAs to use
#   min_mult   multiplier for minimum ("more than") counts
#   censor     FALSE treats minimum counts as exact
#   upper_fun  upper bound for minimum counts (default 2X, as main model; NULL = no bound)
main_bound <- function(s) 2 * s
make_dogs_df <- function(idx = seq_along(cull_size), min_mult = 1, censor = TRUE,
                         upper_fun = main_bound) {
  size <- ifelse(is_min[idx], round(cull_size[idx] * min_mult), cull_size[idx])
  cens <- ifelse(is_min[idx] & censor, "right", "none")
  y2   <- size
  if (!is.null(upper_fun)) {
    up   <- pmax(round(upper_fun(size)), size + 1)
    cens <- ifelse(is_min[idx] & censor, "interval", "none")
    y2   <- ifelse(cens == "interval", up, size)
  }
  data.frame(y = size - 1, y2 = y2 - 1, cens = cens,
             Frequency_12m = events[idx], log2_events = log2(events[idx]))
}

# Refit the cull size model (same settings as 4c)
fit_dogs_bayes <- function(dat, rhs = "Frequency_12m", prior = cull_size_prior) {
  if (rhs == "1") prior <- prior[prior$class == "Intercept", ]   # no slope to put a prior on
  brm(as.formula(paste("y | cens(cens, y2) ~", rhs)),
      data = dat, prior = prior, family = negbinomial(),
      seed = 123, chains = 4, warmup = 1500, iter = 3500, refresh = 0)
}

# Total dogs culled nationally in each simulated year
#   per_event = FALSE: one cull size per LGA, used for all its events (main model)
#   per_event = TRUE : a separate cull size for every event
total_dogs <- function(draws, ev_list, upper_limit = Inf, log2_events = FALSE,
                       per_event = FALSE) {
  n_lga  <- lengths(ev_list)
  n_ev   <- unlist(ev_list)
  d      <- draws[rep(sample.int(nrow(draws), length(ev_list), replace = TRUE), n_lga), ]
  x      <- if (log2_events) log2(n_ev) else n_ev
  mu     <- exp(d$intercept + d$slope * x)
  dogs_lga <- if (per_event) {
    # sum of n_ev independent shifted negative binomials
    rnbinom(length(n_ev), size = d$shape * n_ev, mu = mu * n_ev) + n_ev
  } else {
    pmin(rnbinom(length(n_ev), size = d$shape, mu = mu) + 1, upper_limit) * n_ev
  }
  sim_id <- factor(rep(seq_along(ev_list), n_lga), levels = seq_along(ev_list))
  sapply(split(dogs_lga, sim_id), sum)   # years with no culling LGAs count as 0
}

# One scenario: bootstrap culling probabilities and events, simulate national totals
run_scenario <- function(scenario,
                         draws       = get_dogs_draws(cull_size_model),
                         n_by_type   = c(N_rural, N_urban),
                         p_fun       = function(pr, pu) c(pr, pu),
                         prac_pool   = seq_along(events),
                         event_dist  = "shifted_nb",
                         ev_cap      = Inf,
                         upper_limit = Inf,
                         log2_events = FALSE,
                         per_event   = FALSE,
                         empirical   = FALSE,
                         fixed       = FALSE,
                         B = 1000, sims_per_boot = 100, seed = 123) {
  set.seed(seed)
  if (fixed) {                               # no parameter uncertainty (posterior means)
    draws <- as.data.frame(lapply(draws, mean))
    B <- 1; sims_per_boot <- B_fixed_sims
  }
  totals <- vector("list", B)
  for (b in seq_len(B)) {
    if (fixed) {
      pr <- Pcull_rural; pu <- Pcull_urban; i <- prac_pool
    } else {
      pr <- mean(sample(occ$culled_12m[occ$settlement_type == "rural"], replace = TRUE))
      pu <- mean(sample(occ$culled_12m[occ$settlement_type == "urban"], replace = TRUE))
      i  <- prac_pool[sample.int(length(prac_pool), replace = TRUE)]
    }
    if (!empirical) {
      fe <- tryCatch(suppressWarnings(fit_events_dist(events[i], event_dist)), error = function(e) NULL)
      if (is.null(fe)) next
    }
    n_culling <- Reduce(`+`, Map(function(n, p) rbinom(sims_per_boot, n, p),
                                 n_by_type, p_fun(pr, pu)))
    if (empirical) {                         # observed annual totals, no fitted distributions
      annual <- cull_size[i] * events[i]
      totals[[b]] <- sapply(n_culling, function(k) sum(annual[sample.int(length(annual), k, replace = TRUE)]))
      next
    }
    r_ev    <- r_from_fit(fe, event_dist)
    ev_list <- lapply(n_culling, function(k) if (k == 0) numeric(0) else pmin(r_ev(k), ev_cap))
    totals[[b]] <- total_dogs(draws, ev_list, upper_limit, log2_events, per_event)
  }
  total <- unlist(totals)
  tibble(scenario  = scenario,
         n_boot_ok = sum(!sapply(totals, is.null)),
         median    = median(total),
         lo        = unname(quantile(total, 0.025)),
         hi        = unname(quantile(total, 0.975)),
         pct_dogs  = 100 * median / nat_dog_pop,
         pct_lo    = 100 * lo / nat_dog_pop,
         pct_hi    = 100 * hi / nat_dog_pop)
}

# Save each result as it finishes, and skip ones already done
dir.create("output/table_s5", recursive = TRUE, showWarnings = FALSE)
cached <- function(name, expr) {
  f <- file.path("output/table_s5", paste0(name, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  res <- expr
  saveRDS(res, f)
  res
}
B_fixed_sims <- 100000

## Main estimate -> abstract, Results 3.1.2, Table S5 ------------------------
s5_main <- cached("main", run_scenario("Main model"))
s5_main

## Culling prevalence (no refit) ----------------------------------------------
n_unsampled <- c(N_rural - n_surveyed_rural, N_urban - n_surveyed_urban)
n_surveyed  <- n_surveyed_rural + n_surveyed_urban

s5_flat <- cached("flat", run_scenario(
  "Flat: one culling probability for all LGAs",
  n_by_type = N_districts,
  p_fun = function(pr, pu) (pr * n_surveyed_rural + pu * n_surveyed_urban) / n_surveyed))

unsurveyed <- function(name, f_rural, f_urban) cached(name, run_scenario(
  name,
  n_by_type = c(n_surveyed_rural, n_surveyed_urban, n_unsampled),
  p_fun = function(pr, pu) c(pr, pu, min(1, f_rural * pr), min(1, f_urban * pu))))

s5_unsamp <- bind_rows(
  unsurveyed("Unsurveyed urban LGAs 25% less likely to cull", 1,    0.75),
  unsurveyed("Unsurveyed urban LGAs 50% less likely to cull", 1,    0.5),
  unsurveyed("Unsurveyed LGAs 25% less likely to cull",       0.75, 0.75),
  unsurveyed("Unsurveyed LGAs 50% less likely to cull",       0.5,  0.5),
  unsurveyed("Unsurveyed LGAs 25% more likely to cull",       1.25, 1.25),
  unsurveyed("Unsurveyed LGAs 50% more likely to cull",       1.5,  1.5))

## Event frequency and cull size ----------------------------------------------
# LGA with most dogs culled per year removed (refit)
i_max  <- which.max(cull_size * events)
keep_i <- seq_along(events)[-i_max]
s5_dropmax <- cached("dropmax", run_scenario(
  "LGA with the most dogs culled per year removed",
  draws = get_dogs_draws(fit_dogs_bayes(make_dogs_df(keep_i))),
  prac_pool = keep_i))

# Cull size not related to event frequency (refit, no slope)
s5_indep <- cached("indep", run_scenario(
  "Cull size not related to event frequency",
  draws = get_dogs_draws(fit_dogs_bayes(make_dogs_df(), rhs = "1"))))

# Events and cull size resampled together (no fitted distributions)
s5_emp <- cached("empirical", run_scenario(
  "Events and cull size resampled together (no fitted distributions)", empirical = TRUE))

# One cull size per event rather than per LGA
s5_perevent <- cached("per_event", run_scenario(
  "One cull size per event rather than per LGA", per_event = TRUE))

# Events capped at the highest reported
s5_evcap <- cached("ev_cap", run_scenario(
  sprintf("Events capped at the highest reported (%d per year)", max(events)),
  ev_cap = max(events)))

# Events log2-transformed in cull size model (refit)
s5_log2 <- cached("log2", run_scenario(
  "Events log2-transformed in cull size model",
  draws = get_dogs_draws(fit_dogs_bayes(make_dogs_df(), rhs = "log2_events")),
  log2_events = TRUE))

# Zero-truncated negative binomial for events (main model uses shifted negative binomial)
s5_trnb <- cached("trunc_nb", run_scenario(
  "Zero-truncated negative binomial for events", event_dist = "trunc_nbinom"))

## Minimum ("more than") counts (refit) ---------------------------------------
s5_exact <- cached("min_exact", run_scenario(
  "Minimum counts treated as exact (lower scenario)",
  draws = get_dogs_draws(fit_dogs_bayes(make_dogs_df(censor = FALSE)))))

s5_x15 <- cached("min_x15", run_scenario(
  "Minimum counts between X and 1.5X",
  draws = get_dogs_draws(fit_dogs_bayes(make_dogs_df(upper_fun = function(s) 1.5 * s)))))

s5_x3 <- cached("min_x3", run_scenario(
  "Minimum counts between X and 3X",
  draws = get_dogs_draws(fit_dogs_bayes(make_dogs_df(upper_fun = function(s) 3 * s)))))

# Same upper bound of 300 for every minimum count
s5_x300 <- cached("min_x300", run_scenario(
  "Minimum counts between X and 300",
  draws = get_dogs_draws(fit_dogs_bayes(make_dogs_df(upper_fun = function(s) 300 + 0 * s)))))

s5_unbounded <- cached("min_unbounded", run_scenario(
  "Minimum counts with no upper bound (upper scenario)",
  draws = get_dogs_draws(fit_dogs_bayes(make_dogs_df(upper_fun = NULL)))))

## Dogs per event additionally capped (minimum counts between X and the smaller
## of 2X and the cap; simulated cull sizes capped too) --------------------------
s5_upper <- bind_rows(lapply(c(250, 300, 400), function(ub) cached(paste0("upper", ub), run_scenario(
  paste("Dogs per event capped at", ub),
  draws = get_dogs_draws(fit_dogs_bayes(make_dogs_df(upper_fun = function(s) pmin(2 * s, ub)))),
  upper_limit = ub))))

## Priors (refit) -------------------------------------------------------------
weak_prior <- c(set_prior("normal(0, 5)", class = "Intercept"),
                set_prior("normal(0, 1)", class = "b"))
s5_prior <- cached("weak_prior", run_scenario(
  "Weakly informative priors",
  draws = get_dogs_draws(fit_dogs_bayes(make_dogs_df(), prior = weak_prior))))

## Parameter uncertainty --------------------------------------------------------
s5_fixed <- cached("fixed", run_scenario(
  "Fixed parameters (no parameter uncertainty)", fixed = TRUE))

## Seed check for the main model (Table S5 footnote) ---------------------------
s5_seeds <- bind_rows(lapply(c(1, 42, 999, 2024), function(s) cached(paste0("seed", s),
                                                                     run_scenario(paste("Seed", s), seed = s))))
range(c(s5_main$median, s5_seeds$median))

## Table S5 -------------------------------------------------------------------
table_s5 <- bind_rows(s5_main, s5_flat, s5_unsamp, s5_dropmax, s5_indep, s5_emp,
                      s5_perevent, s5_evcap, s5_log2, s5_trnb, s5_exact, s5_x15,
                      s5_x3, s5_x300, s5_unbounded, s5_upper, s5_prior, s5_fixed) %>%
  mutate(estimate = sprintf("%s (%s-%s)",
                            format(round(median), big.mark = ","),
                            format(round(lo), big.mark = ","),
                            format(round(hi), big.mark = ",")),
         pct = sprintf("%.2f (%.2f-%.2f)", pct_dogs, pct_lo, pct_hi))
print(dplyr::select(table_s5, scenario, n_boot_ok, estimate, pct), n = Inf, width = Inf)
write.csv(dplyr::select(table_s5, scenario, n_boot_ok, estimate, pct),
          "output/table_s5/table_s5.csv", row.names = FALSE)


#### 5. Surveyed vs not surveyed LGAs (Table S3) ####

comp <- councils_2022 %>%
  mutate(group    = ifelse(surveyed == "yes", "Surveyed", "Not surveyed"),
         hdr      = human_pop_2022 / dog_pop_2022,
         zanzibar = zone_tdhs2022 == "Zanzibar")
table(comp$group, comp$urban_rural)

## 5a. Functions
# Standardised mean difference (surveyed minus not surveyed)
smd_cont <- function(x, g) {
  a <- x[g == "Surveyed"]; b <- x[g == "Not surveyed"]
  (mean(a) - mean(b)) / sqrt((var(a) + var(b)) / 2)
}
smd_bin <- function(x, g) {
  p1 <- mean(x[g == "Surveyed"]); p2 <- mean(x[g == "Not surveyed"])
  den <- sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2)
  if (den == 0) NA_real_ else (p1 - p2) / den
}

# Median (IQR); populations rounded to avoid identifying individual LGAs
med_iqr <- function(x, dp = 1, round_to = NA) {
  q <- quantile(x, c(0.5, 0.25, 0.75))
  f <- if (is.na(round_to)) function(v) formatC(v, format = "f", digits = dp)
  else function(v) formatC(round(v / round_to) * round_to, format = "d", big.mark = ",")
  sprintf("%s (%s-%s)", f(q[1]), f(q[2]), f(q[3]))
}

cont_vars <- tribble(
  ~var,                       ~label,                              ~dp, ~round_to, ~log_smd,
  "human_pop_2022",           "Human population, median (IQR)",    0,   1000,      TRUE,
  "dog_pop_2022",             "Dog population, median (IQR)",      0,   100,       TRUE,
  "hdr",                      "Human:dog ratio, median (IQR)",     1,   NA,        FALSE,
  "mean_household_size_2022", "Mean household size, median (IQR)", 1,   NA,        FALSE
)

make_table <- function(dat, stratum) {
  g <- dat$group
  out <- list()
  
  # Population SMDs on the log scale
  for (i in seq_len(nrow(cont_vars))) {
    x  <- dat[[cont_vars$var[i]]]
    xs <- if (cont_vars$log_smd[i]) log(x) else x
    out[[length(out) + 1]] <- tibble(
      Variable       = cont_vars$label[i],
      Surveyed       = med_iqr(x[g == "Surveyed"], cont_vars$dp[i], cont_vars$round_to[i]),
      `Not surveyed` = med_iqr(x[g == "Not surveyed"], cont_vars$dp[i], cont_vars$round_to[i]),
      SMD            = smd_cont(xs, g))
  }
  
  yes_no <- function(label, x) tibble(
    Variable       = label,
    Surveyed       = n_pct(x[g == "Surveyed"]),
    `Not surveyed` = n_pct(x[g == "Not surveyed"]),
    SMD            = smd_bin(x, g))
  
  # Settlement type and zone in the all-LGA table only
  if (n_distinct(dat$urban_rural) > 1) {
    out[[length(out) + 1]] <- yes_no("Urban council, n (%)", dat$urban_rural == "urban")
    for (z in c("North", "East", "South", "West", "Central"))
      out[[length(out) + 1]] <- yes_no(paste0("Zone: ", z, ", n (%)"), dat$zone_grouped == z)
  }
  
  n_reg <- n_distinct(dat$region)
  out[[length(out) + 1]] <- tibble(
    Variable       = "Regions represented, n/N",
    Surveyed       = sprintf("%d/%d", n_distinct(dat$region[g == "Surveyed"]), n_reg),
    `Not surveyed` = sprintf("%d/%d", n_distinct(dat$region[g == "Not surveyed"]), n_reg),
    SMD            = NA_real_)
  
  bind_rows(out) %>%
    mutate(Stratum = sprintf("%s (surveyed n = %d; not surveyed n = %d)", stratum,
                             sum(g == "Surveyed"), sum(g == "Not surveyed")),
           SMD = round(SMD, 2)) %>%
    relocate(Stratum)
}

## 5b. Table
# -> Table S3, Results 3.1.1 (21/31 regions, SMDs)
supp_table <- bind_rows(
  make_table(comp, "All LGAs"),
  make_table(filter(comp, urban_rural == "rural"), "Rural LGAs"),
  make_table(filter(comp, urban_rural == "urban"), "Urban LGAs")
)
print(supp_table, n = Inf)

## 5c. Logistic regression: characteristics associated with being surveyed
# -> Table S3 OR and p columns
# Outcome: surveyed (1) vs not surveyed (0); one model per characteristic
comp <- comp %>%
  mutate(surveyed_bin = as.integer(group == "Surveyed"),
         log2_human   = log2(human_pop_2022),
         log2_dog     = log2(dog_pop_2022),
         urban        = urban_rural == "urban",
         zone_grouped = relevel(factor(zone_grouped), ref = "North"))

s3_predictors <- c(log2_human               = "Human population (per doubling)",
                   log2_dog                 = "Dog population (per doubling)",
                   hdr                      = "Human:dog ratio (per 1-unit increase)",
                   mean_household_size_2022 = "Mean household size (per person)",
                   urban                    = "Urban LGA",
                   zone_grouped             = "Zone (ref: North)")

fit_surveyed <- function(v, dat) {
  m  <- glm(reformulate(v, response = "surveyed_bin"), data = dat, family = binomial)
  m0 <- glm(surveyed_bin ~ 1, data = dat, family = binomial)
  ci <- suppressMessages(exp(confint(m)))[-1, , drop = FALSE]
  tibble(variable = s3_predictors[[v]],
         term     = rownames(ci),
         OR_ci    = sprintf("%.2f (%.2f-%.2f)", exp(coef(m))[-1], ci[, 1], ci[, 2]),
         p_LRT    = signif(anova(m0, m, test = "LRT")$`Pr(>Chi)`[2], 2))
}

s3_glm <- bind_rows(
  bind_rows(lapply(names(s3_predictors), fit_surveyed, dat = comp)) %>%
    mutate(stratum = "All LGAs"),
  bind_rows(lapply(names(s3_predictors)[1:4], fit_surveyed,
                   dat = filter(comp, urban_rural == "rural"))) %>%
    mutate(stratum = "Rural LGAs"),
  bind_rows(lapply(names(s3_predictors)[1:4], fit_surveyed,
                   dat = filter(comp, urban_rural == "urban"))) %>%
    mutate(stratum = "Urban LGAs")
) %>% relocate(stratum)
print(s3_glm, n = Inf)


#### 6. Session info ####
sessionInfo()