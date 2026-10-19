# =============================================================================
# EXPERIMENT 3: Myrmeleon Field Experiment 2025
# Environmental Exploratory Analysis + Bayesian Models per Species
# =============================================================================
# Experiment: day 1 = 22 April 2025, day 71 = 1 July 2025  (n_days = 71)
# Species: Myrmeleon_almohadarum, Myrmeleon_hyalinus
#
# SECTION A — Exploratory environmental analysis
#   - Time series of max daily temp and mean daily RH
#   - Distribution plots (density, boxplot by month/week)
#   - Rolling-mean smooths (7-day, 14-day)
#   - Temp × RH scatter
#   - Summary stats table
#
# SECTION B — Bayesian models
# Structured to follow the same section layout, diagnostics (Rhat +
# trace plots), and summary-table outputs as the experiment 1 and 2
# scripts: B1 Load Data, B2 Prepare Datasets, B3 Priors & MCMC,
# B4 Fit Models, B5 Convergence Diagnostics, B6 Posterior Summaries,
# B7 Posterior Significance Check, B8 Summary Tables, B9 Figures.
#
# Environmental covariates:
#   Hazard model (1):  1-day lag (previous day's max temp). Only temperature
#     is included in the fitted formula below due to correlation with RH.
#   Duration (2-4): species fixed effect only
#   Morphometrics (5a-c): species fixed effect only
#
# Models:
#  1. Larval mortality  — Bernoulli cloglog:
#       s(day_sc, k=4) + s(temp_1d_sc, by=sp, k=3) + sp + (1|ID)
#  2. Days to pupation  — lognormal: sp
#  3. Days pupa to imago  — lognormal: sp
#  4. Total days        — lognormal: sp
#  5a-c. Morphometrics  — Gaussian (raw scale): sp
#
# Required files in the working directory:
#   Data_2025_days_to.csv  (individual summary — event times, morphometrics)
#   data_hobo_2025.csv     (HOBO logger — 15-min temp/RH readings, Section A)
# =============================================================================

# -- 0. Packages ---------------------------------------------------------------
library(brms)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(bayesplot)
library(scales)    # for date_breaks in time series


setwd("LINK TO REPOSITORY")

# Clear stale model caches
file.remove(list.files(getwd(), pattern = "^model_.*exp3\\.rds$", full.names = TRUE))

# set.seed(42)

# =============================================================================
# SECTION A — Exploratory Environmental Analysis
# =============================================================================

# -- A1. Load and parse HOBO data ---------------------------------------------

hobo <- read.csv("data_hobo_2025.csv", fileEncoding = "latin1",
                 check.names = FALSE, stringsAsFactors = FALSE)
hobo$dt   <- as.POSIXct(hobo$Date_Time, format = "%m.%d.%y %I:%M:%S %p")
hobo$date <- as.Date(hobo$dt)

# Daily summaries
daily_env <- hobo %>%
  group_by(date) %>%
  summarise(
    max_temp  = max(Temp_C,  na.rm = TRUE),
    min_temp  = min(Temp_C,  na.rm = TRUE),
    mean_temp = mean(Temp_C, na.rm = TRUE),
    mean_rh   = mean(RH,     na.rm = TRUE),
    max_rh    = max(RH,      na.rm = TRUE),
    min_rh    = min(RH,      na.rm = TRUE),
    n_obs     = sum(!is.na(Temp_C)),
    .groups   = "drop"
  ) %>%
  mutate(
    exp_day   = as.integer(date - as.Date("2025-04-22")) + 1L,
    week      = ceiling(pmax(exp_day, 1) / 7),
    month_lab = format(date, "%B")
  )

# Study-period subset (exp days 1–71)
daily_study <- daily_env %>% filter(exp_day >= 1, exp_day <= 71)

# Rolling means (7-day and 14-day) using zoo::rollmean
# Computed inline with slider::slide_dbl (base-compatible alternative)
roll_mean <- function(x, k) {
  out <- rep(NA_real_, length(x))
  for (i in seq_along(x)) {
    idx <- max(1, i - k + 1):i
    out[i] <- mean(x[idx], na.rm = TRUE)
  }
  out
}

daily_study <- daily_study %>%
  arrange(exp_day) %>%
  mutate(
    roll7_temp  = roll_mean(max_temp, 7),
    roll14_temp = roll_mean(max_temp, 14),
    roll7_rh    = roll_mean(mean_rh,  7),
    roll14_rh   = roll_mean(mean_rh,  14)
  )

cat("=== Environmental Summary (Study Period, Days 1-71) ===\n")
env_summary <- daily_study %>%
  summarise(across(c(max_temp, min_temp, mean_temp, mean_rh, max_rh, min_rh),
                   list(mean = ~round(mean(., na.rm=TRUE), 2),
                        sd   = ~round(sd(., na.rm=TRUE),   2),
                        min  = ~round(min(., na.rm=TRUE),   2),
                        max  = ~round(max(., na.rm=TRUE),   2))))
print(as.data.frame(env_summary))

# -- A2. Exploratory Figures ---------------------------------------------------

pal_env <- c(max_temp = "#D6604D", mean_rh = "#2166AC")

# Figure E1: Time series of daily max temp and mean RH
p_ts_temp <- ggplot(daily_study, aes(x = exp_day)) +
  geom_line(aes(y = max_temp), colour = pal_env["max_temp"], linewidth = 0.5,
            alpha = 0.6) +
  geom_line(aes(y = roll7_temp),  colour = pal_env["max_temp"],
            linewidth = 1.0, linetype = "solid") +
  geom_line(aes(y = roll14_temp), colour = pal_env["max_temp"],
            linewidth = 1.2, linetype = "dashed") +
  labs(x = "Experiment day", y = "Max daily temperature (°C)"
       # subtitle = "Thin: daily; solid: 7-day rolling mean; dashed: 14-day rolling mean"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

p_ts_rh <- ggplot(daily_study, aes(x = exp_day)) +
  geom_line(aes(y = mean_rh), colour = pal_env["mean_rh"], linewidth = 0.5,
            alpha = 0.6) +
  geom_line(aes(y = roll7_rh),  colour = pal_env["mean_rh"],
            linewidth = 1.0, linetype = "solid") +
  geom_line(aes(y = roll14_rh), colour = pal_env["mean_rh"],
            linewidth = 1.2, linetype = "dashed") +
  labs(x = "Experiment day", y = "Mean daily relative humidity (%)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

figE1 <- p_ts_temp / p_ts_rh +
  plot_annotation(title   = "Environmental Conditions during Field Experiment 2025",
                  theme   = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave("plots/FigureE1_Env_TimeSeries.png", figE1, width = 10, height = 7, dpi = 300)


# Figure E2: Weekly boxplots of temp and RH
p_box_temp <- ggplot(daily_study, aes(x = factor(week), y = max_temp)) +
  geom_boxplot(fill = pal_env["max_temp"], alpha = 0.5, outlier.size = 1.5) +
  labs(x = "Experiment week", y = "Max daily temperature (°C)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

p_box_rh <- ggplot(daily_study, aes(x = factor(week), y = mean_rh)) +
  geom_boxplot(fill = pal_env["mean_rh"], alpha = 0.5, outlier.size = 1.5) +
  labs(x = "Experiment week", y = "Mean daily RH (%)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

figE2 <- p_box_temp / p_box_rh 

ggsave("plots/FigureE2_Env_WeeklyBoxplots.png", figE2, width = 10, height = 7, dpi = 300)


# Figure E3: Temp × RH scatter coloured by experiment day
figE3 <- ggplot(daily_study, aes(x = max_temp, y = mean_rh, colour = exp_day)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "loess", colour = "grey30", linewidth = 0.8, se = FALSE) +
  scale_colour_viridis_c(option = "plasma", name = "Experiment Day") +
  labs(x = "Max daily temperature (°C)", y = "Mean daily RH (%)"
       # caption = "Loess smooth shown in grey."
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title   = element_text(face = "bold", size = 14),
        legend.position = "right")

ggsave("plots/FigureE3_Env_TempRH_Scatter.png", figE3, width = 8, height = 5, dpi = 300)

# Figure E4: Temperature range (min-mean-max) ribbon
figE4 <- ggplot(daily_study, aes(x = exp_day)) +
  geom_ribbon(aes(ymin = min_temp, ymax = max_temp),
              fill = pal_env["max_temp"], alpha = 0.2) +
  geom_line(aes(y = mean_temp), colour = pal_env["max_temp"], linewidth = 0.9) +
  labs(x = "Experiment day", y = "Temperature (°C)"
       # caption = "Band: min–max range. Line: daily mean."
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14))
ggsave("plots/FigureE4_Env_TempRange.png", figE4, width = 10, height = 5, dpi = 300)


# =============================================================================
# SECTION B — Bayesian Models
# =============================================================================

# ── B1. Load Data ────────────────────────────────────────────────────────────

n_days <- 71 # length of experiment (22 Apr - 1 Jul 2025)

# --- B1a. 1-day lagged environmental covariates (for the mortality hazard) ---
# For each experiment day d: previous day's max temp (day d-1). Built from
# daily_env, which Section A already computed from the HOBO logger data.
env_lookup <- daily_env %>%
  select(exp_day, max_temp, mean_rh)

env_1d <- tibble(exp_day = 1:n_days) %>%
  mutate(
    temp_1d = map_dbl(exp_day, function(d) {
      v <- env_lookup$max_temp[env_lookup$exp_day == (d - 1)]
      if (length(v) == 0 || all(is.na(v))) NA_real_ else v[1]
    }),
    rh_1d = map_dbl(exp_day, function(d) {
      v <- env_lookup$mean_rh[env_lookup$exp_day == (d - 1)]
      if (length(v) == 0 || all(is.na(v))) NA_real_ else v[1]
    })
  )

# --- B1b. Biological summary data ---------------------------------------------

required_cols <- c("Numb_days_before_death", "Number_days_P", "Number_days_I")

dat_summ <- read.csv("Data_2025_days_to.csv", fileEncoding = "latin1",
                     check.names = FALSE)
names(dat_summ) <- gsub("\n", " ", names(dat_summ))
dat_summ <- dat_summ[, !is.na(names(dat_summ)) & names(dat_summ) != ""]

dat_summ <- dat_summ %>%
  mutate(sp = factor(sp))

sp_levels <- sort(levels(dat_summ$sp))
cat(sprintf("Individuals: %d\n", nrow(dat_summ)))
cat("Species:\n"); print(table(dat_summ$sp))

# ── B2. Prepare Analysis Datasets ─────────────────────────────────────────────

# --- B2a. Larval mortality discrete-time data ---------------------------------
# Build one row per individual-day-at-risk, with the previous day's max temp
# joined in as the hazard model's environmental covariate.
# At risk until death (event=1) or pupation (right-censored at pupa day - 1).
# Individuals that neither died nor pupated are censored at n_days.

mort_events <- dat_summ %>%
  mutate(
    day_death = as.integer(Numb_days_before_death),
    end_day   = case_when(
      !is.na(as.integer(Numb_days_before_death)) ~ as.integer(Numb_days_before_death),
      !is.na(as.integer(Number_days_P))           ~ as.integer(Number_days_P) - 1L,
      TRUE                                         ~ n_days
    )
  )

dt_mortality <- mort_events %>%
  filter(end_day >= 1) %>%
  rowwise() %>%
  reframe(
    ID        = ID,
    sp        = sp,
    day       = seq_len(end_day),
    event     = as.integer(seq_len(end_day) == day_death & !is.na(day_death))
  ) %>%
  mutate(day_scaled = day / n_days) %>%
  left_join(env_1d, by = c("day" = "exp_day")) %>%
  # Scale env covariates; impute rare NAs at lag day 0 with the sample mean
  mutate(
    temp_1d    = if_else(is.na(temp_1d), mean(temp_1d, na.rm = TRUE), temp_1d),
    rh_1d      = if_else(is.na(rh_1d),   mean(rh_1d,   na.rm = TRUE), rh_1d),
    temp_1d_sc = scale(temp_1d)[, 1],
    rh_1d_sc   = scale(rh_1d)[, 1]
  )

cat(sprintf("Mortality: %d rows, %d events\n",
            nrow(dt_mortality), sum(dt_mortality$event)))

# --- B2b. Lognormal duration datasets (species-only, no env covariates) ------
df_pup <- dat_summ %>%
  filter(!is.na(Number_days_P)) %>%
  mutate(days_P = as.integer(Number_days_P))

df_eme <- dat_summ %>%
  filter(!is.na(Number_days_I)) %>%
  mutate(days_I = as.integer(Number_days_I))

df_total <- dat_summ %>%
  filter(!is.na(Number_days_P) & !is.na(Number_days_I)) %>%
  mutate(days_total = as.integer(Number_days_P) + as.integer(Number_days_I))

cat(sprintf("Pupation n=%d | Emergence n=%d | Total n=%d\n",
            nrow(df_pup), nrow(df_eme), nrow(df_total)))

# --- B2c. Morphometric datasets -----------------------------------------------
# Rename columns once for tidiness
dat_summ <- dat_summ %>%
  rename_with(~ "body_length", matches("Body_length_I"))  %>%
  rename_with(~ "wing_length", matches("Length_front_")) %>%
  rename_with(~ "weight_I",    matches("Weight_I"))

df_morph <- dat_summ %>%
  filter(!is.na(body_length) | !is.na(wing_length) | !is.na(weight_I))
cat(sprintf("Morphometrics n=%d\n", nrow(df_morph)))

# ── B3. Priors & MCMC Settings ────────────────────────────────────────────────

n_chains <- 4
n_iter   <- 6000
n_warmup <- 2000

# Mortality hazard model
hazard_priors <- c(
  prior(normal(0, 0.5),  class = b),
  prior(normal(0, 0.5),  class = Intercept),
  prior(normal(0, 1),    class = sds, lb = 0),
  prior(normal(0, 0.25), class = sd,  lb = 0)
)

# Lognormal priors for the duration models (2, 3, 4).
lognorm_priors <- c(
  prior(normal(0, 0.5), class = b),
  prior(normal(0, 0.5), class = Intercept),
  prior(normal(0, 1),   class = sigma, lb = 0)
)


# Cache invalidation helper: remove .rds if it holds no posterior draws
clear_bad_cache <- function(path) {
  if (file.exists(path)) {
    m <- readRDS(path)
    if (is.null(tryCatch(as.array(m), error = function(e) NULL))) {
      message("Removing failed cache: ", path); file.remove(path)
    }
  }
}

# ── B4. Fit Models ────────────────────────────────────────────────────────────

# --- Model 1: Larval mortality (discrete-time hazard) ------------------------
# - s(day_scaled, k=3): temporal baseline hazard (low k: only 12 death events
#   in the whole dataset, so a very flexible smooth isn't identifiable)
# - s(temp_1d_sc, by=sp, k=3): species-specific smooth effect of the previous
#   day's max temperature
# - sp: species main-effect intercept difference
# - (1|ID): individual random intercept

cat("\n--- Model 1: Larval mortality ---\n")

clear_bad_cache("model_Mortality_exp3.rds")

mod_mortality <- brm(
  formula = event ~ s(day_scaled, bs = "cr", k = 3) +
    s(temp_1d_sc, by = sp, bs = "cr", k = 3) +
    sp + (1 | ID),
  data    = dt_mortality,
  family  = bernoulli(link = "cloglog"),
  prior   = hazard_priors,
  chains  = n_chains, iter = n_iter, warmup = n_warmup, cores = 4,
  seed    = 42, init = 0,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  file    = "model_Mortality_exp3"
)

# --- Model 2: Days to pupation (lognormal) -----------------------------------
cat("\n--- Model 2: Days to pupation (lognormal) ---\n")

clear_bad_cache("model_Days_Pupation_exp3.rds")

mod_pup <- brm(
  formula = days_P ~ sp,
  data    = df_pup,
  family  = lognormal(),
  prior   = lognorm_priors,
  chains  = n_chains, iter = n_iter, warmup = n_warmup, cores = 4,
  seed    = 42, init = 0,
  control = list(adapt_delta = 0.95),
  file    = "model_Days_Pupation_exp3"
)

# --- Model 3: Days pupa to imago (lognormal) ---------------------------------
cat("\n--- Model 3: Days pupa to imago (lognormal) ---\n")

clear_bad_cache("model_Days_Emergence_exp3.rds")

mod_eme <- brm(
  formula = days_I ~ sp,
  data    = df_eme,
  family  = lognormal(),
  prior   = lognorm_priors,
  chains  = n_chains, iter = n_iter, warmup = n_warmup, cores = 4,
  seed    = 42, init = 0,
  control = list(adapt_delta = 0.95),
  file    = "model_Days_Emergence_exp3"
)

# --- Model 4: Total days larva to imago (lognormal) -------------------------
cat("\n--- Model 4: Total days to imago (lognormal) ---\n")

clear_bad_cache("model_Days_Total_exp3.rds")

mod_total <- brm(
  formula = days_total ~ sp,
  data    = df_total,
  family  = lognormal(),
  prior   = lognorm_priors,
  chains  = n_chains, iter = n_iter, warmup = n_warmup, cores = 4,
  seed    = 42, init = 0,
  control = list(adapt_delta = 0.95),
  file    = "model_Days_Total_exp3"
)

# --- Models 5a-c: Morphometrics (Gaussian, raw scale) ------------------------
fit_morph <- function(var, label, data) {
  cat(sprintf("\n--- Model 5: %s ---\n", label))
  
  clear_bad_cache(paste0("model_", label, "_exp3.rds"))
  
  d <- data %>% filter(!is.na(.data[[var]]))
  brm(
    formula = reformulate("sp", response = var),
    data    = d,
    family  = gaussian(),
    chains  = n_chains, iter = n_iter, warmup = n_warmup, cores = 4,
    seed    = 42, init = 0,
    control = list(adapt_delta = 0.95),
    file    = paste0("model_", label, "_exp3"),
    file_refit = "on_change"  # refit automatically since the formula/prior
    # changed from the previous z-scored version
  )
}

mod_bl <- fit_morph("body_length", "Body_length",  df_morph)
mod_wl <- fit_morph("wing_length", "Wing_length",  df_morph)
mod_wt <- fit_morph("weight_I",    "Weight_imago", df_morph)

# ── B5. Convergence Diagnostics ───────────────────────────────────────────────

# --- B5a. Rhat summary ---------------------------------------------------------
get_rhats <- function(model, label) {
  s    <- summary(model)
  vals <- c(
    s$fixed[, "Rhat"],
    if (!is.null(s$random))  s$random[[1]][, "Rhat"],
    if (!is.null(s$splines)) s$splines[, "Rhat"]
  )
  vals <- as.numeric(vals)
  cat(sprintf("[%s]  Rhat range: %.4f – %.4f  (all ≤ 1.01: %s)\n",
              label,
              min(vals, na.rm = TRUE), max(vals, na.rm = TRUE),
              ifelse(all(vals <= 1.01, na.rm = TRUE), "YES", "NO")))
  data.frame(model = label, rhat = vals)
}

cat("\n=== Convergence (Rhat) ===\n")
rhat_df <- bind_rows(
  get_rhats(mod_mortality, "1. Mortality"),
  get_rhats(mod_pup,       "2. Days to Pupation (lognormal)"),
  get_rhats(mod_eme,       "3. Days to Emergence (lognormal)"),
  get_rhats(mod_total,     "4. Total Days (lognormal)"),
  get_rhats(mod_bl,        "5a. Body Length"),
  get_rhats(mod_wl,        "5b. Wing Length"),
  get_rhats(mod_wt,        "5c. Weight")
)

# --- B5b. Trace plots ----------------------------------------------------------

make_trace <- function(model, title) {
  mcmc_plot(model,
            variable = c("^b_", "^sd_", "^sds_", "^sigma$"),
            regex    = TRUE,
            type     = "trace") +
    labs(title = title) +
    theme_classic(base_size = 8) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 10),
          axis.text  = element_text(size = 6))
}

cat("\n=== Trace plots ===\n")
trace_list <- list(
  make_trace(mod_mortality, "1. Mortality"),
  make_trace(mod_pup,       "2. Days to Pupation"),
  make_trace(mod_eme,       "3. Days to Emergence"),
  make_trace(mod_total,     "4. Total Days"),
  make_trace(mod_bl,        "5a. Body Length"),
  make_trace(mod_wl,        "5b. Wing Length"),
  make_trace(mod_wt,        "5c. Weight")
)

# ── B6. Posterior Summaries ───────────────────────────────────────────────────

cat("\n=== Posterior Summaries ===\n")
models_list <- list(
  list(mod_mortality, "1. Larval Mortality (hazard)"),
  list(mod_pup,       "2. Days to Pupation (lognormal)"),
  list(mod_eme,       "3. Days Pupa to Imago (lognormal)"),
  list(mod_total,     "4. Total Days to Imago (lognormal)"),
  list(mod_bl,        "5a. Body Length (Gaussian)"),
  list(mod_wl,        "5b. Wing Length (Gaussian)"),
  list(mod_wt,        "5c. Weight (Gaussian)")
)

for (obj in models_list) {
  cat(sprintf("\n========== %s ==========\n", obj[[2]]))
  print(summary(obj[[1]]))
}

# ── B7. Posterior Significance Check ──────────────────────────────────────────

cat("\n=== Fixed-Effect Posteriors (species) ===\n")
report_posteriors <- function(model, label) {
  dr     <- as_draws_df(model)
  b_cols <- grep("^b_", names(dr), value = TRUE)
  b_cols <- b_cols[!grepl("Intercept", b_cols)]
  if (!length(b_cols)) { cat(sprintf("[%s] no fixed effects\n", label)); return() }
  cat(sprintf("\n[%s]\n", label))
  for (col in b_cols) {
    x  <- dr[[col]]
    ci <- quantile(x, c(0.025, 0.975))
    cat(sprintf("  %-45s  median=%7.3f  95%%CI [%7.3f, %7.3f]  P(>0)=%.3f\n",
                col, median(x), ci[1], ci[2], mean(x > 0)))
  }
}

for (obj in models_list)
  report_posteriors(obj[[1]], obj[[2]])

# ── B8. Summary Tables ─────────────────────────────────────────────────────────

# --- B8a. Full parameter table --------------------------------------------------
# Identical helper to experiments 1 and 2 — it just reads whatever blocks
# summary() returns, so it needs no changes for this script's formulas.
extract_summary_table <- function(model, label) {
  s <- summary(model)
  
  pull_block <- function(block, component) {
    if (is.null(block) || nrow(block) == 0) return(NULL)
    df <- as.data.frame(block)
    df$term      <- rownames(block)
    df$component <- component
    df
  }
  
  blocks <- list(
    pull_block(s$fixed,       "population-level (b)"),
    pull_block(s$random[[1]], "group-level SD (ID)"),
    pull_block(s$splines,     "smooth SD (sds)"),
    pull_block(s$spec_pars,   "family-specific")
  )
  
  out <- bind_rows(blocks)
  out$model <- label
  out %>%
    select(model, component, term,
           Estimate, Est.Error, `l-95% CI`, `u-95% CI`, Rhat, Bulk_ESS, Tail_ESS)
}

summary_table <- bind_rows(
  extract_summary_table(mod_mortality, "1. Mortality"),
  extract_summary_table(mod_pup,       "2. Days to Pupation"),
  extract_summary_table(mod_eme,       "3. Days to Emergence"),
  extract_summary_table(mod_total,     "4. Total Days"),
  extract_summary_table(mod_bl,        "5a. Body Length"),
  extract_summary_table(mod_wl,        "5b. Wing Length"),
  extract_summary_table(mod_wt,        "5c. Weight")
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

write_csv(summary_table, "Table_ModelSummaries_exp3.csv")

# --- B8b. Fixed-effect table -----------------------------------------------------
get_fixed_effect_table <- function(model, label) {
  dr     <- as_draws_df(model)
  b_cols <- grep("^b_", names(dr), value = TRUE)
  b_cols <- b_cols[!grepl("Intercept", b_cols)]
  if (!length(b_cols)) return(NULL)
  
  map_dfr(b_cols, function(col) {
    x  <- dr[[col]]
    ci <- quantile(x, c(0.025, 0.975))
    data.frame(
      model      = label,
      term       = col,
      median     = median(x),
      lower_95   = ci[1],
      upper_95   = ci[2],
      p_gt_0     = mean(x > 0)
    )
  })
}

fixed_effect_table <- bind_rows(
  get_fixed_effect_table(mod_mortality, "1. Mortality"),
  get_fixed_effect_table(mod_pup,       "2. Days to Pupation"),
  get_fixed_effect_table(mod_eme,       "3. Days to Emergence"),
  get_fixed_effect_table(mod_total,     "4. Total Days"),
  get_fixed_effect_table(mod_bl,        "5a. Body Length"),
  get_fixed_effect_table(mod_wl,        "5b. Wing Length"),
  get_fixed_effect_table(mod_wt,        "5c. Weight")
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

write_csv(fixed_effect_table, "Table_FixedEffects_exp3.csv")

# ── B9. Figures ────────────────────────────────────────────────────────────────

pal   <- c("Myrmeleon_almohadarum" = "black",
           "Myrmeleon_hyalinus"    = "grey")
xlabs <- c("Myrmeleon_almohadarum" = "M. almohadarum",
           "Myrmeleon_hyalinus"    = "M. hyalinus")

# Helper: posterior species mean ± 95% CI + jitter
plot_sp <- function(model, data, resp_var, ylab, title,
                    back_transform = NULL) {
  nd  <- data.frame(sp = factor(sp_levels, levels = sp_levels))
  pf  <- fitted(model, newdata = nd, summary = TRUE) %>%
    as.data.frame() %>% bind_cols(nd) %>%
    rename(est = Estimate, lo = Q2.5, hi = Q97.5)
  
  if (!is.null(back_transform)) {
    pf  <- pf  %>% mutate(across(c(est, lo, hi), back_transform))
    raw <- data %>% filter(!is.na(.data[[resp_var]]))
  } else {
    raw <- data
  }
  
  ggplot() +
    geom_jitter(data = raw,
                aes(x = sp, y = .data[[resp_var]], colour = sp),
                width = 0.12, alpha = 0.45, size = 1.8) +
    geom_pointrange(data = pf,
                    aes(x = sp, y = est, ymin = lo, ymax = hi, colour = sp),
                    size = 1.0, linewidth = 1.2, shape = 18) +
    scale_colour_manual(values = pal, guide = "none") +
    scale_x_discrete(labels = xlabs) +
    labs(x = NULL, y = ylab, title = title) +
    theme_classic(base_size = 12) +
    theme(axis.text.x  = element_text(face = "italic"),
          plot.title   = element_text(face = "bold", size = 12))
}

# ── Figure 1: Mortality — hazard temporal curve + species x temperature effect
# Predict temporal hazard at mean temperature (sc=0), separately per species
grid_sc <- seq(min(dt_mortality$day_scaled), max(dt_mortality$day_scaled),
               length.out = 150)
nd_mort_sp <- expand.grid(
  day_scaled = grid_sc,
  sp         = factor(sp_levels, levels = sp_levels)
) %>% mutate(temp_1d_sc = 0, ID = NA_character_)

pred_haz <- fitted(mod_mortality, newdata = nd_mort_sp,
                   allow_new_levels = TRUE, re_formula = NA) %>%
  as.data.frame() %>% bind_cols(nd_mort_sp) %>%
  rename(estimate = Estimate, lower = Q2.5, upper = Q97.5) %>%
  mutate(day_orig = day_scaled * n_days, sp = as.character(sp))

p_haz <- ggplot(pred_haz, aes(x = day_orig, y = estimate,
                              colour = sp, fill = sp)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = pal, labels = xlabs, name = "Species") +
  scale_fill_manual(  values = pal, guide = "none") +
  labs(x = "Experiment day", y = "Daily hazard", title = "a) Daily larval mortality hazard") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        legend.text  = element_text(face = "italic"),
        plot.title   = element_text(face = "bold"))

# Marginal temperature effect on mortality: vary temp, hold day at its mean
plot_mort_env <- function(sc_col, orig_col, xlab, panel_title) {
  g  <- seq(min(dt_mortality[[sc_col]], na.rm = TRUE),
            max(dt_mortality[[sc_col]], na.rm = TRUE), length.out = 100)
  nd <- expand.grid(setNames(list(g), sc_col)) %>%
    expand_grid(sp = factor(sp_levels, levels = sp_levels)) %>%
    mutate(day_scaled = mean(dt_mortality$day_scaled), ID = NA_character_)
  pf <- fitted(mod_mortality, newdata = nd,
               allow_new_levels = TRUE, re_formula = NA) %>%
    as.data.frame() %>% bind_cols(nd) %>%
    rename(estimate = Estimate, lower = Q2.5, upper = Q97.5) %>%
    mutate(
      x_orig = .data[[sc_col]] * sd(dt_mortality[[sc_col]], na.rm = TRUE) +
        mean(dt_mortality[[orig_col]], na.rm = TRUE),
      sp     = as.character(sp)
    )
  ggplot(pf, aes(x = x_orig, y = estimate, colour = sp, fill = sp)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = pal, labels = xlabs, name = "Species") +
    scale_fill_manual(  values = pal, guide = "none") +
    labs(x = xlab, y = "Daily hazard", title = panel_title) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom",
          legend.text  = element_text(face = "italic"),
          plot.title   = element_text(face = "bold"))
}

p_mort_t <- plot_mort_env("temp_1d_sc", "temp_1d",
                          "Previous-day max temp (°C)", "b) Temperature effect")

fig1 <- (p_haz) / (p_mort_t) +
  plot_annotation(
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )
ggsave("plots/Figure_Mortality_exp3.png", fig1, width = 11, height = 9, dpi = 300)

# ── Figure 2: Duration (species effect) ───────────────────────────────────────
p2 <- plot_sp(mod_pup,   df_pup,   "days_P",     "Days", "a) Days to pupation")
p3 <- plot_sp(mod_eme,   df_eme,   "days_I",     "Days", "b) Days pupa to imago")
p4 <- plot_sp(mod_total, df_total, "days_total", "Days", "c) Total days larva to imago")

fig2 <- (p2 | p3 | p4) +
  plot_annotation(
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )
ggsave("plots/Figure_Duration_exp3.png", fig2, width = 13, height = 5, dpi = 300)

# ── Figure 3: Morphometrics (species effect) ──────────────────────────────────

p5a <- plot_sp(mod_bl, df_morph, "body_length", "Body Length (mm)",
               "a) Body length")
p5b <- plot_sp(mod_wl, df_morph, "wing_length", "Wing Length (mm)",
               "b) Front wing length")
p5c <- plot_sp(mod_wt, df_morph, "weight_I",    "Weight (g)",
               "c) Adult weight")

fig3 <- (p5a | p5b | p5c) +
  plot_annotation(
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )
ggsave("plots/Figure_Morphometrics_exp3.png", fig3, width = 13, height = 5, dpi = 300)

# ── Figure 4: Posterior predictive checks (all models) ───────────────────────
make_ppc <- function(model, title) {
  pp_check(model, ndraws = 100) +
    labs(title = title) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10))
}

ppc_list <- list(
  make_ppc(mod_mortality, "1. Mortality"),
  make_ppc(mod_pup,       "2. Days to Pupation"),
  make_ppc(mod_eme,       "3. Days to Emergence"),
  make_ppc(mod_total,     "4. Total Days"),
  make_ppc(mod_bl,        "5a. Body Length"),
  make_ppc(mod_wl,        "5b. Wing Length"),
  make_ppc(mod_wt,        "5c. Weight")
)

fig4 <- wrap_plots(ppc_list, ncol = 3) +
  plot_annotation(
    title = "Posterior Predictive Checks",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )
ggsave("plots/Figure_PPC_exp3.png", fig4, width = 13, height = 13, dpi = 300)

# ── Figure 5: Rhat convergence diagnostics ────────────────────────────────────
fig5 <- ggplot(rhat_df, aes(x = rhat, fill = model)) +
  geom_histogram(binwidth = 0.002, colour = "white", alpha = 0.85) +
  geom_vline(xintercept = 1.01, linetype = "dashed", colour = "red") +
  facet_wrap(~model, scales = "free_y", ncol = 3) +
  labs(x       = expression(hat(R)),
       y       = "Count",
       title   = "Model Convergence Diagnostics",
       caption = "Dashed line: Rhat = 1.01 threshold.") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14))
ggsave("plots/Figure_Rhat_exp3.png", fig5, width = 13, height = 8, dpi = 300)

# ── Figure 6: MCMC trace plots ────────────────────────────────────────────────

fig6 <- wrap_plots(trace_list, ncol = 2) +
  plot_annotation(
    title = "MCMC Trace Plots (population-level & variance parameters)",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )
ggsave("plots/Figure_Trace_exp3.png", fig6, width = 20, height = 18, dpi = 300)
