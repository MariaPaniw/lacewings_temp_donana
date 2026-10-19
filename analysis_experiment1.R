# =============================================================================
# EXPERIMENT 1: Myrmeleon Life Cycle Dynamics — Bayesian Models (two species)
# =============================================================================
# Models:
#  1. Larval mortality    — discrete-time hazard (Bernoulli cloglog)
#                          formula: s(day_scaled, by = trt_sp) + sp + (1|ID)
#  2. Days to pupation    — lognormal: treatment * sp
#  3. Days pupa to imago  — lognormal: treatment * sp
#  4. Total days to imago — lognormal: treatment * sp
#  5. Body length, wing length, weight — Gaussian: treatment * sp
#
# Species: Myrmeleon_almohadarum, Myrmeleon_hyalinus
# Treatment: 29°C vs 33°C
# All models include species (sp) as a fixed effect. Model 1 uses a
# treatment:species smooth (trt_sp); models 2-5 use a treatment * sp
# interaction, so each fits a treatment effect, a species effect, and their
# interaction.
#
# Required file in the working directory:
#   Data_2_2024_days_to.csv    (individual summary — event times, morphometrics)
# =============================================================================

# ── 0. Packages ──────────────────────────────────────────────────────────────
library(brms)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(bayesplot)

# set.seed(42)

# ── 1. Load Data ─────────────────────────────────────────────────────────────
setwd("LINK TO REPOSITORY")


# Clear stale model caches
file.remove(list.files(getwd(), pattern = "^model_.*exp1\\.rds$", full.names = TRUE))

n_days   <- 76L # length of experiment

# Summary CSV — event timing and morphometrics
dat_summ <- read.csv("Data_2_2024_days_to.csv",
                     fileEncoding = "latin1", check.names = FALSE)
names(dat_summ)  <- gsub("\n", " ", names(dat_summ))
names(dat_summ)[names(dat_summ) == "t (\xa1C)"] <- "treatment"
names(dat_summ)[names(dat_summ) == "t (¡C)"]    <- "treatment"
dat_summ <- dat_summ[, !is.na(names(dat_summ)) & names(dat_summ) != ""]
dat_summ <- dat_summ %>%
  filter(!is.na(treatment), !is.na(sp)) %>%
  mutate(treatment = factor(as.character(as.integer(treatment)),
                            levels = c("29", "33")),
         sp        = factor(sp))

cat(sprintf("Individuals: %d\n", nrow(dat_summ)))
cat("Species x Treatment:\n"); print(table(dat_summ$sp, dat_summ$treatment))

# ── 2. Prepare Analysis Datasets ──────────────────────────────────────────────

# --- 2a. Larval mortality discrete-time data ----------------------------------
# Build one row per individual-day-at-risk.
# At risk until death (event=1) or pupation (right-censored at pupa day - 1).
# Individuals that neither died nor pupated are censored at n_days (none in this
# dataset, but handled for robustness).

mort_events <- dat_summ %>%
  mutate(
    day_death  = as.integer(Numb_days_before_death),
    day_pupate = as.integer(Number_days_P),
    end_day    = case_when(
      !is.na(as.integer(Numb_days_before_death)) ~ as.integer(Numb_days_before_death),
      !is.na(as.integer(Number_days_P))           ~ as.integer(Number_days_P) - 1L,
      TRUE                                         ~ as.integer(n_days)
    )
  )

dt_mortality <- mort_events %>%
  filter(end_day >= 1) %>%
  rowwise() %>%
  reframe(
    ID        = ID,
    treatment = treatment,
    sp        = sp,
    Stage     = Stage,
    day       = seq_len(end_day),
    event     = as.integer(seq_len(end_day) == day_death & !is.na(day_death))
  ) %>%
  mutate(day_scaled = day / n_days,                          # scale to [0,1]
         trt_sp     = interaction(treatment, sp, drop = TRUE)) # combined for by= smooth

cat(sprintf("Mortality: %d rows, %d events\n",
            nrow(dt_mortality), sum(dt_mortality$event)))

# --- 2b. Lognormal duration datasets ------------------------------------------
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

# --- 2c. Morphometric datasets -------------------------------------------------
# Rename columns once for tidiness
dat_summ <- dat_summ %>%
  rename_with(~ "body_length", matches("Body_length_I"))  %>%
  rename_with(~ "wing_length", matches("Length_front_"))  %>%
  rename_with(~ "weight_I",    matches("Weight_I"))

df_morph <- dat_summ %>%
  filter(!is.na(body_length) | !is.na(wing_length) | !is.na(weight_I))
cat(sprintf("Morphometrics n=%d\n", nrow(df_morph)))

# ── 3. Priors & MCMC Settings ─────────────────────────────────────────────────

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

# Lognormal priors for all duration models (2, 3, 4).
# lognormal() models log(y) ~ Normal, handling right-skewed positive data
# without requiring an explicit dispersion parameter. sigma is on the log scale.
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

# ── 4. Fit Models ─────────────────────────────────────────────────────────────

# --- Model 1: Larval mortality (discrete-time hazard) ------------------------
# Species x treatment interaction modelled via s(day_scaled, by = trt_sp):
# a separate smooth baseline hazard for each treatment x species combination.
# sp also enters as a fixed main effect to capture species intercept differences.
cat("\n--- Model 1: Larval mortality ---\n")

clear_bad_cache("model_Mortality_exp1.rds")

mod_mortality <- brm(
  formula = event ~ s(day_scaled, by = trt_sp, bs = "cr", k = 4) + sp + (1 | ID),
  data    = dt_mortality,
  family  = bernoulli(link = "cloglog"),
  prior   = hazard_priors,
  chains  = n_chains, iter = n_iter, warmup = n_warmup, cores = 4,
  seed    = 42, init = 0,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  file    = "model_Mortality_exp1"
)

# --- Model 2: Days to pupation (lognormal) -----------------------------------
cat("\n--- Model 2: Days to pupation (lognormal) ---\n")

clear_bad_cache("model_Days_Pupation_exp1.rds")

mod_pup <- brm(
  formula = days_P ~ treatment * sp,
  data    = df_pup,
  family  = lognormal(),
  prior   = lognorm_priors,
  chains  = n_chains, iter = n_iter, warmup = n_warmup, cores = 4,
  seed    = 42, init = 0,
  control = list(adapt_delta = 0.95),
  file    = "model_Days_Pupation_exp1"
)

# --- Model 3: Days pupa to imago (lognormal) ---------------------------------
cat("\n--- Model 3: Days pupa to imago (lognormal) ---\n")

clear_bad_cache("model_Days_Emergence_exp1.rds")

mod_eme <- brm(
  formula = days_I ~ treatment * sp,
  data    = df_eme,
  family  = lognormal(),
  prior   = lognorm_priors,
  chains  = n_chains, iter = n_iter, warmup = n_warmup, cores = 4,
  seed    = 42, init = 0,
  control = list(adapt_delta = 0.95),
  file    = "model_Days_Emergence_exp1"
)

# --- Model 4: Total days larva to imago (lognormal) -------------------------
cat("\n--- Model 4: Total days to imago (lognormal) ---\n")

clear_bad_cache("model_Days_Total_exp1.rds")

mod_total <- brm(
  formula = days_total ~ treatment * sp,
  data    = df_total,
  family  = lognormal(),
  prior   = lognorm_priors,
  chains  = n_chains, iter = n_iter, warmup = n_warmup, cores = 4,
  seed    = 42, init = 0,
  control = list(adapt_delta = 0.95),
  file    = "model_Days_Total_exp1"
)

# --- Models 5a-c: Morphometrics (Gaussian) ------------------------
fit_morph <- function(var, label, data) {
  cat(sprintf("\n--- Model 5: %s ---\n", label))
  
  clear_bad_cache(paste0("model_", label, "_exp1.rds"))
  
  d <- data %>% filter(!is.na(.data[[var]]))
  brm(
    formula = reformulate("treatment * sp", response = var),
    data    = d,
    family  = gaussian(),
    chains  = n_chains, iter = n_iter, warmup = n_warmup, cores = 4,
    seed    = 42, init = 0,
    control = list(adapt_delta = 0.95),
    file    = paste0("model_", label, "_exp1"),
    file_refit = "on_change"  # refit automatically since the formula/prior
    # changed from the previous z-scored version
  )
}

mod_bl <- fit_morph("body_length", "Body_length",  df_morph)
mod_wl <- fit_morph("wing_length", "Wing_length",  df_morph)
mod_wt <- fit_morph("weight_I",    "Weight_imago", df_morph)

# ── 5. Convergence Diagnostics ────────────────────────────────────────────────

# --- 5a. Rhat summary ---------------------------------------------------------
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

# --- 5b. Trace plots ----------------------------------------------------------
# Restricted to population-level effects (b_), residual/family sigma,
# and group-level SDs (sd_/sds_) rather than every individual ID-level random
# effect (r_ID...) — plotting a trace panel per larva would produce an
# unreadable wall of panels and adds little beyond the SD of that random
# effect, which is already included.
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

# ── 6. Posterior Summaries ────────────────────────────────────────────────────

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

# ── 7. Posterior Significance Check ──────────────────────────────────────────
# Unlike experiment 2 (a single treatment coefficient per model), every model
# here has multiple population-level terms — treatment, sp, and (for models
# 2-5) their interaction — so all non-intercept b_ coefficients are reported,
# not just one.

cat("\n=== Key Fixed-Effect Posteriors ===\n")
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

# ── 8. Summary Tables ──────────────────────────────────────────────────────────
# Writes two CSVs to the project directory (still the working directory at
# this point — the setwd() to /plots happens in Section 9 below):
#   Table_ModelSummaries_exp1.csv — every estimated parameter, all 7 models
#   Table_FixedEffects_exp1.csv   — every non-intercept b_ coefficient per
#                                    model (treatment, sp, and their
#                                    interaction where present), the tabulated
#                                    version of Section 7's console output

# --- 8a. Full parameter table --------------------------------------------------
# Pulls every parameter block brms reports (population-level effects,
# group-level/smooth SDs, and the family-specific parameter — sigma for the
# lognormal/Gaussian models) into one long-format table with the same
# Estimate/CI/Rhat/ESS columns brms itself uses, so nothing is re-derived.
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

write_csv(summary_table, "Table_ModelSummaries_exp1.csv")

# --- 8b. Fixed-effect table -----------------------------------------------------
# Same numbers as Section 7's console output, captured as a data frame per
# coefficient instead of just printed
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

write_csv(fixed_effect_table, "Table_FixedEffects_exp1.csv")

# ── 9. Figures ────────────────────────────────────────────────────────────────
setwd(".../plots")


pal_trt <- c("29" = "#2A9D8F", "33" = "#B2182B")
trt_labels <- c("29" = "29°C", "33" = "33°C")
sp_labels <- c("Myrmeleon_almohadarum" = "M. almohadarum",
               "Myrmeleon_hyalinus"    = "M. hyalinus")

sp_levels  <- levels(dat_summ$sp)
trt_levels <- levels(dat_summ$treatment)

# Helper: posterior mean ± 95% CI + jitter, faceted by species
plot_treatment <- function(model, data, resp_var, ylab, title,
                           back_transform = NULL) {
  nd  <- expand.grid(
    treatment = factor(trt_levels, levels = trt_levels),
    sp        = factor(sp_levels,  levels = sp_levels)
  )
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
                aes(x = treatment, y = .data[[resp_var]], colour = treatment),
                width = 0.12, alpha = 0.45, size = 1.8) +
    geom_pointrange(data = pf,
                    aes(x = treatment, y = est, ymin = lo, ymax = hi,
                        colour = treatment),
                    size = 1.0, linewidth = 1.2, shape = 18) +
    facet_wrap(~sp, labeller = as_labeller(sp_labels)) +
    scale_colour_manual(values = pal_trt, labels = trt_labels, name = "Treatment",guide = "none") +
    scale_x_discrete(labels = trt_labels) +
    labs(x = "Temperature treatment", y = ylab, title = title) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12),
          strip.text = element_text(face = "italic"))
}

# ── Figure 1: Larval mortality hazard curves (2 treatment colours, faceted by species) ──
scaled_grid <- seq(min(dt_mortality$day_scaled),
                   max(dt_mortality$day_scaled), length.out = 150)
nd_mort <- expand.grid(
  day_scaled = scaled_grid,
  treatment  = factor(trt_levels, levels = trt_levels),
  sp         = factor(sp_levels,  levels = sp_levels)
) %>%
  mutate(trt_sp = interaction(treatment, sp, drop = TRUE),
         ID     = NA_character_)

pred_mort <- fitted(mod_mortality, newdata = nd_mort,
                    allow_new_levels = TRUE, re_formula = NA) %>%
  as.data.frame() %>% bind_cols(nd_mort) %>%
  rename(estimate = Estimate, lower = Q2.5, upper = Q97.5) %>%
  mutate(day_orig = day_scaled * n_days)

p_haz <- ggplot(pred_mort, aes(x = day_orig, y = estimate,
                               colour = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~sp, labeller = as_labeller(sp_labels)) +
  scale_colour_manual(values = pal_trt, labels = trt_labels, name = "Treatment") +
  scale_fill_manual(  values = pal_trt, guide = "none") +
  labs(x = "Experiment day", y = "Daily mortality hazard") +
  theme_classic(base_size = 12) +
  theme(legend.position = c(0.1,0.8), plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "italic"))

fig1 <- p_haz +
  plot_annotation(
    theme   = theme(plot.title = element_text(face = "bold", size = 14))
  )
ggsave("Figure_Mortality_exp1.png", fig1, width = 10, height = 6, dpi = 300)

# ── Figure 2: Lognormal duration models (2, 3, 4) ─────────────────────────────
p2 <- plot_treatment(mod_pup,   df_pup,   "days_P",     "Days",
                     "a) Days to pupation")
p3 <- plot_treatment(mod_eme,   df_eme,   "days_I",     "Days",
                     "b) Days pupa to imago")
p4 <- plot_treatment(mod_total, df_total, "days_total", "Days",
                     "c) Total days larva to imago")

fig2 <- (p2 | p3 | p4) +
  plot_annotation(
    theme   = theme(plot.title = element_text(face = "bold", size = 14))
  )
ggsave("Figure_Duration_exp1.png", fig2, width = 13, height = 5, dpi = 300)

# ── Figure 3: Morphometrics (models 5a–c) ─────────────────────────────────────

p5a <- plot_treatment(mod_bl, df_morph, "body_length", "Body Length (mm)",
                      "a) Body length")
p5b <- plot_treatment(mod_wl, df_morph, "wing_length", "Wing Length (mm)",
                      "b) Front wing length")
p5c <- plot_treatment(mod_wt, df_morph, "weight_I",    "Weight (g)",
                      "c) Adult weight")

fig3 <- (p5a | p5b | p5c) +
  plot_annotation(
    theme   = theme(plot.title = element_text(face = "bold", size = 14))
  )
ggsave("Figure_Morphometrics_exp1.png", fig3, width = 13, height = 5, dpi = 300)

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
ggsave("Figure_PPC_exp1.png", fig4, width = 13, height = 13, dpi = 300)

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
ggsave("Figure_Rhat_exp1.png", fig5, width = 13, height = 8, dpi = 300)

# ── Figure 6: MCMC trace plots ────────────────────────────────────────────────

fig6 <- wrap_plots(trace_list, ncol = 2) +
  plot_annotation(
    title = "MCMC Trace Plots (population-level & variance parameters)",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )
ggsave("Figure_Trace_exp1.png", fig6, width = 14, height = 18, dpi = 300)