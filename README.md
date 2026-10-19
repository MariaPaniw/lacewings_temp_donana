# lacewings_temp_donana

Data and code accompanying the article **"Warmer conditions shorten pupal development in two Mediterranean antlion species"**.

The repository contains individual-level life-history data for two antlion species (Neuroptera: Myrmeleontidae), *Myrmeleon almohadarum* and *Myrmeleon hyalinus*, from Doñana (SW Spain), together with the R code used to fit Bayesian models of larval mortality, development time and adult size. There are three experiments:

| Experiment | Script | Data file | Design | Species | n | Duration |
|---|---|---|---|---|---|---|
| 1 | `analysis_experiment1.R` | `Data_2_2024_days_to.csv` | Laboratory, constant 29 °C vs 33 °C | *M. almohadarum* (25), *M. hyalinus* (151) | 176 | 76 days |
| 2 | `analysis_experiment2.R` | `Data_3_2024_days_to.csv` | Laboratory, constant 27 °C vs 27–39 °C regime | *M. almohadarum* (104) | 104 | 64 days |
| 3 | `analysis_experiment3.R` | `Data_2025_days_to.csv` + `data_hobo_2025.csv` | Field, natural temperature (22 Apr – 1 Jul 2025) | *M. almohadarum* (75), *M. hyalinus* (95) | 170 | 71 days |



---

## Repository structure

```
lacewings_temp_donana/
├── README.md
├── LICENSE                      # GNU GPL v3
├── analysis_experiment1.R       # Exp. 1: 29 vs 33 °C, two species
├── analysis_experiment2.R       # Exp. 2: 27 vs 27–39 °C, M. almohadarum
├── analysis_experiment3.R       # Exp. 3: field 2025, two species + logger data
├── Data_2_2024_days_to.csv      # individual data, experiment 1
├── Data_3_2024_days_to.csv      # individual data, experiment 2
├── Data_2025_days_to.csv        # individual data, experiment 3
└── data_hobo_2025.csv           # HOBO temperature/humidity logger, experiment 3
```

---

## Data

### Individual data files (`Data_*_days_to.csv`)

Each row is one larva followed from the start of the experiment until death or adult emergence (imago). The three files share the same core columns; they differ slightly in names and in whether a temperature-treatment column is present.

| Column | Exp. 1 / 2 name | Exp. 3 name | Description | Units / values |
|---|---|---|---|---|
| Individual ID | `ID` | `ID` | Unique larva identifier. Exp. 1: `D_*` = *M. hyalinus*, `P_*` = *M. almohadarum*; Exp. 2: `E-*`; Exp. 3: `1`–`170` | text |
| Species | `sp` | `sp` | `Myrmeleon_almohadarum` or `Myrmeleon_hyalinus` | text |
| Larval stage at start | `Stage` | `Stage` | Larval instar when the experiment started | `1`, `2`, `3` |
| Larval weight at start | `Weight (g)` (Exp. 2: `Weight  (g)`, two spaces) | `Weight _L (g)` | Larval body mass at the start | g |
| Larval length at start | `Length (mm)` | `Length_L (mm)` | Larval body length at the start | mm |
| Temperature treatment | `t (¡C)` | — (field) | Exp. 1: `29`, `33`; Exp. 2: `27`, `27-39` | °C |
| Day of death | `Numb_days_before_death` | `Numb_days_before_death` | Experiment day on which the individual died; `NA` if it survived to emergence | days |
| Weight at death | `Weight_dead_L_P (g)` | `Weight_ea_L_P (g)` | Mass of the dead larva/pupa | g |
| Stage at death | `Stage of death` | `Stage of death` | `1L`, `2L`, `3L` = 1st–3rd larval instar; `PP` = prepupa; `P` = pupa | code |
| Days to pupation | `Number_days_P` | `Number_days_P` | Days from start of experiment to pupation | days |
| Days pupa → imago | `Number_days_I` | `Number_days_I` | Days from pupation to adult emergence | days |
| Sex | `Gender` | `Gender` | Sex of emerged adult. Exp. 1/2: `m` = male, `w` = female; Exp. 3: `male`, `female` | code |
| Adult body length | `Body_length_I (mm)` | same | Body length of the emerged adult | mm |
| Adult front-wing length | `Length_front_wing_I (mm)` | same | Length of the forewing | mm |
| Adult weight | `Weight_I (g)` | same | Mass of the emerged adult | g |
| Adult condition | `Condition_I` | same | `dev` = fully developed; `underdev` (Exp. 1) / `undev` (Exp. 2) = underdeveloped | code |

`NA` means the event did not occur or was not measured. For example, an individual that died as a larva has `NA` for all pupal and adult columns.

**Encoding quirks**, which the scripts already handle:
- The files are Latin-1 encoded. The degree sign in the treatment header appears as `t (¡C)`.
- Some headers contain line breaks inside quotes, e.g. `"Body_length_I\n(mm)"`. The scripts replace these with spaces.

**Sample summary** (from the files):

| | Exp. 1 | Exp. 2 | Exp. 3 |
|---|---|---|---|
| Individuals | 176 | 104 | 170 |
| Died (with day recorded) | 62 | 15 | 12 |
| Pupated | 134 | 92 | 159 |
| Emerged as adults | 113 | 76 | 150 |

### Logger file (`data_hobo_2025.csv`)

Readings from a HOBO temperature/relative-humidity logger placed at the field site for experiment 3.

| Column | Description | Units / format |
|---|---|---|
| `Log` | Logger record number | integer |
| `Date_Time` | Timestamp | `MM.DD.YY hh:mm:ss AM/PM` |
| `Temp_C` | Air temperature | °C |
| `RH` | Relative humidity | % |

The logger runs from 11 Apr 2025 to 1 Jul 2025 at 15-minute intervals (7,792 records, 3 missing values). Experiment day 1 is 22 Apr 2025, so the first 11 days predate the experiment and are excluded from the analysis.

---

## Code

The three scripts follow the same structure and fit the same set of seven `brms` models. They differ in the predictors used:

| # | Response | Family | Exp. 1 predictors | Exp. 2 predictors | Exp. 3 predictors |
|---|---|---|---|---|---|
| 1 | Larval mortality (daily) | Bernoulli, cloglog link (discrete-time hazard) | `s(day, by = treatment×species)` + species + `(1 \| ID)` | `s(day, by = treatment)` + `(1 \| ID)` | `s(day)` + `s(previous-day max temp, by = species)` + species + `(1 \| ID)` |
| 2 | Days to pupation | Lognormal | treatment × species | treatment | species |
| 3 | Days pupa → imago | Lognormal | treatment × species | treatment | species |
| 4 | Total days larva → imago | Lognormal | treatment × species | treatment | species |
| 5a | Adult body length | Gaussian | treatment × species | treatment | species |
| 5b | Adult front-wing length | Gaussian | treatment × species | treatment | species |
| 5c | Adult weight | Gaussian | treatment × species | treatment | species |

### Script sections

Experiment 3 has an extra **Section A** before the models. It reads the logger file, computes daily max/min/mean temperature and RH, adds 7- and 14-day rolling means, and produces exploratory figures E1–E4. Then all three scripts run these steps:

1. **Load data.** Read the CSV, clean the column names and set the factor levels.
2. **Prepare datasets.**
   - **Mortality data.** The data are expanded to one row per individual per day at risk, with `event = 1` on the day of death. Individuals that pupate are right-censored on the day before pupation. Day is rescaled to [0, 1] for the splines.
   - **Exp. 3 covariate.** The previous day's maximum logger temperature is joined in and standardised.
   - **Other models.** The duration and morphometric subsets keep individuals with non-missing values.
3. **Priors and MCMC settings.**
   - **Priors.** Mortality model: Normal(0, 0.5) on fixed effects and intercept, half-Normal(0, 1) on spline SDs, half-Normal(0, 0.25) on the individual SD. Lognormal duration models: Normal(0, 0.5) on effects and intercept, half-Normal(0, 1) on σ. Gaussian morphometric models: `brms` defaults.
   - **Sampling.** 4 chains, 6,000 iterations each, of which 2,000 are warm-up; `seed = 42`.
4. **Fit models.** Uses `brm()`, with fitted models cached as `model_<name>_exp<N>.rds`.
5. **Convergence diagnostics.** Computes R̂ for all parameters and trace plots for the population-level and variance parameters.
6. **Posterior summaries.** Prints `summary()` for each model.
7. **Posterior effect checks.** Reports the median, 95% CI and P(β > 0) for each non-intercept coefficient.
8. **Summary tables.** Written as CSV files (see below).
9. **Figures.** Written as PNG files (see below).

### Outputs

| Output | Exp. 1 | Exp. 2 | Exp. 3 |
|---|---|---|---|
| Model caches | `model_*_exp1.rds` | `model_*_exp2.rds` | `model_*_exp3.rds` |
| All-parameter table | `Table_ModelSummaries_exp1.csv` | `Table_ModelSummaries_exp2.csv` | `Table_ModelSummaries_exp3.csv` |
| Effect table | `Table_FixedEffects_exp1.csv` | `Table_TreatmentEffects_exp2.csv` | `Table_FixedEffects_exp3.csv` |
| Mortality hazard figure | `Figure_Mortality_exp1.png` | `Figure_Mortality_exp2.png` | `Figure_Mortality_exp3.png` (hazard over time + temperature effect) |
| Development-time figure | `Figure_Duration_exp1.png` | `Figure_Duration_exp2.png` | `Figure_Duration_exp3.png` |
| Adult size figure | `Figure_Morphometrics_exp1.png` | `Figure_Morphometrics_exp2.png` | `Figure_Morphometrics_exp3.png` |
| Posterior predictive checks | `Figure_PPC_exp1.png` | `Figure_PPC_exp2.png` | `Figure_PPC_exp3.png` |
| R̂ histograms | `Figure_Rhat_exp1.png` | `Figure_Rhat_exp2.png` | `Figure_Rhat_exp3.png` |
| Trace plots | `Figure_Trace_exp1.png` | `Figure_Trace_exp2.png` | `Figure_Trace_exp3.png` |
| Environmental figures | — | — | `FigureE1`–`E4_Env_*.png` |

---

## Reproducing the analysis

**Requirements.**
- A recent version of R, with dplyr ≥ 1.1 for `reframe()` and the packages `brms`, `tidyverse`, `patchwork` and `bayesplot`; experiment 3 also needs `scales`.
- `brms` needs a working Stan backend (`rstan` or `cmdstanr`).

```r
install.packages(c("brms", "tidyverse", "patchwork", "bayesplot", "scales"))
```

**Before running a script:**

1. **Set the working directory.** Replace `setwd("LINK TO REPOSITORY")` near the top of the script with the path to your local copy of this repository.
2. **Set the figure folder.**
3. **Expect a full refit.** Each script deletes its own cached `model_*_exp<N>.rds` files at the start. To reuse cached fits instead, comment out the `file.remove(...)` line.

Each script is self-contained and can be run on its own. Fitting is the slow step: seven models with 4 × 6,000 iterations, using 4 cores. The mortality models in particular can take a while.

---

## License

Code and data are released under the [GNU General Public License v3.0](LICENSE).

## Citation

If you use these data or code, please cite the accompanying article: *Warmer conditions shorten pupal development in two Mediterranean antlion species* (full reference to be added on publication).
