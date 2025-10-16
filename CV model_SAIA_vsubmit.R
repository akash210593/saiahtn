## Author: Akash Malhotra
## Last updated: September 2025

rm(list=ls())

library(DiagrammeR)
library(data.table)  # for fast data handling
library(dplyr)
library(tidyr)
library(triangle)
library(extraDistr)
library(progress)

grViz("
digraph MarkovModel {
  rankdir=TB;

  /* IHD Pathway */
  subgraph cluster_IHD {
    label=\"IHD Pathway (First MI Only)\";
    style=filled;
    color=lightgrey;
    node [shape=box, style=filled, fillcolor=white];

    A [label=\"No CVD\"];
    B [label=\"Acute MI (Tunnel)\"];
    C [label=\"Chronic IHD\"];
    D [label=\"Death (IHD/MI)\"];
    E [label=\"Death (Background)\"];

    A -> B [label=\"Acute MI\"];
    A -> A [label=\"Remain\"];
    A -> E [label=\"Background Mortality\"];

    B -> C [label=\"Survive MI\"];
    B -> D [label=\"Direct Death\"];

    C -> C [label=\"Remain\"];
    C -> D [label=\"IHD Death\"];
    C -> E [label=\"Background Mortality\"];
  }

  /* Stroke Pathway */
  subgraph cluster_Stroke {
    label=\"Stroke Pathway (First Stroke Only)\";
    style=filled;
    color=lightblue;
    node [shape=box, style=filled, fillcolor=white];

    F [label=\"No CVD\"];
    G [label=\"Acute Stroke (Tunnel)\"];
    H [label=\"Post-Stroke\"];
    I [label=\"Death (Stroke)\"];
    J [label=\"Death (Background)\"];

    F -> G [label=\"Acute Stroke\"];
    F -> F [label=\"Remain\"];
    F -> J [label=\"Background Mortality\"];

    G -> H [label=\"Survive Stroke\"];
    G -> I [label=\"Direct Death\"];

    H -> H [label=\"Remain\"];
    H -> I [label=\"Stroke Death\"];
    H -> J [label=\"Background Mortality\"];
  }
}
")


# --- Debug Mode Flag ---
replicate_patients <- 1  # Simulate 10x each patient to get whole number events

debug_mode <- FALSE# Set to FALSE when running full 10k simulations

trace_storage <- list()  # To store individual patient state history


# --- Conditional Seed and Draw Setup ---
if (debug_mode) {
  seeds <- c(125634)
  draws <- c(1)
} else {
  set.seed(2025)
  seeds <- sample(100000:999999, 1, replace = FALSE)
  draws <- 1:100
}


transition_probs <- tibble::tibble(
  transition = c(
    "p_mi_direct_death", "p_chronic_ihd_death",
    "p_stroke_direct_death", "p_post_stroke_death",
    "p_background_death"
  ), ##https://pmc.ncbi.nlm.nih.gov/articles/PMC6592597/
  ## Acute MI and stroke lasts 1 month so we could use annual mortlaity as monthly?
  ## General exam we confirmed that these rates are double for PLHIV?
  ## From same paper, for every 5 acute MI outcomes, theres roughly 1 cardiac arrest
  mean = c((0.65), 0.0034, 0.38, 0.0043, 0.0007), ## 0.86% BG mortlaity per year, could adjust to HIV only
  min  = c((0.50), 0.0017, 0.19, 0.00215, 0.00035), ## 50%
  max  = c((0.75), 0.0051, 0.57, 0.00645, 0.00105) ## 150% ## using formula p_monthly = 1 - (1 - p_annual)^(1/12). where range not give, taken extremes close to mean
  #mean = c(0.05, 0.0034, 0.38, 0.0043, 0.0007), ## 0.86% BG mortlaity per year, could adjust to HIV only
  #min  = c(0.05, 0.0034, 0.38, 0.0043, 0.0007), ## 50%
  #max  = c(0.05, 0.0034, 0.38, 0.0043, 0.0007) ## 150% ## using formula p_monthly = 1 - (1 - p_annual)^(1/12). where range not give, taken extremes close to mean
 # For MI mortality https://www.ncbi.nlm.nih.gov/books/NBK459269/ , https://www.ahajournals.org/doi/pdf/10.1161/01.cir.85.6.2100, https://www.ahajournals.org/doi/10.1161/01.cir.96.11.3849 
) ## bg death was too high so divided by 10 further. Need to rationalize --> are people likely of just dying than havign MI/stroke in the first place

# Extract fixed transition probabilities
p_mi_direct_death      <- transition_probs$mean[transition_probs$transition == "p_mi_direct_death"]
p_chronic_ihd_death    <- transition_probs$mean[transition_probs$transition == "p_chronic_ihd_death"]
p_stroke_direct_death  <- transition_probs$mean[transition_probs$transition == "p_stroke_direct_death"]
p_post_stroke_death    <- transition_probs$mean[transition_probs$transition == "p_post_stroke_death"]
p_background_death     <- transition_probs$mean[transition_probs$transition == "p_background_death"]

#mm_red <- 0.1/5 ## 10% reduction in risk per 5mm Hg reduction in bp

# --- Read and Clean the Data ---
df <- read.csv("csv for model.csv")

library(dplyr)

# Step 1: Calculate effbp
df <- df %>%
  group_by(caseid, Phase) %>%
  mutate(
    effbp = case_when(
      Phase == "baseline" ~ first(na.omit(tasist)),
      Phase == "intensive" ~ mean(tasist, na.rm = TRUE),
      Phase == "sustainment" ~ mean(tasist, na.rm = TRUE),
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup()

# Step 2: Get baseline effbp per patient
baseline_effbp <- df %>%
  filter(Phase == "baseline") %>%
  group_by(caseid) %>%
  summarise(baseline_effbp = first(na.omit(effbp))) %>%
  ungroup()


# Step 3: Join baseline effbp back to full dataset
df <- df %>%
  left_join(baseline_effbp, by = "caseid")

# Step 4: Calculate sbp_red = baseline_effbp - effbp
df <- df %>%
  mutate(sbp_red = baseline_effbp - effbp)

# Step 0: Extract baseline cvd_risk for each patient
baseline_risk <- df %>%
  filter(Phase == "baseline") %>%
  group_by(caseid) %>%
  summarise(baseline_cvd_risk = first(na.omit(cvd_risk))*1.5, .groups = "drop") ## 50% higher risk in PLHIV

# Step 1: Join it back to the main df
df <- df %>%
  left_join(baseline_risk, by = "caseid")

# Define percent reductions by phase and arm
reduction_lookup <- data.frame(
  Phase = c("intensive", "sustainment", "intensive", "sustainment"),
  f.arm = c("Control", "Control", "Intervention", "Intervention"),
  reduction_pct = c(0.2083, 0.2930, 0.2946, 0.4027)
)

# Join and apply fixed reduction
df <- df %>%
  left_join(reduction_lookup, by = c("Phase", "f.arm")) %>%
  mutate(
    adj_risk = case_when(
      Phase == "baseline" ~ baseline_cvd_risk,
      TRUE ~ baseline_cvd_risk * (1 - reduction_pct)
    )
  )




## We have now made hte below to be a snapshot, so whaterver the age, sex, smoking status, etc. is at baseline, we consider that
## This is like saying if the same patient happened to start at another phase or arm on that day, how would they have behaved differently
cleaned_df <- df %>%
  group_by(caseid, Phase, f.arm) %>%
  summarise(
    max_cvd_risk = max(adj_risk, na.rm = TRUE),
    total_sup_cost = (first(sup_cost) *
                        ifelse(f.arm[1] == "Control", 4.79, 4.25) / 30) + ifelse(Phase[1] == "baseline", 0.033, 0.12)*1.074, ## 7.4% facility indirect
    sexo = first(sexo),
    age = first(na.omit(age)),
    cluster = first(cluster),
    .groups = "drop"
  ) #4.79, 4.25


#  )





## effectively: screened*eligible*prescribed
## prescribed control baseline: 36%*9%*66% = 
## prescribed intervention baseline: 36%*8%*74% =
## prescribed control sustainment: 100% * 14% *59% =  
## prescribed intervention sustianment: 100% * 11% * 72% = 
# First, get the baseline age for each caseid
baseline_ages <- cleaned_df %>%
  filter(Phase == "baseline") %>%
  select(caseid, baseline_age = age)

# Join the baseline age back to the full dataset and replace the age column
cleaned_df <- cleaned_df %>%
  left_join(baseline_ages, by = "caseid") %>%
  mutate(age = baseline_age) %>%
  select(-baseline_age)






cvd_phase_avg <- cleaned_df %>%
  select(caseid, Phase, f.arm, max_cvd_risk)

# ---- Fast per-person event counts (no traces) ----
# Uses your existing time_horizon_months, and split 55% MI / 45% stroke
time_horizon_months <- 120

per_patient_events <- cvd_phase_avg %>%          # or cvd_phase_avg_expanded if you really need replicates
  mutate(
    # annual risk corresponding to your 10-year baseline risk adj (already applied in max_cvd_risk)
    annual_risk   = 1 - (1 - (max_cvd_risk / 100))^(1/10),
    monthly_risk  = 1 - (1 - annual_risk)^(1/12),
    mi_monthly    = monthly_risk * 0.55, ## 29336, 35,337 IHD incidence
    stroke_monthly= monthly_risk * 0.45
  ) %>%
  transmute(
    caseid, Phase, f.arm,
    acute_mi_count        = time_horizon_months * mi_monthly,
    acute_stroke_count    = time_horizon_months * stroke_monthly,
    mi_deaths_count       = time_horizon_months * mi_monthly  * p_mi_direct_death,
    stroke_deaths_count   = time_horizon_months * stroke_monthly * p_stroke_direct_death
  )

# Per-person means by Phase × Arm
per_person_means <- per_patient_events %>%
  group_by(Phase, f.arm) %>%
  summarise(
    acute_mi_pp       = mean(acute_mi_count, na.rm = TRUE),
    acute_stroke_pp   = mean(acute_stroke_count, na.rm = TRUE),
    mi_deaths_pp      = mean(mi_deaths_count, na.rm = TRUE),
    stroke_deaths_pp  = mean(stroke_deaths_count, na.rm = TRUE),
    n_patients        = n_distinct(caseid),
    .groups = "drop"
  )

# (optional) Per 100k people over the horizon
per_100k <- per_person_means %>%
  mutate(
    acute_mi_per100k      = acute_mi_pp * 100000,
    acute_stroke_per100k  = acute_stroke_pp * 100000,
    mi_deaths_per100k     = mi_deaths_pp * 100000,
    stroke_deaths_per100k = stroke_deaths_pp * 100000
  )

# If you want a single export:
write.csv(per_100k, "per_person_event_counts.csv", row.names = FALSE)

cvd_phase_avg_expanded <- cvd_phase_avg[rep(1:nrow(cvd_phase_avg), each = replicate_patients), ]
cvd_phase_avg_expanded$replicate_id <- rep(1:replicate_patients, times = nrow(cvd_phase_avg))
cvd_phase_avg_expanded$simid <- paste0(cvd_phase_avg_expanded$caseid, "_",
                                       cvd_phase_avg_expanded$Phase, "_",
                                       cvd_phase_avg_expanded$f.arm, "_",
                                       cvd_phase_avg_expanded$replicate_id)

no_cvd_cost_table <- cleaned_df %>%
  group_by(Phase, f.arm) %>%
  summarise(
    mean_cost = mean(total_sup_cost, na.rm = TRUE),
    min_cost = min(total_sup_cost, na.rm = TRUE),
    max_cost = max(total_sup_cost, na.rm = TRUE),
    .groups = "drop"
  )

# Using mean of https://pmc.ncbi.nlm.nih.gov/articles/PMC6592597/ and https://pmc.ncbi.nlm.nih.gov/articles/PMC3973979/ and then reduced
# Mi costs by 32% and stroke costs by 64%. --> maye stroke costs be higher
param_table <- data.frame(
  state = c("Acute_MI", "Chronic_IHD", "MI_Death",  # removed Heart Failure as discussed
            "Acute_Stroke", "Post_Stroke", "Stroke_Death", "Background_Death"),
  cost_mean = c(1.27*547/5, 1.27*12/5, 0, 1.27*317/5, 1.27*14/5, 0, 0), #from Sali's estimates 201 for post MI and 326 for post stroke
  cost_min  = c(1.27*465/5, 1.27*10/5, 0, 1.27*269/5, 1.27*12/5, 0, 0), # 15% up and down
  cost_max  = c(1.27*629/5, 1.27*14/5, 0, 1.27*365/5, 1.27*16/5, 0, 0) # Should be 567 and adjust, and 30 and 47 and adjust
  #Added 27% inpatient indirect costs

  )




disability_weight_table <- data.frame(
 state = c("Acute_MI", "Chronic_IHD", "Acute_Stroke", "Post_Stroke"),
 dw_mean = c(0.439/12, 0.101/12, 0.92/12, 0.266/12),
 dw_min = c(0.405/12, 0.093/12, 0.782/12, 0.228/12),
 dw_max = c(0.477/12, 0.103/12, 0.99/12, 0.295/12)

)

life_table <- read.csv("life_expectancy_by_age.csv")


# Lookup life expectancy by exact age and sex
get_life_expectancy <- function(age, sex) {
  if (sex == 0) {  # Female
    approx(life_table$age, life_table$female, xout = age, rule = 2)$y
  } else {         # Male
    approx(life_table$age, life_table$male, xout = age, rule = 2)$y
  }
}


rbeta_from_bounds <- function(mean, min, max) {
  # Estimate SD assuming min and max represent 95% confidence interval
  sd_est <- (max - min) / (2 * 1.96)
  var_est <- sd_est^2

  # Check bounds
  if (mean <= 0 || mean >= 1 || var_est <= 0) {
    return(mean)
  }

  # Method of moments
  alpha <- ((1 - mean) / var_est - 1 / mean) * mean^2
  beta  <- alpha * (1 / mean - 1)

  # Return sampled beta value
  if (is.finite(alpha) && is.finite(beta) && alpha > 0 && beta > 0) {
    rbeta(1, alpha, beta)
  } else {
    mean  # fallback if invalid
  }
}

# --- Model Settings ---
time_horizon_months <- 120
discount_rate <- 0.03
discount_factor <- (1 / (1 + discount_rate)^(1 / 12))

# --- Storage for all runs ---
all_patient_traces <- list()  # Store patient-level traces

all_results <- list()
progress_bar <- progress_bar$new(total = length(seeds) * length(draws), format = ":percent [:bar] :current/:total")

for (seed in seeds) {
  set.seed(seed)
  for (draw in draws) {
    n_rows <- nrow(cvd_phase_avg_expanded)
    results <- data.frame(
      caseid = cvd_phase_avg_expanded$caseid,
      Phase = cvd_phase_avg_expanded$Phase,
      f.arm = cvd_phase_avg_expanded$f.arm,
      max_cvd_risk = cvd_phase_avg_expanded$max_cvd_risk,
      total_cost = rep(0, n_rows),
      total_ylds = rep(0, n_rows),
      total_ylls = rep(0, n_rows),
      seed = seed,
      draw = draw
    )

    for (i in 1:n_rows) {
      patient <- cvd_phase_avg_expanded[i, ]

      age_val <- cleaned_df %>% filter(caseid == patient$caseid) %>% pull(age) %>% max()
      sex_val <- cleaned_df %>% filter(caseid == patient$caseid) %>% pull(sexo) %>% first()


      annual_risk <- 1 - (1 - patient$max_cvd_risk / 100)^(1 / 10)
      monthly_risk <- 1 - (1 - annual_risk)^(1 / 12)
      monthly_risk_MI <- monthly_risk * 0.55 #IHME incidence is source
      monthly_risk_Stroke <- monthly_risk * 0.45

      phase_costs <- no_cvd_cost_table %>%
        filter(Phase == patient$Phase, f.arm == patient$f.arm)
      no_cvd_cost <- if (nrow(phase_costs) == 0) 0 else {
        rtriangle(1, phase_costs$min_cost, phase_costs$max_cost, phase_costs$mean_cost)
      }

      cost_draws <- param_table %>%
        mutate(cost_draw = mapply(function(min, max, mode) rtriangle(1, min, max, mode),
                                  cost_min, cost_max, cost_mean))

      dw_draws <- disability_weight_table %>%
        mutate(dw_draw = mapply(rbeta_from_bounds, dw_mean, dw_min, dw_max))

      state_costs <- c(No_CVD = no_cvd_cost, setNames(cost_draws$cost_draw, cost_draws$state))
      state_ylds <- c(No_CVD = 0, setNames(dw_draws$dw_draw, dw_draws$state)) ##

      current_state_MI <- "No_CVD"
      current_state_Stroke <- "No_CVD"
      total_cost <- 0
      total_ylds <- 0
      total_ylls <- 0

      mi_state_history <- rep("No_CVD", time_horizon_months + 1)
      stroke_state_history <- rep("No_CVD", time_horizon_months + 1)

      has_died <- FALSE

      for (month in 1:time_horizon_months) {

        # --- Expected Events this Month ---
        expected_mi        <- monthly_risk_MI
        expected_stroke    <- monthly_risk_Stroke
        expected_bg_death  <- p_background_death

        # --- Expected YLDs ---
        ylds_mi_acute      <- expected_mi * state_ylds["Acute_MI"]
        ylds_stroke_acute  <- expected_stroke * state_ylds["Acute_Stroke"]

        # For chronic conditions, assume survivors go on to live in chronic state each month
        ylds_chronic_ihd   <- expected_mi * (1 - p_mi_direct_death) * state_ylds["Chronic_IHD"]
        ylds_post_stroke   <- expected_stroke * (1 - p_stroke_direct_death) * state_ylds["Post_Stroke"]

        ylds_this_month <- ylds_mi_acute + ylds_stroke_acute + ylds_chronic_ihd + ylds_post_stroke
        total_ylds <- total_ylds + (discount_factor^month) * ylds_this_month


        # --- Expected Costs (SAFE LOOKUPS) ---
        cost_mi_acute <- expected_mi * ifelse("Acute_MI" %in% names(state_costs), state_costs[["Acute_MI"]], 0)
        cost_chronic_ihd <- expected_mi * (1 - p_mi_direct_death) * ifelse("Chronic_IHD" %in% names(state_costs), state_costs[["Chronic_IHD"]], 0)

        cost_stroke_acute <- expected_stroke * ifelse("Acute_Stroke" %in% names(state_costs), state_costs[["Acute_Stroke"]], 0)
        cost_post_stroke <- expected_stroke * (1 - p_stroke_direct_death) * ifelse("Post_Stroke" %in% names(state_costs), state_costs[["Post_Stroke"]], 0)


        background_cost <- ifelse("No_CVD" %in% names(state_costs) && !is.na(state_costs[["No_CVD"]]),
                                  state_costs[["No_CVD"]], 0)

        cost_this_month <- background_cost + cost_mi_acute + cost_chronic_ihd + cost_stroke_acute + cost_post_stroke

        total_cost <- total_cost + (discount_factor^month) * cost_this_month



        # --- Expected YLLs ---
        if (!has_died) {
          expected_death_prob <- expected_bg_death +
            expected_mi * p_mi_direct_death +
            expected_stroke * p_stroke_direct_death +
            expected_mi * (1 - p_mi_direct_death) * p_chronic_ihd_death +
            expected_stroke * (1 - p_stroke_direct_death) * p_post_stroke_death

          le <- get_life_expectancy(age_val, sex_val)
          expected_ylls <- expected_death_prob * le
          total_ylls <- total_ylls + (discount_factor^month) * expected_ylls

          # Optional: cap accumulation if 99.9% likely dead
          if (expected_death_prob > 0.999) {
            has_died <- TRUE
          }
        }


        
      }


      results$total_cost[i] <- total_cost
      results$total_ylls[i] <- total_ylls
      results$total_ylds[i] <- total_ylds
      patient_trace <- data.frame(
        Month = 0:time_horizon_months,
        State_MI = mi_state_history,
        State_Stroke = stroke_state_history,
        simid = patient$simid,
        caseid = patient$caseid,
        Phase = patient$Phase,
        f.arm = patient$f.arm
      )




      }

    all_results[[paste0("seed_", seed, "_draw_", draw)]] <- results
    }



    progress_bar$tick()
  }


# Combine and export (optional)
final_results <- bind_rows(all_results)
write.csv(final_results, "simulated_outputs_det.csv", row.names = FALSE)

cat("Simulation complete! Output saved to 'simulated_outputs_det.csv'\n")

summary_by_group <- final_results %>%
  group_by(seed, draw, Phase, f.arm) %>%
  summarise(
    mean_cost = mean(total_cost, na.rm = TRUE),
    mean_ylds = mean(total_ylds, na.rm = TRUE),
    mean_ylls = mean(total_ylls, na.rm = TRUE),
    mean_dalys = mean_ylds + mean_ylls,
    # Stroke Pathway
        .groups = "drop"
  )

## as an example when comparing arms
final_summary <- summary_by_group %>%
   pivot_wider(
    names_from = f.arm,
    values_from = c(mean_cost, mean_dalys)
  ) %>%
  mutate(
    icer = (mean_cost_Intervention - mean_cost_Control) /
      (mean_dalys_Control - mean_dalys_Intervention)  # DALYs averted
  )



# === 1. Between arms (already done, but looped) ===
cat("=== ICER: Intervention vs Control by Phase ===\n")
for (p in unique(final_summary$Phase)) {
  temp <- final_summary %>% filter(Phase == p)

  mean_icer <- mean(temp$icer, na.rm = TRUE)
  lower_icer <- quantile(temp$icer, probs = 0.025, na.rm = TRUE)
  upper_icer <- quantile(temp$icer, probs = 0.975, na.rm = TRUE)

  cat("\nPhase:", p, "\n")
  cat("Mean ICER (Int vs Ctrl):", round(mean_icer, 2), "USD per DALY averted\n")
  cat("95% UI: [", round(lower_icer, 2), ",", round(upper_icer, 2), "]\n")
}


# === 2. Within-arm comparisons across phases ===
cat("\n=== ICERs: Phase Comparisons Within Each Arm ===\n")

# Loop for each arm
for (arm in c("Control", "Intervention")) {
  cat("\n---", arm, "Arm ---\n")

  base <- summary_by_group %>% filter(Phase == "baseline", f.arm == arm)
  intv <- summary_by_group %>% filter(Phase == "intensive", f.arm == arm)
  sust <- summary_by_group %>% filter(Phase == "sustainment", f.arm == arm)

  # Match by seed and draw
  comp_df <- base %>%
    rename(cost_base = mean_cost, daly_base = mean_dalys) %>%
    left_join(intv %>% rename(cost_intv = mean_cost, daly_intv = mean_dalys),
              by = c("seed", "draw")) %>%
    left_join(sust %>% rename(cost_sust = mean_cost, daly_sust = mean_dalys),
              by = c("seed", "draw"))

  # Calculate ICERs
  comp_df <- comp_df %>%
    mutate(
      icer_intensive = (cost_intv - cost_base) / (daly_base - daly_intv),
      icer_sustainment = (cost_sust - cost_base) / (daly_base - daly_sust)
    )

  # Intensive vs Baseline
  mean_i <- mean(comp_df$icer_intensive, na.rm = TRUE)
  lb_i <- quantile(comp_df$icer_intensive, probs = 0.025, na.rm = TRUE)
  ub_i <- quantile(comp_df$icer_intensive, probs = 0.975, na.rm = TRUE)

  cat("Intensive vs Baseline:\n")
  cat("Mean ICER:", round(mean_i, 2), "USD per DALY averted\n")
  cat("95% UI: [", round(lb_i, 2), ",", round(ub_i, 2), "]\n\n")

  # Sustainment vs Baseline
  mean_s <- mean(comp_df$icer_sustainment, na.rm = TRUE)
  lb_s <- quantile(comp_df$icer_sustainment, probs = 0.025, na.rm = TRUE)
  ub_s <- quantile(comp_df$icer_sustainment, probs = 0.975, na.rm = TRUE)

  cat("Sustainment vs Baseline:\n")
  cat("Mean ICER:", round(mean_s, 2), "USD per DALY averted\n")
  cat("95% UI: [", round(lb_s, 2), ",", round(ub_s, 2), "]\n")
}

library(dplyr)


#write.csv(final_results,"final_results287_det.csv")



arm_phase_n <- tibble(
  Phase = rep(c("baseline", "intensive", "sustainment"), each = 2),
  f.arm = rep(c("Control", "Intervention"), 3),
  N = c(151, 183, 151, 183, 151, 183)
)


population_df <- final_results %>%
  group_by(Phase, f.arm) %>%
  summarise(
    N = n_distinct(caseid),
    .groups = "drop"
  )


#https://www.ncbi.nlm.nih.gov/books/NBK459269/
 # https://www.ahajournals.org/doi/pdf/10.1161/01.cir.85.6.2100
# https://www.ahajournals.org/doi/10.1161/01.cir.96.11.3849

