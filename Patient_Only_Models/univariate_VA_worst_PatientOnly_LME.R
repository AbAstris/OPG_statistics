# ==========================================================
# VA (WORST EYE) vs ln(Time Since Diagnosis) — PATIENTS ONLY
# ln(Time Since Diagnosis) Mixed-Effects Model
#
# Description:
#   This script models the relationship between worst-eye
#   visual acuity (VA) and time since diagnosis in a
#   patients-only cohort, using a linear mixed-effects model
#   with a log-transformed time since diagnosis term.
#   Covariates include NF1 status, sex, and log diagnosis age.
#   A random intercept and random slope for ln(time since
#   diagnosis) are included per patient. Only patients with
#   both >= 2 VA measurements and >= 2 FA scans are included.
#   Worst eye is defined as the highest LogMAR value; single-
#   eye measurements are included where one eye is missing.
#
# Inputs:
#   - VA_file.csv : Patient VA data (see column layout below)
#   - FA_file.csv : Used to identify patients with >= 2 FA scans and look up Sex
#   Rename these to match your actual file names.
#
# Outputs:
#   - Console: patient counts and model summary
#   - VA_worst_vs_TimeSinceDx_LnTimeSinceDx.png    : Marginal VA trajectories by NF1 status
#   - Residuals_VA_worst_TimeSinceDx_LnTimeSinceDx.png : Conditional residuals by NF1 status
#
# Requirements:
#   R packages: readr, dplyr, nlme, ggplot2
#   Install via:
#     install.packages(c("readr", "dplyr", "nlme", "ggplot2"))
#
# Usage:
#   Set your working directory to the folder containing your
#   CSV files before running, e.g.:
#     setwd("/path/to/your/data")
#   Or open this script from an RStudio project.
#
# Author: Emily Drabek-Maunder
# Date:   May 5, 2026
# ==========================================================

# ==========================================================
# 0. LIBRARIES
# ==========================================================
library(readr)
library(dplyr)
library(nlme)
library(ggplot2)

# ==========================================================
# 1. LOAD & PREPARE DATA
# ==========================================================

# Input the VA file, shown here as VA_file.csv with column layout:
# ├── OPG ID                  (patient ID — numeric or string)
# ├── Patient ID              (secondary patient ID)
# ├── DOB                     (date of birth)
# ├── Diagnosis date          (date of diagnosis)
# ├── NF1                     (categorical — "Y" or blank/NA)
# ├── Age at VA (yrs)         (age at VA measurement — numeric, 0–20)
# ├── Diagnosis               (age at diagnosis — numeric)
# ├── Years since diagnosis   (time since diagnosis — numeric)
# ├── Date of VA              (date of VA measurement)
# ├── R eye                   (right eye VA — LogMAR)
# ├── L eye                   (left eye VA — LogMAR)
# ├── Both                    (combined VA — LogMAR)
# └── VA_average              (average VA of both eyes — LogMAR)

# Input the FA file, shown here as FA_file.csv, used to identify
# patients with >= 2 FA scans and to look up Sex:
# ├── Subject                 (patient ID — must match OPG ID in VA file)
# ├── Weighted Avg FA Combined (FA measurement)
# ├── Age                     (age in years)
# └── Sex                     (categorical — "F" or "M")

# ── Load Sex from FA file and retain only patients with >= 2 FA scans ─
fa_eligible <- read_csv("FA_JM_updated.csv", show_col_types = FALSE) %>%
  rename(FA = `Weighted Avg FA Combined`) %>%
  mutate(patient_id = as.character(Subject)) %>%
  filter(!is.na(Age), !is.na(FA), Age >= 0, Age <= 20) %>%
  group_by(patient_id) %>%
  filter(n() >= 2) %>%
  ungroup() %>%
  select(patient_id, Sex) %>%
  distinct(patient_id, .keep_all = TRUE)

# ── Load VA data ───────────────────────────────────────────
df <- read_csv("VA_years_JM.csv", show_col_types = FALSE) %>%
  rename(
    ogp_id        = `OPG ID`,
    patient_id_va = `Patient ID`,
    dob           = DOB,
    diag_date     = `Diagnosis date`,
    NF1           = NF1,
    Age           = `Age at VA (yrs)`,
    Diagnosis_age = Diagnosis,
    time_since_dx = `Years since diagnosis`,
    date_va       = `Date of VA`,
    R_eye         = `R eye`,
    L_eye         = `L eye`,
    VA_both       = Both,
    VA_average    = VA_average
  ) %>%
  mutate(
    patient_id = as.character(ogp_id),
    NF1        = ifelse(is.na(NF1) | NF1 == "", "N", NF1),
    Group      = ifelse(NF1 == "Y", "Patient_NF1", "Patient_nonNF1")
  ) %>%
  # Drop rows where both eyes are missing
  filter(!(is.na(R_eye) & is.na(L_eye))) %>%
  # Worst eye = highest LogMAR; handle single-eye missingness
  mutate(
    VA_worst = case_when(
      !is.na(R_eye) & !is.na(L_eye) ~ pmax(R_eye, L_eye),
      !is.na(R_eye) &  is.na(L_eye) ~ R_eye,
       is.na(R_eye) & !is.na(L_eye) ~ L_eye
    )
  ) %>%
  filter(!is.na(Age), !is.na(VA_worst), Age >= 0, Age <= 20) %>%
  filter(!is.na(Diagnosis_age), !is.na(time_since_dx)) %>%
  # Retain only patients with >= 2 FA scans; join Sex
  left_join(fa_eligible, by = "patient_id") %>%
  filter(!is.na(Sex)) %>%
  group_by(patient_id) %>%
  filter(n() >= 2) %>%
  ungroup() %>%
  mutate(
    Group            = factor(Group, levels = c("Patient_nonNF1", "Patient_NF1")),
    Sex              = factor(Sex, levels = c("F", "M")),
    NF1_bin          = as.numeric(Group == "Patient_NF1"),
    ln_diag_age      = log(Diagnosis_age + 1e-3),
    ln_time_since_dx = log(pmax(time_since_dx, 1e-3))
  ) %>%
  select(patient_id, Age, Sex, Group, NF1_bin,
         VA_worst, Diagnosis_age, ln_diag_age,
         time_since_dx, ln_time_since_dx)

cat("Patients (>=2 VA measurements AND >=2 FA scans):", n_distinct(df$patient_id), "\n")
cat("NF1:                                            ", sum(df$NF1_bin == 1 & !duplicated(df$patient_id)), "\n")
cat("non-NF1:                                        ", sum(df$NF1_bin == 0 & !duplicated(df$patient_id)), "\n")

# ==========================================================
# 2. FIT MODEL — ln(Time Since Diagnosis)
# Random intercept + random slope for ln(time since dx)
# Covariates: NF1 status, sex, log diagnosis age
# ==========================================================

ctrl <- lmeControl(maxIter = 5000, msMaxIter = 5000, msMaxEval = 10000)

fit_log <- lme(
  VA_worst ~ NF1_bin * ln_time_since_dx + Sex + ln_diag_age,
  random  = ~ 1 + ln_time_since_dx | patient_id,
  data    = df,
  method  = "ML",
  control = ctrl
)

cat("\n===== Model: ln(Time Since Diagnosis) — Patients Only =====\n")
print(summary(fit_log))

# ==========================================================
# 3. COLOUR PALETTE
# ==========================================================

palette <- c(
  "Patient_nonNF1" = "#2ca02c",   # green
  "Patient_NF1"    = "#ff7f0e"    # orange
)

# ==========================================================
# 4. VA vs TIME SINCE DIAGNOSIS PLOT
# Marginal predictions over time since dx grid
# Sex fixed at F, ln_diag_age fixed at median
# ==========================================================

tdx_grid    <- seq(min(df$time_since_dx, na.rm = TRUE),
                   max(df$time_since_dx, na.rm = TRUE),
                   length.out = 150)
ln_tdx_grid <- log(pmax(tdx_grid, 1e-3))
med_ln_diag <- median(df$ln_diag_age, na.rm = TRUE)

beta       <- fixef(fit_log)
coef_names <- names(beta)
vcov_mat   <- vcov(fit_log)
groups     <- c("Patient_nonNF1", "Patient_NF1")

pred_list <- lapply(groups, function(grp) {

  nf1 <- as.numeric(grp == "Patient_NF1")

  nd <- data.frame(
    time_since_dx    = tdx_grid,
    ln_time_since_dx = ln_tdx_grid,
    ln_diag_age      = med_ln_diag,
    NF1_bin          = nf1,
    Sex              = factor("F", levels = levels(df$Sex)),
    Group            = factor(grp,  levels = levels(df$Group)),
    patient_id       = "new_subject"
  )

  nd$pred <- predict(fit_log, newdata = nd, level = 0)

  # ── Delta method SE ──────────────────────────────────────
  se_vec <- numeric(nrow(nd))
  for (k in seq_len(nrow(nd))) {
    cvec <- setNames(rep(0, length(beta)), coef_names)
    cvec["(Intercept)"]            <- 1
    if ("NF1_bin"                  %in% coef_names) cvec["NF1_bin"]                  <- nf1
    if ("ln_time_since_dx"         %in% coef_names) cvec["ln_time_since_dx"]         <- nd$ln_time_since_dx[k]
    if ("ln_diag_age"              %in% coef_names) cvec["ln_diag_age"]              <- med_ln_diag
    if ("SexM"                     %in% coef_names) cvec["SexM"]                     <- 0
    if ("NF1_bin:ln_time_since_dx" %in% coef_names) cvec["NF1_bin:ln_time_since_dx"] <- nf1 * nd$ln_time_since_dx[k]
    se_vec[k] <- sqrt(as.numeric(t(cvec) %*% vcov_mat %*% cvec))
  }

  nd$se    <- se_vec
  nd$lower <- nd$pred - 1.96 * nd$se
  nd$upper <- nd$pred + 1.96 * nd$se
  nd
})

pred_df <- bind_rows(pred_list) %>%
  mutate(Group = factor(Group, levels = levels(df$Group)))

p_va <- ggplot() +
  geom_point(
    data  = df,
    aes(x = time_since_dx, y = VA_worst, colour = Group),
    alpha = 0.4, size = 1.5
  ) +
  geom_ribbon(
    data  = pred_df,
    aes(x = time_since_dx, ymin = lower, ymax = upper, fill = Group),
    alpha = 0.2
  ) +
  geom_line(
    data      = pred_df,
    aes(x = time_since_dx, y = pred, colour = Group),
    linewidth = 1.1
  ) +
  scale_colour_manual(values = palette,
                      labels = c("Patient (non-NF1)", "Patient (NF1)")) +
  scale_fill_manual(values   = palette,
                    labels   = c("Patient (non-NF1)", "Patient (NF1)")) +
  labs(
    x      = "Time Since Diagnosis (years)",
    y      = "Worst eye VA (LogMAR)",
    colour = "Group",
    fill   = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 14),
    legend.title     = element_text(size = 15),
    axis.text        = element_text(size = 13),
    axis.title       = element_text(size = 15)
  )

print(p_va)

ggsave(
  filename = "VA_worst_vs_TimeSinceDx.png",
  plot     = p_va,
  width    = 8, height = 6, dpi = 300
)
cat("\nSaved: VA_worst_vs_TimeSinceDx.png\n")

# ==========================================================
# 5. RESIDUAL PLOT
# ==========================================================

df$fitted   <- predict(fit_log, newdata = df, level = 1)
df$residual <- df$VA_worst - df$fitted

p_resid <- ggplot(df, aes(x = time_since_dx, y = residual,
                           colour = Group, fill = Group)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(aes(group = Group),
              method    = "loess",
              se        = TRUE,
              linewidth = 1.1,
              alpha     = 0.2) +
  geom_hline(yintercept = 0,
             linetype   = "dashed",
             colour     = "black") +
  scale_colour_manual(values = palette,
                      labels = c("Patient (non-NF1)", "Patient (NF1)")) +
  scale_fill_manual(values   = palette,
                    labels   = c("Patient (non-NF1)", "Patient (NF1)")) +
  labs(
    x      = "Time Since Diagnosis (years)",
    y      = "Residual",
    colour = "Group",
    fill   = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 14),
    legend.title     = element_text(size = 15),
    axis.text        = element_text(size = 13),
    axis.title       = element_text(size = 15)
  )

print(p_resid)

ggsave(
  filename = "Residuals_VA_worst_TimeSinceDx.png",
  plot     = p_resid,
  width    = 8, height = 6, dpi = 300
)
cat("Saved: Residuals_VA_worst_TimeSinceDx.png\n")