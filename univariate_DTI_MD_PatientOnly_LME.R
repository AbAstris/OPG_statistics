# ==========================================================
# MD vs AGE — PATIENTS ONLY — ln(Age) Mixed-Effects Model
#
# Description:
#   This script models the relationship between Mean
#   Diffusivity (MD) in the Optic Radiations and age in a
#   patients-only cohort, using a linear mixed-effects model
#   with a log-transformed age term. Covariates include NF1
#   status, sex, and log diagnosis age. A random intercept
#   and random slope for ln(age) are included per patient.
#   Patients with fewer than 2 scans are excluded.
#
# Inputs:
#   - MD_file.csv : Patient DTI data (see column layout below)
#   Rename this to match your actual file name.
#
# Outputs:
#   - Console: patient counts and model summary
#   - MD_vs_Age_Patients_LogAge.png    : Marginal MD trajectories by NF1 status
#   - Residuals_Patients_LogAge_MD.png : Conditional residuals by NF1 status
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
# Author: [Your name]
# Date:   [Date]
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

# Input the patient DTI file, shown here as MD_file.csv with column layout:
# ├── Subject                    (patient ID — numeric or string)
# ├── Weighted Avg MD Combined   (MD measurement — numeric, e.g. 0.00065)
# ├── Age                        (age in years — numeric, 0–20)
# ├── Sex                        (categorical — "F" or "M")
# ├── NF1                        (categorical — "Y" or blank/NA)
# ├── Diagnosis_age              (age at diagnosis — numeric)
# └── Years_since_diagnosis      (years since diagnosis — numeric)

df <- read_csv("MD_JM_updated.csv", show_col_types = FALSE) %>%
  rename(MD = `Weighted Avg MD Combined`) %>%
  mutate(
    patient_id = as.character(Subject),
    NF1        = ifelse(is.na(NF1) | NF1 == "", "N", NF1),
    Group      = ifelse(NF1 == "Y", "Patient_NF1", "Patient_nonNF1")
  ) %>%
  filter(!is.na(Age), !is.na(MD), Age >= 0, Age <= 20) %>%
  filter(!is.na(Diagnosis_age), !is.na(Years_since_diagnosis)) %>%
  group_by(patient_id) %>%
  filter(n() >= 2) %>%
  ungroup() %>%
  mutate(
    Group       = factor(Group, levels = c("Patient_nonNF1", "Patient_NF1")),
    Sex         = factor(Sex, levels = c("F", "M")),
    NF1_bin     = as.numeric(Group == "Patient_NF1"),
    ln_age      = log(Age           + 1e-3),
    ln_diag_age = log(Diagnosis_age + 1e-3)
  ) %>%
  select(patient_id, Age, Sex, Group, NF1_bin,
         MD, Diagnosis_age, ln_age, ln_diag_age)

cat("Patients (>=2 scans):", n_distinct(df$patient_id), "\n")
cat("NF1:                 ", sum(df$NF1_bin == 1 & !duplicated(df$patient_id)), "\n")
cat("non-NF1:             ", sum(df$NF1_bin == 0 & !duplicated(df$patient_id)), "\n")

# ==========================================================
# 2. FIT MODEL — ln(Age)
# Random intercept + random slope for ln(age) per patient
# Covariates: NF1 status, sex, log diagnosis age
# ==========================================================

ctrl <- lmeControl(maxIter = 5000, msMaxIter = 5000, msMaxEval = 10000)

fit_log <- lme(
  MD ~ NF1_bin * ln_age + Sex + ln_diag_age,
  random  = ~ 1 + ln_age | patient_id,
  data    = df,
  method  = "ML",
  control = ctrl
)

cat("\n===== Model: ln(Age) — Patients Only =====\n")
print(summary(fit_log))

# ==========================================================
# 3. COLOUR PALETTE
# ==========================================================

palette <- c(
  "Patient_nonNF1" = "#2ca02c",   # green
  "Patient_NF1"    = "#ff7f0e"    # orange
)

# ==========================================================
# 4. Y-AXIS SCALING NOTE
# MD values are multiplied by 1e4 in the plot so that axis
# ticks display as plain numbers e.g. 6.0, 7.5 etc.
# The axis label carries the (x 10^-4 mm^2/s) unit explicitly.
# Residuals are multiplied by 1e5 for the same reason.
# ==========================================================

# ==========================================================
# 5. MD vs AGE PLOT
# Marginal predictions over age grid
# Sex fixed at F, ln_diag_age fixed at median
# ==========================================================

age_grid    <- seq(min(df$Age, na.rm = TRUE),
                   max(df$Age, na.rm = TRUE),
                   length.out = 150)
ln_age_grid <- log(age_grid + 1e-3)
med_ln_diag <- median(df$ln_diag_age, na.rm = TRUE)

beta       <- fixef(fit_log)
coef_names <- names(beta)
vcov_mat   <- vcov(fit_log)
groups     <- c("Patient_nonNF1", "Patient_NF1")

pred_list <- lapply(groups, function(grp) {

  nf1 <- as.numeric(grp == "Patient_NF1")

  nd <- data.frame(
    Age         = age_grid,
    ln_age      = ln_age_grid,
    ln_diag_age = med_ln_diag,
    NF1_bin     = nf1,
    Sex         = factor("F", levels = levels(df$Sex)),
    Group       = factor(grp,  levels = levels(df$Group)),
    patient_id  = "new_subject"
  )

  nd$pred <- predict(fit_log, newdata = nd, level = 0)

  # ── Delta method SE ──────────────────────────────────────
  se_vec <- numeric(nrow(nd))
  for (k in seq_len(nrow(nd))) {
    cvec <- setNames(rep(0, length(beta)), coef_names)
    cvec["(Intercept)"]  <- 1
    if ("NF1_bin"        %in% coef_names) cvec["NF1_bin"]        <- nf1
    if ("ln_age"         %in% coef_names) cvec["ln_age"]         <- nd$ln_age[k]
    if ("ln_diag_age"    %in% coef_names) cvec["ln_diag_age"]    <- med_ln_diag
    if ("SexM"           %in% coef_names) cvec["SexM"]           <- 0
    if ("NF1_bin:ln_age" %in% coef_names) cvec["NF1_bin:ln_age"] <- nf1 * nd$ln_age[k]
    se_vec[k] <- sqrt(as.numeric(t(cvec) %*% vcov_mat %*% cvec))
  }

  nd$se    <- se_vec
  nd$lower <- nd$pred - 1.96 * nd$se
  nd$upper <- nd$pred + 1.96 * nd$se
  nd
})

pred_df <- bind_rows(pred_list) %>%
  mutate(Group = factor(Group, levels = levels(df$Group)))

p_md <- ggplot() +
  geom_point(
    data  = df,
    aes(x = Age, y = MD * 1e4, colour = Group),
    alpha = 0.4, size = 1.5
  ) +
  geom_ribbon(
    data  = pred_df,
    aes(x = Age, ymin = lower * 1e4, ymax = upper * 1e4, fill = Group),
    alpha = 0.2
  ) +
  geom_line(
    data      = pred_df,
    aes(x = Age, y = pred * 1e4, colour = Group),
    linewidth = 1.1
  ) +
  scale_colour_manual(values = palette,
                      labels = c("Patient (non-NF1)", "Patient (NF1)")) +
  scale_fill_manual(values   = palette,
                    labels   = c("Patient (non-NF1)", "Patient (NF1)")) +
  labs(
    x      = "Age (years)",
    y      = expression(MD~"in Optic Radiations"~"("*x*" 10"^{-4}~mm^2*"/s)"),
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

print(p_md)

ggsave(
  filename = "MD_vs_Age_Patients_LogAge.png",
  plot     = p_md,
  width    = 8, height = 6, dpi = 300
)
cat("\nSaved: MD_vs_Age_Patients_LogAge.png\n")

# ==========================================================
# 6. RESIDUAL PLOT
# ==========================================================

df$fitted   <- predict(fit_log, newdata = df, level = 1)
df$residual <- df$MD - df$fitted

p_resid <- ggplot(df, aes(x = Age, y = residual * 1e5,
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
    x      = "Age (years)",
    y      = expression("Residual"~"("*x*" 10"^{-5}~mm^2*"/s)"),
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
  filename = "Residuals_Patients_LogAge_MD.png",
  plot     = p_resid,
  width    = 8, height = 6, dpi = 300
)
cat("Saved: Residuals_Patients_LogAge_MD.png\n")