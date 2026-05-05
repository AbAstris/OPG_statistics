# ==========================================================
# AD vs AGE — ln(Age) Mixed-Effects Model
#
# Description:
#   This script models the relationship between Axial
#   Diffusivity (AD) in the Optic Radiations and age using a
#   linear mixed-effects model with a log-transformed age term.
#   It is designed for longitudinal DTI data with two input
#   groups: patients (with and without NF1) and healthy controls.
#
# Inputs:
#   - AD_file.csv         : Patient DTI data (see column layout below)
#   - Control_AD_file.csv : Control DTI data (see column layout below)
#   Rename these to match your actual file names.
#
# Outputs:
#   - Console: subject counts and model summary
#   - AD_vs_Age_LogAge.png  : Marginal AD trajectories per group
#   - Residuals_LogAge_AD.png : Conditional residuals per group
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

# Input the patient DTI file, shown here as AD_file.csv with column layout:
# ├── Subject                    (patient ID — numeric or string)
# ├── Weighted Avg AD Combined   (AD measurement — numeric, e.g. 0.0013)
# ├── Age                        (age in years — numeric, 0–20)
# ├── Sex                        (categorical — "F" or "M")
# └── NF1                        (categorical — "Y" or blank/NA)

# Input the control DTI file, shown here as Control_AD_file.csv
# ├── Date                       (control participant ID — despite the name, used as an ID)
# ├── Weighted Avg AD Combined   (AD measurement — numeric, e.g. 0.0013)
# ├── Age                        (age in years — numeric, 0–20)
# └── Sex                        (categorical — "F" or "M")

# ── Patients ───────────────────────────────────────────────
patients <- read_csv("AD_JM_updated.csv", show_col_types = FALSE) %>%
  rename(AD = `Weighted Avg AD Combined`) %>%
  mutate(
    patient_id = as.character(Subject),
    NF1        = ifelse(is.na(NF1) | NF1 == "", "N", NF1),
    Group      = ifelse(NF1 == "Y", "Patient_NF1", "Patient_nonNF1")
  ) %>%
  filter(!is.na(Age), !is.na(AD), Age >= 0, Age <= 20) %>%
  group_by(patient_id) %>%
  filter(n() >= 2) %>%
  ungroup() %>%
  select(patient_id, Age, Sex, Group, AD)

# ── Controls ───────────────────────────────────────────────
controls <- read_csv("../../control_data/Controls_AD_OR_hdbet.csv",
                     show_col_types = FALSE) %>%
  rename(AD = `Weighted Avg AD Combined`) %>%
  mutate(
    patient_id = as.character(Date),
    Group      = "Control"
  ) %>%
  filter(!is.na(Age), !is.na(AD), Age >= 0, Age <= 20) %>%
  select(patient_id, Age, Sex, Group, AD)

# ── Combine ────────────────────────────────────────────────
df <- bind_rows(patients, controls) %>%
  mutate(
    Group  = factor(Group,
                    levels = c("Control", "Patient_nonNF1", "Patient_NF1")),
    Sex    = factor(Sex, levels = c("F", "M")),
    ln_age = log(Age + 1e-3)
  )

cat("Patients (>=2 scans):", n_distinct(patients$patient_id), "\n")
cat("Controls:            ", n_distinct(controls$patient_id), "\n")
cat("Total subjects:      ", n_distinct(df$patient_id), "\n")

# ==========================================================
# 2. FIT MODEL — Log Age
# ==========================================================

ctrl <- lmeControl(maxIter = 5000, msMaxIter = 5000, msMaxEval = 10000)

fit_log <- lme(
  AD ~ Group * ln_age + Sex,
  random  = ~ 1 | patient_id,
  data    = df,
  method  = "ML",
  control = ctrl
)

cat("\n===== Model: Log Age =====\n")
print(summary(fit_log))

# ==========================================================
# 3. COLOUR PALETTE
# ==========================================================

palette <- c(
  "Control"        = "#1f77b4",   # blue
  "Patient_nonNF1" = "#2ca02c",   # green
  "Patient_NF1"    = "#ff7f0e"    # orange
)

# ==========================================================
# 4. Y-AXIS SCALING NOTE
# AD values are multiplied by 1e3 in the plot so that axis
# ticks display as plain numbers e.g. 12.0, 14.0 etc.
# The axis label carries the (x 10^-3mm^2/s) unit explicitly.
# Residuals are multiplied by 1e5 for the same reason.
# ==========================================================

# ==========================================================
# 5. AD vs AGE PLOT
# ==========================================================

age_grid   <- seq(min(df$Age, na.rm = TRUE),
                  max(df$Age, na.rm = TRUE),
                  length.out = 150)
beta       <- fixef(fit_log)
coef_names <- names(beta)
vcov_mat   <- vcov(fit_log)
groups     <- c("Control", "Patient_nonNF1", "Patient_NF1")

pred_list <- lapply(groups, function(grp) {

  nd <- data.frame(
    Age        = age_grid,
    ln_age     = log(age_grid + 1e-3),
    Group      = factor(grp, levels = levels(df$Group)),
    Sex        = factor("F", levels = levels(df$Sex)),
    patient_id = "new_subject"
  )

  nd$pred <- predict(fit_log, newdata = nd, level = 0)

  # ── Delta method SE ──────────────────────────────────────
  se_vec <- numeric(nrow(nd))
  for (k in seq_len(nrow(nd))) {
    cvec <- setNames(rep(0, length(beta)), coef_names)
    cvec["(Intercept)"] <- 1

    gi <- paste0("Group", grp)
    if (gi %in% coef_names) cvec[gi] <- 1

    if ("ln_age" %in% coef_names) cvec["ln_age"] <- nd$ln_age[k]

    gi_ln <- paste0("Group", grp, ":ln_age")
    if (gi_ln %in% coef_names) cvec[gi_ln] <- nd$ln_age[k]

    se_vec[k] <- sqrt(as.numeric(t(cvec) %*% vcov_mat %*% cvec))
  }

  nd$se    <- se_vec
  nd$lower <- nd$pred - 1.96 * nd$se
  nd$upper <- nd$pred + 1.96 * nd$se
  nd
})

pred_df <- bind_rows(pred_list) %>%
  mutate(Group = factor(Group, levels = levels(df$Group)))

p_ad <- ggplot() +
  geom_point(
    data  = df,
    aes(x = Age, y = AD * 1e3, colour = Group),
    alpha = 0.4, size = 1.5
  ) +
  geom_ribbon(
    data  = pred_df,
    aes(x = Age, ymin = lower * 1e3, ymax = upper * 1e3, fill = Group),
    alpha = 0.2
  ) +
  geom_line(
    data      = pred_df,
    aes(x = Age, y = pred * 1e3, colour = Group),
    linewidth = 1.1
  ) +
  scale_colour_manual(values = palette,
                      labels = c("Control",
                                 "Patient (non-NF1)",
                                 "Patient (NF1)")) +
  scale_fill_manual(values = palette,
                    labels = c("Control",
                               "Patient (non-NF1)",
                               "Patient (NF1)")) +
  labs(
    x      = "Age (years)",
    y      = expression(AD~"in Optic Radiations"~"("*x*" 10"^{-3}~mm^2*"/s)"),
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

print(p_ad)

ggsave(
  filename = "AD_vs_Age_LogAge.png",
  plot     = p_ad,
  width    = 8, height = 6, dpi = 300
)
cat("\nSaved: AD_vs_Age_LogAge.png\n")

# ==========================================================
# 6. RESIDUAL PLOT
# ==========================================================

df$fitted   <- predict(fit_log, newdata = df, level = 1)
df$residual <- df$AD - df$fitted

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
                      labels = c("Control",
                                 "Patient (non-NF1)",
                                 "Patient (NF1)")) +
  scale_fill_manual(values = palette,
                    labels = c("Control",
                               "Patient (non-NF1)",
                               "Patient (NF1)")) +
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
  filename = "Residuals_LogAge_AD.png",
  plot     = p_resid,
  width    = 8, height = 6, dpi = 300
)
cat("Saved: Residuals_LogAge_AD.png\n")
