# ============================================================
# Care-seeking for sexually transmitted infections and
# STI-related symptoms among Nigerian women
# NDHS 2023-24
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(haven)
library(dplyr)
library(survey)

options(survey.lonely.psu = "adjust")


# ============================================================
# 2. LOAD MASTER DATASET
# ============================================================

# The NDHS dataset is not distributed with this analysis code.
# Researchers should obtain the authorised dataset through the
# DHS data-access process and specify its local file path below.

data_path <- "PATH_TO_NDHS_2024_master.rds"

ndhs <- readRDS(data_path)

dim(ndhs)

# ============================================================
# 3. DEFINE STUDY 2 ANALYTICAL POPULATION
# ============================================================

# Women reporting at least one of:
# V763A = STI in last 12 months
# V763B = genital sore/ulcer
# V763C = genital discharge
# V763D = genital itching

study2_eligible_v770 <- with(
  ndhs,
  V763A == 1 |
    V763B == 1 |
    V763C == 1 |
    V763D == 1
)


# Check eligibility
table(
  study2_eligible_v770,
  useNA = "ifany"
)


# Check that all eligible women have valid V770
table(
  study2_eligible_v770,
  is.na(ndhs$V770),
  useNA = "ifany"
)


# Create analytical dataset
study2_data <- ndhs %>%
  filter(study2_eligible_v770)


# Check analytical population
dim(study2_data)


# Check outcome
table(
  study2_data$V770,
  useNA = "ifany"
)

prop.table(
  table(study2_data$V770)
) * 100


# ============================================================
# 4. REPORTED STI-RELATED CONDITIONS/SYMPTOMS
# ============================================================

study2_data$n_reported_conditions <-
  (study2_data$V763A == 1) +
  (study2_data$V763B == 1) +
  (study2_data$V763C == 1) +
  (study2_data$V763D == 1)


table(
  study2_data$n_reported_conditions,
  useNA = "ifany"
)


# Create binary outcome for regression
study2_data$care_sought <- ifelse(
  study2_data$V770 == 1,
  1,
  0
)


# Factor version for categorical burden analysis
study2_data$burden_cat <- factor(
  study2_data$n_reported_conditions,
  levels = 1:4,
  labels = c("1", "2", "3", "4")
)


# ============================================================
# 5. SURVEY DESIGN
# ============================================================

study2_design <- svydesign(
  ids = ~V021,
  strata = ~V022,
  weights = ~I(V005 / 1000000),
  data = study2_data,
  nest = TRUE
)


study2_design


# Overall weighted care-seeking
svymean(
  ~I(V770 == 1),
  study2_design
)


# ============================================================
# 6. TABLE 1
# CHARACTERISTICS AND CARE-SEEKING
# ============================================================


# ---- Education ----

table1_education <- svyby(
  ~I(V770 == 1),
  ~V106,
  study2_design,
  svymean,
  vartype = "ci",
  na.rm = TRUE
)

table1_education


# ---- Wealth ----

table1_wealth <- svyby(
  ~I(V770 == 1),
  ~V190,
  study2_design,
  svymean,
  vartype = "ci",
  na.rm = TRUE
)

table1_wealth


# ---- Residence ----

table1_residence <- svyby(
  ~I(V770 == 1),
  ~V025,
  study2_design,
  svymean,
  vartype = "ci",
  na.rm = TRUE
)

table1_residence


# ---- Region ----

table1_region <- svyby(
  ~I(V770 == 1),
  ~V024,
  study2_design,
  svymean,
  vartype = "ci",
  na.rm = TRUE
)

table1_region


# ---- Age group ----

table1_age <- svyby(
  ~I(V770 == 1),
  ~V013,
  study2_design,
  svymean,
  vartype = "ci",
  na.rm = TRUE
)

table1_age


# ---- Marital/union status ----

table1_marital <- svyby(
  ~I(V770 == 1),
  ~V501,
  study2_design,
  svymean,
  vartype = "ci",
  na.rm = TRUE
)

table1_marital


# ---- Rao-Scott P-values ----

p_education <- svychisq(
  ~V106 + V770,
  study2_design,
  statistic = "F"
)

p_wealth <- svychisq(
  ~V190 + V770,
  study2_design,
  statistic = "F"
)

p_residence <- svychisq(
  ~V025 + V770,
  study2_design,
  statistic = "F"
)

p_region <- svychisq(
  ~V024 + V770,
  study2_design,
  statistic = "F"
)

p_age <- svychisq(
  ~V013 + V770,
  study2_design,
  statistic = "F"
)

p_marital <- svychisq(
  ~V501 + V770,
  study2_design,
  statistic = "F"
)


p_education
p_wealth
p_residence
p_region
p_age
p_marital


# ============================================================
# 7. CRUDE EDUCATION MODEL
# ============================================================

model_edu_crude <- svyglm(
  V770 ~ factor(V106),
  design = study2_design,
  family = quasipoisson(link = "log")
)

summary(model_edu_crude)

exp(
  cbind(
    PR = coef(model_edu_crude),
    confint(model_edu_crude)
  )
)


# ============================================================
# 8. CRUDE WEALTH MODEL
# ============================================================

model_wealth_crude <- svyglm(
  V770 ~ factor(V190),
  design = study2_design,
  family = quasipoisson(link = "log")
)

summary(model_wealth_crude)

exp(
  cbind(
    PR = coef(model_wealth_crude),
    confint(model_wealth_crude)
  )
)


# ============================================================
# 9. FINAL MULTIVARIABLE MODEL
# ============================================================

model_main <- svyglm(
  V770 ~ factor(V106) +
    factor(V190) +
    factor(V013) +
    factor(V025) +
    factor(V024) +
    factor(V501),
  design = study2_design,
  family = quasipoisson(link = "log")
)


summary(model_main)


# Adjusted prevalence ratios
exp(
  cbind(
    aPR = coef(model_main),
    confint(model_main)
  )
)


# ============================================================
# 10. WEALTH TREND ANALYSIS
# ============================================================

model_wealth_trend <- svyglm(
  V770 ~ V190 +
    factor(V106) +
    factor(V013) +
    factor(V025) +
    factor(V024) +
    factor(V501),
  design = study2_design,
  family = quasipoisson(link = "log")
)


summary(model_wealth_trend)


exp(
  cbind(
    aPR = coef(model_wealth_trend),
    confint(model_wealth_trend)
  )
)


# ============================================================
# 11. CONDITION/SYMPTOM BURDEN
# ============================================================


# Weighted care-seeking prevalence by number of conditions

burden_prevalence <- svyby(
  ~I(V770 == 1),
  ~factor(n_reported_conditions),
  study2_design,
  svymean,
  vartype = "ci",
  na.rm = TRUE
)

burden_prevalence


# Rao-Scott test

svychisq(
  ~care_sought + burden_cat,
  study2_design,
  statistic = "F"
)


# ============================================================
# 12. CRUDE BURDEN MODEL
# ============================================================

model_burden_crude <- svyglm(
  care_sought ~ burden_cat,
  design = study2_design,
  family = quasipoisson(link = "log")
)


summary(model_burden_crude)


exp(
  cbind(
    PR = coef(model_burden_crude),
    confint(model_burden_crude)
  )
)


# ============================================================
# 13. ADJUSTED BURDEN MODEL
# ============================================================

model_burden_adj <- svyglm(
  care_sought ~ burden_cat +
    factor(V106) +
    factor(V190) +
    factor(V013) +
    factor(V025) +
    factor(V024) +
    factor(V501),
  design = study2_design,
  family = quasipoisson(link = "log")
)


summary(model_burden_adj)


exp(
  cbind(
    aPR = coef(model_burden_adj),
    confint(model_burden_adj)
  )
)


# ============================================================
# 14. BURDEN TREND MODEL
# ============================================================

model_burden_trend <- svyglm(
  care_sought ~ n_reported_conditions +
    factor(V106) +
    factor(V190) +
    factor(V013) +
    factor(V025) +
    factor(V024) +
    factor(V501),
  design = study2_design,
  family = quasipoisson(link = "log")
)


summary(model_burden_trend)


exp(
  cbind(
    aPR = coef(model_burden_trend),
    confint(model_burden_trend)
  )
)


# ============================================================
# 15. SOURCES OF ADVICE/TREATMENT
# ============================================================

source_vars <- grep(
  "^V770[A-Z]",
  names(study2_data),
  value = TRUE
)


source_summary <- data.frame(
  Variable = source_vars,
  
  Source = sapply(
    study2_data[source_vars],
    function(x) attr(x, "label")
  ),
  
  n = sapply(
    study2_data[source_vars],
    function(x) sum(x == 1, na.rm = TRUE)
  )
)


# Unweighted percentage among those who sought care
source_summary$Percent_unweighted_among_seekers <-
  source_summary$n /
  sum(study2_data$V770 == 1) *
  100


source_summary


# ============================================================
# 16. SURVEY-WEIGHTED SOURCE-OF-CARE ESTIMATES
# ============================================================

seekers_design <- subset(
  study2_design,
  V770 == 1
)


source_weighted <- lapply(
  source_vars,
  function(v) {
    
    f <- as.formula(
      paste0("~", v)
    )
    
    est <- svymean(
      f,
      seekers_design,
      na.rm = TRUE
    )
    
    ci <- confint(est)
    
    data.frame(
      Variable = v,
      n = sum(
        study2_data[study2_data$V770 == 1, v] == 1,
        na.rm = TRUE
      ),
      Weighted_percent =
        as.numeric(coef(est)) * 100,
      Lower_CI =
        as.numeric(ci[, 1]) * 100,
      Upper_CI =
        as.numeric(ci[, 2]) * 100
    )
  }
)


source_weighted <- bind_rows(
  source_weighted
)


source_weighted


# ============================================================
# 17. END OF CORE STUDY 2 ANALYSIS
# ============================================================