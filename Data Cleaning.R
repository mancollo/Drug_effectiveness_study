# ------------------------------------------------------------------------------
# 1. SETUP & LIBRARIES
# ------------------------------------------------------------------------------
library(haven)
library(dplyr)

data_dir <- "D:/DRUG EFFICACY/Dataset/"

# ------------------------------------------------------------------------------
# 2. IMPORT DATASETS
# ------------------------------------------------------------------------------
overall_lab <- read_dta(file.path(data_dir, "Overall_Lab_data.dta"))
child_q     <- read_dta(file.path(data_dir, "Child_questionnaire_April 2021.dta"))

# ------------------------------------------------------------------------------
# 3. PREPARE DEMOGRAPHICS & EXECUTE JOIN WITH DROPS
# ------------------------------------------------------------------------------
# Clean IDs and extract unique demographic profiles
child_demo <- child_q %>%
  mutate(unique_id = trimws(as.character(unique_id))) %>%
  select(unique_id, age, sex) %>%
  distinct(unique_id, .keep_all = TRUE)

# Merge and immediately filter out missing demographics
merged_data <- overall_lab %>%
  mutate(unique_id = trimws(as.character(unique_id))) %>%
  left_join(child_demo, by = "unique_id") %>%
  filter(!is.na(age) & !is.na(sex)) %>% # Drops any row where age or sex is missing
  filter (!is.na(sthinfect) & !is.na(hkinfect) & !is.na(asinfect) & !is.na(ttinfect))


# ------------------------------------------------------------------------------
# 3. DATA CLEANING & AGE CATEGORIZATION
# ------------------------------------------------------------------------------
clean_data <- merged_data %>%
  # Clean text & factors
  mutate(
    unique_id   = trimws(as.character(unique_id)),
    county_code = trimws(as_factor(county_code)),
    school_code = trimws(as_factor(school_code)),
    survey      = trimws(as_factor(survey)),
    sex         = trimws(as_factor(sex))
  ) %>%
  # Filter valid age range and demographic completeness
  filter(!is.na(sex) & !is.na(age)) %>%
  filter(age >= 0 & age <= 18) %>%
  # Create ordered age categories
  mutate(
    age_cat = case_when(
      age <= 5              ~ "<= 5 years", # Fixed label to match factor level below
      age <= 10             ~ "6-10 years",
      age >= 11 & age <= 14 ~ "11-14 years",
      age >= 15 & age <= 18 ~ "15-18 years",
      TRUE                  ~ NA_character_
    ),
    age_cat = factor(age_cat, levels = c("<= 5 years", "6-10 years", "11-14 years", "15-18 years")),
    survey  = factor(survey, levels = c("Baseline", "Follow-up 1", "Follow-up 2", "Follow-up 3"))
  ) %>%
  filter(!is.na(survey) & !is.na(age_cat)) %>%
  # Filter complete infection cases
  filter(!is.na(sthinfect) & !is.na(sthepg)) %>%
  distinct(unique_id, survey, .keep_all = TRUE)

# ------------------------------------------------------------------------------
# 4. FILTER COHORT: KEEP ALL BASELINE + ONLY BASELINE-POSITIVES FOR FOLLOW-UPS
# ------------------------------------------------------------------------------
# Extract unique_ids of children who tested positive for STH (sthinfect == 1) at Baseline
sth_baseline_pos_ids <- clean_data %>%
  filter(survey == "Baseline" & sthinfect == 1) %>%
  pull(unique_id) %>%
  unique()

# Keep ALL Baseline cases (positive & negative) OR Baseline-positive cases in Follow-ups
sth_positive_cohort <- clean_data %>%
  filter(survey == "Baseline" | unique_id %in% sth_baseline_pos_ids)

# ------------------------------------------------------------------------------
# 6. SAVE CLEANED DATASET
# ------------------------------------------------------------------------------
saveRDS(sth_positive_cohort, file.path(data_dir, "STH_Baseline_Cohort.rds"))
write_dta(sth_positive_cohort, file.path(data_dir, "STH_Baseline_Cohort.dta"))

cat("\nCleaned dataset saved successfully to directory!\n")

# ------------------------------------------------------------------------------
# 4. FILTER COHORT: RETAIN ONLY BASELINE-POSITIVE INDIVIDUALS & FOLLOW-UPS
# ------------------------------------------------------------------------------
# Extract unique_ids of children who tested positive for STH (sthinfect == 1) at Baseline
sth_baseline_pos_ids <- clean_data %>%
  filter(survey == "Baseline" & sthinfect == 1) %>%
  pull(unique_id) %>%
  unique()

# Filter dataset to include ONLY baseline-positive individuals across all their survey rounds
sth_longitudinal_cohort <- clean_data %>%
  filter(unique_id %in% sth_baseline_pos_ids)

# ------------------------------------------------------------------------------
# 5. SAVE CLEANED DATASET
# ------------------------------------------------------------------------------
saveRDS(sth_longitudinal_cohort, file.path(data_dir, "STH_Longitudinal_Cohort.rds"))
write_dta(sth_longitudinal_cohort, file.path(data_dir, "STH_Longitudinal_Cohort.dta"))

cat("\nSTH Longitudinal Cohort dataset saved successfully to directory!\n")

