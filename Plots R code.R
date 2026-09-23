# Load required libraries
library(haven)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Read the Stata dataset
df <- read_dta("D:/DRUG EFFICACY/Dataset/Dataset/STH_Longitudinal_Cohort.dta")

# Shared color palette and survey levels
color_palette <- c(
  "Baseline"    = "#0000FF",
  "Follow-up 1" = "#FF0000",
  "Follow-up 2" = "#82E0AA",
  "Follow-up 3" = "#E888E8"
)

# ==============================================================================
# FIGURE 1: PREVALENCE WITH 95% CONFIDENCE INTERVALS (LOWER & UPPER CI) & n
# ==============================================================================

plot_data_prev <- df %>%
  select(Surveycode, Hkinfect, Asinfect, Ttinfect) %>%
  as_factor() %>%
  pivot_longer(
    cols = c(Asinfect, Hkinfect, Ttinfect),
    names_to = "Parasite",
    values_to = "Infected"
  ) %>%
  filter(!is.na(Infected), !is.na(Surveycode)) %>%
  group_by(Parasite, Surveycode) %>%
  summarise(
    Total = n(),
    Positive = sum(Infected == 1 | Infected == "Positive" | Infected == "1" | Infected == "positive", na.rm = TRUE),
    Prevalence = (Positive / Total) * 100,
    # Standard Error for Binomial Proportion
    SE = sqrt((Prevalence / 100) * (1 - (Prevalence / 100)) / Total) * 100,
    # 95% Lower and Upper Confidence Intervals
    Lower_CI = pmax(0, Prevalence - (1.96 * SE)),
    Upper_CI = pmin(100, Prevalence + (1.96 * SE)),
    .groups = "drop"
  ) %>%
  mutate(
    Parasite = case_when(
      Parasite == "Asinfect" ~ "A. lumbricoides",
      Parasite == "Hkinfect" ~ "Hookworm",
      Parasite == "Ttinfect" ~ "T. trichiura"
    ),
    Surveycode = trimws(as.character(Surveycode))
  )

ggplot(plot_data_prev, aes(x = Surveycode, y = Prevalence, fill = Surveycode)) +
  geom_col(width = 0.7) +
  # Lower and Upper 95% CI Error Bars
  geom_errorbar(
    aes(ymin = Lower_CI, ymax = Upper_CI),
    width = 0.2,
    linewidth = 0.5
  ) +
  # Positive sample size label above the upper CI limit
  geom_text(
    aes(y = Upper_CI, label = paste0("n=", Positive)),
    vjust = -0.5,
    size = 3
  ) +
  facet_wrap(~ Parasite, scales = "free_y") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.2))
  ) +
  scale_fill_manual(values = color_palette, name = NULL) +
  labs(
    x = NULL,
    y = "Prevalence (%)"
  ) +
  theme_classic() +
  theme(
    strip.text = element_text(face = "bold.italic", size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 10),
    legend.position = "bottom"
  )

ggsave(
  filename = "Fig1-Prevalence_Scaled.tiff",
  width = 9,
  height = 5,
  dpi = 300,
  compression = "lzw"
)

# ==============================================================================
# FIGURE 2: MEAN INTENSITY WITH 95% CONFIDENCE INTERVALS (LOWER & UPPER CI) & n
# ==============================================================================

plot_data_epg <- df %>%
  select(Surveycode, hkepg, asepg, ttepg) %>%
  as_factor() %>%
  pivot_longer(
    cols = c(asepg, hkepg, ttepg),
    names_to = "Parasite",
    values_to = "EPG"
  ) %>%
  filter(!is.na(EPG), !is.na(Surveycode)) %>%
  group_by(Parasite, Surveycode) %>%
  summarise(
    Mean_EPG = mean(EPG, na.rm = TRUE),
    SE = sd(EPG, na.rm = TRUE) / sqrt(n()),
    # 95% Lower and Upper Confidence Intervals
    Lower_CI = pmax(0, Mean_EPG - (1.96 * SE)),
    Upper_CI = Mean_EPG + (1.96 * SE),
    Positive = sum(EPG > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Parasite = case_when(
      Parasite == "asepg" ~ "A. lumbricoides",
      Parasite == "hkepg" ~ "Hookworm",
      Parasite == "ttepg" ~ "T. trichiura"
    ),
    Surveycode = trimws(as.character(Surveycode))
  )

ggplot(plot_data_epg, aes(x = Surveycode, y = Mean_EPG, fill = Surveycode)) +
  geom_col(width = 0.7) +
  # Lower and Upper 95% CI Error Bars
  geom_errorbar(
    aes(ymin = Lower_CI, ymax = Upper_CI),
    width = 0.2,
    linewidth = 0.5
  ) +
  # Positive sample size label above the upper CI limit
  geom_text(
    aes(y = Upper_CI, label = paste0("n=", Positive)),
    vjust = -0.5,
    size = 3
  ) +
  facet_wrap(~ Parasite, scales = "free_y") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.25))
  ) +
  scale_fill_manual(values = color_palette, name = NULL) +
  labs(
    x = NULL,
    y = "Mean Intensity (EPG)"
  ) +
  theme_classic() +
  theme(
    strip.text = element_text(face = "bold.italic", size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 10),
    legend.position = "bottom"
  )

ggsave(
  filename = "Fig2-Mean_Intensity_Scaled.tiff",
  width = 9,
  height = 5,
  dpi = 300,
  compression = "lzw"
)

# ==============================================================================
# FIGURE 3: INDIVIDUAL EGG COUNT TRAJECTORIES (BASELINE TO FOLLOW-UPS)
# ==============================================================================

long_df <- df %>%
  select(unique_id, Surveycode, hkepg, asepg, ttepg) %>%
  as_factor() %>%
  pivot_longer(
    cols = c(asepg, hkepg, ttepg),
    names_to = "Parasite",
    values_to = "EPG"
  ) %>%
  filter(!is.na(EPG), !is.na(Surveycode)) %>%
  mutate(
    Parasite = case_when(
      Parasite == "asepg" ~ "A) A. lumbricoides",
      Parasite == "hkepg" ~ "B) Hookworm",
      Parasite == "ttepg" ~ "C) T. trichiura"
    ),
    Surveycode = trimws(as.character(Surveycode))
  )

# Extract baseline EPG per subject
baseline_df <- long_df %>%
  filter(tolower(Surveycode) == "baseline") %>%
  select(unique_id, Parasite, Baseline_EPG = EPG)

# Pair baseline directly with each follow-up observation
followup_df <- long_df %>%
  filter(tolower(Surveycode) %in% c("follow-up 1", "follow-up 2", "follow-up 3")) %>%
  inner_join(baseline_df, by = c("unique_id", "Parasite")) %>%
  mutate(
    Comparison = case_when(
      tolower(Surveycode) == "follow-up 1" ~ "baseline to follow-up 1",
      tolower(Surveycode) == "follow-up 2" ~ "baseline to follow-up 2",
      tolower(Surveycode) == "follow-up 3" ~ "baseline to follow-up 3"
    )
  )

segment_data <- followup_df %>%
  pivot_longer(
    cols = c(Baseline_EPG, EPG),
    names_to = "TimePoint",
    values_to = "Egg_Count"
  ) %>%
  mutate(
    X_Label = ifelse(TimePoint == "Baseline_EPG", "baseline", tolower(Surveycode)),
    X_Label = factor(X_Label, levels = c("baseline", "follow-up 1", "follow-up 2", "follow-up 3")),
    Segment_ID = paste(unique_id, Comparison, sep = "_")
  )

ggplot(segment_data, aes(x = X_Label, y = Egg_Count)) +
  geom_line(
    aes(group = Segment_ID, color = Comparison),
    alpha = 0.7,
    linewidth = 0.5
  ) +
  geom_point(color = "black", size = 1.2, alpha = 0.8) +
  facet_wrap(~ Parasite, scales = "free_y", ncol = 2) +
  scale_color_manual(
    values = c(
      "baseline to follow-up 1" = "black",
      "baseline to follow-up 2" = "#FF0000",
      "baseline to follow-up 3" = "#82E0AA"
    ),
    name = NULL
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    x = NULL,
    y = "Individual egg count (epg)"
  ) +
  theme(
    panel.background = element_rect(fill = "#EBEBEB", color = NA),
    panel.grid.major = element_line(color = "white", linewidth = 0.6),
    panel.grid.minor = element_line(color = "white", linewidth = 0.3),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold.italic", size = 12, hjust = 0),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    legend.position = c(0.6, 0.25),
    legend.background = element_blank(),
    legend.key = element_blank()
  )

ggsave(
  filename = "Fig3-Individual_egg_counts.tiff",
  width = 10,
  height = 6,
  dpi = 300,
  compression = "lzw"
)
