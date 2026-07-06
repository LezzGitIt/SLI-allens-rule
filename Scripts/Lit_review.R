## Literature review of methods used to test Allen's rule since 2021.
## Loads and cleans Methods_ecogeographic_rules_5.7.26.xlsx, tabulates the
## analytical approaches used across studies, examines body-size-control
## variable choice, and produces Figures/lit_review.png. Exports the values
## cited inline by Scripts/qmd/Allens_methods_sim.qmd to
## Derived/Rds/lit_review_results.rds.

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(readxl)
library(janitor)

# Lit review formatting ---------------------------------------------------

Lit_review <- read_excel("../Lit_review/Methods_ecogeographic_rules_7.6.26.xlsx",
                         sheet = "Allens") %>%
  clean_names()
Unclear_exclude <- Lit_review %>% filter(include_exclude == "Exclude")
Lr <- Lit_review %>%
  anti_join(Unclear_exclude) %>%
  dplyr::select(
    authors, year, regression_used, class, body_size_control, starts_with("approach")
  ) %>%
  mutate(study = factor(row_number()))

Lr_pivot <- Lr %>%
  pivot_longer(cols = c(starts_with("approach"), -approach_unclear),
               names_to = "approach",
               values_to = "primary_v_secondary") %>%
  mutate(approach = str_remove(approach, "approach_")) %>%
  filter(!is.na(primary_v_secondary) & primary_v_secondary != "M") %>%
  mutate(approach_unclear = ifelse(is.na(approach_unclear), "No", "Yes")) %>%
  filter(approach_unclear == "No")

Lr_per <- Lr_pivot %>%
  summarize(num_studies = n(), .by = c(approach, class)) %>%
  mutate(Tot_studies = sum(num_studies),
         Prop = num_studies / Tot_studies,
         .by = class) %>%
  mutate(Per = Prop * 100,
         approach_label = str_to_sentence(approach),
         approach_label = str_replace_all(approach_label, "_", "\n"),
         approach_label = as.factor(approach_label)) %>%
  mutate(Per = round(Per, 0), Prop = round(Prop, 2)) %>%
  arrange(desc(Per))

Lr_stack <- Lr_pivot %>%
  arrange(class, study) %>%
  mutate(ypos = row_number(), .by = c(approach, class)) %>%
  left_join(Lr_per[, c("approach", "approach_label", "class", "num_studies")])

Lr_per2 <- Lr_per %>%
  group_by(class) %>%
  mutate(approach_facet = fct_reorder(paste(class, approach_label, sep = "___"),
                                      num_studies, .desc = TRUE),
         Per = paste0(Per, "%")) %>%
  ungroup()

Lr_stack2 <- Lr_stack %>%
  left_join(Lr_per2 %>% distinct(class, approach_label, approach_facet),
            by = c("class", "approach_label"))

# Body size control  --------------------------------------------------------
# 46 approaches across 36 studies
Lr %>% filter(approach_raw_values == "P") %>%
  dplyr::select(-regression_used)

## De-duplicate to one row per study before tabulating (body_size_control doesn't vary by approach)
# A study with multiple qualifying approaches would otherwise have its control variables counted once per approach instead of once; the literal string "NA" (typed into some cells) and true NA both mean "no body-size control reported".
Body_size_control_studies <- Lr_pivot %>%
  filter(approach != "raw_values") %>%
  distinct(study, body_size_control) %>%
  mutate(body_size_control = na_if(body_size_control, "NA"))

## 1) How many studies use mass (in any combination, e.g. "Mass, Tarsus") as a body-size control?
Mass_control_tbl <- Body_size_control_studies %>%
  mutate(Mass_control = str_detect(body_size_control, regex("mass", ignore_case = TRUE))) %>%
  tabyl(Mass_control)

Body_size_control_studies 

Lr_pivot %>%
  filter(approach != "raw_values") %>% 
  filter(is.na(body_size_control))

## 2) Overall trait distribution pooled across studies
# Format by splitting multi-control entries by comma, standardize whitespace/capitalization, then classify each individual mention into Mass / Linear / Other.
classify_control <- function(x) {
  x_low <- str_to_lower(x)
  case_when(
    str_detect(x_low, "mass") ~ "Mass",
    str_detect(x_low, "length") | x_low %in% c("wing chord", "tarsus", "keel") ~ "Linear measurement",
    TRUE ~ "Other"
  )
}

Body_size_control_traits <- Body_size_control_studies %>%
  filter(!is.na(body_size_control)) %>%
  separate_longer_delim(body_size_control, delim = ",") %>%
  mutate(body_size_control = str_squish(body_size_control),
         body_size_control = str_to_sentence(body_size_control),
         control_group = classify_control(body_size_control)) %>% 
  filter(body_size_control != "Mass^2")
Trait_group_tbl <- Body_size_control_traits %>% tabyl(control_group)
Trait_group_tbl 

# Figure ------------------------------------------------------------------
# Produce lit review figure
fig_lit_review <- ggplot() +
  geom_col(data = Lr_per2,
           aes(x = approach_facet, y = num_studies),
           fill = "grey85", color = "black") +
  geom_point(data = Lr_stack2,
             aes(x = approach_facet, y = ypos, color = study),
             alpha = 0.7, size = 3) +
  geom_line(data = Lr_stack2,
            aes(x = approach_facet, y = ypos, group = study, color = study),
            alpha = 0.75) +
  geom_text(data = Lr_per2,
            aes(x = approach_facet, y = num_studies, label = Per),
            vjust = -0.7, size = 4.5) +
  facet_wrap(~class, scales = "free_x",
             labeller = labeller(.default = label_value)) +
  scale_x_discrete(labels = function(x) sub(".*___", "", x)) +
  theme(axis.text.x = element_text(size = 9, vjust = 0.58, angle = 55),
        legend.position = "none") +
  labs(x = NULL, y = "Number of studies")
fig_lit_review

ggsave("Figures/lit_review.png", fig_lit_review, bg = "white",
       width = 9, height = 6.5, units = "in", dpi = 300)

# Export ---------------------------------------------------------------------
dir.create("Derived/Rds", showWarnings = FALSE)
saveRDS(
  list(
    Unclear_exclude = Unclear_exclude,
    Lit_review = Lit_review,
    Lr = Lr,
    Lr_per = Lr_per
  ),
  "Derived/Rds/lit_review_results.rds"
)
