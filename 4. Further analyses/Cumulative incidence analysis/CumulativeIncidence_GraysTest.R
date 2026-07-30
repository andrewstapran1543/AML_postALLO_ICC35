renv::install("survminer")
renv::install(c("survminer", "ggplot2"))
library("survminer")
require("survival")
library("readxl")

library(readxl)
library(dplyr)
library(stringr)
library(glue)
library(dplyr)
renv::install("bioc::cmprsk")
library(cmprsk)

cutoff <- 2000  # set your cutoff
landmarked <- FALSE

# ---- reading the dataset -----
d35_input_full <- read_excel(
  "/g/scb2/zaugg/stapran/24_09_FACS_Project_Iter2/DATA/patient_data_excels/d35/FinalPopulationPercentages_ClinicalInfo_3000threshold_2026_06_28.xlsx"
)

d35_input_full_copy <- d35_input_full %>%
  filter(
    OC_d35_PassCrit_FCAnalysis == "Yes" |
      VC_d35_PassCrit_FCAnalysis == "Yes"
  )

# Percent features
PERCENT_FEATURES <- names(d35_input_full_copy)[str_detect(names(d35_input_full_copy), " of ")]
myeloid_endings <- c(" of DC/HLA-DRpos_APC", " of Monocytes", " of Myeloid")
myeloid_endings_and_CD45pos <- c(myeloid_endings, " of CD45pos")

myeloid_child <- PERCENT_FEATURES[
  str_ends(PERCENT_FEATURES, paste0(myeloid_endings, collapse = "|"))
] %>%
  str_split_fixed(" of ", 2) %>%
  .[, 1] %>%
  unique()

myeloid_params_exclude <- PERCENT_FEATURES[
  str_starts(PERCENT_FEATURES, paste0(myeloid_child, collapse = "|")) &
    str_ends(PERCENT_FEATURES, paste0(myeloid_endings_and_CD45pos, collapse = "|"))
]

PERCENT_FEATURES <- setdiff(PERCENT_FEATURES, myeloid_params_exclude)

# Rename LowLevel -> High Resolution
old_names <- PERCENT_FEATURES[str_detect(PERCENT_FEATURES, "LowLevel")]
new_names <- str_replace(old_names,"LowLevel ","High Resolution ")
PERCENT_FEATURES_renamer <- setNames(new_names, old_names)

d35_input_full_copy <- d35_input_full_copy %>%
  rename_with(
    ~ str_replace(.x, "LowLevel ", "High Resolution "),
    all_of(names(PERCENT_FEATURES_renamer))
  )
PERCENT_FEATURES <- str_replace(PERCENT_FEATURES, "LowLevel ", "High Resolution ")

# Ratio features
RATIO_FEATURES <- names(d35_input_full_copy)[str_detect(names(d35_input_full_copy), "-to-")]

# Optional arcsin sqrt transform
# d35_input_full_copy <- d35_input_full_copy %>%
#   mutate(across(
#     all_of(PERCENT_FEATURES),
#     ~ asin(sqrt(pmin(pmax(.x / 100, 0), 1)))
#   ))

OUTCOME_COLS <- list(
  RFS = list(time = "Time_to_Relapse_from_TPL", event = "Relapse_Present_1_NotObserved_0"),
  OS = list(time = "Time_to_Death_from_TPL", event = "Death_Present_1_NotObserved_0"),
  aGVHD = list(time = "Time_to_aGVHD", event = "aGVHD_grading_by_Katja"),
  cGVHD = list(time = "Time_to_cGVHD", event = "cGVHD_grading_by_Katja")
)

d35_input_full_copy <- d35_input_full_copy %>% mutate(FC_Sampling_Days_After_Allo_d35 = as.integer(round(as.numeric(FC_Sampling_Days_After_Allo_d35), 0)))
d35_input_full_copy$Cutoff_Plus_SamplingDay <- cutoff + d35_input_full_copy$FC_Sampling_Days_After_Allo_d35

# version where we do 2000 days post-TPL
if (landmarked == FALSE) {
  apply_cutoff <- function(df, outcome, cutoff) {
    time_col <- OUTCOME_COLS[[outcome]]$time
    event_col <- OUTCOME_COLS[[outcome]]$event
    
    df %>%
      mutate(
        "{time_col}" := as.numeric(.data[[time_col]]),
        "{event_col}" := as.numeric(.data[[event_col]])
      ) %>%
      mutate(
        "{event_col}" := if_else(.data[[time_col]] >= cutoff, 0, .data[[event_col]]),
        "{time_col}"  := if_else(.data[[time_col]] >= cutoff, cutoff, .data[[time_col]])
      )
  }
  for (outcome in names(OUTCOME_COLS)) {
    d35_input_full_copy <- apply_cutoff(d35_input_full_copy, outcome, cutoff)
  }
} else {
  # version where do 2000 days post-sampling
  apply_cutoff <- function(df, outcome, cutoff) {
    time_col <- OUTCOME_COLS[[outcome]]$time
    event_col <- OUTCOME_COLS[[outcome]]$event
    
    df %>%
      mutate(
        "{time_col}" := as.numeric(.data[[time_col]]),
        "{event_col}" := as.numeric(.data[[event_col]])
      ) %>%
      mutate(
        "{event_col}" := if_else(.data[[time_col]] >= .data[["Cutoff_Plus_SamplingDay"]], 0, .data[[event_col]]),
        "{time_col}"  := if_else(.data[[time_col]] >= .data[["Cutoff_Plus_SamplingDay"]], .data[["Cutoff_Plus_SamplingDay"]], .data[[time_col]])
      )
  }
  for (outcome in names(OUTCOME_COLS)) {
    d35_input_full_copy <- apply_cutoff(d35_input_full_copy, outcome, cutoff)
  } 
}



d35_input_full_copy <- d35_input_full_copy %>%
  mutate(
    across(
      c(
        Relapse_Present_1_NotObserved_0,
        aGVHD_grading_by_Katja,
        cGVHD_grading_by_Katja,
        FC_Sampling_Days_After_Allo_d35,
        Time_to_aGVHD,
        Time_to_cGVHD,
        Time_to_Relapse_from_TPL,
        Time_to_Death_from_TPL
      ),
      as.numeric
    )
  ) %>%
  mutate(
    Relapse_Present_1_NotObserved_0 = as.integer(Relapse_Present_1_NotObserved_0),
    aGVHD_grading_by_Katja = as.integer(aGVHD_grading_by_Katja),
    cGVHD_grading_by_Katja = as.integer(cGVHD_grading_by_Katja),
    
    # FC_Sampling_Days_After_Allo_d35 =
    #   as.integer(round(FC_Sampling_Days_After_Allo_d35, 0)),
    
    MAX_FOLLOW_UP = pmax(
      Time_to_aGVHD,
      Time_to_cGVHD,
      Time_to_Relapse_from_TPL,
      Time_to_Death_from_TPL,
      na.rm = TRUE
    )
  ) 

if (landmarked == TRUE) {
  d35_input_full_copy <- d35_input_full_copy %>%
    mutate(
      Time_to_Death_from_TPL =
        Time_to_Death_from_TPL - FC_Sampling_Days_After_Allo_d35,
      
      Time_to_Relapse_from_TPL =
        Time_to_Relapse_from_TPL - FC_Sampling_Days_After_Allo_d35,
      
      Time_to_aGVHD =
        Time_to_aGVHD - FC_Sampling_Days_After_Allo_d35,
      
      Time_to_cGVHD =
        Time_to_cGVHD - FC_Sampling_Days_After_Allo_d35,
      
      MAX_FOLLOW_UP =
        MAX_FOLLOW_UP - FC_Sampling_Days_After_Allo_d35
    )
  x_axis_label <- "Days post alloHCT sampling"
} else {
  x_axis_label <- "Days post alloHCT"
}
# -----

# ----- setting up medians & calculating the scores -----
d35_input_full_copy_С1 <- d35_input_full_copy %>% subset(Cohort == 'OC')
d35_input_full_copy_С2 <- d35_input_full_copy %>% subset(Cohort == 'VC')
d35_input_full_copy_С1$

CD8_parameter_General = "CD8 TEM adv1 GPR56+ CD57+ of CD8 TEM"
CD8_parameter_HighRes = "High Resolution CD8 CCR7- CD45RAmid CD27dim CD57+ GPR56+ of CD45pos"

CD56dimNK_parameter_General = "CD56dim NK CD57+ GPR56+ of Lymphocytes"
CD56dimNK_parameter_HighRes = "High Resolution CD56dim NK CD45RA+ CD57+ GPR56+ of Lymphocytes"

CD4TCM_parameter = "CD4 CCR7mid CD45RA- of T"
CD4TCMNVratio_parameter = "CD4_Naive+CD4_TCM-to-CD4_rest"

medians <- list(
  cohort1 = list(
    CD8_parameter_General = median(as.numeric(d35_input_full_copy_С1[[CD8_parameter_General]])),
    CD8_parameter_HighRes = median(as.numeric(d35_input_full_copy_С1[[CD8_parameter_HighRes]])),
    CD56dimNK_parameter_General = median(as.numeric(d35_input_full_copy_С1[[CD56dimNK_parameter_General]])),
    CD56dimNK_parameter_HighRes = median(as.numeric(d35_input_full_copy_С1[[CD56dimNK_parameter_HighRes]])),
    CD4TCM_parameter = median(as.numeric(d35_input_full_copy_С1[[CD4TCM_parameter]])),
    CD4TCMNVratio_parameter = median(as.numeric(d35_input_full_copy_С1[[CD4TCMNVratio_parameter]]))
  ),
  
  cohort2 = list(
    CD8_parameter_General = median(as.numeric(d35_input_full_copy_С2[[CD8_parameter_General]])),
    CD8_parameter_HighRes = median(as.numeric(d35_input_full_copy_С2[[CD8_parameter_HighRes]])),
    CD56dimNK_parameter_General = median(as.numeric(d35_input_full_copy_С2[[CD56dimNK_parameter_General]])),
    CD56dimNK_parameter_HighRes = median(as.numeric(d35_input_full_copy_С2[[CD56dimNK_parameter_HighRes]])),
    CD4TCM_parameter = median(as.numeric(d35_input_full_copy_С2[[CD4TCM_parameter]])),
    CD4TCMNVratio_parameter = median(as.numeric(d35_input_full_copy_С2[[CD4TCMNVratio_parameter]]))
  )
)
# -----

# ----- score definition: -----
# d35_input_full_copy_С1$CD8_score <- ifelse(as.numeric(d35_input_full_copy_С1[[CD8_parameter_General]]) > medians$cohort1$CD8_parameter_General, 'CD8hi', 'CD8lo')
d35_input_full_copy_С1$CD8_score <- ifelse(as.numeric(d35_input_full_copy_С1[[CD8_parameter_HighRes]]) > medians$cohort1$CD8_parameter_HighRes, 'CD8hi', 'CD8lo')
# d35_input_full_copy_С1$CD4TCM_score <- ifelse(as.numeric(d35_input_full_copy_С1[[CD4TCM_parameter]]) > medians$cohort1$CD4TCM_parameter, 0, 1)
# d35_input_full_copy_С1$CD4TCMNVratio_score <- ifelse(as.numeric(d35_input_full_copy_С1[[CD4TCMNVratio_parameter]]) > medians$cohort1$CD4TCMNVratio_parameter, 0, 1)
# d35_input_full_copy_С1$CD4ScoreOverall <- d35_input_full_copy_С1$CD4TCM_score + d35_input_full_copy_С1$CD4TCMNVratio_score
# d35_input_full_copy_С1$CD4ScoreOverall <- ifelse(d35_input_full_copy_С1$CD4ScoreOverall == 2, 'CD4TCM/NVlo', 'CD4TCM/NVhi')
# d35_input_full_copy_С1$FinalScore <- paste(d35_input_full_copy_С1$CD8_score, d35_input_full_copy_С1$CD4ScoreOverall, sep = " & ")

# d35_input_full_copy_С2$CD8_score <- ifelse(as.numeric(d35_input_full_copy_С2[[CD8_parameter_General]]) > medians$cohort2$CD8_parameter_General, 'CD8hi', 'CD8lo')
d35_input_full_copy_С2$CD8_score <- ifelse(as.numeric(d35_input_full_copy_С2[[CD8_parameter_HighRes]]) > medians$cohort2$CD8_parameter_HighRes, 'CD8hi', 'CD8lo')
# d35_input_full_copy_С2$CD4TCM_score <- ifelse(as.numeric(d35_input_full_copy_С2[[CD4TCM_parameter]]) > medians$cohort2$CD4TCM_parameter, 0, 1)
# d35_input_full_copy_С2$CD4TCMNVratio_score <- ifelse(as.numeric(d35_input_full_copy_С2[[CD4TCMNVratio_parameter]]) > medians$cohort2$CD4TCMNVratio_parameter, 0, 1)
# d35_input_full_copy_С2$CD4ScoreOverall <- d35_input_full_copy_С2$CD4TCM_score + d35_input_full_copy_С2$CD4TCMNVratio_score
# d35_input_full_copy_С2$CD4ScoreOverall <- ifelse(d35_input_full_copy_С2$CD4ScoreOverall == 2, 'CD4TCM/NVlo', 'CD4TCM/NVhi')
# d35_input_full_copy_С2$FinalScore <- paste(d35_input_full_copy_С2$CD8_score, d35_input_full_copy_С2$CD4ScoreOverall, sep = " & ")

# merging_pattern <- list(
#   # 'CD8hi\nCD4TCM/NVlo' = c('CD8hi & CD4TCM/NVlo'),
#   # 'CD8hi\nCD4TCM/NVhi' = c('CD8hi & CD4TCM/NVhi'),
#   # 'CD8lo' = c('CD8lo & CD4TCM/NVlo', 'CD8lo & CD4TCM/NVhi')
#   
#   'GvLhi_GvHDlo' = c('CD8hi & CD4TCM/NVlo'),
#   'Other' = c('CD8lo & CD4TCM/NVlo', 'CD8lo & CD4TCM/NVhi', 'CD8hi & CD4TCM/NVhi')
# )
# 
# lookup <- unlist(
#   lapply(names(merging_pattern), function(name) {
#     setNames(rep(name, length(merging_pattern[[name]])),
#              merging_pattern[[name]])
#   })
# )
# d35_input_full_copy_С1$FinalScore_Grouped <- lookup[as.character(d35_input_full_copy_С1$FinalScore)]
# d35_input_full_copy_С2$FinalScore_Grouped <- lookup[as.character(d35_input_full_copy_С2$FinalScore)]
# -----

# ----- resetting RFS to 0/1/2 -----
df_cr_c1 <- d35_input_full_copy_С1 %>%
  mutate(
    cr_event = case_when(
      Relapse_Present_1_NotObserved_0 == 1 ~ 1,
      Relapse_Present_1_NotObserved_0 == 0 & Death_Present_1_NotObserved_0 == 1 ~ 2,
      TRUE ~ 0
    ),
    cr_time = case_when(
      Relapse_Present_1_NotObserved_0 == 1 ~ Time_to_Relapse_from_TPL,
      Relapse_Present_1_NotObserved_0 == 0 & Death_Present_1_NotObserved_0 == 1 ~ Time_to_Death_from_TPL,
      TRUE ~ Time_to_Relapse_from_TPL
    )
  )

df_cr_c2 <- d35_input_full_copy_С2 %>%
  mutate(
    cr_event = case_when(
      Relapse_Present_1_NotObserved_0 == 1 ~ 1,
      Relapse_Present_1_NotObserved_0 == 0 & Death_Present_1_NotObserved_0 == 1 ~ 2,
      TRUE ~ 0
    ),
    cr_time = case_when(
      Relapse_Present_1_NotObserved_0 == 1 ~ Time_to_Relapse_from_TPL,
      Relapse_Present_1_NotObserved_0 == 0 & Death_Present_1_NotObserved_0 == 1 ~ Time_to_Death_from_TPL,
      TRUE ~ Time_to_Relapse_from_TPL
    )
  )
# -----


ci_c1 <- cuminc(
  ftime = df_cr_c1$Time_to_Relapse_from_TPL,
  fstatus = df_cr_c1$Relapse_Present_1_NotObserved_0,
  group = df_cr_c1$CD8_score,
  cencode = 0
)
print(ci_c1$Tests)

ci_c2 <- cuminc(
  ftime = df_cr_c2$Time_to_Relapse_from_TPL,
  fstatus = df_cr_c2$Relapse_Present_1_NotObserved_0,
  group = df_cr_c2$CD8_score,
  cencode = 0
)
print(ci_c2$Tests)
