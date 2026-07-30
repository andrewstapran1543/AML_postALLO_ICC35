################################
#### Importing the packages ####
################################
renv::install("survminer")
renv::install(c("survminer", "ggplot2"))
library("survminer")
require("survival")
library("readxl")

library(readxl)
library(dplyr)
library(stringr)
library(glue)

######################################
#### Specifying horizoning cutoff ####
######################################
cutoff <- 2000  # set your cutoff





#####################################################################################
#### Reading the dataframe with clinical metadata and immune flow cytometry data ####
#####################################################################################
# path to the full dataset with clinical metadata and flow cytometry data
d35_input_full_path <- ''
d35_input_full <- read_excel(d35_input_full_path)

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
PERCENT_FEATURES_renamer <- PERCENT_FEATURES[str_detect(PERCENT_FEATURES, "LowLevel")]
names(PERCENT_FEATURES_renamer) <- PERCENT_FEATURES_renamer
PERCENT_FEATURES_renamer <- str_replace(PERCENT_FEATURES_renamer, "LowLevel ", "High Resolution ")

d35_input_full_copy <- d35_input_full_copy %>%
  rename_with(
    ~ str_replace(.x, "LowLevel ", "High Resolution "),
    all_of(names(PERCENT_FEATURES_renamer))
  )

PERCENT_FEATURES <- str_replace(PERCENT_FEATURES, "LowLevel ", "High Resolution ")

# Ratio features
RATIO_FEATURES <- names(d35_input_full_copy)[str_detect(names(d35_input_full_copy), "-to-")]
OUTCOME_COLS <- list(
  RFS = list(time = "Time_to_Relapse_from_TPL", event = "Relapse_Present_1_NotObserved_0"),
  OS = list(time = "Time_to_Death_from_TPL", event = "Death_Present_1_NotObserved_0"),
  aGVHD = list(time = "Time_to_aGVHD", event = "aGVHD_grading_by_Katja"),
  cGVHD = list(time = "Time_to_cGVHD", event = "cGVHD_grading_by_Katja")
)

d35_input_full_copy <- d35_input_full_copy %>% mutate(FC_Sampling_Days_After_Allo_d35 = as.integer(round(as.numeric(FC_Sampling_Days_After_Allo_d35), 0)))
d35_input_full_copy$Cutoff_Plus_SamplingDay <- cutoff + d35_input_full_copy$FC_Sampling_Days_After_Allo_d35

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

    MAX_FOLLOW_UP = pmax(
      Time_to_aGVHD,
      Time_to_cGVHD,
      Time_to_Relapse_from_TPL,
      Time_to_Death_from_TPL,
      na.rm = TRUE
    )
  ) %>%
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
# -----




##############################################################################
#### Score: CD8 TEM +/+ of CD8 TEM; CD4 TCM of T; CD4TCM/NV_to_rest ratio ####
#### Showing the KM plots for all patients for this score - OS as endpoint ###
##############################################################################

# ----- setting up medians & calculating the scores -----
d35_input_full_copy_С1 <- d35_input_full_copy %>% subset(Cohort == 'OC')
d35_input_full_copy_С2 <- d35_input_full_copy %>% subset(Cohort == 'VC')

CD8_parameter_General = "CD8 TEM adv1 GPR56+ CD57+ of CD8 TEM"
CD8_parameter_HighRes = "LowLevel CD8 CCR7- CD45RAmid CD27+ CD57+ GPR56+ of CD45pos"

CD56dimNK_parameter_General = "CD56dim NK CD57+ GPR56+ of Lymphocytes"
CD56dimNK_parameter_HighRes = "LowLevel CD56dim NK CD45RA+ CD57+ GPR56+ of Lymphocytes"

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

# score definition:
d35_input_full_copy_С1$CD8_score <- ifelse(as.numeric(d35_input_full_copy_С1[[CD8_parameter_General]]) > medians$cohort1$CD8_parameter_General, 'CD8hi', 'CD8lo')
d35_input_full_copy_С1$CD4TCM_score <- ifelse(as.numeric(d35_input_full_copy_С1[[CD4TCM_parameter]]) > medians$cohort1$CD4TCM_parameter, 0, 1)
d35_input_full_copy_С1$CD4TCMNVratio_score <- ifelse(as.numeric(d35_input_full_copy_С1[[CD4TCMNVratio_parameter]]) > medians$cohort1$CD4TCMNVratio_parameter, 0, 1)
d35_input_full_copy_С1$CD4ScoreOverall <- d35_input_full_copy_С1$CD4TCM_score + d35_input_full_copy_С1$CD4TCMNVratio_score
d35_input_full_copy_С1$CD4ScoreOverall <- ifelse(d35_input_full_copy_С1$CD4ScoreOverall == 2, 'CD4TCM/NVlo', 'CD4TCM/NVhi')
d35_input_full_copy_С1$FinalScore <- paste(d35_input_full_copy_С1$CD8_score, d35_input_full_copy_С1$CD4ScoreOverall, sep = " & ")

d35_input_full_copy_С2$CD8_score <- ifelse(as.numeric(d35_input_full_copy_С2[[CD8_parameter_General]]) > medians$cohort2$CD8_parameter_General, 'CD8hi', 'CD8lo')
d35_input_full_copy_С2$CD4TCM_score <- ifelse(as.numeric(d35_input_full_copy_С2[[CD4TCM_parameter]]) > medians$cohort2$CD4TCM_parameter, 0, 1)
d35_input_full_copy_С2$CD4TCMNVratio_score <- ifelse(as.numeric(d35_input_full_copy_С2[[CD4TCMNVratio_parameter]]) > medians$cohort2$CD4TCMNVratio_parameter, 0, 1)
d35_input_full_copy_С2$CD4ScoreOverall <- d35_input_full_copy_С2$CD4TCM_score + d35_input_full_copy_С2$CD4TCMNVratio_score
d35_input_full_copy_С2$CD4ScoreOverall <- ifelse(d35_input_full_copy_С2$CD4ScoreOverall == 2, 'CD4TCM/NVlo', 'CD4TCM/NVhi')
d35_input_full_copy_С2$FinalScore <- paste(d35_input_full_copy_С2$CD8_score, d35_input_full_copy_С2$CD4ScoreOverall, sep = " & ")
# -----

# ---- merging aGvHD grades ----
merging_pattern <- list(
  # 'CD8hi\nCD4TCM/NVlo' = c('CD8hi & CD4TCM/NVlo'),
  # 'CD8hi\nCD4TCM/NVhi' = c('CD8hi & CD4TCM/NVhi'),
  # 'CD8lo' = c('CD8lo & CD4TCM/NVlo', 'CD8lo & CD4TCM/NVhi')

  'GvLhi_GvHDlo' = c('CD8hi & CD4TCM/NVlo'),
  'Other' = c('CD8lo & CD4TCM/NVlo', 'CD8lo & CD4TCM/NVhi', 'CD8hi & CD4TCM/NVhi')
)

lookup <- unlist(
  lapply(names(merging_pattern), function(name) {
    setNames(rep(name, length(merging_pattern[[name]])),
             merging_pattern[[name]])
  })
)
d35_input_full_copy_С1$FinalScore_Grouped <- lookup[as.character(d35_input_full_copy_С1$FinalScore)]
d35_input_full_copy_С2$FinalScore_Grouped <- lookup[as.character(d35_input_full_copy_С2$FinalScore)]

cohort_of_interest <- d35_input_full_copy_С2
time_col <- OUTCOME_COLS[["OS"]]$time
event_col <- OUTCOME_COLS[["OS"]]$event
fit <- survfit(
  Surv(cohort_of_interest[[time_col]], cohort_of_interest[[event_col]]) ~ FinalScore_Grouped,
  data = cohort_of_interest
)

# specifying the output directory for KM curves:
save_dir <- ''
custom_palette = c("#4DAF4A", "#377EB8", "#A65628", "#E41A1C")
# custom_palette = c("#377EB8", "#4DAF4A", "#A65628", "#E41A1C")

custom_legend_title = "Final score categories:"
custom_legend_labs <- gsub("FinalScore_Grouped=", "", names(fit$strata))
# file_name = 'CD8General_CD4TCM_CD4TCMNVratio_Cohort2.pdf'
file_name = 'CD8General_CD4TCM_CD4TCMNVratio_Cohort2_2curves.pdf'


fontsize_small <- 8
fontsize_big <- 10

# ----- printing the plot (At Risk, Events, Censored) -----
p <- ggsurvplot(
  fit, 
  data = cohort_of_interest, 
  size = 1,
  palette = custom_palette,
  conf.int = FALSE,
  pval = FALSE,
  risk.table = TRUE,
  cumevents = TRUE,
  cumcensor = TRUE,
  # tables.height = 0.15,
  risk.table.height = 0.12,
  cumevents.height = 0.12,
  cumcensor.height = 0.16,
  
  legend.title = custom_legend_title,
  legend.labs = custom_legend_labs,
  ggtheme = theme_bw(),
  fontsize = fontsize_small / ggplot2::.pt,
  break.time.by = 365
)

# shared table style
table_theme <- theme(
  panel.border = element_rect(colour = "grey60", fill = NA, linewidth = 0.3),
  panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
  panel.grid.minor = element_blank()
)

fix_table_theme <- function(x, ylab, show_x_axis = FALSE) {
  x <- x + table_theme + labs(y = ylab)
  
  x$theme <- x$theme %+replace% theme(
    axis.text.y = element_text(size = fontsize_small),
    axis.title.y = element_text(
      size = fontsize_big,
      angle = 90,
      vjust = 0.5,
      margin = margin(r = 2)
    ),
    axis.ticks.y = element_blank(),
    plot.title = element_blank(),
    
    axis.title.x = if (show_x_axis) {
      element_text(
        size = fontsize_big,
        margin = margin(t = 2)
      )
    } else {
      element_blank()
    },
    axis.text.x  = if (show_x_axis) element_text(size = fontsize_small) else element_blank(),
    axis.ticks.x = if (show_x_axis) element_line(linewidth = 0.3) else element_blank()
  )
  
  x
}

p$table <- fix_table_theme(p$table, "At risk", show_x_axis = FALSE)
p$cumevents <- fix_table_theme(p$cumevents, "Events", show_x_axis = FALSE)
p$ncensor.plot <- fix_table_theme(p$ncensor.plot, "Censored", show_x_axis = TRUE) +
  labs(x = "Days post alloHCT sampling")

pval <- surv_pvalue(fit, data = cohort_of_interest)$pval.txt

p$plot <- p$plot +
  labs(
    x = NULL,
    y = "Overall survival"
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      keywidth = unit(0.7, "cm"),
      keyheight = unit(0.5, "cm"),
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  theme(
    legend.position = c(0.35, 0.225),
    legend.text = element_text(size = fontsize_big),
    legend.title = element_text(size = fontsize_big, hjust = 0.5),
    legend.spacing.x = unit(0.2, "cm"),
    legend.margin = margin(3, 3, 3, 3),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
    legend.box.background = element_blank(),
    legend.key = element_rect(fill = "white"),
    
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = fontsize_small),
    axis.title.y = element_text(size = fontsize_big)
  ) +
  annotate(
    "label",
    x = 1400,
    y = 0.2,
    label = pval,
    size = fontsize_big / ggplot2::.pt,
    fill = "white",
    color = "black",
    label.r = unit(0, "lines"),
    linewidth = 0.4
  )
# -----

# ----- printing the plot (At Risk) -----
p <- ggsurvplot(
  fit, 
  data = cohort_of_interest, 
  size = 1,
  palette = custom_palette,
  conf.int = FALSE,
  pval = FALSE,
  risk.table = TRUE,
  # cumevents = TRUE,
  # cumcensor = TRUE,
  # tables.height = 0.15,
  risk.table.height = 0.35,
  # cumevents.height = 0.12,
  # cumcensor.height = 0.16,
  
  legend.title = custom_legend_title,
  legend.labs = custom_legend_labs,
  ggtheme = theme_bw(),
  fontsize = fontsize_small / ggplot2::.pt,
  break.time.by = 500
)

# shared table style
table_theme <- theme(
  panel.border = element_rect(colour = "grey60", fill = NA, linewidth = 0.3),
  panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
  panel.grid.minor = element_blank()
)

fix_table_theme <- function(x, ylab, show_x_axis = FALSE) {
  x <- x + table_theme + labs(y = ylab)
  
  x$theme <- x$theme %+replace% theme(
    axis.text.y = element_text(size = fontsize_small),
    axis.title.y = element_text(
      size = fontsize_big,
      angle = 90,
      vjust = 0.5,
      margin = margin(r = 2)
    ),
    axis.ticks.y = element_blank(),
    plot.title = element_blank(),
    
    axis.title.x = if (show_x_axis) {
      element_text(
        size = fontsize_big,
        margin = margin(t = 2)
      )
    } else {
      element_blank()
    },
    axis.text.x  = if (show_x_axis) element_text(size = fontsize_small) else element_blank(),
    axis.ticks.x = if (show_x_axis) element_line(linewidth = 0.3) else element_blank()
  )
  
  x
}

p$table <- fix_table_theme(p$table, "At risk", show_x_axis = TRUE) +
  labs(x = "Days post alloHCT sampling")
pval <- surv_pvalue(fit, data = cohort_of_interest)$pval.txt

p$plot <- p$plot +
  labs(
    x = NULL,
    y = "Overall survival"
  ) +
  guides(
    color = guide_legend(
      ncol = 1,
      byrow = TRUE,
      keywidth = unit(0.7, "cm"),
      keyheight = unit(0.5, "cm"),
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  theme(
    legend.position = c(0.25, 0.225),
    legend.text = element_text(size = fontsize_big),
    legend.title = element_text(size = fontsize_big, hjust = 0.5),
    legend.spacing.x = unit(0.2, "cm"),
    legend.margin = margin(3, 3, 3, 3),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
    legend.box.background = element_blank(),
    legend.key = element_rect(fill = "white"),
    
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = fontsize_small),
    axis.title.y = element_text(size = fontsize_big)
  ) +
  annotate(
    "label",
    x = 1400,
    y = 0.2,
    label = pval,
    size = fontsize_big / ggplot2::.pt,
    fill = "white",
    color = "black",
    label.r = unit(0, "lines"),
    linewidth = 0.4
  )
# -----

pdf(glue('{save_dir}/{file_name}'), width = 6, height = 5)
print(p, newpage = FALSE)
dev.off()
# -----








##############################################################################
#### Score: CD8 TEM +/+ of CD8 TEM; CD4 TCM of T; CD4TCM/NV_to_rest ratio ####
#### CD8TEM lo patients - KM plots with cGVHD grades with OS as endpoint  ####
##############################################################################
# ---- merging aGvHD grades ----
categories_analysed <- c(0,1)
selected_CD8 <- 'CD8hi'
# 'CD8lo' / 'CD8hi'

# ---- censoring -----
censor_window <- 180
time_col <- OUTCOME_COLS[["OS"]]$time
event_col <- OUTCOME_COLS[["OS"]]$event

d35_input_full_copy_С1_CD8lo <- subset(d35_input_full_copy_С1, CD8_score == selected_CD8)
d35_input_full_copy_С1_CD8lo$FU_first_180days <- ifelse(
  d35_input_full_copy_С1_CD8lo[[time_col]] + d35_input_full_copy_С1_CD8lo$FC_Sampling_Days_After_Allo_d35 < 180,
  'lost_to_FU_within_first_180days', 'sufficient_FU')
d35_input_full_copy_С1_CD8lo <- subset(d35_input_full_copy_С1_CD8lo, (FU_first_180days == 'sufficient_FU') & (cGVHD_grading_by_Katja %in% categories_analysed))



d35_input_full_copy_С2_CD8lo <- subset(d35_input_full_copy_С2, CD8_score == selected_CD8)
d35_input_full_copy_С2_CD8lo$FU_first_180days <- ifelse(
  d35_input_full_copy_С2_CD8lo[[time_col]] + d35_input_full_copy_С2_CD8lo$FC_Sampling_Days_After_Allo_d35 < 180,
  'lost_to_FU_within_first_180days', 'sufficient_FU')
d35_input_full_copy_С2_CD8lo <- subset(d35_input_full_copy_С2_CD8lo, (FU_first_180days == 'sufficient_FU') & (cGVHD_grading_by_Katja %in% categories_analysed))
# -----

merging_pattern <- list(
  # "Grade 0" = c(0),
  # "Grade 1" = c(1),
  # "Grade 2" = c(2)

  # "Grade 0" = c(0),
  # "Grade 1-2" = c(1,2)
  
  "Grade 0" = c(0),
  "Grade 1" = c(1)
)

lookup <- unlist(
  lapply(names(merging_pattern), function(name) {
    setNames(rep(name, length(merging_pattern[[name]])),
             merging_pattern[[name]])
  })
)
d35_input_full_copy_С1_CD8lo$cGvHD_grouped <- lookup[as.character(d35_input_full_copy_С1_CD8lo$cGVHD_grading_by_Katja)]
d35_input_full_copy_С2_CD8lo$cGvHD_grouped <- lookup[as.character(d35_input_full_copy_С2_CD8lo$cGVHD_grading_by_Katja)]

cohort_of_interest <- d35_input_full_copy_С2_CD8lo
fit <- survfit(
  Surv(cohort_of_interest[[time_col]], cohort_of_interest[[event_col]]) ~ cGvHD_grouped,
  data = cohort_of_interest
)

# specifying the output directory for KM curves:
save_dir <- ''
custom_palette = c("#377EB8", "#4DAF4A", "#A65628", "#E41A1C")
custom_legend_title = "cGvHD grades:"
custom_legend_labs <- gsub("cGvHD_grouped=", "", names(fit$strata))
# file_name = 'GvLScore_CD8TEMlo_cGVHD_grades_0_1_Cohort2.pdf'
file_name = 'GvLScore_CD8TEMhi_cGVHD_grades_0_1_Cohort2.pdf'

fontsize_small <- 8
fontsize_big <- 10

# ----- printing the plot (At Risk, Events, Censored) -----
p <- ggsurvplot(
  fit, 
  data = cohort_of_interest, 
  size = 1,
  palette = custom_palette,
  conf.int = FALSE,
  pval = FALSE,
  risk.table = TRUE,
  cumevents = TRUE,
  cumcensor = TRUE,
  # tables.height = 0.15,
  risk.table.height = 0.12,
  cumevents.height = 0.12,
  cumcensor.height = 0.16,
  
  legend.title = custom_legend_title,
  legend.labs = custom_legend_labs,
  ggtheme = theme_bw(),
  fontsize = fontsize_small / ggplot2::.pt,
  break.time.by = 365
)

# shared table style
table_theme <- theme(
  panel.border = element_rect(colour = "grey60", fill = NA, linewidth = 0.3),
  panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
  panel.grid.minor = element_blank()
)

fix_table_theme <- function(x, ylab, show_x_axis = FALSE) {
  x <- x + table_theme + labs(y = ylab)
  
  x$theme <- x$theme %+replace% theme(
    axis.text.y = element_text(size = fontsize_small),
    axis.title.y = element_text(
      size = fontsize_big,
      angle = 90,
      vjust = 0.5,
      margin = margin(r = 2)
    ),
    axis.ticks.y = element_blank(),
    plot.title = element_blank(),
    
    axis.title.x = if (show_x_axis) {
      element_text(
        size = fontsize_big,
        margin = margin(t = 2)
      )
    } else {
      element_blank()
    },
    axis.text.x  = if (show_x_axis) element_text(size = fontsize_small) else element_blank(),
    axis.ticks.x = if (show_x_axis) element_line(linewidth = 0.3) else element_blank()
  )
  
  x
}

p$table <- fix_table_theme(p$table, "At risk", show_x_axis = FALSE)
p$cumevents <- fix_table_theme(p$cumevents, "Events", show_x_axis = FALSE)
p$ncensor.plot <- fix_table_theme(p$ncensor.plot, "Censored", show_x_axis = TRUE) +
  labs(x = "Days post alloHCT sampling")

pval <- surv_pvalue(fit, data = cohort_of_interest)$pval.txt

p$plot <- p$plot +
  labs(
    x = NULL,
    y = "Overall survival"
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      keywidth = unit(0.7, "cm"),
      keyheight = unit(0.5, "cm"),
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  theme(
    legend.position = c(0.35, 0.225),
    legend.text = element_text(size = fontsize_big),
    legend.title = element_text(size = fontsize_big, hjust = 0.5),
    legend.spacing.x = unit(0.2, "cm"),
    legend.margin = margin(3, 3, 3, 3),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
    legend.box.background = element_blank(),
    legend.key = element_rect(fill = "white"),
    
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = fontsize_small),
    axis.title.y = element_text(size = fontsize_big)
  ) +
  annotate(
    "label",
    x = 1400,
    y = 0.2,
    label = pval,
    size = fontsize_big / ggplot2::.pt,
    fill = "white",
    color = "black",
    label.r = unit(0, "lines"),
    linewidth = 0.4
  )
# -----

# ----- printing the plot (At Risk) -----
p <- ggsurvplot(
  fit, 
  data = cohort_of_interest, 
  size = 1,
  palette = custom_palette,
  conf.int = FALSE,
  pval = FALSE,
  risk.table = TRUE,
  # cumevents = TRUE,
  # cumcensor = TRUE,
  # tables.height = 0.15,
  risk.table.height = 0.35,
  # cumevents.height = 0.12,
  # cumcensor.height = 0.16,
  
  legend.title = custom_legend_title,
  legend.labs = custom_legend_labs,
  ggtheme = theme_bw(),
  fontsize = fontsize_small / ggplot2::.pt,
  break.time.by = 500
)

# shared table style
table_theme <- theme(
  panel.border = element_rect(colour = "grey60", fill = NA, linewidth = 0.3),
  panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
  panel.grid.minor = element_blank()
)

fix_table_theme <- function(x, ylab, show_x_axis = FALSE) {
  x <- x + table_theme + labs(y = ylab)
  
  x$theme <- x$theme %+replace% theme(
    axis.text.y = element_text(size = fontsize_small),
    axis.title.y = element_text(
      size = fontsize_big,
      angle = 90,
      vjust = 0.5,
      margin = margin(r = 2)
    ),
    axis.ticks.y = element_blank(),
    plot.title = element_blank(),
    
    axis.title.x = if (show_x_axis) {
      element_text(
        size = fontsize_big,
        margin = margin(t = 2)
      )
    } else {
      element_blank()
    },
    axis.text.x  = if (show_x_axis) element_text(size = fontsize_small) else element_blank(),
    axis.ticks.x = if (show_x_axis) element_line(linewidth = 0.3) else element_blank()
  )
  
  x
}

p$table <- fix_table_theme(p$table, "At risk", show_x_axis = TRUE) +
  labs(x = "Days post alloHCT sampling")
pval <- surv_pvalue(fit, data = cohort_of_interest)$pval.txt

p$plot <- p$plot +
  labs(
    x = NULL,
    y = "Overall survival"
  ) +
  guides(
    color = guide_legend(
      ncol = 1,
      byrow = TRUE,
      keywidth = unit(0.7, "cm"),
      keyheight = unit(0.5, "cm"),
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  theme(
    legend.position = c(0.25, 0.225),
    legend.text = element_text(size = fontsize_big),
    legend.title = element_text(size = fontsize_big, hjust = 0.5),
    legend.spacing.x = unit(0.2, "cm"),
    legend.margin = margin(3, 3, 3, 3),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
    legend.box.background = element_blank(),
    legend.key = element_rect(fill = "white"),
    
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = fontsize_small),
    axis.title.y = element_text(size = fontsize_big)
  ) +
  annotate(
    "label",
    x = 1400,
    y = 0.2,
    label = pval,
    size = fontsize_big / ggplot2::.pt,
    fill = "white",
    color = "black",
    label.r = unit(0, "lines"),
    linewidth = 0.4
  )
# -----

pdf(glue('{save_dir}/{file_name}'), width = 6, height = 5)
print(p, newpage = FALSE)
dev.off()
# -----











###################################################################################################
#### Score: CD8 TEM +/+ of CD8 TEM; CD4 TCM of T / CD4 TCM/NV to rest; CD56dim NK +/+ general  ####
####          Showing the KM plots for all patients for this score - OS as endpoint            ####
###################################################################################################

# ----- setting up medians & calculating the scores -----
d35_input_full_copy_С1 <- d35_input_full_copy %>% subset(Cohort == 'OC')
d35_input_full_copy_С2 <- d35_input_full_copy %>% subset(Cohort == 'VC')

CD8_parameter_General = "CD8 TEM adv1 GPR56+ CD57+ of CD8 TEM"
CD8_parameter_HighRes = "LowLevel CD8 CCR7- CD45RAmid CD27+ CD57+ GPR56+ of CD45pos"

CD56dimNK_parameter_General = "CD56dim NK CD57+ GPR56+ of Lymphocytes"
CD56dimNK_parameter_HighRes = "LowLevel CD56dim NK CD45RA+ CD57+ GPR56+ of Lymphocytes"

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

# score definition:
d35_input_full_copy_С1$CD8_score <- ifelse(as.numeric(d35_input_full_copy_С1[[CD8_parameter_General]]) > medians$cohort1$CD8_parameter_General, 'CD8hi', 'CD8lo')
d35_input_full_copy_С1$CD4_score <- ifelse(as.numeric(d35_input_full_copy_С1[[CD4TCM_parameter]]) > medians$cohort1$CD4TCM_parameter, 0, 1)
# d35_input_full_copy_С1$CD4_score <- ifelse(as.numeric(d35_input_full_copy_С1[[CD4TCMNVratio_parameter]]) > medians$cohort1$CD4TCMNVratio_parameter, 0, 1)
# d35_input_full_copy_С1$CD56dimNK_general_score <- ifelse(as.numeric(d35_input_full_copy_С1[[CD56dimNK_parameter_General]]) > medians$cohort1$CD56dimNK_parameter_General, 1, 0)
# d35_input_full_copy_С1$GvHD_Score <- d35_input_full_copy_С1$CD4_score + d35_input_full_copy_С1$CD56dimNK_general_score
d35_input_full_copy_С1$CD56dimNK_highres_score <- ifelse(as.numeric(d35_input_full_copy_С1[[CD56dimNK_parameter_HighRes]]) > medians$cohort1$CD56dimNK_parameter_HighRes, 1, 0)
d35_input_full_copy_С1$GvHD_Score <- d35_input_full_copy_С1$CD4_score + d35_input_full_copy_С1$CD56dimNK_highres_score
d35_input_full_copy_С1$GvHD_Score <- ifelse(d35_input_full_copy_С1$GvHD_Score == 2, 'GvHDlo', 'GvHDhi')
d35_input_full_copy_С1$FinalScore <- paste(d35_input_full_copy_С1$CD8_score, d35_input_full_copy_С1$GvHD_Score, sep = " & ")

d35_input_full_copy_С2$CD8_score <- ifelse(as.numeric(d35_input_full_copy_С2[[CD8_parameter_General]]) > medians$cohort2$CD8_parameter_General, 'CD8hi', 'CD8lo')
d35_input_full_copy_С2$CD4_score <- ifelse(as.numeric(d35_input_full_copy_С2[[CD4TCM_parameter]]) > medians$cohort2$CD4TCM_parameter, 0, 1)
# d35_input_full_copy_С2$CD4_score <- ifelse(as.numeric(d35_input_full_copy_С2[[CD4TCMNVratio_parameter]]) > medians$cohort2$CD4TCMNVratio_parameter, 0, 1)
# d35_input_full_copy_С2$CD56dimNK_general_score <- ifelse(as.numeric(d35_input_full_copy_С2[[CD56dimNK_parameter_General]]) > medians$cohort2$CD56dimNK_parameter_General, 1, 0)
# d35_input_full_copy_С2$GvHD_Score <- d35_input_full_copy_С2$CD4_score + d35_input_full_copy_С2$CD56dimNK_general_score
d35_input_full_copy_С2$CD56dimNK_highres_score <- ifelse(as.numeric(d35_input_full_copy_С2[[CD56dimNK_parameter_HighRes]]) > medians$cohort2$CD56dimNK_parameter_HighRes, 1, 0)
d35_input_full_copy_С2$GvHD_Score <- d35_input_full_copy_С2$CD4_score + d35_input_full_copy_С2$CD56dimNK_highres_score
d35_input_full_copy_С2$GvHD_Score <- ifelse(d35_input_full_copy_С2$GvHD_Score == 2, 'GvHDlo', 'GvHDhi')
d35_input_full_copy_С2$FinalScore <- paste(d35_input_full_copy_С2$CD8_score, d35_input_full_copy_С2$GvHD_Score, sep = " & ")
# -----

# ---- merging aGvHD grades ----
merging_pattern <- list(
  'GvLhi_GvHDlo' = c("CD8hi & GvHDlo"),
  'Other' = c("CD8hi & GvHDhi", "CD8lo & GvHDlo", "CD8lo & GvHDhi")
)

lookup <- unlist(
  lapply(names(merging_pattern), function(name) {
    setNames(rep(name, length(merging_pattern[[name]])),
             merging_pattern[[name]])
  })
)
d35_input_full_copy_С1$FinalScore_Grouped <- lookup[as.character(d35_input_full_copy_С1$FinalScore)]
d35_input_full_copy_С2$FinalScore_Grouped <- lookup[as.character(d35_input_full_copy_С2$FinalScore)]

cohort_of_interest <- d35_input_full_copy_С2
time_col <- OUTCOME_COLS[["OS"]]$time
event_col <- OUTCOME_COLS[["OS"]]$event
fit <- survfit(
  Surv(cohort_of_interest[[time_col]], cohort_of_interest[[event_col]]) ~ FinalScore_Grouped,
  data = cohort_of_interest
)

# specifying the output directory for KM curves:
save_dir <- ''
custom_palette = c("#4DAF4A", "#377EB8", "#A65628", "#E41A1C")
custom_legend_title = "Final score categories:"
custom_legend_labs <- gsub("FinalScore_Grouped=", "", names(fit$strata))
# file_name = 'CD8General_CD4TCMNVratio_HighResCD56dimNK_Cohort2.pdf'
file_name = 'CD8General_CD4TCM_HighResCD56dimNK_Cohort2.pdf'


fontsize_small <- 8
fontsize_big <- 10

# ----- printing the plot (At Risk, Events, Censored) -----
p <- ggsurvplot(
  fit, 
  data = cohort_of_interest, 
  size = 1,
  palette = custom_palette,
  conf.int = FALSE,
  pval = FALSE,
  risk.table = TRUE,
  cumevents = TRUE,
  cumcensor = TRUE,
  # tables.height = 0.15,
  risk.table.height = 0.12,
  cumevents.height = 0.12,
  cumcensor.height = 0.16,
  
  legend.title = custom_legend_title,
  legend.labs = custom_legend_labs,
  ggtheme = theme_bw(),
  fontsize = fontsize_small / ggplot2::.pt,
  break.time.by = 365
)

# shared table style
table_theme <- theme(
  panel.border = element_rect(colour = "grey60", fill = NA, linewidth = 0.3),
  panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
  panel.grid.minor = element_blank()
)

fix_table_theme <- function(x, ylab, show_x_axis = FALSE) {
  x <- x + table_theme + labs(y = ylab)
  
  x$theme <- x$theme %+replace% theme(
    axis.text.y = element_text(size = fontsize_small),
    axis.title.y = element_text(
      size = fontsize_big,
      angle = 90,
      vjust = 0.5,
      margin = margin(r = 2)
    ),
    axis.ticks.y = element_blank(),
    plot.title = element_blank(),
    
    axis.title.x = if (show_x_axis) {
      element_text(
        size = fontsize_big,
        margin = margin(t = 2)
      )
    } else {
      element_blank()
    },
    axis.text.x  = if (show_x_axis) element_text(size = fontsize_small) else element_blank(),
    axis.ticks.x = if (show_x_axis) element_line(linewidth = 0.3) else element_blank()
  )
  
  x
}

p$table <- fix_table_theme(p$table, "At risk", show_x_axis = FALSE)
p$cumevents <- fix_table_theme(p$cumevents, "Events", show_x_axis = FALSE)
p$ncensor.plot <- fix_table_theme(p$ncensor.plot, "Censored", show_x_axis = TRUE) +
  labs(x = "Days post alloHCT sampling")

pval <- surv_pvalue(fit, data = cohort_of_interest)$pval.txt

p$plot <- p$plot +
  labs(
    x = NULL,
    y = "Overall survival"
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      keywidth = unit(0.7, "cm"),
      keyheight = unit(0.5, "cm"),
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  theme(
    legend.position = c(0.35, 0.225),
    legend.text = element_text(size = fontsize_big),
    legend.title = element_text(size = fontsize_big, hjust = 0.5),
    legend.spacing.x = unit(0.2, "cm"),
    legend.margin = margin(3, 3, 3, 3),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
    legend.box.background = element_blank(),
    legend.key = element_rect(fill = "white"),
    
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = fontsize_small),
    axis.title.y = element_text(size = fontsize_big)
  ) +
  annotate(
    "label",
    x = 1400,
    y = 0.2,
    label = pval,
    size = fontsize_big / ggplot2::.pt,
    fill = "white",
    color = "black",
    label.r = unit(0, "lines"),
    linewidth = 0.4
  )
# -----

# ----- printing the plot (At Risk) -----
p <- ggsurvplot(
  fit, 
  data = cohort_of_interest, 
  size = 1,
  palette = custom_palette,
  conf.int = FALSE,
  pval = FALSE,
  risk.table = TRUE,
  # cumevents = TRUE,
  # cumcensor = TRUE,
  # tables.height = 0.15,
  risk.table.height = 0.35,
  # cumevents.height = 0.12,
  # cumcensor.height = 0.16,
  
  legend.title = custom_legend_title,
  legend.labs = custom_legend_labs,
  ggtheme = theme_bw(),
  fontsize = fontsize_small / ggplot2::.pt,
  break.time.by = 500
)

# shared table style
table_theme <- theme(
  panel.border = element_rect(colour = "grey60", fill = NA, linewidth = 0.3),
  panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
  panel.grid.minor = element_blank()
)

fix_table_theme <- function(x, ylab, show_x_axis = FALSE) {
  x <- x + table_theme + labs(y = ylab)
  
  x$theme <- x$theme %+replace% theme(
    axis.text.y = element_text(size = fontsize_small),
    axis.title.y = element_text(
      size = fontsize_big,
      angle = 90,
      vjust = 0.5,
      margin = margin(r = 2)
    ),
    axis.ticks.y = element_blank(),
    plot.title = element_blank(),
    
    axis.title.x = if (show_x_axis) {
      element_text(
        size = fontsize_big,
        margin = margin(t = 2)
      )
    } else {
      element_blank()
    },
    axis.text.x  = if (show_x_axis) element_text(size = fontsize_small) else element_blank(),
    axis.ticks.x = if (show_x_axis) element_line(linewidth = 0.3) else element_blank()
  )
  
  x
}

p$table <- fix_table_theme(p$table, "At risk", show_x_axis = TRUE) +
  labs(x = "Days post alloHCT sampling")
pval <- surv_pvalue(fit, data = cohort_of_interest)$pval.txt

p$plot <- p$plot +
  labs(
    x = NULL,
    y = "Overall survival"
  ) +
  guides(
    color = guide_legend(
      ncol = 1,
      byrow = TRUE,
      keywidth = unit(0.7, "cm"),
      keyheight = unit(0.5, "cm"),
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  theme(
    legend.position = c(0.25, 0.225),
    legend.text = element_text(size = fontsize_big),
    legend.title = element_text(size = fontsize_big, hjust = 0.5),
    legend.spacing.x = unit(0.2, "cm"),
    legend.margin = margin(3, 3, 3, 3),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
    legend.box.background = element_blank(),
    legend.key = element_rect(fill = "white"),
    
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = fontsize_small),
    axis.title.y = element_text(size = fontsize_big)
  ) +
  annotate(
    "label",
    x = 1400,
    y = 0.2,
    label = pval,
    size = fontsize_big / ggplot2::.pt,
    fill = "white",
    color = "black",
    label.r = unit(0, "lines"),
    linewidth = 0.4
  )
# -----

pdf(glue('{save_dir}/{file_name}'), width = 6, height = 5)
print(p, newpage = FALSE)
dev.off()
# -----







