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


##########################################################
#### Specifying horizoning cutoff and break intervals ####
##########################################################
cutoff <- 750  # set your cutoff
landmarked <- FALSE
break_time_by <- 150

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





############################################
#### Introducing extra metadata columns ####
############################################
d35_input_full_copy$Cohort <- ifelse(d35_input_full_copy$Cohort == 'OC', 'Cohort1', 'Cohort2')
d35_input_full_copy$ELN_Score_AdvGenetics <- paste0('ELN', d35_input_full_copy$ELN_Score)

# ----- monosomal KT / complex KT / chr5 aberration -----
cols_to_clear <- c("monosomal_KT_1yes", "complexKT_1yes", "Monosomy5_Del5")
d35_input_full_copy <- d35_input_full_copy[complete.cases(d35_input_full_copy[cols_to_clear]) & apply(d35_input_full_copy[cols_to_clear] != "NA", 1, all),]
d35_input_full_copy$`Monosomal KT | Complex KT | Chr5 del` <- ifelse(
  (d35_input_full_copy$monosomal_KT_1yes == 1) | (d35_input_full_copy$complexKT_1yes == 1) | (d35_input_full_copy$Monosomy5_Del5 == 1),
  'Yes', 'No'
)
# -----

# ----- ELN3 + MDS-related mutations in Cohort 2 -----
d35_input_full_copy$ELN_Score_AdvGenetics <- ifelse(
  (d35_input_full_copy$ELN_Score_AdvGenetics == 'ELN3') & (d35_input_full_copy$Cohort == 'Cohort2'),
  ifelse(d35_input_full_copy$`MDS-related molecular_mutations_0_1_2is2ormore` %in% c(1,2), 'ELN3 with\nMDS mutation', 'ELN3 w/o\nMDS mutation'),
  d35_input_full_copy$ELN_Score_AdvGenetics
)
# -----

# ----- ELN3 + MDS-related mutations + monosomal KT in Cohort 2 -----
d35_input_full_copy$ELN_Score_AdvGenetics_MDSandMonKT <- ifelse(
  (d35_input_full_copy$ELN_Score_AdvGenetics == 'ELN3') & (d35_input_full_copy$Cohort == 'Cohort2') & (d35_input_full_copy$`MDS-related molecular_mutations_0_1_2is2ormore` %in% c(0)) & (d35_input_full_copy$monosomal_KT_1yes == 1),
  'ELN3 w/o MDS mut\nwith mon.KT', d35_input_full_copy$ELN_Score_AdvGenetics
)
# -----


# ----- ELN3 + monosomal KT -----
cols_to_clear <- c("monosomal_KT_1yes")
d35_input_full_copy <- d35_input_full_copy[complete.cases(d35_input_full_copy[cols_to_clear]) & apply(d35_input_full_copy[cols_to_clear] != "NA", 1, all),]
d35_input_full_copy$Monosomal_KT <- ifelse(d35_input_full_copy$monosomal_KT_1yes == 1,'Yes', 'No')
d35_input_full_copy$ELN_Score_Monosomal_KT <- ifelse(
  (d35_input_full_copy$ELN_Score_AdvGenetics == 'ELN3'),
  ifelse(d35_input_full_copy$monosomal_KT_1yes == 1, 'ELN3 with\nmonosomal KT', 'ELN3 w/o\nmonosomal KT'),
  d35_input_full_copy$ELN_Score_AdvGenetics
)
d35_input_full_copy$ELN_Score_Monosomal_KT_3vsothers <- ifelse(
  (d35_input_full_copy$ELN_Score_AdvGenetics == 'ELN3') & (d35_input_full_copy$monosomal_KT_1yes == 1),
  'ELN3 with\nmonosomal KT', 'Other'
)
# -----





########################################################################
#### Kaplan-Meier curves: splitting by clinical covariate of choice ####
########################################################################
cohort_selected <- "Cohort1"
parameter <- 'Monosomal KT | Complex KT | Chr5 del'
# 'ELN_Score'
# 'Monosomal KT | Complex KT | Chr5 del'
# 'ELN_Score_AdvGenetics'
# 'Monosomal_KT'
# 'ELN_Score_Monosomal_KT'

# plotting two cohorts:
cohort_of_interest <- subset(d35_input_full_copy, Cohort == cohort_selected)
time_col <- OUTCOME_COLS[["OS"]]$time
event_col <- OUTCOME_COLS[["OS"]]$event

cohort_of_interest <- cohort_of_interest[!is.na(cohort_of_interest[[parameter]]),]
fit <- survfit(
  Surv(cohort_of_interest[[time_col]], cohort_of_interest[[event_col]]) ~ cohort_of_interest[[parameter]],
  data = cohort_of_interest
)

# specificying the output directory for the KM plots:
save_dir <- ''
custom_palette = c("#4DAF4A", "#377EB8", "#A65628", "#E41A1C")
custom_legend_title = "Monosomal KT\nComplex KT\nChr5 del"
# custom_legend_title = "ELN Score + MDS-related\nadverse genetics"
# custom_legend_title = "ELN Score\nMDS/mon.KT"
# custom_legend_title = "ELN Score\nMonosomal KT"
# custom_legend_title = "ELN Score (3) &\nMonosomal KT"
# custom_legend_title = "Monosomal KT"
custom_legend_labs <- gsub("cohort_of_interest[[parameter]]=", "", names(fit$strata),fixed = TRUE)

param_name <- gsub(" ", "_", parameter, fixed = TRUE)
param_name <- gsub("|", "or", param_name, fixed = TRUE)
file_name = glue("{param_name}_{cohort_selected}_OS.pdf")

fontsize_small <- 8
fontsize_big <- 10

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
  risk.table.height = 0.3,
  # cumevents.height = 0.12,
  # cumcensor.height = 0.16,
  
  legend.title = custom_legend_title,
  legend.labs = custom_legend_labs,
  ggtheme = theme_bw(),
  fontsize = fontsize_small / ggplot2::.pt,
  break.time.by = break_time_by
  # censor.shape = 124,   # vertical tick
  # censor.size = 2       # smaller than default
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
  labs(x = x_axis_label)
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
    x = cutoff * 0.7,
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







###########################################################################
#### Kaplan-Meier curves: plotting the overall survival of two cohorts ####
###########################################################################
cohort_of_interest <- d35_input_full_copy
cohort_of_interest$Cohort <- ifelse(cohort_of_interest$Cohort == 'OC', 'Cohort 1', 'Cohort 2')
time_col <- OUTCOME_COLS[["OS"]]$time
event_col <- OUTCOME_COLS[["OS"]]$event
fit <- survfit(
  Surv(cohort_of_interest[[time_col]], cohort_of_interest[[event_col]]) ~ Cohort,
  data = cohort_of_interest
)

# specificying the output directory for the KM plots:
save_dir <- ''
custom_palette = c("#4DAF4A", "#377EB8", "#A65628", "#E41A1C")
custom_legend_title = "Patients"
custom_legend_labs <- gsub("Cohort=", "", names(fit$strata))
file_name = 'Cohort1&2_SimpleKM.pdf'


fontsize_small <- 8
fontsize_big <- 10

# ----- printing the plot (At Risk) -----
p <- ggsurvplot(
  fit, 
  data = cohort_of_interest, 
  size = 1,
  palette = custom_palette,
  conf.int = FALSE,
  pval = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.3,
  
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
  labs(x = x_axis_label)
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






########################################################################
#### Kaplan-Meier curves: splitting the cohort by final aGVHD grade ####
########################################################################
cohort_selected <- "Cohort2"
plot_type <- "AtRisk"
# 'AtRisk_Events_Censored' / 'AtRisk'

merging_pattern <- list(
  # "Grade 0" = c(0),
  # "Grade 1" = c(1),
  # "Grade 2" = c(2),
  # "Grades 3-4" = c(3, 4)
   
  "Grade 0" = c(0),
  "Grade 1-2" = c(1,2),
  "Grades 3-4" = c(3, 4)
  
  # "Grade 0-1-2" = c(0,1,2),
  # "Grades 3-4" = c(3, 4)
)

# plotting the KM curves:
# ---- ||| -----
lookup <- unlist(
  lapply(names(merging_pattern), function(name) {
    setNames(rep(name, length(merging_pattern[[name]])),
             merging_pattern[[name]])
  })
)
d35_input_full_copy$aGvHD_grouped <- lookup[as.character(d35_input_full_copy$aGVHD_grading_by_Katja)]

d35_input_full_copy_С1 <- d35_input_full_copy %>% subset(Cohort == 'OC')
d35_input_full_copy_С2 <- d35_input_full_copy %>% subset(Cohort == 'VC')

if (cohort_selected == 'Cohort1') {
  cohort_of_interest <- d35_input_full_copy_С1 
} else {
  cohort_of_interest <- d35_input_full_copy_С2 
}
time_col <- OUTCOME_COLS[["OS"]]$time
event_col <- OUTCOME_COLS[["OS"]]$event
fit <- survfit(
  Surv(cohort_of_interest[[time_col]], cohort_of_interest[[event_col]]) ~ aGvHD_grouped,
  data = cohort_of_interest
)

if (plot_type == 'AtRisk') {
  # specificying the output directory for the KM plots (with AtRisk count table):
  save_dir <- ''
} else {
  # specificying the output directory for the KM plots (with AtRisk, Events, Censored count tables):
  save_dir <- ''
}
custom_palette = c("#4DAF4A", "#377EB8", "#A65628", "#E41A1C")
custom_legend_title = "aGvHD grades"
custom_legend_labs <- gsub("aGvHD_grouped=", "", names(fit$strata))

group_string <- merging_pattern |>
  lapply(paste0, collapse = "") |>
  unlist(use.names = FALSE) |>
  paste(collapse = "_")

file_name = glue('aGvHD_grades_{group_string}_OS_{cohort_selected}_cutat{cutoff}days.pdf')

fontsize_small <- 8
fontsize_big <- 10

if (plot_type == 'AtRisk') {
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
    risk.table.height = 0.3,
    # cumevents.height = 0.12,
    # cumcensor.height = 0.16,
    
    legend.title = custom_legend_title,
    legend.labs = custom_legend_labs,
    ggtheme = theme_bw(),
    fontsize = fontsize_small / ggplot2::.pt,
    break.time.by = break_time_by
    # censor.shape = 124,   # vertical tick
    # censor.size = 2       # smaller than default
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
    labs(x = x_axis_label)
  # scale_x_continuous(
  #   breaks = c(0, 500, 1000, 1500, 1900),
  #   labels = c("0", "500", "1000", "1500","1900")
  # )
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
      x = cutoff * 0.7,
      y = 0.2,
      label = pval,
      size = fontsize_big / ggplot2::.pt,
      fill = "white",
      color = "black",
      label.r = unit(0, "lines"),
      linewidth = 0.4
    )
  # -----
} else {
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
    risk.table.height = 0.16,
    cumevents.height = 0.16,
    cumcensor.height = 0.2,
    
    legend.title = custom_legend_title,
    legend.labs = custom_legend_labs,
    ggtheme = theme_bw(),
    fontsize = fontsize_small / ggplot2::.pt,
    break.time.by = break_time_by
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
    labs(x = x_axis_label)
  
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
      legend.position = c(0.2, 0.25),
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
      x = cutoff * 0.7,
      y = 0.2,
      label = pval,
      size = fontsize_big / ggplot2::.pt,
      fill = "white",
      color = "black",
      label.r = unit(0, "lines"),
      linewidth = 0.4
    )
  # -----
}

pdf(glue('{save_dir}/{file_name}'), width = 6, height = 5)
print(p, newpage = FALSE)
dev.off()
# -----







########################################################################
#### Kaplan-Meier curves: splitting the cohort by final cGVHD grade ####
########################################################################
cohort_selected <- "Cohort2"
plot_type <- "AtRisk_Events_Censored"
# 'AtRisk_Events_Censored' / 'AtRisk'

merging_pattern <- list(
  "Grade 0" = c(0),
  "Grade 1" = c(1),
  "Grade 2" = c(2)
)

# ---- plotting the KM curves: ----
lookup <- unlist(
  lapply(names(merging_pattern), function(name) {
    setNames(rep(name, length(merging_pattern[[name]])),
             merging_pattern[[name]])
  })
)
d35_input_full_copy$cGvHD_grouped <- lookup[as.character(d35_input_full_copy$cGVHD_grading_by_Katja)]

d35_input_full_copy_С1 <- d35_input_full_copy %>% subset(Cohort == 'OC')
d35_input_full_copy_С2 <- d35_input_full_copy %>% subset(Cohort == 'VC')

if (cohort_selected == 'Cohort1') {
  cohort_of_interest <- d35_input_full_copy_С1 
} else {
  cohort_of_interest <- d35_input_full_copy_С2 
}
time_col <- OUTCOME_COLS[["OS"]]$time
event_col <- OUTCOME_COLS[["OS"]]$event
fit <- survfit(
  Surv(cohort_of_interest[[time_col]], cohort_of_interest[[event_col]]) ~ cGvHD_grouped,
  data = cohort_of_interest
)

if (plot_type == 'AtRisk') {
  # specificying the output directory for the KM plots (with AtRisk count table):
  save_dir <- ''
} else {
  # specificying the output directory for the KM plots (with AtRisk, Events, Censored count tables):
  save_dir <- ''
}
custom_palette = c("#4DAF4A", "#377EB8", "#A65628", "#E41A1C")
custom_legend_title = "cGvHD grades"
custom_legend_labs <- gsub("cGvHD_grouped=", "", names(fit$strata))

group_string <- merging_pattern |>
  lapply(paste0, collapse = "") |>
  unlist(use.names = FALSE) |>
  paste(collapse = "_")

file_name = glue('cGvHD_grades_{group_string}_OS_{cohort_selected}.pdf')


fontsize_small <- 8
fontsize_big <- 10
if (plot_type == 'AtRisk') {
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
    risk.table.height = 0.3,
    # cumevents.height = 0.12,
    # cumcensor.height = 0.16,
    
    legend.title = custom_legend_title,
    legend.labs = custom_legend_labs,
    ggtheme = theme_bw(),
    fontsize = fontsize_small / ggplot2::.pt,
    break.time.by = break_time_by
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
    labs(x = x_axis_label)
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
      x = cutoff * 0.7,
      y = 0.2,
      label = pval,
      size = fontsize_big / ggplot2::.pt,
      fill = "white",
      color = "black",
      label.r = unit(0, "lines"),
      linewidth = 0.4
    )
  # ----- 
} else {
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
    risk.table.height = 0.16,
    cumevents.height = 0.16,
    cumcensor.height = 0.2,
    
    legend.title = custom_legend_title,
    legend.labs = custom_legend_labs,
    ggtheme = theme_bw(),
    fontsize = fontsize_small / ggplot2::.pt,
    break.time.by = break_time_by
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
    labs(x = x_axis_label)
  
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
      x = cutoff * 0.7,
      y = 0.2,
      label = pval,
      size = fontsize_big / ggplot2::.pt,
      fill = "white",
      color = "black",
      label.r = unit(0, "lines"),
      linewidth = 0.4
    )
  # -----
}

pdf(glue('{save_dir}/{file_name}'), width = 6, height = 5)
print(p, newpage = FALSE)
dev.off()
# -----
