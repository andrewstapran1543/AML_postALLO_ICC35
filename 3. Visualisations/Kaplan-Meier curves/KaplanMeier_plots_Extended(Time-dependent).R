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
library(tidyr)
library(ggplot2)
library(cowplot)




##########################################################
#### Specifying horizoning cutoff and break intervals ####
##########################################################
cutoff <- 2000  # set your cutoff
landmarked <- TRUE
break_time_by <- 100




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
d35_input_full_copy_С1 <- d35_input_full_copy %>% subset(Cohort == 'OC')
d35_input_full_copy_С2 <- d35_input_full_copy %>% subset(Cohort == 'VC')

#############################################
#### cGVHD & aGVHD categories - from TPL ####
####   Simon Makuch curves for cGVHD     ####
#############################################

cohort_selected <- 'Cohort2'
# 'Cohort1' / 'Cohort2'
GVHD_type <- 'aGVHD'
# 'aGVHD' / 'cGVHD'


# ---- merging GVHD groups ----
if (GVHD_type == 'aGVHD') {
  merge_col <- 'aGVHD_grading_by_Katja'
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
} else {
  merge_col <- 'cGVHD_grading_by_Katja'
  merging_pattern <- list(
    "Grade 0" = c(0),
    "Grade 1" = c(1),
    "Grade 2" = c(2)
  )
}

lookup <- unlist(
  lapply(names(merging_pattern), function(name) {
    setNames(rep(name, length(merging_pattern[[name]])),
             merging_pattern[[name]])
  })
)
d35_input_full_copy_С1$GVHD_grouped <- lookup[as.character(d35_input_full_copy_С1[[merge_col]])]
d35_input_full_copy_С2$GVHD_grouped <- lookup[as.character(d35_input_full_copy_С2[[merge_col]])]
# -----


# ----- reshaping the dataframe -----

if (GVHD_type == 'cGVHD') {
  # ----- Cohot 1: -----
  tv_rows <- list()
  for(i in seq_len(nrow(d35_input_full_copy_С1))){
    
    patient <- d35_input_full_copy_С1[i, ]
    
    # No cGVHD or cGVHD after death/censoring
    if(
      is.na(patient$Time_to_cGVHD) |
      patient$cGVHD_grading_by_Katja == 0 |
      patient$Time_to_cGVHD >= patient$Time_to_Death_from_TPL
    ){
      
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = 0,
        stop = patient$Time_to_Death_from_TPL,
        event = patient$Death_Present_1_NotObserved_0,
        score = patient$GVHD_grouped,
        state = "cGVHD0",
        gvhd1 = 0,
        gvhd2 = 0
      )
      
    } else {
      
      # before cGVHD
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = 0,
        stop = patient$Time_to_cGVHD,
        event = 0,
        score = patient$GVHD_grouped,
        state = "cGVHD0",
        gvhd1 = 0,
        gvhd2 = 0
      )
      
      # after cGVHD
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = patient$Time_to_cGVHD,
        stop = patient$Time_to_Death_from_TPL,
        event = patient$Death_Present_1_NotObserved_0,
        score = patient$GVHD_grouped,
        state = ifelse(
          patient$cGVHD_grading_by_Katja == 1,
          "cGVHD1",
          "cGVHD2"
        ),
        gvhd1 = as.integer(patient$cGVHD_grading_by_Katja == 1),
        gvhd2 = as.integer(patient$cGVHD_grading_by_Katja == 2)
      )
    }
  }
  
  d35_input_full_copy_С1_tv_df <- bind_rows(tv_rows)
  
  d35_input_full_copy_С1_tv_df$state <- factor(
    d35_input_full_copy_С1_tv_df$state,
    levels = c("cGVHD0","cGVHD1","cGVHD2")
  )
  rownames(d35_input_full_copy_С1_tv_df) <- NULL
  # -----
  
  # ----- Cohort 2: -----
  tv_rows <- list()
  for(i in seq_len(nrow(d35_input_full_copy_С2))){
    
    patient <- d35_input_full_copy_С2[i, ]
    
    # No cGVHD or cGVHD after death/censoring
    if(
      is.na(patient$Time_to_cGVHD) |
      patient$cGVHD_grading_by_Katja == 0 |
      patient$Time_to_cGVHD >= patient$Time_to_Death_from_TPL
    ){
      
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = 0,
        stop = patient$Time_to_Death_from_TPL,
        event = patient$Death_Present_1_NotObserved_0,
        score = patient$GVHD_grouped,
        state = "cGVHD0",
        gvhd1 = 0,
        gvhd2 = 0
      )
      
    } else {
      
      # before cGVHD
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = 0,
        stop = patient$Time_to_cGVHD,
        event = 0,
        score = patient$GVHD_grouped,
        state = "cGVHD0",
        gvhd1 = 0,
        gvhd2 = 0
      )
      
      # after cGVHD
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = patient$Time_to_cGVHD,
        stop = patient$Time_to_Death_from_TPL,
        event = patient$Death_Present_1_NotObserved_0,
        score = patient$GVHD_grouped,
        state = ifelse(
          patient$cGVHD_grading_by_Katja == 1,
          "cGVHD1",
          "cGVHD2"
        ),
        gvhd1 = as.integer(patient$cGVHD_grading_by_Katja == 1),
        gvhd2 = as.integer(patient$cGVHD_grading_by_Katja == 2)
      )
    }
  }
  
  d35_input_full_copy_С2_tv_df <- bind_rows(tv_rows)
  
  d35_input_full_copy_С2_tv_df$state <- factor(
    d35_input_full_copy_С2_tv_df$state,
    levels = c("cGVHD0","cGVHD1","cGVHD2")
  )
  rownames(d35_input_full_copy_С2_tv_df) <- NULL
  # -----
} else {
  # ----- Cohot 1: -----
  tv_rows <- list()
  for(i in seq_len(nrow(d35_input_full_copy_С1))){
    
    patient <- d35_input_full_copy_С1[i, ]
    
    # No aGVHD or aGVHD after death/censoring
    if(
      is.na(patient$Time_to_aGVHD) |
      patient$aGVHD_grading_by_Katja == 0 |
      patient$Time_to_aGVHD >= patient$Time_to_Death_from_TPL
    ){
      
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = 0,
        stop = patient$Time_to_Death_from_TPL,
        event = patient$Death_Present_1_NotObserved_0,
        score = patient$GVHD_grouped,
        state = "aGVHD0",
        gvhd12 = 0,
        gvhd34 = 0
      )
      
    } else {
      
      # before aGVHD
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = 0,
        stop = patient$Time_to_aGVHD,
        event = 0,
        score = patient$GVHD_grouped,
        state = "aGVHD0",
        gvhd12 = 0,
        gvhd34 = 0
      )
      
      # after aGVHD
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = patient$Time_to_aGVHD,
        stop = patient$Time_to_Death_from_TPL,
        event = patient$Death_Present_1_NotObserved_0,
        score = patient$GVHD_grouped,
        state = ifelse(
          patient$aGVHD_grading_by_Katja == 1,
          "aGVHD12",
          "aGVHD34"
        ),
        gvhd12 = as.integer(patient$aGVHD_grading_by_Katja %in% c(1,2)),
        gvhd34 = as.integer(patient$aGVHD_grading_by_Katja %in% c(3,4))
      )
    }
  }
  
  d35_input_full_copy_С1_tv_df <- bind_rows(tv_rows)
  
  d35_input_full_copy_С1_tv_df$state <- factor(
    d35_input_full_copy_С1_tv_df$state,
    levels = c("aGVHD0","aGVHD12","aGVHD34")
  )
  rownames(d35_input_full_copy_С1_tv_df) <- NULL
  # -----
  
  # ----- Cohort 2: -----
  tv_rows <- list()
  for(i in seq_len(nrow(d35_input_full_copy_С2))){
    
    patient <- d35_input_full_copy_С2[i, ]
    
    # No aGVHD or aGVHD after death/censoring
    if(
      is.na(patient$Time_to_aGVHD) |
      patient$aGVHD_grading_by_Katja == 0 |
      patient$Time_to_aGVHD >= patient$Time_to_Death_from_TPL
    ){
      
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = 0,
        stop = patient$Time_to_Death_from_TPL,
        event = patient$Death_Present_1_NotObserved_0,
        score = patient$GVHD_grouped,
        state = "aGVHD0",
        gvhd12 = 0,
        gvhd34 = 0
      )
      
    } else {
      
      # before aGVHD
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = 0,
        stop = patient$Time_to_aGVHD,
        event = 0,
        score = patient$GVHD_grouped,
        state = "aGVHD0",
        gvhd12 = 0,
        gvhd34 = 0
      )
      
      # after aGVHD
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = patient$Time_to_aGVHD,
        stop = patient$Time_to_Death_from_TPL,
        event = patient$Death_Present_1_NotObserved_0,
        score = patient$GVHD_grouped,
        state = ifelse(
          patient$aGVHD_grading_by_Katja == 1,
          "aGVHD12",
          "aGVHD34"
        ),
        gvhd12 = as.integer(patient$aGVHD_grading_by_Katja %in% c(1,2)),
        gvhd34 = as.integer(patient$aGVHD_grading_by_Katja %in% c(3,4))
      )
    }
  }
  
  d35_input_full_copy_С2_tv_df <- bind_rows(tv_rows)
  
  d35_input_full_copy_С2_tv_df$state <- factor(
    d35_input_full_copy_С2_tv_df$state,
    levels = c("aGVHD0","aGVHD12","aGVHD34")
  )
  rownames(d35_input_full_copy_С2_tv_df) <- NULL
  # -----
}

if (cohort_selected == 'Cohort1') {
  main_df <- d35_input_full_copy_С1_tv_df
} else if (cohort_selected == 'Cohort2') {
  main_df <- d35_input_full_copy_С2_tv_df
}

fit <- survfit(Surv(start, stop, event) ~ state,data = main_df)

# specify the output directory for extended KM plots:
save_dir <- ''
custom_palette = c("#377EB8", "#4DAF4A",  "#A65628", "#E41A1C")

custom_legend_title = glue("{GVHD_type} states:")
custom_legend_labs <- gsub("state=", "", names(fit$strata))

group_string <- merging_pattern |>
  lapply(paste0, collapse = "") |>
  unlist(use.names = FALSE) |>
  paste(collapse = "_")
if (landmarked == TRUE) {
  landmarked_add_on <- '_landmarked_d35'
} else {
  landmarked_add_on <- '_no_landmarking'
}
file_name = glue('{cohort_selected}_{GVHD_type}_{group_string}_OS{landmarked_add_on}_cutat{cutoff}days.pdf')


fontsize_small <- 8
fontsize_big <- 10

width_selected <- 6
height_selected <- 5

# ---- main survival plot, without default risk table ----
p <- ggsurvplot(
  fit,
  data = main_df,
  size = 1,
  palette = custom_palette,
  conf.int = FALSE,
  pval = FALSE,
  risk.table = FALSE,
  legend.title = custom_legend_title,
  legend.labs = custom_legend_labs,
  ggtheme = theme_bw(),
  fontsize = fontsize_small / ggplot2::.pt,
  break.time.by = break_time_by
)

# ---- time-dependent Cox p-values ----
cox_state <- coxph(
  Surv(start, stop, event) ~ state,
  data = main_df
)

cox_sum <- summary(cox_state)

if (GVHD_type == 'cGVHD') {
  p_gvhd1 <- cox_sum$coefficients["statecGVHD1", "Pr(>|z|)"]
  p_gvhd2 <- cox_sum$coefficients["statecGVHD2", "Pr(>|z|)"] 
  label_states <- paste0(
    "Time-dependent Cox:\n",
    "cGvHD (1 vs 0): p = ", signif(p_gvhd1, 2), "\n",
    "cGVHD (2 vs 0): p = ", signif(p_gvhd2, 2)
  )
} else {
  p_gvhd12 <- cox_sum$coefficients["stateaGVHD12", "Pr(>|z|)"]
  p_gvhd34 <- cox_sum$coefficients["stateaGVHD34", "Pr(>|z|)"]
  label_states <- paste0(
    "Time-dependent Cox:\n",
    "aGvHD (1/2 vs 0): p = ", signif(p_gvhd12, 2), "\n",
    "aGVHD (3/4 vs 0): p = ", signif(p_gvhd34, 2)
  )
}


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
    legend.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.4
    ),
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
    x = cutoff * 0.5,
    y = 0.2,
    label = label_states,
    size = fontsize_big / ggplot2::.pt,
    fill = "white",
    color = "black",
    label.r = unit(0, "lines"),
    linewidth = 0.4,
    hjust = 0
  )

x_pad <- 25
x_limits_shared <- c(-x_pad, cutoff + x_pad)
time_breaks <- seq(0, cutoff, by = break_time_by)
x_scale_shared <- scale_x_continuous(
  breaks = time_breaks,
  limits = x_limits_shared,
  expand = expansion(mult = c(0, 0))
)
coord_surv <- coord_cartesian(
  xlim = x_limits_shared,
  expand = FALSE,
  clip = "on"
)

p$plot <- p$plot +
  x_scale_shared +
  coord_surv +
  scale_y_continuous(
    limits = c(-0.05, 1.05),
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0, 0))
  )
# ---- custom state-occupancy table ----
time_breaks <- seq(0, cutoff, by = break_time_by)
if (GVHD_type == 'cGVHD') {
  state_levels <- c("cGVHD0", "cGVHD1", "cGVHD2") 
} else {
  state_levels <- c("aGVHD0", "aGVHD12", "aGVHD34") 
}

risk_df <- lapply(time_breaks, function(t) {
  main_df %>%
    filter(start <= t, stop > t) %>%
    count(state) %>%
    complete(
      state = factor(state_levels, levels = state_levels),
      fill = list(n = 0)
    ) %>%
    mutate(time = t)
}) %>%
  bind_rows() %>%
  mutate(
    state = factor(state, levels = rev(state_levels))
  )

table_theme <- theme(
  panel.border = element_rect(
    colour = "grey60",
    fill = NA,
    linewidth = 0.3
  ),
  panel.grid.major = element_line(
    linewidth = 0.25,
    colour = "grey90"
  ),
  panel.grid.minor = element_blank()
)

coord_table <- coord_cartesian(
  xlim = x_limits_shared,
  clip = "on"
)
risk_table_plot <- ggplot(risk_df, aes(x = time, y = state, label = n)) +
  geom_text(size = fontsize_small / ggplot2::.pt) +
  x_scale_shared +
  coord_table + 
  scale_y_discrete(drop = FALSE) +
  labs(
    x = x_axis_label,
    y = "Patients in state"
  ) +
  theme_bw() +
  table_theme +
  theme(
    axis.text.y = element_text(size = fontsize_small),
    axis.text.x = element_text(size = fontsize_small),
    axis.title.x = element_text(size = fontsize_big, margin = margin(t = 2)),
    axis.title.y = element_text(size = fontsize_big, angle = 90, vjust = 0.5),
    axis.ticks.y = element_blank()
  )

# ---- combine plot and custom table ----
final_plot <- plot_grid(
  p$plot,
  risk_table_plot,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(0.72, 0.28)
)
# ---- save ----
pdf(glue('{save_dir}/{file_name}'), width = width_selected, height = height_selected)
print(final_plot, newpage = FALSE)
dev.off()
# -----







#########################################################
#### Score components vs cGVHD > landmarked analysis ####
####    Simon Makuch curves for cGVHD grades         ####
#########################################################

# ----- setting up medians & calculating the scores -----
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






##############################################################################
#### Score: CD8 TEM +/+ of CD8 TEM; CD4 TCM of T; CD4TCM/NV_to_rest ratio ####
####               Simon Makuch curves for cGVHD grades                   ####
##############################################################################

# ---- merging score groups ----
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
# -----


# ----- reshaping the dataframe -----
# Cohot 1:
tv_rows <- list()
for(i in seq_len(nrow(d35_input_full_copy_С1))){
  
  patient <- d35_input_full_copy_С1[i, ]
  
  # No cGVHD or cGVHD after death/censoring
  if(
    is.na(patient$Time_to_cGVHD) |
    patient$cGVHD_grading_by_Katja == 0 |
    patient$Time_to_cGVHD >= patient$Time_to_Death_from_TPL
  ){
    
    tv_rows[[length(tv_rows)+1]] <- data.frame(
      id = patient$ISH_ID_TPL_date_MERGED,
      start = 0,
      stop = patient$Time_to_Death_from_TPL,
      event = patient$Death_Present_1_NotObserved_0,
      score = patient$FinalScore_Grouped,
      state = "cGVHD0",
      gvhd1 = 0,
      gvhd2 = 0
    )
    
  } else {
    
    # before cGVHD
    tv_rows[[length(tv_rows)+1]] <- data.frame(
      id = patient$ISH_ID_TPL_date_MERGED,
      start = 0,
      stop = patient$Time_to_cGVHD,
      event = 0,
      score = patient$FinalScore_Grouped,
      state = "cGVHD0",
      gvhd1 = 0,
      gvhd2 = 0
    )
    
    # after cGVHD
    tv_rows[[length(tv_rows)+1]] <- data.frame(
      id = patient$ISH_ID_TPL_date_MERGED,
      start = patient$Time_to_cGVHD,
      stop = patient$Time_to_Death_from_TPL,
      event = patient$Death_Present_1_NotObserved_0,
      score = patient$FinalScore_Grouped,
      state = ifelse(
        patient$cGVHD_grading_by_Katja == 1,
        "cGVHD1",
        "cGVHD2"
      ),
      gvhd1 = as.integer(patient$cGVHD_grading_by_Katja == 1),
      gvhd2 = as.integer(patient$cGVHD_grading_by_Katja == 2)
    )
  }
}

d35_input_full_copy_С1_tv_df <- bind_rows(tv_rows)

d35_input_full_copy_С1_tv_df$state <- factor(
  d35_input_full_copy_С1_tv_df$state,
  levels = c("cGVHD0","cGVHD1","cGVHD2")
)
rownames(d35_input_full_copy_С1_tv_df) <- NULL


# Cohort 2:
tv_rows <- list()
for(i in seq_len(nrow(d35_input_full_copy_С2))){
  
  patient <- d35_input_full_copy_С2[i, ]
  
  # No cGVHD or cGVHD after death/censoring
  if(
    is.na(patient$Time_to_cGVHD) |
    patient$cGVHD_grading_by_Katja == 0 |
    patient$Time_to_cGVHD >= patient$Time_to_Death_from_TPL
  ){
    
    tv_rows[[length(tv_rows)+1]] <- data.frame(
      id = patient$ISH_ID_TPL_date_MERGED,
      start = 0,
      stop = patient$Time_to_Death_from_TPL,
      event = patient$Death_Present_1_NotObserved_0,
      score = patient$FinalScore_Grouped,
      state = "cGVHD0",
      gvhd1 = 0,
      gvhd2 = 0
    )
    
  } else {
    
    # before cGVHD
    tv_rows[[length(tv_rows)+1]] <- data.frame(
      id = patient$ISH_ID_TPL_date_MERGED,
      start = 0,
      stop = patient$Time_to_cGVHD,
      event = 0,
      score = patient$FinalScore_Grouped,
      state = "cGVHD0",
      gvhd1 = 0,
      gvhd2 = 0
    )
    
    # after cGVHD
    tv_rows[[length(tv_rows)+1]] <- data.frame(
      id = patient$ISH_ID_TPL_date_MERGED,
      start = patient$Time_to_cGVHD,
      stop = patient$Time_to_Death_from_TPL,
      event = patient$Death_Present_1_NotObserved_0,
      score = patient$FinalScore_Grouped,
      state = ifelse(
        patient$cGVHD_grading_by_Katja == 1,
        "cGVHD1",
        "cGVHD2"
      ),
      gvhd1 = as.integer(patient$cGVHD_grading_by_Katja == 1),
      gvhd2 = as.integer(patient$cGVHD_grading_by_Katja == 2)
    )
  }
}

d35_input_full_copy_С2_tv_df <- bind_rows(tv_rows)

d35_input_full_copy_С2_tv_df$state <- factor(
  d35_input_full_copy_С2_tv_df$state,
  levels = c("cGVHD0","cGVHD1","cGVHD2")
)
rownames(d35_input_full_copy_С2_tv_df) <- NULL
# -----

score_group_selected <- c('GvLhi_GvHDlo', 'Other')
# c('GvLhi_GvHDlo') / c('Other') / c('GvLhi_GvHDlo', 'Other')
d35_input_full_copy_С1_tv_df_SELECTED <- subset(d35_input_full_copy_С1_tv_df, score %in% score_group_selected) 
d35_input_full_copy_С2_tv_df_SELECTED <- subset(d35_input_full_copy_С2_tv_df, score %in% score_group_selected)

cohort_selected <- 'Cohort2'
# 'Cohort1' / 'Cohort2'
if (cohort_selected == 'Cohort1') {
  main_df <- d35_input_full_copy_С1_tv_df_SELECTED
} else if (cohort_selected == 'Cohort2') {
  main_df <- d35_input_full_copy_С2_tv_df_SELECTED
}

fit <- survfit(Surv(start, stop, event) ~ state,data = main_df)
# ggsurvplot(fit,data = main_df,risk.table = TRUE,conf.int = FALSE)

# specify the output directory for extendede KM plots:
save_dir <- ''
custom_palette = c("#377EB8", "#4DAF4A",  "#A65628", "#E41A1C")

custom_legend_title = "cGvHD states:"
custom_legend_labs <- gsub("state=", "", names(fit$strata))
score_group_selected_name <- paste(score_group_selected, collapse = '*')
file_name = glue('{cohort_selected}_CD8TEMCD4TCMCD4NVTCMRATIO_FinalScoreGroup{score_group_selected_name}.pdf')


fontsize_small <- 8
fontsize_big <- 10

width_selected <- 6
height_selected <- 5


# ---- main survival plot, without default risk table ----
p <- ggsurvplot(
  fit,
  data = main_df,
  size = 1,
  palette = custom_palette,
  conf.int = FALSE,
  pval = FALSE,
  risk.table = FALSE,
  legend.title = custom_legend_title,
  legend.labs = custom_legend_labs,
  ggtheme = theme_bw(),
  fontsize = fontsize_small / ggplot2::.pt,
  break.time.by = 500
)

# ---- time-dependent Cox p-values ----
cox_state <- coxph(
  Surv(start, stop, event) ~ state,
  data = main_df
)

cox_sum <- summary(cox_state)

p_gvhd1 <- cox_sum$coefficients["statecGVHD1", "Pr(>|z|)"]
p_gvhd2 <- cox_sum$coefficients["statecGVHD2", "Pr(>|z|)"]

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
    legend.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.4
    ),
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
    x = 1000,
    y = 0.2,
    label = paste0(
      "Time-dependent Cox:\n",
      "cGvHD (1 vs 0): p = ", signif(p_gvhd1, 2), "\n",
      "cGVHD (2 vs 0): p = ", signif(p_gvhd2, 2)
    ),
    size = fontsize_big / ggplot2::.pt,
    fill = "white",
    color = "black",
    label.r = unit(0, "lines"),
    linewidth = 0.4,
    hjust = 0
  )

x_limits <- c(0, 2000)
x_scale_shared <- scale_x_continuous(
  breaks = time_breaks,
  limits = x_limits,
  expand = expansion(mult = c(0.04, 0.04))
)
p$plot <- p$plot + x_scale_shared

# ---- custom state-occupancy table ----
time_breaks <- seq(0, 2000, by = 500)
state_levels <- c("cGVHD0", "cGVHD1", "cGVHD2")

risk_df <- lapply(time_breaks, function(t) {
  main_df %>%
    filter(start <= t, stop > t) %>%
    count(state) %>%
    complete(
      state = factor(state_levels, levels = state_levels),
      fill = list(n = 0)
    ) %>%
    mutate(time = t)
}) %>%
  bind_rows() %>%
  mutate(
    state = factor(state, levels = rev(state_levels))
  )

table_theme <- theme(
  panel.border = element_rect(
    colour = "grey60",
    fill = NA,
    linewidth = 0.3
  ),
  panel.grid.major = element_line(
    linewidth = 0.25,
    colour = "grey90"
  ),
  panel.grid.minor = element_blank()
)

risk_table_plot <- ggplot(risk_df, aes(x = time, y = state, label = n)) +
  geom_text(size = fontsize_small / ggplot2::.pt) +
  x_scale_shared +
  scale_y_discrete(drop = FALSE) +
  labs(
    x = "Days post alloHCT sampling",
    y = "Patients in state"
  ) +
  theme_bw() +
  table_theme +
  theme(
    axis.text.y = element_text(size = fontsize_small),
    axis.text.x = element_text(size = fontsize_small),
    axis.title.x = element_text(size = fontsize_big, margin = margin(t = 2)),
    axis.title.y = element_text(size = fontsize_big, angle = 90, vjust = 0.5),
    axis.ticks.y = element_blank()
  )

# ---- combine plot and custom table ----
final_plot <- plot_grid(
  p$plot,
  risk_table_plot,
  ncol = 1,
  align = "v",
  axis = "tblr",
  rel_heights = c(0.72, 0.28)
)

# ---- save ----
pdf(glue('{save_dir}/{file_name}'), width = width_selected, height = height_selected)
print(final_plot, newpage = FALSE)
dev.off()
# -----










##############################################
####   CD8 TEM +/+ of CD8 TEM: Hi vs Lo   ####
#### Simon Makuch curves for cGVHD grades ####
##############################################

cGVHD012_grouping <- '0_12'
# '0_1_2' / '0_12'

if (cGVHD012_grouping == '0_1_2') {
  # ----- reshaping the dataframe -----
  # Cohot 1:
  tv_rows <- list()
  for(i in seq_len(nrow(d35_input_full_copy_С1))){
    
    patient <- d35_input_full_copy_С1[i, ]
    
    # No cGVHD or cGVHD after death/censoring
    if(
      is.na(patient$Time_to_cGVHD) |
      patient$cGVHD_grading_by_Katja == 0 |
      patient$Time_to_cGVHD >= patient$Time_to_Death_from_TPL
    ){
      
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = 0,
        stop = patient$Time_to_Death_from_TPL,
        event = patient$Death_Present_1_NotObserved_0,
        score = patient$CD8_score,
        state = "cGVHD0",
        gvhd1 = 0,
        gvhd2 = 0
      )
      
    } else {
      
      # before cGVHD
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = 0,
        stop = patient$Time_to_cGVHD,
        event = 0,
        score = patient$CD8_score,
        state = "cGVHD0",
        gvhd1 = 0,
        gvhd2 = 0
      )
      
      # after cGVHD
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = patient$Time_to_cGVHD,
        stop = patient$Time_to_Death_from_TPL,
        event = patient$Death_Present_1_NotObserved_0,
        score = patient$CD8_score,
        state = ifelse(
          patient$cGVHD_grading_by_Katja == 1,
          "cGVHD1",
          "cGVHD2"
        ),
        gvhd1 = as.integer(patient$cGVHD_grading_by_Katja == 1),
        gvhd2 = as.integer(patient$cGVHD_grading_by_Katja == 2)
      )
    }
  }
  
  d35_input_full_copy_С1_tv_df <- bind_rows(tv_rows)
  
  d35_input_full_copy_С1_tv_df$state <- factor(
    d35_input_full_copy_С1_tv_df$state,
    levels = c("cGVHD0","cGVHD1","cGVHD2")
  )
  rownames(d35_input_full_copy_С1_tv_df) <- NULL
  
  
  # Cohort 2:
  tv_rows <- list()
  for(i in seq_len(nrow(d35_input_full_copy_С2))){
    
    patient <- d35_input_full_copy_С2[i, ]
    
    # No cGVHD or cGVHD after death/censoring
    if(
      is.na(patient$Time_to_cGVHD) |
      patient$cGVHD_grading_by_Katja == 0 |
      patient$Time_to_cGVHD >= patient$Time_to_Death_from_TPL
    ){
      
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = 0,
        stop = patient$Time_to_Death_from_TPL,
        event = patient$Death_Present_1_NotObserved_0,
        score = patient$CD8_score,
        state = "cGVHD0",
        gvhd1 = 0,
        gvhd2 = 0
      )
      
    } else {
      
      # before cGVHD
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = 0,
        stop = patient$Time_to_cGVHD,
        event = 0,
        score = patient$CD8_score,
        state = "cGVHD0",
        gvhd1 = 0,
        gvhd2 = 0
      )
      
      # after cGVHD
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = patient$Time_to_cGVHD,
        stop = patient$Time_to_Death_from_TPL,
        event = patient$Death_Present_1_NotObserved_0,
        score = patient$CD8_score,
        state = ifelse(
          patient$cGVHD_grading_by_Katja == 1,
          "cGVHD1",
          "cGVHD2"
        ),
        gvhd1 = as.integer(patient$cGVHD_grading_by_Katja == 1),
        gvhd2 = as.integer(patient$cGVHD_grading_by_Katja == 2)
      )
    }
  }
  
  d35_input_full_copy_С2_tv_df <- bind_rows(tv_rows)
  
  d35_input_full_copy_С2_tv_df$state <- factor(
    d35_input_full_copy_С2_tv_df$state,
    levels = c("cGVHD0","cGVHD1","cGVHD2")
  )
  rownames(d35_input_full_copy_С2_tv_df) <- NULL
  # -----
} else if (cGVHD012_grouping == '0_12') {
  # ----- reshaping the dataframe (cGVHD0 vs cGVHD1/2) -----
  # Cohot 1:
  tv_rows <- list()
  for(i in seq_len(nrow(d35_input_full_copy_С1))){
    
    patient <- d35_input_full_copy_С1[i, ]
    
    # No cGVHD or cGVHD after death/censoring
    if(
      is.na(patient$Time_to_cGVHD) |
      patient$cGVHD_grading_by_Katja == 0 |
      patient$Time_to_cGVHD >= patient$Time_to_Death_from_TPL
    ){
      
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = 0,
        stop = patient$Time_to_Death_from_TPL,
        event = patient$Death_Present_1_NotObserved_0,
        score = patient$CD8_score,
        state = "cGVHD0",
        gvhd12 = 0
      )
      
    } else {
      
      # before cGVHD
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = 0,
        stop = patient$Time_to_cGVHD,
        event = 0,
        score = patient$CD8_score,
        state = "cGVHD0",
        gvhd12 = 0
      )
      
      # after cGVHD
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = patient$Time_to_cGVHD,
        stop = patient$Time_to_Death_from_TPL,
        event = patient$Death_Present_1_NotObserved_0,
        score = patient$CD8_score,
        state = "cGVHD1/2",
        gvhd12 = 1
      )
    }
  }
  
  d35_input_full_copy_С1_tv_df <- bind_rows(tv_rows)
  
  d35_input_full_copy_С1_tv_df$state <- factor(
    d35_input_full_copy_С1_tv_df$state,
    levels = c("cGVHD0","cGVHD1/2")
  )
  rownames(d35_input_full_copy_С1_tv_df) <- NULL
  
  
  # Cohort 2:
  tv_rows <- list()
  for(i in seq_len(nrow(d35_input_full_copy_С2))){
    
    patient <- d35_input_full_copy_С2[i, ]
    
    # No cGVHD or cGVHD after death/censoring
    if(
      is.na(patient$Time_to_cGVHD) |
      patient$cGVHD_grading_by_Katja == 0 |
      patient$Time_to_cGVHD >= patient$Time_to_Death_from_TPL
    ){
      
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = 0,
        stop = patient$Time_to_Death_from_TPL,
        event = patient$Death_Present_1_NotObserved_0,
        score = patient$CD8_score,
        state = "cGVHD0",
        gvhd12 = 0
      )
      
    } else {
      
      # before cGVHD
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = 0,
        stop = patient$Time_to_cGVHD,
        event = 0,
        score = patient$CD8_score,
        state = "cGVHD0",
        gvhd12 = 0
      )
      
      # after cGVHD
      tv_rows[[length(tv_rows)+1]] <- data.frame(
        id = patient$ISH_ID_TPL_date_MERGED,
        start = patient$Time_to_cGVHD,
        stop = patient$Time_to_Death_from_TPL,
        event = patient$Death_Present_1_NotObserved_0,
        score = patient$CD8_score,
        state = "cGVHD1/2",
        gvhd12 = 1
      )
    }
  }
  
  d35_input_full_copy_С2_tv_df <- bind_rows(tv_rows)
  
  d35_input_full_copy_С2_tv_df$state <- factor(
    d35_input_full_copy_С2_tv_df$state,
    levels = c("cGVHD0","cGVHD1/2")
  )
  rownames(d35_input_full_copy_С2_tv_df) <- NULL
  # -----
}

score_group_selected <- c('CD8lo')
# c('CD8hi') / c('CD8lo') / c('CD8hi', 'CD8lo')
d35_input_full_copy_С1_tv_df_SELECTED <- subset(d35_input_full_copy_С1_tv_df, score %in% score_group_selected) 
d35_input_full_copy_С2_tv_df_SELECTED <- subset(d35_input_full_copy_С2_tv_df, score %in% score_group_selected)

cohort_selected <- 'Cohort2'
# 'Cohort1' / 'Cohort2'
if (cohort_selected == 'Cohort1') {
  main_df <- d35_input_full_copy_С1_tv_df_SELECTED
} else if (cohort_selected == 'Cohort2') {
  main_df <- d35_input_full_copy_С2_tv_df_SELECTED
}

fit <- survfit(Surv(start, stop, event) ~ state,data = main_df)
# ggsurvplot(fit,data = main_df,risk.table = TRUE,conf.int = FALSE)

# specify the output directory for extendede KM plots:
save_dir <- ''
custom_palette = c("#377EB8", "#4DAF4A",  "#A65628", "#E41A1C")

custom_legend_title = "cGvHD states:"
custom_legend_labs <- gsub("state=", "", names(fit$strata))
score_group_selected_name <- paste(score_group_selected, collapse = '*')
file_name = glue('{cohort_selected}_CD8TEMparam_{score_group_selected_name}_cGVHD{cGVHD012_grouping}_grouped.pdf')


fontsize_small <- 8
fontsize_big <- 10

width_selected <- 6
height_selected <- 5


# ---- main survival plot, without default risk table ----
p <- ggsurvplot(
  fit,
  data = main_df,
  size = 1,
  palette = custom_palette,
  conf.int = FALSE,
  pval = FALSE,
  risk.table = FALSE,
  legend.title = custom_legend_title,
  legend.labs = custom_legend_labs,
  ggtheme = theme_bw(),
  fontsize = fontsize_small / ggplot2::.pt,
  break.time.by = 500
)

# ---- time-dependent Cox p-values ----
cox_state <- coxph(
  Surv(start, stop, event) ~ state,
  data = main_df
)

cox_sum <- summary(cox_state)
if (cGVHD012_grouping == '0_1_2') {
  p_gvhd1 <- cox_sum$coefficients["statecGVHD1", "Pr(>|z|)"]
  p_gvhd2 <- cox_sum$coefficients["statecGVHD2", "Pr(>|z|)"] 
  label_plot <- paste0(
    "Time-dependent Cox:\n",
    "cGvHD (1 vs 0): p = ", signif(p_gvhd1, 2), "\n",
    "cGVHD (2 vs 0): p = ", signif(p_gvhd2, 2)
  )
  state_levels <- c("cGVHD0", "cGVHD1", "cGVHD2")
} else if (cGVHD012_grouping == '0_12') {
  p_gvhd12 <- cox_sum$coefficients["statecGVHD1/2", "Pr(>|z|)"]
  label_plot <- paste0(
    "Time-dependent Cox:\n",
    "cGvHD (1/2 vs 0): p = ", signif(p_gvhd12, 2)
  )
  state_levels <- c("cGVHD0", "cGVHD1/2")
}

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
    legend.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.4
    ),
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
    x = 1000,
    y = 0.2,
    label = label_plot,
    size = fontsize_big / ggplot2::.pt,
    fill = "white",
    color = "black",
    label.r = unit(0, "lines"),
    linewidth = 0.4,
    hjust = 0
  )

x_limits <- c(0, 2000)
x_scale_shared <- scale_x_continuous(
  breaks = time_breaks,
  limits = x_limits,
  expand = expansion(mult = c(0.04, 0.04))
)
p$plot <- p$plot + x_scale_shared

# ---- custom state-occupancy table ----
time_breaks <- seq(0, 2000, by = 500)
# state_levels <- c("cGVHD0", "cGVHD1", "cGVHD2")

risk_df <- lapply(time_breaks, function(t) {
  main_df %>%
    filter(start <= t, stop > t) %>%
    count(state) %>%
    complete(
      state = factor(state_levels, levels = state_levels),
      fill = list(n = 0)
    ) %>%
    mutate(time = t)
}) %>%
  bind_rows() %>%
  mutate(
    state = factor(state, levels = rev(state_levels))
  )

table_theme <- theme(
  panel.border = element_rect(
    colour = "grey60",
    fill = NA,
    linewidth = 0.3
  ),
  panel.grid.major = element_line(
    linewidth = 0.25,
    colour = "grey90"
  ),
  panel.grid.minor = element_blank()
)

risk_table_plot <- ggplot(risk_df, aes(x = time, y = state, label = n)) +
  geom_text(size = fontsize_small / ggplot2::.pt) +
  x_scale_shared +
  scale_y_discrete(drop = FALSE) +
  labs(
    x = "Days post alloHCT sampling",
    y = "Patients in state"
  ) +
  theme_bw() +
  table_theme +
  theme(
    axis.text.y = element_text(size = fontsize_small),
    axis.text.x = element_text(size = fontsize_small),
    axis.title.x = element_text(size = fontsize_big, margin = margin(t = 2)),
    axis.title.y = element_text(size = fontsize_big, angle = 90, vjust = 0.5),
    axis.ticks.y = element_blank()
  )

# ---- combine plot and custom table ----
final_plot <- plot_grid(
  p$plot,
  risk_table_plot,
  ncol = 1,
  align = "v",
  axis = "tblr",
  rel_heights = c(0.72, 0.28)
)

# ---- save ----
pdf(glue('{save_dir}/{file_name}'), width = width_selected, height = height_selected)
print(final_plot, newpage = FALSE)
dev.off()
# -----





##############################################
####   CD4 parameters: Hi vs Lo           ####
#### Simon Makuch curves for cGVHD grades ####
##############################################

criteria_fulfilled <- c(2)
d35_input_full_copy_С1$CD4ScoreOverall <- d35_input_full_copy_С1$CD4TCM_score + d35_input_full_copy_С1$CD4TCMNVratio_score
d35_input_full_copy_С1$CD4ScoreOverall <- ifelse(d35_input_full_copy_С1$CD4ScoreOverall %in% criteria_fulfilled, 'CD4ParamLo', 'CD4ParamHi')
d35_input_full_copy_С2$CD4ScoreOverall <- d35_input_full_copy_С2$CD4TCM_score + d35_input_full_copy_С2$CD4TCMNVratio_score
d35_input_full_copy_С2$CD4ScoreOverall <- ifelse(d35_input_full_copy_С2$CD4ScoreOverall %in% criteria_fulfilled, 'CD4ParamLo', 'CD4ParamHi')

# ----- reshaping the dataframe -----
# Cohot 1:
tv_rows <- list()
for(i in seq_len(nrow(d35_input_full_copy_С1))){
  
  patient <- d35_input_full_copy_С1[i, ]
  
  # No cGVHD or cGVHD after death/censoring
  if(
    is.na(patient$Time_to_cGVHD) |
    patient$cGVHD_grading_by_Katja == 0 |
    patient$Time_to_cGVHD >= patient$Time_to_Death_from_TPL
  ){
    
    tv_rows[[length(tv_rows)+1]] <- data.frame(
      id = patient$ISH_ID_TPL_date_MERGED,
      start = 0,
      stop = patient$Time_to_Death_from_TPL,
      event = patient$Death_Present_1_NotObserved_0,
      score = patient$CD4ScoreOverall,
      state = "cGVHD0",
      gvhd1 = 0,
      gvhd2 = 0
    )
    
  } else {
    
    # before cGVHD
    tv_rows[[length(tv_rows)+1]] <- data.frame(
      id = patient$ISH_ID_TPL_date_MERGED,
      start = 0,
      stop = patient$Time_to_cGVHD,
      event = 0,
      score = patient$CD4ScoreOverall,
      state = "cGVHD0",
      gvhd1 = 0,
      gvhd2 = 0
    )
    
    # after cGVHD
    tv_rows[[length(tv_rows)+1]] <- data.frame(
      id = patient$ISH_ID_TPL_date_MERGED,
      start = patient$Time_to_cGVHD,
      stop = patient$Time_to_Death_from_TPL,
      event = patient$Death_Present_1_NotObserved_0,
      score = patient$CD4ScoreOverall,
      state = ifelse(
        patient$cGVHD_grading_by_Katja == 1,
        "cGVHD1",
        "cGVHD2"
      ),
      gvhd1 = as.integer(patient$cGVHD_grading_by_Katja == 1),
      gvhd2 = as.integer(patient$cGVHD_grading_by_Katja == 2)
    )
  }
}

d35_input_full_copy_С1_tv_df <- bind_rows(tv_rows)

d35_input_full_copy_С1_tv_df$state <- factor(
  d35_input_full_copy_С1_tv_df$state,
  levels = c("cGVHD0","cGVHD1","cGVHD2")
)
rownames(d35_input_full_copy_С1_tv_df) <- NULL


# Cohort 2:
tv_rows <- list()
for(i in seq_len(nrow(d35_input_full_copy_С2))){
  
  patient <- d35_input_full_copy_С2[i, ]
  
  # No cGVHD or cGVHD after death/censoring
  if(
    is.na(patient$Time_to_cGVHD) |
    patient$cGVHD_grading_by_Katja == 0 |
    patient$Time_to_cGVHD >= patient$Time_to_Death_from_TPL
  ){
    
    tv_rows[[length(tv_rows)+1]] <- data.frame(
      id = patient$ISH_ID_TPL_date_MERGED,
      start = 0,
      stop = patient$Time_to_Death_from_TPL,
      event = patient$Death_Present_1_NotObserved_0,
      score = patient$CD4ScoreOverall,
      state = "cGVHD0",
      gvhd1 = 0,
      gvhd2 = 0
    )
    
  } else {
    
    # before cGVHD
    tv_rows[[length(tv_rows)+1]] <- data.frame(
      id = patient$ISH_ID_TPL_date_MERGED,
      start = 0,
      stop = patient$Time_to_cGVHD,
      event = 0,
      score = patient$CD4ScoreOverall,
      state = "cGVHD0",
      gvhd1 = 0,
      gvhd2 = 0
    )
    
    # after cGVHD
    tv_rows[[length(tv_rows)+1]] <- data.frame(
      id = patient$ISH_ID_TPL_date_MERGED,
      start = patient$Time_to_cGVHD,
      stop = patient$Time_to_Death_from_TPL,
      event = patient$Death_Present_1_NotObserved_0,
      score = patient$CD4ScoreOverall,
      state = ifelse(
        patient$cGVHD_grading_by_Katja == 1,
        "cGVHD1",
        "cGVHD2"
      ),
      gvhd1 = as.integer(patient$cGVHD_grading_by_Katja == 1),
      gvhd2 = as.integer(patient$cGVHD_grading_by_Katja == 2)
    )
  }
}

d35_input_full_copy_С2_tv_df <- bind_rows(tv_rows)

d35_input_full_copy_С2_tv_df$state <- factor(
  d35_input_full_copy_С2_tv_df$state,
  levels = c("cGVHD0","cGVHD1","cGVHD2")
)
rownames(d35_input_full_copy_С2_tv_df) <- NULL
# -----

score_group_selected <- c('CD4ParamLo')
# c('CD4ParamHi') / c('CD4ParamLo') / c('CD4ParamHi', 'CD4ParamLo')
d35_input_full_copy_С1_tv_df_SELECTED <- subset(d35_input_full_copy_С1_tv_df, score %in% score_group_selected) 
d35_input_full_copy_С2_tv_df_SELECTED <- subset(d35_input_full_copy_С2_tv_df, score %in% score_group_selected)

cohort_selected <- 'Cohort2'
# 'Cohort1' / 'Cohort2'
if (cohort_selected == 'Cohort1') {
  main_df <- d35_input_full_copy_С1_tv_df_SELECTED
} else if (cohort_selected == 'Cohort2') {
  main_df <- d35_input_full_copy_С2_tv_df_SELECTED
}

fit <- survfit(Surv(start, stop, event) ~ state,data = main_df)
# ggsurvplot(fit,data = main_df,risk.table = TRUE,conf.int = FALSE)

# specify the output directory for extendede KM plots:
save_dir <- ''
custom_palette = c("#377EB8", "#4DAF4A",  "#A65628", "#E41A1C")

custom_legend_title = "cGvHD states:"
custom_legend_labs <- gsub("state=", "", names(fit$strata))
score_group_selected_name <- paste(score_group_selected, collapse = '*')
CD4_crit_fulfilled <- paste(criteria_fulfilled, collapse = '*')
file_name = glue('{cohort_selected}_CD4param_{score_group_selected_name}_critfulfilled_{CD4_crit_fulfilled}.pdf')


fontsize_small <- 8
fontsize_big <- 10

width_selected <- 6
height_selected <- 5


# ---- main survival plot, without default risk table ----
p <- ggsurvplot(
  fit,
  data = main_df,
  size = 1,
  palette = custom_palette,
  conf.int = FALSE,
  pval = FALSE,
  risk.table = FALSE,
  legend.title = custom_legend_title,
  legend.labs = custom_legend_labs,
  ggtheme = theme_bw(),
  fontsize = fontsize_small / ggplot2::.pt,
  break.time.by = 500
)

# ---- time-dependent Cox p-values ----
cox_state <- coxph(
  Surv(start, stop, event) ~ state,
  data = main_df
)

cox_sum <- summary(cox_state)

p_gvhd1 <- cox_sum$coefficients["statecGVHD1", "Pr(>|z|)"]
p_gvhd2 <- cox_sum$coefficients["statecGVHD2", "Pr(>|z|)"]

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
    legend.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.4
    ),
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
    x = 1000,
    y = 0.2,
    label = paste0(
      "Time-dependent Cox:\n",
      "cGvHD (1 vs 0): p = ", signif(p_gvhd1, 2), "\n",
      "cGVHD (2 vs 0): p = ", signif(p_gvhd2, 2)
    ),
    size = fontsize_big / ggplot2::.pt,
    fill = "white",
    color = "black",
    label.r = unit(0, "lines"),
    linewidth = 0.4,
    hjust = 0
  )

x_limits <- c(0, 2000)
x_scale_shared <- scale_x_continuous(
  breaks = time_breaks,
  limits = x_limits,
  expand = expansion(mult = c(0.04, 0.04))
)
p$plot <- p$plot + x_scale_shared

# ---- custom state-occupancy table ----
time_breaks <- seq(0, 2000, by = 500)
state_levels <- c("cGVHD0", "cGVHD1", "cGVHD2")

risk_df <- lapply(time_breaks, function(t) {
  main_df %>%
    filter(start <= t, stop > t) %>%
    count(state) %>%
    complete(
      state = factor(state_levels, levels = state_levels),
      fill = list(n = 0)
    ) %>%
    mutate(time = t)
}) %>%
  bind_rows() %>%
  mutate(
    state = factor(state, levels = rev(state_levels))
  )

table_theme <- theme(
  panel.border = element_rect(
    colour = "grey60",
    fill = NA,
    linewidth = 0.3
  ),
  panel.grid.major = element_line(
    linewidth = 0.25,
    colour = "grey90"
  ),
  panel.grid.minor = element_blank()
)

risk_table_plot <- ggplot(risk_df, aes(x = time, y = state, label = n)) +
  geom_text(size = fontsize_small / ggplot2::.pt) +
  x_scale_shared +
  scale_y_discrete(drop = FALSE) +
  labs(
    x = "Days post alloHCT sampling",
    y = "Patients in state"
  ) +
  theme_bw() +
  table_theme +
  theme(
    axis.text.y = element_text(size = fontsize_small),
    axis.text.x = element_text(size = fontsize_small),
    axis.title.x = element_text(size = fontsize_big, margin = margin(t = 2)),
    axis.title.y = element_text(size = fontsize_big, angle = 90, vjust = 0.5),
    axis.ticks.y = element_blank()
  )

# ---- combine plot and custom table ----
final_plot <- plot_grid(
  p$plot,
  risk_table_plot,
  ncol = 1,
  align = "v",
  axis = "tblr",
  rel_heights = c(0.72, 0.28)
)

# ---- save ----
pdf(glue('{save_dir}/{file_name}'), width = width_selected, height = height_selected)
print(final_plot, newpage = FALSE)
dev.off()
# -----











