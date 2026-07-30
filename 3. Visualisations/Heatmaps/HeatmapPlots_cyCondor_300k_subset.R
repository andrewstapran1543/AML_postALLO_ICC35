# loading the packages:
library(cyCONDOR)
library(qs)
library(glue)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggridges)
library(ggrastr)
library(Hmisc)
library(Polychrome)
library(ComplexHeatmap)
library(circlize)


# ----- d35 dataset -----
path_to_the_cyCondor_dataset_VISsubset_full_anntations <- ''
main_dataset_d35 <- qread(path_to_the_cyCondor_dataset_VISsubset_full_anntations)
# -----

scale <- TRUE
cohort_to_vis = 'Cohort 2, day 35'
# 'both' / 'Cohort 1, day 35' / 'Cohort 2, day 35'
chosen_color_scale <- "standard"
# 'green_col' / 'purple_col' / 'standard'
exclusion_markers_list <- c("FSC-A","FSC-H","FSC-W","SSC-A","SSC-H","SSC-W")
heatmap_design <- 'pops_rows_&_markers_cols'
# 'pops_rows_&_markers_cols' / 'markers_rows_&_pops_cols'
population_split <- 'NoSplitInAnnotation'
# 'SplitInAnnotation' / 'NoSplitInAnnotation'
visualize_details <- 'CD8_CD4_full_names'
# 'GeneralAnnotation' / 'CD8_CD4_more_details' / 'CD8_CD4_full_names'
sub_superscripts <- 'PlainAIVersion'
# 'AlreadyInR' / 'PlainAIVersion'

if (visualize_details == 'GeneralAnnotation') {
  vis_add_on <- 'GenAnnot'
} else if (visualize_details == 'CD8_CD4_more_details') {
  vis_add_on <- 'CD8CD4details'
} else if (visualize_details == 'CD8_CD4_full_names') {
  vis_add_on <- 'CD8CD4detailsfull'
}

# ----- d35 dict -----
if (visualize_details == 'GeneralAnnotation') {
  parent1_population_dict_d35 <- c(
    "ALL_res45_14"="Exclude",
    "ALL_res45_5"="Exclude",
    "ALL_res45_45"="B",
    "ALL_res45_37"="CD56bright NK",
    "ALL_res45_38"="CD56bright NK",
    "ALL_res45_44"="CD56bright NK",
    "ALL_res45_40"="CD56bright NK",
    "ALL_res45_43"="CD56bright NK",
    "ALL_res55_47"="CD56bright NK",
    "ALL_res55_52"="CD56bright NK",
    "ALL_res45_34"="CD56dim NK",
    "CD3pos_res60_16"="CD56dim NK",
    "CD3pos_res60_29"="CD56dim NK",
    "CD3pos_res60_37"="CD56dim NK",
    "CD3pos_res60_54"="CD56dim NK",
    "CD3pos_res70_40"="CD56dim NK",
    "ALL_res45_35"="CD56dim NK",
    "ALL_res55_42"="CD56dim NK",
    "ALL_res55_49"="CD56dim NK",
    "CD3pos_res60_44"="CD3+ CD56dim+",
    "CD3pos_res60_60"="CD3+ CD56dim+",
    "CD3pos_res70_49"="CD3+ CD56dim+",
    "CD3pos_res70_57"="CD3+ CD56dim+",
    "CD3pos_res60_21"="CD3+ CD56dim+",
    "CD3pos_res60_27"="CD3+ CD56dim+",
    "CD3pos_res60_41"="CD3+ CD56dim+",
    "CD3pos_res60_14"="CD3+ CD56dim+",
    "CD3pos_res60_19"="CD3+ CD56dim+",
    "CD3pos_res60_20"="CD3+ CD56dim+",
    "CD3pos_res60_22"="CD3+ CD56dim+",
    "CD3pos_res60_35"="CD3+ CD56dim+",
    "CD3pos_res60_51"="CD3+ CD56dim+",
    "CD3pos_res70_1"="CD4 TCM",
    "CD3pos_res70_8"="CD4 TCM",
    "CD3pos_res70_9"="CD4 TCM",
    "CD3pos_res60_23"="CD4 Naive",
    "CD3pos_res60_13"="CD4 TEM",
    "CD3pos_res60_12"="CD4 TEM",
    "CD3pos_res60_17"="CD4 TEM",
    "CD3pos_res60_18"="CD4 TEM",
    "CD3pos_res70_14"="CD4 TEM",
    "CD3pos_res60_3"="CD4 TEM",
    "CD3pos_res60_8"="CD4 TEM",
    "CD3pos_res60_2"="CD4 TEM",
    "CD3pos_res60_48"="CD8 TEMRA",
    "CD3pos_res60_38"="CD8 TEMRA",
    "CD3pos_res60_28"="CD8 TEMRA",
    "CD3pos_res60_52"="CD8 TEMRA",
    "CD3pos_res60_59"="CD8 TEMRA",
    "CD3pos_res60_45"="CD8 Naive",
    "CD3pos_res60_25"="CD8 TEM-adv0",
    "CD3pos_res70_51"="CD8 TEM-adv0",
    "CD3pos_res60_31"="CD8 TEM-adv0",
    "CD3pos_res60_46"="CD8 TEM-adv0",
    "CD3pos_res60_24"="CD8 TEM-adv1",
    "CD3pos_res70_58"="CD8 TEM-adv1",
    "CD3pos_res70_43"="CD8 TEM-adv1",
    "CD3pos_res70_36"="CD8 TEM-adv1",
    "CD3pos_res60_15"="CD8 TEM-adv1",
    "CD3pos_res60_26"="CD8 TEM-adv1",
    "CD3pos_res60_33"="CD8 TEM-adv1",
    "CD3pos_res60_34"="CD8 TEM-adv1",
    "CD3pos_res60_42"="CD8 TEM-adv1",
    "CD3pos_res60_58"="CD8 TEM-adv1",
    "CD3pos_res60_47"="CD8 TEM-adv1",
    "CD3pos_res60_57"="CD8 TEM-adv1",
    "CD3pos_res60_40"="CD8 TEM-adv1",
    "CD3pos_res60_56"="CD8 TEM-adv1",
    "CD3pos_res60_43"="CD8 TEM-adv2",
    "CD3pos_res60_49"="CD8 TEM-adv2",
    "CD3pos_res60_55"="CD8 TEM-adv2",
    "CD3pos_res60_53"="Double Negative T",
    "CD3pos_res70_56"="Double Negative T",
    "CD3pos_res60_30"="Double Negative T",
    "CD3pos_res60_9"="gdT Vd1+",
    "CD3pos_res60_4"="gdT Vd1+",
    "CD3pos_res60_5"="gdT Vd2+ Vg9+",
    "CD3pos_res60_11"="gdT Vd2+ Vg9+",
    "CD3pos_res60_6"="gdT Vd2+ Vg9+",
    "CD3pos_res60_7"="gdT Vd2+ Vg9+",
    "CD3pos_res60_10"="gdT Vd2+ Vg9+",
    "ALL_res45_41"="DCs / HLA-DR+ APCs",
    "ALL_res55_44"="DCs / HLA-DR+ APCs",
    "ALL_res55_39"="DCs / HLA-DR+ APCs",
    "ALL_res45_7"="DCs / HLA-DR+ APCs",
    "ALL_res45_17"="Monocytes",
    "ALL_res55_12"="Monocytes",
    "ALL_res55_14"="Monocytes",
    "ALL_res55_16"="Monocytes",
    "ALL_res55_20"="Monocytes",
    "ALL_res55_21"="Monocytes",
    "ALL_res55_26"="Monocytes",
    "ALL_res55_8"="Monocytes",
    "ALL_res45_27"="Monocytes",
    "ALL_res55_15"="Monocytes",
    "ALL_res55_23"="Monocytes",
    "ALL_res55_45"="Monocytes",
    "ALL_res55_6"="Monocytes",
    "ALL_res45_42"="Monocytes",
    "ALL_res55_13"="Monocytes",
    "ALL_res55_4"="Monocytes"
  )
} else if (visualize_details == 'CD8_CD4_more_details') {
  parent1_population_dict_d35 <- c(
    "ALL_res45_14"="Exclude",
    "ALL_res45_5"="Exclude",
    "ALL_res45_45"="B",
    "ALL_res45_37"="CD56bright NK CD45RA+",
    "ALL_res45_38"="CD56bright NK CD45RA+ CD8mid GPR56+",
    "ALL_res45_44"="CD56bright NK CD45RA+ CD8mid GPR56+",
    "ALL_res45_40"="CD56bright NK CD45RA+ GPR56+",
    "ALL_res45_43"="CD56bright NK CD45RA+ GPR56+",
    "ALL_res55_47"="CD56bright NK CD45RA+ GPR56+",
    "ALL_res55_52"="CD56bright NK CD45RA+ GPR56+",
    "ALL_res45_34"="CD56dim NK CD45RA+ CD57+ CD8+ GPR56+",
    "CD3pos_res60_16"="CD56dim NK CD45RA+ CD57+ CD8+ GPR56+",
    "CD3pos_res60_29"="CD56dim NK CD45RA+ CD57+ GPR56+",
    "CD3pos_res60_37"="CD56dim NK CD45RA+ CD57+ GPR56+",
    "CD3pos_res60_54"="CD56dim NK CD45RA+ CD57+ GPR56+",
    "CD3pos_res70_40"="CD56dim NK CD45RA+ CD57+ GPR56+",
    "ALL_res45_35"="CD56dim NK CD45RA+ GPR56+",
    "ALL_res55_42"="CD56dim NK CD45RA+ GPR56+",
    "ALL_res55_49"="CD56dim NK CD45RA+ GPR56+",
    "CD3pos_res60_44"="CD3+ CD56dim CD8- CD45RA+ CD57+ GPR56+",
    "CD3pos_res60_60"="CD3+ CD56dim CD8- CD45RA+ CD57+ GPR56+",
    "CD3pos_res70_49"="CD3+ CD56dim CD8- CD45RA+ CD57+ GPR56+",
    "CD3pos_res70_57"="CD3+ CD56dim CD8- CD45RA+ CD57+ GPR56+",
    "CD3pos_res60_21"="CD3+ CD56dim CD8+ CD45RA+ CD57+ GPR56+",
    "CD3pos_res60_27"="CD3+ CD56dim CD8+ CD45RA+ CD57+ GPR56+",
    "CD3pos_res60_41"="CD3+ CD56dim CD8+ CD45RA+ CD57+ GPR56+",
    "CD3pos_res60_14"="CD3+ CD56dim CD8+ CD45RAmid CD57+ GPR56+",
    "CD3pos_res60_19"="CD3+ CD56dim CD8+ CD45RAmid CD57+ GPR56+",
    "CD3pos_res60_20"="CD3+ CD56dim CD8+ CD45RAmid CD57+ GPR56+",
    "CD3pos_res60_22"="CD3+ CD56dim CD8dim CD45RA+ CD57+ GPR56+",
    "CD3pos_res60_35"="CD3+ CD56dim CD8dim CD45RA+ CD57+ GPR56+",
    "CD3pos_res60_51"="CD3+ CD56dim CD8dim CD45RA+ CD57+ GPR56+",
    "CD3pos_res70_1"="CD4 TCM CD27+ CD29+",
    "CD3pos_res70_8"="CD4 TCM CD27+ CD29+",
    "CD3pos_res70_9"="CD4 TCM CD27+ CD29+ HLA-DR+",
    "CD3pos_res60_23"="CD4 Naive",
    "CD3pos_res70_14"="CD4 TEM CD27hi CD29+",
    "CD3pos_res60_12"="CD4 TEM CD27hi CD29+ CD57+",
    "CD3pos_res60_17"="CD4 TEM CD27hi CD29+ CD57+",
    "CD3pos_res60_18"="CD4 TEM CD27hi CD29+ CD57+",
    "CD3pos_res60_13"="CD4 TEM CD27dim CD29+ CD57+ GPR56+",
    "CD3pos_res60_2"="CD4 TEM CD27lo CD29+",
    "CD3pos_res60_3"="CD4 TEM CD27lo CD29+ GPR56+",
    "CD3pos_res60_8"="CD4 TEM CD27lo CD29+ GPR56+",
    "CD3pos_res60_48"="CD8 TEMRA CD27dim CD57+ GPR56+",
    "CD3pos_res60_38"="CD8 TEMRA CD27+ GPR56+",
    "CD3pos_res60_28"="CD8 TEMRA CD57+ GPR56+",
    "CD3pos_res60_52"="CD8 TEMRA CD57+ GPR56+",
    "CD3pos_res60_59"="CD8 TEMRA CD57+ GPR56+",
    "CD3pos_res60_45"="CD8 Naive",
    "CD3pos_res60_31"="CD8 TEM-adv0 CD27+",
    "CD3pos_res70_51"="CD8 TEM-adv0 CD27+",
    "CD3pos_res60_25"="CD8 TEM-adv0 CD27+ CD57+",
    "CD3pos_res60_46"="CD8dim TEM-adv0 GPR56+",
    "CD3pos_res60_24"="CD8 TEM-adv1",
    "CD3pos_res70_58"="CD8 TEM-adv1 CD27+",
    "CD3pos_res70_43"="CD8 TEM-adv1 CD27+ CD56dim",
    "CD3pos_res70_36"="CD8 TEM-adv1 CD27+ CD57+",
    "CD3pos_res60_15"="CD8 TEM-adv1 CD27+ CD57+ GPR56+",
    "CD3pos_res60_26"="CD8 TEM-adv1 CD27+ CD57+ GPR56+",
    "CD3pos_res60_33"="CD8 TEM-adv1 CD27+ CD57+ GPR56+",
    "CD3pos_res60_34"="CD8 TEM-adv1 CD27+ CD57+ GPR56+",
    "CD3pos_res60_42"="CD8 TEM-adv1 CD27+ CD57+ GPR56+",
    "CD3pos_res60_40"="CD8 TEM-adv1 CD27+ GPR56+",
    "CD3pos_res60_47"="CD8 TEM-adv1 CD27+ GPR56+",
    "CD3pos_res60_58"="CD8 TEM-adv1 CD57+ GPR56+",
    "CD3pos_res60_56"="CD8 TEM-adv1 GPR56+",
    "CD3pos_res60_57"="CD8 TEM-adv1 GPR56+",
    "CD3pos_res60_43"="CD8 TEM-adv2 CD57+ GPR56+",
    "CD3pos_res60_49"="CD8 TEM-adv2 CD57+ GPR56+",
    "CD3pos_res60_55"="CD8 TEM-adv2 CD57+ GPR56+",
    "CD3pos_res60_53"="Double Negative T",
    "CD3pos_res70_56"="Double Negative T",
    "CD3pos_res60_30"="Double Negative T",
    "CD3pos_res60_9"="gdT Vd1+ TEMRA",
    "CD3pos_res60_4"="gdT Vd1+ Naive",
    "CD3pos_res60_5"="gdT Vd2+ Vg9+ TEM",
    "CD3pos_res60_10"="gdT Vd2+ Vg9+ TEM",
    "CD3pos_res60_11"="gdT Vd2+ Vg9+ TEM",
    "CD3pos_res60_6"="gdT Vd2+ Vg9+ TEM",
    "CD3pos_res60_7"="gdT Vd2+ Vg9+ TEM",
    "ALL_res45_41"="DCs / HLA-DR+ APCs",
    "ALL_res55_44"="DCs / HLA-DR+ APCs",
    "ALL_res45_7"="DCs / HLA-DR+ APCs",
    "ALL_res55_39"="DCs / HLA-DR+ APCs",
    "ALL_res45_17"="Monocytes",
    "ALL_res55_12"="Monocytes",
    "ALL_res55_14"="Monocytes",
    "ALL_res55_16"="Monocytes",
    "ALL_res55_20"="Monocytes",
    "ALL_res55_21"="Monocytes",
    "ALL_res55_26"="Monocytes",
    "ALL_res55_8"="Monocytes",
    "ALL_res45_27"="Monocytes",
    "ALL_res55_15"="Monocytes",
    "ALL_res55_23"="Monocytes",
    "ALL_res55_45"="Monocytes",
    "ALL_res55_6"="Monocytes",
    "ALL_res45_42"="Monocytes",
    "ALL_res55_13"="Monocytes",
    "ALL_res55_4"="Monocytes"
  )
} else if (visualize_details == 'CD8_CD4_full_names') {
  parent1_population_dict_d35 <- c(
    "ALL_res45_14"="Exclude",
    "ALL_res45_5"="Exclude",
    "ALL_res45_45"="B",
    "ALL_res45_37"="CD56bright NK CD8- GPR56-",
    "ALL_res45_38"="CD56bright NK CD8mid GPR56+",
    "ALL_res45_44"="CD56bright NK CD8mid GPR56+",
    "ALL_res45_40"="CD56bright NK CD8- GPR56+",
    "ALL_res45_43"="CD56bright NK CD8- GPR56+",
    "ALL_res55_47"="CD56bright NK CD8- GPR56+",
    "ALL_res55_52"="CD56bright NK CD8- GPR56+",
    
    "ALL_res45_34"="CD56dim NK CD57+ CD8+ GPR56+",
    "CD3pos_res60_16"="CD56dim NK CD57+ CD8+ GPR56+",
    "CD3pos_res60_29"="CD56dim NK CD57+ CD8- GPR56+",
    "CD3pos_res60_37"="CD56dim NK CD57+ CD8- GPR56+",
    "CD3pos_res60_54"="CD56dim NK CD57+ CD8- GPR56+",
    "CD3pos_res70_40"="CD56dim NK CD57+ CD8- GPR56+",
    "ALL_res45_35"="CD56dim NK CD57- CD8- GPR56+",
    "ALL_res55_42"="CD56dim NK CD57- CD8- GPR56+",
    "ALL_res55_49"="CD56dim NK CD57- CD8- GPR56+",
    
    "CD3pos_res60_44"="CD3+ CD56dim CD8- CD45RA+ CD57+ GPR56+",
    "CD3pos_res60_60"="CD3+ CD56dim CD8- CD45RA+ CD57+ GPR56+",
    "CD3pos_res70_49"="CD3+ CD56dim CD8- CD45RA+ CD57+ GPR56+",
    "CD3pos_res70_57"="CD3+ CD56dim CD8- CD45RA+ CD57+ GPR56+",
    "CD3pos_res60_21"="CD3+ CD56dim CD8+ CD45RA+ CD57+ GPR56+",
    "CD3pos_res60_27"="CD3+ CD56dim CD8+ CD45RA+ CD57+ GPR56+",
    "CD3pos_res60_41"="CD3+ CD56dim CD8+ CD45RA+ CD57+ GPR56+",
    "CD3pos_res60_14"="CD3+ CD56dim CD8+ CD45RAmid CD57+ GPR56+",
    "CD3pos_res60_19"="CD3+ CD56dim CD8+ CD45RAmid CD57+ GPR56+",
    "CD3pos_res60_20"="CD3+ CD56dim CD8+ CD45RAmid CD57+ GPR56+",
    "CD3pos_res60_22"="CD3+ CD56dim CD8dim CD45RA+ CD57+ GPR56+",
    "CD3pos_res60_35"="CD3+ CD56dim CD8dim CD45RA+ CD57+ GPR56+",
    "CD3pos_res60_51"="CD3+ CD56dim CD8dim CD45RA+ CD57+ GPR56+",
    
    "CD3pos_res70_1"="CD4 TCM CD27+ CD29+ HLA-DR-",
    "CD3pos_res70_8"="CD4 TCM CD27+ CD29+ HLA-DR-",
    "CD3pos_res70_9"="CD4 TCM CD27+ CD29+ HLA-DR+",
    
    "CD3pos_res60_23"="CD4 Naive",
    
    "CD3pos_res70_14"="CD4 TEM CD27dim CD29+ CD57- GPR56-",
    "CD3pos_res60_12"="CD4 TEM CD27dim CD29+ CD57+ GPR56-",
    "CD3pos_res60_17"="CD4 TEM CD27dim CD29+ CD57+ GPR56-",
    "CD3pos_res60_18"="CD4 TEM CD27dim CD29+ CD57+ GPR56-",
    "CD3pos_res60_13"="CD4 TEM CD27dim CD29+ CD57+ GPR56+",
    "CD3pos_res60_2"="CD4 TEM CD27- CD29+ CD57- GPR56-",
    "CD3pos_res60_3"="CD4 TEM CD27- CD29+ CD57- GPR56+",
    "CD3pos_res60_8"="CD4 TEM CD27- CD29+ CD57- GPR56+",
    
    "CD3pos_res60_48"="CD8 TEMRA CD27dim CD57+ GPR56+",
    "CD3pos_res60_38"="CD8 TEMRA CD27dim CD57- GPR56+",
    "CD3pos_res60_28"="CD8 TEMRA CD27- CD57+ GPR56+",
    "CD3pos_res60_52"="CD8 TEMRA CD27- CD57+ GPR56+",
    "CD3pos_res60_59"="CD8 TEMRA CD27- CD57+ GPR56+",
    
    "CD3pos_res60_45"="CD8 Naive",
    
    "CD3pos_res60_31"="CD8 TEM-adv0 CD27+ CD57- GPR56-",
    "CD3pos_res70_51"="CD8 TEM-adv0 CD27+ CD57- GPR56-",
    "CD3pos_res60_25"="CD8 TEM-adv0 CD27+ CD57+ GPR56-",
    "CD3pos_res60_46"="CD8dim TEM-adv0 CD27dim CD57- GPR56+",
    
    "CD3pos_res60_24"="CD8 TEM-adv1 CD27- CD56- CD57- GPR56-",
    "CD3pos_res70_58"="CD8 TEM-adv1 CD27+ CD56- CD57- GPR56-",
    "CD3pos_res70_43"="CD8 TEM-adv1 CD27+ CD56dim CD57- GPR56-",
    "CD3pos_res70_36"="CD8 TEM-adv1 CD27dim CD56- CD57+ GPR56-",
    "CD3pos_res60_15"="CD8 TEM-adv1 CD27dim CD56- CD57+ GPR56+",
    "CD3pos_res60_26"="CD8 TEM-adv1 CD27dim CD56- CD57+ GPR56+",
    "CD3pos_res60_33"="CD8 TEM-adv1 CD27dim CD56- CD57+ GPR56+",
    "CD3pos_res60_34"="CD8 TEM-adv1 CD27dim CD56- CD57+ GPR56+",
    "CD3pos_res60_42"="CD8 TEM-adv1 CD27dim CD56- CD57+ GPR56+",
    "CD3pos_res60_40"="CD8 TEM-adv1 CD27dim CD56- CD57- GPR56+",
    "CD3pos_res60_47"="CD8 TEM-adv1 CD27dim CD56- CD57- GPR56+",
    "CD3pos_res60_58"="CD8 TEM-adv1 CD27- CD56- CD57+ GPR56+",
    "CD3pos_res60_56"="CD8 TEM-adv1 CD27- CD56- CD57- GPR56+",
    "CD3pos_res60_57"="CD8 TEM-adv1 CD27- CD56- CD57- GPR56+",
    
    "CD3pos_res60_43"="CD8 TEM-adv2 CD57+ GPR56+",
    "CD3pos_res60_49"="CD8 TEM-adv2 CD57+ GPR56+",
    "CD3pos_res60_55"="CD8 TEM-adv2 CD57+ GPR56+",
    
    "CD3pos_res60_53"="DN TEMRA CD57+ GPR56+",
    "CD3pos_res70_56"="DN TEMRA CD57+ GPR56+",
    "CD3pos_res60_30"="DN TEM CD45RAmid CD57- GPR56-",

    "CD3pos_res60_9"="gdT Vd1+ TEMRA CD27- CD57+ GPR56+",
    "CD3pos_res60_4"="gdT Vd1+ Naive CD27+ CD57- GPR56-",
    
    "CD3pos_res60_5"="gdT Vd2+ Vg9+ TEM CD45RAmid CD56dim GPR56+",
    "CD3pos_res60_10"="gdT Vd2+ Vg9+ TEM CD45RAmid/+ CD56dim CD57+ GPR56+",
    "CD3pos_res60_11"="gdT Vd2+ Vg9+ TEM CD45RAmid/+ CD56dim CD57+ GPR56+",
    "CD3pos_res60_6"="gdT Vd2+ Vg9+ TEM CD45RAmid/+ CD56dim CD57+ GPR56+",
    "CD3pos_res60_7"="gdT Vd2+ Vg9+ TEM CD45RAmid/+ CD56- CD57+ GPR56+",
    
    "ALL_res45_41"="DCs / HLA-DR+ APCs",
    "ALL_res55_44"="DCs / HLA-DR+ APCs",
    "ALL_res45_7"="DCs / HLA-DR+ APCs",
    "ALL_res55_39"="DCs / HLA-DR+ APCs",
    
    "ALL_res45_17"="Monocytes",
    "ALL_res55_12"="Monocytes",
    "ALL_res55_14"="Monocytes",
    "ALL_res55_16"="Monocytes",
    "ALL_res55_20"="Monocytes",
    "ALL_res55_21"="Monocytes",
    "ALL_res55_26"="Monocytes",
    "ALL_res55_8"="Monocytes",
    "ALL_res45_27"="Monocytes",
    "ALL_res55_15"="Monocytes",
    "ALL_res55_23"="Monocytes",
    "ALL_res55_45"="Monocytes",
    "ALL_res55_6"="Monocytes",
    "ALL_res45_42"="Monocytes",
    "ALL_res55_13"="Monocytes",
    "ALL_res55_4"="Monocytes"
  )
}
# ------

# ----- Building a heatmap for d35 -----
technical_channels <- c("FSC-A","FSC-H","FSC-W","SSC-A","SSC-H","SSC-W")
if (length(setdiff(exclusion_markers_list, technical_channels)) > 0) {
  excl_add_on <- paste("EXCL", paste(setdiff(exclusion_markers_list, technical_channels), collapse = "_"), sep = "_")
} else {
  excl_add_on <- 'ALL_MARKERS'
}

if (cohort_to_vis == 'Cohort 1, day 35') {
  coh_add_on <- 'OCd35'
} else if (cohort_to_vis == 'Cohort 2, day 35') {
  coh_add_on <- 'VCd35'
} else {
  coh_add_on <- 'BothCoh'
}

# defining markers
# -----
markers_all_d35_dict = list(
  "FSC-A" = "technical",
  "FSC-H" = "technical",
  "FSC-W" = "technical",
  "SSC-A" = "technical",
  "SSC-H" = "technical",
  "SSC-W" = "technical",
  "CD45" = "pan\nPBMC",
  "CD14" = "Myeloid/\nAPC",
  "CD303" = "Myeloid/\nAPC",
  "HLA-DR" = "Myeloid/\nAPC",
  "CD19" = "B",
  "CD56" = "NK",
  "CD3" = "T cell\ntype marker",
  "CD4" = "T cell\ntype marker",
  "CD8" = "T cell\ntype marker",
  "gdTCR" = "T cell\ntype marker",
  "TCRVd1" = "T cell\ntype marker",
  "TCRVd2" = "T cell\ntype marker",
  "TCRVg9" = "T cell\ntype marker",
  "CD45RA" = "T cell\nmaturation\nmarker",
  "CCR7" = "T cell\nmaturation\nmarker",
  "CD27" = "T cell\nmaturation\nmarker",
  "GPR56" = "State\nmarker",
  "CD57" = "State\nmarker",
  "CD29" = "State\nmarker",
  "CXCR4" = "State\nmarker"
)
# -----
markers_all_d35 <- names(markers_all_d35_dict)

# ----- defining broad cats -----
if (visualize_details == 'GeneralAnnotation') {
  populations_order_dict <- list(
    "Monocytes" = "Mono",
    "DCs / HLA-DR+ APCs" = "APC",
    "B" = "APC",
    "CD56dim NK" = "NK",
    "CD56bright NK" = "NK",
    "CD3+ CD56dim+" = "CD3+ CD56+",
    "CD4 Naive" = "CD4 T",
    "CD4 TCM" = "CD4 T",
    "CD4 TEM" = "CD4 T",
    "CD8 Naive" = "CD8 T",
    "CD8 TEM-adv0" = "CD8 T",
    "CD8 TEM-adv1" = "CD8 T",
    "CD8 TEM-adv2" = "CD8 T",
    "CD8 TEMRA" = "CD8 T",
    "gdT Vd1+" = "gdT",
    "gdT Vd2+ Vg9+" = "gdT",
    "Double Negative T" = "Other T"
  )
} else if (visualize_details == 'CD8_CD4_more_details') {
  populations_order_dict <- list(
    "Monocytes"="Monocytes",
    "DCs / HLA-DR+ APCs"="APC",
    "B"="APC",
    "CD56dim NK CD45RA+ CD57+ CD8+ GPR56+"="NK",
    "CD56dim NK CD45RA+ CD57+ GPR56+"="NK",
    "CD56dim NK CD45RA+ GPR56+"="NK",
    "CD56bright NK CD45RA+"="NK",
    "CD56bright NK CD45RA+ CD8mid GPR56+"="NK",
    "CD56bright NK CD45RA+ GPR56+"="NK",
    "CD3+ CD56dim CD8- CD45RA+ CD57+ GPR56+"="CD3+ CD56dim",
    "CD3+ CD56dim CD8+ CD45RA+ CD57+ GPR56+"="CD3+ CD56dim",
    "CD3+ CD56dim CD8+ CD45RAmid CD57+ GPR56+"="CD3+ CD56dim",
    "CD3+ CD56dim CD8dim CD45RA+ CD57+ GPR56+"="CD3+ CD56dim",
    "CD4 Naive"="CD4",
    "CD4 TCM CD27+ CD29+"="CD4",
    "CD4 TCM CD27+ CD29+ HLA-DR+"="CD4",
    "CD4 TEM CD27hi CD29+"="CD4",
    "CD4 TEM CD27hi CD29+ CD57+"="CD4",
    "CD4 TEM CD27dim CD29+ CD57+ GPR56+"="CD4",
    "CD4 TEM CD27lo CD29+"="CD4",
    "CD4 TEM CD27lo CD29+ GPR56+"="CD4",
    "CD8 Naive"="CD8",
    "CD8dim TEM-adv0 GPR56+"="CD8",
    "CD8 TEM-adv0 CD27+"="CD8",
    "CD8 TEM-adv0 CD27+ CD57+"="CD8",
    "CD8 TEM-adv1"="CD8",
    "CD8 TEM-adv1 CD27+"="CD8",
    "CD8 TEM-adv1 CD27+ CD56dim"="CD8",
    "CD8 TEM-adv1 CD27+ CD57+"="CD8",
    "CD8 TEM-adv1 CD27+ CD57+ GPR56+"="CD8",
    "CD8 TEM-adv1 CD27+ GPR56+"="CD8",
    "CD8 TEM-adv1 CD57+ GPR56+"="CD8",
    "CD8 TEM-adv1 GPR56+"="CD8",
    "CD8 TEM-adv2 CD57+ GPR56+"="CD8",
    "CD8 TEMRA CD27dim CD57+ GPR56+"="CD8",
    "CD8 TEMRA CD27+ GPR56+"="CD8",
    "CD8 TEMRA CD57+ GPR56+"="CD8",
    "gdT Vd1+ Naive"="gdT",
    "gdT Vd1+ TEMRA"="gdT",
    "gdT Vd2+ Vg9+ TEM"="gdT",
    "Double Negative T"="Other T"
  )
} else if (visualize_details == 'CD8_CD4_full_names') {
  populations_order_dict <- list(
    "Monocytes"="Monocytes",
    "DCs / HLA-DR+ APCs"="APC",
    "B"="APC",
    
    "CD56dim NK CD57+ CD8+ GPR56+"="NK",
    "CD56dim NK CD57+ CD8- GPR56+"="NK",
    "CD56dim NK CD57- CD8- GPR56+"="NK",
    
    "CD56bright NK CD8- GPR56-"="NK",
    "CD56bright NK CD8mid GPR56+"="NK",
    "CD56bright NK CD8- GPR56+"="NK",
    
    "CD3+ CD56dim CD8- CD45RA+ CD57+ GPR56+"="CD3+ CD56dim",
    "CD3+ CD56dim CD8+ CD45RA+ CD57+ GPR56+"="CD3+ CD56dim",
    "CD3+ CD56dim CD8+ CD45RAmid CD57+ GPR56+"="CD3+ CD56dim",
    "CD3+ CD56dim CD8dim CD45RA+ CD57+ GPR56+"="CD3+ CD56dim",
    
    "CD4 Naive"="CD4",
    
    "CD4 TCM CD27+ CD29+ HLA-DR-"="CD4",
    "CD4 TCM CD27+ CD29+ HLA-DR+"="CD4",
    
    "CD4 TEM CD27dim CD29+ CD57- GPR56-"="CD4",
    "CD4 TEM CD27dim CD29+ CD57+ GPR56-"="CD4",
    "CD4 TEM CD27dim CD29+ CD57+ GPR56+"="CD4",
    "CD4 TEM CD27- CD29+ CD57- GPR56-"="CD4",
    "CD4 TEM CD27- CD29+ CD57- GPR56+"="CD4",
    
    "CD8 Naive"="CD8",
    
    "CD8dim TEM-adv0 CD27dim CD57- GPR56+"="CD8",
    "CD8 TEM-adv0 CD27+ CD57- GPR56-"="CD8",
    "CD8 TEM-adv0 CD27+ CD57+ GPR56-"="CD8",
    
    "CD8 TEM-adv1 CD27- CD56- CD57- GPR56-"="CD8",
    "CD8 TEM-adv1 CD27+ CD56- CD57- GPR56-"="CD8",
    "CD8 TEM-adv1 CD27+ CD56dim CD57- GPR56-"="CD8",
    "CD8 TEM-adv1 CD27dim CD56- CD57+ GPR56-"="CD8",
    "CD8 TEM-adv1 CD27dim CD56- CD57+ GPR56+"="CD8",
    "CD8 TEM-adv1 CD27dim CD56- CD57- GPR56+"="CD8",
    "CD8 TEM-adv1 CD27- CD56- CD57+ GPR56+"="CD8",
    "CD8 TEM-adv1 CD27- CD56- CD57- GPR56+"="CD8",
    
    "CD8 TEM-adv2 CD57+ GPR56+"="CD8",
    
    "CD8 TEMRA CD27dim CD57+ GPR56+"="CD8",
    "CD8 TEMRA CD27dim CD57- GPR56+"="CD8",
    "CD8 TEMRA CD27- CD57+ GPR56+"="CD8",
    
    "DN TEMRA CD57+ GPR56+"="Other T",
    "DN TEMRA CD57+ GPR56+"="Other T",
    "DN TEM CD45RAmid CD57- GPR56-"="Other T",
    "gdT Vd1+ TEMRA CD27- CD57+ GPR56+"="gdT",
    "gdT Vd1+ Naive CD27+ CD57- GPR56-"="gdT",
    "gdT Vd2+ Vg9+ TEM CD45RAmid CD56dim GPR56+"="gdT",
    "gdT Vd2+ Vg9+ TEM CD45RAmid/+ CD56dim CD57+ GPR56+"="gdT",
    "gdT Vd2+ Vg9+ TEM CD45RAmid/+ CD56dim CD57+ GPR56+"="gdT",
    "gdT Vd2+ Vg9+ TEM CD45RAmid/+ CD56dim CD57+ GPR56+"="gdT",
    "gdT Vd2+ Vg9+ TEM CD45RAmid/+ CD56- CD57+ GPR56+"="gdT"
  )
}
# -----


populations_order <- names(populations_order_dict)
broad_cats <- c(paste(unlist(populations_order_dict)))
broad_cats <- factor(broad_cats, levels = unique(broad_cats))


# retrieving the expression matrix
cluster_map <- parent1_population_dict_d35
set_of_markers <- markers_all_d35[! markers_all_d35 %in% exclusion_markers_list]
marker_states <- paste(unlist(markers_all_d35_dict[set_of_markers]))
marker_states <- factor(marker_states, levels = unique(marker_states))
df_35 <- as.data.frame(main_dataset_d35$expr$orig[, set_of_markers, drop = FALSE])
df_35$cluster <- c(main_dataset_d35$clustering$Final_FlowSOM_Annotation_v2026$FinalCluster) 
df_35$cohort <- c(main_dataset_d35$anno$cell_anno$batch_coh)
df_35$cell_type <- recode(df_35$cluster,!!!cluster_map)
df_35$cohort <- ifelse(df_35$cohort == 'orig', 'Cohort 1, day 35', 'Cohort 2, day 35')
df_35$cell_type <- ifelse(df_35$cell_type == 'CD4 Naïve', 'CD4 Naive', ifelse(df_35$cell_type == 'CD8 Naïve', 'CD8 Naive', as.character(df_35$cell_type)))


# subsetting to a specific cohort:
if (cohort_to_vis == 'both') {
  df_35 <- df_35
} else if (cohort_to_vis == 'Cohort 1, day 35') {
  df_35 <- subset(df_35, cohort == 'Cohort 1, day 35')
} else if (cohort_to_vis == 'Cohort 2, day 35') {
  df_35 <- subset(df_35, cohort == 'Cohort 2, day 35')
}

# calculating medians per cell type:
median_matrix_d35 <- df_35 %>%
  group_by(cell_type) %>%
  summarise(across(-c(cluster, cohort), median, na.rm = TRUE))
median_matrix_d35 <- subset(median_matrix_d35, cell_type != 'Exclude')

# converting to matrix & scaling the median expression:
mat_d35 <- median_matrix_d35 %>%
  column_to_rownames("cell_type") %>%
  as.matrix()
if (scale) {
  # mat_d35 <- scale(mat_d35)
  
  mat_d35 <- apply(mat_d35, 2, function(x) {
    x <- as.numeric(scale(x))
    x <- pmax(pmin(x, 3), -3)
    x
  })
  
  rownames(mat_d35) <- median_matrix_d35$cell_type
  add_on <- 'Scaled'
} else {
  add_on <- 'Unscaled'
}

# rows as markers? or columns as markers?
mat_d35 <- if (heatmap_design == 'markers_rows_&_pops_cols') mat_d35[set_of_markers, populations_order] else mat_d35[populations_order, set_of_markers]

# choosing the color scheme
if (chosen_color_scale == "green_col") {
  color_scale <- c("#f2f2f2", "#d9ead9", "#b7d7b7", "#5aaE61", "#0b3d2e") 
} else if (chosen_color_scale == "purple_col") {
  color_scale <- c("#f2f2f2", "#ead6e8", "#d9a6cf", "#c05ab3", "#b30086")
} else {
  color_scale <- c("#3d85c6", "#f2f2f2", "#cd4025")
}
min_val <- round(range(mat_d35, na.rm = TRUE)[1])
max_val <- round(range(mat_d35, na.rm = TRUE)[2])
at_positions <- seq(min_val,max_val,(max_val - min_val) / (length(color_scale) - 1))
col_fun = colorRamp2(at_positions, color_scale)

# building a heatmap:
# specifying the output folder or the heatmap:
output_folder <- ''
pdf(glue('{output_folder}/d35_heatmap_broad_populations_{coh_add_on}_{add_on}_{excl_add_on}_{chosen_color_scale}_{vis_add_on}_SubSuper{sub_superscripts}.pdf'),
    width = 1.2 * ncol(mat_d35),
    height = 0.7 * nrow(mat_d35),
    family = "Helvetica")

# ----- relabelling the populations nicely -----
if (visualize_details == 'GeneralAnnotation') {
  subset_labels <- c(
    "CD56dim NK" = "CD56^{dim}~NK",
    "CD56bright NK" = "CD56^{bright}~NK",
    "DCs / HLA-DR+ APCs" = "DCs~'/'~'HLA-DR'^'+'~APCs",
    "CD3+ CD56dim+" = "CD3^'+'~CD56^{dim}",
    "CD4 Naive"   = "CD4~T['NV']",
    "CD4 TEM"     = "CD4~T['EM']",
    "CD4 TCM"     = "CD4~T['CM']",
    "CD8 Naive"   = "CD8~T['NV']",
    "CD8 TEM-adv0" = "CD8~T['EM-adv0']",
    "CD8 TEM-adv1" = "CD8~T['EM-adv1']",
    "CD8 TEM-adv2" = "CD8~T['EM-adv2']",
    "CD8 TEMRA"   = "CD8~T['EMRA']",
    "gdT Vd1+" = "gamma*delta~T~V*delta*'1'^'+'",
    "gdT Vd2+ Vg9+" = "gamma*delta~T~V*delta*'2'^'+'~V*gamma*'9'^'+'"
  )
  
} else if (visualize_details == 'CD8_CD4_more_details') {
  subset_labels <- c(
    "DCs / HLA-DR+ APCs" = "DCs~'/'~'HLA-DR'^'+'~APCs",
    
    "CD56dim NK CD45RA+ CD57+ CD8+ GPR56+" = "CD56^{dim}~NK~CD45RA^'+'~CD57^'+'~CD8^'+'~GPR56^'+'",
    "CD56dim NK CD45RA+ CD57+ GPR56+" = "CD56^{dim}~NK~CD45RA^'+'~CD57^'+'~GPR56^'+'",
    "CD56dim NK CD45RA+ GPR56+" = "CD56^{dim}~NK~CD45RA^'+'~GPR56^'+'",
    "CD56bright NK CD45RA+" = "CD56^{bright}~NK~CD45RA^'+'",
    "CD56bright NK CD45RA+ CD8mid GPR56+" ="CD56^{bright}~NK~CD45RA^'+'~CD8^{mid}~GPR56^'+'",
    "CD56bright NK CD45RA+ GPR56+" = "CD56^{bright}~NK~CD45RA^'+'~GPR56^'+'",
    
    "CD3+ CD56dim CD8- CD45RA+ CD57+ GPR56+" ="CD3^'+'~CD56^{dim}~CD8^'-'~CD45RA^'+'~CD57^'+'~GPR56^'+'",
    "CD3+ CD56dim CD8+ CD45RA+ CD57+ GPR56+" ="CD3^'+'~CD56^{dim}~CD8^'+'~CD45RA^'+'~CD57^'+'~GPR56^'+'",
    "CD3+ CD56dim CD8+ CD45RAmid CD57+ GPR56+" ="CD3^'+'~CD56^{dim}~CD8^'+'~CD45RA^{mid}~CD57^'+'~GPR56^'+'",
    "CD3+ CD56dim CD8dim CD45RA+ CD57+ GPR56+" ="CD3^'+'~CD56^{dim}~CD8^{dim}~CD45RA^'+'~CD57^'+'~GPR56^'+'",
    
    "CD4 Naive" ="CD4~T['NV']",
    "CD4 TCM CD27+ CD29+" = "CD4~T['CM']~CD27^'+'~CD29^'+'",
    "CD4 TCM CD27+ CD29+ HLA-DR+" = "CD4~T['CM']~CD27^'+'~CD29^'+'~HLA-DR^'+'",
    "CD4 TEM CD27hi CD29+" = "CD4~T['EM']~CD27^{hi}~CD29^'+'",
    "CD4 TEM CD27hi CD29+ CD57+" = "CD4~T['EM']~CD27^{hi}~CD29^'+'~CD57^'+'",
    "CD4 TEM CD27dim CD29+ CD57+ GPR56+" = "CD4~T['EM']~CD27^{dim}~CD29^'+'~CD57^'+'~GPR56^'+'",
    "CD4 TEM CD27lo CD29+" = "CD4~T['EM']~CD27^{lo}~CD29^'+'",
    "CD4 TEM CD27lo CD29+ GPR56+" = "CD4~T['EM']~CD27^{lo}~CD29^'+'~GPR56^'+'",
    
    "CD8 Naive" ="CD8~T['NV']",
    "CD8dim TEM-adv0 GPR56+" = "CD8^{dim}~T['EM-adv0']~GPR56^'+'",
    "CD8 TEM-adv0 CD27+" = "CD8~T['EM-adv0']~CD27^'+'",
    "CD8 TEM-adv0 CD27+ CD57+" = "CD8~T['EM-adv0']~CD27^'+'~CD57^'+'",
    "CD8 TEM-adv1" = "CD8~T['EM-adv1']",
    "CD8 TEM-adv1 CD27+" = "CD8~T['EM-adv1']~CD27^'+'",
    "CD8 TEM-adv1 CD27+ CD56dim" = "CD8~T['EM-adv1']~CD27^'+'~CD56^{dim}",
    "CD8 TEM-adv1 CD27+ CD57+" = "CD8~T['EM-adv1']~CD27^'+'~CD57^'+'",
    "CD8 TEM-adv1 CD27+ CD57+ GPR56+" = "CD8~T['EM-adv1']~CD27^'+'~CD57^'+'~GPR56^'+'",
    "CD8 TEM-adv1 CD27+ GPR56+" = "CD8~T['EM-adv1']~CD27^'+'~GPR56^'+'",
    "CD8 TEM-adv1 CD57+ GPR56+" = "CD8~T['EM-adv1']~CD57^'+'~GPR56^'+'",
    "CD8 TEM-adv1 GPR56+" = "CD8~T['EM-adv1']~GPR56^'+'",
    "CD8 TEM-adv2 CD57+ GPR56+" = "CD8~T['EM-adv2']~CD57^'+'~GPR56^'+'",
    "CD8 TEMRA CD27dim CD57+ GPR56+" = "CD8~T['EMRA']~CD27^{dim}~CD57^'+'~GPR56^'+'",
    "CD8 TEMRA CD27+ GPR56+" = "CD8~T['EMRA']~CD27^'+'~GPR56^'+'",
    "CD8 TEMRA CD57+ GPR56+" = "CD8~T['EMRA']~CD57^'+'~GPR56^'+'",
    
    
    "gdT Vd1+ Naive" = "gamma*delta~T~V*delta*'1'^'+'~T['NV']",
    "gdT Vd1+ TEMRA" = "gamma*delta~T~V*delta*'1'^'+'~T['EMRA']",
    "gdT Vd2+ Vg9+ TEM" = "gamma*delta~T~V*delta*'2'^'+'~V*gamma*'9'^'+'~T['EM']"
  )
} else if (visualize_details == 'CD8_CD4_full_names') {
  subset_labels <- c(
    "DCs / HLA-DR+ APCs" = "DCs~'/'~'HLA-DR'^'+'~APCs",
    
    "CD56dim NK CD57+ CD8+ GPR56+" = "CD56^{dim}~NK~CD57^'+'~CD8^'+'~GPR56^'+'",
    "CD56dim NK CD57+ CD8- GPR56+" = "CD56^{dim}~NK~CD57^'+'~CD8^'-'~GPR56^'+'",
    "CD56dim NK CD57- CD8- GPR56+" = "CD56^{dim}~NK~CD57^'-'~CD8^'-'~GPR56^'+'",
    
    "CD56bright NK CD8- GPR56-" = "CD56^{bright}~NK~CD8^'-'~GPR56^'-'",
    "CD56bright NK CD8mid GPR56+" ="CD56^{bright}~NK~CD8^{mid}~GPR56^'+'",
    "CD56bright NK CD8- GPR56+" = "CD56^{bright}~NK~CD8^'-'~GPR56^'+'",
    
    "CD3+ CD56dim CD8- CD45RA+ CD57+ GPR56+" ="CD3^'+'~CD56^{dim}~CD8^'-'~CD45RA^'+'~CD57^'+'~GPR56^'+'",
    "CD3+ CD56dim CD8+ CD45RA+ CD57+ GPR56+" ="CD3^'+'~CD56^{dim}~CD8^'+'~CD45RA^'+'~CD57^'+'~GPR56^'+'",
    "CD3+ CD56dim CD8+ CD45RAmid CD57+ GPR56+" ="CD3^'+'~CD56^{dim}~CD8^'+'~CD45RA^{mid}~CD57^'+'~GPR56^'+'",
    "CD3+ CD56dim CD8dim CD45RA+ CD57+ GPR56+" ="CD3^'+'~CD56^{dim}~CD8^{dim}~CD45RA^'+'~CD57^'+'~GPR56^'+'",
    
    "CD4 Naive" ="CD4~T['NV']",
    
    "CD4 TCM CD27+ CD29+ HLA-DR-" = "CD4~T['CM']~CD27^'+'~CD29^'+'~HLA-DR^'-'",
    "CD4 TCM CD27+ CD29+ HLA-DR+" = "CD4~T['CM']~CD27^'+'~CD29^'+'~HLA-DR^'+'",
    
    "CD4 TEM CD27dim CD29+ CD57- GPR56-" = "CD4~T['EM']~CD27^{dim}~CD29^'+'~CD57^'-'~GPR56^'-'",
    "CD4 TEM CD27dim CD29+ CD57+ GPR56-" = "CD4~T['EM']~CD27^{dim}~CD29^'+'~CD57^'+'~GPR56^'-'",
    "CD4 TEM CD27dim CD29+ CD57+ GPR56+" = "CD4~T['EM']~CD27^{dim}~CD29^'+'~CD57^'+'~GPR56^'+'",
    "CD4 TEM CD27- CD29+ CD57- GPR56-" = "CD4~T['EM']~CD27^'-'~CD29^'+'~CD57^'-'~GPR56^'-'",
    "CD4 TEM CD27- CD29+ CD57- GPR56+" = "CD4~T['EM']~CD27^'-'~CD29^'+'~CD57^'-'~GPR56^'+'",
    
    "CD8 Naive" ="CD8~T['NV']",
    
    "CD8dim TEM-adv0 CD27dim CD57- GPR56+" = "CD8^{dim}~T['EM-adv0']~CD27^{dim}~CD57^'-'~GPR56^'+'",
    "CD8 TEM-adv0 CD27+ CD57- GPR56-" = "CD8~T['EM-adv0']~CD27^'+'~CD57^'-'~GPR56^'-'",
    "CD8 TEM-adv0 CD27+ CD57+ GPR56-" = "CD8~T['EM-adv0']~CD27^'+'~CD57^'+'~GPR56^'-'",
    
    "CD8 TEM-adv1 CD27- CD56- CD57- GPR56-" = "CD8~T['EM-adv1']~CD27^'-'~CD56^'-'~CD57^'-'~GPR56^'-'",
    "CD8 TEM-adv1 CD27+ CD56- CD57- GPR56-" = "CD8~T['EM-adv1']~CD27^'+'~CD56^'-'~CD57^'-'~GPR56^'-'",
    "CD8 TEM-adv1 CD27+ CD56dim CD57- GPR56-" = "CD8~T['EM-adv1']~CD27^'+'~CD56^{dim}~CD57^'-'~GPR56^'-'",
    "CD8 TEM-adv1 CD27dim CD56- CD57+ GPR56-" = "CD8~T['EM-adv1']~CD27^{dim}~CD56^'-'~CD57^'+'~GPR56^'-'",
    "CD8 TEM-adv1 CD27dim CD56- CD57+ GPR56+" = "CD8~T['EM-adv1']~CD27^{dim}~CD56^'-'~CD57^'+'~GPR56^'+'",
    "CD8 TEM-adv1 CD27dim CD56- CD57- GPR56+" = "CD8~T['EM-adv1']~CD27^{dim}~CD56^'-'~CD57^'-'~GPR56^'+'",
    "CD8 TEM-adv1 CD27- CD56- CD57+ GPR56+" = "CD8~T['EM-adv1']~CD27^'-'~CD56^'-'~CD57^'+'~GPR56^'+'",
    "CD8 TEM-adv1 CD27- CD56- CD57- GPR56+" = "CD8~T['EM-adv1']~CD27^'-'~CD56^'-'~CD57^'-'~GPR56^'+'",
    
    "CD8 TEM-adv2 CD57+ GPR56+" = "CD8~T['EM-adv2']~CD57^'+'~GPR56^'+'",
    
    "CD8 TEMRA CD27dim CD57+ GPR56+" = "CD8~T['EMRA']~CD27^{dim}~CD57^'+'~GPR56^'+'",
    "CD8 TEMRA CD27dim CD57- GPR56+" = "CD8~T['EMRA']~CD27^{dim}~CD57^'-'~GPR56^'+'",
    "CD8 TEMRA CD27- CD57+ GPR56+" = "CD8~T['EMRA']~CD27^'-'~CD57^'+'~GPR56^'+'",
    
    
    "gdT Vd1+ Naive" = "gamma*delta~T~V*delta*'1'^'+'~T['NV']",
    "gdT Vd1+ TEMRA" = "gamma*delta~T~V*delta*'1'^'+'~T['EMRA']",
    "gdT Vd2+ Vg9+ TEM" = "gamma*delta~T~V*delta*'2'^'+'~V*gamma*'9'^'+'~T['EM']"
  )
}
# -----


if (heatmap_design == 'markers_rows_&_pops_cols') {
  cell_types <- colnames(mat_d35)
  cell_types_expr <- ifelse(cell_types %in% names(subset_labels),subset_labels[cell_types],paste0("'", cell_types, "'"))
} else if (heatmap_design == 'pops_rows_&_markers_cols') {
  cell_types <- rownames(mat_d35)
  cell_types_expr <- ifelse(cell_types %in% names(subset_labels),subset_labels[cell_types],paste0("'", cell_types, "'"))
}

if (sub_superscripts == 'AlreadyInR') {
  used_cell_types <- parse(text = cell_types_expr)
} else if (sub_superscripts == 'PlainAIVersion') {
  used_cell_types <- cell_types
}

heatmap_plt <- Heatmap(mat_d35,
                       name = glue("{add_on} median expression"),
                       cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       
                       column_order = if (heatmap_design == 'markers_rows_&_pops_cols') populations_order else set_of_markers,
                       column_labels = if (heatmap_design == 'markers_rows_&_pops_cols') used_cell_types else set_of_markers,
                       
                       row_order    = if (heatmap_design == 'pops_rows_&_markers_cols') populations_order else set_of_markers,
                       row_labels = if (heatmap_design == 'pops_rows_&_markers_cols') used_cell_types else set_of_markers,
                       
                       heatmap_legend_param = list(
                         at = seq(min_val,max_val, 1),
                         title_gp = gpar(fontsize = 25, fontfamily = "Helvetica"),
                         labels_gp = gpar(fontsize = 25, fontfamily = "Helvetica"),
                         direction = 'horizontal',
                         # title_position = "leftcenter-rot",
                         title_position = "topcenter",
                         title_gap = unit(6, "cm"),
                         # legend_height = unit(12, "cm"),
                         grid_width = unit(4, "cm"),
                         legend_width = unit(12, "cm"),
                         title = "Scaled median normalized MFI"
                       ),
                       row_names_side = "left",
                       col = col_fun,
                       rect_gp = gpar(col = "white", lwd = 1),
                       border_gp = gpar(lwd = 2),
                       
                       row_names_gp = gpar(fontsize = 25, fontfamily = "Helvetica"),
                       column_names_gp = gpar(fontsize = 25, fontfamily = "Helvetica"),
                       row_title_gp = gpar(fontsize = 25, fontfamily = "Helvetica"),
                       column_title_gp = gpar(fontsize = 25, fontfamily = "Helvetica"),
                       
                       cluster_column_slices = FALSE,
                       # column_split = col_split_list,
                       column_split = if (heatmap_design == 'markers_rows_&_pops_cols') broad_cats else marker_states,
                       column_title = if ((heatmap_design == 'markers_rows_&_pops_cols') & (population_split == 'NoSplitInAnnotation')) NULL else "%s",
                       # column_gap = if (heatmap_design == 'markers_rows_&_pops_cols') unit(1.5, "mm") else NULL,
                       column_gap = unit(2, "mm"),
                       
                       cluster_row_slices = FALSE,
                       # row_split = row_split_list,
                       row_split = if (heatmap_design == 'pops_rows_&_markers_cols') broad_cats else marker_states,
                       row_title = if ((heatmap_design == 'pops_rows_&_markers_cols') & (population_split == 'NoSplitInAnnotation')) NULL else "%s",
                       column_names_rot = 45,
                       # row_gap = if (heatmap_design == 'pops_rows_&_markers_cols') unit(1.5, "mm") else NULL,
                       row_gap = unit(2, "mm"),
                       
                       border = TRUE
                       )
draw(heatmap_plt, heatmap_legend_side = "bottom", padding = unit(c(8, 200, 8, 8), "mm"))
dev.off()
# -----

