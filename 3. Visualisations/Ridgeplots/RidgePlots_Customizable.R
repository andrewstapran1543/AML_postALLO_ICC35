################################
#### Importing the packages ####
################################

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
library(dplyr)
library(tidyverse)
library(Hmisc)
library(glue)
library(ggridges)

##############################
#### reading the dataset: ####
##############################
# path to the subsetted cyCondor dataset with full final annotations:
main_dataset_d35_path <- ''
main_dataset_d35 <- qread(main_dataset_d35_path)


#############################################
#### choosing the annotation resolution: ####
#############################################
visualize_details <- 'CD8_CD4_full_names'
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
    "CD3pos_res60_13"="CD4 TEM CD27hi CD29+ CD57+ GPR56+",
    "CD3pos_res60_2"="CD4 TEM CD27lo CD29+",
    "CD3pos_res60_3"="CD4 TEM CD27lo CD29+ GPR56+",
    "CD3pos_res60_8"="CD4 TEM CD27lo CD29+ GPR56+",
    "CD3pos_res60_48"="CD8 TEMRA CD27+ CD57+ GPR56+",
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
    "CD3pos_res60_53"="Double Negative T TEMRA CD57+ GPR56+",
    "CD3pos_res70_56"="Double Negative T TEMRA CD57+ GPR56+",
    "CD3pos_res60_30"="Double Negative T TEM CD57- GPR56-",
    "CD3pos_res60_9"="gdT Vd1+ TEMRA CD27- CD57+ GPR56+",
    "CD3pos_res60_4"="gdT Vd1+ Naïve CD27+ CD57- GPR56-",
    "CD3pos_res60_5"="gdT Vd2+ Vg9+ TEM CD56dim CD57- GPR56+",
    "CD3pos_res60_10"="gdT Vd2+ Vg9+ TEM CD56dim CD57+ GPR56+",
    "CD3pos_res60_11"="gdT Vd2+ Vg9+ TEM CD56dim CD57+ GPR56+",
    "CD3pos_res60_6"="gdT Vd2+ Vg9+ TEM CD56dim CD57+ GPR56+",
    "CD3pos_res60_7"="gdT Vd2+ Vg9+ TEM CD56- CD57+ GPR56+",
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


##############################################
#### specifying the markers to visualize: ####
##############################################
all_markers <- c("CD45", "CD14", "CD303","HLA-DR", "CD19",
                 "CD56", "CD3", "CD4", "CD8", "gdTCR", "TCRVd1", "TCRVd2", "TCRVg9",
                 "CD45RA", "CCR7", "CD27", "GPR56", "CD57", "CD29", "CXCR4")


##########################################
#### reshaping the expression matrix: ####
##########################################
cluster_map <- parent1_population_dict_d35
set_of_markers <- all_markers
df_35 <- as.data.frame(main_dataset_d35$expr$orig[, set_of_markers, drop = FALSE])
df_35$cluster <- c(main_dataset_d35$clustering$Final_FlowSOM_Annotation_v2026$FinalCluster)
df_35$cohort <- c(main_dataset_d35$anno$cell_anno$batch_coh)
df_35_long <- df_35 %>% pivot_longer(cols = -c(cluster, cohort), names_to = "marker", values_to = "expression")
df_35_long$subset <- recode(df_35_long$cluster,!!!cluster_map)
df_35_long <- subset(df_35_long, subset %nin% c('Exclude'))
df_35_long$marker <- factor(df_35_long$marker, levels = set_of_markers)
df_35_long$cohort <- ifelse(df_35_long$cohort == 'orig', 'Cohort 1', 'Cohort 2')
df_35_long$subset <- ifelse(df_35_long$subset == 'CD4 Naïve', 'CD4 Naive', ifelse(df_35_long$subset == 'CD8 Naïve', 'CD8 Naive', as.character(df_35_long$subset)))
df_total <- df_35_long



################################################################
#### which populations and cohorts to show? in which order? ####
################################################################
pops_to_show <- 'all'
# 'shared' / 'all' / 'selected'

# ----- ordering the populations: -----
population_levels <- c(
  "CD4 Naive",
  "CD4 TCM CD27+ CD29+ HLA-DR-",
  "CD4 TCM CD27+ CD29+ HLA-DR+",
  "CD4 TEM CD27dim CD29+ CD57- GPR56-",
  "CD4 TEM CD27dim CD29+ CD57+ GPR56-",
  "CD4 TEM CD27dim CD29+ CD57+ GPR56+",
  "CD4 TEM CD27- CD29+ CD57- GPR56-",
  "CD4 TEM CD27- CD29+ CD57- GPR56+",
  "CD8 Naive",
  "CD8dim TEM-adv0 CD27dim CD57- GPR56+",
  "CD8 TEM-adv0 CD27+ CD57- GPR56-",
  "CD8 TEM-adv0 CD27+ CD57+ GPR56-",
  "CD8 TEM-adv1 CD27- CD56- CD57- GPR56-",
  "CD8 TEM-adv1 CD27+ CD56- CD57- GPR56-",
  "CD8 TEM-adv1 CD27+ CD56dim CD57- GPR56-",
  "CD8 TEM-adv1 CD27dim CD56- CD57+ GPR56-",
  "CD8 TEM-adv1 CD27dim CD56- CD57+ GPR56+",
  "CD8 TEM-adv1 CD27dim CD56- CD57- GPR56+",
  "CD8 TEM-adv1 CD27- CD56- CD57+ GPR56+",
  "CD8 TEM-adv1 CD27- CD56- CD57- GPR56+",
  "CD8 TEM-adv2 CD57+ GPR56+",
  "CD8 TEMRA CD27dim CD57+ GPR56+",
  "CD8 TEMRA CD27dim CD57- GPR56+",
  "CD8 TEMRA CD27- CD57+ GPR56+",
  "gdT Vd1+ TEMRA CD27- CD57+ GPR56+",
  "gdT Vd1+ Naïve CD27+ CD57- GPR56-",
  "gdT Vd2+ Vg9+ TEM CD56dim CD57- GPR56+",
  "gdT Vd2+ Vg9+ TEM CD56dim CD57+ GPR56+",
  "gdT Vd2+ Vg9+ TEM CD56- CD57+ GPR56+",
  "Double Negative T TEMRA CD57+ GPR56+",
  "Double Negative T TEM CD57- GPR56-",
  "CD3+ CD56dim CD8- CD45RA+ CD57+ GPR56+",
  "CD3+ CD56dim CD8+ CD45RA+ CD57+ GPR56+",
  "CD3+ CD56dim CD8+ CD45RAmid CD57+ GPR56+",
  "CD3+ CD56dim CD8dim CD45RA+ CD57+ GPR56+",
  "CD56dim NK CD57+ CD8+ GPR56+",
  "CD56dim NK CD57+ CD8- GPR56+",
  "CD56dim NK CD57- CD8- GPR56+",
  "CD56bright NK CD8- GPR56-",
  "CD56bright NK CD8mid GPR56+",
  "CD56bright NK CD8- GPR56+",
  "B",
  "DCs / HLA-DR+ APCs",
  "Monocytes"
)
# -----

if (pops_to_show == 'all') {
  population_levels <- population_levels
} else if (pops_to_show == 'shared') {
  shared_list_pops = intersect(parent1_population_dict_d35,parent1_population_dict_d100)
  population_levels <- intersect(population_levels, shared_list_pops)
} else if (pops_to_show == 'selected') {
  # population_levels <- c(
  #   "CD4 Naive","CD4 TCM","CD4 TEM",
  #   "CD8 Naive","CD8dim TEM-adv0","CD8 TEM-adv0","CD8 TEM-adv1","CD8 TEM-adv2","CD8 TEMRA"
  # )
  population_levels <- c("CD8 Naive","CD8 TEM-adv0","CD8 TEM-adv1","CD8 TEM-adv2","CD8 TEMRA")
  # population_levels <- c(
  #   "CD8 Naive",
  # 
  #   "CD8 TEM-adv0 CD27+ CD57- GPR56-",
  #   "CD8 TEM-adv0 CD27+ CD57+ GPR56-",
  #   "CD8dim TEM-adv0 CD27- CD57- GPR56+",
  # 
  #   "CD8 TEM-adv1 CD27- CD56- CD57- GPR56-",
  #   "CD8 TEM-adv1 CD27+ CD56- CD57- GPR56-",
  #   "CD8 TEM-adv1 CD27+ CD56dim CD57- GPR56-",
  #   "CD8 TEM-adv1 CD27+ CD56- CD57+ GPR56-",
  #   "CD8 TEM-adv1 CD27+ CD56- CD57+ GPR56+",
  #   "CD8 TEM-adv1 CD27+ CD56- CD57- GPR56+",
  #   "CD8 TEM-adv1 CD27- CD56- CD57+ GPR56+",
  #   "CD8 TEM-adv1 CD27- CD56- CD57- GPR56+",
  # 
  #   "CD8 TEM-adv2 CD57+ GPR56+",
  # 
  #   "CD8 TEMRA CD27+ CD57+ GPR56+",
  #   "CD8 TEMRA CD27+ CD57- GPR56+",
  #   "CD8 TEMRA CD27- CD57+ GPR56+"
  # )
}

df_visualise = subset(df_total, subset %in% population_levels)
df_visualise$subset <- factor(df_visualise$subset, levels = rev(population_levels))


cohorts_levels <- c("Cohort 1", "Cohort 2")
df_visualise = subset(df_visualise, cohort %in% cohorts_levels)

# markers_levels <- c("CD14", "CD56", "CD3", "CD4", "CD8", "CD45RA", "CCR7")
# markers_levels <- c("CD45RA", "CCR7")
# markers_levels <- c("CD27", "CD57", "GPR56")
markers_levels <- c("CD45", "CD14", "CD303","HLA-DR", "CD19",
                    "CD56", "CD3", "CD4", "CD8", "gdTCR", "TCRVd1", "TCRVd2", "TCRVg9",
                    "CD45RA", "CCR7", "CD27", "GPR56", "CD57", "CD29", "CXCR4")
df_visualise$marker <- factor(df_visualise$marker, levels = markers_levels)
df_visualise$cohort <- factor(df_visualise$cohort,levels = c("Cohort 2", "Cohort 1"))   # adapt to your exact values




################################
#### Building the ridgeplot ####
################################
level_selected <- 'Parent1'
markers_path <- paste(markers_levels, collapse = "_")
ncol_selected = length(markers_levels)


# Extra piece of code for splitting the ridgeplot in two parts - for Suppl. Fig1
# df_visualise$lineage <- ifelse(grepl("^CD4", df_visualise$subset),"CD4 states","CD8 states")
# vline_df <- data.frame(
#   marker = c("CD45RA", "CD45RA","CD45RA", "CCR7"),
#   xintercept = c(2, 3, 3.75, 2.75)
# )
# label_df <- data.frame(
#   marker  = c("CCR7", "CCR7",
#               "CD45RA", "CD45RA", "CD45RA", "CD45RA"),
#   lineage = c("CD4 states", "CD4 states",
#               "CD4 states", "CD4 states", "CD4 states", "CD4 states"),
#   x       = c(1.7, 3.7,
#               1.0, 2.5, 3.4, 4.2),
#   label   = c("neg", "pos",
#               "neg", "mid", "mid/+", "pos"))
# df_visualise$marker <- factor(df_visualise$marker, levels = c("CD45RA", "CCR7"))
# vline_df$marker     <- factor(vline_df$marker,     levels = c("CD45RA", "CCR7"))
# label_df$marker     <- factor(label_df$marker,     levels = c("CD45RA", "CCR7"))



# specifying the output folder to ridgeplots:
output_folder <- ''
pdf(glue('{output_folder}/CustomPopulationsPaper_June2026_v2_fullalpha_ColorScheme2.pdf'),
    width =  1.5 * length(unique(df_visualise$marker)),
    height = 0.4 * length(unique(df_visualise$subset)))
plt <- ggplot(df_visualise, aes(x = expression, y = subset, fill = cohort)) +
  geom_density_ridges(
    aes(height = after_stat(ndensity)),
    stat = "density_ridges",
    scale = 0.7,
    alpha = 1,
    color = "black",
    linewidth = 0.6
    # rel_min_height = 0.0001
  ) + 
  # facet_grid(lineage ~ marker,scales = "free_y",space = "free_y", switch = "y") +
  facet_wrap(~ marker, ncol = length(markers_levels)) +
  labs(x = "",y = "",title = "") +
  scale_y_discrete(
    # labels = function(x) parse(text = subset_labels[x]),
    expand = expansion(add = c(0.7, 1.5))
  ) +
  coord_cartesian(clip = "off") + 
  scale_fill_manual(
    values = c(
      "Cohort 1" = "#a8c5df", # "#a8c5df" / "#1f77b4"
      "Cohort 2" = "#d17575" # "#d17575" / "#d62728"
    ),
    breaks = c("Cohort 1", "Cohort 2"),
    labels = c("Cohort 1      ", "Cohort 2")
  ) + 
  guides(fill = guide_legend(nrow = 1)) +
  theme_minimal() + 
  theme(
    strip.text.x = element_text(size = 15),
    strip.text.y = element_blank(),
    strip.background.y = element_blank(),
    
    panel.spacing.x = unit(1, "lines"),
    panel.spacing.y = unit(0.5, "lines"),
    # strip.text = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.key.width = unit(0.6, "cm"),
    legend.key.height = unit(0.45, "cm"),
    legend.spacing.x = unit(20, "cm"),
    
    legend.box.spacing = unit(-0.25, "cm"),
    axis.text = element_text(colour = "black"),
    text = element_text(color = "black"),
    axis.text.x = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 15, color = "black"),
    plot.margin = margin(1, 1, 1, 1),
) + scale_x_continuous(
  breaks = c(0, 2, 4)
)
  # geom_vline(
  #   data = vline_df,
  #   aes(xintercept = xintercept),
  #   linetype = "dashed",
  #   color = "darkred",
  #   linewidth = 0.7,
  #   inherit.aes = FALSE
  # ) + 
  # geom_label(
  #   data = label_df,
  #   aes(x = x, y = Inf, label = label),
  #   inherit.aes = FALSE,
  #   vjust = 1.2,
  #   size = 3,
  #   label.size = 0.25,
  #   # fill = "white",
  #   color = "black",
  #   label.padding = unit(0.08, "lines")
  # )
print(plt)
dev.off()




