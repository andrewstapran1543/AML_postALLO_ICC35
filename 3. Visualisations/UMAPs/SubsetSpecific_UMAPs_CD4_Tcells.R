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
library(viridis)
source('~/NewRidgePlottingFunction.R')
source('~/Custom_CyCONDOR_Functions.R')

# ----- reading the dataset: -----
# path to the subsetted cyCondor dataset with full final annotations:
main_dataset_d35_path <- ''
main_dataset_d35 <- qread('')

# ----- subsetting to CD8 cells only: -----
main_dataset_d35_cd4 <- filter_fcd(fcd = main_dataset_d35, cell_ids = rownames(main_dataset_d35$expr$orig)[main_dataset_d35$clustering$Final_FlowSOM_Annotation_v2026$Parent3Name == "CD4"])

# ----- running the UMAP on the CD8 subset: -----
main_dataset_d35_cd4 <- runUMAP(fcd = main_dataset_d35_cd4, 
                                input_type = "pca", 
                                data_slot = "HARMONY_PCA_norm",
                                prefix = 'CD4_SubClust',
                                seed = 91)

subcluster_or_total <- 'subcluster'
if (subcluster_or_total == 'subcluster') {
  working_dataset <- main_dataset_d35_cd4
  umap_name <- "CD4_SubClust_pca_HARMONY_PCA_norm"
  vis_name_add_on <- "CD4_subsets_"

} else if (subcluster_or_total == 'total_dataset') {
  working_dataset <- main_dataset_d35
  umap_name <- "250610_FinalAnnot_pca_HARMONY_PCA_norm"
  vis_name_add_on <- "FullUMAP_CD4_subsets_"
}
working_dataset$anno$cell_anno$batch_coh_FULL_NAME <- ifelse(working_dataset$anno$cell_anno$batch_coh == 'orig', 'Cohort 1', 'Cohort 2')


# ----- Visualizing the populations on the UMAP: -----
# specifying the output folder to CD4-specific UMAPs:
output_folder <- ''
level_selected <- 'FullPopulationName_NO_EXTRA_ANNOT'
CD4_TCM <- c("CD4 CCR7mid CD45RA- CD27+ CD29+","CD4 CCR7mid CD45RA- CD27+ CD29+ HLA-DR+")
CD4_TCM_HLADRpos <- c("CD4 CCR7mid CD45RA- CD27+ CD29+ HLA-DR+")
CD4_TCM_HLADRneg <- c("CD4 CCR7mid CD45RA- CD27+ CD29+")
CD4_NV <- c("CD4 Naive CD27+")
CD4_TEM <- c("CD4 TEM CD27hi CD29+","CD4 TEM CD27hi CD29+ CD57+","CD4 TEM CD27hi CD29+ CD57+ GPR56+","CD4 TEM CD27lo CD29+","CD4 TEM CD27lo CD29+ GPR56+")

x <- working_dataset$clustering$Final_FlowSOM_Annotation_v2026[[level_selected]]

# vis_name <- glue('{vis_name_add_on}Fig3_CD4_TCM_HighRes_GGPLOT')
# working_dataset$clustering$Final_FlowSOM_Annotation_v2026$ClusterToVis <-
#   factor(ifelse(x %in% CD4_TCM,'CD4 TCM', 'Other CD4 T'),levels = c("Other CD4 T",'CD4 TCM'))
# my_palette <- c('CD4 TCM' = "#C27BA0", "Other CD4 T" = "#ededed")

# vis_name <- glue('{vis_name_add_on}Fig3_CD4_TCM_SubPops_HighRes_GGPLOT')
# working_dataset$clustering$Final_FlowSOM_Annotation_v2026$ClusterToVis <-
#   factor(ifelse(x %in% CD4_TCM_HLADRpos,'CD4 TCM CD27+ CD29+ HLA-DR+',
#                 ifelse(x %in% CD4_TCM_HLADRneg, 'CD4 TCM CD27+ CD29+ HLA-DR-',
#                        'Other CD4 T')), levels = c("Other CD4 T",'CD4 TCM CD27+ CD29+ HLA-DR+', 'CD4 TCM CD27+ CD29+ HLA-DR-'))
# my_palette <- c('CD4 TCM CD27+ CD29+ HLA-DR+' = "#93C47D", 'CD4 TCM CD27+ CD29+ HLA-DR-' = "#B4A7D6", "Other CD4 T" = "#ededed")
# 
# vis_name <- glue('{vis_name_add_on}Fig3_CD4_NV_HighRes_GGPLOT')
# working_dataset$clustering$Final_FlowSOM_Annotation_v2026$ClusterToVis <-
#   factor(ifelse(x %in% CD4_NV,'CD4 NV', 'Other CD4 T'),levels = c("Other CD4 T",'CD4 NV'))
# my_palette <- c('CD4 NV' = "#93C47D", "Other CD4 T" = "#ededed")

vis_name <- glue('{vis_name_add_on}Fig3_CD4_NV_CD4_TCM_CD4_TEM_HighRes_GGPLOT')
working_dataset$clustering$Final_FlowSOM_Annotation_v2026$ClusterToVis <-
  factor(ifelse(x %in% CD4_NV,'CD4 NV',
                ifelse(x %in% CD4_TCM, 'CD4 TCM',ifelse(x %in% CD4_TEM, 'CD4 TEM', 'Other'))),levels = c('CD4 NV', 'CD4 TCM', "CD4 TEM", "Other"))
my_palette <- c('CD4 NV' = "#d29eb9", 'CD4 TCM' = "#b884c6","CD4 TEM" = "#6b51ad", "Other" = "#ededed")






# ----- ggplot 2 version of the UMAP -----
pdf(file = glue("{output_folder}{vis_name}.pdf"), width = 3.5, height = 8, family = "Helvetica")
umap_df <- data.frame(
  UMAP1 = working_dataset$umap[[umap_name]][,1],
  UMAP2 = working_dataset$umap[[umap_name]][,2],
  ClusterToVis = working_dataset$clustering$Final_FlowSOM_Annotation_v2026$ClusterToVis,
  Cohort = working_dataset$anno$cell_anno$batch_coh_FULL_NAME   # or your cohort column
)
umap_df$ClusterToVis <- as.character(umap_df$ClusterToVis)
umap_df$ClusterToVis <- factor(umap_df$ClusterToVis,levels = names(my_palette))
umap_df <- umap_df[order(umap_df$ClusterToVis), ]
set.seed(91)
umap_df <- umap_df[sample(seq_len(nrow(umap_df))), ]

# annotation positions per facet/cohort
arrow_df <- umap_df %>%
  dplyr::group_by(Cohort) %>%
  dplyr::summarise(
    x_min = min(UMAP1, na.rm = TRUE),
    x_max = max(UMAP1, na.rm = TRUE),
    y_min = min(UMAP2, na.rm = TRUE),
    y_max = max(UMAP2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    x0 = x_min,
    y0 = y_min,
    arrow_len_x = 0.14 * (x_max - x_min),
    arrow_len_y = 0.14 * (y_max - y_min),
    
    x_umap1 = x0 + arrow_len_x / 2,
    y_umap1 = y0 - 0.04 * (y_max - y_min),
    
    x_umap2 = x0 - 0.04 * (x_max - x_min),
    y_umap2 = y0 + arrow_len_y / 2
  )

arrow_df <- arrow_df %>%
  dplyr::mutate(
    x0 = x_min - 0.04 * (x_max - x_min),
    y0 = y_min - 0.04 * (y_max - y_min),
    x_umap1 = x0 + arrow_len_x / 2,
    y_umap1 = y0 - 0.06 * (y_max - y_min),
    x_umap2 = x0 - 0.06 * (x_max - x_min),
    y_umap2 = y0 + arrow_len_y / 2
  )

plot1 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = ClusterToVis)) +
  geom_point(size = 0.2, alpha = 1) +
  # ggrastr::rasterize(
  #   geom_point(size = 0.2, alpha = 0.7),
  #   dpi = 300
  # ) + 
  facet_wrap(~ Cohort, ncol = 1) +
  
  geom_segment(
    data = arrow_df,
    aes(x = x0, y = y0, xend = x0 + arrow_len_x, yend = y0),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.22, "cm")),
    linewidth = 0.8,
    color = "black"
  ) +
  geom_segment(
    data = arrow_df,
    aes(x = x0, y = y0, xend = x0, yend = y0 + arrow_len_y),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.22, "cm")),
    linewidth = 0.8,
    color = "black"
  ) +
  geom_text(
    data = arrow_df,
    aes(x = x_umap1, y = y_umap1, label = "UMAP 1"),
    inherit.aes = FALSE,
    family = "Helvetica",
    fontface = "bold",
    size = 3
  ) +
  geom_text(
    data = arrow_df,
    aes(x = x_umap2, y = y_umap2, label = "UMAP 2"),
    inherit.aes = FALSE,
    angle = 90,
    family = "Helvetica",
    fontface = "bold",
    size = 3
  ) +
  scale_color_manual(values = my_palette) +
  labs(color = NULL) +
  theme_void() +
  theme(
    text = element_text(family = "Helvetica"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(2, "cm")
  ) +
  guides(colour = guide_legend(
    ncol = 2,
    override.aes = list(size = 6, alpha = 1)
  ))
print(plot1)
dev.off()
# -----




