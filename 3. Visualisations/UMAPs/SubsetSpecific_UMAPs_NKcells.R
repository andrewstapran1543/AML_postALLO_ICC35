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
main_dataset_d35_nk <- filter_fcd(fcd = main_dataset_d35, cell_ids = rownames(main_dataset_d35$expr$orig)[main_dataset_d35$clustering$Final_FlowSOM_Annotation_v2026$Parent3Name %in% c("CD56bright NK", "CD56dim NK")])

# ----- running the UMAP on the CD8 subset: -----
main_dataset_d35_nk <- runUMAP(fcd = main_dataset_d35_nk, 
                                input_type = "pca", 
                                data_slot = "HARMONY_PCA_norm",
                                prefix = 'NK_SubClust',
                                seed = 91)


subcluster_or_total <- 'subcluster'
if (subcluster_or_total == 'subcluster') {
  working_dataset <- main_dataset_d35_nk
  umap_name <- "NK_SubClust_pca_HARMONY_PCA_norm"
  vis_name_add_on <- "NK_subsets_"
  
} else if (subcluster_or_total == 'total_dataset') {
  working_dataset <- main_dataset_d35
  umap_name <- "250610_FinalAnnot_pca_HARMONY_PCA_norm"
  vis_name_add_on <- "FullUMAP_NK_subsets_"
}
working_dataset$anno$cell_anno$batch_coh_FULL_NAME <- ifelse(working_dataset$anno$cell_anno$batch_coh == 'orig', 'Cohort 1', 'Cohort 2')


# ----- Visualizing the populations on the UMAP: -----
# specifying the output folder to NK-specific UMAPs:
output_folder <- ''
level_selected <- 'FullPopulationName_NO_EXTRA_ANNOT'

NK_dim <- c("CD56dim NK CD45RA+ GPR56+", "CD56dim NK CD45RA+ CD57+ GPR56+", "CD56dim NK CD45RA+ CD57+ CD8+ GPR56+")
NK_bright <- c("CD56bright NK CD45RA+", "CD56bright NK CD45RA+ GPR56+", "CD56bright NK CD45RA+ CD8mid GPR56+")
NK_dim_GPR56pos_CD57pos <- c("CD56dim NK CD45RA+ CD57+ GPR56+", "CD56dim NK CD45RA+ CD57+ CD8+ GPR56+")
NK_dim_GPR56pos_CD57pos_CD45RApos_HighRes <- c("CD56dim NK CD45RA+ CD57+ GPR56+")
NK_dim_GPR56pos_CD57pos_CD45RApos_CD8pos_HighRes <- c("CD56dim NK CD45RA+ CD57+ CD8+ GPR56+")
x <- working_dataset$clustering$Final_FlowSOM_Annotation_v2026[[level_selected]]

# vis_name <- glue('{vis_name_add_on}Fig3_NK_dim')
# working_dataset$clustering$Final_FlowSOM_Annotation_v2026$ClusterToVis <-
#   factor(ifelse(x %in% NK_dim,'CD56dim NK', 'Other NK'),levels = c("Other NK",'CD56dim NK'))
# my_palette <- c('CD56dim NK' = "#93C47D", "Other NK" = "#ededed")
# 
# vis_name <- glue('{vis_name_add_on}Fig3_NK_bright')
# working_dataset$clustering$Final_FlowSOM_Annotation_v2026$ClusterToVis <-
#   factor(ifelse(x %in% NK_bright,'CD56bright NK', 'Other NK'),levels = c("Other NK",'CD56bright NK'))
# my_palette <- c('CD56bright NK' = "#93C47D", "Other NK" = "#ededed")
# 
# vis_name <- glue('{vis_name_add_on}Fig3_NK_dim_bright')
# working_dataset$clustering$Final_FlowSOM_Annotation_v2026$ClusterToVis <-
#   factor(ifelse(x %in% NK_bright,'CD56bright NK',ifelse(x %in% NK_dim, 'CD56dim NK', 'Other NK')),
#          levels = c("Other NK",'CD56bright NK', 'CD56dim NK'))
# my_palette <- c('CD56bright NK' = "#377EB8", 'CD56dim NK' = "#93C47D", "Other NK" = "#ededed")
# 
# vis_name <- glue('{vis_name_add_on}Fig3_NK_dim_GPR56pos_CD57pos')
# working_dataset$clustering$Final_FlowSOM_Annotation_v2026$ClusterToVis <-
#   factor(ifelse(x %in% NK_dim_GPR56pos_CD57pos,'CD56dim NK CD57+ GPR56+', 'Other NK'),levels = c("Other NK", 'CD56dim NK CD57+ GPR56+'))
# my_palette <- c('CD56dim NK CD57+ GPR56+' = "#93C47D", "Other NK" = "#ededed")
# 
# vis_name <- glue('{vis_name_add_on}Fig3_NK_dim_GPR56pos_CD57pos_CD45RApos_HighRes')
# working_dataset$clustering$Final_FlowSOM_Annotation_v2026$ClusterToVis <-
#   factor(ifelse(x %in% NK_dim_GPR56pos_CD57pos_CD45RApos_HighRes,'CD56dim NK CD45RA+ CD57+ GPR56+', 'Other NK'),levels = c("Other NK", 'CD56dim NK CD45RA+ CD57+ GPR56+'))
# my_palette <- c('CD56dim NK CD45RA+ CD57+ GPR56+' = "#93C47D", "Other NK" = "#ededed")
# 
# vis_name <- glue('{vis_name_add_on}Fig3_NK_dim_GPR56pos_CD57pos_CD45RApos_CD8pos_HighRes')
# working_dataset$clustering$Final_FlowSOM_Annotation_v2026$ClusterToVis <-
#   factor(ifelse(x %in% NK_dim_GPR56pos_CD57pos_CD45RApos_CD8pos_HighRes,"CD56dim NK CD45RA+ CD57+ CD8+ GPR56+", 'Other NK'),levels = c("Other NK", "CD56dim NK CD45RA+ CD57+ CD8+ GPR56+"))
# my_palette <- c("CD56dim NK CD45RA+ CD57+ CD8+ GPR56+" = "#93C47D", "Other NK" = "#ededed")
# 
# vis_name <- glue('{vis_name_add_on}Fig3_NK_dim_GPR56pos_CD57pos_subsets')
# working_dataset$clustering$Final_FlowSOM_Annotation_v2026$ClusterToVis <-
#   factor(ifelse(x %in% NK_dim_GPR56pos_CD57pos_CD45RApos_HighRes,'CD56dim NK CD45RA+ CD57+ GPR56+',
#                 ifelse(x %in% NK_dim_GPR56pos_CD57pos_CD45RApos_CD8pos_HighRes, "CD56dim NK CD45RA+ CD57+ CD8+ GPR56+", 'Other NK')),
#          levels = c("Other NK", "CD56dim NK CD45RA+ CD57+ CD8+ GPR56+", 'CD56dim NK CD45RA+ CD57+ GPR56+'))
# my_palette <- c('CD56dim NK CD45RA+ CD57+ GPR56+' = "#93C47D", "CD56dim NK CD45RA+ CD57+ CD8+ GPR56+" = "#377EB8", "Other NK" = "#ededed")

vis_name <- glue('{vis_name_add_on}Fig3_NK_dim_HighResAndOthers_NK_bright')
working_dataset$clustering$Final_FlowSOM_Annotation_v2026$ClusterToVis <-
  factor(ifelse(x %in% NK_dim_GPR56pos_CD57pos_CD45RApos_HighRes,'CD56dim NK CD45RA+ CD57+ GPR56+',
                ifelse(x %in% NK_dim, "Other CD56dim NK",
                       ifelse(x %in% NK_bright, 'CD56bright NK', "Other"))),
         levels = c("Other CD56dim NK", "CD56bright NK", 'CD56dim NK CD45RA+ CD57+ GPR56+', 'Other'))
my_palette <- c('CD56dim NK CD45RA+ CD57+ GPR56+' = "#6b51ad", "Other CD56dim NK" = "#B4A7D6", "CD56bright NK" = "#d29eb9", "Other" = "#ededed")



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
    ncol = 1,
    override.aes = list(size = 6, alpha = 1)
  ))
print(plot1)
dev.off()
# -----




