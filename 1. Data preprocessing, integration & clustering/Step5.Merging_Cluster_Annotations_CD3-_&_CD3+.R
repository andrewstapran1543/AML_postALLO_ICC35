####################################
#### Activating the environment ####
####################################
working_directory_path <- "~/FACS_FlowSum_Seurat_Analysis"
setwd(working_directory_path)
environment_path <- "~/FACS_FlowSum_Seurat_Analysis/renv/library/linux-rocky-8.8/R-4.4/x86_64-pc-linux-gnu"
.libPaths(environment_path)
renv::activate()


###############################
#### Loading the libraries ####
###############################
library(cyCONDOR)
library(qs)
library(glue)
library(dplyr)


# The purpose of this part is to fuse the full dataset clustering and CD3+ specific clustering
# Optimal resolution was for 1) full object FlowSOM clustering and 2) CD3+ specific FlowSOM clustering was selected using clustree package
# For each cluster that splits at the resolution following the optimal resolution - we have checked whether the resulting clusters are similar/different in phenotype. The resolution after optimal was only selected if the resulting clusters differ from each other in phenotype. E.g., resolution 60 (nclusters = 60) is optimal for CD3+ specific clustering; cluster 35 splits into clusters 23 and 24 at resolution 70. If clusters 23 and 24 are similar in phenotype - we simply keep them as cluster 35 at resolution 60. Otherwise - we label then as cluster 23 and cluster 24 at resolution 70


#####################################################
#### Reading the full and CD3+ specific datasets ####
#####################################################
# Path to the full cyCondor object (all cells with initial FlowSOM clustering)
all_cells_path <- ''
all_cells <- qread(all_cells_path)
# Path to CD3+ Lymphocytes subset of the full cyCondor object (with CD3+ specific clustering)
t_cells_path <- ''
t_cells <- qread(t_cells_path)

# Saving path to the full final cyCondor object
all_cells_save_path <- ''
# Saving path to the small subset of the full final cyCondor object (visualisation purposes)
all_cells_vis_save_path <- ''


########################################################################################
#### Step 1.1. Reading the clustering annotations for the full dataset (CD3- cells) ####
########################################################################################
# Optimal resolution for the full dataset: nclusters = 45
# The next level to check whether any of the res45 clusters split meaningfully: nclusters = 55
# Specifiying the full dataset clustering slot prefix:
full_dataset_clustering_slot <- ''
# Extracting the resolution 45 and 55 cluster_ID information:
all_cells_annotation_table <- data.frame(
  res_45 = all_cells$clustering[[glue("{full_dataset_clustering_slot}_k_45")]][["FlowSOM"]],
  res_55 = all_cells$clustering[[glue("{full_dataset_clustering_slot}_k_55")]][["FlowSOM"]],
  row.names = rownames(all_cells$clustering[[glue("{full_dataset_clustering_slot}_k_45")]])
)
# Which clusters at resolution 45 split into meaningfully (phenotypically) different clusters at resolution 55:
clusters_at_res45_splitting_into_res55 <- c(4,11,12,8,13,28,6,33,36,39)
all_cells_annotation_table$Final_Cluster_ID <- ifelse(all_cells_annotation_table$res_45 %in% clusters_at_res45_splitting_into_res55,
                                                      paste0('ALL_res55_',as.numeric(all_cells_annotation_table$res_55)),
                                                      paste0('ALL_res45_',as.numeric(all_cells_annotation_table$res_45)))
all_cells$clustering$Final_ALLCells_Clustering_res45and55 <- all_cells_annotation_table
all_cells_annotation_table$res_45 <- NULL
all_cells_annotation_table$res_55 <- NULL
all_cells_annotation_table <- all_cells_annotation_table %>% mutate(cell_ID = rownames(.))


###########################################################################
#### Step 1.2. Reading the clustering annotations for the CD3+ dataset ####
###########################################################################
# Optimal resolution for the full dataset: nclusters = 60
# The next level to check whether any of the res60 clusters split meaningfully: nclusters = 70
# Specifiying the full dataset clustering slot prefix:
cd3pos_dataset_clustering_slot <- ''
# Extracting the resolution 60 and 70 cluster_ID information:
cd3pos_cells_annotation_table <- data.frame(
  res_60 = t_cells$clustering[[glue("{cd3pos_dataset_clustering_slot}_k_60")]][["FlowSOM"]],
  res_70 = t_cells$clustering[[glue("{cd3pos_dataset_clustering_slot}_k_70")]][["FlowSOM"]],
  row.names = rownames(t_cells$clustering[[glue("{cd3pos_dataset_clustering_slot}_k_60")]])
)
# Which clusters at resolution 60 split into meaningfully (phenotypically) different clusters at resolution 70:
cd3pos_clusters_at_res60_splitting_into_res70 <- c(1,39,32,36,50)
cd3pos_cells_annotation_table$Final_Cluster_ID <- ifelse(cd3pos_cells_annotation_table$res_60 %in% cd3pos_clusters_at_res60_splitting_into_res70,
                                                         paste0('CD3pos_res70_',as.numeric(cd3pos_cells_annotation_table$res_70)),
                                                         paste0('CD3pos_res60_',as.numeric(cd3pos_cells_annotation_table$res_60)))
t_cells$clustering$Final_CD3posCells_Clustering_res60and70 <- cd3pos_cells_annotation_table
cd3pos_cells_annotation_table$res_60 <- NULL
cd3pos_cells_annotation_table$res_70 <- NULL
cd3pos_cells_annotation_table <- cd3pos_cells_annotation_table %>% mutate(cell_ID = rownames(.))


############################################################################
#### Step 2. Merging the cluster_ID information for CD3- and CD3+ cells ####
############################################################################
# Getting the CD3- cluster_ID information
all_cells_annotation_table <- data.frame(
  Final_Cluster_ID = all_cells$clustering$Final_ALLCells_Clustering_res45and55$Final_Cluster_ID,
  row.names = rownames(all_cells$clustering$FlowSOM_250513_ALL_CELLS_pca_HARMONY_PCA_norm_k_45)
)
all_cells_annotation_table <- all_cells_annotation_table %>% mutate(cell_ID = rownames(.))

# Merging the CD3- cluster_ID information with CD3+ cluster_IDs
merged_df <- all_cells_annotation_table %>%
  left_join(cd3pos_cells_annotation_table, by = "cell_ID") %>%
  mutate(FinalCluster = coalesce(Final_Cluster_ID.y, Final_Cluster_ID.x)) %>%
  select(-Final_Cluster_ID.x, -Final_Cluster_ID.y)
rownames(merged_df) <- merged_df$cell_ID
merged_df <- merged_df %>% select(-cell_ID)

# Recording the final merged cluster_ID information into the full dataset
all_cells$clustering$Final_FlowSOM_Annotation <- merged_df
all_cells$clustering$Final_FlowSOM_Annotation$FinalCluster <- as.factor(all_cells$clustering$Final_FlowSOM_Annotation$FinalCluster)
desired_order <- levels(all_cells$clustering$Final_FlowSOM_Annotation$FinalCluster)

# Saving the full dataset with final cluster_ID annotations
qsave(all_cells, all_cells_save_path)





