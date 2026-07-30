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
library(lisi)
library(ggplot2)
library(ggrastr)
library(readr)


######################################################################################################
#### Step 0. Reading the cyCondor subset (to avoid crashing - recommended to read smaller subset) ####
######################################################################################################
path_to_cyCondor_dataset <- ''
condor_dataset_subset <- qread(path_to_cyCondor_dataset)
unintegrated_pca_slot_name <- ''
# usually something like {unintegrated_pca_name}_orig
integrated_pca_slot_name <- ''
# usually something like {integrated_pca_name}_norm

###########################################################################
#### Step 1. Checking the integration efficiency on the entire dataset ####
###########################################################################
# Integration efficiency calculated using LISI score on UMAP1 and UMAP2
# Selecting a parameter on which the integration efficiency will be estimated
# Can be cohort / experimental batch etc. (max LISI score = maximum number of levels in the selected variable)
chosen_column = "batch_coh"

# Extracting UMAP & medadata columns for unintegrated PCA based UMAP 
unintegrated_UMAP <- cbind(condor_dataset_subset$umap[[unintegrated_pca_slot_name]], condor_dataset_subset$anno$cell_anno)
# calculate LISI score for unintegrated PCA based UMAP
LISI_unintegrated <- compute_lisi(unintegrated_UMAP[,c(1,2)], unintegrated_UMAP, c(chosen_column))
colnames(LISI_unintegrated) <- "lisi"
# combine LISI score with cell annotation
lisi_matrix_unintegrated <- cbind(unintegrated_UMAP, LISI_unintegrated)
lisi_matrix_unintegrated$type <- "unintegrated"

# Extracting UMAP & medadata columns for UMAP calculated on Harmong-integrated PCA
harmony_integrated_UMAP <- cbind(condor_dataset_subset$umap[[integrated_pca_slot_name]], condor_dataset_subset$anno$cell_anno)
# calculate LISI score for batch variable 'exp'
LISI_harmony_integrated <- compute_lisi(harmony_integrated_UMAP[,c(1,2)], harmony_integrated_UMAP, c(chosen_column))
colnames(LISI_harmony_integrated) <- "lisi"
# combine LISI score with cell annotation
lisi_matrix_harmony_integrated <- cbind(harmony_integrated_UMAP, LISI_harmony_integrated)
lisi_matrix_harmony_integrated$type <- "harmony_integrated"

# combine unintegrated and Harmony-integrated LISI matrices
lisi_matrix_full <- rbind(lisi_matrix_harmony_integrated, lisi_matrix_unintegrated)
lisi_matrix_full$type <- factor(lisi_mat$type, levels = c("unintegrated", "harmony_integrated"))

# Visualizing the integration efficiency using LISI score
# The closer to maximum value - the more efficient is the integration
p <- ggplot(data = lisi_matrix_full, aes(y = lisi, x = type, fill = type)) +
  geom_jitter_rast(alpha = 0.01, scale =0.5) +
  geom_violin(alpha = 0.8) +
  scale_fill_manual(values= c("#1C75BC", "#BE1E2D"))+
  theme_bw() +
  theme(aspect.ratio = 2, panel.grid = element_blank(),
        text= element_text(size=16, color= "black")) +
  xlab("")+
  ylab("LISI score")

# Saving the LISI csv table and pdf plot
LISI_csv_table_save_path <- ''
write_csv(lisi_matrix_full, LISI_csv_table_save_path)
LISI_pdf_plot <- ''
pdf(LISI_pdf_plot, width = 12, height = 7)
print(p)
dev.off()

