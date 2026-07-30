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
library(readr)

##################################################
#### Setting up the paths to cyCondor objects ####
##################################################
# Saving path to the main object with all cells
saving_path <- ''
# Saving path to the smaller subsetted object (used for UMAPs & visualisations)
saving_path_small <- ''


#######################################################
#### Step 1. Reading the CSV files into the object ####
#######################################################
folder_with_filtered_csv_files <- ''
annotation_table_path <- ''
# Reading the filtered CSV files into the cyCondor object
# max_cell limit should be set depending on what is the maximum number of cells per file
condor_dataset <- prep_fcd(data_path = folder_with_filtered_csv_files,
                        max_cell = 2000000,
                        useCSV = TRUE,
                        transformation = "auto_logi",
                        remove_param = c("Time"),
                        anno_table = annotation_table_path,
                        filename_col = "filename",
                        seed = 91,
                        verbose = TRUE)
# Saving the main cyCondor object
qsave(condor_dataset, saving_path)


#####################################################
#### Step 2.1 Running PCA on unintegrated object ####
#####################################################
# Excluding the following columns from PCA analysis: FSC-A/H/W, SSC-A/H/W
condor_dataset <- runPCA(fcd = condor_dataset,
                      data_slot = "orig",
                      discard = T,
                      markers = c("FSC-A","FSC-H","FSC-W","SSC-A","SSC-H","SSC-W"),
                      seed = 91)
# Saving the main cyCondor object
qsave(condor_dataset, saving_path)


####################################################################################################
#### Step 2.2 Adding TOTAL_BATCH column - combining the cohort & experimental batch information ####
####################################################################################################
condor_dataset$anno$cell_anno$batch_exp <- as.factor(condor_dataset$anno$cell_anno$batch_exp)
condor_dataset$anno$cell_anno$TOTAL_BATCH <- paste(condor_dataset$anno$cell_anno$batch_coh, condor_dataset$anno$cell_anno$batch_exp, sep = "_")
# Saving the main cyCondor object
qsave(condor_dataset, saving_path)


########################################################################################
#### Step 3. Using Harmony algorithm to integrate FCS data on TOTAL_BATCH parameter ####
########################################################################################
# Running the integration using the harmonize_PCA() function
# The FCS files were normalised using batch controls beforehand in FlowJo - so integrating per cohort is more important than integrating per batch
prefix_harmony_pca <- ''
condor_dataset <- harmonize_PCA(fcd = condor_dataset,
                                batch_var = c("TOTAL_BATCH"),
                                data_slot = "orig",
                                prefix = glue("{prefix_harmony_pca}"),
                                seed = 91)
# Saving the main cyCondor object
qsave(condor_dataset, saving_path)


###############################################################################
#### Step 4. Running FlowSOM clustering on integrated Principal Components ####
###############################################################################
# This is the main clustering step. Further the CD3+ T cells would be taken aside and clustered separately
# Running the FlowSOM clustering at different resolutions for CD3+ Lymphocytes
# What would be the prefix for clustering slot in the cyCondor object?
prefix_selected <- ''
for (nclusters in c(30,40,50,60,70)) {
  condor_dataset <- runFlowSOM(fcd = condor_dataset,
                               input_type = "pca",
                               data_slot = glue("{prefix_harmony_pca}_norm"),
                               nClusters = nclusters,
                               seed = 91,
                               ret_model = FALSE,
                               prefix = prefix_selected)
}
# Saving the main cyCondor object
qsave(condor_dataset, saving_path)


