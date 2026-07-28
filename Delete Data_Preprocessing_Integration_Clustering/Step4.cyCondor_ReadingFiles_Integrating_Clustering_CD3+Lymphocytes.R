####################################
#### Activating the environment ####
####################################
working_directory_path <- "/g/zaugg/stapran/PhD_Projects/FACS_FlowSum_Seurat_Analysis"
setwd(working_directory_path)
environment_path <- "/g/zaugg/stapran/PhD_Projects/FACS_FlowSum_Seurat_Analysis/renv/library/linux-rocky-8.8/R-4.4/x86_64-pc-linux-gnu"
.libPaths(environment_path)
renv::activate()


###############################
#### Loading the libraries ####
###############################
library(cyCONDOR)
library(qs)


##################################################
#### Setting up the paths to cyCondor objects ####
##################################################
# Main dataset (created in Step 2 - with main FlowSOM clustering)
main_full_dataset_path <- ''
# Main CD3+ dataset (subset to CD3+ Lymphocytes from the big dataset)
saving_path <- ''
# Small subset of the CD3+ dataset (for visualisation)
saving_path_small <- ''


##########################################
#### Step 1. Reading the main dataset ####
##########################################
condor_dataset <- qread(main_full_dataset_path)


#############################################
#### Step 2. Subsetting to CD3+ clusters ####
#############################################
# CD3+ clusters were selected based on ridgeplot-based visualisation of CD3 expression across clusters in the main dataset
cd3pos_clusters_res45 <- c('1','2','3','9','10','15','16','18','19','20','21','22','23','24','25','26','29','30','31','32')
# subsetting the main dataset to only CD3+ clusters
cd3pos_cells_full <- filter_fcd(fcd = condor_dataset,
                                cell_ids = rownames(condor_dataset$expr$orig)[condor_dataset$clustering$FlowSOM_250513_ALL_CELLS_pca_HARMONY_PCA_norm_k_45$FlowSOM
                                                                              %in% cd3pos_clusters_res45])
# removing the main dataset to free up RAM for the operations on CD3+ dataset
rm(condor_dataset)
# saving the main CD3+ dataset
qsave(cd3pos_cells_full, saving_path)


#########################################################
#### Step 3. Running PCA on unintegrated CD3+ object ####
#########################################################
# Since we have subsetted the main dataset to CD3+ clusters - we have to run the integration again. That is why we first run unintegrated PCA
# Excluding the following columns from PCA analysis: FSC-A/H/W, SSC-A/H/W
cd3pos_prefix_unintegrated_PCA <- ''
cd3pos_cells_full <- runPCA(fcd = cd3pos_cells_full,
                            data_slot = "orig",
                            discard = T,
                            markers = c("FSC-A","FSC-H","FSC-W","SSC-A","SSC-H","SSC-W"),
                            seed = 91,
                            prefix = cd3pos_prefix_unintegrated_PCA)
# saving the main CD3+ dataset
qsave(cd3pos_cells_full, saving_path)


##############################################################################################
#### Step 3.1. Using Harmony algorithm to integrate CD3+ dataset on TOTAL_BATCH parameter ####
##############################################################################################
# Running the integration using the harmonize_PCA() function
cd3pos_prefix_HARMintegrated_PCA <- ''
cd3pos_cells_full <- harmonize_PCA(fcd = cd3pos_cells_full,
                                   batch_var = c("TOTAL_BATCH"),
                                   data_slot = glue("{cd3pos_prefix_unintegrated_PCA}_orig"),
                                   prefix = cd3pos_prefix_HARMintegrated_PCA,
                                   seed = 91)
# saving the main CD3+ dataset
qsave(cd3pos_cells_full, saving_path)


###############################################################################
#### Step 4. Running FlowSOM clustering on integrated Principal Components ####
###############################################################################
# This is the second clustering step specifically for CD3+ Lymphocytes. Now we have information about both: initial FlowSOM clustering & CD3+ specific clustering
# Running the FlowSOM clustering at different resolutions for CD3+ Lymphocytes
# What would be the prefix for CD3+ clustering slot in the cyCondor object?
cd3pos_prefix_selected <- ''
for (nclusters in c(30,40,50,60,70)) {
  cd3pos_cells_full <- runFlowSOM(fcd = cd3pos_cells_full,
                                  input_type = "pca",
                                  data_slot = glue("{cd3pos_prefix_HARMintegrated_PCA}_norm"),
                                  nClusters = nclusters,
                                  seed = 91,
                                  ret_model = FALSE,
                                  prefix = cd3pos_prefix_selected) 
}
# saving the main CD3+ dataset
qsave(cd3pos_cells_full, saving_path)