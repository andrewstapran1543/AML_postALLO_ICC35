####################################
#### Activating the environment ####
####################################
working_directory_path <- "~/FACS_FlowSum_Seurat_Analysis"
setwd(working_directory_path)
environment_path <- "~/FACS_FlowSum_Seurat_Analysis/renv/library/linux-rocky-8.8/R-4.4/x86_64-pc-linux-gnu"
.libPaths(environment_path)
renv::activate()


# Importantly! After every modification to the main object - to carry over the changes to the smaller visualisation subset we have to subset the big object again and run the UMAP on harmony integreated PCA slot:


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
#### Step 1. Reading the main object ####
#######################################################
main_cycondor_dataset <- qread(saving_path)


###############################################################
#### Step 2. Subsetting the main object to a small dataset ####
###############################################################
# Subsetting the main dataset to a smaller object 
# In our case we are subsetting by equal number of cells per cohort (specified in subset_parameter)
subset_parameter <- ''
subset_parameter_size <- 150000
condor_dataset_subset <- subset_fcd_byparam(main_cycondor_dataset,
                                            param = subset_parameter,
                                            size = subset_parameter_size,
                                            seed = 91)
# Saving the subsetted cyCondor object
qsave(condor_dataset_subset, saving_path_small)


######################################################################################
#### Step 3. Running UMAP on unintegrated and integrated PCAs - for visualisation ####
######################################################################################
# Running the UMAP for unintegrated PCA:
unintegrated_pca_slot <- ''
# usually something like {unintegrated_PCA_prefix}_orig
condor_dataset_subset <- runUMAP(fcd = condor_dataset_subset,
                          input_type = "pca",
                          data_slot = unintegrated_pca_slot,
                          seed = 91)

# Running the UMAP for PCA integrated on TOTAL_BATCH (cohort + experimental batch) using Harmony algorithm:
integrated_pca_slot <- ''
# usually something like {unintegrated_PCA_prefix}_norm
condor_dataset_subset <- runUMAP(fcd = condor_dataset_subset,
                          input_type = "pca",
                          data_slot = integrated_pca_slot,
                          seed = 91,
                          ret_model = T)
# Saving the subsetted cyCondor object
qsave(condor_dataset_subset, saving_path_small)


