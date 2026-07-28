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
library(dplyr)
library(readr)
library(dataPreparation)


######################################################
#### Step 1. Reading CSV files into one dataframe ####
######################################################
# Pre-gate CD45+ live cells in FlowJo
# Export the files as channel values from FlowJo in CSV format
# Store the files in one directory
data_dir <- ""
# Get all CSV file names
csv_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
# Read all CSVs, adding a filename column
merged_data <- csv_files %>%
  lapply(function(file) {
    print(basename(file))
    df <- read_csv(file)
    df$filename <- basename(file)
    return(df)
  }) %>%
  bind_rows()


########################################################################
#### Step 2. Removing 0.01 percentile from each end for each marker ####
########################################################################
# Excluding the FSC / SSC / Filename / ViabilityDye from channel list where winsorisation would happen
cluster.cols <- colnames(merged_data)[!grepl( "FSC|SSC|Time|File|L_D|dead|Comp|iability|filename",colnames(merged_data))]
# Remove the 0.01% top and bottom percentiles
merged_data <- remove_percentile_outlier(data_set = merged_data, cols = cluster.cols, percentile = 0.01)
merged_data_csv_path <- ''
# Writing the merged & filtered CSV table
write_csv(merged_data, merged_data_csv_path)



###############################################################
#### Step 3. Splitting the big table into individual files ####
###############################################################
output_dir <- ""
# Splitting and recording individual files as CSV
split(merged_data, merged_data$filename) %>%
  lapply(function(df) {
    file_name <- unique(df$filename)
    print(file_name)
    df <- select(df, -filename)
    write_csv(df, file.path(output_dir, file_name))
  })

