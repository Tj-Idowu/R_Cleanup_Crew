# Load necessary libraries
library(data.table)
library(zlibbioc)

# Read a single file and return its contents with file name
read_file <- function(file_path) {
  content <- fread(cmd=paste0("zcat ", file_path))
  file_name <- basename(file_path)
  # Add file name as a column
  content[, file := file_name]
    return(content)
}

# Read all files in the directory
read_files <- function(directory) {
  # Get list of files
  files <- list.files(directory, pattern = "\\.gz$", full.names = TRUE)
  contents <- lapply(files, read_file)
  # rbind files
  all_data <- rbindlist(contents)
  return(all_data)
}

# Determine directory path
directory_path <- "Simulated_Metagenomes/Tools_Results/zero/kma"
all_data <- read_files(directory_path)

head(all_data)

# Export dataframe as csv
write.csv(all_data, "all_kma_Simulated_Wastewater.csv", row.names=FALSE)
