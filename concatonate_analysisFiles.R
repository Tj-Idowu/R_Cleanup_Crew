# Load necessary libraries
library(data.table)
#library(pbapply) # for progress bar during file reading
library(zlibbioc) # for reading .gz files

# Function to read a single file and return its contents with file name
read_file <- function(file_path) {
  # Read file
  content <- fread(cmd=paste0("zcat ", file_path))
  
  # Extract file name
  file_name <- basename(file_path)
  
  # Add file name as a column
  content[, file := file_name]
  
  # Return content
  return(content)
}

# Function to read all files in a directory
read_files <- function(directory) {
  # Get list of files
  files <- list.files(directory, pattern = "\\.gz$", full.names = TRUE)
  
  # Read files in parallel with progress bar
  contents <- lapply(files, read_file)
  
  # Concatenate contents into a single dataframe
  all_data <- rbindlist(contents)
  
  return(all_data)
}

# Replace 'directory_path' with the path to your folder containing .gz files
directory_path <- "/scratch/12355656/Benchmarking/Simulated_Metagenomes/Tools_Results/zero/kma"
all_data <- read_files(directory_path)

# Print the first few rows of the combined dataframe
print(head(all_data))

# Export dataframe as csv
write.csv(all_data, "/scratch/12355656/Benchmarking/Simulated_Metagenomes/results/all_kma_Simulated_Wastewater.csv", row.names=FALSE)
