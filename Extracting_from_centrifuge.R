# Load libraries
library(readr)
library(stringr)
library(tidyr)
library(dplyr)
library(data.table)

# Set working directory
setwd("/Benchmarking/Simulated_Metagenomes/results/centrifuge")

# Load in .csv sample files
F3<- read_csv("Bacillus_F3_classification_centrifuge.csv")

# Print the first few rows of the dataframe
print(head(F3))

# Remove ".tsv" from file column
F3$file<- gsub('.tsv','', F3$file)

# Pathogen sequence
# To select all the data based on one sequence source from one column. Chromosome or plasmid of bacteria but not whole bacteria genome.
# Selected sample sequences
anthracis_pl1<-F3[grep("NC_007530.2*", F3$Read_ID), ]

# To select all the data for the whole bacterial genome with no regard for the reference.
# All sample sequences detected
all_anthracis<- F3[grep("^NC_007530.2*|^NC_007322.2*|^NC_007323.2*", F3$readID), ]

# To select all the data that mapped to the specific reference and count them
# All positives
AP <- F3[grep("1392", F3$taxID), ]
AP_count<- as.data.frame(table(AP$file))
colnames(AP_count)<- c("File_name", "Read_count")
AP_count$File_name<- as.character(AP_count$File_name)

# To select all the data from "All sample" that is mapped to the correct reference
# True positives
TP <- all_anthracis[grep("1392", F3$taxID), ]
TP_count<- as.data.frame(table(TP$file))
colnames(TP_count)<- c("File_name", "Read_count")
TP_count$File_name<- as.character(TP_count$File_name)

# To select all the data from "All sample" that is not mapped to the reference
# False negative
FN <- all_anthracis[!grepl("1392", all_anthracis$taxID), ]
FN_count<- as.data.frame(table(FN$file))
colnames(FN_count)<- c("File_name", "Read_count")
FN_count$File_name<- as.character(FN_count$File_name)

# To select all the reads that map to the reference but are not the sample sequence
# False positive
FP<- AP[!grepl("^NC_007530.2*|^NC_007322.2*|^NC_007323.2*", AP$readID), ]
FP_count<- as.data.frame(table(FP$file))
colnames(FP_count)<- c("File_name", "Read_count")
FP_count$File_name<- as.character(FP_count$File_name)

# Decoy reads that map to the pathogen reference
D_mapping<- AP[grep("^NC_005957.1*|^NC_006578.1*", AP$readID), ]
D_count<- as.data.frame(table(D_mapping$file))
colnames(D_count)<- c("File_name", "Read_count")
D_count$File_name<- as.character(D_count$File_name)

# Change table format for all positive data
# Define empty dataframe
AP_table<- data.frame(
  row_names = c("r1-0","r1000-1","r100-1","r10-1","r1-1","r1-10","r1-100","r1-1000","r0-1"),
  stringsAsFactors = FALSE)

# Extract column identifiers from filled_df's file_name column
col_names_filled <- as.numeric(sapply(strsplit(AP_count$File_name, "_"), "[[", 5))

# Add column names to the empty dataframe
AP_table <- cbind(AP_table, matrix(NA, nrow = nrow(AP_table), ncol = max(col_names_filled)))

# Rename the columns of the empty dataframe
colnames(AP_table)[2:ncol(AP_table)] <- as.character(1:max(col_names_filled))

# Extract row identifiers from AP_count's file_name column
row_names_filled <- sapply(strsplit(AP_count$File_name, "_"), "[[", 4)

# Match row and column identifiers with AP_table's row and column names
matching_rows <- match(row_names_filled, AP_table$row_names)
matching_cols <- match(col_names_filled, colnames(AP_table)[-1])

# Fill in the values from filled_df into the corresponding cells of AP_table
for (i in 1:length(matching_rows)) {
  if (!is.na(matching_rows[i]) && !is.na(matching_cols[i])) {
    AP_table[matching_rows[i], matching_cols[i] + 1] <- AP_count$Read_count[i]
  }
}
colnames(AP_table)[1]<- "Conditions"
AP_table[is.na(AP_table)] <- 0

# Write out into a .csv file
write_csv(AP_table, "Bacillus_F3_centrifuge_allPositives.csv", row.names=FALSE)

#############################################

# Change table format for all positive data for the true positives
# Make sure that all the variables have been changed to the appropriate ones
TP_table<- data.frame(
  row_names = c("r1-0","r1000-1","r100-1","r10-1","r1-1","r1-10","r1-100","r1-1000","r0-1"),
  stringsAsFactors = FALSE)
col_names_filled <- as.numeric(sapply(strsplit(TP_count$File_name, "_"), "[[", 5))
TP_table <- cbind(TP_table, matrix(NA, nrow = nrow(TP_table), ncol = max(col_names_filled)))
colnames(TP_table)[2:ncol(TP_table)] <- as.character(1:max(col_names_filled))
row_names_filled <- sapply(strsplit(TP_count$File_name, "_"), "[[", 4)
matching_rows <- match(row_names_filled, TP_table$row_names)
matching_cols <- match(col_names_filled, colnames(TP_table)[-1])
for (i in 1:length(matching_rows)) {
  if (!is.na(matching_rows[i]) && !is.na(matching_cols[i])) {
    TP_table[matching_rows[i], matching_cols[i] + 1] <- TP_count$Read_count[i]
  }
}
colnames(TP_table)[1]<- "Conditions"
TP_table[is.na(TP_table)] <- 0
write_csv(TP_table, "Bacillus_F3_centrifuge_truePositives.csv", row.names=FALSE)

# False positive table
FP_table<- data.frame(
  row_names = c("r1-0","r1000-1","r100-1","r10-1","r1-1","r1-10","r1-100","r1-1000","r0-1"),
  stringsAsFactors = FALSE)
col_names_filled <- as.numeric(sapply(strsplit(FP_count$File_name, "_"), "[[", 5))
FP_table <- cbind(FP_table, matrix(NA, nrow = nrow(FP_table), ncol = max(col_names_filled)))
colnames(FP_table)[2:ncol(FP_table)] <- as.character(1:max(col_names_filled))
row_names_filled <- sapply(strsplit(FP_count$File_name, "_"), "[[", 4)
matching_rows <- match(row_names_filled, FP_table$row_names)
matching_cols <- match(col_names_filled, colnames(FP_table)[-1])
for (i in 1:length(matching_rows)) {
  if (!is.na(matching_rows[i]) && !is.na(matching_cols[i])) {
    FP_table[matching_rows[i], matching_cols[i] + 1] <- FP_count$Read_count[i]
  }
}
colnames(FP_table)[1]<- "Conditions"
FP_table[is.na(FP_table)] <- 0
write_csv(FP_table, "Bacillus_F3_centrifuge_falsePositives.csv", row.names=FALSE)

# False negative table
FN_table<- data.frame(
  row_names = c("r1-0","r1000-1","r100-1","r10-1","r1-1","r1-10","r1-100","r1-1000","r0-1"),
  stringsAsFactors = FALSE)
col_names_filled <- as.numeric(sapply(strsplit(FN_count$File_name, "_"), "[[", 5))
FN_table <- cbind(FN_table, matrix(NA, nrow = nrow(FN_table), ncol = max(col_names_filled)))
colnames(FN_table)[2:ncol(FN_table)] <- as.character(1:max(col_names_filled))
row_names_filled <- sapply(strsplit(FN_count$File_name, "_"), "[[", 4)
matching_rows <- match(row_names_filled, FN_table$row_names)
matching_cols <- match(col_names_filled, colnames(FN_table)[-1])
for (i in 1:length(matching_rows)) {
  if (!is.na(matching_rows[i]) && !is.na(matching_cols[i])) {
    FN_table[matching_rows[i], matching_cols[i] + 1] <- FN_count$Read_count[i]
  }
}
colnames(FN_table)[1]<- "Conditions"
FN_table[is.na(FN_table)] <- 0
write_csv(FN_table, "Bacillus_F3_centrifuge_falseNegatives.csv", row.names=FALSE)

# Decoy mapping
Decoy_table<- data.frame(
  row_names = c("r1-0","r1000-1","r100-1","r10-1","r1-1","r1-10","r1-100","r1-1000","r0-1"),
  stringsAsFactors = FALSE)
col_names_filled <- as.numeric(sapply(strsplit(D_count$File_name, "_"), "[[", 5))
Decoy_table <- cbind(Decoy_table, matrix(NA, nrow = nrow(Decoy_table), ncol = max(col_names_filled)))
colnames(Decoy_table)[2:ncol(Decoy_table)] <- as.character(1:max(col_names_filled))
row_names_filled <- sapply(strsplit(D_count$File_name, "_"), "[[", 4)
matching_rows <- match(row_names_filled, Decoy_table$row_names)
matching_cols <- match(col_names_filled, colnames(Decoy_table)[-1])
for (i in 1:length(matching_rows)) {
  if (!is.na(matching_rows[i]) && !is.na(matching_cols[i])) {
    Decoy_table[matching_rows[i], matching_cols[i] + 1] <- D_count$Read_count[i]
  }
}
colnames(Decoy_table)[1]<- "Conditions"
Decoy_table[is.na(Decoy_table)] <- 0
write_csv(Decoy_table, "Bacillus_F3_centrifuge_DecoytoPathRef.csv", row.names=FALSE)
