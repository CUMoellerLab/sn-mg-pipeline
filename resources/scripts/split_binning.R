#!/Library/Frameworks/R.framework/Versions/Current/Resources/bin/ Rscript

args <- commandArgs(trailingOnly = TRUE)
n_groups <- as.numeric(args[1])
binning_file <- args[2]

if (n_groups < 1) {
  message("'n_groups' parameter must be greater than 1.")
}

# n_groups <- 10
# binning_file <- "/Users/danielsprockett/Software/sn-mg-pipeline/resources/scripts/binning_test.txt"

binning_df <- read.table(file = binning_file, sep = "\t", header = TRUE)
sample_list <- binning_df$Sample
n_sample <- length(sample_list)
message(n_sample, " samples found in sample file ", basename(binning_file), ".")

Read_Groups_index <- which(binning_df$Read_Groups == "A")
n_Read_Groups <- length(Read_Groups_index)
Read_Groups_names <- binning_df[Read_Groups_index, "Sample"]
message(n_Read_Groups, " of those samples are in the 'read' binning group, and will be mapped against every sample. ")
message(paste(Read_Groups_names, "\n"))

message(paste0("Splitting file ", basename(binning_file), " into ", n_groups,  " even groups."))

binning_file_name <- paste0(tools::file_path_sans_ext(binning_file), "_subset_", 1, ".txt")
message(basename(binning_file_name), " contains the 'read' files mapped against their own contigs.")
reads_df <- binning_df[Read_Groups_index, , drop = FALSE]
write.table(reads_df, file = binning_file_name, quote = FALSE, sep = '\t', row.names = FALSE)
reads_mqpping_df <- reads_df
reads_mqpping_df$Contig_Groups <- ""

message("Subsequent files contains those same 'read' files mapped against all other files.")

contig_names <- sample_list[!sample_list %in% Read_Groups_names]
groups_list <- split(contig_names, sort(seq_along(contig_names) %% (n_groups - 1)))

for (i in seq_along(groups_list)) {
  igroup <- groups_list[i]
  
  isample <- binning_df[binning_df$Sample %in% igroup[[1]], , drop = FALSE]
  isample_reads <- rbind(reads_mqpping_df, isample)
  binning_file_name <- paste0(tools::file_path_sans_ext(binning_file), "_subset_", (i + 1), ".txt")
  write.table(isample_reads, file = binning_file_name, quote = FALSE, sep = '\t', row.names = FALSE)

}


