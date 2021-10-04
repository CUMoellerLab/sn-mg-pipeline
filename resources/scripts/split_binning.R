#!/Library/Frameworks/R.framework/Versions/Current/Resources/bin/ Rscript

args <- commandArgs(trailingOnly = TRUE)
n_groups <- as.numeric(args[1])
binning_file <- args[2]

# n_groups <- 10
# binning_file <- "/Users/danielsprockett/Software/sn-mg-pipeline/resources/scripts/binning_test.txt"

binning_df <- read.table(file = binning_file, header = TRUE)
sample_list <- binning_df$Sample
n_sample <- length(sample_list)
message(n_sample, " samples found in sample file ", basename(sample_file), ".")

message(paste0("Splitting file ", basename(binning_file), " into ", n_groups,  " even groups."))

groups_list <- split(sample_list, sort(seq_along(sample_list) %% n_groups))

for (i in seq_along(groups_list)) {
  igroup <- groups_list[i]
  
  isample <- sample_df[sample_df$Sample %in% igroup[[1]], , drop = FALSE]
  sample_file_name <- paste0(tools::file_path_sans_ext(sample_file), "_subset_", i, ".txt")
  write.table(isample, file = sample_file_name, quote = FALSE, sep = '\t', row.names = FALSE)

}


