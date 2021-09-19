#!/Library/Frameworks/R.framework/Versions/Current/Resources/bin/ Rscript

args <- commandArgs(trailingOnly = TRUE)
n_groups <- as.numeric(args[1])
sample_file <- args[2]
unit_file <- args[3]

# n_groups <- 10
# sample_file <- "/Users/danielsprockett/Software/sn-mg-pipeline/resources/scripts/sample_test.txt"
# unit_file <- "/Users/danielsprockett/Software/sn-mg-pipeline/resources/scripts/unit_test.txt"

sample_df <- read.table(file = sample_file, header = TRUE)
sample_list <- sample_df$Sample
n_sample <- length(sample_list)
message(n_sample, " samples found in sample file ", basename(sample_file), ".")

unit_df <- read.table(file = unit_file, header = TRUE)
n_unit <- length(unit_df$Sample)
message(n_unit, " samples found in unit file ", basename(unit_file),  ".")

if (n_sample > n_unit) {
  warning("Not all samples in " , basename(sample_file), " were not found in ", basename(unit_file))
}

if (n_sample < n_unit) {
  warning(basename(unit_file), " is being subset to only the samples found in " , basename(sample_file), ".")
  units_df <- units_df[units_df$Sample %in% sample_list, ]
}

message(paste0("Splitting files ", basename(sample_file), " and " , basename(unit_file), " into ", n_groups,  " even groups."))

groups_list <- split(sample_list, sort(seq_along(sample_list) %% n_groups))

for (i in seq_along(groups_list)) {
  igroup <- groups_list[i]
  
  isample <- sample_df[sample_df$Sample %in% igroup[[1]], , drop = FALSE]
  sample_file_name <- paste0(tools::file_path_sans_ext(sample_file), "_subset_", i, ".txt")
  write.table(isample, file = sample_file_name, quote = FALSE, sep = '\t', row.names = FALSE)
  
  iunit <- unit_df[unit_df$Sample %in% igroup[[1]], , drop = FALSE]
  unit_file_name <- paste0(tools::file_path_sans_ext(unit_file), "_subset_", i, ".txt")
  write.table(iunit, file = unit_file_name, quote = FALSE, sep = '\t', row.names = FALSE)

}


