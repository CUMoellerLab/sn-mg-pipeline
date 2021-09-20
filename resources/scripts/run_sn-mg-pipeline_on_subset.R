#!/Library/Frameworks/R.framework/Versions/Current/Resources/bin/ Rscript

username <- "dds232"
server <- "cbsumoeller.biohpc.cornell.edu"
rsync_fp <- "/workdir/Sprockett/Projects/CU03_Liddell2020/sn-mg-pipeline/output/"

args <- commandArgs(trailingOnly = TRUE)
sample_file <- args[1]
unit_file <- args[2]

sample_df <- read.table(file = sample_file, header = TRUE)
sample_list <- sample_df$Sample
n_sample <- length(sample_list)
message(n_sample, " samples found in sample file ", basename(sample_file), ".")

unit_df <- read.table(file = unit_file, header = TRUE)
unit_list <- unit_df$Sample
n_unit <- length(unit_df$Sample)
message(n_unit, " samples found in unit file ", basename(unit_file),  ".")

sample_to_keep <- unit_list[unit_list %in% sample_list]

if (length(sample_to_keep) == 0) {
  warning("No samples in " , basename(sample_file), " were not found in ", basename(unit_file), ". Please check your files.")
}

if (n_sample > n_unit) {
  warning("Not all samples in " , basename(sample_file), " were not found in ", basename(unit_file))
}

if (n_sample < n_unit) {
  warning(basename(unit_file), " is being subset to only the samples found in " , basename(sample_file), ".")
  unit_df <- unit_df[unit_df$Sample %in% sample_list, ]
}

dir.create("sn-mg-pipeline/raw_data/")
scp_list <- c(rbind(unit_df$R1, unit_df$R2))
message("Transferring ", length(scp_list), " fastq files to sn-mg-pipeline/raw_data/ using scp.")

for (i in seq_along(scp_list)) {
  i_file <- scp_list[i]
  i_cmd <- paste0("scp ", username, "@", server, ":", i_file, " sn-mg-pipeline/raw_data/")
  system(i_cmd)
}
file_list <- list.files("sn-mg-pipeline/raw_data/")
if (length(scp_list) == length(file_list)) {
  message("All ", length(scp_list), " files downloaded correctly.")
}

message("Moving samples.txt and updated units.txt to sn-mg-pipeline/resources/config/")
sample_cp_cmd <- paste0("cp ", sample_file, " sn-mg-pipeline/resources/config/samples.txt")
system(sample_cp_cmd)

unit_R1 <- paste0( "raw_data/", basename(unit_df$R1))
unit_R2 <- paste0( "raw_data/", basename(unit_df$R2))

updated_units_df <- data.frame(Sample = unit_df$Sample,
                               Unit = unit_df$Unit,
                               R1 = unit_R1,
                               R2 = unit_R2)

write.table(updated_units_df, file = "sn-mg-pipeline/resources/config/units.txt", quote = FALSE, sep = '\t', row.names = FALSE)

setwd("sn-mg-pipeline")

message("Begin running sn-mg-pipeline.")
snakemake_cmd <- "snakemake -c all --use-conda --conda-prefix ~/snakemake_envs -k"
system(snakemake_cmd)

message("Begin using rsync to transfer output/ to ", rsync_fp)
rsync_cmd <- paste0("rsync -a output/ ", server, ":", rsync_fp)
system(rsync_cmd)

endres_cmd <- "/programs/bin/labutils/endres.pl"
system(endres_cmd)


