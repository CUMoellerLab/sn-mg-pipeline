#!/Library/Frameworks/R.framework/Versions/Current/Resources/bin/ Rscript

username <- "dds232"
server <- "cbsumoeller.biohpc.cornell.edu"
rsync_fp <- "/workdir/Sprockett/Projects/CU03_Liddell2020/sn-mg-pipeline/output/"

args <- commandArgs(trailingOnly = TRUE)
binning_file <- args[1]

# binning_file <- "/Users/danielsprockett/Software/sn-mg-pipeline/resources/scripts/binning_test_subset_2.txt"

binning_df <- read.table(file = binning_file, sep = "\t", header = TRUE)
sample_list <- binning_df$Sample
n_sample <- length(sample_list)

contig_indexes <- which(binning_df$Contig_Groups == "A")
contig_fps_from <- paste0("/workdir/Sprockett/Projects/CU03_Liddell2020/sn-mg-pipeline/", binning_df[contig_indexes, "Contigs"])

read_indexes <- which(binning_df$Read_Groups == "A")
read_sample_ids <- binning_df[read_indexes, "Sample"]
read_R1_fps <- paste0( "/workdir/Sprockett/Projects/CU03_Liddell2020/sn-mg-pipeline/output/qc/host_filter/nonhost/", read_sample_ids, ".R1.fastq.gz")
read_R2_fps <- paste0( "/workdir/Sprockett/Projects/CU03_Liddell2020/sn-mg-pipeline/output/qc/host_filter/nonhost/", read_sample_ids, ".R2.fastq.gz")
read_fps_from <- c(rbind(read_R1_fps, read_R2_fps))
read_fps_to <- "sn-mg-pipeline/output/qc/host_filter/nonhost/"

message(n_sample, " samples found in binning file ", basename(binning_file), ".")
message("    of those, ", length(read_indexes), " are in the 'read' group and ", length(contig_indexes), " are in the 'contig' group.")

dir.create(read_fps_to, recursive = TRUE)
message("Transferring ", length(read_fps_from), " R1 and R2 fastq files to ", read_fps_to , " using scp.")

read_tar_fp <- "/workdir/Sprockett/Projects/CU03_Liddell2020/sn-mg-pipeline/output/qc/host_filter/nonhost/prototype_reads_for_binning.tar.gz"
scp_cmd <- paste0("scp ", username, "@", server, ":", read_tar_fp, " ", read_fps_to)
system(scp_cmd)
untar_cmd <- "tar -xvf sn-mg-pipeline/output/qc/host_filter/nonhost/prototype_reads_for_binning.tar.gz -C sn-mg-pipeline/output/qc/host_filter/nonhost/"
system(untar_cmd)

contig_fps_to <- "sn-mg-pipeline/output/assemble/megahit/"
dir.create(contig_fps_to, recursive = TRUE)
message("Transferring ", length(contig_fps_from), " fasta contig files to ", contig_fps_to , " using scp.")

for (i in seq_along(contig_fps_from)) {
  i_file <- contig_fps_from[i]
  i_cmd <- paste0("scp ", username, "@", server, ":", i_file, " ", contig_fps_to)
  system(i_cmd)
}
file_list <- list.files(contig_fps_to)
if (length(contig_fps_from) == length(file_list)) {
  message("All ", length(contig_fps_from), " files downloaded correctly.")
}

message("Copying ", basename(binning_file), " to sn-mg-pipeline/resources/config/" )
binning_cp_cmd <- paste0("cp ", basename(binning_file), " sn-mg-pipeline/resources/config/binning.txt")
system(binning_cp_cmd)

setwd("sn-mg-pipeline")

message("Begin running sn-mg-pipeline module Snakefile-bin.")
snakemake_cmd <- "snakemake -s Snakefile-bin -c all --use-conda --conda-prefix ~/snakemake_envs -k"
system(snakemake_cmd)

message("Begin using rsync to transfer output/ to ", rsync_fp)
rsync_cmd <- paste0("rsync -a output/ ", server, ":", rsync_fp)
system(rsync_cmd)

endres_cmd <- "/programs/bin/labutils/endres.pl"
system(endres_cmd)
