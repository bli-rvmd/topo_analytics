library(dplyr)
library(tools)

####
# run options
####
data_dir <- "./data/temp"
allele_table_txt_file <- "Alleles_frequency_table_around_sgRNA_GAAGCAGAAATACCACGACC.txt"
cols <- c("aligned_seq", "ref_seq", "unedited", "n_del", "n_ins", "n_mut", "n_reads", "perc_reads")

# allele_table_txt_file <- "Alleles_frequency_table.txt"
# cols <- c("aligned_seq", "ref_seq", "ref_name", "read_status", "n_del", "n_ins", "n_mut", "n_reads", "perc_reads")

# allele_table_txt_file <- "TOPO_smad4_Alleles_frequency_table_around_sgRNA_GATGTGTCATAGACAAGGTG.txt"
# cols <- c("aligned_seq", "ref_seq", "unedited", "n_del", "n_ins", "n_mut", "n_reads", "perc_reads") 

# %_reads cutoff
perc_reads_cutoff <- 0.5

# number of duplication
n_dups <- 100


####
# preprocessing
####
df_allele_tbl <- read.table(file.path(data_dir, allele_table_txt_file), header = T, sep = "\t", check.names = F, row.names = NULL)
colnames(df_allele_tbl) <- cols

df_allele_tbl <- df_allele_tbl %>%
  filter(perc_reads >= perc_reads_cutoff) 


####
# conversion
####

fastq <- c()

for (idx in 1:nrow(df_allele_tbl)) {
  
  s <- df_allele_tbl$aligned_seq[idx]
  
  # remove '-' 
  s <- gsub("-", "", s)
  
  for (i in 1:n_dups) {
    
    # generate fastq header line for s
    fastq <- c(fastq, paste0("@AB00001:123:12ABCDEFH:1:1:1:", (idx-1) * n_dups + i, " 1:N:0:AAAAAAAA+GGGGGGGG"))
    
    fastq <- c(fastq, s)
    
    fastq <- c(fastq, "+")
    
    # generate quality score line for s
    tmp <- gsub("N", "#", s)
    tmp <- gsub("[^#]", "I", tmp)
    
    fastq <- c(fastq, tmp)
    
  }
  
}

####
# save output
####
output_file <- paste0("fastq_", file_path_sans_ext(allele_table_txt_file), ".fastq")

writeLines(fastq, con = file.path(data_dir, output_file))


