#### This script converts aligned sequences of variants from allelic frequency table to VCF format 


library(Biostrings)
library(VariantAnnotation)
library(GenomicRanges)
library(biomaRt)
library(dplyr)

## load GRCm39 mouse genome assembly 2020
library(BSgenome.Mmusculus.UCSC.mm39)
GRCm39 <- BSgenome.Mmusculus.UCSC.mm39

####
# Test reading in an allelic table file

### FIXME! This script only works with single guide RNA for now! 

## procedure of converting a allelic freq table to vcf format
# step 1 - consolidate seqs in the txt file grouping up aligned seqs of same variants
# step 2 - remove seqs of non-contiguous indels / snps that are deemed of subpar topo sequencing quality
# step 3 - diff between aligned_sequence and reference_sequence, build logic for identifying genomic coordinate of gene edit for insertion, deletion, and point mutation, respectively
# step 4 - convert logic from previous step into vcf format, then output
####

## utility functions
# find matches of a seq in GRCm39
find.seq.match.grcm39 <- function(dna_seq) {
  
  # reverse complement of dna_seq
  rev_dna_seq <- Biostrings::reverseComplement(Biostrings::DNAString(dna_seq))
  
  # find matches of reference seq in GRCm39
  res_seq_query <- data.frame()
  
  chrs <- GenomeInfoDb::seqnames(GRCm39)
  
  for (chr in chrs[!grepl("_", chrs)]) {
    
    chr_sequence <- GRCm39[[chr]]
    
    print(paste0(idx, ": ", chr))
    
    # find matches of chr_sequence
    match_positions <- Biostrings::matchPattern(dna_seq, chr_sequence)
    
    if (length(match_positions) > 0) {
      print(paste("Found match on chromosome", chr))
      print(match_positions)
      
      res_seq_query <- rbind(res_seq_query, c(chr, start(match_positions), end(match_positions)))
    }
    
    # find for reverse complement of dna_seq
    match_positions_rev <- Biostrings::matchPattern(rev_dna_seq, chr_sequence)
    
    if (length(match_positions_rev) > 0) {
      print(paste("Found match of reverse complement on chromosome", chr))
      print(match_positions_rev)
      
      res_seq_query <- rbind(res_seq_query, c(chr, start(match_positions_rev), end(match_positions_rev)))
    }
    
  }
  
  colnames(res_seq_query) <- c("chr", "start", "end") 
  
  return (res_seq_query)
  # 
  # # stop if there are multiple matches of a reference sequence
  # if (nrow(res_seq_query) != 1) {
  #   
  #   print(res_seq_query)
  #   stop(paste0("Check row ", idx, " of sequence that it has multiple matches in the GRCm39 genome!\n"))
  # }
  # 
}


# find canonical transcript of a mouse gene
find.mouse.gene.canonical.transcript <- function(gene_symbol) {
  
  # Set up Ensembl connection for the mouse genome
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl") ## corresponding to library(BSgenome.Mmusculus.UCSC.mm39) # Mouse genome assembly (mm39 2020)
  
  ## find canonical transcript of a gene
  results <- getBM(
    attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_transcript_id", "transcript_is_canonical"),
    filters = "external_gene_name",
    values = gene_symbol,
    mart = ensembl
  )
  
  canonical_transcript <- results %>%
    dplyr::filter(transcript_is_canonical == 1) %>%
    dplyr::pull(ensembl_transcript_id)
  
  return (canonical_transcript)
}


## step 0 - read in a test file and QC
# min frequency alleles to keep (by default 0.5%) 
min_perc_reads <- 0.5 


# allelic_freq_table_txt <- "./data/temp/topo_pcr_60_cic_2D_Alleles_frequency_table_around_sgRNA_GAAGCAGAAATACCACGACC.txt" # pcr_60 cic 2D (large delections and point mutation)
# allelic_freq_table_txt <- "./data/temp/topo_pcr_58_cdkn2a_2D_Alleles_frequency_table_around_sgRNA_CGGTGCAGATTCGAACTGCG.txt" # pcr_58 cdkn2a 2D (small deletion and small insertion)
# Cellecta Order# 101888 Smad4 KO
# allelic_freq_table_txt <- "./data/temp/Cellecta_smad4_KO_Alleles_frequency_table_around_sgRNA_TGTCACCATACAGAGAACAT.txt"
# Cellecta Order# 101888 Trp53 KO
# allelic_freq_table_txt <- "./data/temp/Cellecta_trp53_Alleles_frequency_table_around_sgRNA_GACCCTGTCACCGAGACCCC.txt"
# Cellecta Order# 101888 Cdkn2a KO
allelic_freq_table_txt <- "./data/temp/Cellecta_cdkn2a_Alleles_frequency_table_around_sgRNA_GTCGAGCGGCAGGCGACCCC.txt"

# change colnames of allele freq table txt file
df_af <- read.delim(allelic_freq_table_txt, header = T, sep = "\t")
colnames(df_af)[(ncol(df_af) - 1):ncol(df_af)] <- c("n_reads", "perc_reads")

# get sgRNA seq
sgRNA_seq <- gsub(".*_sgRNA_([A-Z]+)\\.txt$", "\\1", allelic_freq_table_txt)

# remove reads that are below min_perc_reads
df_af <- df_af %>%
  dplyr::filter(perc_reads >= min_perc_reads)

## steps 1, 2, 3, 4 - process each row in df_af

df_res <- do.call(rbind, lapply(1:nrow(df_af), function(idx) {
  
  ## define ref and mutated seqs
  mutated <- DNAString(df_af$Aligned_Sequence[idx])
  reference <- DNAString(df_af$Reference_Sequence[idx])
  
  print(paste0("Handling variant ", idx, ": "))
  
  ## detect if there are non-contiguous indels / snps contained in aligned sequence or in reference sequence
  # skip the sequence if true
  if (grepl("-[^-]+-", as.character(mutated)) | grepl("-[^-]+-", as.character(reference))) {
    # next 
    print("Skipping - there are non-contiguous indels / snps contained in the aligned sequence or reference sequence")
    return (c())
  }
  
  ## detect if there are multiple types of mutations occurring in editing region (by default +/- bp of expected cut site)
  ## skip the sequence if true
  if (sum(df_af[idx, c("n_deleted", "n_inserted", "n_mutated")] != 0) != 1) {
    # next
    print("Skipping - unedited seq or there are multiple types of mutations occurring in the editing region")
    return (c())
  }
  
  # number of reads
  n_reads <- df_af[idx, "n_reads"]
  
  ## if there's deletion in mutated seq
  if (df_af[idx, "n_deleted"] != 0) {
    
    del_size <- df_af[idx, "n_deleted"]
    
    # find position of deletion 
    first_index <- which(unlist(strsplit(as.character(mutated), "")) == "-")[1]
    
    # # upstream_seq
    # upstream_seq <- mutated[1:(first_index - 1)]

    # find match of reference seq
    res_seq_query <- find.seq.match.grcm39(reference)
    
    # stop if there are multiple matches of a reference sequence
    if (nrow(res_seq_query) != 1) {
      
      print(res_seq_query)
      stop(paste0("Check row ", idx, " of reference sequence that it has multiple matches in the GRCm39 genome!\n"))
      
    }
    
    # return position of the deletion in VCF format if first_index of deletion in mutated seq is > 1
    if (first_index > 1) {
      
      position <- as.integer(res_seq_query$start) + first_index - 2
      
    }
    
    # return position of the deletion in VCF format if first_index of deletion in mutated seq == 1
    if (first_index == 1) {
      
      n_minus_signs <- sum(unlist(strsplit(as.character(mutated), "")) == "-")
      
      position <- as.integer(res_seq_query$start) - (del_size - n_minus_signs) - 1
      
    }
    
    v_res <- c(
      res_seq_query$chr, #CHROM
      position, # POS
      paste0("del-", del_size), # ID
      as.character(Biostrings::getSeq(GRCm39, # REF by getting seq from GRCm39
                                      res_seq_query$chr, 
                                      position, 
                                      position + del_size)), 
      as.character(Biostrings::getSeq(GRCm39, # ALT by getting seq from GRCm39
                                      res_seq_query$chr, 
                                      position, 
                                      position)), 
      ".", # QUAL
      ".", # FILTER
      n_reads # INFO
    )
    
  }
  
  
  ## if there's an insertion in mutated seq - and there's supposed to be '-' signs in reference seq
  if (df_af[idx, "n_inserted"] != 0) {
    
    ins_size <- df_af[idx, "n_inserted"]
 
    # find location of insertion
    first_index <- which(unlist(strsplit(as.character(reference), "")) == "-")[1]
    
    # find match of ref seq
    cleaned_ref <- gsub("-", "", as.character(reference))
    res_seq_query <- find.seq.match.grcm39(cleaned_ref)
    
    # stop if there are multiple matches of a reference seq
    if (nrow(res_seq_query) != 1) {
      
      print(res_seq_query)
      stop(paste0("Check row ", idx, " of reference sequence that it has multiple matches in the GRCm39 genome!\n"))
      
    }
    
    # return position of the insertion found in VCF format if first_index of deletion in reference seq is > 1
    if (first_index > 1) {
      
      position <- as.integer(res_seq_query$start) + first_index - 2
      
    }
    
    # return position of the insertion found in VCF format if first_index of deletion in reference seq == 1
    if (first_index == 1) {
      
      position <- as.integer(res_seq_start$start) - ins_size - 1
      
    }
    
    v_res <- c(
      res_seq_query$chr, # CHROM
      position, # POS
      paste0("ins+", ins_size), # ID
      as.character(Biostrings::getSeq(GRCm39,  # REF 
                                      res_seq_query$chr, 
                                      position, 
                                      position)), 
      as.character(Biostrings::getSeq(GRCm39,  # ALT 
                                      res_seq_query$chr, 
                                      position, 
                                      position + ins_size)), 
      ".", # QUAL
      ".", # FILTER
      n_reads # INFO
    )
       
  }
  
  
  ## if there's a SNP in mutated seq 
  if (df_af[idx, "n_mutated"] != 0) {
    
    # find index of sgRNA in reference seq
    sg_idx_in_ref <- Biostrings::matchPattern(sgRNA_seq, reference)
    
    n_offset_bp <- 3 # N bps beyond sgRNA both 3' and 5' ends to be considered as quantifying region, set 3 as default
    
    start_idx <- start(sg_idx_in_ref) - 3
    end_idx <- end(sg_idx_in_ref) + 3
    
    # obtain mutated and ref seqs around sgRNA
    sg_mut_seq <- substr(as.character(mutated), start_idx, end_idx)
    sg_ref_seq <- substr(as.character(reference), start_idx, end_idx)
    
    # find indices of mutations while ignoring positions of "N"'s
    sg_mut_chars <- strsplit(sg_mut_seq, "")[[1]]
    sg_ref_chars <- strsplit(sg_ref_seq, "")[[1]]
    
    diff_pos <- which(sg_mut_chars != sg_ref_chars & sg_mut_chars != "N" & sg_ref_chars != "N")
    
    if (length(diff_pos) != 1) {
      
      stop(paste0("Check row ", idx, " of reference sequence that it has multiple SNPs around sgRNA region!\n"))
      
    }
    
    # find coordinate of sg_ref_seq in GRCm39
    res_seq_query <- find.seq.match.grcm39(sg_ref_seq)
    position <- as.integer(res_seq_query$start) + diff_pos - 1
    
    v_res <- c(
      res_seq_query$chr, #CHROM
      position, # POS
      paste0(sg_ref_chars[diff_pos], "->", sg_mut_chars[diff_pos]), # ID
      as.character(Biostrings::getSeq(GRCm39, # REF
                                      res_seq_query$chr, 
                                      position, 
                                      position)), 
      sg_mut_chars[diff_pos], # ALT
      ".", # QUAL
      ".",  # FILTER
      n_reads # INFO
    )
    
  }
 
  return (v_res)
  
}))


## rename df_res
df_res <- as.data.frame(df_res)
colnames(df_res) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

## merge rows of identical mutations
df_res_merged <- df_res %>%
  dplyr::group_by(across(-INFO)) %>%
  dplyr::summarize(INFO = sum(as.integer(INFO)), .groups = "drop")
  

## group up counts of unedited seqs and update INFO field 
df_af_unedited <- df_af %>%
  dplyr::filter(Unedited == "True")

n_unedited_reads <- sum(as.integer(df_af_unedited$n_reads))
n_edited_reads <- sum(as.integer(df_res_merged$INFO))
n_total_reads <- n_unedited_reads + n_edited_reads

df_res_output <- df_res_merged %>%
  dplyr::mutate(
    INFO = paste0(as.integer(INFO), ",", n_total_reads, ",", round(as.integer(INFO) / n_total_reads, 3))
  )

## output df_res_output in VCF format
write.table(df_res_output, 
            file = file.path("./data/temp", sub("\\.txt$", ".vcf", paste0("VCF_", basename(allelic_freq_table_txt)))), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)


##################### TEST BELOW #######################

stop()

## query reference seq
res_query <- data.frame()

for (chr in seqnames(GRCm39)) {
  chr_sequence <- GRCm39[[chr]]  # Get the chromosome sequence
  
  # Find matches (note to remove '-' from reference)
  match_positions <- matchPattern(gsub('-', '', reference), chr_sequence)
  
  if (length(match_positions) > 0) {
    print(paste("Found match on chromosome", chr))
    print(match_positions)
    
    res_query <- rbind(res_query, c(chr, start(match_positions), end(match_positions)))
  }
}

colnames(res_query) <- c("chr", "start", "end")


## find diffs
mismatch_pos <- which(reference != mutated)






