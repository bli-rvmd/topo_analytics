#################
# Begin setting runtime parameters
#################

fwd_primer <- toupper("acctccagcgtattctggtagtccagagtggatg")
rev_primer <- toupper("gcactcttaacagctgagccacgtct")

#################
# End setting runtime parameters
#################

# Install the necessary packages if not already installed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# # Install BSgenome and mouse genome data package
# BiocManager::install("BSgenome")
# BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

# Load the necessary libraries
library(BSgenome)
# library(BSgenome.Mmusculus.UCSC.mm10) # Mouse genome assembly (mm10)
library(BSgenome.Mmusculus.UCSC.mm39) # Mouse genome assembly (mm39 2020)
library(Biostrings) 
library(rtracklayer)

# Define the mouse genome object
GRCm39 <- BSgenome.Mmusculus.UCSC.mm39

# utility function
# find matches of a seq in GRCm39
find.seq.match.grcm39 <- function(dna_seq) {
  
  # reverse complement of dna_seq
  rev_dna_seq <- Biostrings::reverseComplement(Biostrings::DNAString(dna_seq))
  
  # find matches of reference seq in GRCm39
  res_seq_query <- data.frame()
  
  chrs <- GenomeInfoDb::seqnames(GRCm39)
  
  for (chr in chrs[!grepl("_", chrs)]) {
    
    chr_sequence <- GRCm39[[chr]]
    
    # print(paste0(idx, ": ", chr))
    print(chr)
    
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

# find matches of fwd and reverse primers
fwd_primer_query <- find.seq.match.grcm39(fwd_primer)
rev_primer_query <- find.seq.match.grcm39(rev_primer)

chr_info <- unique(fwd_primer_query$chr, rev_primer_query$chr)
start_info <- min(as.numeric(c(fwd_primer_query[c("start", "end")], rev_primer_query[c("start", "end")])))
end_info <- max(as.numeric(c(fwd_primer_query[c("start", "end")], rev_primer_query[c("start", "end")])))

amplicon_sequence <- Biostrings::getSeq(GRCm39, names = chr_info, start = start_info, end = end_info)


# Display the reference sequence
print(as.character(amplicon_sequence))
