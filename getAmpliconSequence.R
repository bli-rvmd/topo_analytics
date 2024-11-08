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
mouse_genome <- BSgenome.Mmusculus.UCSC.mm39

# Specify the genomic coordinates 
# cdkn2a primers 17 and 18
# chromosome <- "chr4"   # Specify the chromosome
# start_pos <- 89276463  # Starting position
# end_pos <- 89277265     # Ending position

# cic primers 1 and 2
chromosome <- "chr7"
start_pos <- 24985221
end_pos <- 24985332

# Retrieve the reference sequence for the specified region
reference_sequence <- getSeq(mouse_genome, names = chromosome, start = start_pos, end = end_pos)

# Display the reference sequence
print(as.character(reference_sequence))

# Optionally, you can write the sequence to a file
# writeXStringSet(reference_sequence, filepath = "mouse_reference_sequence.fasta")


#########
query_seq <- DNAString("GGCCAGCATGGAGTAAGACC")
hits <- vmatchPattern(query_seq, mouse_genome)



#########################
library(VariantAnnotation)

chain <- import.chain("mm10ToMm39.over.chain.gz")
