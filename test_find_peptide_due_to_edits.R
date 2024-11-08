library(biomaRt)
library(ensembldb)
library(GenomicRanges)
library(Biostrings)


# Set up Ensembl connection for the mouse genome
ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl") ## corresponding to library(BSgenome.Mmusculus.UCSC.mm39) # Mouse genome assembly (mm39 2020)


sgRNA <- DNAString("GAAGCAGAAATACCACGACC")

# Query using the mutation coordinates
chromosome <- "7"  # replace with your chromosome
position <- "24985211"  # replace with your position
gene_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description'),
                   filters = c('chromosome_name', 'start', 'end'),
                   values = list(chromosome, position, position),
                   mart = ensembl)

gene_info


## retrieve ref gene seq
gene_id <- gene_info$ensembl_gene_id[1]  # assuming first result is your gene
cds_seq <- biomaRt::getSequence(id = gene_id, type = "ensembl_gene_id", seqType = "cdna", mart = ensembl)

## cds info
cds_info <- biomaRt::getBM(attributes = c('ensembl_transcript_id', 'cds_start', 'cds_end', 'strand'),
                           filters = 'ensembl_gene_id', values = gene_id, mart = ensembl)


## get cdna sequence from transcript id
uniq_transcript_ids <- unique(cds_info$ensembl_transcript_id)
# transcript id #1 as an example
tmp_id <- uniq_transcript_ids[1]
# cdna of transcript
tmp_id_cdna <- biomaRt::getSequence(id = tmp_id, type = "ensembl_transcript_id", 
                                   seqType = "cdna", mart = ensembl)
# genomic sequence of transcript
tmp_id_seq <- biomaRt::getSequence(id = tmp_id, type = "ensembl_transcript_id", 
                                   seqType = "gene_exon_intron", mart = ensembl)

# coding sequence of transcript
tmp_id_coding <- biomaRt::getSequence(id = tmp_id, type = "ensembl_transcript_id", 
                                      seqType = "coding", mart = ensembl)
# translate the coding seq
Biostrings::translate(DNAString(tmp_id_coding$coding))

# peptide sequence of transcript
tmp_id_protein_seq <- biomaRt::getSequence(id = tmp_id, type = "ensembl_transcript_id", 
                                           seqType = "peptide", mart = ensembl)

# find if sgRNA matches in coding sequence
tmp_match <- Biostrings::matchPattern(sgRNA, DNAString(tmp_id_coding$coding))


# View the sequence
tmp_cdna <- DNAString(cds_seq$cdna[14])

matchPattern(sgRNA, tmp_cdna)

Biostrings::translate(tmp_cdna)


### temp

for (idx in 1:length(cds_seq$cdna)) {
  
  tmp_cdna <- DNAString(cds_seq$cdna[idx])
  
  print(as.character(Biostrings::translate(tmp_cdna)))
  
}

tmp_cdna <- DNAString(substr(cds_seq$cdna[12], start = 65, stop = nchar(cds_seq$cdna[12])))
Biostrings::translate(tmp_cdna)

                   