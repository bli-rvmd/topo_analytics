#### this script merges information from VCF and VEP outputs, resulting in gene-transcript-specific quantification of editing and functional impact

#########
# Begin setting runtime parameters
#########

vcf_file <- "/Users/bli/Docker/20241024_ElimBio_Order_864502_TOPO_Smad4_3D/pcr_14_topo_21_M13_F/CRISPResso_on_pcr_14_filtered_seqs_all/VCF_Alleles_frequency_table_around_sgRNA_GATGTGTCATAGACAAGGTG.vcf"

vep_file <- "/Users/bli/Docker/20241024_ElimBio_Order_864502_TOPO_Smad4_3D/pcr_14_topo_21_M13_F/CRISPResso_on_pcr_14_filtered_seqs_all/VEP_Alleles_frequency_table_around_sgRNA_GATGTGTCATAGACAAGGTG.txt"

gene_symbol <- "Smad4"

output_summary_xlsx_file <- "/Users/bli/Docker/20241024_ElimBio_Order_864502_TOPO_Smad4_3D/pcr_14_topo_21_M13_F/CRISPResso_on_pcr_14_filtered_seqs_all/topo_21_sequencing_validation_analysis_report.xlsx"

#########
# End setting runtime parameters
#########

################################ !!! DO NOT CHANGE ANYTHING BELOW !!! ####################################

library(dplyr)
library(tidyr)
library(openxlsx)

## step 0 - read in an VCF file and its associated VEP result output file
# vcf_file <- "./data/temp/VCF_Cellecta_cdkn2a_Alleles_frequency_table_around_sgRNA_GTCGAGCGGCAGGCGACCCC.vcf"
# vep_file <- "./data/temp/Cellecta_cdkn2a_KO_VEP.txt"
# gene_symbol <- "Cdkn2a"

# vcf_file <- "./data/temp/VCF_topo_pcr_60_cic_2D_Alleles_frequency_table_around_sgRNA_GAAGCAGAAATACCACGACC.vcf"
# vep_file <- "./data/temp/cic_2D_VEP_from_VCF.txt"
# gene_symbol <- "Cic"

# vcf_file <- "./data/temp/VCF_topo_pcr_58_cdkn2a_2D_Alleles_frequency_table_around_sgRNA_CGGTGCAGATTCGAACTGCG.vcf"
# vep_file <- "./data/temp/cdkn2a_2D_VEP_from_VCF.txt"
# gene_symbol <- "Cdkn2a"


## step 1 - load files and preprocessing

# vcf
df_vcf_raw <- read.delim(vcf_file, header = T, sep = "\t")
# rename first col
colnames(df_vcf_raw)[1] <- "CHROM"
# split last col into multiple info fields
df_vcf <- df_vcf_raw %>%
  tidyr::separate(INFO, into = c("n_reads", "n_total", "ratio"), sep = ",") %>% 
  dplyr::select(CHROM, # retain informative fields
    ID, 
    n_reads, 
    n_total, 
    ratio
  )

# vep
df_vep_raw <- read.delim(vep_file, header = T, sep = "\t")
# rename first col
colnames(df_vep_raw)[1] <- "ID"
df_vep <- df_vep_raw %>%
  dplyr::filter(tolower(SYMBOL) == tolower(gene_symbol)) %>% # filter out genes that are not of gene_symbol
  dplyr::filter(grepl("protein_coding", BIOTYPE, ignore.case = TRUE)) %>% # filter out transcripts that are not of protein_coding related
  dplyr::select(ID,  # retain informative fields
                Location, 
                Allele, 
                REF_ALLELE, 
                Consequence, 
                IMPACT, 
                SYMBOL, 
                Feature, 
                BIOTYPE, 
                EXON, 
                INTRON, 
                cDNA_position, 
                CDS_position, 
                Protein_position, 
                Amino_acids, 
                Codons)


## step 2 - merge df_vcf and df_vep by ID field
# df <- df_vcf %>%
#   dplyr::left_join(df_vep, by = "ID") %>%
#   dplyr::mutate(ratio = as.numeric(ratio)) %>%
#   dplyr::group_by(Feature) %>% 
#   dplyr::mutate(
#     sum_ratio = paste0(sum(ratio) * 100, "%"),
#     sum_ratio_HIGH = paste0(sum(ratio[IMPACT == "HIGH"]) * 100, "%"), 
#     sum_ratio_HIGH_MODERATE = paste0(sum(ratio[IMPACT %in% c("HIGH", "MODERATE")]) * 100, "%"), 
#     sum_ratio_MODIFIER = paste0(sum(ratio[IMPACT == "MODIFIER"]) * 100, "%")
#   ) 

df <- df_vcf %>%
  dplyr::left_join(df_vep, by = "ID") %>% # joining df_vcf with df_vep
  dplyr::mutate(ratio = as.numeric(ratio)) %>% #
  dplyr::mutate(Percent_Edited = paste0(ratio * 100, "%")) %>% # raname ratio to percent_edited %>%
  dplyr::rename(TranscriptID = Feature) %>% # rename feature to transcriptID
  dplyr::rename(VariantID = ID) %>% # rename ID to variantID
  dplyr::relocate(TranscriptID, VariantID, Percent_Edited, IMPACT, .before = tidyr::everything()) %>% # relocate positions of key columns to very beginning
  dplyr::arrange(TranscriptID) # sort by TranscriptID


df_summary <- df %>%
  dplyr::group_by(TranscriptID) %>%
  dplyr::summarise(
    sum_ratio = paste0(sum(ratio) * 100, "%"), 
    sum_ratio_HIGH = paste0(sum(ratio[IMPACT == "HIGH"]) * 100, "%"),
    sum_ratio_MODERATE = paste0(sum(ratio[IMPACT == "MODERATE"]) * 100, "%"), 
    sum_ratio_HIGH_OR_MODERATE = paste0(sum(ratio[IMPACT %in% c("HIGH", "MODERATE")]) * 100, "%"), 
    sum_ratio_MODIFIER = paste0(sum(ratio[IMPACT == "MODIFIER"]) * 100, "%")
  )

## rename sum_ratio... into %_edited...
colnames(df_summary) <- gsub("^sum\\_ratio", "Percent\\_Edited_Impact", colnames(df_summary))
colnames(df_summary) <- gsub("\\_Impact$", "", colnames(df_summary))

## cap percentage of df_summary
cap_percent <- function(column) {
  # Remove the '%' sign and convert to numeric
  numeric_values <- as.numeric(sub("%", "", column))
  # Cap values at 100
  capped_values <- pmin(numeric_values, 100)
  # Convert back to strings with '%'
  return(paste0(capped_values, "%"))
}

df_summary[-1] <- lapply(df_summary[-1], cap_percent)

## output df and df_summary to sequencing_validation_analysis_report
# create a new workbook
wb <- createWorkbook()

# add worksheets to workbook
addWorksheet(wb, "sequencing_validation_analysis")
addWorksheet(wb, "summary")

# write data to worksheets
writeData(wb, sheet = "sequencing_validation_analysis", df)
writeData(wb, sheet = "summary", df_summary)

# save workbook
saveWorkbook(wb, output_summary_xlsx_file, overwrite = TRUE)

