## 
library(ggplot2)
library(dplyr)
library(RColorBrewer) # For better color palettes

## 
data_dir <- "/Users/bli/Docker/20240919_KPTC_ElimBio_order_861172_TOPO_cdkn2a_cic_nutlin"

list_exp <- list(
  
  "pcr_37" = "Cdkn2a_P13_NoNutlin", 
  "pcr_39" = "Cdkn2a_P13_Nutlin",
  "pcr_41" = "Cdkn2a_P15_Nutlin", 
  "pcr_43" = "Cic_P13_NoNutlin", 
  "pcr_45" = "Cic_P13_Nutlin", 
  "pcr_47" = "Cic_P15_Nutlin"
  
)

## read in allelic table for each pcr run

list_af_tbl <- list()
for (x in names(list_exp)) {
  
  file_x <- list.files(file.path(data_dir, x, paste0("CRISPResso_on_", x, "_filtered_seqs_all")), 
                       "^Alleles_frequency_table.*\\.txt$", full.names = T)
  
  df_x <- read.table(file_x, header = T, sep = "\t", check.names = F, row.names = NULL)
  colnames(df_x) <- c("aligned_seq", "ref_Seq", "unedited", "n_deleted", "n_inserted", "n_mutated", "n_reads", "perc_reads")
  
  # 
  list_af_tbl[[x]] <- df_x
}

## For cdkn2a

dat_pcr_37 <- data.frame(
  edit = c("-61", "+1", "ref", "-398", "-207", "-216", "-7", "-6", "-120"), 
  perc = c(22.5, 35, 10, 5, 5, 5, 5, 7.5, 2.5)
)
dat_pcr_37$target <- "Cdkn2a"
dat_pcr_37$condition <- "P13_No_Nutlin"


dat_pcr_39 <- data.frame(
  edit = c("-398", "-207", "+1", "-61", "-7", "ref", "-8", "-6"), 
  perc = c(54.8, 24.8, 9.5, 2.4, 2.4, 2.4, 2.4, 2.4)
)
dat_pcr_39$target <- "Cdkn2a"
dat_pcr_39$condition <- "P13_Nutlin"


dat_pcr_41 <- data.frame(
  edit = c("-207", "-216", "-6", "-398", "-7", "ref", "+1"), 
  perc = c(40.9, 22.7, 13.6, 4.5, 4.5, 4.5, 4.5)
)
dat_pcr_41$target <- "Cdkn2a"
dat_pcr_41$condition <- "P15_Nutlin"


dat_cdkn2a <- rbind(dat_pcr_37, rbind(dat_pcr_39, dat_pcr_41))

dat_cdkn2a$edit <- as.factor(dat_cdkn2a$edit)

all_combinations_cdkn2a <- expand.grid(
  edit = unique(dat_cdkn2a$edit),
  condition = unique(dat_cdkn2a$condition)
)

# Merge the original data with the full set of combinations
complete_data_cdkn2a <- merge(all_combinations_cdkn2a, dat_cdkn2a, by = c("edit", "condition"), all.x = TRUE)

# Replace missing 'perc' values with 0
complete_data_cdkn2a$perc[is.na(complete_data_cdkn2a$perc)] <- 0

# plot

ggplot(complete_data_cdkn2a, aes(x = condition, y = perc, group = edit, color = edit, shape = edit)) +
  geom_line() +
  geom_point(size = 2.5) +
  theme_bw() +
  labs(x = "Condition", y = "Percentage", color = "Edit", shape = "Edit") +
  scale_color_brewer(palette = "Paired") +  # Use a more discernible color palette
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),   # Increase x-axis text size
        axis.text.y = element_text(size = 12),                          # Increase y-axis text size
        axis.title = element_text(size = 14),                           # Increase axis title size
        legend.title = element_text(size = 13),                         # Increase legend title size
        legend.text = element_text(size = 12)) +
  ggtitle("KPTC Cdkn2a Editing") + 
  scale_shape_manual(values = 1:length(unique(complete_data_cdkn2a$edit))) 


## For Cic

dat_pcr_43 <- data.frame(
  edit = c("-3", "-1", "ref", "-23", "-111", "-4", "-139", "S1", "-5"), 
  perc = c(25.6, 12.8, 30.8, 10.3, 7.7, 5.1, 2.6, 2.6, 2.6)
)
dat_pcr_43$target <- "Cic"
dat_pcr_43$condition <- "P13_No_Nutlin"


dat_pcr_45 <- data.frame(
  edit = c("ref", "-111", "S1", '-3', "-5"), 
  perc = c(81.4, 9.3, 4.7, 2.3, 2.3)
)
dat_pcr_45$target <- "Cic"
dat_pcr_45$condition <- "P13_Nutlin"


dat_pcr_47 <- data.frame(
  edit = c("ref", "-4", "-23", "-2", "-3", "-1", "-40", "-111", "S1"), 
  perc = c(35.5, 19.3, 9.7, 9.7, 6.5, 6.5, 3.2, 3.2, 6.4)
)
dat_pcr_47$target <- "Cic"
dat_pcr_47$condition <- "P15_Nutlin"


dat_cic <- rbind(dat_pcr_43, rbind(dat_pcr_45, dat_pcr_47))

dat_cic$edit <- as.factor(dat_cic$edit)

all_combinations_cic <- expand.grid(
  edit = unique(dat_cic$edit),
  condition = unique(dat_cic$condition)
)

# Merge the original data with the full set of combinations
complete_data_cic <- merge(all_combinations_cic, dat_cic, by = c("edit", "condition"), all.x = TRUE)

# Replace missing 'perc' values with 0
complete_data_cic$perc[is.na(complete_data_cic$perc)] <- 0

# plot

ggplot(complete_data_cic, aes(x = condition, y = perc, group = edit, color = edit, shape = edit)) +
  geom_line() +
  geom_point(size = 2.5) +
  theme_bw() +
  labs(x = "Condition", y = "Percentage", color = "Edit", shape = "Edit") +
  scale_color_brewer(palette = "Paired") +  # Use a more discernible color palette
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),   # Increase x-axis text size
        axis.text.y = element_text(size = 12),                          # Increase y-axis text size
        axis.title = element_text(size = 14),                           # Increase axis title size
        legend.title = element_text(size = 13),                         # Increase legend title size
        legend.text = element_text(size = 12)) +
  ggtitle("KPTC Cic Editing") + 
  scale_shape_manual(values = 1:length(unique(complete_data_cic$edit))) 
