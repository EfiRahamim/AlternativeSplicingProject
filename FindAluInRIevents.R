library(argparse)
# CLI Arguments
parser <- ArgumentParser(description='Filter RI events from PSI-Sigma filtered results and look for Alu within them.\n')
parser$add_argument("-i", action="store", dest="PSI_Sigma_file", help="PSI-Sigma filtered results file (.csv)")
parser$add_argument("-o", action="store", dest="output_dir", help="Output directory. Files of RI events and Alu regions would be created.")
parser$add_argument("-l1", action="store", dest="lable1", help="Group 1 lable")
parser$add_argument("-l1", action="store", dest="lable2", help="Group 2 lable")
user_args <- parser$parse_args()
stopifnot(!is.null(user_args$PSI_Sigma_file) && !is.null(user_args$lable1) && !is.null(user_args$lable2))
print(user_args)

# Arguments assignment
PSI_Sigma_file <- user_args$PSI_Sigma_file
output_dir <- user_args$output_dir
lable1 <- user_args$lable1
lable2 <- user_args$lable2

# DEBUG Arguments
PSI_Sigma_file <- "/private10/Projects/Efi/CRG/SF3B1_WT/PSI-Sigma/Melanoma/GencodeGTF/Results/SplicingEventsFiltered-DMSO_vs_dCEMM1-PSI20_Pvalue0.05_FDR0.05.csv"
output_dir <- "/private10/Projects/Efi/CRG/SF3B1_WT/RI_EditingIndex/Melanoma/DMSO_vs_dCEMM1/"
lable1 <- "DMSO"
lable2 <- "dCEMM1"


# read results file
results <- read.csv(PSI_Sigma_file)
# creat BED file of Intron Retention events for each group
groupA_ri_exons <- results[(results$Event.Type=='IR' | results$Event.Type=='IR (overlapping region)') & results$ΔPSI.... < 0, "Target.Exon"]
groupB_ri_exons <- results[(results$Event.Type=='IR' | results$Event.Type=='IR (overlapping region)') & results$ΔPSI.... > 0, "Target.Exon"]
# split by : and -
groupA_BED_ri_exons <- strsplit(groupA_ri_exons, "[:-]")
groupB_BED_ri_exons <- strsplit(groupB_ri_exons, "[:-]")
# combine character vector into matrix
groupA_BED_ri_exons_matrix <- do.call(rbind, groupA_BED_ri_exons)
groupB_BED_ri_exons_matrix <- do.call(rbind, groupB_BED_ri_exons)
# conver to data frame and remove duplicates
groupA_BED_ri_exons_df <- unique(as.data.frame(groupA_BED_ri_exons_matrix))
groupB_BED_ri_exons_df <- unique(as.data.frame(groupB_BED_ri_exons_matrix))
# change colnames
colnames(groupA_BED_ri_exons_df) <- c("chr", "start", "end")
colnames(groupB_BED_ri_exons_df) <- c("chr", "start", "end")

groupA_BED_ri_exons_df$start <- as.numeric(as.character(groupA_BED_ri_exons_df$start))
groupA_BED_ri_exons_df$end <- as.numeric(as.character(groupA_BED_ri_exons_df$end))
groupB_BED_ri_exons_df$start <- as.numeric(as.character(groupB_BED_ri_exons_df$start))
groupB_BED_ri_exons_df$end <- as.numeric(as.character(groupB_BED_ri_exons_df$end))

sorted_groupA_BED_ri_exons_df <- groupA_BED_ri_exons_df[order(groupA_BED_ri_exons_df$chr, groupA_BED_ri_exons_df$start, groupA_BED_ri_exons_df$end), ]
write.table(sorted_BED_ri_exons_df, file = file.path(output_dir,"IR_events.sorted.bed"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE) # Write the sorted matrix into a BED file
