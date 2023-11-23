library(argparse)
# CLI Arguments
parser <- ArgumentParser(description='Create PCA plot of splicing analysis, based on PSI scores. Works with PSI-Sigma tool results.')
parser$add_argument("-i", action="store", dest="PSI_Sigma_dir", help="Super directory of PSI-Sigma result (Where sub-directories of comparisons are located)")
parser$add_argument("-dir_name", action="store", dest="output_dir_name", help="The name of the actual directories of the results")
parser$add_argument("-file_name", action="store", dest="output_file_name", help="The name of PSI-Sigma results file. default: 'PSI-Sigma_r10_ir3.sorted.txt' ", default="PSI-Sigma_r10_ir3.sorted.txt")
parser$add_argument("-info", action="store", dest="group_info_file", help="Group information file in tab-delimetered format. Must have heared: 'Sample\tGroup'. Order of samples must much the order specified for the PSI-Sigma tool in the bam.txt files.", default="PSI-Sigma_r10_ir3.sorted.txt")
user_args <- parser$parse_args()
stopifnot(!is.null(user_args$PSI_Sigma_dir) && !is.null(user_args$output_dir) && !is.null(user_args$output_dir_name) && !is.null(user_args$group_info_file))
print(user_args)

## DEBUG Arguments
group_info_file <-"/private10/Projects/Efi/CRG/SF3B1_WT/PSI-Sigma/UM/Input/GroupInfo.txt"
PSI_Sigma_dir <- "/private10/Projects/Efi/CRG/SF3B1_WT/PSI-Sigma/UM/"
output_dir_name <- "Output_copy"
output_file_name <- "PSI-Sigma_r10_ir3.sorted.txt"

## Arguments assignment
group_info_file <- user_args$group_info_file
PSI_Sigma_dir <- user_args$PSI_Sigma_dir
output_dir_name <- user_args$output_dir_name
output_file_name <- user_args$output_file_name

library(dplyr)
library(tidyverse)

group_info <- read.csv(group_info_file, sep='\t')
cols_to_keep <- c("Event.Region", "Gene.Symbol", "Target.Exon", "Event.Type", "Reference.Transcript")
psi_df <- list()
output_dirs <- dir(PSI_Sigma_dir, recursive = T, pattern = output_dir_name, full.names = T, include.dirs = TRUE)
for (dir in output_dirs){
  comparison <- basename(dirname((dir)))
  group_a <- str_match(comparison, "(\\w+)_vs_(\\w+)")[2]
  group_b <- str_match(comparison, "(\\w+)_vs_(\\w+)")[3]
  result_file <- list.files(dir, pattern= output_file_name, full.names=T)
  df <- read.csv(result_file, sep='\t')
  splited_n <- as.data.frame(do.call(rbind,strsplit(df$N.Values, "|", fixed = T)))
  splited_t <- as.data.frame(do.call(rbind,strsplit(df$T.Values, "|", fixed = T)))
  splited <- cbind(splited_n, splited_t)
  splited[splited=="na"] <- 0
  splited <- as.data.frame(sapply(splited, as.numeric))
  samples_a <- subset(group_info, Group==group_a)$Sample
  samples_b <- subset(group_info, Group==group_b)$Sample
  colnames(splited) <- c(samples_a, samples_b)
  psi_df[[comparison]] <- cbind(df[,cols_to_keep], splited)
}
psi_unified <- psi_df %>% reduce(inner_join, by=cols_to_keep)
t_psi_unified <- t(psi_unified[,-c(1:5)])
t_psi_unified <- t_psi_unified[, apply(t_psi_unified, 2, function(x) sd(x) != 0)]
# Check if there are any columns left
if (ncol(t_psi_unified) == 0) {
  stop("No columns with non-zero variance. Unable to perform PCA.")
}
pca_result <- prcomp(unique(t_psi_unified), scale=T, center=T)
pca_data <- as.data.frame(pca_result$x)
pca_data <- unique(pca_data)
rownames(pca_data) <- sub("\\.x|\\.y", "", rownames(pca_data))
pca_data$Sample <- rownames(pca_data)
pca_data$Sample <- sub("\\.x|\\.y", "", pca_data$Sample)
pca_data <- merge(pca_data, group_info, by = "Sample")
unified.pca.plot <- autoplot(pca_result, 
                             data = pca_data, 
                             colour = 'Group') + 
  ggtitle("PCA Plot of Splicing Analysis (based on PSI scores)")+ 
  geom_point(size = 6) + 
  aes(color = Group)
unified.pca.plot
ggsave(unified.pca.plot, filename = "PCA_plot.png",path =PSI_Sigma_dir ,width = 10, height =6, dpi=300 )
