library(argparse)
# CLI Arguments
parser <- ArgumentParser(description='Filter PSI-Sigma initial results.\n The filter is according to: P_value, FDR and ΔPSI value.')
parser$add_argument("-i", action="store", dest="PSI_Sigma_dir", help="Super directory of PSI-Sigma result (Where sub-directories of comparisons are located)")
parser$add_argument("-o", action="store", dest="output_dir", help="Directory to store the output.")
parser$add_argument("-dir_name", action="store", dest="output_dir_name", required=FALSE, help="The name of the actual directories of the results")
parser$add_argument("-file_name", action="store", dest="output_file_name", help="The name of PSI-Sigma results file. default: 'PSI-Sigma_r10_ir3.sorted.txt' ", default="PSI-Sigma_r10_ir3.sorted.txt")
parser$add_argument("-psi", action="store", dest="delta_PSI", help="ΔPSI (in precentages) for filtering. default: 20.", type="integer", default=20)
parser$add_argument("-pval", action="store", dest="p_value", help="P-value for filtering. default: 0.05", type="double", default=0.05)
parser$add_argument("-fdr", action="store", dest="fdr", help="FDR for filtering. default: 0.05", type="double", default=0.05)
parser$add_argument("--novelSS", action="store_true", dest="novelSS", help="PSI-Sigma results include novel transcripts. Plot will be generated for results with and without novel transcripts.")
parser$add_argument("-gene_prefix", action="store", dest="gene_prefix", help="Prefix of novel genes. default: MSTRG", default="MSTRG")
parser$add_argument("-ncol", action="store", dest="ncol", help="Number of columns to plot in the volcano plots. default: 3", type="integer", default=3)
user_args <- parser$parse_args()
stopifnot(!is.null(user_args$PSI_Sigma_dir) && !is.null(user_args$output_dir))
print(user_args)

# # DEBUG Arguments
PSI_Sigma_dir <- "/private10/Projects/Efi/CRG/GBM/PSI-Sigma/NonSotredGTF/"
output_dir <- "/private10/Projects/Efi/CRG/GBM/PSI-Sigma/NonSotredGTF/"
output_dir_name <- NULL
output_file_name <- "PSI-Sigma_r10_ir3.sorted.txt"
delta_PSI = 20
p_value = 0.05
fdr = 0.05
ncol_plot <- 3
novelSS = T
gene_prefix = "MSTRG"


# Arguments assignment
PSI_Sigma_dir <- user_args$PSI_Sigma_dir
output_dir <- user_args$output_dir
output_dir_name <- user_args$output_dir_name
output_file_name <- user_args$output_file_name
delta_PSI <- user_args$delta_PSI
p_value <- user_args$p_value
fdr <- user_args$fdr
novelSS <- user_args$novelSS
gene_prefix <- user_args$gene_prefix
ncol_plot <- user_args$ncol

library(dplyr)
library(ggplot2)

# read results files into list of data frames
if (!is.null(output_dir_name) ){
  output_dirs <- dir(PSI_Sigma_dir, recursive = T, pattern = output_dir_name, full.names = T, include.dirs = TRUE)
} else {
  output_dirs <- list.dirs(PSI_Sigma_dir, recursive = F, full.names = T)
}
full_df_list <- list()
filtered_df_list <- list()
comparisons <- c()
for (dir in output_dirs){
  if (!is.null(output_dir_name) ){
    comparison <- basename(dirname(dir)) # get comparison groups
  } else {
    comparison <- basename(dir) # get comparison groups
  }
  
  comparisons <- append(comparisons, comparison) # add comparison group to vector
  result_file <- list.files(dir, pattern= output_file_name, full.names=T) # get result file
  df <- read.csv(result_file, sep='\t') # read result file
  full_df_list[[comparison]] <- df
  filtered_df_list[[comparison]] <- subset(df,abs(ΔPSI....)>= delta_PSI & T.test.p.value < p_value & FDR..BH. < fdr )
}

merged_results <- data.frame(matrix(nrow = 0, ncol = 15))
colnames(merged_results) <- c(colnames(full_df_list[[comparisons[1]]]), "Comparison")
merged_results_filtered <- data.frame(matrix(nrow = 0, ncol = 15))
colnames(merged_results_filtered) <- c(colnames(filtered_df_list[[comparisons[1]]]), "Comparison")

for (comparison in comparisons){
  # write filtered results to csv file
  filtered_results_file <- file.path(output_dir, paste0("SplicingEventsFiltered-",comparison,"-PSI",delta_PSI,"_Pvalue", p_value,"_FDR", fdr,".csv"))
  write.csv(filtered_df_list[[comparison]], file=filtered_results_file, row.names = F) # write filtered results file
  full_df_list[[comparison]]$Comparison <- comparison
  filtered_df_list[[comparison]]$Comparison <- comparison
  merged_results <- rbind(merged_results, full_df_list[[comparison]])
  merged_results_filtered <- rbind(merged_results_filtered, filtered_df_list[[comparison]])
}

# create volcano plots
volcano_plot <- ggplot(merged_results, aes(x = ΔPSI...., y = -log10(T.test.p.value))) +
  geom_point(aes(color = ifelse(ΔPSI.... >= delta_PSI & T.test.p.value < p_value, "Inclusion",
                                 ifelse(ΔPSI.... <= -delta_PSI & T.test.p.value < p_value, "Exclusion", "Not Significant")))) +
  facet_wrap(~Comparison, ncol = ncol_plot) + 
  labs(x = "ΔPSI", y = "-log10(P-Value)", subtitle = paste0("Thresholds: |ΔPSI| > ",delta_PSI ,"%, ", "P-Value < ",p_value)) +
  #ggtitle(paste("ΔPSI:", comparison)) +
  ggtitle("Volcano plots of differential splicing events.") +
  scale_color_manual(values = c("blue", "red", "gray"))+
  guides(color = guide_legend(title = NULL))+
  geom_vline(xintercept = c(-delta_PSI, delta_PSI)) +
  geom_hline(yintercept = -log10(p_value))

volcano_plot
ggsave(filename="VolcanoPlots.png", plot = volcano_plot, path=output_dir, width = 10, height = 6, dpi = 300)

# create bar plots of the significant splicing events
merged_results_filtered$Comparison <- factor(merged_results_filtered$Comparison, levels = comparisons) # reorder the comparisons lables
filtered_splicing_event_barplot <- ggplot(merged_results_filtered, aes(y=Comparison, fill=Event.Type))+
  geom_bar(position = "stack")+
  geom_text(stat = "count", aes(label = after_stat(count)), position=position_stack(0.5)) +
  #theme_minimal()+
  labs(title = "PSI-Sigma results: Differential Splicing Events",
       subtitle = paste0("Thresholds: |ΔPSI| > ",delta_PSI ,"%, ", "P-Value < ",p_value), y = "Count")
filtered_splicing_event_barplot
ggsave(filename="SignificantSplicingEvents.png", plot = filtered_splicing_event_barplot, path=output_dir, width = 10, height = 6, dpi = 300)

# create plots without novel transcripts
if (novelSS){
  gene_prefix_pattern <- paste0("^", gene_prefix)
  volcano_plot <- ggplot(subset(merged_results, !grepl(gene_prefix_pattern, Gene.Symbol)), 
                         aes(x = ΔPSI...., y = -log10(T.test.p.value))) +
    geom_point(aes(color = ifelse(ΔPSI.... >= delta_PSI & T.test.p.value < p_value, "Inclusion",
                                   ifelse(ΔPSI.... <= -delta_PSI & T.test.p.value < p_value, "Exclusion", "Not Significant")))) +
    facet_wrap(~Comparison,ncol = ncol_plot) + 
    labs(x = "ΔPSI", y = "-log10(P-Value)", subtitle = paste0("Thresholds: |ΔPSI| > ",delta_PSI ,"%, ", "P-Value < ",p_value)) +
    #ggtitle(paste("ΔPSI:", comparison)) +
    ggtitle("Volcano plots of differential splicing events (without novel transcripts)") +
    scale_color_manual(values = c("blue", "red", "gray"))+
    guides(color = guide_legend(title = NULL))+
    geom_vline(xintercept = c(-delta_PSI, delta_PSI)) +
    geom_hline(yintercept = -log10(p_value))
  
  volcano_plot
  ggsave(filename="VolcanoPlots_noNovelSS.png", plot = volcano_plot, path=output_dir, width = 10, height = 6, dpi = 300)
  
  filtered_splicing_event_barplot <- ggplot(subset(merged_results_filtered, !grepl(gene_prefix_pattern, Gene.Symbol)), 
                                                   aes(y=Comparison, fill=Event.Type))+
    geom_bar(position = "stack")+
    geom_text(stat = "count", aes(label = after_stat(count)), position=position_stack(0.5)) +
    #theme_minimal()+
    labs(title = "PSI-Sigma results: Differential Splicing Events (without novel transcripts)",
         subtitle = paste0("Thresholds: |ΔPSI| > ",delta_PSI ,"%, ", "P-Value < ",p_value), y = "Count")
  filtered_splicing_event_barplot
  ggsave(filename="SignificantSplicingEvents_noNovelSS.png", plot = filtered_splicing_event_barplot, path=output_dir, width = 10, height = 6, dpi = 300)
}
