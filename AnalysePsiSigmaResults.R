library(argparse)
# CLI Arguments
parser <- ArgumentParser(description='Filter PSI-Sigma initial results.\n The filter is according to: P_value, FDR and ΔPSI value. Also calculate TPM of gene expression for each transcript.')
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
parser$add_argument("--filter_by_TM", action="store_true", dest="filter_tm", help="Filter results by Trans-Membrane (TM) proteins")
parser$add_argument("-tm", action="store", dest="tm_table", help="Path of TransMembrane Domains table (UniProt). Defualt: /private10/Projects/Efi/General/transmembrane_Nov23.csv", default="/private10/Projects/Efi/General/transmembrane_Nov23.csv")
parser$add_argument("-salmon", action="store", dest="salmon_dir", help="Directory of salmon results.")
parser$add_argument("-group_info", action="store", dest="group_info_file", help="Tab-delimetered file of 'Sample' and 'Group'. Header required.")

user_args <- parser$parse_args()
stopifnot(!is.null(user_args$PSI_Sigma_dir) && !is.null(user_args$output_dir))
print(user_args)

# # DEBUG Arguments
# PSI_Sigma_dir <- "/private10/Projects/Efi/AML/PSI-Sigma/All/NoTreatments/"
# output_dir <- "/private10/Projects/Efi/AML/PSI-Sigma/All/NoTreatments/TM_Results/"
# output_dir_name <- NULL
# output_file_name <- "PSI-Sigma_r10_ir3.sorted.txt"
# delta_PSI = 20
# p_value = 0.05
# fdr = 0.05
# ncol_plot <- 3
# novelSS = F
# gene_prefix = "MSTRG"
# filter_tm <- T
# tm_table <- "/private10/Projects/Efi/General/transmembrane_Nov23.csv"
# salmon_dir <- "/private10/Projects/Efi/AML/Salmon_gencode_v28/"
# salmon_suffix = ".quant.genes.sf"
# group_info_file <- "/private10/Projects/Efi/AML/Salmon_gencode_v28/NoTreatments-Info.txt"


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
filter_tm <- user_args$filter_tm
tm_table <- user_args$tm_table
salmon_dir <- user_args$salmon_dir
salmon_suffix <- ".quant.genes.sf"
group_info_file <- user_args$group_info_file

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(VennDiagram)


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
  result_file <- list.files(dir, pattern= output_file_name, full.names=T) # get result file
  if (length(result_file) == 0){
    next
  }
  if (!is.null(output_dir_name) ){
    comparison <- basename(dirname(dir)) # get comparison groups
  } else {
    comparison <- basename(dir) # get comparison groups
  }
  comparisons <- append(comparisons, comparison) # add comparison group to vector
  df <- read.csv(result_file, sep='\t') # read result file
  # filter by TM proteins
  if (filter_tm){
    # load TM table
    tm <- read.csv(tm_table)
    # merge AS results with TM table according to Gene Synbol
    df <- df[df$Gene.Symbol %in% tm$GeneSymbol,]
    print(paste0("Merging with TM table: ", nrow(df), " results."))
  }
  full_df_list[[comparison]] <- df
  filtered_df_list[[comparison]] <- subset(df,abs(ΔPSI....)>= delta_PSI & T.test.p.value < p_value & FDR..BH. < fdr )
}

merged_results <- data.frame(matrix(nrow = 0, ncol = 15))
colnames(merged_results) <- c(colnames(full_df_list[[comparisons[1]]]), "Comparison")
merged_results_filtered <- data.frame(matrix(nrow = 0, ncol = 15))
colnames(merged_results_filtered) <- c(colnames(filtered_df_list[[comparisons[1]]]), "Comparison")

# read Gene Expression data
salmon_files <- list.files(path=salmon_dir, pattern = salmon_suffix, full.names = T, recursive = T )
col_names <- c("Name")
merged_tpm <- data.frame(matrix(nrow = 0, ncol = length(col_names)))
colnames(merged_tpm) <- col_names
# read GroupInfo data 
group_info <- read.csv(group_info_file, sep='\t', header=T)
groups <- unique(group_info$Group)
for (file in salmon_files){
  sample <- gsub(salmon_suffix, "",basename(file))
  if (!(sample %in% group_info$Sample)){
    next
  }
  df <- read.csv(file, sep ='\t', header = T)[,c(1,4)]
  colnames(df)[2] <- sample
  if (nrow(merged_tpm) == 0){
    merged_tpm <- df
  } else {
    merged_tpm <- merge(merged_tpm, df, by = "Name")
  }
}
for (group in groups){
  samples <- subset(group_info, Group==group)$Sample
  if (length(samples) == 1) { # handle case of 1 sample in a group
    merged_tpm[[paste0("TPM_mean_",group)]] <- merged_tpm[,samples]
    merged_tpm[[paste0("TPM_std_",group)]] <- 0
  } else {
    merged_tpm[[paste0("TPM_mean_",group)]] <- rowMeans(merged_tpm[,samples])
    merged_tpm[[paste0("TPM_std_",group)]] <- apply(merged_tpm[,samples], 1, sd)
  }
}
merged_tpm$FixedTranscript <- sub("\\..*", "", merged_tpm$Name)
target_exons_list <- list() # for Vann diagram
for (comparison in comparisons){
  #add TPM value for each transcript in each group
  TPM_mean_groupA <- paste0("TPM_mean_", strsplit(comparison, "_vs_")[[1]][1])
  TPM_mean_groupB <- paste0("TPM_mean_", strsplit(comparison, "_vs_")[[1]][2])
  TPM_std_groupA <- paste0("TPM_std_", strsplit(comparison, "_vs_")[[1]][1])
  TPM_std_groupB <- paste0("TPM_std_", strsplit(comparison, "_vs_")[[1]][2])
  
  filtered_df_list[[comparison]]$FixedTranscript <- gsub("Ex\\.|TSS\\.|Ex\\.TSS\\.", "", filtered_df_list[[comparison]]$Reference.Transcript)
  filtered_df_list[[comparison]]$FixedTranscript <- sub("\\..*", "", filtered_df_list[[comparison]]$FixedTranscript)
  # filtered_df_list[[comparison]] <- merge(filtered_df_list[[comparison]],
  #                                         merged_tpm[,c("FixedTranscript",TPM_mean_groupA, TPM_mean_groupB,TPM_std_groupA,TPM_std_groupB)],
  #                                         by.x='FixedTranscript', by.y='FixedTranscript',
  #                                         all.x = T)
  filtered_df_list[[comparison]] <- merge(filtered_df_list[[comparison]],
                                          merged_tpm[,c("Name",TPM_mean_groupA, TPM_mean_groupB,TPM_std_groupA,TPM_std_groupB)],
                                          by.x='Gene.Symbol', by.y='Name',
                                          all.x = T)
  # write filtered results to csv file
  filtered_results_file <- file.path(output_dir, paste0("SplicingEventsFiltered-",comparison,"-PSI",delta_PSI,"_Pvalue", p_value,"_FDR", fdr,".csv"))
  write.csv(filtered_df_list[[comparison]], file=filtered_results_file, row.names = F) # write filtered results file
  full_df_list[[comparison]]$Comparison <- comparison
  filtered_df_list[[comparison]]$Comparison <- comparison
  merged_results <- rbind(merged_results, full_df_list[[comparison]])
  merged_results_filtered <- rbind(merged_results_filtered, filtered_df_list[[comparison]]%>%select(-TPM_mean_groupA, -TPM_mean_groupB, -TPM_std_groupA, -TPM_std_groupB))
  target_exons_list[[comparison]] <- paste0(filtered_df_list[[comparison]]$Target.Exon,
                                            filtered_df_list[[comparison]]$FixedTranscript,
                                            filtered_df_list[[comparison]]$Event.Type)
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
if (length(filtered_df_list) == 1 ){
  ordered_event_types <- names(sort(table(merged_results_filtered$Event.Type), decreasing = TRUE))
  merged_results_filtered$Event.Type <- factor(merged_results_filtered$Event.Type, levels = ordered_event_types)
  filtered_splicing_event_barplot <- ggplot(merged_results_filtered, aes(x=Event.Type,fill=Event.Type))+
    geom_bar(position = "stack")+
    geom_text(stat = "count", aes(label = after_stat(count)), position=position_stack(0.5)) +
    #theme_minimal()+
    labs(title = "PSI-Sigma results: Differential Splicing Events",
         subtitle = "|ΔPSI| > 20%, P-value & FDR < 0.05", y = "Count")
  filtered_splicing_event_barplot
  ggsave(filename="SignificantSplicingEvents.png", plot = filtered_splicing_event_barplot, path=output_dir, width = 10, height = 6, dpi = 300)
  
} else {
  merged_results_filtered$Comparison <- factor(merged_results_filtered$Comparison, levels = comparisons) # reorder the comparisons lables
  filtered_splicing_event_barplot <- ggplot(merged_results_filtered, aes(y=Comparison, fill=Event.Type))+
    geom_bar(position = "stack")+
    geom_text(stat = "count", aes(label = after_stat(count)), position=position_stack(0.5)) +
    #theme_minimal()+
    labs(title = "PSI-Sigma results: Differential Splicing Events",
         subtitle = paste0("Thresholds: |ΔPSI| > ",delta_PSI ,"%, ", "P-Value < ",p_value), y = "Count")
  filtered_splicing_event_barplot
  ggsave(filename="SignificantSplicingEvents.png", plot = filtered_splicing_event_barplot, path=output_dir, width = 10, height = 6, dpi = 300)
}
if (length(filtered_df_list) >= 2 ){
  # create Vann diagram of splicing events interesections
  num_groups <- length(target_exons_list)
  circle_colors <- rainbow(num_groups) # Generate colors dynamically based on the number of groups
  palette <- brewer.pal(num_groups, "Dark2")
  plot_path = file.path(output_dir, "VennDiagram.png")
  venn.diagram(target_exons_list, category.names = comparisons, 
               filename=plot_path,
               disable.logging=T,
               force.unique = T,
               height=2000, width=3000, resolution=250,
               col = palette,
               fill = c(alpha(palette, 0.3)),
               main = "Division of significant splicing events among different comparisons",
               main.cex=1.6,#,
               sub = paste0("Thresholds: |ΔPSI| > ",delta_PSI ,"%, ", "P-Value < ",p_value)
  )
}
# creat BED file of Intron Retention events
# groupA_ri_exons <- merged_results_filtered[(merged_results_filtered$Event.Type=='IR' | merged_results_filtered$Event.Type=='IR (overlapping region)') & merged_results_filtered$ΔPSI.... > 0, "Target.Exon"]
# BED_ri_exons <- strsplit(groupA_ri_exons, "[:-]")
# BED_ri_exons_matrix <- do.call(rbind, BED_ri_exons)
# BED_ri_exons_df <- as.data.frame(BED_ri_exons_matrix)
# BED_ri_exons_df <- unique(BED_ri_exons_df)
# colnames(BED_ri_exons_df) <- c("chr", "start", "end")
# BED_ri_exons_df$start <- as.numeric(as.character(BED_ri_exons_df$start))
# BED_ri_exons_df$end <- as.numeric(as.character(BED_ri_exons_df$end))
# sorted_BED_ri_exons_df <- BED_ri_exons_df[order(BED_ri_exons_df$chr, BED_ri_exons_df$start, BED_ri_exons_df$end), ]
# write.table(sorted_BED_ri_exons_df, file = file.path(output_dir,"IR_events.sorted.bed"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE) # Write the sorted matrix into a BED file

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
