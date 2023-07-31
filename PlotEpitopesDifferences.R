# Step 1: Load the required libraries
library(ggplot2)
library(tidyr)
library(dplyr)

SB_result_path <- "/private5/Projects/Efi/AS/AML/X204SC22102645-Z01-F005/rMATS/Mock_Indisulam_6h/A3SS/StrongBinders_All.csv"
output_dir <- "/private5/Projects/Efi/AS/AML/X204SC22102645-Z01-F005/rMATS/Mock_Indisulam_6h/A3SS/"
group1 <- "Mock"
group2 <- "Indisulam"
as <- "A3SS"

if (as=='SE') {
  AS_type <- "Exon Skipping"
  color<-"steelblue"
} else if (as == 'A5SS') {
  AS_type <- "Alternative 5' Splice Site"
  color<-"goldenrod"
} else if (as == 'A3SS') {
  AS_type <- "Alternative 3' Splice Site"
  color<-"orchid"
} else if (as == 'RI') {
  AS_type <- "Intron retention"
  color<-"plum"
}

# Step 2: Read the CSV file into a data frame
gene_data <- read.csv(SB_result_path)

# Step 3: Calculate the sum of positive values in each row and add it as a new column
gene_data <- gene_data %>%
  rowwise() %>%
  mutate(Sum_Positive_Values = sum(c_across(starts_with("HLA."))[c_across(starts_with("HLA.")) > 0]))

# Sort the data frame by "Sum_Positive_Values" in descending order and keep only the top 25 genes
top_25_genes <- gene_data %>%
  arrange(desc(Sum_Positive_Values)) %>%
  distinct(GeneSymbol, .keep_all = TRUE) %>%
  head(2)

# Reorder genes based on "Sum_Positive_Values"
top_25_genes$GeneSymbol <- reorder(top_25_genes$GeneSymbol, -top_25_genes$Sum_Positive_Values)

# Create a bar plot for the top 25 genes in descending order of "Sum_Positive_Values"
plot <- ggplot(top_25_genes, aes(x = GeneSymbol, y = Sum_Positive_Values)) +
  geom_bar(stat = "identity", fill = color) +
  labs(x = "Genes", y = "Sum of new epitopes", title=paste0("Top Genes with Highest Sum of New Epitopes*\n",group1, " vs. " ,group2, " - ", AS_type ),
       caption = "*HLA Types that were checked: HLA-A0101,HLA-A0201,HLA-A0301,HLA-A2402,HLA-A2601,HLA-B0702,HLA-B0801,HLA-B1501,HLA-B2705,HLA-B3901,HLA-B4001,HLA-B5801")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

plot
ggsave(plot, path = output_dir, filename = "EpitopesDifferences.png",bg=NULL, width = 10, height = 6, dpi = 300)



















# # Step 3: Reshape the data from wide to long format
# gene_data_long <- gather(gene_data, HLA_allele, Value, starts_with("HLA."))
# gene_data_long_filtered <- gene_data_long[gene_data_long$Value > 0, ]
# gene_data_long_filtered$GeneSymbol <- factor(gene_data_long_filtered$GeneSymbol, levels = unique(gene_data_long_filtered$GeneSymbol[order(-gene_data_long_filtered$Value)]))
# 
# # Step 6: Group the data by HLA allele and gene, and arrange the data in descending order by Value within each group
# gene_data_long_sorted <- gene_data_long_filtered %>%
#   group_by(HLA_allele, GeneSymbol) %>%
#   arrange(desc(Value)) %>%
#   ungroup()
# 
# # Step 7: Select the top 25 genes for each HLA allele
# top_25_data <- gene_data_long_sorted %>%
#   group_by(HLA_allele) %>%
#   top_n(10)
# 
# # Step 4: Create the plot
# plot <- ggplot(top_25_data, aes(x = GeneSymbol, y = Value, fill = HLA_allele)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(x = "Gene", y = "HLA Allele Value", fill = "HLA Allele",
#        title = paste0("Epitopes difference between ", group1, " and ",group2, " in ", AS_type)) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         plot.title = element_text(hjust = 0.5, face = "bold"),
#         plot.background = element_rect(fill = "white"), complete = TRUE) +
#   facet_wrap(~ HLA_allele, ncol = 2,scales = "free_x")
# 
# # save the plot
# #ggsave(plot, path = output_dir, filename = "EpitopesDifferences.png",bg=NULL, width = 10, height = 6, dpi = 300)
# plot