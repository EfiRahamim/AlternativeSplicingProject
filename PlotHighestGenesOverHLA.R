# Load necessary libraries
library(ggplot2)
library(dplyr)

# Read the CSV file into a data frame
data <- read.csv("/private10/Projects/Efi/AML_secondbatch/test/SE_TestForNetMHC/Treatment_StrongBinders_All.csv")  # Replace with your actual CSV file name

# Initialize an empty list to store plots
plots <- list()

# Iterate through each HLA allele column
hla_columns <- colnames(data)[grepl("^HLA", colnames(data))]
for (hla_col in hla_columns) {
  # Sort data by strong binders value in descending order and select top 25 genes
  top_genes <- data %>%
    arrange(desc(!!sym(hla_col))) %>%
    group_by(GeneSymbol) %>%
    slice(1) %>%
    ungroup() %>%
    arrange(desc(!!sym(hla_col))) %>%
    top_n(25, wt = !!sym(hla_col))
  
  # Create a bar plot for the top genes in the current HLA allele
  plot <- ggplot(top_genes, aes(x = reorder(GeneSymbol, -!!sym(hla_col)), y = !!sym(hla_col))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = hla_col, x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Append the plot to the list of plots
  plots[[hla_col]] <- plot
}

# Arrange and display the plots using gridExtra package
library(gridExtra)
grid.arrange(grobs = plots, ncol = 4, top = "Top genes with highest strong binders value - Treatment group", bottom ="GeneSymbol", left="Strong Binders")  # Change ncol as per your preference
