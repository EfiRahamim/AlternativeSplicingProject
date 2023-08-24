library(tidyverse)
data <- read.csv("/private10/Projects/Efi/AML/rMATS/netMHCcombined.csv")
data_long <- data %>%
  pivot_longer(cols = starts_with("HLA"), names_to = "HLA", values_to = "StrongBinders") #%>%
  #pivot_longer(cols = c("SplicingType", "Group"), names_to = "Category", values_to = "Value")
treatments_desired_order <- c("SF mutations", "Mock_6h", "Indisulam", "Pladienolide-B", "Mock_18h", "5-Azacytidine", "FB23-2")
data_long$Group <- factor(data_long$Group, levels = desired_order)
# plot distribution of AS type in each treatment for each HLA allel
plot <- ggplot(data_long, aes(x = Group, y = StrongBinders, fill = SplicingType)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ HLA, ncol = 3) +
  labs(x = "Group", y = "Strong Binders") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plot)


# plot distribution of rank in each treatment for each HLA allel
plot <- ggplot(data_long, aes(x = HLA, y = StrongBinders, fill = as.factor(Rank))) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Group, ncol = 2) +
  labs(x = "HLA Allele", y = "Strong Binders") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c("dodgerblue4", "dodgerblue"))
print(plot)
