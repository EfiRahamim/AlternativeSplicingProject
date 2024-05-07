library(VennDiagram)
library(dplyr)
data <- read.csv("/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/NovelStrongBindingEpitopes_noDups.csv")
out_dir <- "/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/"
# data_long <- data %>%
#   #pivot_longer(cols = starts_with("HLA"), names_to = "HLA_BindingScore", values_to = "Value", values_drop_na = TRUE)
#   pivot_longer(cols = ends_with("_Rank"), names_to = "HLA_Rank", values_to = "Rank", values_drop_na = T) %>%
#   pivot_longer(cols = ends_with("_nM"), names_to = "HLA_nM", values_to = "nM", values_drop_na = T) %>%
#   mutate(HLA = gsub("_Rank", "", HLA_Rank)) %>%
#   select(-HLA_Rank, -HLA_nM) %>%
#   relocate(HLA, .after = 4)
data_long <- data
#treatments_desired_order <- c("No Treatment", "Mock (6h)", "Indisulam", "Pladienolide-B", "Mock (18h)", "5-Azacytidine", "FB23-2")
treatments_desired_order <- c("NoTreatmentNoSF","NoTreatmentSF","Mock6","Indisulam","PladB","Mock18","5Aza","FB23-2")
data_long$Group <- factor(data_long$Group, levels = treatments_desired_order)

# Temporary: filter out A3SS|A5SS events since sequences may not be corrected
# data_long <- data_long %>%
#   filter(!grepl('A5SS|A3SS',Splicing.Event))

# 1. Plot count of each treatment in each HLA type, divided by Splicing Event
HLA_treatment_AStype <- data_long %>%
  group_by(Group, HLA, Splicing.Event) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = Group, y = count, fill = Splicing.Event)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~HLA, ncol = 3) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title = "Novel Neo-Epitopes From Treatments", y = "Novel Neo-Epitopes Count")
  #theme_minimal()
print(HLA_treatment_AStype)
ggsave(HLA_treatment_AStype, 
       path = "/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/", 
       filename = "NovelNEinTreatments.png",
       bg=NULL, width = 10, height = 6, dpi = 300)

# 2. Plot count of each HLA type in each group
#treatments_desired_order <- c("Mock (6h)", "Mock (18h)", "Indisulam", "5-Azacytidine", "Pladienolide-B","FB23-2", "SF Mutations")
treatments_desired_order <- c("NoTreatmentNoSF","NoTreatmentSF","Mock6","Mock18","Indisulam","5Aza","PladB","FB23-2")
data_long$Group <- factor(data_long$Group, levels = treatments_desired_order)
group_colors <- c("gold", "firebrick","dodgerblue", "deeppink" , "darkviolet","darksalmon","darkorange", "darkgreen", "darkblue", "cyan", "coral", "cadetblue", "red")
Treatment_HLA <- data_long %>%
  group_by(HLA, Group) %>%
  summarise(count = n()) %>%
  arrange(Group, desc(count)) %>%  # Arrange within each Group
  ungroup() %>%
  group_by(Group) %>%
  mutate(HLA = factor(HLA, levels = unique(HLA))) %>% 
  ggplot(aes(x = HLA, y = count, fill = HLA)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Group, ncol = 2) +
  labs(title = "Novel Neo-Epitopes Distribution over HLA alleles", y = "Novel Neo-Epitopes Count")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = group_colors)
print(Treatment_HLA)
ggsave(Treatment_HLA, 
       path = "/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/", 
       filename = "Treatment_grid_descending.png",
       bg=NULL, width = 10, height = 6, dpi = 300)

# plot rank distribution
rank_plot <- data_long %>%
  mutate(rank_category = case_when(
    Rank < 0.1 ~ "<0.1",
    Rank <= 0.5 ~ "0.1-0.5",
    TRUE ~ ">0.5"
  )) %>%
  group_by(Group, HLA, rank_category) %>%
  summarise(count = n()) %>%
  arrange(Group, desc(count)) %>%  # Arrange within each Group
  ungroup() %>%
  group_by(Group) %>%
  mutate(HLA = factor(HLA, levels = unique(HLA))) %>%
  ggplot(aes(x = HLA, y = count, fill = rank_category)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Group, ncol = 2) +
  labs(x = "HLA Allele", y = "Strong Binders", title = "%Rank Distribution") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c("dodgerblue4", "dodgerblue"))
print(rank_plot)
ggsave(rank_plot, 
       path = "/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/", 
       filename = "Rank_plot.png",bg=NULL, width = 10, height = 6, dpi = 300)

# plot nM distribution
nM_plot <- data_long %>%
  mutate(nM_category = case_when(
    nM < 10 ~ "<10nM",
    nM <= 30 ~ "10nM-30nM",
    TRUE ~ "30nM-50nM"
  )) %>%
  group_by(Group, HLA, nM_category) %>%
  summarise(count = n()) %>%
  arrange(Group, desc(count)) %>%  # Arrange within each Group
  ungroup() %>%
  group_by(Group) %>%
  mutate(HLA = factor(HLA, levels = unique(HLA))) %>%
  ggplot(aes(x = HLA, y = count, fill = nM_category)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Group, ncol = 2) +
  labs(x = "HLA Allele", y = "Strong Binders", title = "Affinity (nM) Distribution") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c("brown4", "chocolate","burlywood"))
print(nM_plot)
ggsave(nM_plot, 
       path = "/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/", 
       filename = "Affinity_nM_plot.png",bg=NULL, width = 10, height = 6, dpi = 300)

# plot UpSet plot
library(UpSetR)
# Count the number of occurrences of each peptide in each group
df_counts <- data_long %>%
  group_by(Group, Peptide) %>%
  summarise(count = as.integer(n() > 0)) %>%
  ungroup()

# Reshape the data frame to have groups as columns and counts as values
df_reshaped <- df_counts %>%
  pivot_wider(names_from = Group, values_from = count, values_fill = 0)

df_reshaped <- as.data.frame(df_reshaped)

str(df_reshaped)
upset_plot <- upset(df_reshaped[,-1], 
                    sets=names(df_reshaped)[-1],
                    main.bar.color = "#4e79a7",
                    order.by = "freq",
                    empty.intersections = "on")
# save the upset plot
pdf(file="/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/UpSetPlot_AllTreatments.pdf", onefile=FALSE) # or other device
upset_plot
dev.off()

# check for peptides that appear in many groups and many HLA allels
# Count the number of unique groups and HLA's for each peptide 
candidate_peptides <- data_long %>%
  filter(!grepl('A5SS|A3SS', Splicing.Event)) %>%
  rowwise()%>%
  mutate(Gene=strsplit(ID,"_")[[1]][1])%>%
  group_by(Peptide, Gene, Splicing.Event) %>%
  summarize(Group_Count = n_distinct(Group),
            Groups = paste(unique(Group), collapse = ", "),
            HLA_Count = n_distinct(HLA),
            HLAs = paste(unique(HLA), collapse = ", "),
            Avg.Rank = mean(Rank),
            Avg.nM = mean(nM),
            Avg.PSI_PeptideSource = mean(Avg.PSI_PeptideSource),
            Avg.PSI_OtherGroups = mean(Avg.PSI_OtherGroup),
            Avg.TPM_PeptideSource = mean(Avg.TPM_PeptideSource),
            Avg.TPM_OtherGroups = mean(Avg.TPM_OtherGroup),
            SplicingIndex = paste(unique(SplicingIndex), collapse = "/"))

# write candidate peptides file
write.csv(candidate_peptides,
          file = paste0(out_dir, "PeptidesCandidates.csv"),
          row.names = F)
