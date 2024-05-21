library(VennDiagram)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)

data_long <- read.csv("/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/18-Hours-Treatments/NovelStrongBindingEpitopes_noDups.csv")
out_dir <- '/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/18-Hours-Treatments/'


#treatments_desired_order <- c("No Treatment", "Mock (6h)", "Indisulam", "Pladienolide-B", "Mock (18h)", "5-Azacytidine", "FB23-2")
#treatments_desired_order <- c("NoTreatmentNoSF","NoTreatmentSF","Mock6","Indisulam","PladB","Mock18","5Aza","FB23-2")
#treatments_desired_order <- c("DMSO", "H3B8800")
#treatments_desired_order <- c('DMSO','Indisulam', 'PladienolideB')
#treatments_desired_order <- c('DMSO','dCEMM1','Indisulam', 'PladienolideB')
#treatments_desired_order <- c("NoTreatmentSF","Mock6","Indisulam","PladB")
treatments_desired_order <- c("NoTreatmentSF","Mock18","5Aza","FB23-2")
controls <- c("NoTreatmentSF","Mock18")
data_long$Group <- factor(data_long$Group, levels = treatments_desired_order)

# add Annotated/Non-Annotated column
data_long$Annotated <- ifelse(grepl('Ex.', data_long$Reference.Transcript), 'Non-Annotated', 'Annotated')
data_long$Annotated <- factor(data_long$Annotated, levels = c('Non-Annotated', 'Annotated'))
# modify transcripts
data_long$FixedTranscript <- gsub('Ex.|TSS.','',data_long$Reference.Transcript)

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
       path = out_dir, 
       filename = "NovelNEinTreatments.png",
       bg=NULL, width = 10, height = 6, dpi = 300)

# 2. Plot count of each HLA type in each group
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
       path = out_dir, 
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
       path = out_dir, 
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
       path = out_dir, 
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
pdf(file=paste0(out_dir,"/UpSetPlot_AllTreatments.pdf"), onefile=FALSE) # or other device
upset_plot
dev.off()

# check for peptides that appear in many groups and many HLA allels
# Count the number of unique groups and HLA's for each peptide 
candidate_peptides <- data_long %>%
  rowwise()%>%
  mutate(Gene=strsplit(ID,"_")[[1]][1])%>%
  group_by(Peptide, Gene, Splicing.Event, Annotated) %>%
  summarize(Group_Count = n_distinct(Group),
            Groups = paste(unique(Group), collapse = ", "),
            HLA_Count = n_distinct(HLA),
            HLAs = paste(unique(HLA), collapse = ", "),
            #Transcript_Counts = n_distinct(FixedTranscript),
            #Transcripts = paste(unique(FixedTranscript), collapse = ", "),
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

# count novel NE in case groups by Annotated dividing
annotated_plot <- data_long %>%
  filter(!Group %in% controls) %>%
  group_by(Peptide, Group, Annotated) %>%
  summarise(count = n()) %>%
  #mutate(g_annotated=paste(Group,'-',Annotated))%>%
  ggplot(aes(x=Group, fill = Annotated))+
  geom_bar()+
  scale_fill_brewer(palette = "Paired")+
  labs(title='Novel Neo-Epitopes Counts',
       subtitle = 'Thresholds: %Rank < 0.5, Affinity (nM) < 50',
       fill = 'Splice Junction Type')
annotated_plot
ggsave(annotated_plot, 
       path = out_dir, 
       filename = "AnnotatedNeoEpitopesCounts_plot.png",bg=NULL, width = 10, height = 6, dpi = 300)


# count novel NE in all groups
counts_plot <- ggplot(data_long, aes(x=Group, fill=Group))+
         geom_bar()+
  labs(title='Novel Neo-Epitopes Counts',
       subtitle = 'Thresholds: %Rank < 0.5, Affinity (nM) < 50')
counts_plot
ggsave(counts_plot, 
       path = out_dir, 
       filename = "NovelNeoEpitopesCounts_plot.png",bg=NULL, width = 10, height = 6, dpi = 300)

### create plots of Rank and nM comparisons
summary_data <- data_long %>%
  group_by(Group) %>%
  summarize(Rank_mean=mean(Rank),
            Rank_std = sd(Rank),
            Rank_se=sd(Rank) / sqrt(length(Rank)),
            nM_mean=mean(nM),
            nM_std = sd(nM),
            nM_se=sd(nM) / sqrt(length(nM)))
# calc p-values for Rank
Rank_p_test <- data_long %>% 
  compare_means(Rank~Group, data=.)
Rank_p_test$label <- paste(Rank_p_test$group1,"-", Rank_p_test$group2,":",Rank_p_test$p.signif,"(",Rank_p_test$method,")")
sign_label <- paste(Rank_p_test$label, collapse = '\n')
# compare Rank between groups
rank_plot <- ggplot(summary_data) +
  geom_bar( aes(x=Group, y=Rank_mean), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=Group, ymin=Rank_mean-Rank_std, ymax=Rank_mean+Rank_std), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  geom_text(aes(x = max(as.numeric(Group)), y = max(Rank_mean+Rank_std)+0.1), label = sign_label, hjust = 1, vjust = 1) +
  #geom_text(aes(x = mean(as.numeric(Group)), y = max(Rank_mean+Rank_std)), label = sign_label, hjust = 1, vjust = 1) +
  labs(title = 'Error Bars of Rank')
  
# calc p-values for nM
nM_p_test <- data_long %>% 
  compare_means(nM~Group, data=.)
nM_p_test$label <- paste(nM_p_test$group1,"-", nM_p_test$group2,":",nM_p_test$p.signif,"(",nM_p_test$method,")")
sign_label <- paste(nM_p_test$label, collapse = '\n')
nM_plot <- ggplot(summary_data) +
  geom_bar( aes(x=Group, y=nM_mean), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=Group, ymin=nM_mean-nM_std, ymax=nM_mean+nM_std), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  geom_text(aes(x = max(as.numeric(Group)), y = max(nM_mean+nM_std)+10), label = sign_label, hjust = 1, vjust = 1) +
  #geom_text(aes(x = mean(as.numeric(Group)), y = max(nM_mean+nM_std)+10), label = sign_label, hjust = 1, vjust = 1) +
  labs(title = 'Error Bars of Affinity (nM)')

g <- grid.arrange(rank_plot, nM_plot, ncol=2)
ggsave(g, 
       path = out_dir, 
       filename = "RankAffinityMeansCompare_plot_STD.png",bg=NULL, width = 10, height = 6, dpi = 300)


