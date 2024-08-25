library(VennDiagram)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)

# Output directory for plotting
out_dir <- "/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/PerPatient/IDMO1343/"

# Define the control and treatments groups
treatments_desired_order <- c("Mock6","Indisulam", "Madrasin", "H3B-8800")
controls <- c("Mock6")
treatments_desired_order <- c("Mock6","Indisulam")

treatments_desired_order <- c("Mock18","5Aza")
controls <- c("Mock18")

treatments_desired_order <- c("NoTreatmentSF","Mock18","FB23-2", "5Aza")
controls <- c("NoTreatmentSF","Mock18")

# find the neoepitopes file and read it
neo_epitopes_file = list.files(out_dir, pattern = 'NovelStrongBindingEpitopes_noDups.csv', full.names = TRUE)
data_long <- read.csv(neo_epitopes_file)
# treat groups names as factor
data_long$Group <- factor(data_long$Group, levels = treatments_desired_order)
# filter out peptides of control groups
data_long <- filter(data_long, !(Group %in% controls))


# add Annotated/Non-Annotated column
data_long$Annotated <- ifelse(grepl('Ex.', data_long$Reference.Transcript), 'Non-Annotated', 'Annotated')
data_long$Annotated <- factor(data_long$Annotated, levels = c('Non-Annotated', 'Annotated'))
# modify transcripts
data_long$FixedTranscript <- gsub('Ex.|TSS.','',data_long$Reference.Transcript)

# 1. Plot count of each treatment in each HLA type, divided by Splicing Event
HLA_treatment_AStype <- data_long %>%
  #group_by(Group, HLA, Splicing.Event) %>%
  mutate(Splicing.Event = ifelse(Splicing.Event == 'IR_OLR', 'IR',
                                 ifelse(Splicing.Event == 'TSS_A5SS', 'A5SS',
                                        ifelse(Splicing.Event == 'TSS_A3SS', 'A3SS', Splicing.Event)))) %>%
  group_by(Peptide, Group, Splicing.Event) %>%
  summarise(count = n()) %>%
  #ggplot(aes(x = Group, y = count, fill = Splicing.Event)) +
  #geom_bar(stat = "identity", position = "dodge") +
  #facet_wrap(~HLA, ncol = 3) +
  #theme_bw()+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  #labs(title = "Novel Neo-Epitopes From Splicing Events", y = "Novel Neo-Epitopes Count")
  #theme_minimal()
  ggplot(aes(x=Group, fill=Splicing.Event))+
  geom_bar(position="fill")+
  #scale_fill_brewer(palette = "Paired")+
  labs(title='Novel Neo-Epitopes Divided by Splicing Events',
       subtitle = 'Thresholds: %Rank < 0.5, Affinity (nM) < 50',
       fill = 'Splicing Event',
       y='Propotion')+
  geom_text(aes(label=paste0(signif(..count.. / tapply(..count.., ..x.., sum)[as.character(..x..)], digits=3)*100,'%')),
    stat="count",
    position=position_fill(vjust=0.5),
    size=3)

print(HLA_treatment_AStype)
ggsave(HLA_treatment_AStype, 
       path = out_dir, 
       filename = "NovelNEBySplicingEvents.png",
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
  facet_wrap(~Group, ncol = 3) +
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
  facet_wrap(~ Group, ncol = 3) +
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
  facet_wrap(~ Group, ncol = 3) +
  labs(x = "HLA Allele", y = "Strong Binders", title = "Affinity (nM) Distribution") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c("brown4", "chocolate","burlywood"))
print(nM_plot)
ggsave(nM_plot, 
       path = out_dir, 
       filename = "Affinity_nM_plot.png",bg=NULL, width = 10, height = 6, dpi = 300)

# plot UpSet plot
if (length(setdiff(treatments_desired_order, controls)) > 1){
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
  upset_plot <- upset(df_reshaped[,-1], 
                      sets=names(df_reshaped)[-1],
                      main.bar.color = "#4e79a7",
                      order.by = "freq",
                      empty.intersections = "on",
                      set_size.show=T)
                      #set_size.angles=45,
                      #set_size.numbers_size=8)
  pdf(file=paste0(out_dir,"/UpSetPlot_AllTreatments.pdf"), onefile=FALSE,width=12, height=10) # or other device
  upset_plot
  dev.off()
}

# check for peptides that appear in many groups and many HLA allels
# Count the number of unique groups and HLA's for each peptide 
candidate_peptides <- data_long %>%
  rowwise()%>%
  mutate(Gene=strsplit(ID,"_")[[1]][1])%>%
  
  
  #### CHECK!!#####
  #group_by(Peptide, Gene, Splicing.Event,Exon.Type,Annotated)%>%
  group_by(Peptide, Gene, Splicing.Event, Target.Exon, Exon.Type,Annotated)%>%
  #### CHECK!!#####
  
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
  geom_bar(position="fill")+
  scale_fill_brewer(palette = "Paired")+
  labs(title='Novel Neo-Epitopes Counts',
       subtitle = 'Thresholds: %Rank < 0.5, Affinity (nM) < 50',
       fill = 'Splice Junction Type',
       y='Propotion')+
  geom_text(
    aes(label=paste0(signif(..count.. / tapply(..count.., ..x.., sum)[as.character(..x..)], digits=3)*100,'%')),
    stat="count",
    position=position_fill(vjust=0.5),
    size=4)
  #geom_text(stat = "count", aes(label = after_stat(count)), position=position_stack(0.5))
annotated_plot
ggsave(annotated_plot, 
       path = out_dir, 
       filename = "AnnotatedNeoEpitopesCounts_plot.png",bg=NULL, width = 10, height = 6, dpi = 300)
