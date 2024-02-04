# filter true TM events
dir_path <- "/private10/Projects/Efi/AML/PSI-Sigma/SRSF2/6-Hours-Treatments/TM_Results/SplicingEventsFiles/NoTreatment_vs_PladB/"
comparison <- "NoTreatment_vs_PladB"
results_file <- paste0(dir_path, comparison,"_analyzed.csv")
candidate_TM_events <- paste0(dir_path,"/SplicingEventsFiles/Candidate_TM_events.txt")
results_df <- read.csv(results_file)
colnames(results_df)[10] <- 'dPSI'
colnames(results_df)[11] <- 'T.test.p.value'
TM_events_df <- unique(read.csv(candidate_TM_events, header = F))
results_filtered <- merge(results_df, TM_events_df, by.x = 'X', by.y = 'V1')
write.csv(results_filtered, row.names = F, 
          file = paste0(dir_path,comparison,"_TM_candidates.csv"))


# check for mutual events in 6-hours-treatments 

# read files of filtered TM events
NTNSF_NT_df <- read.csv("/private10/Projects/Efi/AML/PSI-Sigma/All/NoTreatments/TM_Results/NoTreatmentNoSF_vs_NoTreatmentSF/NoTreatmentNoSF_vs_NoTreatmentSF_TM_candidates.csv")
NT_pladB_df <- read.csv("/private10/Projects/Efi/AML/PSI-Sigma/SRSF2/6-Hours-Treatments/TM_Results/SplicingEventsFiles/NoTreatment_vs_PladB/NoTreatment_vs_PladB_TM_candidates.csv")
NT_indisulam_df <- read.csv("/private10/Projects/Efi/AML/PSI-Sigma/SRSF2/6-Hours-Treatments/TM_Results/SplicingEventsFiles/NoTreatment_vs_Indisulam/NoTreatment_vs_Indisulam_TM_candidates.csv")
Mock6_pladB_df <- read.csv("/private10/Projects/Efi/AML/PSI-Sigma/SRSF2/6-Hours-Treatments/TM_Results/SplicingEventsFiles/Mock6_vs_PladB/Mock6_vs_PladB_TM_candidates.csv")
Mock6_indisulam_df <- read.csv("/private10/Projects/Efi/AML/PSI-Sigma/SRSF2/6-Hours-Treatments/TM_Results/SplicingEventsFiles/Mock6_vs_Indisulam/Mock6_vs_Indisulam_TM_candidates.csv")
# mutate data frames for "Inclusion/Exclusion" column
NTNSF_NT_df <- NTNSF_NT_df %>%
  mutate(`Inc/Exc` = ifelse(dPSI > 0, "Inclusion", "Exclusion"))
NT_pladB_df <- NT_pladB_df %>%
  mutate(`Inc/Exc` = ifelse(dPSI > 0, "Inclusion", "Exclusion"))
NT_indisulam_df <- NT_indisulam_df %>%
  mutate(`Inc/Exc` = ifelse(dPSI > 0, "Inclusion", "Exclusion"))
Mock6_pladB_df <- Mock6_pladB_df %>%
  mutate(`Inc/Exc` = ifelse(dPSI > 0, "Inclusion", "Exclusion"))
Mock6_indisulam_df <- Mock6_indisulam_df %>%
  mutate(`Inc/Exc` = ifelse(dPSI > 0, "Inclusion", "Exclusion"))
# subset data frames
cols_to_keep <- c("X", "Gene.Symbol", "Event.Region","Target.Exon","Event.Type","Exon.Type","FixedTranscript","dPSI","Inc/Exc","InclusionAAseq", "ExclusionAAseq")
NTNSF_NT_df <- NTNSF_NT_df[,c(cols_to_keep,colnames(NTNSF_NT_df)[17], colnames(NTNSF_NT_df)[18] )]
NT_pladB_df <- NT_pladB_df[,c(cols_to_keep,colnames(NT_pladB_df)[17], colnames(NT_pladB_df)[18] )]
NT_indisulam_df <- NT_indisulam_df[,c(cols_to_keep,colnames(NT_indisulam_df)[17], colnames(NT_indisulam_df)[18] )]
Mock6_pladB_df <- Mock6_pladB_df[,c(cols_to_keep,colnames(Mock6_pladB_df)[17], colnames(Mock6_pladB_df)[18] )]
Mock6_indisulam_df <- Mock6_indisulam_df[,c(cols_to_keep,colnames(Mock6_indisulam_df)[17], colnames(Mock6_indisulam_df)[18] )]
# define column to merge by
col_to_merge <- c('FixedTranscript', 'Event.Region', 'Gene.Symbol', 'Target.Exon', 'Event.Type', 'Inc/Exc')
# merge
library(tidyverse)

# merge all data frames
list_df_all <- list(NT_pladB_df, NT_indisulam_df,Mock6_pladB_df,Mock6_indisulam_df)
merged_all <- list_df_all %>% reduce(inner_join, by=col_to_merge)
#merged_df <- merge(pladB_df, indisulam_df, by <- col_to_merge)
write.csv(merged_all, row.names = F, file = "/private10/Projects/Efi/AML/PSI-Sigma/All/6-Hours-Treatments/TM_Results/MutualEventsCandidates/6HtreatmentsAll_mutualEvents.csv")

# merge only NoTreatments
list_df_NT <- list(NT_pladB_df, NT_indisulam_df)
merged_NT <- list_df_NT %>% reduce(inner_join, by=col_to_merge)
write.csv(merged_NT, row.names = F, file = "/private10/Projects/Efi/AML/PSI-Sigma/SRSF2/6-Hours-Treatments/TM_Results/SplicingEventsFiles/MutualEventsCandidates/NoTreatments_mutualEvents.csv")

# merge only Mock6
list_df_mock6 <- list(Mock6_pladB_df,Mock6_indisulam_df)
merged_mock6 <- list_df_mock6 %>% reduce(inner_join, by=col_to_merge)
write.csv(merged_mock6, row.names = F, file = "/private10/Projects/Efi/AML/PSI-Sigma/All/6-Hours-Treatments/TM_Results/MutualEventsCandidates/Mock6_mutualEvents.csv")

# merge only indisulam
list_df_indisulam <- list(NT_indisulam_df,Mock6_indisulam_df)
merged_indisulam <- list_df_indisulam %>% reduce(inner_join, by=col_to_merge)
write.csv(merged_indisulam, row.names = F, file = "/private10/Projects/Efi/AML/PSI-Sigma/SRSF2/6-Hours-Treatments/TM_Results/SplicingEventsFiles/MutualEventsCandidates/Indisulam_mutualEvents.csv")

# merge only PladB
list_df_pladB <- list(Mock6_pladB_df,NT_pladB_df)
merged_pladB <- list_df_pladB %>% reduce(inner_join, by=col_to_merge)
write.csv(merged_pladB, row.names = F, file = "/private10/Projects/Efi/AML/PSI-Sigma/SRSF2/6-Hours-Treatments/TM_Results/SplicingEventsFiles/MutualEventsCandidates/PladB_mutualEvents.csv")

# count results
counts <- data.frame(Intersection = c('All', 
                                 'Indisulam+PladB (vs. NoTreatments)',
                                 'Indisulam+PladB (vs. Mock)', 
                                 'Indisulam (vs. Mock/No Treatment)',
                                 'PladB (vs. Mock/No Treatment)'),
                     Mutual_Events = c(length(merged_all$Gene.Symbol),
                                       length(merged_NT$Gene.Symbol),
                                       length(merged_mock6$Gene.Symbol),
                                       length(merged_indisulam$Gene.Symbol),
                                       length(merged_pladB$Gene.Symbol)))
library(ggplot2)
intersect_plot<-ggplot(counts, aes(x = Intersection, y = Mutual_Events, fill = Intersection)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Mutual_Events), vjust = -0.5, size = 3) +
  labs(title = "Mutual splicing events among different comparisons - only in SRSF2mut samples", y = "Mutual Events") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),# angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
intersect_plot
ggsave(intersect_plot, filename = "/private10/Projects/Efi/AML/PSI-Sigma/SRSF2/6-Hours-Treatments/TM_Results/TM_IntersectionEvents.png",width = 10, height = 6, dpi = 300 )

# check for mutual events in 18-hours-treatments (only NoTreatments vs treatments)

# filter true TM events
results_file <- "/private10/Projects/Efi/AML/PSI-Sigma/All/18-Hours-Treatments/TM_Results/NoTreatmentSF_vs_5Aza/NoTreatmentSF_vs_5Aza_analyzed.csv"
candidate_TM_events <- "/private10/Projects/Efi/AML/PSI-Sigma/All/18-Hours-Treatments/TM_Results/NoTreatmentSF_vs_5Aza/SplicingEventsFiles/Candidate_TM_events.txt"
results_df <- read.csv(results_file)
colnames(results_df)[10] <- 'dPSI'
colnames(results_df)[11] <- 'T.test.p.value'
TM_events_df <- unique(read.csv(candidate_TM_events, header = F))
results_filtered <- merge(results_df, TM_events_df, by.x = 'X', by.y = 'V1')
write.csv(results_filtered, row.names = F, 
          file = "/private10/Projects/Efi/AML/PSI-Sigma/All/18-Hours-Treatments/TM_Results/NoTreatmentSF_vs_5Aza/NoTreatmentSF_vs_5Aza_TM_candidates.csv")

# read files of filtered TM events
NT_5Aza_df <- read.csv("/private10/Projects/Efi/AML/PSI-Sigma/All/18-Hours-Treatments/TM_Results/NoTreatmentSF_vs_5Aza/NoTreatmentSF_vs_5Aza_TM_candidates.csv")
NT_FB23_df <- read.csv("/private10/Projects/Efi/AML/PSI-Sigma/All/18-Hours-Treatments/TM_Results/NoTreatmentSF_vs_FB23-2/NoTreatmentSF_vs_FB23-2_TM_candidates.csv")
# mutate data frames for "Inclusion/Exclusion" column
NT_5Aza_df <- NT_5Aza_df %>%
  mutate(`Inc/Exc` = ifelse(dPSI > 0, "Inclusion", "Exclusion"))
NT_FB23_df <- NT_FB23_df %>%
  mutate(`Inc/Exc` = ifelse(dPSI > 0, "Inclusion", "Exclusion"))
# subset data frames
cols_to_keep <- c("X", "Gene.Symbol", "Event.Region","Target.Exon","Event.Type","Exon.Type","FixedTranscript","dPSI","Inc/Exc","InclusionAAseq", "ExclusionAAseq")
NT_5Aza_df <- NT_5Aza_df[,c(cols_to_keep,colnames(NT_5Aza_df)[17], colnames(NT_5Aza_df)[18] )]
NT_FB23_df <- NT_FB23_df[,c(cols_to_keep,colnames(NT_FB23_df)[17], colnames(NT_FB23_df)[18] )]
# define column to merge by
col_to_merge <- c('FixedTranscript', 'Event.Region', 'Gene.Symbol', 'Target.Exon', 'Event.Type', 'Inc/Exc')
# merge
library(tidyverse)

# merge all data frames
list_df <- list(NT_5Aza_df, NT_FB23_df)
merged <- list_df %>% reduce(inner_join, by=col_to_merge)
#merged_df <- merge(pladB_df, indisulam_df, by <- col_to_merge)
write.csv(merged, row.names = F, file = "/private10/Projects/Efi/AML/PSI-Sigma/All/18-Hours-Treatments/TM_Results/MutualEventsCandidates/NoTreatments_mutualEvents.csv")

