results_file <- "/private10/Projects/Efi/AML/PSI-Sigma/SRSF2/6-Hours-Treatments/TM_Results/SplicingEventsFiles/NoTreatment_vs_PladB/NoTreatment_vs_PladB_analyzed.csv"
true_TM_events_file <- "/private10/Projects/Efi/AML/PSI-Sigma/SRSF2/6-Hours-Treatments/TM_Results/SplicingEventsFiles/NoTreatment_vs_PladB/SplicingEventsFiles/True_TM_events.txt"
results_df <- read.csv(results_file)
TM_events_df <- unique(read.csv(true_TM_events_file, header = F))
results_filtered <- merge(results_df, TM_events_df, by.x = 'X', by.y = 'V1')
write.csv(results_filtered, row.names = F, file = "/private10/Projects/Efi/AML/PSI-Sigma/SRSF2/6-Hours-Treatments/TM_Results/SplicingEventsFiles/NoTreatment_vs_PladB/NoTreatment_vs_PladB_TMevents.csv")


results_file <- "/private10/Projects/Efi/AML/PSI-Sigma/SRSF2/6-Hours-Treatments/TM_Results/SplicingEventsFiles/NoTreatment_vs_Indisulam/NoTreatment_vs_Indisulam_analyzed.csv"
true_TM_events_file <- "/private10/Projects/Efi/AML/PSI-Sigma/SRSF2/6-Hours-Treatments/TM_Results/SplicingEventsFiles/NoTreatment_vs_Indisulam/SplicingEventsFiles/True_TM_events.txt"
results_df <- read.csv(results_file)
TM_events_df <- unique(read.csv(true_TM_events_file, header = F))
results_filtered <- merge(results_df, TM_events_df, by.x = 'X', by.y = 'V1')
write.csv(results_filtered, row.names = F, file = "/private10/Projects/Efi/AML/PSI-Sigma/SRSF2/6-Hours-Treatments/TM_Results/SplicingEventsFiles/NoTreatment_vs_Indisulam/NoTreatment_vs_Indisulam_TMevents.csv")


pladB_df <- read.csv("/private10/Projects/Efi/AML/PSI-Sigma/SRSF2/6-Hours-Treatments/TM_Results/SplicingEventsFiles/NoTreatment_vs_PladB/NoTreatment_vs_PladB_TMevents.csv")
indisulam_df <- read.csv("/private10/Projects/Efi/AML/PSI-Sigma/SRSF2/6-Hours-Treatments/TM_Results/SplicingEventsFiles/NoTreatment_vs_Indisulam/NoTreatment_vs_Indisulam_TMevents.csv")
pladB_df <- pladB_df %>%
  mutate(`Inc/Exc` = ifelse(dPSI > 0, "Inclusion", "Exclusion"))
indisulam_df <- indisulam_df %>%
  mutate(`Inc/Exc` = ifelse(dPSI > 0, "Inclusion", "Exclusion"))
col_to_merge <- c('FixedTranscript', 'Event.Region', 'Gene.Symbol', 'Target.Exon', 'Event.Type', 'Inc/Exc')
merged_df <- merge(pladB_df, indisulam_df, by <- col_to_merge)
write.csv(merged_df, row.names = F, file = "/private10/Projects/Efi/AML/PSI-Sigma/SRSF2/6-Hours-Treatments/TM_Results/SplicingEventsFiles/NoTreatments_mutualEvents.csv")
