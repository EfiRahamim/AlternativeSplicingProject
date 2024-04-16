import pandas as pd

file_paths = ["/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/NoTreatments/NoTreatmentNoSF_vs_NoTreatmentSF/NoTreatmentNoSF_NoTreatmentSF_merged_rank0.5_aff50.0.csv",
"/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/6-Hours-Treatments/Mock6_vs_Indisulam/Mock6_Indisulam_merged_rank0.5_aff50.0.csv",
"/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/6-Hours-Treatments/Mock6_vs_PladB/Mock6_PladB_merged_rank0.5_aff50.0.csv",
"/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/6-Hours-Treatments/NoTreatmentSF_vs_Indisulam/NoTreatmentSF_Indisulam_merged_rank0.5_aff50.0.csv",
"/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/6-Hours-Treatments/NoTreatmentSF_vs_PladB/NoTreatmentSF_PladB_merged_rank0.5_aff50.0.csv",
"/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/18-Hours-Treatments/Mock18_vs_5Aza/Mock18_5Aza_merged_rank0.5_aff50.0.csv",
"/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/18-Hours-Treatments/Mock18_vs_FB23-2/Mock18_FB23-2_merged_rank0.5_aff50.0.csv",
"/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/18-Hours-Treatments/NoTreatmentSF_vs_5Aza/NoTreatmentSF_5Aza_merged_rank0.5_aff50.0.csv",
"/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/18-Hours-Treatments/NoTreatmentSF_vs_FB23-2/NoTreatmentSF_FB23-2_merged_rank0.5_aff50.0.csv"]  

dataframes = [pd.read_csv(file) for file in file_paths]
combined_dataframe = pd.concat(dataframes, ignore_index=True)
combined_dataframe_noDups = combined_dataframe.drop_duplicates()

data = combined_dataframe_noDups
df_rank = data.melt(id_vars=['Group','Splicing Event','Peptide','ID', 'Avg.PSI_PeptideSource','Avg.PSI_OtherGroup'], 
                value_vars= data.columns[data.columns.str.endswith('_Rank')],
                var_name='HLA_Rank', 
                value_name='Rank').dropna(subset=['Rank'])
df_nM = data.melt(id_vars=['Group','Splicing Event','Peptide','ID', 'Avg.PSI_PeptideSource','Avg.PSI_OtherGroup'], 
                value_vars= data.columns[data.columns.str.endswith('_nM')],
                var_name='HLA_nM', 
                value_name='nM').dropna(subset=['nM'])

df_rank['HLA'] = df_rank['HLA_Rank'].str.split('_').str[0]
df_nM['HLA'] = df_nM['HLA_nM'].str.split('_').str[0] 

df_rank.drop(columns=['HLA_Rank'], inplace=True)
df_nM.drop(columns=['HLA_nM'], inplace=True)

# Rearrange columns in each melted DataFrame
df_rank = df_rank[['Group','Splicing Event','Peptide','ID','Avg.PSI_PeptideSource','Avg.PSI_OtherGroup', 'HLA', 'Rank']]
df_nM = df_nM[['Group','Splicing Event','Peptide','ID','Avg.PSI_PeptideSource','Avg.PSI_OtherGroup', 'HLA', 'nM']]

df_merged = pd.merge(df_rank, df_nM, on=['Group','Splicing Event','Peptide','ID','Avg.PSI_PeptideSource','Avg.PSI_OtherGroup', 'HLA'], how='inner')

# Identify peptides in both control and treatments
control_peptides = set(df_merged[(df_merged['Group'] == 'Mock6') | (df_merged['Group'] == 'NoTreatmentSF')]['Peptide'])
treatment_peptides = set(df_merged[(df_merged['Group'] == 'Indisulam') | (df_merged['Group'] == 'PladB')]['Peptide'])

# Peptides in both control and treatments
common_peptides = control_peptides.intersection(treatment_peptides)

# Filter the data frame
#filtered_df = df[~(((df['Group'] == 'Mock6') | (df['Group'] == 'NoTreatmentSF')) & (df['Peptide'].isin(common_peptides)))]
filtered_df = df_merged[~(df_merged['Peptide'].isin(common_peptides))]

# Identify peptides in both control and treatments
control_peptides = set(df_merged[(df_merged['Group'] == 'Mock18') | (df_merged['Group'] == 'NoTreatmentSF')]['Peptide'])
treatment_peptides = set(df_merged[(df_merged['Group'] == '5Aza') | (df_merged['Group'] == 'FB23-2')]['Peptide'])

# Peptides in both control and treatments
common_peptides = control_peptides.intersection(treatment_peptides)
filtered_df = filtered_df[~(filtered_df['Peptide'].isin(common_peptides))]


#filtered_df.to_csv('/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/merge_test.csv', index=False)
filtered_df.to_csv("/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/NovelStrongBindingEpitopes_noDups.csv", index=False)
#combined_dataframe_noDups.to_csv("/private10/Projects/Efi/AML/SplicingAnalysis_March2024/SplicingEvents/_forNEanalysis/NovelStrongBindingEpitopes_noDups.csv", index=False)
