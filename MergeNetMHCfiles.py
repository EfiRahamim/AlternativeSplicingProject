import pandas as pd

file_paths = ["/private10/Projects/Efi/AML/NeoEpitopesAnalyze/5-Azacytidine_FB23-2_merged.csv",
"/private10/Projects/Efi/AML/NeoEpitopesAnalyze/Control_SFmutations_5-Azacytidine_merged.csv",
"/private10/Projects/Efi/AML/NeoEpitopesAnalyze/Control_SFmutations_FB23-2_merged.csv",
"/private10/Projects/Efi/AML/NeoEpitopesAnalyze/Control_SFmutations_Indisulam_merged.csv",
"/private10/Projects/Efi/AML/NeoEpitopesAnalyze/Control_SFmutations_Mock_18h_merged.csv",
"/private10/Projects/Efi/AML/NeoEpitopesAnalyze/Control_SFmutations_Mock_6h_merged.csv",
"/private10/Projects/Efi/AML/NeoEpitopesAnalyze/Control_SFmutations_Pladienolide-B_merged.csv",
"/private10/Projects/Efi/AML/NeoEpitopesAnalyze/Indisulam_Pladienolide-B_merged.csv",
"/private10/Projects/Efi/AML/NeoEpitopesAnalyze/Mock_18h_5-Azacytidine_merged.csv",
"/private10/Projects/Efi/AML/NeoEpitopesAnalyze/Mock_18h_FB23-2_merged.csv",
"/private10/Projects/Efi/AML/NeoEpitopesAnalyze/Mock_6h_Indisulam_merged.csv",
"/private10/Projects/Efi/AML/NeoEpitopesAnalyze/Mock_6h_Pladienolide-B_merged.csv"
]  

dataframes = [pd.read_csv(file) for file in file_paths]
combined_dataframe = pd.concat(dataframes, ignore_index=True)
combined_dataframe_noDups = combined_dataframe.drop_duplicates()

combined_dataframe.to_csv("/private10/Projects/Efi/AML/rMATS/NovelStrongBindingEpitopes.csv", index=False)
