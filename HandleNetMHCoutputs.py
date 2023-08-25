import pandas as pd
import os, glob
# input directory and group names
in_dir = "/private10/Projects/Efi/AML/rMATS/Control_PladienolideB/"
group1_name = "Pladienolide-B"
group2_name = "Control_SFmutations"
# list of HLA types
HLA_types = ["HLA-A0101","HLA-A0201","HLA-A0301","HLA-A2402","HLA-A2601","HLA-B0702","HLA-B0801","HLA-B1501","HLA-B2705","HLA-B3901","HLA-B4001","HLA-B5801"]
# header of data frame
unified_col = ['Pos', 'Peptide', 'ID',
               'HLA-A0101_nM','HLA-A0101_Rank','HLA-A0101_Core',
                'HLA-A0201_nM','HLA-A0201_Rank','HLA-A0201_Core',
                 'HLA-A0301_nM','HLA-A0301_Rank','HLA-A0301_Core',
                  'HLA-A2402_nM','HLA-A2402_Rank','HLA-A2402_Core',
                   'HLA-A2601_nM','HLA-A2601_Rank','HLA-A2601_Core',
                    'HLA-B0702_nM','HLA-B0702_Rank','HLA-B0702_Core',
                     'HLA-B0801_nM','HLA-B0801_Rank','HLA-B0801_Core',
                      'HLA-B1501_nM','HLA-B1501_Rank','HLA-B1501_Core',
                       'HLA-B2705_nM','HLA-B2705_Rank','HLA-B2705_Core',
                        'HLA-B3901_nM','HLA-B3901_Rank','HLA-B3901_Core',
                         'HLA-B4001_nM','HLA-B4001_Rank','HLA-B4001_Core',
                          'HLA-B5801_nM','HLA-B5801_Rank','HLA-B5801_Core',
                           "H_Avg_Ranks", "N_binders"]
# get all netMHC excel files of group1 and group2
group1_files_pattern = os.path.join(in_dir, "*", "*", "*", "*", "netMHC*", f"*{group1_name}.xls")
group2_files_pattern = os.path.join(in_dir, "*", "*", "*", "*", "netMHC*", f"*{group2_name}.xls")
group1_files = glob.glob(group1_files_pattern)
group2_files = glob.glob(group2_files_pattern)
group1_merged = pd.DataFrame()
group2_merged = pd.DataFrame()
# filter each file by rank < 0.5 and merge them all (for each group)
for group1_file, group2_file in zip(group1_files, group2_files):
    # read file and change header
    df_1 = pd.read_csv(group1_file, sep='\t', header = 1, index_col=False)
    df_2 = pd.read_csv(group2_file, sep='\t', header = 1, index_col=False)
    df_1.columns = unified_col
    df_2.columns = unified_col
    # subset data frame by strong binding (sb) value (rank < 0.5)
    df_1_filteredSB = df_1[(df_1['HLA-A0101_Rank'] <= 0.5) | 
               (df_1['HLA-A0201_Rank'] <= 0.5) |
               (df_1['HLA-A0301_Rank'] <= 0.5) |
               (df_1['HLA-A2402_Rank'] <= 0.5) |
               (df_1['HLA-A2601_Rank'] <= 0.5) |
               (df_1['HLA-B0702_Rank'] <= 0.5) |
               (df_1['HLA-B0801_Rank'] <= 0.5) |
               (df_1['HLA-B1501_Rank'] <= 0.5) |
               (df_1['HLA-B2705_Rank'] <= 0.5) |
               (df_1['HLA-B3901_Rank'] <= 0.5) |
               (df_1['HLA-B4001_Rank'] <= 0.5) |
               (df_1['HLA-B5801_Rank'] <= 0.5)]
    df_2_filteredSB = df_2[(df_2['HLA-A0101_Rank'] <= 0.5) | 
            (df_2['HLA-A0201_Rank'] <= 0.5) |
            (df_2['HLA-A0301_Rank'] <= 0.5) |
            (df_2['HLA-A2402_Rank'] <= 0.5) |
            (df_2['HLA-A2601_Rank'] <= 0.5) |
            (df_2['HLA-B0702_Rank'] <= 0.5) |
            (df_2['HLA-B0801_Rank'] <= 0.5) |
            (df_2['HLA-B1501_Rank'] <= 0.5) |
            (df_2['HLA-B2705_Rank'] <= 0.5) |
            (df_2['HLA-B3901_Rank'] <= 0.5) |
            (df_2['HLA-B4001_Rank'] <= 0.5) |
            (df_2['HLA-B5801_Rank'] <= 0.5)]
    if group1_merged.empty:
        group1_merged = df_1_filteredSB
    elif group2_merged.empty:
        group2_merged = df_2_filteredSB
    else:
        group1_merged = pd.concat([group1_merged, df_1_filteredSB], ignore_index=True)
        group2_merged = pd.concat([group2_merged, df_2_filteredSB], ignore_index=True)
# seperate each data frame by the HLA alleles and make symmetric difference between groups
merged_all = pd.DataFrame()
for HLA_type in HLA_types:
    # filter data frame 1 by rank value and current HLA allele
    sb_group1 = group1_merged.loc[group1_merged[f'{HLA_type}_Rank'] <=0.5, ["Pos", "Peptide", "ID",f'{HLA_type}_Rank']]
    sb_group1.insert(0, "Group", group1_name)
    # filter data frame 2 by rank value and current HLA allele
    sb_group2 = group2_merged.loc[group2_merged[f'{HLA_type}_Rank'] <=0.5, ["Pos", "Peptide", "ID",f'{HLA_type}_Rank']]
    sb_group2.insert(0, "Group", group2_name)
    # find the symetric difference between the data frames
    symmetric_diff_df = pd.concat([sb_group1[~sb_group1['Peptide'].isin(sb_group2['Peptide'])], sb_group2[~sb_group2['Peptide'].isin(sb_group1['Peptide'])]])
    
    if merged_all.empty:
        merged_all = symmetric_diff_df
    else:
        merged_all = pd.concat([merged_all, symmetric_diff_df], ignore_index=True)

group_counts = merged_all.groupby('Group').size()
print(group_counts)

# netMHC output
group_1_output = "/private10/Projects/Efi/AML/rMATS/Control_PladienolideB/SE/ENSG00000282218.1_RP1-179P9.3/260616/ENST00000052569.10/netMHC_Rank0.5/netMHC_exclusionAA_Pladienolide-B.xls"
group_2_output = "/private10/Projects/Efi/AML/rMATS/Control_PladienolideB/SE/ENSG00000282218.1_RP1-179P9.3/260616/ENST00000052569.10/netMHC_Rank0.5/netMHC_inclusionAA_Control_SFmutations.xls"
group1_name = "Pladienolide-B"
group2_name = "Control_SFmutations"
# read files
df_1 = pd.read_csv(group_1_output, sep='\t', header = 1, index_col=False)
df_2 = pd.read_csv(group_2_output, sep='\t', header = 1, index_col=False)
df_1.columns = unified_col
df_2.columns = unified_col
# subset data frames by strong binding (sb) value
sb_df_1 = df_1[(df_1['HLA-A0101_Rank'] <= 0.5) | 
               (df_1['HLA-A0201_Rank'] <= 0.5) |
               (df_1['HLA-A0301_Rank'] <= 0.5) |
               (df_1['HLA-A2402_Rank'] <= 0.5) |
               (df_1['HLA-A2601_Rank'] <= 0.5) |
               (df_1['HLA-B0702_Rank'] <= 0.5) |
               (df_1['HLA-B0801_Rank'] <= 0.5) |
               (df_1['HLA-B1501_Rank'] <= 0.5) |
               (df_1['HLA-B2705_Rank'] <= 0.5) |
               (df_1['HLA-B3901_Rank'] <= 0.5) |
               (df_1['HLA-B4001_Rank'] <= 0.5) |
               (df_1['HLA-B5801_Rank'] <= 0.5)]
merged_all = pd.DataFrame()
for HLA_type in HLA_types:
    # filter data frame 1 by rank value
    sb_group1 = df_1.loc[df_1[f'{HLA_type}_Rank'] <=0.5, ["Pos", "Peptide", "ID",f'{HLA_type}_nM',f'{HLA_type}_Rank',f'{HLA_type}_Core' ]]
    sb_group1.insert(0, "Group", group1_name)
    # filter data frame 1 by rank value
    sb_group2 = df_2.loc[df_2[f'{HLA_type}_Rank'] <=0.5, ["Pos", "Peptide", "ID",f'{HLA_type}_nM',f'{HLA_type}_Rank',f'{HLA_type}_Core' ]]
    sb_group2.insert(0, "Group", group2_name)
    # find the symetric difference between the data frames
    symmetric_diff_df = pd.concat([sb_group1[~sb_group1['Peptide'].isin(sb_group2['Peptide'])], sb_group2[~sb_group2['Peptide'].isin(sb_group2['Peptide'])]])
    if merged_all.empty:
        merged_all = symmetric_diff_df
    else:
        merged_all = pd.concat([merged_all, symmetric_diff_df])

print()
