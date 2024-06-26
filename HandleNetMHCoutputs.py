import pandas as pd
import os, glob,argparse
import regex as re
# CLI arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Filter NetMHC output by StrongBinding values and merge outputs together to one file. Also filtering out mutual peptides that appear in both group.\nNOTE: For Avg. PSI calculation make sure the file '*_analyzed.csv' is located in the input directory!.")
parser.add_argument("-i", action='store', dest='in_dir', required=True, help="Input directory of comparisons groups")
parser.add_argument("-o", action='store', dest='output_dir', required=True, help="Output directory for merged file to be written")
parser.add_argument("-l1", action='store', dest='group1_name', required=True, help="Name of group 1")
parser.add_argument("-l2", action='store', dest='group2_name', required=True, help="Name of group 2")
parser.add_argument("--rank", action='store', dest='rank', required=False, type=float, default=0.5, help="Rank value for filtering (default: 0.5)")
parser.add_argument("--aff", action='store', dest='aff', required=False, type=float, default=float('inf'), help="Affinity (nM) value for filtering (default: None")
parser.add_argument('--PSIsigma', action='store_true', dest='psi_sigma', help='Set for PSI-Sigma tool results.')
user_args = parser.parse_args()

# Parse CLI arguments
in_dir = user_args.in_dir
output_dir = user_args.output_dir
group1_name = user_args.group1_name
group2_name = user_args.group2_name
rank = user_args.rank
aff = user_args.aff

# DEBUG Arguments
# in_dir = "/private10/Projects/Efi/AML/rMATS/Mock6h_Indisulam/"
# output_dir = "/private10/Projects/Efi/AML/NeoEpitopesAnalyze/"
# group1_name = "Mock_6h"
# group2_name = "Indisulam"

# list of HLA types
HLA_types = ["HLA-A0101","HLA-A0201","HLA-A0301","HLA-A1101", "HLA-A2402","HLA-A2601","HLA-B0702","HLA-B0801","HLA-B1501","HLA-B2705","HLA-B3901","HLA-B4001","HLA-B5801"]
# header of data frame
unified_col = ['Pos', 'Peptide', 'ID',
               'HLA-A0101_nM','HLA-A0101_Rank','HLA-A0101_Core',
                'HLA-A0201_nM','HLA-A0201_Rank','HLA-A0201_Core',
                 'HLA-A0301_nM','HLA-A0301_Rank','HLA-A0301_Core',
                  'HLA-A1101_nM','HLA-A1101_Rank','HLA-A1101_Core',
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
if user_args.psi_sigma:
    group1_files_pattern = os.path.join(in_dir,'SplicingEventsFiles', "*", "*", "*", "*", "netMHC*", f"*{group1_name}.xls")
    group2_files_pattern = os.path.join(in_dir,'SplicingEventsFiles', "*", "*", "*", "*", "netMHC*", f"*{group2_name}.xls")
else:
    group1_files_pattern = os.path.join(in_dir, "*", "*", "*", "*", "netMHC*", f"*{group1_name}.xls")
    group2_files_pattern = os.path.join(in_dir, "*", "*", "*", "*", "netMHC*", f"*{group2_name}.xls")
group1_files = glob.glob(group1_files_pattern)
group2_files = glob.glob(group2_files_pattern)
group1_merged = pd.DataFrame()
group2_merged = pd.DataFrame()
# get and read the 'analyzed.csv' file for adding the Avg.PSI values
analyzed_file_pattern = os.path.join(in_dir, '*_analyzed.csv') # check for the '*_azanlyzed.csv' file
matching_files = glob.glob(analyzed_file_pattern)
if not matching_files:
    raise FileNotFoundError(f"No file matching the pattern {analyzed_file_pattern}")
analyzed_file = matching_files[0]
analyzed_df = pd.read_csv(analyzed_file)
# filter each file by rank and affinity score and merge them all (for each group)
for group1_file, group2_file in zip(group1_files, group2_files):
    # get type of splicing event
    if user_args.psi_sigma:
        splicing_event = group1_file.strip("/").split("/")[-4]
    else:
        splicing_event_pattern = fr"{in_dir}(\w+)/"
        splicing_event = re.search(splicing_event_pattern, group1_file).group(1)
    # read file and change header
    try:
        df_1 = pd.read_csv(group1_file, sep='\t', header = 1, index_col=False)
        df_2 = pd.read_csv(group2_file, sep='\t', header = 1, index_col=False)
    except Exception as e:
        print(f"Error occurred: {e}")
        print(f"Files: {group1_file}, {group2_file}\nMaybe files are empty. Skipping.")        
        continue
    df_1.columns = unified_col
    df_2.columns = unified_col
    df_1.insert(0, "Splicing Event", splicing_event)
    df_2.insert(0, "Splicing Event", splicing_event)
    # subset data frame by strong binding (sb) value: rank and affinity scores
    columns_to_check_rank = ['HLA-A0101_Rank', 'HLA-A0201_Rank', 'HLA-A0301_Rank', 'HLA-A1101_Rank',
                         'HLA-A2402_Rank', 'HLA-A2601_Rank', 'HLA-B0702_Rank', 'HLA-B0801_Rank',
                         'HLA-B1501_Rank', 'HLA-B2705_Rank', 'HLA-B3901_Rank', 'HLA-B4001_Rank',
                         'HLA-B5801_Rank']
    columns_to_check_nM = ['HLA-A0101_nM', 'HLA-A0201_nM', 'HLA-A0301_nM', 'HLA-A1101_nM',
                        'HLA-A2402_nM', 'HLA-A2601_nM', 'HLA-B0702_nM', 'HLA-B0801_nM',
                        'HLA-B1501_nM', 'HLA-B2705_nM', 'HLA-B3901_nM', 'HLA-B4001_nM',
                        'HLA-B5801_nM']
    df_1_filteredSB = df_1[(df_1[columns_to_check_rank] <= rank).any(axis=1) &
                        (df_1[columns_to_check_nM] <= aff).any(axis=1)]
    df_2_filteredSB = df_2[(df_2[columns_to_check_rank] <= rank).any(axis=1) &
                        (df_2[columns_to_check_nM] <= aff).any(axis=1)]
    # Add the splicing index, Reference Transcript, PSI and TPM values of the specific splicing event in each group
    if user_args.psi_sigma:
        # add the splicing index
        splicing_index = int(group1_file.strip("/").split("/")[-3]) # get the index of the splicng event according to the 'analyzed' file.
        df_1_filteredSB.insert(1, 'SplicingIndex', splicing_index) # splicing index of event
        df_2_filteredSB.insert(1, 'SplicingIndex', splicing_index) # splicing index of event
        # add the refernece transcript
        ref_transcript = analyzed_df.loc[splicing_index, 'Reference.Transcript'] # get the reference transcript of the splicing event
        df_1_filteredSB.insert(1, 'Reference.Transcript', ref_transcript) # reference transcript of event
        df_2_filteredSB.insert(1, 'Reference.Transcript', ref_transcript) # reference transcript of event
        # add the Avg.PSI values of the groups
        df_1_filteredSB.insert(1, 'Avg.PSI_PeptideSource', analyzed_df.loc[splicing_index, f'Avg.PSI_{group1_name}']) # PSI of groupA of splicing (peptide source)
        df_1_filteredSB.insert(2, 'Avg.PSI_OtherGroup', analyzed_df.loc[splicing_index, f'Avg.PSI_{group2_name}']) # PSI of groupB of splicing (not peptide source)
        df_2_filteredSB.insert(1, 'Avg.PSI_PeptideSource', analyzed_df.loc[splicing_index, f'Avg.PSI_{group2_name}']) # PSI of groupA of splicing (peptide source)
        df_2_filteredSB.insert(2, 'Avg.PSI_OtherGroup', analyzed_df.loc[splicing_index, f'Avg.PSI_{group1_name}']) # PSI of groupB of splicing (not peptide source)
        # add the TPM values of the groups
        df_1_filteredSB.insert(1, 'Avg.TPM_PeptideSource', analyzed_df.loc[splicing_index, f'TPM_mean_{group1_name}']) # TPM of groupA of splicing (peptide source)
        df_1_filteredSB.insert(2, 'Avg.TPM_OtherGroup', analyzed_df.loc[splicing_index, f'TPM_mean_{group2_name}']) # TPM of groupB of splicing (not peptide source)
        df_2_filteredSB.insert(1, 'Avg.TPM_PeptideSource', analyzed_df.loc[splicing_index, f'TPM_mean_{group2_name}']) # TPM of groupA of splicing (peptide source)
        df_2_filteredSB.insert(2, 'Avg.TPM_OtherGroup', analyzed_df.loc[splicing_index, f'TPM_mean_{group1_name}']) # TPM of groupB of splicing (not peptide source)
    
    if not group1_merged.empty and not group2_merged.empty:
        group1_merged = pd.concat([group1_merged, df_1_filteredSB], ignore_index=True)
        group2_merged = pd.concat([group2_merged, df_2_filteredSB], ignore_index=True)
    else:
        if group1_merged.empty:
            group1_merged = df_1_filteredSB
        if group2_merged.empty:
            group2_merged = df_2_filteredSB

# seperate each data frame by the HLA alleles and make symmetric difference between groups
merged_all = pd.DataFrame()
for HLA_type in HLA_types:
    cols_to_keep = ["Splicing Event", "Peptide", "ID",'SplicingIndex','Reference.Transcript', 'Avg.PSI_PeptideSource','Avg.PSI_OtherGroup','Avg.TPM_PeptideSource','Avg.TPM_OtherGroup', f'{HLA_type}_Rank', f'{HLA_type}_nM']
    # filter data frame 1 by rank and affinity values at current HLA allele
    sb_group1 = group1_merged.loc[(group1_merged[f'{HLA_type}_Rank'] <=rank) & (group1_merged[f'{HLA_type}_nM'] <= aff), cols_to_keep]
    sb_group1.insert(0, "Group", group1_name)
    # filter data frame 2 by rank and affinity values at current HLA allele
    sb_group2 = group2_merged.loc[(group2_merged[f'{HLA_type}_Rank'] <=rank) & (group2_merged[f'{HLA_type}_nM'] <= aff), cols_to_keep]
    sb_group2.insert(0, "Group", group2_name)
    # find the symetric difference between the data frames
    symmetric_diff_df = pd.concat([sb_group1[~sb_group1['Peptide'].isin(sb_group2['Peptide'])], sb_group2[~sb_group2['Peptide'].isin(sb_group1['Peptide'])]])
    # merge the symetric difference data frame to he previous ones
    if merged_all.empty:
        merged_all = symmetric_diff_df
    else:
        merged_all = pd.concat([merged_all, symmetric_diff_df], ignore_index=True)
# drop duplicates rows
merged_all.drop_duplicates()
output_file = os.path.join(output_dir, f"{group1_name}_{group2_name}_merged_rank{rank}_aff{aff}.csv")
merged_all.to_csv(output_file, index=False)
print("Done.")
