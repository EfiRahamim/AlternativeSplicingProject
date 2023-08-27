# run the trans-membrane prediction step on existing genes directories
# 1. Create PDB file for each sequence
# 2. run MembraneFold of PDB file

import csv, os, glob, subprocess

rmats_output = "/private10/Projects/Efi/AML/rMATS/Mock6h_Pladienolide/SE/SE_filtered_byTM.csv"
genes_dir = "/private10/Projects/Efi/AML/tests/SE_TestForOmegaFold/"
def file_to_list(rmats_output_file):
    # open file for reading and extract the rows from the file
    with open(rmats_output_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = list(reader)
    return rows

# read rMATS output file
rows = file_to_list(rmats_output)
directories_pattern = os.path.join(genes_dir,"*")
directories = glob.glob(directories_pattern)
pathes = [os.path.abspath(dir) for dir in directories if os.path.isdir(dir)]
exist_dirs = []
# find exist directories of TM genes
for row in rows:
    dir_path = os.path.join(genes_dir, f"{row['GeneID']}_{row['geneSymbol']}")
    if dir_path in pathes:
        exist_dirs.append(dir_path)
print(f"Directories of trans-membrane genes: {exist_dirs}")
for gene_dir in exist_dirs:
    transcripts_dir_pattern = os.path.join(gene_dir, "*","*")
    directories = glob.glob(transcripts_dir_pattern)
    transcripts_dirs = [os.path.abspath(dir) for dir in directories if os.path.isdir(dir)]
    for transcript_dir in transcripts_dirs:
        membranefold_dir = os.path.join(transcript_dir, "MembraneFold")
        if not os.path.isdir(membranefold_dir):
            os.mkdir(membranefold_dir)
        os.chdir(membranefold_dir)
        files = glob.glob(os.path.join(transcript_dir, "*"))
        print(files)
        for file in files:
            if "inclusionAA" in file or "exclusionAA" in file:
                print(f"Running OmegaFold in {file}..")
                omegafold_command = ['omegafold', file, membranefold_dir]
                try:
                    output = subprocess.check_output(omegafold_command, universal_newlines=True)
                    print(output)
                except:
                    print(f"Error in running omega fold command for {file}")
        






