# run OmegaFold on inclusion and exclusion AA sequences of trans-membrane genes
import csv, os, glob, subprocess, multiprocessing, argparse
# CLI arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Run OmegaFold on inclusion and exclusion amino-acid FASTA files on trans-membrane genes. The script uses the csv file of filtered rMATS results and check for existing directories with the trans membrane genes.")
parser.add_argument("-i", action='store', dest='input_dir', required=True, help="Input directory on groups comparisons")
user_args = parser.parse_args()
input_dir = user_args.input_dir

def file_to_list(rmats_output_file):
    # open file for reading and extract the rows from the file
    with open(rmats_output_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = list(reader)
    return rows

def run_OmegaFold(fasta_file):
    omegafold_command = ['omegafold', fasta_file, "."]
    print(f"Running OmegaFold on {fasta_file}..")
    try:
        output = subprocess.check_output(omegafold_command, universal_newlines=True)
        print(output)
    except:
        print(f"Error in running omega fold command for {fasta_file}")

def get_rMATS_filtered_file(AS_dir):
    file = glob.glob(os.path.join(AS_dir, "*_filtered_byTM.csv"))[0]
    return file

def findExistGenes(rows, genes_dir):
    id_directories_pattern = os.path.join(genes_dir,"*","*")
    id_directories = glob.glob(id_directories_pattern)
    id_pathes = [os.path.abspath(dir) for dir in id_directories if os.path.isdir(dir)]
    exist_dirs = []
    # find exist directories of TM genes (be their ID number of splicing event)
    for row in rows:
        id_dir_path = os.path.join(genes_dir, f"{row['GeneID']}_{row['geneSymbol']}", row['ID'])
        if id_dir_path in id_pathes:
            exist_dirs.append(id_dir_path)
    print(f"Directories of trans-membrane genes: {exist_dirs}")
    return exist_dirs


# find AS directories
AS_dirs = [] 
AS_dirs_pattern = ['SE', 'A5SS', 'A3SS', 'MXE', 'RI']
for as_type in AS_dirs_pattern:
    AS_dirs.append(glob.glob(os.path.join(input_dir, as_type))[0])
for as_dir in AS_dirs:
    rmats_output = get_rMATS_filtered_file(as_dir)
    # read rMATS output file
    rows = file_to_list(rmats_output)
    # find existing genes directory
    exist_dirs = findExistGenes(rows, as_dir)
    # run over each gene
    for id_dir in exist_dirs:
        # find transcripts directories
        transcripts_dir_pattern = os.path.join(id_dir, "*")
        transcripts_dirs = glob.glob(transcripts_dir_pattern)
        # run over each transcript
        for transcript_dir in transcripts_dirs:
            # OmegaFold directory creation
            omegafold_dir = os.path.join(transcript_dir, "OmegaFold")
            if not os.path.isdir(omegafold_dir):
                os.mkdir(omegafold_dir)
            os.chdir(omegafold_dir)
            # get AA fasta files (inclusion and exclusion)
            AA_fasta_files = glob.glob(os.path.join(transcript_dir, "*AA*.fasta"))
            # run OmegaFold in parallel on both FASTA files
            pool = multiprocessing.Pool(processes=2)
            results = pool.map(run_OmegaFold, AA_fasta_files)
            pool.close()
            pool.join()

        






