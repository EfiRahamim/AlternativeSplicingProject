import os
import biolib
from Bio import Entrez, SeqIO
from Bio import Seq
import shutil
import pandas as pd
from pybedtools import BedTool

### parsing arguments
# def parse_args():
#     parser = argparse.ArgumentParser(
#         description=('add a prefix to each file as it is copied'
#                      ' into a destination directory'))
#     parser.add_argument('prefix', help='the prefix to add to each file')
#     parser.add_argument('dest_dir',
#                         help='the destination directory to copy files to')
#     parser.add_argument('files', nargs='+', help='the files to copy')
#     return parser.parse_args()



### global variables

# rMATS outputs - wanted columns for BED files
SE_wanted_col = [3,5,6,4]
RI_wanted_col = [3,5,6,4]
MXE_1_wanted_col = [3,5,6,4]
MXE_2_wanted_col = [3,7,8,4]
A5SS_wanted_col = [3,] # TBD
A3SS_wanted_col = [] # TBD

### REQUIERD ARGUMENTS FROM COMMAND LINE ###
results_dir = "/private5/Projects/Efi/AS/pipeline_tests/testRun_results"
genome_fasta_file = "/private/dropbox/Genomes/Human/hg38/hg38.fa"
# list of genes that has significant exon(s)
#path_of_genes_list = "/private5/Projects/Efi/AS/pipeline_tests/gene_list.txt"
path_of_genes_list = "/private5/Projects/Efi/AS/TCGA-BRCA/test/output/SE.MATS.JCEC.significant_genes.txt"
#rmats_output = "/private5/Projects/Efi/Carcinoma_HCRN/rmats/output/SE.MATS.JC.significant.txt"
rmats_output = "/private5/Projects/Efi/AS/TCGA-BRCA/test/output/SE.MATS.JCEC.significant.txt"
# event = SE/RI/MXE/A5SS/A3SS

# lists to hold TM and not TM genes
TM_genes =[]
notTM_genes=[]
# list of genes with AS exon that does not devide by 3
notDevideByThree =[]




"""
This function predict trans-membrane (TM) domains for proteins, using DeepTMHMM tool.
Input: directory of DeepTMHMM results.
output: If TM domain exist - TRUE
        Else - FALSE + delete directory of results
"""
def check_TM(result_dir_path):
    # check if TM domain(s) exist
    result_file_path = os.path.join(result_dir_path,"deeptmhmm_results.md")
    if not os.path.exists(result_file_path):
        print(f'File {result_file_path} does not exist.')
        return
    print("check for TMR's.")
    # open result file
    with open(result_file_path, "r") as result_file:
        for line in result_file:
            if "Number of predicted TMRs: " in line:
                tmr_count_str = line.split("Number of predicted TMRs: ")[-1].strip()
                print(f"found {tmr_count_str} TMR's.")
                break
    # check for TMR's
    if tmr_count_str.isdigit():
        tmr_count = int(tmr_count_str)
        if not tmr_count > 0:
            return False
        else:
            return True

def get_gene_NA_seq(gene_id,gene_dir, type):
    # Creat NA sequence in FASTA
    print("Fetching: ", gene_id,"...")
    handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="fasta_cds_na", retmode="text")
    record_na = handle.read()
    na_fasta_path = os.path.join(gene_dir,gene_id+"_fullNA_"+type+".fasta")
    with open(na_fasta_path, 'w') as na_fasta_file:
        na_fasta_file.write(record_na)
    print("Fetching ", gene, " Completed.")
    return na_fasta_path

def get_gene_AA_seq(na_fasta_path,gene_id, gene_dir,type, spliced=False):
    # Create AA sequence in FASTA
    print("Translating: ", gene_id, "...")
    dna_seq = SeqIO.read(na_fasta_path, "fasta")
    aa_seq = dna_seq.seq.translate()
    record_aa = ">" + dna_seq.id + "\n" + str(aa_seq) + "\n"
    if spliced:
        aa_fasta_path = os.path.join(gene_dir,gene_id+"_splicedAA_"+type+".fasta")
    else:
        aa_fasta_path = os.path.join(gene_dir,gene_id+"_fullAA_"+type+".fasta")
    with open(aa_fasta_path, 'w') as aa_fasta_file:
        aa_fasta_file.write(record_aa)
    print("Translating ",gene, " Completed.")
    return aa_fasta_path

def make_dirs():
    # create directory for all NA sequences
    if not os.path.exists("./Full_NA_Seq"):
        os.mkdir("./Full_NA_Seq")
    # create directory for AA sequences
    if not os.path.exists("./Full_AA_Seq"):
        os.mkdir("./Full_AA_Seq")
    # create directory for DeepTMHMM results
    if not os.path.exists("./DeepTMHMM"):
        os.mkdir("./DeepTMHMM")
    # create directory for spliced NA sequences
    if not os.path.exists("./Spliced_NA_Seq"):
        os.mkdir("./Spliced_NA_Seq")
    # create directory for spliced AA sequences
    if not os.path.exists("./Spliced_AA_Seq"):
        os.mkdir("./Spliced_AA_Seq")
    
def run_DeepTMHMM(deeptmhmm, aa_path, gene_id, gene_dir, type,spliced=False):
    # create command for the tool
    cmd = '--fasta '+aa_path
    print("Sending command to deeptmhmm.")
    # run DeepTMHMM
    deeptmhmm_job = deeptmhmm.cli(args=cmd)
    # create result dir according to gene ID and save results
    if spliced:
        result_dir_path = os.path.join(gene_dir,"DeepTMHMM/"+gene_id+"_spliced_"+type)
    else:
        result_dir_path = os.path.join(gene_dir,"DeepTMHMM/"+gene_id+"_full_"+type)
    print(f"Saving results in {result_dir_path}")
    deeptmhmm_job.save_files(result_dir_path)
    return result_dir_path

def create_bed_file(rmats_file_path, wanted_col, gene, gene_dir, AS_type):
    # read rMATS results file
    rmats_df = pd.read_csv(rmats_file_path, sep="\t")
    # select subset of data, according to the wanted columns list (given) and wanted gene (given)
    sub_rmats_df = rmats_df.loc[rmats_df["GeneID"]==gene]
    sub_rmats_df = sub_rmats_df.iloc[:,wanted_col]
    # write output to '.bed' file
    output_path = os.path.join(gene_dir, f"{AS_type}_AS_exonsRegions.bed")
    print(f"Writing AS Exons Regions to:{output_path}")
    with open (output_path, 'w') as output_bed_file:
        sub_rmats_df.to_csv(output_path, sep = '\t', header=False, index=False)
    return output_path

def create_spliced_seq(na_fasta_file, bed_file_path, gene, gene_dir, type):
    print(f"Creating spliced NA sequence for {gene}...")
    # get the sequence of the AS exon
    bed_file = BedTool(bed_file_path)
    fasta_file = BedTool(genome_fasta_file)
    sequences = bed_file.sequence(fi=fasta_file)
    exon_fasta = SeqIO.read(sequences.seqfn, 'fasta')
    with open(bed_file_path, 'r') as b:
        line = b.readline()
        strand = line[-2] 
    # check if exon on positive or negative strand
    if strand == '+':
        exon_seq = exon_fasta.seq
    elif strand == "-":
        exon_seq = exon_fasta.reverse_complement().seq
    else:
        print(f"Strand of {gene} is not recognized. Skipping.")
        return
    # check if exon length is devided by 3 (realted to ORF consequences)
    if len(exon_seq)%3==0:
        print(f"Exon length of {gene} is devided by 3. Not hurting ORF.")
        full_seq_fasta = SeqIO.read(na_fasta_path, "fasta")
        full_seq = full_seq_fasta.seq
        exon_seq = exon_seq.upper() # make sure all the sequence is in upper case (relevant for repeatings sequences)
        spliced_seq = full_seq.replace(exon_seq,'')
        # make sure the splicing succeed
        if len(spliced_seq) == len(full_seq):
            print(f"Splicing was not done correctly. Skipping.")
            return
        record_spliced_na = ">" + full_seq_fasta.id + "\n" + str(spliced_seq) + "\n"
        na_spliced_fasta_path = os.path.join(gene_dir,gene+"_splicedNA_"+type+".fasta")
        with open(na_spliced_fasta_path, 'w') as na_spliced_fasta_file:
            na_spliced_fasta_file.write(record_spliced_na)
        print(f"Splicing of {gene} has completed. Sequence has been saved at:{na_spliced_fasta_path}")
        return na_spliced_fasta_path
    else:
        print(f"AS exon of gene {gene} is not devided by 3. Skipping.")
        # delete gene directory
        print(f"Deleting {gene_dir}...")
        shutil.rmtree(gene_dir, onerror=None)
        #notDevideByThree.append(gene)
        with open ("NotDevideBy3Exons.txt", 'a') as f:
            f.write(gene+'\n')
        return   




Entrez.email = 'efirahamim@gmail.com'
genes_list=[]
with open(path_of_genes_list, 'r') as genes_file:
    for line in genes_file:
        values=line.strip().split('\t')
        genes_list.append([values[0], float(values[1])])
# remove double quets from the gene ID's
genes_list = [[gene[0].replace('"',''), gene[1]] for gene in genes_list]
# load DeepTMHMM tool
deeptmhmm = biolib.load('DTU/DeepTMHMM')
#make_dirs()
# create FASTA file of NA and AA sequence for each gene that has TM domain
for gene_and_score in genes_list:
    gene = gene_and_score[0]
    # set the 'type' list - normal and tumor. If IncAvgLevel is positive - [normal, tumor]. else - [tumor, normal]
    if gene_and_score[1] > 0:
        type = ['normal','tumor']
    elif gene_and_score[1] < 0:
        type = ['tumor', 'normal']
    else:
        print(f"No significant IncAvgLevel in {gene}. Skipping.")
        continue
    # create directory for gene
    gene_dir = os.path.join(results_dir, gene)
    # check if current gene was already checked in previous running
    if os.path.isdir(gene_dir):
        print(f"Gene {gene} was already checked.\nResults can be found in {gene_dir}.\nSkipping to next gene.")
        continue
    else:
        os.mkdir(gene_dir)
    # get nucleic acid sequence (FASTA)
    na_fasta_path = get_gene_NA_seq(gene, gene_dir,type[0])
    #na_fasta_path = "/private5/Projects/Efi/AS/pipeline_tests/NM_001077628.3/NM_001077628.3_fullNA.fasta"
    # get amino acid sequence (FASTA)
    aa_fasta_path = get_gene_AA_seq(na_fasta_path, gene, gene_dir,type[0])
    #aa_fasta_path = "/private5/Projects/Efi/AS/pipeline_tests/NM_001077628.3/NM_001077628.3_fullAA.fasta"
    # Check if gene has TM domain(s) - using DeepTMHMM tool
    deeptmhmm_result_path = run_DeepTMHMM(deeptmhmm,aa_fasta_path,gene, gene_dir,type[0])
    #deeptmhmm_result_path = "/private5/Projects/Efi/AS/pipeline_tests/NM_001077628.3/DeepTMHMM"
    is_TM = check_TM(deeptmhmm_result_path)
    if is_TM:
        TM_genes.append(gene)
        # make bed file of AS exons regions for this genes
        bed_file = create_bed_file(rmats_output,SE_wanted_col,gene, gene_dir)
        # spliced the AS exon from the full NA sequence
        na_spliced_seq_path = create_spliced_seq(na_fasta_path, bed_file, gene, gene_dir,type[1])
        # continue only if exon is devided by three
        if na_spliced_seq_path != None:
            # create the spliced AA sequence
            aa_spliced_seq_path = get_gene_AA_seq(na_spliced_seq_path, gene, gene_dir,type[1] ,spliced=True)
            # run DeepTMHMM on the spliced protein sequence
            spliced_prot_deeptmhmm = run_DeepTMHMM(deeptmhmm,aa_spliced_seq_path,gene,gene_dir, type[1],spliced=True )
    else:
        #notTM_genes.append(gene)
        with open ("NoTMGenes.txt", 'a') as f:
            f.write(gene + '\n')
        # delete gene directory
        print(f"Deleting {gene_dir}...")
        shutil.rmtree(gene_dir, onerror=None)
    print("############\n\n############")
# write genes withput TM domains
# with open ("NoTMGenes.txt", 'w') as f:
#     f.write('\n'.join(notTM_genes))
# write genes with AS exon that does not devide by 3
# with open ("NotDevideBy3Exons.txt", 'w') as f:
#     f.write('\n'.join(notDevideByThree))








