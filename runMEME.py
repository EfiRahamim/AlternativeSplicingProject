# must activate meme conda env: /home/alu/rahamie4/anaconda3/envs/meme
# 1. read splicing events result File
# 2. get sequence of every target exon +- 250 bases
# 3. run the sequence in meme
# 4. save results in a disired folder
import pandas as pd
import requests, sys, os, subprocess, multiprocessing, time

server="https://rest.ensembl.org"
results_file = "/private10/Projects/Efi/ArielBashari/PSI-Sigma_gencodeGTF/Results_withTPM/SplicingEventsFiltered-Scrambled_vs_4repeats-PSI20_Pvalue0.05_FDR0.05.csv"
output_dir = "/private10/Projects/Efi/ArielBashari/MEME/Scrambled_vs_decoy4/"
bases_to_add = 250

def getSequence(row):
  # define strand as 1/-1
  if row['Strand'] == '+':
    strand = '1'
  elif row['Strand'] == '-':
    strand = '-1'
  # get sequence by region
  chr = row['Target.Exon'].split(":")[0]
  start = str(int(row['Target.Exon'].split(":")[1].split("-")[0]) - bases_to_add)
  end = str(int(row['Target.Exon'].split(":")[1].split("-")[1]) + bases_to_add)
  regions = f"{chr}:{start}..{end}:{strand}"
  ext_region = f"/sequence/region/human/{regions}?"
  #print(f"Region seq request: {server+ext_region}")
  r_r = requests.get(server+ext_region, headers={ "Content-Type" : "text/x-fasta"}) 
  if not r_r.ok:
    print(f"Error in getting sequence by regions: {server+ext_region}. Exit")
    r_r.raise_for_status()
    sys.exit
    #return None
  seq = r_r.text.upper()
  return seq
def run_MEME(seq_file, target_dir):
  command = ["meme",
            seq_file,
            "-oc",
            target_dir,
            "-dna",
            "-mod",
            "anr",
            "-nmotifs",
            "3",
            "-nostatus"]
  try:
    output = subprocess.check_output(command, universal_newlines=True)
  except:
    print("Error in running MEME on subprocess. Exit.")
def run_event(index, row):
  seq = getSequence(row)
  gene_dir = os.path.join(output_dir, row['Gene.Symbol'])
  if not os.path.isdir(gene_dir):
    os.mkdir(gene_dir)
  transcript_dir = os.path.join(gene_dir, row['FixedTranscript'])
  if not os.path.isdir(transcript_dir):
    os.mkdir(transcript_dir)
  eventType_dir = os.path.join(transcript_dir, row['Event.Type'])
  if not os.path.isdir(eventType_dir):
    os.mkdir(eventType_dir)
  target_exon = (row['Target.Exon'].replace(":", "_")).replace("-", "_")
  target_dir = os.path.join(eventType_dir, target_exon)
  if not os.path.isdir(target_dir):
    os.mkdir(target_dir)
#  seq_file = os.path.join(target_dir, "sequence.fasta")
#   with open (seq_file, 'w') as f:
#     f.write(seq)
  seq_file = os.path.join(output_dir,"sequences.fasta")
  with open (seq_file, 'a') as f:
    f.write(seq)
  #run_MEME(seq_file, target_dir)

start_time = time.time()
results = pd.read_csv(results_file, index_col=False) # read PSI-Sigma results file
pool = multiprocessing.Pool(processes=5)
pool.starmap(run_event,results.iterrows())
pool.close()
pool.join()
end_time = time.time()
total_time = end_time-start_time
print(f"Execution time: ~{int(total_time)/60:.2f} minutes (~{int(total_time)/3600:.2f} hours).")


