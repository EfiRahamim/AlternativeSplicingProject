# must activate meme conda env: /home/alu/rahamie4/anaconda3/envs/meme
# 1. read splicing events result File
# 2. get sequence of every target exon +- 250 bases
# 3. run the sequence in meme
# 4. save results in a disired folder
import pandas as pd
import requests, sys, os, subprocess, multiprocessing, time

server="https://rest.ensembl.org"
results_file = "/private10/Projects/Efi/ArielBashari/PSI-Sigma_gencodeGTF/Results/Scrambled_vs_4repeats_analyzed.csv"
output_dir = "/private10/Projects/Efi/ArielBashari/Motifs/Scrambled_vs_decoy4/XSTREME/"
bases_to_add = 250
onlySES=True
splitIncExc=True

def getSequence(row, bases_to_add):
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
  r_r = requests.get(server+ext_region, headers={ "Content-Type" : "text/plain"}) 
  if not r_r.ok:
    print(f"Error in getting sequence by regions: {server+ext_region}. Exit")
    r_r.raise_for_status()
    sys.exit
    #return None
  seq = r_r.text.upper()
  #full_seq = f">{row['Event.Region']}_{row['Gene.Symbol']}_{row['Event.Type']}_dPSI:{row['dPSI']}\n"+seq+'\n'
  return seq
def run_MEME(seq_file, target_dir):
  # command for conda
  conda_command = ["meme",
            seq_file,
            "-oc",
            target_dir,
            "-dna",
            "-mod",
            "anr",
            "-nmotifs",
            "3",
            "-nostatus"]
  # command for docker
  docker_command = ["docker", "run" ,"--rm", 
             "-v" ,"/private10/Projects/Efi/:/private10/Projects/Efi/:rw", 
             "-u" ,"$(id -u ${USER}):$(id -g ${USER})", 
             "memesuite/memesuite:latest", "xstreme", 
             "--p", seq_file, 
             "--oc", target_dir,
             "-m", "/opt/meme/share/meme-5.5.5/db/motif_databases/RNA/Ray2013_rbp_Homo_sapiens.meme", 
             "--dna2rna", "--minw", "6", "--maxw", "15"]
  try:
    output = subprocess.check_output(docker_command, universal_newlines=True)
  except:
    print("Error in running MEME on subprocess. Exit.")
def run_event(index, row):
  if onlySES and row['Event.Type'] != 'SES':
      return
  global target_exons_list
  if row['Target.Exon'] in target_exons_list:
    return
  else:
    target_exons_list.append(row['Target.Exon'])
  seq = getSequence(row, bases_to_add=250)
  #full_seq = f">{row['Event.Region']}_{row['Gene.Symbol']}_{row['Event.Type']}_dPSI:{row['dPSI']}\n"+seq+'\n'
  seq_record = f">{row['Event.Region']}_{row['Gene.Symbol']}_{row['Event.Type']}_dPSI:{row['dPSI']}\n"

  if splitIncExc and row['dPSI'] > 0: #Inclusion form
    inc_seq_file = os.path.join(output_dir,"sequences_Inc.fasta")
    full_seq = seq_record+seq+'\n'
    with open (inc_seq_file, 'a') as f:
      f.write(full_seq)
      return
  elif splitIncExc and row['dPSI'] < 0: #Exclusion form
    exon = getSequence(row, bases_to_add=0)
    exc_seq = seq.replace(exon, "")
    exc_seq_file = os.path.join(output_dir,"sequences_Exc.fasta")
    full_seq = seq_record+exc_seq+'\n'
    with open (exc_seq_file, 'a') as f:
      f.write(full_seq)
      return

  # gene_dir = os.path.join(output_dir, row['Gene.Symbol'])
  # if not os.path.isdir(gene_dir):
  #   os.mkdir(gene_dir)
  # transcript_dir = os.path.join(gene_dir, row['Reference.Transcript'])
  # if not os.path.isdir(transcript_dir):
  #   os.mkdir(transcript_dir)
  # eventType_dir = os.path.join(transcript_dir, row['Event.Type'])
  # if not os.path.isdir(eventType_dir):
  #   os.mkdir(eventType_dir)
  # target_exon = (row['Target.Exon'].replace(":", "_")).replace("-", "_")
  # target_dir = os.path.join(eventType_dir, target_exon)
  # if not os.path.isdir(target_dir):
  #   os.mkdir(target_dir)
  # seq_file = os.path.join(target_dir, "sequence.fasta")
  # with open (seq_file, 'w') as f:
  #   f.write(seq)
  
  seq_file = os.path.join(output_dir,"sequences_all.fasta")
  full_seq = seq_record+seq+'\n'
  with open (seq_file, 'a') as f:
    f.write(full_seq)
  # run_MEME(seq_file, target_dir)

start_time = time.time()
target_exons_list = []
results = pd.read_csv(results_file, index_col=False) # read PSI-Sigma results file
pool = multiprocessing.Pool(processes=5)
pool.starmap(run_event,results.iterrows())
pool.close()
pool.join()
end_time = time.time()
total_time = end_time-start_time
print(f"Execution time: ~{int(total_time)/60:.2f} minutes (~{int(total_time)/3600:.2f} hours).")


