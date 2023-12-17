# 1. Read the results table
# 2. for each transcript - find the mRNA and CDS of the sequence according to the GTF file
# 4. look for the target exon within the sequence
# 5. Count the exons that has been found and the exons that hasent.
# 6. Create inclusion and exclusion amino acid sequences for cases of alternative exo within the CDS

import csv, requests, multiprocessing, sys, os, subprocess, re, time, argparse
import pandas as pd
from Bio import SeqIO, Seq
from Bio.Seq import Seq
from io import StringIO

# CLI arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Analyzing transcripts of PSI-Sigma results: Creates the inclusion and exclusion amino acids sequences for each relevant transcript. Results will be saved in the giving output directory.\n")
parser.add_argument("-i", action='store', dest='input_file', required=True, help="Input file: PSI-Sigma filtered csv file")
parser.add_argument("-o", action='store', dest='output_dir', required=True, help="Output directory for the splicing events sequences")
parser.add_argument("-gtf_file", action='store', dest='gtf_file', help="GTF file that has been used in PSI-Sigma analysis.", default="/private10/Projects/Efi/General/gencode.v28.annotation.gtf")
parser.add_argument("-gtf_map", action='store', dest='gtf_map', required=True, help="GTF map file from PSI-Sigma (located in PSI-Sigma results directory).")
parser.add_argument("-l1", action='store', dest='group_A', required=True, help="Lable of group 1 (first group in PSI-sigma analysis)")
parser.add_argument("-l2", action='store', dest='group_B', required=True, help="Lable of group 2 (second group in PSI-sigma analysis)")
user_args = parser.parse_args()

input_file = user_args.input_file
output_dir = user_args.output_dir
gtf_file = user_args.gtf_file
gtf_map = user_args.gtf_map
group_A = user_args.group_A
group_B = user_args.group_B

# DEBUG Arguments
# input_file = "/private10/Projects/Efi/CRG/GBM/PSI-Sigma/GencodeGTF/TM_Genes/SplicingEventsFiltered-DMSO_vs_H3B8800-PSI20_Pvalue0.05_FDR0.05.csv"
# group_A = 'DMSO'
# group_B = 'H3B880'
# output_dir = "/private10/Projects/Efi/CRG/GBM/PSI-Sigma/GencodeGTF/TM_Genes/"
# gtf_map = "/private10/Projects/Efi/CRG/GBM/PSI-Sigma/GencodeGTF/DMSO_vs_H3B8800/gencode.v28.annotation.gtf.mapping.txt"
# gtf_file = "/private10/Projects/Efi/General/gencode.v28.annotation.gtf"

def getExonByRow(row):
  # define strand as 1/-1
  if row['Strand'] == '+':
    strand = '1'
  elif row['Strand'] == '-':
    strand = '-1'
  # get sequence by region
  regions = row['Target.Exon'].replace("-","..")
  regions = f"{regions}:{strand}"
  ext_region = f"/sequence/region/human/{regions}?"
  #print(f"Region seq request: {server+ext_region}")
  r_r = requests.get(server+ext_region, headers={ "Content-Type" : "text/plain"}) 
  if not r_r.ok:
    print(f"Error in getting exon by regions: {server+ext_region}. Skipping.")
    r_r.raise_for_status()
    sys.exit
    #return None
  exon = r_r.text.upper()
  return exon
def getExonByRegions(chr, start, end, strand_sign):
  # define strand as 1/-1
  if strand_sign == '+':
    strand = '1'
  elif strand_sign == '-':
    strand = '-1'
  # get sequence by region
  regions = f"{chr}:{start}..{end}:{strand}"
  ext_region = f"/sequence/region/human/{regions}?"
  #print(f"Region seq request: {server+ext_region}")
  r_r = requests.get(server+ext_region, headers={ "Content-Type" : "text/plain"}) 
  if not r_r.ok:
    print(f"Error in getting exon by regions: {server+ext_region}. Skipping.")
    r_r.raise_for_status()
    sys.exit
    # print(f"Error in getting exon by regions. Skipping.")
    # return None
  exon = r_r.text.upper()
  return exon
def parse_gtf(filename):
  start_time = time.time()
  print("Parsing GTF file..")
  columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'] # columns of interest
  data = []
  with open(filename, 'r') as file:
      for line in file:
          if line.startswith('#'):
              continue  # Skip comment lines
          parts = line.strip().split('\t') # split GTF fields
          parts[8]=parts[8].replace('"', "") 
          attributes = dict(item.strip().split(' ') for item in parts[8].split(';') if item.strip()) # split 'attributes' field
          parts[8] = attributes  # Replace 'attribute' field with a dictionary
          data.append(parts)
  df = pd.DataFrame(data, columns=columns)
  end_time = time.time()
  total_time = end_time - start_time
  print(f"Done parsing GTF file (~{int(total_time)} sec).")
  return df
def getSeqFromGTF (transcript_id, chr):
  filtered_gtf_df = gtf_df[(gtf_df['seqname'] == chr) & (gtf_df['feature'] == 'exon')] # filter GTF data frame
  transcript_found = False # flag
  transcript_seq = ""
  for index, row in filtered_gtf_df.iterrows(): # run over the filtered GTF data frame
    if row['attribute']['transcript_id'] == transcript_id and not transcript_found: # trascript was found
      transcript_found = True
      print(f"Transcript {transcript_id} was found in the GTF file.")
    if row['attribute']['transcript_id'] == transcript_id and transcript_found: # build sequence
      chr, start, end, strand = row['seqname'], row['start'], row['end'], row['strand']
      #print(f"Exon: {start} - {end}")
      exon = getExonByRegions(chr, start, end, strand)
      if strand == '+':
        transcript_seq = transcript_seq+exon
      elif strand == '-':
        transcript_seq = exon+transcript_seq
      continue
    if row['attribute']['transcript_id'] != transcript_id and transcript_found: # end of transcript
      break
  print(f"Transcript {transcript_id}: {len(transcript_seq)} nucleotides.")
  return transcript_seq
def getCDSFromGTF (transcript_id, chr):
  cds_found = False
  cds_seq = ""
  filtered_gtf_df = gtf_df[(gtf_df['seqname'] == chr) &
                           ((gtf_df['feature'] == 'CDS') | (gtf_df['feature'] == 'stop_codon'))] # filter GTF data frame
  if filtered_gtf_df.shape[0] == 0: # check if it's a coding transcript
    print(f"No CDS for transcript {transcript_id}. Skipping.")
    return cds_seq
  for index, row in filtered_gtf_df.iterrows(): # run over the filtered GTF data frame
    if row['feature'] == 'CDS' and row['attribute']['transcript_id'] == transcript_id and not cds_found: # CDS was found
      cds_found = True
      print(f"CDS of {transcript_id} was found in the GTF file.")
    if row['feature'] == 'CDS' and cds_found: # build the CDS sequence
      chr, start, end, strand = row['seqname'], row['start'], row['end'], row['strand']
      #print(f"CDS: {start} - {end}")
      exon = getExonByRegions(chr, start, end, strand)
      if strand == '+':
        cds_seq = cds_seq+exon
      elif strand == '-':
        cds_seq = exon+cds_seq
      continue
    if row['feature'] == 'stop_codon' and cds_found: # reached the stop codon of the CDS
      chr, start, end, strand = row['seqname'], row['start'], row['end'], row['strand']
      #print(f"Stop codon: {start} - {end}")
      exon = getExonByRegions(chr, start, end, strand)
      if strand == '+':
        cds_seq = cds_seq+exon
      elif strand == '-':
        cds_seq = exon+cds_seq
      break # end of CDS
  print(f"CDS of {transcript_id}: {len(cds_seq)} nucleotides.")
  return cds_seq
def getAminoAcidSeq(inclusion_seq_object, exclusion_seq_object, transcript_id):
  # create AA sequences
  try:
      AA_inclusion_seq = inclusion_seq_object.translate(cds=True) 
  except:
    print(f"Error in translating inclusion seq of transcript {transcript_id}. Skipping.")
    return None, None
  try:
      AA_exclusion_seq = exclusion_seq_object.translate(to_stop=True) # translation is terminated at the first in frame stop codon
  except:
    print(f"Error in translating exclusion seq of transcript {transcript_id}. Skipping.")
    return None, None
  return AA_inclusion_seq, AA_exclusion_seq 
def createDir (output_dir, gene_name, transcript_id, splicing_event, index):
  # modify splicing event name to be used in file system paths
  splicing_event = splicing_event.replace("|", "_") # for 'TSS|' cases
  if splicing_event == 'IR (overlapping region)':
    splicing_event = 'IR_OLR'
  # super-directory of files
  files_dir = os.path.join(output_dir,'SplicingEventsFiles')
  if not os.path.isdir(files_dir):
    os.mkdir(files_dir)
  # directory of gene
  gene_dir = os.path.join(files_dir,gene_name)
  if not os.path.isdir(gene_dir):
      os.mkdir(gene_dir)
  # directory of transcript
  transcript_dir = os.path.join(gene_dir,transcript_id)
  if not os.path.isdir(transcript_dir):
    os.mkdir(transcript_dir)
  # directory of splicing event
  splicing_event_dir = os.path.join(transcript_dir,splicing_event)
  if not os.path.isdir(splicing_event_dir):
    os.mkdir(splicing_event_dir)
  # directory of event ID
  eventID_dir = os.path.join(splicing_event_dir,str(index))
  if not os.path.isdir(eventID_dir):
    os.mkdir(eventID_dir)
  return eventID_dir
def SaveSeqAsFastaFiles(AA_inclusion_seq, AA_exclusion_seq, files_dir, row):
  # define groups to isoforms
  if row['dPSI'] > 0 :
    inclusion_group = group_B
    exclusion_group = group_A
  else:
    inclusion_group = group_A
    exclusion_group = group_B
  # files pathes
  inclusion_path = os.path.join(files_dir, 'Inclusion_'+inclusion_group+'.fasta')
  exclusion_path = os.path.join(files_dir, 'Exclusion_'+exclusion_group+'.fasta')
  # files content
  inclusion_record = ">"+row['Gene.Symbol']+"_"+row['Reference.Transcript']+'_'+row['Event.Type']+'_Inclusion_'+inclusion_group+'\n'+str(AA_inclusion_seq)
  exclusion_record = ">"+row['Gene.Symbol']+"_"+row['Reference.Transcript']+'_'+row['Event.Type']+'_Exclusion_'+exclusion_group+'\n'+str(AA_exclusion_seq)
  # write files
  with open(inclusion_path, 'w') as f:
    f.write(inclusion_record)
  with open(exclusion_path, 'w') as f:
    f.write(exclusion_record)
def run_with_GTF(index,row):
  global transcripts_dict, cds_dict
  # remove non-required prefixes
  transcript_id = row['Reference.Transcript'].replace("Ex.", "")
  transcript_id = transcript_id.replace("TSS.", "")
  chr = row['Event.Region'].split(":")[0]
  # get transcript seq
  if transcript_id in transcripts_dict.keys(): # use prepared sequence
    transcript_seq = transcripts_dict[transcript_id]
    print(f"Trascript {transcript_id} was already analyzed. Using the prepared sequence.")  
  else : # build the sequence
    transcript_seq = getSeqFromGTF(transcript_id, chr)
    transcripts_dict[transcript_id] = transcript_seq
  if transcript_seq is None:
    row['Transcript found?'] = 'no'
    return row
  else :
    row['Transcript found?'] = 'yes'
  # check if exon in the transcript (mRNA)
  exon_seq = getExonByRow(row)
  if exon_seq in transcript_seq:
    row['Exon in transcript?'] = 'yes'
  else:
    row['Exon in transcript?'] = 'no'
    return row
  # check if exon in the CDS - using the GTF file
  if transcript_id in cds_dict.keys(): # use prepared CDS sequence
    cds_seq = cds_dict[transcript_id]
    print(f"CDS of trascript {transcript_id} was already analyzed. Using the prepared sequence.")
  else: # build the CDS sequence
    cds_seq = getCDSFromGTF(transcript_id, chr)
    cds_dict[transcript_id] = cds_seq
  if len(cds_seq) == 0: # Non-coding transcript
    row['Exon in CDS?'] = 'Non-coding transcript'
    return row
  if exon_seq in cds_seq:
    row['Exon in CDS?'] = 'yes'
    if (len(exon_seq) % 3 == 0): # check if exon is devided by 3
      row['Exon devide by 3?'] = 'yes'
    else:
      row['Exon devide by 3?'] = 'no'
    # get inclusion and exclusion sequences
    inclusion_seq = cds_seq
    exclusion_seq = cds_seq.replace(exon_seq, "") # remove alternative exon
    # translate RNA sequences into protein
    AA_inclusion_seq, AA_exclusion_seq = getAminoAcidSeq(Seq(inclusion_seq), Seq(exclusion_seq), transcript_id)
    if AA_inclusion_seq == None or AA_exclusion_seq == None:
        return row
    else:
      files_dir = createDir(output_dir, row['Gene.Symbol'], transcript_id, row['Event.Type'], index)
      SaveSeqAsFastaFiles(AA_inclusion_seq, AA_exclusion_seq, files_dir, row)
      row['InclusionAAseq'] = str(AA_inclusion_seq)
      row['ExclusionAAseq'] = str(AA_exclusion_seq)
  else:
    row['Exon in CDS?'] = 'no'
  return row

start_time = time.time()
server="https://rest.ensembl.org"

# Pre-processing of the data
results = pd.read_csv(input_file, index_col=False) # read PSI-Sigma results file
gene_map = pd.read_csv(gtf_map, delimiter='\t', index_col=False) # read GTF map file
gene_map.columns = ['Transcript', 'Gene.Symbol', 'Strand']
gene_map.drop('Transcript', axis=1, inplace=True)
gene_map.drop_duplicates(inplace=True)
merged_results = pd.merge(results, gene_map, on='Gene.Symbol', how='inner') # add strand to results
merged_results = merged_results.rename(columns={merged_results.columns[8]: 'dPSI'}) # rename 'delta PSI' column 
# add columns of exon location
merged_results['Transcript found?'] = 'NA'
merged_results['Exon in transcript?'] = 'NA'
merged_results['Exon in CDS?'] = 'NA'
merged_results['Exon devide by 3?'] = 'NA'
merged_results['InclusionAAseq'] = 'NA'
merged_results['ExclusionAAseq'] = 'NA'

# global dictionaries
transcripts_dict = {}
cds_dict = {}
gtf_df = parse_gtf(gtf_file) # parse the GTF file into shared data frame
# run analyze in parallel
pool = multiprocessing.Pool(processes=5) 
updated_rows = pool.starmap(run_with_GTF,merged_results.iterrows())
pool.close()
pool.join()
# write updated results file
merged_results_updated = pd.DataFrame(updated_rows)
updated_results_file = os.path.join(output_dir, group_A+'_vs_'+group_B+'_analyzed.csv')
merged_results_updated.to_csv(updated_results_file, index=True)
# print splicing events locations to the screen
print(f"Found transcripts: {(merged_results_updated['Transcript found?']=='yes').sum()} out of {merged_results_updated.shape[0]}")
print(f"Exon in transcripts: {(merged_results_updated['Exon in transcript?']=='yes').sum()} out of {(merged_results_updated['Transcript found?']=='yes').sum()}")
print(f"Exon in CDS: {(merged_results_updated['Exon in CDS?']=='yes').sum()} out of {(merged_results_updated['Exon in transcript?']=='yes').sum()}")
print(f"Exon devide by 3: {(merged_results_updated['Exon devide by 3?']=='yes').sum()} out of {(merged_results_updated['Exon in CDS?']=='yes').sum()}")
end_time = time.time()
total_time = end_time-start_time
print(f"Execution time: ~{int(total_time)/60:.2f} minutes (~{int(total_time)/3600:.2f} hours).")

