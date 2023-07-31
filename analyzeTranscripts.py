from Bio import SeqIO, Seq
import os, re
import requests, sys, csv, argparse, subprocess, multiprocessing
from io import StringIO
from Bio.Seq import Seq

# CLI arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Analyzing transcripts of rMATS results: Creating full and spliced nucleic and amino acids sequences for each transcript. Results will be saved in the giving output directory.\nNOTE: NOT SUITABLE FOR MXE EVENTS (YET)!!!")
parser.add_argument("-i", action='store', dest='input', required=True, help="Input file: rMATS filtered csv file")
parser.add_argument("-type", action='store', dest='as_type', required=True, help="Type of Splicing Event - SE/A5SS/A3SS/MXE/RI.", choices=["SE","A3SS","A5SS","MXE","RI"])
# parser.add_argument("-a5ss", action='store', dest='A5SS_output', required=False, help="Input file: A5SS filtered csv file")
# parser.add_argument("-a3ss", action='store', dest='A3SS_output', required=False, help="Input file: A3SS filtered csv file")
# parser.add_argument("-mxe", action='store', dest='MXE_output', required=False, help="Input file: MXE filtered csv file")
# parser.add_argument("-ri", action='store', dest='RI_output', required=False, help="Input file: RI filtered csv file")
parser.add_argument('-o', action='store', dest='output_dir', default=os.getcwd(), help="Output directory")
parser.add_argument('-l1', action='store', dest='label1', help="Lable of first group in the analyze (e.g. Normal)")
parser.add_argument('-l2', action='store', dest='label2', help="Lable of second group in the analyze (e.g. Tumor)")
parser.add_argument('--mulby3', action='store_true', dest='mulby3', help='If set - filterring out spliced transcripts that were not multipule of 3.')
user_args = parser.parse_args()

def defineType(row):
    # define type
  if float(row['IncLevelDifference']) > 0:
    #return ['Normal', 'Tumor']
    return [user_args.label1, user_args.label2]
  elif float(row['IncLevelDifference']) < 0:
    #return ['Tumor', 'Normal']
    return [user_args.label2, user_args.label1]
  else:
    print(f"Gene {row['GeneID']} has 'IncLevelDifference' equal to zero. Skipping.")
    with open ("NoType.txt", 'a') as f:
      f.write(row['GeneID'] + '\n')
    return None

def getExon(row):
  # define strand as 1/-1
  if row['strand'] == '+':
    strand = '1'
  elif row['strand'] == '-':
    strand = '-1'
  # get sequence by region
  regions = f"{row['chr']}:{row['exonStart_1base']}..{row['exonEnd']}:{strand}"
  ext_region = f"/sequence/region/human/{regions}?"
  print(f"Region seq request: {ext_region}")
  r_r = requests.get(server+ext_region, headers={ "Content-Type" : "text/plain"}) 
  if not r_r.ok:
    #r_r.raise_for_status()
    #sys.exit
    print(f"Error in getting exon by regions. Skipping.")
    with open ("NoExon.txt", 'a') as f:
      f.write(row['GeneID'] + '\n')
    return None
  exon = r_r.text.upper()
  return exon

def getSeq(row, transcript, as_type):
  # get sequence by trascript ID
  # transcripts = row['optional transcripts'].split(";")
  # found = False
  # for transcript in transcripts:
  #   #transcript = transcript.replace(".","/.")
  ext_transcript = f"/sequence/id/{transcript.replace('.','/.')}?type=cds"
  print(f"transcript sequence request: {server+ext_transcript}")
  r_t = requests.get(server+ext_transcript, headers={ "Content-Type" : "text/x-fasta"})
  if not r_t.ok:
    #r_t.raise_for_status()
    #sys.exit
    print(f"Error in: {server+ext_transcript}")
    return None, None
  # edit the "\n" in the sequence - keep only first occurance of "\n"
  seq = r_t.text.replace("\n", "")
  match_start_seq = re.search(r'\.\d+([a-zA-Z])', seq) # find the end index of the transcript name - ends with a dot and 1 or 2 numbers after
  seq_start_index = match_start_seq.start(1)
  seq = seq[:seq_start_index]+"\n"+seq[seq_start_index:]
  #print("seq before inclusion: ", seq)
  # handle case of RI: full seq will include the retained intron
  if as_type == 'RI':
    # check that upstream/downstream exons appears only once in the seq
    if seq.count(row['upstreamExonSeq']) != 1:
      print("Warning in RI event: upstreamExonSeq is not found or has more than one occurrence. Trying by downstreamExonSeq..")
      if seq.count(row['downstreamExonSeq']) != 1:
        print("downstreamExonSeq is not found or has more than one occurrence. Skipping to next transcript.")
        return None, None
      else: 
        # create the full seq by adding the retained intron before the downstreamExon
        seq = seq.replace(row['downstreamExonSeq'], row['RetainedIntronSeq']+row['downstreamExonSeq'])
        print("Creating the inclusion sequence according to the downstreamExon")    
    else:
      # create the full seq by adding the retainded intron after the upstreamExon
      seq = seq.replace(row['upstreamExonSeq'], row['upstreamExonSeq']+row['RetainedIntronSeq'])
      print("Creating the inclusion sequence according to the uptreamExon")
    #print("seq after inclusion: ", seq)
    #print("length of retained intron: ", len(row['RetainedIntronSeq']))
  # make FASTA record for the seuence
  fasta_str = StringIO(seq)
  # Parse the contents of the StringIO object as a FASTA file
  records = list(SeqIO.parse(fasta_str, "fasta"))
  # if exon in records[0].seq:
  #   print(f"Gene: {row['GeneID']} Matched transcript: {transcript}")
  #   found = True
  #   #break
  return seq, records[0]
  # else:
  #   continue
  # if not found:
  #   print(f"No matching transcript were found for {row['GeneID']}")
  #   with open ("NoSeq.txt", 'a') as f:
  #     f.write(row['GeneID'] + '\n')
  #   return None, None, None

def spliceSeq(fullSeqRecord, exon):
  spliced_seq = str(fullSeqRecord.seq).replace(exon, "")
  if len(spliced_seq) == len(exon):
    return None
  return spliced_seq
  
def find_spliced_peptide(full, spliced):
  # find the first change, from the begining of the sequence
  i = 0
  j = 0
  while i < len(full) and j < len(spliced):
    if full[i] == spliced[j]:
      i = i+1
      j = j+1
    else:
      break
  start = i
  # find the first change, from the end of the sequence
  i = len(full)-1
  j=len(spliced)-1
  while i >= 0 and j >= 0:
    if full[i] == spliced[j]:
      i = i-1
      j = j-1
    else:
      break
  end = i
  #print(f"{start+1}:{end+1}\n{full[start:end]}")
  return full[start:end], start+1, end+1

def getAminoAcidSeq(full_seq_record, spliced_seq):
  # create AA sequences
  try:
    if (not user_args.mulby3) or (user_args.as_type == 'RI'):
      AA_full_seq = full_seq_record.seq.translate(stop_symbol="") # retainded intron can affects the ORF, causing short isoform
    else:
      AA_full_seq = full_seq_record.seq.translate(cds=True) 
  except:
    print("Error in translating full RNA seq. Skipping")
    return None, None
  AA_spliced_seq = Seq(spliced_seq)
  try:
    if not user_args.mulby3: # if spliced seq is not a multiple of 3 - short peptide can be formes
      AA_spliced_seq = AA_spliced_seq.translate(stop_symbol="")  
    else:
      AA_spliced_seq = AA_spliced_seq.translate(cds=True)
  except:
    print("Error in translating spliced RNA seq. Skipping.")
    return None, None
  return AA_full_seq, AA_spliced_seq 

def run_one_row(row):
  # check if current gene has matching transcript(s)
  transcripts = row['optional transcripts'].split(';')
  if len(transcripts) == 0 or transcripts[0] == '':
    print(f"Gene {row['GeneID']} has no matching transcripts. Skipping to next gene.\n")
    #continue
    return
  else:
    # define type of full/spliced sequences
    type = defineType(row)
    if type is None:
      #continue
      return
    # get AS exon sequence
    if user_args.as_type in ['SE', 'A5SS', 'A3SS']:
      exon = row['AlternativeExonSeq']
    elif user_args.as_type == 'RI':
      exon = row['RetainedIntronSeq']
    elif user_args.as_type == 'MXE':
      exon = row['AlternativeExonSeq1']
      exon_2 = row['AlternativeExonSeq2']
    if exon is None:
      print(f"Error in retreving exon sequence from the input file. Skipping to next gene.")
      #continue
      return
    # run over each transcript
    for transcript in transcripts:
      print(f"Entering row ID {row['ID']}, GeneID {row['GeneID']}, Transcript {transcript}")    
      # get full sequence of current transcript
      NA_seq, full_seq_record = getSeq(row, transcript, user_args.as_type)
      if full_seq_record is None:
        continue
      # create the spliced sequence
      spliced_seq = spliceSeq(full_seq_record, exon)
      if spliced_seq is None:
        print(f"Error in splicing transcript {transcript} of gene {row['GeneID']}. Skipping")
        with open('noSplicing.txt', 'a') as f:
          f.write(transcript + '\n')
        continue
      # check if spliced sequence is a multiple of 3 (only if flag '--mulby3' passed by the user)
      if (user_args.mulby3) and (len(spliced_seq) % 3 != 0):
        print(f"Length of spliced seq is {len(spliced_seq)}. not multiple of 3. Skipping.")
        with open ("noMulByThree.txt", 'a') as f:
          f.write(transcript + '\n')
        continue
      # get Amino Acid sequences
      AA_full_seq, AA_spliced_seq = getAminoAcidSeq(full_seq_record, spliced_seq)
      if AA_full_seq == None or AA_spliced_seq == None:
        continue
      # create directories
      gene_dir = os.path.join(user_args.output_dir,row['GeneID']+'_'+row['geneSymbol'])
      if not os.path.isdir(gene_dir):
        os.mkdir(gene_dir)
      id_dir = os.path.join(gene_dir, row['ID']) 
      if not os.path.isdir(id_dir):
        os.mkdir(id_dir)
      transcript_dir = os.path.join(id_dir, transcript)
      # check if current transcript was already checked in previous running
      if os.path.isdir(transcript_dir):
        print(f"Transcript {transcript} of gene {row['GeneID']} was already checked.\nResults can be found in {transcript_dir}.\nSkipping to next transcript.")
        continue
      # save NA sequences
      os.mkdir(transcript_dir)
      na_fasta_path = os.path.join(transcript_dir,transcript+"_fullNA_"+type[0]+".fasta")
      with open(na_fasta_path, 'w') as na_fasta_file:
        na_fasta_file.write(NA_seq)
      spliced_seq_record = ">"+row['GeneID']+"_"+transcript+"_"+type[1]+'\n'+str(spliced_seq)+'\n'
      spliced_na_fasta_path = os.path.join(transcript_dir,transcript+"_splicedNA_"+type[1]+".fasta")
      with open(spliced_na_fasta_path, 'w') as spliced_na_fasta_file:
        spliced_na_fasta_file.write(spliced_seq_record)
      # save AA sequences
      AA_full_seq_record = ">"+row['GeneID']+"_"+transcript+"_"+type[0]+'\n'+str(AA_full_seq)+'\n'
      AA_spliced_seq_record = ">"+row['GeneID']+"_"+transcript+"_"+type[1]+'\n'+str(AA_spliced_seq)+'\n'
      AA_full_seq_path = os.path.join(transcript_dir,transcript+'_fullAA_'+type[0]+'.fasta')
      AA_spliced_seq_path = os.path.join(transcript_dir,transcript+'_splicedAA_'+type[1]+'.fasta')
      with open(AA_full_seq_path, 'w') as f:
        f.write(AA_full_seq_record)
      with open(AA_spliced_seq_path, 'w') as f:
        f.write(AA_spliced_seq_record)
      # save AA spliced peptide
      spliced_peptide, start,end = find_spliced_peptide(str(AA_full_seq),str(AA_spliced_seq))
      spliced_peptide_record = "Spliced peptide coordinates: "+str(start)+"-"+str(end)+"\n"+str(spliced_peptide)
      spliced_peptide_path=os.path.join(transcript_dir,"spliced_peptide.txt")
      with open(spliced_peptide_path, 'w') as f:
        f.write(spliced_peptide_record)


def run_analyse(input_file):
  # run over the input file (e.g. filtered rMATS output file)
  with open(input_file, 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    rows = list(reader)
  pool = multiprocessing.Pool(processes=10) 
  pool.map(run_one_row,rows)
  pool.close()
  pool.join()
    # # run analys for each row (e.g. each AS event)
    # for row in rows:
    #   # check if current gene has matching transcript(s)
    #   transcripts = row['optional transcripts'].split(';')
    #   if len(transcripts) == 0 or transcripts[0] == '':
    #     print(f"Gene {row['GeneID']} has no matching transcripts. Skipping to next gene.\n")
    #     continue
    #   else:
    #     # define type of full/spliced sequences
    #     type = defineType(row)
    #     if type is None:
    #       continue
    #     # get AS exon sequence
    #     if user_args.as_type in ['SE', 'A5SS', 'A3SS']:
    #       exon = row['AlternativeExonSeq']
    #     elif user_args.as_type == 'RI':
    #       exon = row['RetainedIntronSeq']
    #     elif user_args.as_type == 'MXE':
    #       exon = row['AlternativeExonSeq1']
    #       exon_2 = row['AlternativeExonSeq2']
    #     if exon is None:
    #       print(f"Error in retreving exon sequence from the input file. Skipping to next gene.")
    #       continue
    #     # run over each transcript
    #     for transcript in transcripts:
    #       print(f"Entering row ID {row['ID']}, GeneID {row['GeneID']}, Transcript {transcript}")    
    #       # get full sequence of current transcript
    #       NA_seq, full_seq_record = getSeq(row, transcript, user_args.as_type)
    #       if full_seq_record is None:
    #         continue
    #       # create the spliced sequence
    #       spliced_seq = spliceSeq(full_seq_record, exon)
    #       if spliced_seq is None:
    #         print(f"Error in splicing transcript {transcript} of gene {row['GeneID']}. Skipping")
    #         with open('noSplicing.txt', 'a') as f:
    #           f.write(transcript + '\n')
    #         continue
    #       if len(spliced_seq) % 3 != 0:
    #         print(f"Length of spliced seq is {len(spliced_seq)}. not multiple of 3. Skipping.")
    #         with open ("noMulByThree.txt", 'a') as f:
    #           f.write(transcript + '\n')
    #         continue
    #       # get Amino Acid sequences
    #       AA_full_seq, AA_spliced_seq = getAminoAcidSeq(full_seq_record, spliced_seq)
    #       if AA_full_seq == None or AA_spliced_seq == None:
    #         continue
    #       # create directories
    #       gene_dir = os.path.join(user_args.output_dir,row['GeneID']+'_'+row['geneSymbol'])
    #       if not os.path.isdir(gene_dir):
    #         os.mkdir(gene_dir)
    #       id_dir = os.path.join(gene_dir, row['ID']) 
    #       if not os.path.isdir(id_dir):
    #         os.mkdir(id_dir)
    #       transcript_dir = os.path.join(id_dir, transcript)
    #       # check if current transcript was already checked in previous running
    #       if os.path.isdir(transcript_dir):
    #         print(f"Transcript {transcript} of gene {row['GeneID']} was already checked.\nResults can be found in {transcript_dir}.\nSkipping to next transcript.")
    #         continue
    #       # save NA sequences
    #       os.mkdir(transcript_dir)
    #       na_fasta_path = os.path.join(transcript_dir,transcript+"_fullNA_"+type[0]+".fasta")
    #       with open(na_fasta_path, 'w') as na_fasta_file:
    #         na_fasta_file.write(NA_seq)
    #       spliced_seq_record = ">"+row['GeneID']+"_"+transcript+"_"+type[1]+'\n'+str(spliced_seq)+'\n'
    #       spliced_na_fasta_path = os.path.join(transcript_dir,transcript+"_splicedNA_"+type[1]+".fasta")
    #       with open(spliced_na_fasta_path, 'w') as spliced_na_fasta_file:
    #         spliced_na_fasta_file.write(spliced_seq_record)
    #       # save AA sequences
    #       AA_full_seq_record = ">"+row['GeneID']+"_"+transcript+"_"+type[0]+'\n'+str(AA_full_seq)+'\n'
    #       AA_spliced_seq_record = ">"+row['GeneID']+"_"+transcript+"_"+type[1]+'\n'+str(AA_spliced_seq)+'\n'
    #       AA_full_seq_path = os.path.join(transcript_dir,transcript+'_fullAA_'+type[0]+'.fasta')
    #       AA_spliced_seq_path = os.path.join(transcript_dir,transcript+'_splicedAA_'+type[1]+'.fasta')
    #       with open(AA_full_seq_path, 'w') as f:
    #         f.write(AA_full_seq_record)
    #       with open(AA_spliced_seq_path, 'w') as f:
    #         f.write(AA_spliced_seq_record)
    #       # save AA spliced peptide
    #       spliced_peptide, start,end = find_spliced_peptide(str(AA_full_seq),str(AA_spliced_seq))
    #       spliced_peptide_record = "Spliced peptide coordinates: "+str(start)+"-"+str(end)+"\n"+str(spliced_peptide)
    #       spliced_peptide_path=os.path.join(transcript_dir,"spliced_peptide.txt")
    #       with open(spliced_peptide_path, 'w') as f:
    #         f.write(spliced_peptide_record)

# server URL prefix for REST API requests
server = "https://rest.ensembl.org"
type=[]
noType=[]
noRegion =[]
noSeq=[]

print(f"Start analyze transcripts for {user_args.as_type}.")
os.chdir(user_args.output_dir)
run_analyse(user_args.input)
print(f"Done analyze transcripts for {user_args.as_type}.")







# # run over SE output file
# with open(user_args.SE_output, 'r') as csvfile:
#   reader = csv.DictReader(csvfile)
#   rows = list(reader)
#   # run over each row - each AS event
#   for row in rows:
#     # check if current gene has matching transcript(s)
#     transcripts = row['optional transcripts'].split(';')
#     if len(transcripts) == 0 or transcripts[0] == '':
#       print(f"Gene {row['GeneID']} has no matching transcripts. Skipping to next gene.\n")
#       continue
#     else:
#       # define type of full/spliced sequences - Normal or Tumor
#       type = defineType(row)
#       if type is None:
#         continue
#       # get AS exon sequence
#       exon = getExon(row)
#       if exon is None:
#         continue    
#       # run over each transcript
#       for transcript in transcripts:
#         print(f"Entering row ID {row['ID']}, GeneID {row['GeneID']}, Transcript {transcript}")    
#         # get full sequence of current transcript
#         NA_seq, full_seq_record = getSeq(row, transcript, exon)
#         if full_seq_record is None:
#           continue    
#         # create the spliced sequence
#         spliced_seq = spliceSeq(full_seq_record, exon)
#         if spliced_seq is None:
#           print(f"Error in splicing transcript {transcript} of gene {row['GeneID']}. Skipping")
#           with open('noSplicing.txt', 'a') as f:
#             f.write(transcript + '\n')
#           continue
#         if len(spliced_seq) % 3 != 0:
#           print(f"Length of spliced seq is {len(spliced_seq)}. not multiple of 3. Skipping.")
#           with open ("noMulByThree.txt", 'a') as f:
#             f.write(transcript + '\n')
#           continue
#         # create AA sequences
#         try:
#           AA_full_seq = full_seq_record.seq.translate(cds=True)
#         except:
#           print("Error in translating full RNA seq. Skipping")
#           continue
#         AA_spliced_seq = Seq(spliced_seq)
#         try:
#           AA_spliced_seq = AA_spliced_seq.translate(cds=True)
#         except:
#           print("Error in translating spliced RNA seq. Skipping.")
#           continue
#         # create directories
#         gene_dir = os.path.join(user_args.output_dir,row['GeneID']+'_'+row['geneSymbol'])
#         if not os.path.isdir(gene_dir):
#           os.mkdir(gene_dir)
#         id_dir = os.path.join(gene_dir, row['ID']) 
#         if not os.path.isdir(id_dir):
#           os.mkdir(id_dir)
#         transcript_dir = os.path.join(id_dir, transcript)
#         # check if current transcript was already checked in previous running
#         if os.path.isdir(transcript_dir):
#           print(f"Transcript {transcript} of gene {row['GeneID']} was already checked.\nResults can be found in {transcript_dir}.\nSkipping to next transcript.")
#           continue
#         # save NA sequences
#         os.mkdir(transcript_dir)
#         na_fasta_path = os.path.join(transcript_dir,transcript+"_fullNA_"+type[0]+".fasta")
#         with open(na_fasta_path, 'w') as na_fasta_file:
#           na_fasta_file.write(NA_seq)
#         spliced_seq_record = ">"+row['GeneID']+"_"+transcript+"_"+type[1]+'\n'+str(spliced_seq)+'\n'
#         spliced_na_fasta_path = os.path.join(transcript_dir,transcript+"_splicedNA_"+type[1]+".fasta")
#         with open(spliced_na_fasta_path, 'w') as spliced_na_fasta_file:
#           spliced_na_fasta_file.write(spliced_seq_record)
#         # save AA sequences
#         AA_full_seq_record = ">"+row['GeneID']+"_"+transcript+"_"+type[0]+'\n'+str(AA_full_seq)+'\n'
#         AA_spliced_seq_record = ">"+row['GeneID']+"_"+transcript+"_"+type[1]+'\n'+str(AA_spliced_seq)+'\n'
#         AA_full_seq_path = os.path.join(transcript_dir,transcript+'_fullAA_'+type[0]+'.fasta')
#         AA_spliced_seq_path = os.path.join(transcript_dir,transcript+'_splicedAA_'+type[1]+'.fasta')
#         with open(AA_full_seq_path, 'w') as f:
#           f.write(AA_full_seq_record)
#         with open(AA_spliced_seq_path, 'w') as f:
#           f.write(AA_spliced_seq_record)
#         # save AA spliced peptide
#         spliced_peptide, start,end = find_spliced_peptide(str(AA_full_seq),str(AA_spliced_seq))
#         spliced_peptide_record = "Spliced peptide coordinates: "+str(start)+"-"+str(end)+"\n"+str(spliced_peptide)
#         spliced_peptide_path=os.path.join(transcript_dir,"spliced_peptide.txt")
#         with open(spliced_peptide_path, 'w') as f:
#           f.write(spliced_peptide_record)
# print("Done analyze transcripts.")