import argparse
import csv
import subprocess, multiprocessing
import requests
import re
from io import StringIO
from Bio.Seq import Seq
from Bio import SeqIO

# CLI arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Matching optional transcripts for each Alternative Splicing event, based on the SE_to_isoform.py script of rMATS developers (took from GitHub: https://github.com/Xinglab/rmats-turbo/files/7237006/se_to_isoforms.zip). The matching transcripts will be added to the end of the table, in a new column. Detecting also novel transcripts. Conda environment ORFfinder must be activated!")
parser.add_argument("-se", action='store', dest='SE_output', required=False, help="Input file: SE filtered csv file")
parser.add_argument("-a5ss", action='store', dest='A5SS_output', required=False, help="Input file: A5SS filtered csv file")
parser.add_argument("-a3ss", action='store', dest='A3SS_output', required=False, help="Input file: A3SS filtered csv file")
parser.add_argument("-mxe", action='store', dest='MXE_output', required=False, help="Input file: MXE filtered csv file")
parser.add_argument("-ri", action='store', dest='RI_output', required=False, help="Input file: RI filtered csv file")
parser.add_argument('-gtf', action='store', dest='gtf', help="GTF file that was used for running STAR and rMATS. If detecting novel transcript, GTF should contain novel transcripts as well (can be obtained from StringTie)")
parser.add_argument('-script', action='store', dest='script', default="/private5/Projects/Efi/AS/rmats-turbo/se_to_isoforms.py", help='The \'SE_to_isoform.py\' script')
#SE_output = "/private5/Projects/Efi/AS/TCGA-BRCA/gencode_v36/significant_results/SE/SE_filtered.csv"
#gtf = "/private5/Projects/Efi/AS/gencode.v36.annotation.gtf"

user_args = parser.parse_args()

server = "https://rest.ensembl.org"

# Function to run the gtf_to_transcript.py script and parse the output
def run_gtf_to_transcript(chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE):
    # create the shell command
    command = ["python", user_args.script , "--chr", str(chr), "--strand", strand, "--exon-start", str(exonStart_0base),
               "--exon-end", str(exonEnd), "--upstream-start", str(upstreamES), "--upstream-end", str(upstreamEE),
               "--downstream-start", str(downstreamES), "--downstream-end", str(downstreamEE), "--gtf", user_args.gtf]
    # run the command
    output = subprocess.check_output(command, universal_newlines=True)
    transcripts = []
    # split output of command to transcripts ID's
    for line in output.strip().split('\n'):
        if line.startswith("transcript_id"):
            # find the transcript_id using regex
            transcript_pattern = r"transcript_id: \"(\w+\.\d+\.*\d*)\","
            matche = re.search(transcript_pattern, line)
            transcript_id = matche.group(1)
            transcripts.append(transcript_id)
    return transcripts

def getExon(chr, strand_sign, exonStart, exonEnd, geneID, as_type):
  # define strand as 1/-1
  if strand_sign == '+':
    strand = '1'
  elif strand_sign == '-':
    strand = '-1'
  # get sequence by region - 1base (!!!)
  regions = f"{chr}:{exonStart}..{exonEnd}:{strand}"
  ext_region = f"/sequence/region/human/{regions}?"
  print(f"Region seq request: {server+ext_region}")
  r_r = requests.get(server+ext_region, headers={ "Content-Type" : "text/plain"}) 
  if not r_r.ok:
    #r_r.raise_for_status()
    #sys.exit
    print(f"Error in getting exon by regions. Skipping.")
    with open (f"NoExon_{as_type}.txt", 'a') as f:
      f.write(geneID + '\n')
    return None
  exon = r_r.text.upper()
  return exon

def findMatchingTranscripts(row, transcripts, exon, as_type):
  # get sequence by trascript ID
  # transcripts = row['optional transcripts'].split(";")
  found = False
  matching_transcripts = []
  # run over each transcript
  for transcript in transcripts:
    # create and run URL 
    transcript = transcript.replace(".","/.")
    ext_transcript = f"/sequence/id/{transcript}?type=cds"
    print(f"Transcript sequence request: {server+ext_transcript}")
    r_t = requests.get(server+ext_transcript, headers={ "Content-Type" : "text/x-fasta"})
    if not r_t.ok:
      #r_t.raise_for_status()
      #sys.exit
      print(f"Error in: {server+ext_transcript}")
      continue
    fasta_str = StringIO(r_t.text)
    # Parse the contents of the StringIO object as a FASTA file
    records = list(SeqIO.parse(fasta_str, "fasta"))
    if exon != None and exon in records[0].seq:
      found = True
      transcript = transcript.replace("/.",".")
      print(f"Gene: {row['GeneID']} Matched transcript: {transcript}")
      matching_transcripts.append(transcript)
      #break
      #return r_t.text, records[0], transcript
    else:
      continue
  if not found: # check if at least one matching transcript was found
    print(f"No matching transcript were found for {row['GeneID']}")
    with open (f"NoTranscripts_{as_type}.txt", 'a') as f:
      f.write(row['GeneID'] + '\n')
    #return None, None, None
  return matching_transcripts

def get_command_args(row, as_type, mxe_second_exon = False):
  # define command parameters according to each AS type
  if as_type == 'SE':
    chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE = row['chr'], row['strand'], row['exonStart_0base'], row['exonEnd'], row['upstreamES'], row['upstreamEE'], row['downstreamES'], row['downstreamEE']
  elif (as_type == 'A5SS' and row['strand'] == '+') or (as_type=='A3SS' and row['strand'] == '-'):
    chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE = row['chr'], row['strand'], row['shortEE'], row['longExonEnd'], row['shortES'], row['shortEE'], row['flankingES'], row['flankingEE']
  elif (as_type == 'A5SS' and row['strand'] == '-') or (as_type=='A3SS' and row['strand'] == '+'):
    chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE = row['chr'], row['strand'], row['longExonStart_0base'], row['shortES'], row['flankingES'], row['flankingEE'], row['shortES'], row['shortEE']
  elif as_type == 'RI':
    chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE = row['chr'], row['strand'], row['riExonStart_0base'], row['riExonEnd'], row['upstreamES'], row['upstreamEE'], row['downstreamES'], row['downstreamEE']
  elif as_type == 'MXE' and not mxe_second_exon:
    chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE = row['chr'], row['strand'], row['X1stExonStart_0base'], row['X1stExonEnd'], row['upstreamES'], row['upstreamEE'], row['downstreamES'], row['downstreamEE']
  elif as_type == 'MXE' and mxe_second_exon:
    chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE = row['chr'], row['strand'], row['X2ndExonStart_0base'], row['X2nsExonEnd'], row['upstreamES'], row['upstreamEE'], row['downstreamES'], row['downstreamEE']
  else:
    print(f"Error in define type of splicing - {as_type} in function 'get_command_args'. skipping.")
    return None
  return chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE

def file_to_list(rmats_output_file):
  # open file for reading and extract the rows from the file
  with open(rmats_output_file, 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    rows = list(reader)
  return rows

def define_exon_coordinates(row, as_type, mxe_second_exon = False, ri_second_exon = False):
  # define alternative exon coordinates for each AS type
  if as_type == 'SE':
    exonStart, exonEnd = row['exonStart_1base'], row['exonEnd']
  elif as_type == 'RI' and not ri_second_exon:
    exonStart, exonEnd = row['riExonStart_1base'], row['upstreamEE']
  elif as_type == 'RI' and ri_second_exon:
    exonStart, exonEnd = row['downstreamES_1base'], row['downstreamEE']
  elif (as_type  == 'A5SS' and row['strand'] == '+') or (as_type == 'A3SS' and row['strand'] == '-'):
    exonStart, exonEnd = row['shortEE_1base'], row['longExonEnd']
  elif (as_type  == 'A5SS' and row['strand'] == '-') or (as_type == 'A3SS' and row['strand'] == '+'):
    exonStart, exonEnd = row['longExonStart_1base'], row['shortES']
  elif as_type == 'MXE' and not mxe_second_exon:
    exonStart, exonEnd = row['X1stExonStart_1base'], row['X1stExonEnd']
  elif as_type == 'MXE' and mxe_second_exon:
    exonStart, exonEnd = row['X2ndExonStart_1base'], row['X2ndExonEnd']
  else:
    print(f"Error in define type of splicing - {as_type} in function 'define_exon_coordinates'. skipping.")
    return None
  return exonStart, exonEnd

def run_one_row(row, as_type):
  print(f"Running on ID: {row['ID']}, GeneID: {row['GeneID']}")
  # define command arguments
  chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE = get_command_args(row, as_type)
  # run command and fine transcripts
  transcripts = run_gtf_to_transcript(chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE)
  print(f"Optional transcripts: {transcripts}.")
  # find matching transcript according to the AS type
  if as_type in ['SE', 'A5SS', 'A3SS']:
    exonStart, exonEnd = define_exon_coordinates(row, as_type)
    exon_seq = getExon(row['chr'], row['strand'], exonStart, exonEnd, row['GeneID'], as_type)
    row['AlternativeExonSeq'] = exon_seq
    print(f"Looking for matching transcripts for GeneID {row['GeneID']} in {as_type} type.") 
    matching_transcripts = findMatchingTranscripts(row, transcripts, exon_seq, as_type)
  elif as_type == 'RI':
    exonStart1, exonEnd1 = define_exon_coordinates(row, as_type, ri_second_exon=False)
    exonStart2, exonEnd2 = define_exon_coordinates(row,as_type, ri_second_exon=True)
    exon1_seq = getExon(row['chr'], row['strand'], exonStart1, exonEnd1, row['GeneID'], as_type)
    exon2_seq = getExon(row['chr'], row['strand'], exonStart2, exonEnd2, row['GeneID'], as_type)
    retained_intron = getExon(row['chr'], row['strand'], row['upstreamEE_1base'], row['downstreamES'], row['GeneID'], as_type)
    row['upstreamExonSeq'], row['downstreamExonSeq'], row['RetainedIntronSeq'] = exon1_seq, exon2_seq, retained_intron
    print(f"Looking for matching transcripts for GeneID {row['GeneID']} in {as_type} type.")
    optional_transcripts1 = findMatchingTranscripts(row, transcripts, exon1_seq, as_type)
    optional_transcripts2 = findMatchingTranscripts(row, transcripts, exon2_seq, as_type)
    matching_transcripts = list(set(optional_transcripts1).intersection(optional_transcripts2))
  elif as_type == 'MXE':
    exonStart1, exonEnd1 = define_exon_coordinates(row, as_type, mxe_second_exon=False)
    exonStart2, exonEnd2 = define_exon_coordinates(row,as_type, mxe_second_exon=True)
    exon1_seq = getExon(row['chr'], row['strand'], exonStart1, exonEnd1, row['GeneID'], as_type)
    exon2_seq = getExon(row['chr'], row['strand'], exonStart2, exonEnd2, row['GeneID'], as_type)
    row['AlternativeExonSeq1'] = exon1_seq
    row['AlternativeExonSeq2'] = exon2_seq
    print(f"Looking for matching transcripts for GeneID {row['GeneID']} in {as_type} type.")
    optional_transcripts1 = findMatchingTranscripts(row, transcripts, exon1_seq, as_type)
    optional_transcripts2 = findMatchingTranscripts(row, transcripts, exon2_seq, as_type)
    uniqe_matching_transcripts = list(set(optional_transcripts1).symmetric_difference(optional_transcripts2))
    both_exon_transcripts = list(set(optional_transcripts1).intersection(optional_transcripts2))
    matching_transcripts_group1 = [transcript for transcript in uniqe_matching_transcripts if transcript in optional_transcripts1]
    matching_transcripts_group2 = [transcript for transcript in uniqe_matching_transcripts if transcript in optional_transcripts2]
    print(f"Transcripts mathcing both exons: {both_exon_transcripts}")
    print(f"Unique transcripts: {uniqe_matching_transcripts}")
    print(f"Transcript for first exon: {matching_transcripts_group1}")
    print(f"Transcript for second exon: {matching_transcripts_group2}")
    if len(matching_transcripts_group1) == 0 or len(matching_transcripts_group2) == 0:
      print("No matching transcript for exon1 or exon2. Returning transcripts that contain both exons.")
      matching_transcripts = both_exon_transcripts
    else:
      matching_transcripts = uniqe_matching_transcripts
  row['optional transcripts'] = ';'.join(matching_transcripts)
  print(f"Matching transcripts: {matching_transcripts}")
  return row

def run_analyse(rmats_output_file, as_type):
  rows = file_to_list(rmats_output_file) # read the file
  #print(rows)
  print(f"---------------------------------Start Proccesing {as_type}---------------------------------")
  pool = multiprocessing.Pool(processes=10) 
  rows_update = pool.starmap(run_one_row,[(row,as_type) for row in rows])
  pool.close()
  pool.join()
  # for row in rows: # run over each row in the file
  #   print(f"Running on ID: {row['ID']}, GeneID: {row['GeneID']}")
  #   # define command arguments
  #   chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE = get_command_args(row, as_type)
  #   # run command and fine transcripts
  #   transcripts = run_gtf_to_transcript(chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE)
  #   # find matching transcript according to the AS type
  #   if as_type in ['SE', 'A5SS', 'A3SS']:
  #     exonStart, exonEnd = define_exon_coordinates(row, as_type)
  #     exon_seq = getExon(row['chr'], row['strand'], exonStart, exonEnd, row['GeneID'])
  #     row['AlternativeExonSeq'] = exon_seq 
  #     matching_transcripts = findMatchingTranscripts(row, transcripts, exon_seq)
  #   elif as_type == 'RI':
  #     exonStart1, exonEnd1 = define_exon_coordinates(row, as_type, ri_second_exon=False)
  #     exonStart2, exonEnd2 = define_exon_coordinates(row,as_type, ri_second_exon=True)
  #     exon1_seq = getExon(row['chr'], row['strand'], exonStart1, exonEnd1, row['GeneID'])
  #     exon2_seq = getExon(row['chr'], row['strand'], exonStart2, exonEnd2, row['GeneID'])
  #     retained_intron = getExon(row['chr'], row['strand'], row['upstreamEE_1base'], row['downstreamES'], row['GeneID'])
  #     row['upstreamExonSeq'], row['downstreamExonSeq'], row['RetainedIntronSeq'] = exon1_seq, exon2_seq, retained_intron
  #     optional_transcripts1 = findMatchingTranscripts(row, transcripts, exon1_seq)
  #     optional_transcripts2 = findMatchingTranscripts(row, transcripts, exon2_seq)
  #     matching_transcripts = list(set(optional_transcripts1).intersection(optional_transcripts2))
  #   elif as_type == 'MXE':
  #     exonStart1, exonEnd1 = define_exon_coordinates(row, as_type, mxe_second_exon=False)
  #     exonStart2, exonEnd2 = define_exon_coordinates(row,as_type, mxe_second_exon=True)
  #     exon1_seq = getExon(row['chr'], row['strand'], exonStart1, exonEnd1, row['GeneID'])
  #     exon2_seq = getExon(row['chr'], row['strand'], exonStart2, exonEnd2, row['GeneID'])
  #     row['AlternativeExonSeq1'] = exon1_seq
  #     row['AlternativeExonSeq2'] = exon2_seq
  #     optional_transcripts1 = findMatchingTranscripts(row, transcripts, exon1_seq)
  #     optional_transcripts2 = findMatchingTranscripts(row, transcripts, exon2_seq)
  #     uniqe_matching_transcripts = list(set(optional_transcripts1).symmetric_difference(optional_transcripts2))
  #     both_exon_transcripts = list(set(optional_transcripts1).intersection(optional_transcripts2))
  #     matching_transcripts_group1 = [transcript for transcript in uniqe_matching_transcripts if transcript in optional_transcripts1]
  #     matching_transcripts_group2 = [transcript for transcript in uniqe_matching_transcripts if transcript in optional_transcripts2]
  #     print(f"Transcripts mathcing both exons: {both_exon_transcripts}")
  #     print(f"Unique transcripts: {uniqe_matching_transcripts}")
  #     print(f"Transcript for first exon: {matching_transcripts_group1}")
  #     print(f"Transcript for second exon: {matching_transcripts_group2}")
  #     if len(matching_transcripts_group1) == 0 or len(matching_transcripts_group2) == 0:
  #       print("No matching transcript for exon1 or exon2. Returning transcripts that contain both exons.")
  #       matching_transcripts = both_exon_transcripts
  #     else:
  #       matching_transcripts = uniqe_matching_transcripts
  #   row['optional transcripts'] = ';'.join(matching_transcripts)
  #   print(f"Matching transcripts: {matching_transcripts}")
  print(f"---------------------------------Done Proccesing {as_type}---------------------------------")
  return rows_update
  
def update_csv_file(original_csv_file, updated_rows, as_type):
  # write to the csv file the new column with the matching transcripts
  with open (original_csv_file, 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    if as_type in ['SE', 'A5SS', 'A3SS']:
      fieldnames = reader.fieldnames + ['AlternativeExonSeq','optional transcripts']
    elif as_type == 'RI':
      fieldnames = reader.fieldnames + ['upstreamExonSeq', 'downstreamExonSeq', 'RetainedIntronSeq', 'optional transcripts']
    elif as_type == 'MXE':
      fieldnames = reader.fieldnames + ['AlternativeExonSeq1', 'AlternativeExonSeq2', 'optional transcripts']
  with open(original_csv_file, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(updated_rows)

as_type_dict = {user_args.SE_output:'SE',
                user_args.A5SS_output:'A5SS',
                user_args.A3SS_output:'A3SS',
                user_args.MXE_output:'MXE',
                user_args.RI_output:'RI'}

for as_event in [as_input for as_input in [user_args.SE_output,user_args.A5SS_output,user_args.A3SS_output,user_args.MXE_output,user_args.RI_output] if as_input]:
  updates_rows = run_analyse(as_event, as_type_dict[as_event])
  update_csv_file(as_event,updates_rows,as_type_dict[as_event])

# if user_args.SE_output:
#   updates_rows = run_analyse(user_args.SE_output, 'SE')
#   update_csv_file(user_args.A5SS_output, updates_rows )
# if user_args.A5SS_output:
#   updates_rows = run_analyse(user_args.A5SS_output, 'A5SS')
#   update_csv_file(user_args.A5SS_output, updates_rows )
# # if user_args.A3SS_output:
# #   updates_rows = run_analyse(user_args.A3SS_output, 'A3SS')
# #   update_csv_file(user_args.A3SS_output, updates_rows, 'A3SS')
# if user_args.MXE_output:
#   updates_rows = run_analyse(user_args.MXE_output, 'MXE')
#   update_csv_file(user_args.MXE_output, updates_rows)
# if user_args.RI_output:
#   updates_rows = run_analyse(user_args.RI_output, 'RI')
#   update_csv_file(user_args.RI_output, updates_rows)


   
  
#def run_se_to_isoform(rmats_output, as_type):
   

# Read the CSV file and add the "optional transcripts" column
# with open(user_args.SE_output, 'r') as csvfile:
#     reader = csv.DictReader(csvfile)
#     rows = list(reader)
#     for row in rows:
#         print(f"Running on ID: {row['ID']}, GeneID: {row['GeneID']}")
#         transcripts = run_gtf_to_transcript(row['chr'], row['strand'], row['exonStart_0base'], row['exonEnd'],
#                                             row['upstreamES'], row['upstreamEE'], row['downstreamES'], row['downstreamEE'],
#                                             user_args.gtf)
#         exon = getExon(row)
#         matching_transcripts = findMatchingTranscripts(row, transcripts, exon)
#         row['optional transcripts'] = ';'.join(matching_transcripts)
#         print(f"Transcripts: {matching_transcripts}")
# # Write the updated CSV file
# fieldnames = reader.fieldnames + ['optional transcripts']
# with open(user_args.SE_output, 'w', newline='') as csvfile:
#     writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
#     writer.writeheader()
#     writer.writerows(rows)

print("Done proccesing transcripts matching.")
