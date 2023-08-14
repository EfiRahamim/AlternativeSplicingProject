from Bio import SeqIO, Seq
import os, re
import requests, sys, csv, argparse, subprocess, multiprocessing
from io import StringIO
from Bio.Seq import Seq

# CLI arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Analyzing transcripts of rMATS results: Creating full and spliced nucleic and amino acids sequences for each transcript. Results will be saved in the giving output directory.\nNOTE: NOT SUITABLE FOR MXE EVENTS (YET)!!!")
parser.add_argument("-i", action='store', dest='input', required=True, help="Input file: rMATS filtered csv file")
parser.add_argument("-type", action='store', dest='as_type', required=True, help="Type of Splicing Event - SE/A5SS/A3SS/MXE/RI.", choices=["SE","A3SS","A5SS","MXE","RI"])
parser.add_argument("-gtf", action='store', dest='gtf', required=True, help="GTF file for handeling novel transcripts. Should be GTF file from StringTie output.")
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

def getExonByRegions(chr, strand_sign, exonStart, exonEnd):
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
    print(f"Error in getting exon by regions.")
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
  # make FASTA record for the sequence
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

def getAminoAcidSeq(inclusion_seq_object, exclusion_seq_object):
  # create AA sequences
  try:
    if (not user_args.mulby3) or (user_args.as_type == 'RI'):
      AA_inclusion_seq = inclusion_seq_object.translate(stop_symbol="") # retainded intron can affects the ORF, causing short isoform
    else:
      AA_inclusion_seq = inclusion_seq_object.translate(cds=True) 
  except:
    print("Error in translating inclusion seq. Skipping")
    return None, None
  #AA_spliced_seq = Seq(exclusion_seq_record)
  try:
    if not user_args.mulby3: # if spliced seq is not a multiple of 3 - short peptide can be formed
      AA_exclusion_seq = exclusion_seq_object.translate(stop_symbol="")  
    else:
      AA_exclusion_seq = exclusion_seq_object.translate(cds=True)
  except:
    print("Error in translating exclusion seq. Skipping.")
    return None, None
  return AA_inclusion_seq, AA_exclusion_seq 

# function to get the sequence of a novel transcript using the GTF file
def getNovelTranscriptFasta(novel_transcript_id):
  with open(user_args.gtf, 'r') as gtf:
    transcript_found = False
    novel_transcript_seq = ""
    for line in gtf:
      if line.startswith("#"):
        continue
      line_fields = line.strip().split('\t') # split the line by tab delimeter and get the dields of the line
      feature = line_fields[2] # feature type: transcript or exon
      attributes = line_fields[8] # attributes of feature (last column in the line)
      transcript_id_pattern = r"transcript_id \"(\w+\.\d+\.*\d*)\"" # regex to capture the transcript_id within the attributes
      matche = re.search(transcript_id_pattern, attributes) # search for transcript_id string in the attributes
      if matche == None:
        print(f"Error in capturing transcript_id attribute in {attributes}")
        continue
      transcript_id = matche.group(1) # capture the transcript
      if feature == "transcript" and transcript_id == novel_transcript_id:
        transcript_found = True
        print(f"Novel transcript {novel_transcript_id} was found in the GTF file")
        continue
      if feature == "exon" and transcript_found:
        chr, strand, exonStart, exonEnd = line_fields[0], line_fields[6],line_fields[3], line_fields[4]
        exon = getExon(chr, strand, exonStart, exonEnd, "NoGene", "NoAStype") # get the exon
        if strand == "+":
          novel_transcript_seq = novel_transcript_seq+exon # positive strand: add the exon to the sequence, at the end
        elif strand == "-":
          novel_transcript_seq = exon+novel_transcript_seq # negative strand: add the exon to the sequence, at the beggining
        print(f"Adding exon: {exonStart}-{exonEnd}")
        continue
      if feature == "transcript" and transcript_found:
        print(f"Done creating sequence of novel transcript {novel_transcript_id}")
        break
  # check if any transcript was found
  if not transcript_found:
    return None
  # save the transcript in a fasta file format 
  if not os.path.isdir(os.path.join(os.getcwd(), "NovelTranscrips")):
    os.mkdir(os.path.join(os.getcwd(), "NovelTranscrips"))
  novel_transcript_fasta = "> " + novel_transcript_id + "\n" + novel_transcript_seq
  output_file = os.path.join(os.getcwd(),"NovelTranscrips", f'{novel_transcript_id}_transcript.fasta')
  print(f"Writing temporary FASTA file of novel transcript to: {output_file}. File will be deleted after checking.") 
  with open(output_file, 'w') as fasta_file:
    fasta_file.write(novel_transcript_fasta)
  return output_file

# function to find the most possible ORF of the novel transcript
def findORF(novel_transcript_fasta):
  # use the ORFfinder tool (installed via conda: /home/alu/rahamie4/anaconda3/envs/ORFfinder)
  command = ["ORFfinder", "-in", novel_transcript_fasta, "-s", "0", "-n", "TRUE","-outfmt", "1"]
  try:
    output = subprocess.check_output(command, universal_newlines=True)
  except:
    print("Error in running ORFfinder on subprocess. Exit.")
    exit
  # parse the output of ORFfinder into FAST file
  ORFfinder_output = StringIO(output)
  records = list(SeqIO.parse(ORFfinder_output, "fasta"))
  # find the longest ORF
  ORF_lengths = [len(record.seq) for record in records]
  longest_ORF_index = ORF_lengths.index(max(ORF_lengths))
  longest_ORF_record = records[longest_ORF_index]
  # check if most possible ORF is on the positive strand (more accurate) or on the negative strand (less accurate and need to be checked)
  onPositiveStrand = checkORFStrand(longest_ORF_record.id) 
  if onPositiveStrand == None:
    return None
  if not onPositiveStrand:
    print("Most possible ORF that has been found is not on the positive strand.")
    return None
  # return the sequence of the most possible ORF in the given novel transcript sequence
  return longest_ORF_record.seq

# function to check the strand of the ORF
def checkORFStrand(orf_fasta_id):
  pattern = r":(\d+)-(\d+)$"
  match = re.search(pattern, orf_fasta_id)
  if match == None:
    print(f"Error in retreving ORF position of ORF ID: {orf_fasta_id}")
    return None
  start_pos = int(match.group(1))
  end_pos = int(match.group(2))
  if start_pos < end_pos:
    return True
  else:
    return False

# function to handle case of novel transcript  
def HandleNovelTranscriptCase(novel_transcript_id):
  # get fasta file of the novel transcript
  novel_transcript_file = getNovelTranscriptFasta(novel_transcript_id)
  if novel_transcript_file ==  None:
    print(f"Error in creating novel transcript sequence for transcript: {novel_transcript_id}. Skipping.")
    return
  # find the most possible ORF
  novel_transcript_ORF = findORF(novel_transcript_file)
  # remove FASTA file of novel transcript
  #os.remove(novel_transcript_file)
  # check the returned ORF 
  if novel_transcript_ORF == None:
    return None
  # return the ORF of the novel transcript
  return novel_transcript_ORF

# function to handle case of annotated transcript
def HandleAnnotatedTranscriptCase(annotated_transcript_id):
  # create URL command for retreving transcript sequence from Ensemble
  ext_transcript = f"/sequence/id/{annotated_transcript_id.replace('.','/.')}?type=cds"
  print(f"Annotated transcript sequence request: {server+ext_transcript}")
  r_t = requests.get(server+ext_transcript, headers={ "Content-Type" : "text/x-fasta"})
  if not r_t.ok:
    #r_t.raise_for_status()
    #sys.exit
    print(f"Error in: {server+ext_transcript}")
    return None
  # parse the output as FASTA text
  fasta_str = StringIO(r_t.text)
  # Parse the contents of the StringIO object as a FASTA file
  records = list(SeqIO.parse(fasta_str, "fasta"))
  # return the annotated transcript sequence
  return records[0].seq

# function to create the inclusion sequence from the transcript
def getInclusionExclusionSeq(transcript_seq, row, as_type, transcript_id):
  # The idea: look for the upstream and downstream sequences in the transcript and join the alternative exon between them. works both if the current transcript is the inclusion or exlusion version
  
  # get the upstream, downstream and alternative exon sequence
  if as_type in ['SE', 'MXE']:
    upstream_seq = getExonByRegions(row['chr'], row['strand'],str(int(row['upstreamES'])+1), row['upstreamEE'])
    downstream_seq = getExonByRegions(row['chr'], row['strand'],str(int(row['downstreamES'])+1), row['downstreamEE'])
    if as_type == 'SE':
      alternative_exon_seq = row['AlternativeExonSeq']
    elif as_type == 'MXE':
      alternative_exon_seq1 = row['AlternativeExonSeq1']
      alternative_exon_seq2 = row['AlternativeExonSeq2']
  elif as_type == 'A5SS':
    alternative_exon_seq = row['AlternativeExonSeq']
    upstream_seq = getExonByRegions(row['chr'], row['strand'], exonStart=row['shortES_1base'], exonEnd=row['shortEE'])
    downstream_seq = getExonByRegions(row['chr'], row['strand'], exonStart=str(int(row['flankingES'])+1), exonEnd=row['flankingEE'])
  elif as_type == 'A3SS':
    alternative_exon_seq = row['AlternativeExonSeq']
    upstream_seq = getExonByRegions(row['chr'], row['strand'], exonStart=str(int(row['flankingES'])+1), exonEnd=row['flankingEE'])
    downstream_seq = getExonByRegions(row['chr'], row['strand'], exonStart=row['shortES_1base'], exonEnd=row['shortEE'])
  elif as_type == 'RI':
    alternative_exon_seq = row['RetainedIntronSeq']
    upstream_seq = row['upstreamExonSeq']
    downstream_seq = row['downstreamExonSeq']
  
  # look for the upstream&downstream sequences and create the inclusion and exclusion sequences
  if as_type in ['SE', 'RI', 'MXE']:
    if row['strand'] == '+':
      # pattern to look from the upstream+downstream sequence in the positive strand
      up_and_down_seq_pattern = fr"({upstream_seq}.*{downstream_seq})"
      matche = re.search(up_and_down_seq_pattern,transcript_seq)
      if matche == None:
        print(f"Error in finding upstream and downstream exons in transcript: {transcript_id}, Gene: {row['GeneID']}, AS Type: {as_type}")
        return None
      up_and_down_seq = matche.group(0)
      if as_type == 'MXE':
        # MXE type in forward strand: the inclusion form includes the 1st exon and skips the 2nd exon
        inclusion_seq = transcript_seq.replace(up_and_down_seq, upstream_seq+alternative_exon_seq1+downstream_seq)
        exclusion_seq = transcript_seq.replace(up_and_down_seq, upstream_seq+alternative_exon_seq2+downstream_seq)
      else:
        # other cases: SE & RI - inclusion form include the exon and exclusion exclude the exon
        inclusion_seq = transcript_seq.replace(up_and_down_seq, upstream_seq+alternative_exon_seq+downstream_seq)
        exclusion_seq = transcript_seq.replace(up_and_down_seq, upstream_seq+downstream_seq)
    elif row['strand'] == '-':
      # pattern to look from the upstream+downstream sequence in the negative strand
      up_and_down_seq_pattern = fr"({downstream_seq}.*{upstream_seq})"
      matche = re.search(up_and_down_seq_pattern,transcript_seq)
      if matche == None:
        print(f"Error in finding upstream and downstream exons in transcript: {transcript_id}, Gene: {row['GeneID']}, AS Type: {as_type}")
        return None, None
      up_and_down_seq = matche.group(0)
      if as_type == 'MXE':
        # MXE type in negative strand: the inclusion form includes the 2st exon and skips the 1nd exon
        inclusion_seq = transcript_seq.replace(up_and_down_seq, upstream_seq+alternative_exon_seq2+downstream_seq)
        exclusion_seq = transcript_seq.replace(up_and_down_seq, upstream_seq+alternative_exon_seq1+downstream_seq)
      else:
        # other cases: SE & RI - inclusion form include the exon and exclusion exclude the exon
        inclusion_seq = transcript_seq.replace(up_and_down_seq, downstream_seq+alternative_exon_seq+upstream_seq)
        exclusion_seq = transcript_seq.replace(up_and_down_seq, downstream_seq+upstream_seq)
  elif as_type in ['A5SS', 'A3SS']:
    # A5SS/A3SS cases: all sequences are designed as positive
    inclusion_seq = transcript_seq.replace(up_and_down_seq, upstream_seq+alternative_exon_seq+downstream_seq)
    exclusion_seq = transcript_seq.replace(up_and_down_seq, upstream_seq+downstream_seq)
  
  # check inclusion and exclusion sequences not same length
  if len(inclusion_seq) == len(exclusion_seq):
    print(f"Length of inclusion and exlusion sequences are the same. Skipping")
    return None, None
  else:
    return inclusion_seq, exclusion_seq


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
      # check wether is an annotated transcript or noval transcript
      if not transcript.startswith("ENST"):
        # novel transcript
        transcript_seq = HandleNovelTranscriptCase(transcript)
      else:
        # annotated transcript
        transcript_seq = HandleAnnotatedTranscriptCase(transcript)
      # check validity of transcript seq
      if transcript_seq == None:
        print(f"Error in retreving transcript sequences of {transcript}. Skipping.")
        with open('NoTranscripSeq.txt', 'a') as f:
          f.write(transcript + '\n')
        continue
      # create the full and spliced sequenced from the transcript
      inclusion_seq, exclusion_seq = getInclusionExclusionSeq(transcript_seq, row, user_args.as_type)
      if inclusion_seq == None or exclusion_seq == None:
        print(f"Error in getting inclusion and/or exclusion sequences. Skipping.")
        with open('NoInclusionExclusionSeq.txt', 'a') as f:
          f.write(transcript + '\n')
          continue
      # # get full sequence of current transcript
      # NA_seq, full_seq_record = getSeq(row, transcript, user_args.as_type)
      # if full_seq_record is None:
      #   continue
      # # create the spliced sequence
      # spliced_seq = spliceSeq(full_seq_record, exon)
      # if spliced_seq is None:
      #   print(f"Error in splicing transcript {transcript} of gene {row['GeneID']}. Skipping")
      #   with open('noSplicing.txt', 'a') as f:
      #     f.write(transcript + '\n')
      #   continue
      # check if inclusion/exclusion sequences are multiple of 3 (only if flag '--mulby3' passed by the user)
      if (user_args.mulby3) and ((len(inclusion_seq) % 3 != 0) or (len(exclusion_seq) % 3 != 0)):
        print(f"Length of Inclusion/Exclusion sequences are not multiple of 3. Skipping.")
        with open ("noMulByThree.txt", 'a') as f:
          f.write(transcript + '\n')
        continue
      # get Amino Acid sequences
      inclusion_seq_object = Seq(inclusion_seq)
      exclusion_seq_object = Seq(exclusion_seq)
      AA_inclusion_seq, AA_exclusion_seq = getAminoAcidSeq(inclusion_seq_object, exclusion_seq_object)
      if AA_inclusion_seq == None or AA_exclusion_seq == None:
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
      inclusion_seq_record = ">"+row['GeneID']+"_"+transcript+"_"+type[0]+'\n'+str(inclusion_seq)+'\n'
      inclusion_na_fasta_path = os.path.join(transcript_dir,transcript+"_inclusionNA_"+type[0]+".fasta")
      with open(inclusion_na_fasta_path, 'w') as inclusion_na_fasta_file:
        inclusion_na_fasta_file.write(inclusion_seq_record)
      exclusion_seq_record = ">"+row['GeneID']+"_"+transcript+"_"+type[1]+'\n'+str(exclusion_seq)+'\n'
      exclusion_na_fasta_path = os.path.join(transcript_dir,transcript+"_exclusionNA_"+type[1]+".fasta")
      with open(exclusion_na_fasta_path, 'w') as exclusion_na_fasta_file:
        exclusion_na_fasta_file.write(exclusion_seq_record)
      # save AA sequences
      AA_inclusion_seq_record = ">"+row['GeneID']+"_"+transcript+"_"+type[0]+'\n'+str(AA_inclusion_seq)+'\n'
      AA_exclusion_seq_record = ">"+row['GeneID']+"_"+transcript+"_"+type[1]+'\n'+str(AA_exclusion_seq)+'\n'
      AA_inclusion_seq_path = os.path.join(transcript_dir,transcript+'_inclusionAA_'+type[0]+'.fasta')
      AA_exclusion_seq_path = os.path.join(transcript_dir,transcript+'_exclusionAA_'+type[1]+'.fasta')
      with open(AA_inclusion_seq_path, 'w') as f:
        f.write(AA_inclusion_seq_record)
      with open(AA_exclusion_seq_path, 'w') as f:
        f.write(AA_exclusion_seq_record)
      # save AA spliced (alternative) peptide
      spliced_peptide, start,end = find_spliced_peptide(str(AA_inclusion_seq),str(AA_exclusion_seq))
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