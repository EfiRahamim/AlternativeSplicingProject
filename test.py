import re, requests, subprocess, os
from io import StringIO
from Bio import SeqIO

# 1. read the GTF file and find the entry with the required transcript
# 2. read the next entries specifing the exons of the required transcript and achieve each exon sequence by using the getExon function.
# 3. concatenate all the exons into one sequence
# 4. define the possible ORF's in the sequence
# 5. Maybe: Choose the longest ORF as the most possible sequence
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

# function to get the sequence of a novel transcript using the GTF file
def getNovelTranscriptFasta(novel_transcript_id):
  with open(gtf_file, 'r') as gtf:
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



server = "https://rest.ensembl.org"
gtf_file="/private10/Projects/Efi/AML/StringTie/Control_SF_Mutations/Control.gtf"
print(HandleNovelTranscriptCase("ENST00000361445.8"))



