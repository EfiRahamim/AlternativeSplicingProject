import re, requests

server = "https://rest.ensembl.org"

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
def getNovelTranscript(novel_transcript_id):
  # 1. read the GTF file and find the entry with the required transcript
  # 2. read the next entries specifing the exons of the required transcript and achieve each exon sequence by using the getExon function.
  # 3. concatenate all the exons into one sequence
  # 4. define the possible ORF's in the sequence
  # 5. Maybe: Choose the longest ORF as the most possible sequence
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
        print(f"Error in capturing transcript_id attribute in attributes {attributes}")
        continue
      transcript_id = matche.group(1) # capture the transcript
      if feature == "transcript" and transcript_id == novel_transcript_id:
        transcript_found = True
        print(f"Transcript {novel_transcript_id} was found in the GTF file")
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
        print("Done creating sequence")
        break
    return novel_transcript_seq

gtf_file="/private10/Projects/Efi/AML/StringTie/Control_SF_Mutations/Control.gtf"
novel_transcript_seq = getNovelTranscript("ENST00000428771.6")
print(novel_transcript_seq)


