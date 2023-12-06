import pandas as pd
def parse_gtf(filename):
  columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
  data = []
  with open(filename, 'r') as file:
      for line in file:
          if line.startswith('#'):
              continue  # Skip comment lines
          parts = line.strip().split('\t')
          parts[8]=parts[8].replace('"', "")
          attributes = dict(item.strip().split(' ') for item in parts[8].split(';') if item.strip())
          parts[8] = attributes  # Replace 'attribute' field with a dictionary
          data.append(parts)
  df = pd.DataFrame(data, columns=columns)
  return df
gtf_df = parse_gtf("/private10/Projects/Efi/General/test.gtf")
for index, row in gtf_df.iterrows():
    print(row['attribute']['gene_id'])