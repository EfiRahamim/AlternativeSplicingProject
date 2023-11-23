import os, glob
#import pandas as pd
# Specify the file path
file_path = "/private10/Projects/Efi/CRG/GBM/StringTie/GroupInfo.txt"
files_to_merge = []
for root, _, files in os.walk("/private10/Projects/Efi/CRG/GBM/StringTie/"):
    for file in files:
        if file.endswith(".gtf"):
            files_to_merge.append(os.path.abspath(os.path.join(root, file)))
print(files_to_merge)
command = ["stringtie --merge",
            "-G",
            "-o",
            "-i"]
command.extend(files_to_merge)
print(command)
