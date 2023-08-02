import os
import pandas as pd

# Replace 'your_file.csv' with the actual path to your CSV file.
csv_file_path = '/private5/Projects/Efi/lostfound_private5_rahamie4.csv'

# Read the CSV file into a DataFrame.
df = pd.read_csv(csv_file_path)

# Assuming the "FullPath" column is named "FullPath" in the CSV file.
full_path_column = 'FullPath'

# Loop through each row and display the content of the corresponding path.
for index, row in df.iterrows():
    full_path = row[full_path_column]
    print(f"Content of Path '{full_path}':")
    
    try:
        # List the content of the directory using os.listdir.
        content = os.listdir(full_path)
        
        # Display the content (list of files and directories).
        print(content)
    except FileNotFoundError:
        print(f"Directory not found: '{full_path}'")
    except PermissionError:
        print(f"Permission denied for directory: '{full_path}'")
    except Exception as e:
        print(f"An error occurred while accessing '{full_path}': {e}")
    
    print("------------------------------------------------------")

