import sys, os, glob, multiprocessing, argparse, subprocess, re, csv

# CLI arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Running netMHC on inclusion & exclusion sequences of each transcript. Results will be saved in the directory of each transcript from the giving output directory. Using the 'netMHC' tool located in: /private/common/Software/netMHC-4.0/Linux_x86_64/bin/netMHC. HLA allels that are being used for this analysis: HLA-A0101,HLA-A0201,HLA-A0301,HLA-A1101,HLA-A2402,HLA-A2601,HLA-B0702,HLA-B0801,HLA-B1501,HLA-B2705,HLA-B3901,HLA-B4001,HLA-B5801. The output of this script is the amount of strong binders that were found in each HLA allel for each group.")
parser.add_argument("-i", action='store', dest='input_dir', required=True, help="Input directory of genes directories")
parser.add_argument("-l1", action='store', dest='lable1', required=True, help="Lable of first group (same order as in the splicing analysis).")
parser.add_argument("-l2", action='store', dest='lable2', required=True, help="Lable of second group (same order as in the splicing analysis).")
parser.add_argument("-rank", action='store', dest='rank', required=False, default='0.5', help="Threshold for high binding peptides (Rank)")
parser.add_argument("-as_type", action='store', dest='as_type', required=False,choices=['SE', 'A5SS', 'A3SS', 'MXE', 'RI'] ,help="Type of splicing event. Required only for rMATS results")
parser.add_argument('--PSIsigma', action='store_true', dest='psi_sigma', help='Set for PSI-Sigma tool results.')
user_args = parser.parse_args()


# get list of absolut pathes of transcripts directories
def getPaths(input_dir):
    #directories_pattern = input_dir+"*/*/*"
    if user_args.psi_sigma:
        directories_pattern = os.path.join(input_dir,"*","*","*","*")
    else:
        directories_pattern = os.path.join(input_dir,"*","*","*")
    directories = glob.glob(directories_pattern)
    pathes = [os.path.abspath(dir) for dir in directories if os.path.isdir(dir)] #and dir.find("sashimiplots") == -1]
    if not pathes: # check if pathes exist
        print("Error in getting pathes of transcripts. Maybe no transcripts exist?")
        sys.exit()
    return pathes

# get list of absolut paths files of current transcript directory
def get_absolute_file_paths(directory):
    file_paths = []
    for root, _, files in os.walk(directory):
        for file in files:
            file_path = os.path.abspath(os.path.join(root, file))
            file_paths.append(file_path)
    return file_paths

# define the type of the inclusion and exclusion sequences - control or case
def getGroupsForms(group1Seq_file):
    groups={"1":"", "2":""}
    if ("Inclusion" in group1Seq_file) or ("inclusion" in group1Seq_file) :
        groups['1'] = "Inclusion"
        groups['2'] = "Exclusion"
    elif ("Exclusion" in group1Seq_file) or ("exclusion" in group1Seq_file) :
        groups['1'] = "Exclusion"
        groups['2'] = "Inclusion"
    else:
        groups=None
    return groups

def get_HLA_strongBinders(netMHC_output,hla_type):
    # # find the index of the line in the output that reffers to the strong binders of this allel
    # index=netMHC_output.find(f"Allele {hla_type}. Number of high binders")
    # if index != -1 and netMHC_output[index+41].isnumeric():
    #     return int(netMHC_output[index+41]) # return the value of string binders
    # else:
    #     print(f"Error in getting 'Strong Binders' value from {hla_type}. Setting to (-1).")
    #     return -1
    
    # find the strong binders value by using regex
    pattern=f"Allele {hla_type}. Number of high binders " + r"(\d+)."
    matche = re.search(pattern, netMHC_output)
    # check if the regex captured the strong binders value
    if matche != None and matche.group(1).isnumeric():
        return int(matche.group(1))
    else:
        print(f"Warning: Could not find 'Strong Binders' value from {hla_type}. Setting to zero (0).")
        return 0
 
def get_GeneSymbol_transcriptID(path):
    # this function gets a path and extract the Gene Symbol and TranscriptID from it
    # define regex
    pattern = r"/(\w+\.\d+\.*\d*\.*\d*_(\w+-*\w+.*\w+))/(\d+)/(\w+\.\d+\.*\d*\.*\d*)/*$"
    # capture the occurances
    match = re.search(pattern, path)
    if match:
        geneSymbol = match.group(2)
        transcriptID = match.group(4)
        return geneSymbol, transcriptID
    else:
        print(f"Error in finding GeneSymbol and TranscriptID in path {path}")
        return None, None

# calculate the difference of SB values of each HLA allele between the two dictionaries
def calculate_SB_difference(sb_dict_1, sb_dict_2):
    sb_difference_dict={}
    # the difference will be between the first dict to the second dict.
    # check if first dict if the normal group. if yes, sign will be 1, else - sign will be (-1) to change the difference value respectivly 
    if sb_dict_1['Group'] == user_args.lable1:
        sign=1
    elif sb_dict_1['Group'] == user_args.lable2:
        sign=-1
    for key in sb_dict_1.keys():
        if key.startswith("HLA"):
            sb_difference_dict[key] = (sb_dict_2[key]-sb_dict_1[key])*sign
    return sb_difference_dict

# run the netMHC command
def run_netmhc(fasta_file, netMHC_dir, form, group):
    # create file name for netMHC output
    if form.lower() == "exclusion":
        output_file = f"netMHC_exclusionAA_{group}.xls"
    elif form.lower() == "inclusion":
        output_file = f"netMHC_inclusionAA_{group}.xls"
    # create output file path
    output_path = os.path.join(netMHC_dir,output_file)
    # create shell command for netMHC
    command =["/private/common/Software/netMHC/netMHC-4.0/Linux_x86_64/bin/netMHC",
              "-hlalist", "/private/common/Software/netMHC/netMHC-4.0/data/allelelist",
              "-syn", "/private/common/Software/netMHC/netMHC-4.0/Linux_x86_64/data/synlists/%s.synlist",
              "-thrfmt", "/private/common/Software/netMHC/netMHC-4.0/threshold/%s.thr",
              "-rdir", "/private/common/Software/netMHC/netMHC-4.0/Linux_x86_64",
              "-version", "/private/common/Software/netMHC/netMHC-4.0/Linux_x86_64/data/version ",
              "-rth", user_args.rank,
              "-a", "HLA-A0101,HLA-A0201,HLA-A0301,HLA-A1101,HLA-A2402,HLA-A2601,HLA-B0702,HLA-B0801,HLA-B1501,HLA-B2705,HLA-B3901,HLA-B4001,HLA-B5801",
              "-l", "8,9,10,11",
              "-f",fasta_file,
              "-tdir", netMHC_dir,
              "-xls", 
              "-xlsfile", output_path]
    # run the command
    #print('Net-MHC command to run:', ' '.join(map(str, [str(item) if isinstance(item, float) else item for item in command])))
    try:
        #print('Net-MHC ommand to run:', ' '.join(command))
        output = subprocess.check_output(command, universal_newlines=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print("Error in running netMHC on subprocess. Exit.")
        print(e.returncode)
        print(e.output)
    # save the log of the proccess
    netMHC_log = os.path.join(netMHC_dir, "netMHC_Log.txt")
    with open (netMHC_log, 'w') as netMHClog:
        netMHClog.write(output)
    # set dictionary of HLA types and their string binders that were found
    sb_HLA_dict = {"HLA-A0101":0, "HLA-A0201":0, "HLA-A0301":0,"HLA-A1101":0, "HLA-A2402":0, "HLA-A2601":0, "HLA-B0702":0, "HLA-B0801":0, "HLA-B1501":0, "HLA-B2705":0,"HLA-B3901":0, "HLA-B4001":0, "HLA-B5801":0}
    # find how many strong binders were found for each HLA type
    for hla_type in sb_HLA_dict.keys():
        sb_HLA_dict[hla_type]=get_HLA_strongBinders(output,hla_type)
    return sb_HLA_dict # return dictionary with updates strong binders value for each allele

# save results in a csv file
def save_results_to_csv(list_of_dicts, output_dir, filename):
    #dicts_list=[dict_sb_1,dict_sb_2,difference_dict]
    output_file = os.path.join(output_dir,filename)
    with open (output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile,fieldnames=list_of_dicts[0].keys())
        writer.writeheader()
        writer.writerows(list_of_dicts)

def get_dicts_from_exist_file(file_path):
    with open(file_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = list(reader)
        for row in rows:
            # if row['Group'] == user_args.lable2+"-"+user_args.lable1:
            #     return row
            # if "-"+user_args.lable1 in row['Group']:
            #     return row
            if row['Group'] == user_args.lable1:
                group1_dict = dict(row)
            elif row['Group'] == user_args.lable2:
                group2_dict = dict(row)
            else:
                return None, None
    return group1_dict, group2_dict

# run the analyze steps on the current transcript directory
def runAnalyze(transcript_dir):
    # get GeneSymbol and TranscriptID of current path
    if user_args.psi_sigma:
        as_type = transcript_dir.strip("/").split("/")[-2]
        transcriptID = transcript_dir.strip("/").split("/")[-3]
        geneSymbol = transcript_dir.strip("/").split("/")[-4] 
    else:
        geneSymbol, transcriptID = get_GeneSymbol_transcriptID(transcript_dir)
        as_type = user_args.as_type
    if geneSymbol == None or transcriptID == None:
        return
    # create 'netMHC' results directory
    netMHC_dir = os.path.join(transcript_dir,f"netMHC_Rank{user_args.rank}")
    # check if netMHC was already ran on this directory
    if os.path.isdir(netMHC_dir):
        exist_files = glob.glob(os.path.join(netMHC_dir, '*.xls'))
        if len(exist_files) == 2:
            print(f"NetMHC was already ran on {netMHC_dir}. Skipping.")
            return
    if not os.path.isdir(netMHC_dir):
        os.mkdir(netMHC_dir)
    # Get all files in the directory
    files = get_absolute_file_paths(transcript_dir)
    # Search for group1 AA sequence fasta file
    if user_args.psi_sigma:
        group1Seq_file = next((os.path.abspath(file) for file in files if file.endswith(f'_{user_args.lable1}.fasta')), None)
    else:
        group1Seq_file = next((os.path.abspath(file) for file in files if file.endswith('.fasta') and f'AA_{user_args.lable1}' in file), None)
    if group1Seq_file is None:
        print(f"Error in locating {user_args.lable1} fasta file in {transcript_dir}. Skipping.")
        return
    groups=getGroupsForms(group1Seq_file) # define the forms (inclusion/exclusion) of the groups (lable1/lable2) 
    if groups is None:
        print(f"Error in define forms of groups in {transcript_dir}. Skipping.")
        return
    # run netMHC command on group1 sequence
    group1_seq_sb_dict = run_netmhc(group1Seq_file, netMHC_dir,groups['1'], user_args.lable1)
    # add GeneSymbol,TranscriptID and Form keys at the beggining of the dictionary
    group1_dict = {"GeneSymbol": geneSymbol, "TranscriptID": transcriptID, "Group": user_args.lable1,"SplicingType":as_type, "Form": groups['1'], "Rank": user_args.rank}
    group1_dict.update(group1_seq_sb_dict)
    # Search for group2 AA sequence fasta file
    if user_args.psi_sigma:
        group2Seq_file = next((os.path.abspath(file) for file in files if file.endswith(f'_{user_args.lable2}.fasta')), None)
    else:
        group2Seq_file = next((os.path.abspath(file) for file in files if file.endswith('.fasta') and f'AA_{user_args.lable2}' in file), None)
    # run netmHC command on group2 sequence
    group2_seq_sb_dict = run_netmhc(group2Seq_file, netMHC_dir, groups['2'], user_args.lable2)
    # add GeneSymbol,TranscriptID and Form keys at the beggining of the dictionary
    group2_dict = {"GeneSymbol": geneSymbol, "TranscriptID": transcriptID, "Group": user_args.lable2, "SplicingType":as_type, "Form": groups['2'], "Rank": user_args.rank}
    group2_dict.update(group2_seq_sb_dict)

    return group1_dict, group2_dict


# global args
list_of_dicts_group1 = [] 
list_of_dicts_group2 = []
list_of_dicts = [] 

if __name__ == '__main__':
    # get absolute paths of transcripts directories
    pathes = getPaths(user_args.input_dir)
    # run netMHC in parallel
    pool = multiprocessing.Pool(processes=10) 
    results = pool.map(runAnalyze, pathes)
    pool.close()
    pool.join()
    print("Done proccessing netMHC on samples.")
