import sys, os, glob, multiprocessing, argparse, subprocess, re, csv

# CLI arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Running netMHC on inclusion & exclusion sequences of each transcript. Results will be saved in the directory of each transcript from the giving output directory. Using the 'netMHC' tool located in: /private/common/Software/netMHC-4.0/Linux_x86_64/bin/netMHC. HLA allels that are being used for this analysis: HLA-A0101,HLA-A0201,HLA-A0301,HLA-A2402,HLA-A2601,HLA-B0702,HLA-B0801,HLA-B1501,HLA-B2705,HLA-B3901,HLA-B4001,HLA-B5801. The output of this script is the amount of strong binders that were found in each HLA allel for each group.")
parser.add_argument("-i", action='store', dest='input_dir', required=True, help="Input directory of genes directories")
parser.add_argument("-l1", action='store', dest='lable1', required=True, help="Lable of first group (same order as in the rMATS analysis).")
parser.add_argument("-l2", action='store', dest='lable2', required=True, help="Lable of second group (same order as in the rMATS analysis).")
user_args = parser.parse_args()


# get list of absolut pathes of transcripts directories
def getPaths(input_dir):
    #directories_pattern = input_dir+"*/*/*"
    directories_pattern = os.path.join(input_dir,"*","*","*")
    directories = glob.glob(directories_pattern)
    pathes = [os.path.abspath(dir) for dir in directories if os.path.isdir(dir)] #and dir.find("sashimiplots") == -1]
    if pathes is None:
        print("Error in getting pathes of transcripts.")
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
    if "inclusion" in group1Seq_file:
        groups['1'] = "inclusion"
        groups['2'] = "exclusion"
    elif "exclusion" in group1Seq_file:
        groups['1'] = "exclusion"
        groups['2'] = "inclusion"
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
        print(f"Error in getting 'Strong Binders' value from {hla_type}. Setting to (-1).")
        return -1
 
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
    if form == "exclusion":
        output_file = f"netMHC_exclusionAA_{group}.xls"
    elif form == "inclusion":
        output_file = f"netMHC_inclusionAA_{group}.xls"
    # create output file path
    output_path = os.path.join(netMHC_dir,output_file)
    # create shell command for netMHC
    command =["/private/common/Software/netMHC-4.0/Linux_x86_64/bin/netMHC",
              "-hlalist", "/private/common/Software/netMHC-4.0/data/allelelist",
              "-syn", "/private/common/Software/netMHC-4.0/Linux_x86_64/data/synlists/%s.synlist",
              "-thrfmt", "/private/common/Software/netMHC-4.0/threshold/%s.thr",
              "-rdir", "/private/common/Software/netMHC-4.0/Linux_x86_64",
              "-version", "/private/common/Software/netMHC-4.0/Linux_x86_64/data/version ",
              "-a", "HLA-A0101,HLA-A0201,HLA-A0301,HLA-A2402,HLA-A2601,HLA-B0702,HLA-B0801,HLA-B1501,HLA-B2705,HLA-B3901,HLA-B4001,HLA-B5801",
              "-l", "8,9,10,11",
              "-f",fasta_file,
              "-tdir", netMHC_dir,
              "-xls", 
              "-xlsfile", output_path]
    # run the command
    try:
        output = subprocess.check_output(command, universal_newlines=True)
    except:
        print("Error in running netMHC on subprocess. Exit.")
    # save the log of the proccess
    netMHC_log = os.path.join(netMHC_dir, "netMHC_Log.txt")
    with open (netMHC_log, 'w') as netMHClog:
        netMHClog.write(output)
    # set dictionary of HLA types and their string binders that were found
    sb_HLA_dict = {"HLA-A0101":0, "HLA-A0201":0, "HLA-A0301":0, "HLA-A2402":0, "HLA-A2601":0, "HLA-B0702":0, "HLA-B0801":0, "HLA-B1501":0, "HLA-B2705":0,"HLA-B3901":0, "HLA-B4001":0, "HLA-B5801":0}
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
    # assign lists as global so they can be modified within the function
    #global list_of_dicts_group1, list_of_dicts_group2
    #print(f"Analyzing transcript: {transcript_dir}")
    # get GeneSymbol and TranscriptID of current path
    geneSymbol, transcriptID = get_GeneSymbol_transcriptID(transcript_dir)
    if geneSymbol == None or transcriptID == None:
        return
    # create 'netMHC' results directory
    netMHC_dir = os.path.join(transcript_dir,"netMHC")
    if os.path.isdir(netMHC_dir):
        #print(f"Transcript {transcript_dir} was already checked. Results can be found at {netMHC_dir}.")
        #return
        for file in os.listdir(netMHC_dir):
            if file.endswith(".csv"):
                os.chdir(netMHC_dir)
                group1_dict, group2_dict = get_dicts_from_exist_file(os.path.abspath(file))
                #print(difference_dict)
                os.chdir(user_args.input_dir)
                if group1_dict is None or group2_dict is None:
                    print(f"Could not find group1 and group2 dicts from exist directory {netMHC_dir}")
                    return
                return group1_dict, group2_dict
                #list_of_dicts_group1.append(group1_dict)
                #list_of_dicts_group2.append(group2_dict)
                #print(list_of_dicts_group1, list_of_dicts_group2)
                return
    os.mkdir(netMHC_dir)
    # Get all files in the directory
    #files = os.listdir(transcript_dir)
    files = get_absolute_file_paths(transcript_dir)
    #print(f"Files found in {transcript_dir}:{files}")
    # Search for group1 AA sequence fasta file
    group1Seq_file = next((os.path.abspath(file) for file in files if file.endswith('.fasta') and f'AA_{user_args.lable1}' in file), None)
    groups=getGroupsForms(group1Seq_file) # define the forms (inclusion/exclusion) of the groups (lable1/lable2) 
    if groups is None:
        print(f"Error in define forms of groups in {transcript_dir}. Skipping.")
        return
    # run netMHC command on group1 sequence
    group1_seq_sb_dict = run_netmhc(group1Seq_file, netMHC_dir,groups['1'], user_args.lable1)
    # add GeneSymbol,TranscriptID and Form keys at the beggining of the dictionary
    group1_dict = {"GeneSymbol": geneSymbol, "TranscriptID": transcriptID, "Group": user_args.lable1, "Form": groups['1']}
    group1_dict.update(group1_seq_sb_dict)
    # Search for group2 AA sequence fasta file
    group2Seq_file = next((os.path.abspath(file) for file in files if file.endswith('.fasta') and f'AA_{user_args.lable2}' in file), None)
    # run netmHC command on group2 sequence
    group2_seq_sb_dict = run_netmhc(group2Seq_file, netMHC_dir, groups['2'], user_args.lable2)
    # add GeneSymbol,TranscriptID and Form keys at the beggining of the dictionary
    group2_dict = {"GeneSymbol": geneSymbol, "TranscriptID": transcriptID, "Group": user_args.lable2, "Form": groups['2']}
    group2_dict.update(group2_seq_sb_dict)

    # calculate the difference between SB in each HLA allele
    # differences_dict = calculate_SB_difference(update_inclusion_dict, update_exclusion_dict)
    # update_difference_dict = {"GeneSymbol": geneSymbol, "TranscriptID": transcriptID, "Group":user_args.lable2+"-"+user_args.lable1}
    # update_difference_dict.update(differences_dict)
    # save_results_to_csv([update_inclusion_dict, update_exclusion_dict,update_difference_dict], netMHC_dir, filename=geneSymbol+"_"+transcriptID+"_"+"StrongBinders.csv")
    # return update_difference_dict
    
    # Without calculating the diffrences - just  save the HLA dict's of the groups in one csv file and return the dicts
    save_results_to_csv([group1_dict, group2_dict], netMHC_dir, filename=geneSymbol+"_"+transcriptID+"_"+"StrongBinders.csv")
    return group1_dict, group2_dict
    #list_of_dicts_group1.append(group1_dict)
    #list_of_dicts_group2.append(group2_dict)
    # if groups['inclusion'] == user_args.lable2:
    #     return group1_dict
    # elif groups['exclusion'] == user_args.lable2:
    #     return group2_dict


    #print( update_inclusion_dict,"\n",update_exclusion_dict,"\n",update_difference_dict)


# global args
list_of_dicts_group1 = [] 
list_of_dicts_group2 = []  

if __name__ == '__main__':
    # get absolute paths of transcripts directories
    pathes = getPaths(user_args.input_dir)
    #print(f"Directories that has been found in {user_args.input_dir}:\n{pathes}")
    pool = multiprocessing.Pool(processes=10) 
    #list_of_differences_dict = pool.map(runAnalyze, pathes)
    #list_of_dicts_group1, list_of_dicts_group2 = pool.map(runAnalyze, pathes)
    #pool.map(runAnalyze, pathes)
    results = pool.map(runAnalyze, pathes)
    pool.close()
    pool.join()
    #print(results)
    for group1_dict, group2_dict in results:
        list_of_dicts_group1.append(group1_dict)
        list_of_dicts_group2.append(group2_dict)
    #print(list_of_dicts_group1)
    #print(list_of_dicts_group2)
    #print(list_of_dicts_group1, list_of_dicts_group2)
    #print("list of dicts: ", list_of_differences_dict)
    save_results_to_csv(list_of_dicts_group1, user_args.input_dir, filename=f"{user_args.lable1}_StrongBinders_All.csv")
    save_results_to_csv(list_of_dicts_group2, user_args.input_dir, filename=f"{user_args.lable2}_StrongBinders_All.csv")
    print("Done proccessing netMHC on samples.")
