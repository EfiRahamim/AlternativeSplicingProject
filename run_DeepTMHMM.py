import os,sys, glob, shutil, multiprocessing, subprocess, argparse, biolib
biolib.utils.STREAM_STDOUT = False


# CLI arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Running DeepTMHMM on full & spliced sequences of each transcript. Results will be saved in the directory of each transcript from the giving output directory.")
parser.add_argument("-i", action='store', dest='input_dir', required=True, help="Input directory of genes directories")
parser.add_argument("-l1", action='store', dest='lable1', required=True, help="Lable of first group (same order as in the rMATS analysis)")
parser.add_argument("-l2", action='store', dest='lable2', required=True, help="Lable of second group (same order as in the rMATS analysis)")
parser.add_argument('--PSIsigma', action='store_true', dest='psi_sigma', help='Set for PSI-Sigma tool results.')
user_args = parser.parse_args()

# activate conda environment: pybedtools (located in "/home/alu/rahamie4/anaconda3/envs/pybedtools/")
#subprocess.run(["conda", "activate", "pybedtools"], shell=True, executable="/home/alu/rahamie4/anaconda3/envs/pybedtools/bin/")
# check conda env ws activated correctly
#if not "CONDA_PREFIX" in os.environ or not os.environ["CONDA_PREFIX"]:
#  raise RuntimeError("Conda environemt \'pybedtools\' was not activated correctly")

# global variables
#input_dir = "/private5/Projects/Efi/AS/pipeline_tests/tmp_results/"
true_TM_events = [] # list of true TM splicing events ID
# get list of absolut pathes of transcripts directories
def getPaths(input_dir):
    #directories_pattern = input_dir+"*/*/*"
    if user_args.psi_sigma:
        directories_pattern = os.path.join(input_dir,"*","*","*","*")
    else:
        directories_pattern = os.path.join(input_dir,"*","*","*")
    directories = glob.glob(directories_pattern)
    pathes = [os.path.abspath(dir) for dir in directories if os.path.isdir(dir) and dir.find("sashimiplots") == -1]
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

# run deepTMHMM
def run_DeepTMHMM(deeptmhmm, aa_path, tmhmm_dir, type,spliced=False):
    # create result dir according to gene ID and save results
    if spliced:
        result_dir_path = os.path.join(tmhmm_dir,"Exclusion_"+type)
        result_prefix = 'Exclusion'
    else:
        result_dir_path = os.path.join(tmhmm_dir,"Inclusion_"+type)
        result_prefix = 'Inclusion'
    # create command for the tool
    cmd = '--fasta '+aa_path
    print("Sending command to deeptmhmm.")
    # run DeepTMHMM
    deeptmhmm_job = deeptmhmm.cli(args=cmd, machine='local',result_prefix='my_help_testr' )
    print(f"Saving results in {result_dir_path}")
    deeptmhmm_job.save_files(result_dir_path)
    return result_dir_path

# define type order - Normal/Tumor
def getType(inclusionSeq_file):
    type=[]
    if user_args.lable1 in inclusionSeq_file:
        type = [user_args.lable1, user_args.lable2]
    elif user_args.lable2 in inclusionSeq_file:
        type = [user_args.lable2, user_args.lable1]
    else:
        type=None
    return type

# check if transcript has TM domain
def check_TM(result_dir_path):
    # check if TM domain(s) exist
    result_file_path = os.path.join(result_dir_path,"deeptmhmm_results.md")
    if not os.path.exists(result_file_path):
        print(f'File {result_file_path} does not exist.')
        return
    print("check for TMR's.")
    # open result file
    with open(result_file_path, "r") as result_file:
        for line in result_file:
            if "Number of predicted TMRs: " in line:
                tmr_count_str = line.split("Number of predicted TMRs: ")[-1].strip()
                print(f"found {tmr_count_str} TMR's.")
                break
    # check for TMR's
    if tmr_count_str.isdigit():
        tmr_count = int(tmr_count_str)
        if not tmr_count > 0:
            return False
        else:
            return True

# run the analyze steps on the current transcript directory
def runAnalyze(event_dir):
    global true_TM_events
    print(f"Analyzing event: {event_dir}")
    type=[]
    # create 'deepTMHMM' results directory
    tmhmm_dir = os.path.join(event_dir,"deepTMHMM")
    if os.path.isdir(tmhmm_dir) and len(os.listdir(tmhmm_dir)) != 0:
        print(f"Event {event_dir} was already checked. Results can be found at {tmhmm_dir}.")
        if len(os.listdir(tmhmm_dir)) == 2:
            true_TM_events.append(os.path.basename(event_dir))
            with open ('True_TM_events.txt', 'a') as f:
                f.write(os.path.basename(event_dir)+'\n')
        elif len(os.listdir(tmhmm_dir)) == 1:
            with open('noTMdomain.txt', 'a') as f:
                f.write(event_dir + '\n')
        return
    if not os.path.isdir(tmhmm_dir):
        os.mkdir(tmhmm_dir)
    # Get all files in the directory
    #files = os.listdir(transcript_dir)
    files = get_absolute_file_paths(event_dir)
    print(f"Files found in {event_dir}:{files}")
    # Search for inclusion AA sequence fasta file
    if user_args.psi_sigma:
        inclusionSeq_file = next((os.path.abspath(file) for file in files if file.endswith('.fasta') and 'Inclusion' in file), None)
    else:
        inclusionSeq_file = next((os.path.abspath(file) for file in files if file.endswith('.fasta') and 'inclusionAA' in file), None)
    type=getType(inclusionSeq_file)
    if type is None:
        print(f"Error in define types of {event_dir}. Skipping.")
        return
    # load DeepTMHMM tool
    deeptmhmm = biolib.load('DTU/DeepTMHMM')
    # run deepTMHMM on inclusion sequence
    tmhmm_inclusionSeq_path = run_DeepTMHMM(deeptmhmm,inclusionSeq_file,tmhmm_dir,type[0], spliced=False)
    # check if inclusion transcript has TM domain
    if not check_TM(tmhmm_inclusionSeq_path):
        print(f"No TM domains has been found in {event_dir}.")
        with open('noTMdomain.txt', 'a') as f:
            f.write(event_dir + '\n')
        # delete transcript directory
        #print(f"Deleting directory {event_dir}...")
        #shutil.rmtree(event_dir, onerror=None)
        return
    # Search for spliced AA sequence fasta file
    if user_args.psi_sigma:
        exclusionSeq_file = next((os.path.abspath(file) for file in files if file.endswith('.fasta') and 'Exclusion' in file), None)
    else:
        exclusionSeq_file = next((os.path.abspath(file) for file in files if file.endswith('.fasta') and 'exclusionAA' in file), None)
    # run deepTMHMM on spliced sequence
    tmhmm_splicedSeq_path = run_DeepTMHMM(deeptmhmm,exclusionSeq_file,tmhmm_dir,type[1],spliced=True)
    true_TM_events.append(os.path.basename(event_dir))
    with open ('True_TM_events.txt', 'a') as f:
        f.write(os.path.basename(event_dir)+'\n')




if __name__ == '__main__':
    os.chdir(user_args.input_dir)
    # load DeepTMHMM tool
    #deeptmhmm = biolib.load('DTU/DeepTMHMM')
    # get absolute paths of transcripts directories
    pathes = getPaths(user_args.input_dir)
    print(f"Directories that has been found in {user_args.input_dir}:\n{pathes}")
    pool = multiprocessing.Pool(processes=2) # deepTMHMM is very costly. Not recommanded to run in parallel
    pool.map(runAnalyze, pathes)
    pool.close()
    pool.join()
    print("Done proccessing DeepTMHMM on samples.")
    print(f"Splicing events with true TM domains: {true_TM_events}")
    # with open ('True_TM_events.txt', 'w') as f:
    #     f.write('\n'.join(true_TM_events))




