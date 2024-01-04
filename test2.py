import biolib, time
biolib.utils.STREAM_STDOUT = True # Stream progress from app in real time
aa_path = "/private10/Projects/Efi/CRG/GBM/PSI-Sigma/GencodeGTF/TM_Genes/SplicingEventsFiles/AGPAT5/ENST00000285518.10/MES/151/Inclusion_DMSO.fasta"
out_dir = "/private10/Projects/Efi/CRG/GBM/PSI-Sigma/GencodeGTF/TM_Genes/SplicingEventsFiles/AGPAT5/ENST00000285518.10/MES/151/deepTMHMM/"
cmd = '--fasta '+ aa_path
start_time =time.time()
print("Loading biolib job...")
deeptmhmm = biolib.load('DTU/DeepTMHMM')
end_load_time = time.time()
print(f"Biolib job loaded succesfully! {end_load_time-start_time:.2f} sec.")
deeptmhmm_job = deeptmhmm.cli(args=cmd, machine='local') # Blocks until done
deeptmhmm_job.save_files(out_dir) # Saves all results to `result` dir
