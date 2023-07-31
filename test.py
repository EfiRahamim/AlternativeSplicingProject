import biolib
biolib.utils.STREAM_STDOUT = True
membranefold = biolib.load('KU/MembraneFold')
membranefold_job = membranefold.cli(args='--file /private5/Projects/Efi/AS/TCGA-BRCA/gencode_v36/pipeline_results/ENSG00000139880.19_CDH24/100472/ENST00000397359.7/ENST00000397359.7_splicedAA_Normal.fasta', machine='local') # Blocks until done
membranefold_job.save_files('/private5/Projects/Efi/AS/TCGA-BRCA/gencode_v36/pipeline_results/ENSG00000139880.19_CDH24/100472/MembraneFold/') # Saves all results to `folder_name` dir
