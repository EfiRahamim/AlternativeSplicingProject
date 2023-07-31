##### Alternative-Transmembrane Exons Project Pipeline

### Step 2: Filtering significant results

# mkdir significant_results
# cd significant_results
# # filter SE.MATS/RI.MATS/A5SS/A3SS outputs
# awk -F"\t" 'NR == 1 || ($19 < 0.05 && $20 < 0.1) {print}' SE.MATS_txt_file > SE.MATS_significant.txt 
# # filter MXE.MATS outputs
# awk -F"\t" 'NR == 1 || ($21 < 0.05 && $22 < 0.1) {print}' MXE.MATS_txt_file > MXE.MATS_significant.txt

# ### Step 3: Creating AA sequence for each gene
# # create list of genes that has significant alternative exon(s) - only the genes ID
# awk 'NR>1 {print $2}' "/private5/Projects/Efi/Carcinoma_HCRN/rmats/output/SE.MATS.JC.significant.txt" | sort -u > /private5/Projects/Efi/Carcinoma_HCRN/rmats/output/significant_genes.txt
# # create list of genes that has significant alternative exon(s) - genes ID and IncLevelAvg value (positive - exon appears more in normal than in tumor. negative - the opposite)
# awk 'BEGIN {OFS="\t"} NR>1 && $23!=0 {print $2, $23}' "/private5/Projects/Efi/Carcinoma_HCRN/rmats/output/SE.MATS.JC.significant.txt" | sort -u > /private5/Projects/Efi/Carcinoma_HCRN/rmats/output/significant_genes.txt
# # generate AA sequence for each gene

SE_OUTPUT=/private5/Projects/Efi/AS/TCGA-BRCA/gencode_v36/output/SE.MATS.JC.txt
OUTPUT_DIR=/private5/Projects/Efi/AS/TCGA-BRCA/gencode_v36/pipeline_results_wide/

Rscript /private5/Projects/Efi/AS/scripts/filter_SE_results.R -i $SE_OUTPUT -o $OUTPUT_DIR

SE_FILTERED="$OUTPUT_DIR/SE_filtered.csv"
python3 /private5/Projects/Efi/AS/scripts/rMATS_to_transcripts.py -i $SE_FILTERED

python3 /private5/Projects/Efi/AS/scripts/analyze_transcripts.py -i $SE_FILTERED -o $OUTPUT_DIR

conda activate pybedtools
python3 /private5/Projects/Efi/AS/scripts/run_DeepTMHMM.py -i $OUTPUT_DIR
