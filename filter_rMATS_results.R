library(argparse)

# CLI Arguments
parser <- ArgumentParser(description='Filter rMATS initial results.\n The filter is according to: P_value, FDR, TransMembrane Genes and IncLevelDifference value')
parser$add_argument("-i", action="store", dest="input", help="Directory of rMATS result files ({SE/A5SS/A3SS/MXE/RI}.MATS.JC.txt)")
parser$add_argument("-o", action="store", dest="output", help="Output directory to write the filtered table. Default: working directory.", default = getwd() )
#parser$add_argument("-type", action="store", dest="type", help="Type of Splicing Event - SE/A5SS/A3SS/MXE/RI.", choices=c("SE","A3SS","A5SS","MXE","RI"))
parser$add_argument("--filter_by_TM", action="store_true", dest="filter_tm", help="Filter results by Trans-Membrane (TM) proteins")
parser$add_argument("-tm", action="store", dest="tm", help="Path of TransMembrane Domains table (UniProt). Defualt: /private10/Projects/Efi/General/transmembrane_Aug23.csv", default="/private10/Projects/Efi/General/transmembrane_Aug23.csv")
parser$add_argument("-pvalue", action="store", dest="pvalue", help="P-Value threshold. Defualt: 0.05", type="double", default=0.05)
parser$add_argument("-fdr", action="store", dest="fdr", help="FDR threshold. Default: 0.05", type="double", default=0.05)
parser$add_argument("-ild", action="store", dest="ild", help="Absolute Inclusion Level Difference threshold value. Default: 0.1", type="double", default=0.1)
user_args <- parser$parse_args()
stopifnot(!is.null(user_args$input)) # halt execution if input file not giving
print(user_args)

result_files <- list.files(user_args$input, pattern = ".MATS.JC.txt", full.names = T)
print(result_files)
for (file in result_files){
    raw_results <- read.csv(file, sep = '\t')
    filtered_results <- subset(raw_results, PValue != 0 & PValue < user_args$pvalue & FDR != 0 & FDR < user_args$fdr & abs(IncLevelDifference) >= user_args$ild)
    print(file)
    print(paste0("Filtering by P-Value ,FDR and IncLevelDifference: ",nrow(filtered_results), " results."))
    # filter by TM proteins
    if (user_args$filter_tm){
        # load TM table
        tm <- read.csv(user_args$tm)
        # merge AS results with TM table according to Gene Synbol
        filtered_results <- filtered_results[filtered_results$geneSymbol %in% tm$GeneSymbol,]
        print(paste0("Merging with TM table: ", nrow(filtered_results), " results."))
    }
    # define splicing type and change 0base to 1base accordingly
    as_type = NULL
    if (grepl("SE", file, ignore.case = TRUE)){
        as_type = "SE"
        filtered_results$exonStart_1base <- filtered_results$exonStart_0base+1
    } else if (grepl("A3SS|A5SS", file, ignore.case = TRUE)){
        filtered_results$longExonStart_1base <- filtered_results$longExonStart_0base+1
        filtered_results$shortES_1base <- filtered_results$shortES+1
        filtered_results$shortEE_1base <- filtered_results$shortEE+1
        if (grepl("A3SS", file, ignore.case = TRUE)){as_type="A3SS"}
        else {as_type="A5SS"}
    } else if (grepl("MXE", file, ignore.case = TRUE)){
        as_type = "MXE"
        filtered_results$X1stExonStart_1base <- filtered_results$X1stExonStart_0base+1
        filtered_results$X2ndExonStart_1base <- filtered_results$X2ndExonStart_0base+1
    } else if (grepl("RI", file, ignore.case = TRUE)){
        as_type = "RI"
        filtered_results$riExonStart_1base <- filtered_results$riExonStart_0base+1
        filtered_results$downstreamES_1base <- filtered_results$downstreamES +1
        filtered_results$upstreamEE_1base <- filtered_results$upstreamEE + 1
    }

    if (user_args$filter_tm){
        output_file_path = paste0(user_args$output,"/",as_type,"_filtered_byTM","_ild",user_args$ild,"_pv",user_args$pvalue,"_fdr",user_args$fdr,".csv")
    } else {
        output_file_path = paste0(user_args$output,"/",as_type,"_filtered","_ild",user_args$ild,"_pv",user_args$pvalue,"_fdr",user_args$fdr,".csv")
    }
    print(paste("Writing file:", output_file_path))
    write.csv(filtered_results, output_file_path, row.names=FALSE)
}



# # --------OLD VERSION---------
# # load SE results (raw results)
# # as <- read.csv("/private5/Projects/Efi/AS/TCGA-BRCA/gencode_v36/significant_results/SE/SE.MATS.JCEC_significant.txt", sep = '\t')
# raw_SE <- read.csv(user_args$input, sep = '\t')
# # filter p_value & FDR < 0.05 & |IncLevelDifference| >= ild
# filtered_SE <- subset(raw_SE, PValue != 0 & PValue < user_args$pvalue & FDR != 0 & FDR < user_args$fdr & abs(IncLevelDifference) >= user_args$ild)
# print(paste0("Filtering by P-Value ,FDR and IncLevelDifference: ",nrow(filtered_SE), " results."))
# 
# # filter by TM proteins
# if (user_args$filter_tm){
#     # load TM table
#     tm <- read.csv(user_args$tm)
#     # merge AS results with TM table according to Gene Synbol
#     filtered_SE <- filtered_SE[filtered_SE$geneSymbol %in% tm$GeneSymbol,]
#     print(paste0("Merging with TM table: ", nrow(filtered_SE), " results."))
# }
# 
# # match UniProtID to each Gene Symbol
# # filtered_SE$uniProtID <- NA
# # for (i in seq_len(nrow(filtered_SE))) {
# #   filtered_SE$uniProtID[i] <- tm$uniProtId[tm$GeneSymbol == filtered_SE$geneSymbol[i]]
# # }
# 
# # change 0base to 1base
# if (user_args$type == "SE"){
#     filtered_SE$exonStart_1base <- filtered_SE$exonStart_0base+1
# } else if (user_args$type %in% c("A3SS","A5SS")){
#     filtered_SE$longExonStart_1base <- filtered_SE$longExonStart_0base+1
#     filtered_SE$shortES_1base <- filtered_SE$shortES+1
#     filtered_SE$shortEE_1base <- filtered_SE$shortEE+1
# } else if (user_args$type == "MXE"){
#     filtered_SE$X1stExonStart_1base <- filtered_SE$X1stExonStart_0base+1
#     filtered_SE$X2ndExonStart_1base <- filtered_SE$X2ndExonStart_0base+1
# } else if (user_args$type == "RI"){
#     filtered_SE$riExonStart_1base <- filtered_SE$riExonStart_0base+1
#     filtered_SE$downstreamES_1base <- filtered_SE$downstreamES +1
#     filtered_SE$upstreamEE_1base <- filtered_SE$upstreamEE + 1
# }
# 
# 
# # filter out exons not divided by 3 (according to zero based Coordinate Systems: https://www.biostars.org/p/84686/)
# #filtered_SE$exonDevideBy3 <- ((filtered_SE$exonEnd-filtered_SE$exonStart_0base)%%3)
# #filtered_SE_significant <- subset(filtered_SE, exonDevideBy3 == 0)
# #print(paste0("Filtering by exon length: ", nrow(filtered_SE_significant), " results."))
# filtered_SE_significant <- filtered_SE
# # write results as csv table
# if (user_args$filter_tm){
#     output_file_path = paste0(user_args$output,"/",user_args$type,"_filtered_byTM.csv")
# } else {
#     output_file_path = paste0(user_args$output,"/",user_args$type,"_filtered.csv")
# }
# #write.csv(filtered_SE_significant, output_file_path, row.names=FALSE)
# 
