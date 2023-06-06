

####------------------------------------------------####
#                                                      #
#  This is a script to easily generate an RData list   #
#                                                      #
####------------------------------------------------####



####----User Input----####

#will accept .gmt, .tsv, or .txt
GeneSet_file <- '~/R/data/GeneSetData/Misc/ERstress_Genes_v2.txt'
header.gs <- TRUE

#must use .RData extension
OutFile_PathAndName <- '~/R/data/GeneSetData/Misc/ERstress_Genes.RData'






bioCpacks <- "clusterProfiler"
installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
if (any(installed_packages_BIOC == FALSE)) {
  BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
}
invisible(lapply(bioCpacks, library, character.only = TRUE))

####----Read File----####

if (file_ext(GeneSet_file) == "gmt") {
  GeneSet <- read.gmt(GeneSet_file)
}
if (file_ext(GeneSet_file) == "tsv" || file_ext(GeneSet_file) == "txt") {
  GeneSet <- read.delim(GeneSet_file, header = header.gs, sep = '\t')
}

colnames(GeneSet) <- c("term","gene")


####----Generate RData List----####

gsDataList <- list()
for (i in unique(GeneSet[,1])){
  gsDataList[[i]] <- GeneSet[GeneSet[,1] == i,]$gene
}

save(gsDataList, file = OutFile_PathAndName)











