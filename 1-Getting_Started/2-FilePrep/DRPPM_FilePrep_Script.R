####----PATH - SURVEYORS File Preparation----####
# author: Alyssa Obermayer (alyssa.obermayer@moffitt.org)

# Documentation
#' R script to ensure user files are in proper format
#'
#' @param Expression_File - A file path in string format to a users expression matrix dataframe.
#' @param Expression_Delim - A string denoting the delimiter in regular expression of the expression matrix.
#' @param Meta_File - A file path in string format to a users meta/clinical dataframe.
#' @param Meta_Delim - A string denoting the delimiter in regular expression of the meta/clinical dataframe.
#' @param Survival_Time_Columns - A singular string or vector of strings of the column names in the meta/clinical data that represent survival time.
#' @param Survival_Time_Units - A string of either "Days", "Months", or "Years" denoting the unit of time the survival column is in.
#' @param Survival_ID_Columns - A singular string or vector of strings of the column names in the meta/clinical data that represent survival event ID.
#' @param Outfile_Path - A string denoting the file path where the user would like the output files to be when completed
#' @returns A properly formatted expression, clinical, and clinical parameter files for use in the PATH-SURVEYORS suite of tools.


####---- User Input ----####

Expression_File <- ""
Expression_Delim <- "\t"

Meta_File <- ""
Meta_Delim <- "\t"

Survival_Time_Columns <- c("OS.time","EFS.time")
Survival_Time_Units <- "Days"
Survival_ID_Columns <- c("OS","EFS")

Outfile_Path <- "~/R/"





####---- Load Packages ----####

library(readr)
library(dplyr)
library(tools)
library(stringr)


####---- Read in Files ----####

expr <- as.data.frame(read_delim(Expression_File, delim = Expression_Delim, col_names = T))
meta <- as.data.frame(read_delim(Meta_File, delim = Meta_Delim, col_names = T))


####---- Check Outfile Path ----####

if (dir.exists(Outfile_Path)) {
  # Add a forward slash if missing from the end of path
  last_char <- str_sub(Outfile_Path,-1,-1)
  if (last_char != "/") {
    Outfile_Path <- paste(Outfile_Path,"/",sep = "")
  }
}
if (!dir.exists(Outfile_Path)) {
  Outfile_Path <- getwd()
  Outfile_Path <- paste(Outfile_Path,"/",sep = "")
}


####---- Match Sample Names ----####

## Clean Names
colnames(expr) <- gsub("[[:punct:]]",".",colnames(expr))
meta[,1] <- gsub("[[:punct:]]",".",meta[,1])

## Check sample names match
exprSamp <- colnames(expr)[-1]
metaSamp <- meta[,1]
SampSame <- intersect(exprSamp,metaSamp)

## Sample Text Alert
if (length(SampSame) != length(colnames(expr)[-1]) | length(SampSame) != length(rownames(meta))) {
  SampCheck_line1 <- print("Uneven number of samples or mismatching names between expression and meta file.")
  SampCheck_line2 <- print(paste("Number of matching sample names:",length(SampSame)))
  SampCheck_line2 <- print("New files based on similar samples")
  SampCheck_lines <- paste(SampCheck_line1,SampCheck_line2,SampCheck_line3,sep="\n")
} else {
  SampCheck_lines <- "All sample names match between expression and meta data."
}
print(SampCheck_lines)

GeneCol <- colnames(expr)[1]
expr <- expr[,c(GeneCol,SampSame)]
meta <- meta[which(meta[,1] %in% SampSame),]


####---- Clean Expression Data ----####

## Check that expression columns are numeric
expr2 <- expr
isChar <- unname(which(sapply(expr2, function(x) is.character(x))))
isChar <-  isChar[-1]
if (length(isChar) > 0) {
  expr2[isChar] <- sapply(expr2[isChar],as.numeric)
}

## Check for duplicated genes
expr3 <- expr2
colnames(expr3)[1] <- "Symbol"
if (TRUE %in% duplicated(expr3[,1])) {
  DupGene_Line <- "Duplicated genes found. Summarizing to gene with highest average expression."
  expr3 <- expr3 %>%
    group_by(Symbol) %>%
    summarise_all(max)
} else {
  DupGene_Line <- "No duplicated genes found."
}
print(DupGene_Line)
colnames(expr3)[1] <- GeneCol

## Filter out rows with all 0
expr_0 <- length(which(rowSums(expr3[,c(2:ncol(expr3))], na.rm = T) == 0))
if (expr_0 > 0) {
  expr4 <- expr3[which(rowSums(expr3[,c(2:ncol(expr3))], na.rm = T) > 0),]
  Expr0_Line <- paste0(expr_0," genes found expression level of zero across all samples.")
} else {
  Expr0_Line <- "All genes have an average expression level above zero."
}
print(Expr0_Line)


####---- Clean Meta Data ----####

## Check survival time units
if (Survival_Time_Units != "Days") {
  TimeUnit_Line <- "Survival time column(s) converted to days."
  print(TimeUnit_Line)
  if (Survival_Time_Units == "Months") {
    for (i in Survival_Time_Columns) {
      meta[,i] <- meta[,i] * 30.4375
    }
  }
  if (Survival_Time_Units == "Years") {
    for (i in Survival_Time_Columns) {
      meta[,i] <- meta[,i] * 365.25
    }
  }
}

## Remove column white space
colnames(meta) <- gsub(" ","_",colnames(meta))


####---- Generate Meta Parameter File ----####

metaCols <- colnames(meta)
SampNameCol <- colnames(meta)[1]
if (length(Survival_Time_Columns) > 0 && length(Survival_ID_Columns) > 0) {
  metacol_feature <- metaCols[-1]
  metacol_feature <- metacol_feature[!metacol_feature %in% c(Survival_Time_Columns,Survival_ID_Columns)]
  MetaParam1 <- data.frame(Clinical_Column_Name = SampNameCol,
                           Clinical_Column_Type = "SampleName")
  MetaParam2 <- data.frame(Clinical_Column_Name = Survival_Time_Columns,
                           Clinical_Column_Type = "SurvivalTime")
  MetaParam3 <- data.frame(Clinical_Column_Name = Survival_ID_Columns,
                           Clinical_Column_Type = "SurvivalID")
  MetaParam4 <- data.frame(Clinical_Column_Name = metacol_feature,
                           Clinical_Column_Type = "Feature")
  MetaParam <- rbind(MetaParam1,MetaParam2,MetaParam3,MetaParam4)
}



####---- Write reformatted data to files ----####

## Expression
expr_pathfile <- tools::file_path_sans_ext(Expression_File)
expr_fileName <- paste0(Outfile_Path,basename(expr_pathfile),"_Cleaned.txt")
write_delim(expr4,expr_fileName, delim = '\t')

## Meta
meta_pathfile <- tools::file_path_sans_ext(Meta_File)
meta_fileName <- paste0(Outfile_Path,basename(meta_pathfile),"_Cleaned.txt")
write_delim(meta,meta_fileName, delim = '\t')

## Meta Parameters
meta_pathfile <- tools::file_path_sans_ext(Meta_File)
param_fileName <- paste0(Outfile_Path,basename(meta_pathfile),"_Parameters.txt")
write_delim(MetaParam,param_fileName, delim = '\t', col_names = F)


print("~~~ Complete ~~~")




