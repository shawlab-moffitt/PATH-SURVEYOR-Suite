####----PATH - SURVEYOR Immune Deconvolution----####
# author: Alyssa Obermayer (alyssa.obermayer@moffitt.org)

# Documentation
#' R script to calculate immune deconvolution scores for use in PATH-SURVEYOR applications
#'
#' @param ProjectName - A string of the users desired project name.
#' @param Expression_Matrix_File - A file path in string format to a users expression matrix dataframe.
#' @param Meta_Data_File - A file path in string format to a users meta/clinical dataframe.
#' @param Meta_Data_Param_File - (Optional) A file path in string format to a users meta/clinical data parameters dataframe.
#' @param Output_Path - A string denoting the file path where the user would like the output files to be when completed.
#' @param quantiseq - Boolean - TRUE to run the method - FALSE to not run the method
#' @param mcp_counter - Boolean - TRUE to run the method - FALSE to not run the method
#' @param xcell - Boolean - TRUE to run the method - FALSE to not run the method
#' @param epic - Boolean - TRUE to run the method - FALSE to not run the method
#' @param abis - Boolean - TRUE to run the method - FALSE to not run the method
#' @param estimate - Boolean - TRUE to run the method - FALSE to not run the method
#' @param cibersort - (Optional) Boolean - TRUE to run the method - FALSE to not run the method
#' @param cibersort_abs - (Optional) Boolean - TRUE to run the method - FALSE to not run the method
#' @param CIBERSORT_Script - (Optional) "Path/To/CIBERSORT.R"
#' @param LM22_File - (Optional) "Path/To/LM22.txt"
#' @returns An immune deconvolution score matrix, and a meta/clinical data and parameters file with appeneded columns of immune deconvolution scores that can be loaded into the PATHJ-SURVEYORS app.


####----User Input----####

ProjectName <- "PANICI_Melanome_VanAllen"
Expression_Matrix_File <- "Test_Input_Data/Expression_Data.zip"
Clinical_Data_File <- "Test_Input_Data/Clinical_Data.txt"
Clinical_Data_Param_File <- "Test_Input_Data/Clinical_Parameters.txt"
Output_Path <- "Test_Output_Data/"

quantiseq <- TRUE
mcp_counter <- TRUE
xcell <- TRUE
epic <- TRUE
abis <- TRUE
estimate <- TRUE

##--CIBERSORT Option --##
#-Requires License and specific files from https://cibersortx.stanford.edu/csdownload.php
#-Requires CIBERSORT.R and LM22.txt files

cibersort <- FALSE
cibersort_abs <- FALSE

# CIBERSORT.R path and file name
CIBERSORT_Script <- "Path/To/CIBERSORT.R"
# LM22 path and file name
LM22_File <- "Path/To/LM22.txt"





####----Check and Install Packages----####

## Check if Immune deconvolution package is installed
immudecon <- "immunedeconv"
immudecon_check <- immudecon %in% rownames(installed.packages())
if (immudecon_check == TRUE) {
  library(immunedeconv)
}
if (immudecon_check == FALSE) {
  print("Immune Deconvolution package not found, please install")
  print("Please run 'remotes::install_github('omnideconv/immunedeconv')'")
}
packages <- c("dplyr","readr","stringr")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))



####----Read Input----####

##--Expression--##
expr <- as.data.frame(read_delim(Expression_Matrix_File,delim = '\t', col_names = T))
rownames(expr) <- expr[,1]
expr <- expr[,-1]
colnames(expr) <- gsub("[[:punct:]]",".",colnames(expr))

##--Meta--##
if (file.exists(Clinical_Data_File)) {
  meta <- as.data.frame(read_delim(Clinical_Data_File, delim = '\t', col_names = T))
  meta[,1] <- gsub("[[:punct:]]",".",meta[,1])
  
}
colnames(meta)[1] <- "SampleName"
if (file.exists(Clinical_Data_Param_File)) {
  meta_param <- as.data.frame(read_delim(Clinical_Data_Param_File,delim = '\t', col_names = F))
}


##--CIBERSORT Dependecies--##
if (cibersort == TRUE | cibersort_abs == TRUE) {
  if (file.exists(CIBERSORT_Script) & file.exists(LM22_File)) {
    set_cibersort_binary(CIBERSORT_Script)
    set_cibersort_mat(LM22_File)
  }
  else if (!file.exists(CIBERSORT_Script) & !file.exists(LM22_File)) {
    cibersort <- FALSE
    cibersort_abs <- FALSE
  }
}

##--Output Path--##
if (dir.exists(Output_Path)) {
  # Add a forward slash if missing from the end of path
  last_char <- str_sub(Output_Path,-1,-1)
  if (last_char != "/") {
    Output_Path <- paste(Output_Path,"/",sep = "")
  }
}
if (!dir.exists(Output_Path)) {
  Output_Path <- getwd()
  Output_Path <- paste(Output_Path,"/",sep = "")
}

##--Deconvolution Method Parameters--##
deconv_params <- c(quantiseq = quantiseq,cibersort = cibersort,cibersort_abs = cibersort_abs,
                   mcp_counter = mcp_counter,xcell = xcell,epic = epic,abis = abis,estimate = estimate)



####----Immune Deconvolution----####

deconv_res_df <- data.frame(SampleName = colnames(expr))

for (i in 1:length(deconv_params)) {
  if (deconv_params[i] == TRUE) {
    deconv_method <- names(deconv_params)[i]
    deconv_res <- as.data.frame(deconvolute(expr, deconv_method))
    rownames(deconv_res) <- deconv_res[,1]
    deconv_res <- deconv_res[,-1]
    deconv_res <- as.data.frame(t(deconv_res))
    colnames(deconv_res) <- paste(gsub(" ","_",colnames(deconv_res)),deconv_method,"PreProcessedScore",sep = "_")
    deconv_res$SampleName <- rownames(deconv_res)
    deconv_res_df <- merge(deconv_res_df,deconv_res, by = "SampleName")
  }
}



####----Merge With Existing Data----####

##--Meta--##
if (file.exists(Clinical_Data_File)) {
  meta2 <- merge(meta,deconv_res_df, by = "SampleName")
  write_delim(meta2,paste(Output_Path,ProjectName,"_Clinical_ImmDeconv.txt", sep = ""),delim = '\t')
}

##--Meta Param--##
deconv_cols2 <- colnames(deconv_res_df)[c(2:ncol(deconv_res_df))]
deconv_meta_params <- data.frame(deconv_cols2,"Feature")
if (file.exists(Clinical_Data_Param_File)) {
  colnames(deconv_meta_params) <- colnames(meta_param)
  meta_param2 <- rbind(meta_param,deconv_meta_params)
  write_delim(meta_param2,paste(Output_Path,ProjectName,"_Clinical_Params_ImmDeconv.txt", sep = ""),delim = '\t', col_names = F)
}
if (!file.exists(Clinical_Data_Param_File)) {
  meta_param2 <- deconv_meta_params
  write_delim(meta_param2,paste(Output_Path,ProjectName,"_Clinical_Params_ImmDeconv.txt", sep = ""),delim = '\t', col_names = F)
}

##--Immune Deconvolution Scores--##

write_delim(deconv_res_df,paste(Output_Path,ProjectName,"_ImmuneDeconvolution.txt", sep = ""),delim = '\t')




