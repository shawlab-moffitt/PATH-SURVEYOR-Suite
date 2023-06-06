

####----Command Line Arguments----####

args=commandArgs(trailingOnly = T)

Param_File <- args[1]


############################ Copyright 2022 Moffitt Cancer Center ############################
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##############################################################################################


####----Packages----####

# Base Packages
packages <- c("gtsummary","tidyr","dplyr","DT","ggpubr","tibble","stringr","tools",
              "survival","readr","survminer","gridExtra","BiocManager","matrixStats")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))
# Bioconductor packages
bioCpacks <- c("GSVA","clusterProfiler")
installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
if (any(installed_packages_BIOC == FALSE)) {
  BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
}
invisible(lapply(bioCpacks, library, character.only = TRUE))


####----Read in Parameter File----####

params <- as.data.frame(read_delim(Param_File,delim = "\t", col_names = F))

## Project Name
Project_Name <- params[which(params[,1] == "Project_Name"),2]
Project_Name <- gsub(" ","_",Project_Name)
Project_Name <- gsub("[[:punct:]]","_",Project_Name)

## Expression File
Expression_Matrix_file <- params[which(params[,1] == "Expression_Matrix_file"),2]

## Meta File
Meta_Data_File <- params[which(params[,1] == "Meta_Data_File"),2]

## Gene Set File
Gene_Set_File <- params[which(params[,1] == "Gene_Set_File"),2]

## Gene Set Name
Gene_Set_Name <- params[which(params[,1] == "Gene_Set_Name"),2]
Gene_Set_Name <- gsub(" ",".",Gene_Set_Name)
Gene_Set_Name <- gsub("[[:punct:]]",".",Gene_Set_Name)
if (!exists("Gene_Set_Name")) {
  Gene_Set_Name <- NA
}

## Output File Path
Output_File_Path <- params[which(params[,1] == "Output_File_Path"),2]
if (!exists("Output_File_Path")) {
  Output_File_Path <- NA
}
if (is.na(Output_File_Path)) {
  Output_File_Path <- getwd()
}
# Add a forward slash if missing from the end of path
last_char <- str_sub(Output_File_Path,-1,-1)
if (last_char != "/") {
  Output_File_Path <- paste(Output_File_Path,"/",sep = "")
}

## Survival Time ID
Survival_Time <- params[which(params[,1] == "Survival_Time_Label"),2]

## Survival ID ID
Survival_ID <- params[which(params[,1] == "Survival_ID_Label"),2]

## Rank By Genes Option
Rank_Genes_Choice <- params[which(params[,1] == "Rank_Genes"),2]

## Covariate Column Name
Covariate_Column_Label <- params[which(params[,1] == "Covariate_Column_Label"),2]

## Covariate Reference
Covariate_Reference <- params[which(params[,1] == "Covariate_Reference"),2]


####----Read in Files----####

print("####----Loading in Files----####")
ST <- Sys.time()
##--Expression--##
expr <- as.data.frame(read_delim(Expression_Matrix_file,delim = '\t', col_names = T))
colnames(expr)[1] <- "Symbol"
if (TRUE %in% duplicated(expr[,1])) {
  expr <- expr %>%
    group_by(Symbol) %>%
    summarise_all(max)
}
isChar <- unname(which(sapply(expr, function(x) is.character(x))))
isChar <-  isChar[-1]
if (length(isChar) > 0) {
  expr[isChar] <- sapply(expr[isChar],as.numeric)
}
expr <- as.data.frame(expr)
rownames(expr) <- expr[,1]
expr <- expr[,-1]
colnames(expr) <- gsub("[[:punct:]]",".",colnames(expr))
expr <- expr[which(rowSums(expr, na.rm = T) > 0),]
expr <- as.matrix(expr)

##--Meta--##
meta <- as.data.frame(read_delim(Meta_Data_File,delim = '\t',col_names = T))
colnames(meta)[1] <- "SampleName"
meta$SampleName <- gsub("[[:punct:]]",".",meta$SampleName)

## Check that all samples have survival data and covariate data if needed
meta <- meta[!is.na(meta[,Survival_Time]),]
meta <- meta[!is.na(meta[,Survival_ID]),]

SampSame <- intersect(colnames(expr),meta$SampleName)
expr <- expr[,SampSame]
meta <- meta[which(meta$SampleName %in% SampSame),]

## Check if user provided covariate
if (!exists("Covariate_Column_Label")) {
  Covariate_Column_Label <- NA
  Covariate_Reference <- NA
}
if (!is.na(Covariate_Column_Label)) {
  meta <- meta[!is.na(meta[,Covariate_Column_Label]),]
}


##--Gene Set--##
# If user loads GMT File
if (!is.na(Gene_Set_File)) {
  
  if (tools::file_ext(Gene_Set_File) == "gmt") {
    GeneSet <- read.gmt(Gene_Set_File)
    colnames(GeneSet) <- c("term","gene")
    GeneSetList <- list()
    for (i in unique(GeneSet[,1])){
      GeneSetList[[i]] <- GeneSet[GeneSet[,1] == i,]$gene
    }
    #GeneSet_list[[Gene_Set_Name]] <- GeneSetList
  }
  # If user loads TSV/TXT file
  if (tools::file_ext(Gene_Set_File) == "tsv" || tools::file_ext(Gene_Set_File) == "txt") {
    GeneSet <- read.delim(Gene_Set_File, header = T, sep = '\t')
    colnames(GeneSet) <- c("term","gene")
    GeneSetList <- list()
    for (i in unique(GeneSet[,1])){
      GeneSetList[[i]] <- GeneSet[GeneSet[,1] == i,]$gene
    }
    #GeneSet_list[[Gene_Set_Name]] <- GeneSetList
  }
  # If user loads RData list file
  if (tools::file_ext(Gene_Set_File) == "RData") {
    loadRData <- function(fileName){
      #loads an RData file, and returns it
      load(fileName)
      get(ls()[ls() != "fileName"])
    }
    GeneSetList <- loadRData(Gene_Set_File)
    #GeneSet_list[[Gene_Set_Name]] <- GeneSetList
  }
  
}



## High-Low function
highlow = function(mat) {
  new_mat = mat;
  new_mat[mat > quantile(as.numeric(mat),na.rm = T)[3]] = "High_AboveMedian";
  new_mat[mat <= quantile(as.numeric(mat),na.rm = T)[3]] = "Low_BelowMedian";
  return (new_mat)
}


####----Pathway Ranking----####

if (Rank_Genes_Choice == "FALSE") {

  ####----ssGSEA----####
  print("####----Performing ssGSEA Analysis----####")
  ## Calculate ssGSEA
  gs_name <- Gene_Set_Name                    #Gene set name
  gs <- GeneSetList                           #Gene Set contents
  ssgsea <- gsva(expr,gs,method = "ssgsea")   #Outputs sample names as column names
  ssgsea_t <- as.data.frame(t(ssgsea))        #Transpose to make column names gs names
  ssgsea_t$SampleName <- rownames(ssgsea_t)   #Change gs names to row names
  ssgsea_t <- ssgsea_t %>%                    #Make rownames first column to output to file
    relocate(SampleName)
  write_delim(as.data.frame(ssgsea_t),paste(Output_File_Path,Project_Name,"_",gs_name,"_ssGSEA_score.txt", sep = ""), delim = '\t')
  
  ## Generate Summary Columns
  ssgsea <- as.data.frame(ssgsea)
  col_to_calc <- ncol(ssgsea)
  ssgsea$Median_Score <- rowMedians(as.matrix(ssgsea[,c(1:col_to_calc)]))
  ssgsea$Mean_Score <- rowMeans(as.matrix(ssgsea[,c(1:col_to_calc)]))
  ssgsea$Min_Score <- apply(ssgsea[,c(1:col_to_calc)], 1, FUN = min)
  ssgsea$Max_Score <- apply(ssgsea[,c(1:col_to_calc)], 1, FUN = max)
  
  ## Extract Summary Columns
  ssgsea_summ <- ssgsea %>%
    select(Median_Score,Mean_Score,Min_Score,Max_Score)
  rownames(ssgsea_summ) <- gsub("[[:punct:]]",".",rownames(ssgsea_summ))
  rownames(ssgsea_summ) <- gsub(" ",".",rownames(ssgsea_summ))
  ssgsea_summ$Feature <- rownames(ssgsea_summ)
  ssgsea_summ <- ssgsea_summ %>%
    relocate(Feature)
  
  
  ####----Perform Median Cutpoint----####
  print("####----Performing Median Cutpoint----####")
  ##--Median Cutpoint of ssGSEA--##
  gs_name <- Gene_Set_Name
  # Place sample names as row names
  rownames(ssgsea_t) <- ssgsea_t[,1]
  ssgsea_t <- ssgsea_t[,-1]
  ## loop through each gene set to get BIN
  ## Perform High/Low function
  new_col_list <- c()
  for (j in colnames(ssgsea_t)) {
    ssgsea_t[,paste(j,'.MedianCutP',sep = "")] <- highlow(ssgsea_t[, which(colnames(ssgsea_t) == j)])
    new_col_list <- c(new_col_list,paste(j,'.MedianCutP',sep = ""))
  }
  ## Reformat
  ssGSEA_BIN <- ssgsea_t[,new_col_list]
  ssGSEA_BIN$SampleName <- rownames(ssGSEA_BIN)
  ssGSEA_BIN <- ssGSEA_BIN %>%
    relocate(SampleName)
  write_delim(ssGSEA_BIN,paste(Output_File_Path,Project_Name,"_",gs_name,"_MedianCutP.txt", sep = ""), delim = '\t')
  
  
  ####----CoxH Ranking----#### 
  print("####----Performing Coxh Analysis----####")
  gs_name <- Gene_Set_Name
  ## Clean Gene Set names
  colnames(ssGSEA_BIN) <- gsub("[[:punct:]]",".",colnames(ssGSEA_BIN))
  colnames(ssGSEA_BIN) <- gsub(" ",".",colnames(ssGSEA_BIN))
  ## Merge with Meta data
  ssGSEA_BIN_meta <- merge(meta,ssGSEA_BIN,by = "SampleName",all.y = T)
  ## Binary columns to loop through
  calc_cols <- grep("..MedianCutP$", colnames(ssGSEA_BIN_meta), value = T)
  
  ##--Generate Headers--##
  ProjName_Line <- paste("## Project Name:",Project_Name)
  Survival_Line <- paste("## Survival/Event Variables:",Survival_Time,"///",Survival_ID)
  Samples_Line <- paste("## Number of Samples:",nrow(meta))
  expr_file <- paste("## Expression Matrix File Name:",basename(Expression_Matrix_file))
  meta_file <- paste("## Meta Data File Name:",basename(Meta_Data_File))
  header <- c("variable","var_label","var_type","reference_row","row_type","header_row",
              "N_obs","N_event","N","coefficients_type","coefficients_label","Criteria","term",
              "var_class","var_nlevels","contrasts","contrasts_type","n_obs","n_event",
              "exposure","Hazard_Ratio","std.error","statistic","nevent","conf.low",
              "conf.high","ci","High_AboveMedian_Pval","Concordance","Likelihood_Ratio_Pval","Wald_Test_Pval",
              "Logrank_Test_Pval","Chisq","df","HazardAssumption_Global_Pval","Median_Score","Mean_Score","Min_Score","Max_Score")
  header_add <- c("variable","var_label","var_type","reference_row","row_type","header_row",
              "N_obs","N_event","N","coefficients_type","coefficients_label","Criteria","term",
              "var_class","var_nlevels","contrasts","contrasts_type","n_obs","n_event",
              "exposure","Hazard_Ratio","std.error","statistic","nevent","conf.low",
              "conf.high","ci","High_AboveMedian_Pval","Concordance","Likelihood_Ratio_Pval","Wald_Test_Pval",
              "Logrank_Test_Pval","Chisq_Global","df_Global","HazardAssumption_Pval_Global","Chisq_Feature","df_Feature","HazardAssumption_Pval_Feature"
              ,"Chisq_covariate","df_Covariate","HazardAssumption_Pval_Covariate","Median_Score","Mean_Score","Min_Score","Max_Score")
  header_int <- c("variable","var_label","var_type","reference_row","row_type","header_row",
                  "N_obs","N_event","N","coefficients_type","coefficients_label","Criteria","term",
                  "var_class","var_nlevels","contrasts","contrasts_type","n_obs","n_event",
                  "exposure","Hazard_Ratio","std.error","statistic","nevent","conf.low",
                  "conf.high","ci","High_AboveMedian_Pval","Concordance","Likelihood_Ratio_Pval","Wald_Test_Pval",
                  "Logrank_Test_Pval","Chisq_Global","df_Global","HazardAssumption_Pval_Global","Chisq_Feature","df_Feature","HazardAssumption_Pval_Feature",
                  "Chisq_Covariate","df_Covariate","HazardAssumption_Pval_Covariate","Chisq_FeatureCovariate","df_FeatureCovariate","HazardAssumption_Pval_FeatureCovariate",
                  "Median_Score","Mean_Score","Min_Score","Max_Score")
  
  
  ####----Univariate CoxH Analysis----####
  ## If no covariate being analyzed
  if (is.na(Covariate_Column_Label)) {
    ## Add Additional Headers
    geneset_file <- paste("## Gene Set File Name:",basename(Gene_Set_File))
    line_geneset <- paste("## Gene Set:",gs_name,"Of",length(GeneSetList),"Gene Sets")
    line_analysis <- paste("## ssGSEA was used to calculate pathway scores and group samples Above or Below median ssGSEA score per gene set")
    line_analysis2 <- paste("## Pathways were ranked by Coxph Likelihood_Ratio_Pval")
    line_analysis3 <- paste("## Hazard Ratios > 1 represent pathway scores associated with High-Risk")
    param_lines <- paste(ProjName_Line,expr_file,meta_file,geneset_file,Survival_Line,Samples_Line,line_geneset,line_analysis,line_analysis2,line_analysis3,sep = "\n")
    file_made <- paste(Output_File_Path,Project_Name,"_",gs_name,"_coxh.txt", sep = "")
    write(param_lines,file = file_made, append = T, sep = '\t', ncolumns = 39)
    
    ##--CoxH Loop--##
    out_df <- as.data.frame(read_delim(file_made,delim = '\t', col_names = F))
    out_df[c(2:39)] <- NA
    out_df <- rbind(out_df,header)
    write(header,file = file_made, append = T, sep = '\t', ncolumns = 39)
    k <- 1
    for (j in calc_cols) {
      meta_ssgsea_sub <- ssGSEA_BIN_meta[,c(Survival_Time,Survival_ID,j)]
      meta_ssgsea_sub[,j] <- as.factor(meta_ssgsea_sub[,j])
      colnames(meta_ssgsea_sub)[which(colnames(meta_ssgsea_sub) == Survival_Time)] <- "time"
      colnames(meta_ssgsea_sub)[which(colnames(meta_ssgsea_sub) == Survival_ID)] <- "ID"
      colnames(meta_ssgsea_sub)[which(colnames(meta_ssgsea_sub) == j)] <- "Feature"
      
      if (length(levels(as.factor(meta_ssgsea_sub[,"Feature"]))) > 1) {
      
        meta_ssgsea_sub[,"Feature"] <- relevel(meta_ssgsea_sub[,"Feature"], ref = "Low_BelowMedian")
        
        tab <- coxph(as.formula(paste("Surv(time,ID) ~ Feature",sep = "")),
                     data = meta_ssgsea_sub)
        temp_tab <- tab %>% 
          gtsummary::tbl_regression(exp = TRUE) %>%
          as_gt()
        temp_tab_df <- as.data.frame(temp_tab)
        ## Get additional Pvals
        out <- capture.output(summary(tab))
        con_line <- grep("^Concordance=",out,value = T)
        con_v <- str_split(con_line," ")[[1]][2]
        lik_line <- grep("^Likelihood ratio test=",out,value = T)
        lik_p <- str_split(lik_line,"=")[[1]][3]
        wal_line <- grep("^Wald test",out,value = T)
        wal_p <- str_split(wal_line,"=")[[1]][3]
        sco_line <- grep("^Score ",out,value = T)
        sco_p <- str_split(sco_line,"=")[[1]][3]
        if (min(table(meta_ssgsea_sub$Feature)) > 1 & !any(grepl("Inf",temp_tab_df$ci)) & temp_tab_df[3,21] != "") {
          zph <- capture.output(cox.zph(tab))
          HA_line <- zph[3]
          chi_n <- str_split(HA_line,"\\s+")[[1]][2]
          dof_n <- str_split(HA_line,"\\s+")[[1]][3]
          pva_n <- str_split(HA_line,"\\s+")[[1]][4]
        } else {
          chi_n <- NA
          dof_n <- NA
          pva_n <- NA
        }
        temp_tab_vect <- as.character(c(temp_tab_df[3,]))
        temp_tab_vect <- c(temp_tab_vect,con_v,lik_p,wal_p,sco_p,chi_n,dof_n,pva_n)
        ## Add summary values
        variable <- file_path_sans_ext(j)
        var_scores <- unname(unlist(c(ssgsea_summ[variable,c(2:5)])))
        temp_tab_vect <- c(temp_tab_vect,var_scores)
        temp_tab_vect[1] <- variable
        temp_tab_vect[2] <- j
        temp_tab_vect[12] <- "High_AboveMedian"
        temp_tab_vect[13] <- paste0(variable,"High_AboveMedian")
        out_df <- rbind(out_df,temp_tab_vect)
        write(temp_tab_vect,file = file_made, append = T, sep = '\t', ncolumns = 39)
        if (k %% 100 == 0) {
          print(paste("Features Processed:",k))
        }
        k <- k+1
        
      }
    }
    
    ##--Cal Adj Pvals Univar--##
    out_df_top <- out_df[grep("##",out_df[,1]),]
    out_df_top[,c(40:42)] <- NA
    tab_df3 <- out_df[grep("##",out_df[,1],invert = T),]
    colnames(tab_df3) <- header
    tab_df3 <- tab_df3[-1,]
    Likelihood_Ratio_adjPval_BH <- p.adjust(as.numeric(tab_df3$Likelihood_Ratio_Pval), method = "BH")
    Wald_Test_adjPval_BH <- p.adjust(as.numeric(tab_df3$Wald_Test_Pval), method = "BH")
    Logrank_Test_adjPval_BH <- p.adjust(as.numeric(tab_df3$Logrank_Test_Pval), method = "BH")
    tab_df3 <- cbind(tab_df3,Likelihood_Ratio_adjPval_BH,Wald_Test_adjPval_BH,Logrank_Test_adjPval_BH)
    tab_df3$Likelihood_Ratio_Pval <- as.numeric(tab_df3$Likelihood_Ratio_Pval)
    tab_df3_ordered <- tab_df3[order(tab_df3$Likelihood_Ratio_Pval, decreasing = F),]
    tab_df3_ordered$Likelihood_Ratio_Pval <- formatC(tab_df3_ordered$Likelihood_Ratio_Pval)
    tab_df3_ordered$Wald_Test_Pval <- formatC(as.numeric(tab_df3_ordered$Wald_Test_Pval))
    tab_df3_ordered$Logrank_Test_Pval <- formatC(as.numeric(tab_df3_ordered$Logrank_Test_Pval))
    ## Reorder Columns and write final table
    tab_df3_ordered <- tab_df3_ordered %>%
      relocate(variable,Hazard_Ratio,ci,High_AboveMedian_Pval,Concordance,Likelihood_Ratio_Pval,Wald_Test_Pval,
               Logrank_Test_Pval,Likelihood_Ratio_adjPval_BH,Wald_Test_adjPval_BH,Logrank_Test_adjPval_BH,HazardAssumption_Global_Pval,
               Chisq,df,Criteria,Median_Score,Mean_Score,Min_Score,Max_Score)
    new_header <- colnames(tab_df3_ordered)
    colnames(tab_df3_ordered) <- colnames(out_df_top)
    tab_df3_ordered <- rbind(out_df_top,new_header,tab_df3_ordered)
    new_file <- gsub(".txt","_ranked.txt",file_made)
    write_delim(tab_df3_ordered, new_file, delim = '\t', col_names = F, na = "")
    print("CoxH Ranking Complete")
  }
  
  ####----Bivariate CoxH Analysis----####
  if (!is.na(Covariate_Column_Label)) {
    
    ##--Generate Covariate Specific Headers--##
    ## Covariate Specific Lines
    line_covar <- paste("## Covariate and Reference Variable:",Covariate_Column_Label,"Referencing",Covariate_Reference)
    covar_brk <- as.data.frame(table(ssGSEA_BIN_meta[,Covariate_Column_Label]))
    covars <- c()
    for (k in 1:nrow(covar_brk)) {
      var <- as.character(covar_brk[k,1])
      num <- as.character(covar_brk[k,2])
      varnum <- c("### Variable:",var,"(","N =",num,")")
      covars <- c(covars,paste(varnum,collapse = " "))
    }
    line_covar2 <- paste("## Covariate Breakdown:\n",paste(covars,collapse = "\n"),sep = "")
    
    ## Geneset Lines
    geneset_file <- paste("## Gene Set File Name:",basename(Gene_Set_File))
    line_geneset <- paste("## Gene Set:",gs_name,"Of",length(GeneSetList),"Gene Sets")
    line_analysis_int <- paste("## ssGSEA was used to calculate pathway scores and group samples Above or Below median ssGSEA score per gene set then ran through interactive Cox Proportional Hazard Analysis with the covariate")
    line_analysis_add <- paste("## ssGSEA was used to calculate pathway scores and group samples Above or Below median ssGSEA score per gene set then ran through additive Cox Proportional Hazard Analysis with the covariate")
    line_analysis2 <- paste("## Pathways were ranked by Coxph Likelihood_Ratio_Pval")
    line_analysis3 <- paste("## Hazard Ratios > 1 represent pathway scores associated with High-Risk")
    param_lines_add <- paste(ProjName_Line,expr_file,meta_file,geneset_file,Survival_Line,Samples_Line,line_geneset,line_analysis_add,line_analysis2,line_analysis3,line_covar,line_covar2,sep = "\n")
    param_lines_int <- paste(ProjName_Line,expr_file,meta_file,geneset_file,Survival_Line,Samples_Line,line_geneset,line_analysis_int,line_analysis2,line_analysis3,line_covar,line_covar2,sep = "\n")
    file_made_add <- paste(Output_File_Path,Project_Name,"_",gs_name,"_",Covariate_Column_Label,"_Additive_coxh.txt", sep = "")
    file_made_int <- paste(Output_File_Path,Project_Name,"_",gs_name,"_",Covariate_Column_Label,"_Interactive_coxh.txt", sep = "")
    write(param_lines_add,file = file_made_add, append = T, sep = '\t', ncolumns = 45)
    write(param_lines_int,file = file_made_int, append = T, sep = '\t', ncolumns = 48)
    
    ##--CoxH Loop--##
    out_df_add <- as.data.frame(read_delim(file_made_add,delim = '\t', col_names = F))
    out_df_add[c(2:45)] <- NA
    out_df_add <- rbind(out_df_add,header_add)
    out_df_int <- as.data.frame(read_delim(file_made_int,delim = '\t', col_names = F))
    out_df_int[c(2:48)] <- NA
    out_df_int <- rbind(out_df_int,header_int)
    write(header,file = file_made_add, append = T, sep = '\t', ncolumns = 45)
    write(header,file = file_made_int, append = T, sep = '\t', ncolumns = 48)
    k <- 1
    for (j in calc_cols) {
      meta_ssgsea_sub <- ssGSEA_BIN_meta[,c(Survival_Time,Survival_ID,j,Covariate_Column_Label)]
      # factor hi/lo column
      meta_ssgsea_sub[,j] <- as.factor(meta_ssgsea_sub[,j])
      meta_ssgsea_sub[,Covariate_Column_Label] <- as.factor(meta_ssgsea_sub[,Covariate_Column_Label])
      colnames(meta_ssgsea_sub)[which(colnames(meta_ssgsea_sub) == Survival_Time)] <- "time"
      colnames(meta_ssgsea_sub)[which(colnames(meta_ssgsea_sub) == Survival_ID)] <- "ID"
      colnames(meta_ssgsea_sub)[which(colnames(meta_ssgsea_sub) == j)] <- "Feature"
      colnames(meta_ssgsea_sub)[which(colnames(meta_ssgsea_sub) == Covariate_Column_Label)] <- "Covariate"
      # factor reference column
      if (is.na(Covariate_Reference)) {
        current_ref <- levels(meta_ssgsea_sub[,"Covariate"])[1]
        print(paste("No Covariate Reference Given. Using ",current_ref," as Coxh reference variable.",sep = ""))
      } else {
        meta_ssgsea_sub[,"Covariate"] <- relevel(meta_ssgsea_sub[,"Covariate"], ref = Covariate_Reference)
      }
      if (length(levels(as.factor(meta_ssgsea_sub[,"Feature"]))) > 1) {
        
        meta_ssgsea_sub[,"Feature"] <- relevel(meta_ssgsea_sub[,"Feature"], ref = "Low_BelowMedian")
        
        ##--Additive--##
        tab_add <- coxph(as.formula(paste("Surv(time,ID) ~ Feature + Covariate",sep = "")),
                         data = meta_ssgsea_sub)
        temp_tab_add <- tab_add %>% 
          gtsummary::tbl_regression(exp = TRUE) %>%
          as_gt()
        temp_tab_add_df <- as.data.frame(temp_tab_add)
        out_add <- capture.output(summary(tab_add))
        con_line_add <- grep("^Concordance=",out_add,value = T)
        con_v_add <- str_split(con_line_add," ")[[1]][2]
        lik_line_add <- grep("^Likelihood ratio test=",out_add,value = T)
        lik_p_add <- str_split(lik_line_add,"=")[[1]][3]
        wal_line_add <- grep("^Wald test",out_add,value = T)
        wal_p_add <- str_split(wal_line_add,"=")[[1]][3]
        sco_line_add <- grep("^Score ",out_add,value = T)
        sco_p_add <- str_split(sco_line_add,"=")[[1]][3]
        if (min(table(meta_ssgsea_sub$Feature)) > 1 & !any(grepl("Inf",temp_tab_df$ci)) & temp_tab_df[3,21] != "") {
          xph_add <- capture.output(cox.zph(tab_add))
          HA_Feat_line_add <- xph_add[2]
          chi_Feat_n_add <- str_split(HA_Feat_line_add,"\\s+")[[1]][2]
          dof_Feat_n_add <- str_split(HA_Feat_line_add,"\\s+")[[1]][3]
          pva_Feat_n_add <- str_split(HA_Feat_line_add,"\\s+")[[1]][4]
          HA_Covar_line_add <- xph_add[3]
          chi_Covar_n_add <- str_split(HA_Covar_line_add,"\\s+")[[1]][2]
          dof_Covar_n_add <- str_split(HA_Covar_line_add,"\\s+")[[1]][3]
          pva_Covar_n_add <- str_split(HA_Covar_line_add,"\\s+")[[1]][4]
          HA_Global_line_add <- xph_add[4]
          chi_Global_n_add <- str_split(HA_Global_line_add,"\\s+")[[1]][2]
          dof_Global_n_add <- str_split(HA_Global_line_add,"\\s+")[[1]][3]
          pva_Global_n_add <- str_split(HA_Global_line_add,"\\s+")[[1]][4]
        } else {
          chi_Feat_n_add <- NA
          dof_Feat_n_add <- NA
          pva_Feat_n_add <- NA
          chi_Covar_n_add <- NA
          dof_Covar_n_add <- NA
          pva_Covar_n_add <- NA
          chi_Global_n_add <- NA
          dof_Global_n_add <- NA
          pva_Global_n_add <- NA
        }
        
        temp_tab_vect_add <- as.character(c(temp_tab_add_df[3,]))
        temp_tab_vect_add <- c(temp_tab_vect_add,con_v_add,lik_p_add,wal_p_add,sco_p_add,
                               chi_Feat_n_add,dof_Feat_n_add,pva_Feat_n_add,
                               chi_Covar_n_add,dof_Covar_n_add,pva_Covar_n_add,
                               chi_Global_n_add,dof_Global_n_add,pva_Global_n_add)
        ## Add summary values
        variable <- file_path_sans_ext(j)
        var_scores <- unname(unlist(c(ssgsea_summ[variable,c(2:5)])))
        temp_tab_vect_add <- c(temp_tab_vect_add,var_scores)
        temp_tab_vect_add[1] <- variable
        temp_tab_vect_add[2] <- j
        temp_tab_vect_add[12] <- "High_AboveMedian"
        temp_tab_vect_add[13] <- paste0(variable,"High_AboveMedian")
        out_df_add <- rbind(out_df_add,temp_tab_vect_add)
        write(temp_tab_vect_add,file = file_made_add, append = T, sep = '\t', ncolumns = 39)
        
        ##--Interactive--##
        tab_int <- coxph(as.formula(paste("Surv(time,ID) ~ Feature * Covariate",sep = "")),
                         data = meta_ssgsea_sub)
        temp_tab_int <- tab_int %>% 
          gtsummary::tbl_regression(exp = TRUE) %>%
          as_gt()
        temp_tab_int_df <- as.data.frame(temp_tab_int)
        out_int <- capture.output(summary(tab_int))
        con_line_int <- grep("^Concordance=",out_int,value = T)
        con_v_int <- str_split(con_line_int," ")[[1]][2]
        lik_line_int <- grep("^Likelihood ratio test=",out_int,value = T)
        lik_p_int <- str_split(lik_line_int,"=")[[1]][3]
        wal_line_int <- grep("^Wald test",out_int,value = T)
        wal_p_int <- str_split(wal_line_int,"=")[[1]][3]
        sco_line_int <- grep("^Score ",out_int,value = T)
        sco_p_int <- str_split(sco_line_int,"=")[[1]][3]
        if (min(table(meta_ssgsea_sub$Feature)) > 1 & !any(grepl("Inf",temp_tab_df$ci)) & temp_tab_df[3,21] != "") {
          xph_int <- capture.output(cox.zph(tab_int))
          HA_Feat_line_int <- xph_int[2]
          chi_Feat_n_int <- str_split(HA_Feat_line_int,"\\s+")[[1]][2]
          dof_Feat_n_int <- str_split(HA_Feat_line_int,"\\s+")[[1]][3]
          pva_Feat_n_int <- str_split(HA_Feat_line_int,"\\s+")[[1]][4]
          HA_Covar_line_int <- xph_int[3]
          chi_Covar_n_int <- str_split(HA_Covar_line_int,"\\s+")[[1]][2]
          dof_Covar_n_int <- str_split(HA_Covar_line_int,"\\s+")[[1]][3]
          pva_Covar_n_int <- str_split(HA_Covar_line_int,"\\s+")[[1]][4]
          HA_FeatCovar_line_int <- xph_int[4]
          chi_FeatCovar_n_int <- str_split(HA_FeatCovar_line_int,"\\s+")[[1]][2]
          dof_FeatCovar_n_int <- str_split(HA_FeatCovar_line_int,"\\s+")[[1]][3]
          pva_FeatCovar_n_int <- str_split(HA_FeatCovar_line_int,"\\s+")[[1]][4]
          HA_Global_line_int <- xph_int[5]
          chi_Global_n_int <- str_split(HA_Global_line_int,"\\s+")[[1]][2]
          dof_Global_n_int <- str_split(HA_Global_line_int,"\\s+")[[1]][3]
          pva_Global_n_int <- str_split(HA_Global_line_int,"\\s+")[[1]][4]
        } else {
          chi_Feat_n_int <- NA
          dof_Feat_n_int <- NA
          pva_Feat_n_int <- NA
          chi_Covar_n_int <- NA
          dof_Covar_n_int <- NA
          pva_Covar_n_int <- NA
          chi_Global_n_int <- NA
          dof_Global_n_int <- NA
          pva_Global_n_int <- NA
        }
        
        temp_tab_vect_int <- as.character(c(temp_tab_int_df[3,]))
        temp_tab_vect_int <- c(temp_tab_vect_int,con_v_int,lik_p_int,wal_p_int,sco_p_int,
                               chi_Feat_n_int,dof_Feat_n_int,pva_Feat_n_int,
                               chi_Covar_n_int,dof_Covar_n_int,pva_Covar_n_int,
                               chi_FeatCovar_n_int,dof_FeatCovar_n_int,pva_FeatCovar_n_int,
                               chi_Global_n_int,dof_Global_n_int,pva_Global_n_int)
        ## Add summary values
        variable <- file_path_sans_ext(j)
        var_scores <- unname(unlist(c(ssgsea_summ[variable,c(2:5)])))
        temp_tab_vect_int <- c(temp_tab_vect_int,var_scores)
        temp_tab_vect_int[1] <- variable
        temp_tab_vect_int[2] <- j
        temp_tab_vect_int[12] <- "High_AboveMedian"
        temp_tab_vect_int[13] <- paste0(variable,"High_AboveMedian")
        out_df_int <- rbind(out_df_int,temp_tab_vect_int)
        write(temp_tab_vect_int,file = file_made_int, append = T, sep = '\t', ncolumns = 48)
        if (k %% 100 == 0) {
          print(paste("Features Processed:",k))
        }
        k <- k+1
      }
    }
        
        
    ##--Calc Adj Pval Add--##
    out_df_add_top <- out_df_add[grep("##",out_df_add[,1]),]
    out_df_add_top[,c(46:48)] <- NA
    tab_df_add <- out_df_add[grep("##",out_df_add[,1],invert = T),]
    colnames(tab_df_add) <- header_add
    tab_df_add <- tab_df_add[-1,]
    Likelihood_Ratio_adjPval_BH <- p.adjust(as.numeric(tab_df_add$Likelihood_Ratio_Pval), method = "BH")
    Wald_Test_adjPval_BH <- p.adjust(as.numeric(tab_df_add$Wald_Test_Pval), method = "BH")
    Logrank_Test_adjPval_BH <- p.adjust(as.numeric(tab_df_add$Logrank_Test_Pval), method = "BH")
    tab_df_add <- cbind(tab_df_add,Likelihood_Ratio_adjPval_BH,Wald_Test_adjPval_BH,Logrank_Test_adjPval_BH)
    tab_df_add$Likelihood_Ratio_Pval <- as.numeric(tab_df_add$Likelihood_Ratio_Pval)
    tab_df_add_ordered <- tab_df_add[order(tab_df_add$Likelihood_Ratio_Pval, decreasing = F),]
    tab_df_add_ordered$Likelihood_Ratio_Pval <- formatC(tab_df_add_ordered$Likelihood_Ratio_Pval)
    tab_df_add_ordered$Wald_Test_Pval <- formatC(as.numeric(tab_df_add_ordered$Wald_Test_Pval))
    tab_df_add_ordered$Logrank_Test_Pval <- formatC(as.numeric(tab_df_add_ordered$Logrank_Test_Pval))
    tab_df_add_ordered <- tab_df_add_ordered %>%
      relocate(variable,Hazard_Ratio,ci,High_AboveMedian_Pval,Concordance,Likelihood_Ratio_Pval,Wald_Test_Pval,
               Logrank_Test_Pval,Likelihood_Ratio_adjPval_BH,Wald_Test_adjPval_BH,Logrank_Test_adjPval_BH,
               HazardAssumption_Pval_Global,HazardAssumption_Pval_Feature,HazardAssumption_Pval_Covariate,Criteria,
               Median_Score,Mean_Score,Min_Score,Max_Score)
    new_header_add <- colnames(tab_df_add_ordered)
    colnames(tab_df_add_ordered) <- colnames(out_df_add_top)
    tab_df_add_ordered <- rbind(out_df_add_top,new_header_add,tab_df_add_ordered)
    new_file_add <- gsub(".txt","_ranked.txt",file_made_add)
    write_delim(tab_df_add_ordered, new_file_add, delim = '\t', col_names = F, na = "")
    
    
    ##--Calc Adj Pval Int--##
    out_df_int_top <- out_df_int[grep("##",out_df_int[,1]),]
    out_df_int_top[,c(49:51)] <- NA
    tab_df_int <- out_df_int[grep("##",out_df_int[,1],invert = T),]
    colnames(tab_df_int) <- header_int
    tab_df_int <- tab_df_int[-1,]
    Likelihood_Ratio_adjPval_BH <- p.adjust(as.numeric(tab_df_int$Likelihood_Ratio_Pval), method = "BH")
    Wald_Test_adjPval_BH <- p.adjust(as.numeric(tab_df_int$Wald_Test_Pval), method = "BH")
    Logrank_Test_adjPval_BH <- p.adjust(as.numeric(tab_df_int$Logrank_Test_Pval), method = "BH")
    tab_df_int <- cbind(tab_df_int,Likelihood_Ratio_adjPval_BH,Wald_Test_adjPval_BH,Logrank_Test_adjPval_BH)
    tab_df_int$Likelihood_Ratio_Pval <- as.numeric(tab_df_int$Likelihood_Ratio_Pval)
    tab_df_int_ordered <- tab_df_int[order(tab_df_int$Likelihood_Ratio_Pval, decreasing = F),]
    tab_df_int_ordered$Likelihood_Ratio_Pval <- formatC(tab_df_int_ordered$Likelihood_Ratio_Pval)
    tab_df_int_ordered$Wald_Test_Pval <- formatC(as.numeric(tab_df_int_ordered$Wald_Test_Pval))
    tab_df_int_ordered$Logrank_Test_Pval <- formatC(as.numeric(tab_df_int_ordered$Logrank_Test_Pval))
    tab_df_int_ordered <- tab_df_int_ordered %>%
      relocate(variable,Hazard_Ratio,ci,High_AboveMedian_Pval,Concordance,Likelihood_Ratio_Pval,Wald_Test_Pval,
               Logrank_Test_Pval,Likelihood_Ratio_adjPval_BH,Wald_Test_adjPval_BH,Logrank_Test_adjPval_BH,
               HazardAssumption_Pval_Global,HazardAssumption_Pval_FeatureCovariate,HazardAssumption_Pval_Feature,HazardAssumption_Pval_Covariate,
               Criteria,Median_Score,Mean_Score,Min_Score,Max_Score)
    new_header_int <- colnames(tab_df_int_ordered)
    colnames(tab_df_int_ordered) <- colnames(out_df_int_top)
    tab_df_int_ordered <- rbind(out_df_int_top,new_header_int,tab_df_int_ordered)
    new_file_int <- gsub(".txt","_ranked.txt",file_made_int)
    write_delim(tab_df_int_ordered, new_file_int, delim = '\t', col_names = F, na = "")
    print("CoxH Ranking Complete")
    
  }
  
}



####----Gene Ranking----####

if (Rank_Genes_Choice == "TRUE") {
  
  ####----Perform Median Cutpoint----####
  print("####----Performing Median Cutpoint----####")
  ##--Median Cutpoint of ssGSEA--##
  gs_name <- "Genes"
  ## Expr summary columns
  expr_df <- as.data.frame(expr)
  expr_df$Median_Score <- rowMedians(as.matrix(expr_df[,c(1:ncol(expr_df))]), na.rm = T)
  expr_df$Mean_Score <- rowMeans(as.matrix(expr_df[,c(1:ncol(expr_df))]), na.rm = T)
  expr_df$Min_Score <- apply(expr_df[,c(1:ncol(expr_df))], 1, FUN = min)
  expr_df$Max_Score <- apply(expr_df[,c(1:ncol(expr_df))], 1, FUN = max)
  expr_summ <- expr_df %>%
    select(Median_Score,Mean_Score,Min_Score,Max_Score)
  rownames(expr_summ) <- gsub("[[:punct:]]",".",rownames(expr_summ))
  rownames(expr_summ) <- gsub(" ",".",rownames(expr_summ))
  expr_summ$Feature <- rownames(expr_summ)
  expr_summ <- expr_summ %>%
    relocate(Feature)
  ## Transpose expression data
  expr_t <- as.data.frame(t(expr))
  new_col_list <- c()
  for (i in colnames(expr_t)) {
    expr_t[,paste(i,'.MedianCutP',sep = "")] <- highlow(expr_t[, which(colnames(expr_t) == i)])
    new_col_list <- c(new_col_list,paste(i,'.MedianCutP',sep = ""))
  }
  ## Reformat
  expr_BIN <- expr_t[,new_col_list]
  expr_BIN$SampleName <- rownames(expr_BIN)
  expr_BIN <- expr_BIN %>%
    relocate(SampleName)
  ## Write to file
  write_delim(expr_BIN,paste(Output_File_Path,Project_Name,"_Genes_MedianCutP.txt", sep = ""), delim = '\t')
  
  
  ####----CoxH Ranking----#### 
  print("####----Performing Coxh Analysis----####")
  gs_name <- "Genes"
  ## Clean Gene Set names
  colnames(expr_BIN) <- gsub("[[:punct:]]",".",colnames(expr_BIN))
  colnames(expr_BIN) <- gsub(" ",".",colnames(expr_BIN))
  ## Merge with Meta data
  expr_BIN_meta <- merge(meta,expr_BIN,by = "SampleName",all.y = T)
  ## Binary columns to loop through
  calc_cols <- grep("..MedianCutP$", colnames(expr_BIN_meta), value = T)
  
  ##--Generate Headers--##
  ProjName_Line <- paste("## Project Name:",Project_Name)
  Survival_Line <- paste("## Survival/Event Variables:",Survival_Time,"///",Survival_ID)
  Samples_Line <- paste("## Number of Samples:",nrow(meta))
  expr_file <- paste("## Expression Matrix File Name:",basename(Expression_Matrix_file))
  meta_file <- paste("## Meta Data File Name:",basename(Meta_Data_File))
  header <- c("variable","var_label","var_type","reference_row","row_type","header_row",
              "N_obs","N_event","N","coefficients_type","coefficients_label","Criteria","term",
              "var_class","var_nlevels","contrasts","contrasts_type","n_obs","n_event",
              "exposure","Hazard_Ratio","std.error","statistic","nevent","conf.low",
              "conf.high","ci","High_AboveMedian_Pval","Concordance","Likelihood_Ratio_Pval","Wald_Test_Pval",
              "Logrank_Test_Pval","Chisq","df","HazardAssumption_Global_Pval","Median_Score","Mean_Score","Min_Score","Max_Score")
  header_add <- c("variable","var_label","var_type","reference_row","row_type","header_row",
                  "N_obs","N_event","N","coefficients_type","coefficients_label","Criteria","term",
                  "var_class","var_nlevels","contrasts","contrasts_type","n_obs","n_event",
                  "exposure","Hazard_Ratio","std.error","statistic","nevent","conf.low",
                  "conf.high","ci","High_AboveMedian_Pval","Concordance","Likelihood_Ratio_Pval","Wald_Test_Pval",
                  "Logrank_Test_Pval","Chisq_Global","df_Global","HazardAssumption_Pval_Global","Chisq_Feature","df_Feature","HazardAssumption_Pval_Feature"
                  ,"Chisq_covariate","df_Covariate","HazardAssumption_Pval_Covariate","Median_Score","Mean_Score","Min_Score","Max_Score")
  header_int <- c("variable","var_label","var_type","reference_row","row_type","header_row",
                  "N_obs","N_event","N","coefficients_type","coefficients_label","Criteria","term",
                  "var_class","var_nlevels","contrasts","contrasts_type","n_obs","n_event",
                  "exposure","Hazard_Ratio","std.error","statistic","nevent","conf.low",
                  "conf.high","ci","High_AboveMedian_Pval","Concordance","Likelihood_Ratio_Pval","Wald_Test_Pval",
                  "Logrank_Test_Pval","Chisq_Global","df_Global","HazardAssumption_Pval_Global","Chisq_Feature","df_Feature","HazardAssumption_Pval_Feature",
                  "Chisq_Covariate","df_Covariate","HazardAssumption_Pval_Covariate","Chisq_FeatureCovariate","df_FeatureCovariate","HazardAssumption_Pval_FeatureCovariate",
                  "Median_Score","Mean_Score","Min_Score","Max_Score")
  
  
  ####----Univariate CoxH Analysis----####
  ## If no covariate being analyzed
  if (is.na(Covariate_Column_Label)) {
    ## Add Additional Headers
    line_analysis <- paste("## Samples were grouped to Above or Below median based on raw gene expression then ran through Cox Proportional Hazard Analysis")
    line_analysis2 <- paste("## Genes were ranked by Coxph Likelihood_Ratio_Pval")
    line_analysis3 <- paste("## Hazard Ratios > 1 represent high (above median) raw gene expression associated with High-Risk")
    param_lines <- paste(ProjName_Line,expr_file,meta_file,Survival_Line,Samples_Line,line_analysis,line_analysis2,line_analysis3,sep = "\n")
    file_made <- paste(Output_File_Path,Project_Name,"_",gs_name,"_coxh.txt", sep = "")
    write(param_lines,file = file_made, append = T, sep = '\t', ncolumns = 39)
    
    ##--CoxH Loop--##
    out_df <- as.data.frame(read_delim(file_made,delim = '\t', col_names = F))
    out_df[c(2:39)] <- NA
    out_df <- rbind(out_df,header)
    write(header,file = file_made, append = T, sep = '\t', ncolumns = 39)
    k <- 1
    for (j in calc_cols) {
      meta_expr_sub <- expr_BIN_meta[,c(Survival_Time,Survival_ID,j)]
      meta_expr_sub[,j] <- as.factor(meta_expr_sub[,j])
      colnames(meta_expr_sub)[which(colnames(meta_expr_sub) == Survival_Time)] <- "time"
      colnames(meta_expr_sub)[which(colnames(meta_expr_sub) == Survival_ID)] <- "ID"
      colnames(meta_expr_sub)[which(colnames(meta_expr_sub) == j)] <- "Feature"
      
      if (length(levels(as.factor(meta_expr_sub[,"Feature"]))) > 1) {
        
        meta_expr_sub[,"Feature"] <- relevel(meta_expr_sub[,"Feature"], ref = "Low_BelowMedian")
        
        tab <- coxph(as.formula(paste("Surv(time,ID) ~ Feature",sep = "")),
                     data = meta_expr_sub)
        temp_tab <- tab %>% 
          gtsummary::tbl_regression(exp = TRUE) %>%
          as_gt()
        temp_tab_df <- as.data.frame(temp_tab)
        ## Get additional Pvals
        out <- capture.output(summary(tab))
        con_line <- grep("^Concordance=",out,value = T)
        con_v <- str_split(con_line," ")[[1]][2]
        lik_line <- grep("^Likelihood ratio test=",out,value = T)
        lik_p <- str_split(lik_line,"=")[[1]][3]
        wal_line <- grep("^Wald test",out,value = T)
        wal_p <- str_split(wal_line,"=")[[1]][3]
        sco_line <- grep("^Score ",out,value = T)
        sco_p <- str_split(sco_line,"=")[[1]][3]
        if (min(table(meta_expr_sub$Feature)) > 1 | !any(grepl("Inf",temp_tab_df$ci)) & temp_tab_df[3,21] != "") {
          zph <- capture.output(cox.zph(tab))
          HA_line <- zph[3]
          chi_n <- str_split(HA_line,"\\s+")[[1]][2]
          dof_n <- str_split(HA_line,"\\s+")[[1]][3]
          pva_n <- str_split(HA_line,"\\s+")[[1]][4]
        } else {
          chi_n <- NA
          dof_n <- NA
          pva_n <- NA
        }
        temp_tab_vect <- as.character(c(temp_tab_df[3,]))
        temp_tab_vect <- c(temp_tab_vect,con_v,lik_p,wal_p,sco_p,chi_n,dof_n,pva_n)
        ## Add summary values
        variable <- file_path_sans_ext(j)
        var_scores <- unname(unlist(c(expr_summ[variable,c(2:5)])))
        temp_tab_vect <- c(temp_tab_vect,var_scores)
        temp_tab_vect[1] <- variable
        temp_tab_vect[2] <- j
        temp_tab_vect[12] <- "High_AboveMedian"
        temp_tab_vect[13] <- paste0(variable,"High_AboveMedian")
        out_df <- rbind(out_df,temp_tab_vect)
        write(temp_tab_vect,file = file_made, append = T, sep = '\t', ncolumns = 39)
        if (k %% 100 == 0) {
          print(paste("Features Processed:",k))
        }
        k <- k+1
        
      }
    }
    
    ##--Cal Adj Pvals Univar--##
    out_df_top <- out_df[grep("##",out_df[,1]),]
    out_df_top[,c(40:42)] <- NA
    tab_df3 <- out_df[grep("##",out_df[,1],invert = T),]
    colnames(tab_df3) <- header
    tab_df3 <- tab_df3[-1,]
    Likelihood_Ratio_adjPval_BH <- p.adjust(as.numeric(tab_df3$Likelihood_Ratio_Pval), method = "BH")
    Wald_Test_adjPval_BH <- p.adjust(as.numeric(tab_df3$Wald_Test_Pval), method = "BH")
    Logrank_Test_adjPval_BH <- p.adjust(as.numeric(tab_df3$Logrank_Test_Pval), method = "BH")
    tab_df3 <- cbind(tab_df3,Likelihood_Ratio_adjPval_BH,Wald_Test_adjPval_BH,Logrank_Test_adjPval_BH)
    tab_df3$Likelihood_Ratio_Pval <- as.numeric(tab_df3$Likelihood_Ratio_Pval)
    tab_df3_ordered <- tab_df3[order(tab_df3$Likelihood_Ratio_Pval, decreasing = F),]
    tab_df3_ordered$Likelihood_Ratio_Pval <- formatC(tab_df3_ordered$Likelihood_Ratio_Pval)
    tab_df3_ordered$Wald_Test_Pval <- formatC(as.numeric(tab_df3_ordered$Wald_Test_Pval))
    tab_df3_ordered$Logrank_Test_Pval <- formatC(as.numeric(tab_df3_ordered$Logrank_Test_Pval))
    ## Reorder Columns and write final table
    tab_df3_ordered <- tab_df3_ordered %>%
      relocate(variable,Hazard_Ratio,ci,High_AboveMedian_Pval,Concordance,Likelihood_Ratio_Pval,Wald_Test_Pval,
               Logrank_Test_Pval,Likelihood_Ratio_adjPval_BH,Wald_Test_adjPval_BH,Logrank_Test_adjPval_BH,HazardAssumption_Global_Pval,
               Chisq,df,Criteria,Median_Score,Mean_Score,Min_Score,Max_Score)
    new_header <- colnames(tab_df3_ordered)
    colnames(tab_df3_ordered) <- colnames(out_df_top)
    tab_df3_ordered <- rbind(out_df_top,new_header,tab_df3_ordered)
    new_file <- gsub(".txt","_ranked.txt",file_made)
    write_delim(tab_df3_ordered, new_file, delim = '\t', col_names = F, na = "")
    print("CoxH Ranking Complete")
  }
  
  ####----Bivariate CoxH Analysis----####
  if (!is.na(Covariate_Column_Label)) {
    
    ##--Generate Covariate Specific Headers--##
    ## Covariate Specific Lines
    line_covar <- paste("## Covariate and Reference Variable:",Covariate_Column_Label,"Referencing",Covariate_Reference)
    covar_brk <- as.data.frame(table(expr_BIN_meta[,Covariate_Column_Label]))
    covars <- c()
    for (k in 1:nrow(covar_brk)) {
      var <- as.character(covar_brk[k,1])
      num <- as.character(covar_brk[k,2])
      varnum <- c("### Variable:",var,"(","N =",num,")")
      covars <- c(covars,paste(varnum,collapse = " "))
    }
    line_covar2 <- paste("## Covariate Breakdown:\n",paste(covars,collapse = "\n"),sep = "")
    
    ## Geneset Lines
    line_analysis_int <- paste("## Samples were grouped to Above or Below median based on raw gene expression then ran through interactive Cox Proportional Hazard Analysis with the covariate")
    line_analysis_add <- paste("## Samples were grouped to Above or Below median based on raw gene expression then ran through additive Cox Proportional Hazard Analysis with the covariate")
    line_analysis2 <- paste("## Pathways were ranked by Coxph Likelihood_Ratio_Pval")
    line_analysis3 <- paste("## Hazard Ratios > 1 represent pathway scores associated with High-Risk")
    param_lines_add <- paste(ProjName_Line,expr_file,meta_file,Survival_Line,Samples_Line,line_analysis_add,line_analysis2,line_analysis3,line_covar,line_covar2,sep = "\n")
    param_lines_int <- paste(ProjName_Line,expr_file,meta_file,Survival_Line,Samples_Line,line_analysis_int,line_analysis2,line_analysis3,line_covar,line_covar2,sep = "\n")
    file_made_add <- paste(Output_File_Path,Project_Name,"_",gs_name,"_",Covariate_Column_Label,"_Additive_coxh.txt", sep = "")
    file_made_int <- paste(Output_File_Path,Project_Name,"_",gs_name,"_",Covariate_Column_Label,"_Interactive_coxh.txt", sep = "")
    write(param_lines_add,file = file_made_add, append = T, sep = '\t', ncolumns = 45)
    write(param_lines_int,file = file_made_int, append = T, sep = '\t', ncolumns = 48)
    
    ##--CoxH Loop--##
    out_df_add <- as.data.frame(read_delim(file_made_add,delim = '\t', col_names = F))
    out_df_add[c(2:45)] <- NA
    out_df_add <- rbind(out_df_add,header_add)
    out_df_int <- as.data.frame(read_delim(file_made_int,delim = '\t', col_names = F))
    out_df_int[c(2:48)] <- NA
    out_df_int <- rbind(out_df_int,header_int)
    write(header,file = file_made_add, append = T, sep = '\t', ncolumns = 45)
    write(header,file = file_made_int, append = T, sep = '\t', ncolumns = 48)
    k <- 1
    for (j in calc_cols) {
      meta_ssgsea_sub <- expr_BIN_meta[,c(Survival_Time,Survival_ID,j,Covariate_Column_Label)]
      # factor hi/lo column
      meta_ssgsea_sub[,j] <- as.factor(meta_ssgsea_sub[,j])
      meta_ssgsea_sub[,Covariate_Column_Label] <- as.factor(meta_ssgsea_sub[,Covariate_Column_Label])
      colnames(meta_ssgsea_sub)[which(colnames(meta_ssgsea_sub) == Survival_Time)] <- "time"
      colnames(meta_ssgsea_sub)[which(colnames(meta_ssgsea_sub) == Survival_ID)] <- "ID"
      colnames(meta_ssgsea_sub)[which(colnames(meta_ssgsea_sub) == j)] <- "Feature"
      colnames(meta_ssgsea_sub)[which(colnames(meta_ssgsea_sub) == Covariate_Column_Label)] <- "Covariate"
      # factor reference column
      if (is.na(Covariate_Reference)) {
        current_ref <- levels(meta_ssgsea_sub[,"Covariate"])[1]
        print(paste("No Covariate Reference Given. Using ",current_ref," as Coxh reference variable.",sep = ""))
      } else {
        meta_ssgsea_sub[,"Covariate"] <- relevel(meta_ssgsea_sub[,"Covariate"], ref = Covariate_Reference)
      }
      if (length(levels(as.factor(meta_ssgsea_sub[,"Feature"]))) > 1) {
        
        meta_ssgsea_sub[,"Feature"] <- relevel(meta_ssgsea_sub[,"Feature"], ref = "Low_BelowMedian")
        
        ##--Additive--##
        tab_add <- coxph(as.formula(paste("Surv(time,ID) ~ Feature + Covariate",sep = "")),
                         data = meta_ssgsea_sub)
        temp_tab_add <- tab_add %>% 
          gtsummary::tbl_regression(exp = TRUE) %>%
          as_gt()
        temp_tab_add_df <- as.data.frame(temp_tab_add)
        out_add <- capture.output(summary(tab_add))
        con_line_add <- grep("^Concordance=",out_add,value = T)
        con_v_add <- str_split(con_line_add," ")[[1]][2]
        lik_line_add <- grep("^Likelihood ratio test=",out_add,value = T)
        lik_p_add <- str_split(lik_line_add,"=")[[1]][3]
        wal_line_add <- grep("^Wald test",out_add,value = T)
        wal_p_add <- str_split(wal_line_add,"=")[[1]][3]
        sco_line_add <- grep("^Score ",out_add,value = T)
        sco_p_add <- str_split(sco_line_add,"=")[[1]][3]
        if (min(table(meta_ssgsea_sub$Feature)) > 1 | !any(grepl("Inf",temp_tab_add_df$ci)) & temp_tab_df[3,21] != "") {
          xph_add <- capture.output(cox.zph(tab_add))
          HA_Feat_line_add <- xph_add[2]
          chi_Feat_n_add <- str_split(HA_Feat_line_add,"\\s+")[[1]][2]
          dof_Feat_n_add <- str_split(HA_Feat_line_add,"\\s+")[[1]][3]
          pva_Feat_n_add <- str_split(HA_Feat_line_add,"\\s+")[[1]][4]
          HA_Covar_line_add <- xph_add[3]
          chi_Covar_n_add <- str_split(HA_Covar_line_add,"\\s+")[[1]][2]
          dof_Covar_n_add <- str_split(HA_Covar_line_add,"\\s+")[[1]][3]
          pva_Covar_n_add <- str_split(HA_Covar_line_add,"\\s+")[[1]][4]
          HA_Global_line_add <- xph_add[4]
          chi_Global_n_add <- str_split(HA_Global_line_add,"\\s+")[[1]][2]
          dof_Global_n_add <- str_split(HA_Global_line_add,"\\s+")[[1]][3]
          pva_Global_n_add <- str_split(HA_Global_line_add,"\\s+")[[1]][4]
        } else {
          chi_Feat_n_add <- NA
          dof_Feat_n_add <- NA
          pva_Feat_n_add <- NA
          chi_Covar_n_add <- NA
          dof_Covar_n_add <- NA
          pva_Covar_n_add <- NA
          chi_Global_n_add <- NA
          dof_Global_n_add <- NA
          pva_Global_n_add <- NA
        }
        
        temp_tab_vect_add <- as.character(c(temp_tab_add_df[3,]))
        temp_tab_vect_add <- c(temp_tab_vect_add,con_v_add,lik_p_add,wal_p_add,sco_p_add,
                               chi_Feat_n_add,dof_Feat_n_add,pva_Feat_n_add,
                               chi_Covar_n_add,dof_Covar_n_add,pva_Covar_n_add,
                               chi_Global_n_add,dof_Global_n_add,pva_Global_n_add)
        ## Add summary values
        variable <- file_path_sans_ext(j)
        var_scores <- unname(unlist(c(expr_summ[variable,c(2:5)])))
        temp_tab_vect_add <- c(temp_tab_vect_add,var_scores)
        temp_tab_vect_add[1] <- variable
        temp_tab_vect_add[2] <- j
        temp_tab_vect_add[12] <- "High_AboveMedian"
        temp_tab_vect_add[13] <- paste0(variable,"High_AboveMedian")
        out_df_add <- rbind(out_df_add,temp_tab_vect_add)
        write(temp_tab_vect_add,file = file_made_add, append = T, sep = '\t', ncolumns = 39)
        
        ##--Interactive--##
        tab_int <- coxph(as.formula(paste("Surv(time,ID) ~ Feature * Covariate",sep = "")),
                         data = meta_ssgsea_sub)
        temp_tab_int <- tab_int %>% 
          gtsummary::tbl_regression(exp = TRUE) %>%
          as_gt()
        temp_tab_int_df <- as.data.frame(temp_tab_int)
        out_int <- capture.output(summary(tab_int))
        con_line_int <- grep("^Concordance=",out_int,value = T)
        con_v_int <- str_split(con_line_int," ")[[1]][2]
        lik_line_int <- grep("^Likelihood ratio test=",out_int,value = T)
        lik_p_int <- str_split(lik_line_int,"=")[[1]][3]
        wal_line_int <- grep("^Wald test",out_int,value = T)
        wal_p_int <- str_split(wal_line_int,"=")[[1]][3]
        sco_line_int <- grep("^Score ",out_int,value = T)
        sco_p_int <- str_split(sco_line_int,"=")[[1]][3]
        if (min(table(meta_ssgsea_sub$Feature)) > 1 | !any(grepl("Inf",temp_tab_int_df$ci)) & temp_tab_df[3,21] != "") {
          xph_int <- capture.output(cox.zph(tab_int))
          HA_Feat_line_int <- xph_int[2]
          chi_Feat_n_int <- str_split(HA_Feat_line_int,"\\s+")[[1]][2]
          dof_Feat_n_int <- str_split(HA_Feat_line_int,"\\s+")[[1]][3]
          pva_Feat_n_int <- str_split(HA_Feat_line_int,"\\s+")[[1]][4]
          HA_Covar_line_int <- xph_int[3]
          chi_Covar_n_int <- str_split(HA_Covar_line_int,"\\s+")[[1]][2]
          dof_Covar_n_int <- str_split(HA_Covar_line_int,"\\s+")[[1]][3]
          pva_Covar_n_int <- str_split(HA_Covar_line_int,"\\s+")[[1]][4]
          HA_FeatCovar_line_int <- xph_int[4]
          chi_FeatCovar_n_int <- str_split(HA_FeatCovar_line_int,"\\s+")[[1]][2]
          dof_FeatCovar_n_int <- str_split(HA_FeatCovar_line_int,"\\s+")[[1]][3]
          pva_FeatCovar_n_int <- str_split(HA_FeatCovar_line_int,"\\s+")[[1]][4]
          HA_Global_line_int <- xph_int[5]
          chi_Global_n_int <- str_split(HA_Global_line_int,"\\s+")[[1]][2]
          dof_Global_n_int <- str_split(HA_Global_line_int,"\\s+")[[1]][3]
          pva_Global_n_int <- str_split(HA_Global_line_int,"\\s+")[[1]][4]
        } else {
          chi_Feat_n_int <- NA
          dof_Feat_n_int <- NA
          pva_Feat_n_int <- NA
          chi_Covar_n_int <- NA
          dof_Covar_n_int <- NA
          pva_Covar_n_int <- NA
          chi_Global_n_int <- NA
          dof_Global_n_int <- NA
          pva_Global_n_int <- NA
          chi_FeatCovar_n_int <- NA
          dof_FeatCovar_n_int <- NA
          pva_FeatCovar_n_int <- NA
        }
        
        temp_tab_vect_int <- as.character(c(temp_tab_int_df[3,]))
        temp_tab_vect_int <- c(temp_tab_vect_int,con_v_int,lik_p_int,wal_p_int,sco_p_int,
                               chi_Feat_n_int,dof_Feat_n_int,pva_Feat_n_int,
                               chi_Covar_n_int,dof_Covar_n_int,pva_Covar_n_int,
                               chi_FeatCovar_n_int,dof_FeatCovar_n_int,pva_FeatCovar_n_int,
                               chi_Global_n_int,dof_Global_n_int,pva_Global_n_int)
        ## Add summary values
        variable <- file_path_sans_ext(j)
        var_scores <- unname(unlist(c(expr_summ[variable,c(2:5)])))
        temp_tab_vect_int <- c(temp_tab_vect_int,var_scores)
        temp_tab_vect_int[1] <- variable
        temp_tab_vect_int[2] <- j
        temp_tab_vect_int[12] <- "High_AboveMedian"
        temp_tab_vect_int[13] <- paste0(variable,"High_AboveMedian")
        out_df_int <- rbind(out_df_int,temp_tab_vect_int)
        write(temp_tab_vect_int,file = file_made_int, append = T, sep = '\t', ncolumns = 48)
        if (k %% 100 == 0) {
          print(paste("Features Processed:",k))
        }
        k <- k+1
      }
    }
    
    
    ##--Calc Adj Pval Add--##
    out_df_add_top <- out_df_add[grep("##",out_df_add[,1]),]
    out_df_add_top[,c(46:48)] <- NA
    tab_df_add <- out_df_add[grep("##",out_df_add[,1],invert = T),]
    colnames(tab_df_add) <- header_add
    tab_df_add <- tab_df_add[-1,]
    Likelihood_Ratio_adjPval_BH <- p.adjust(as.numeric(tab_df_add$Likelihood_Ratio_Pval), method = "BH")
    Wald_Test_adjPval_BH <- p.adjust(as.numeric(tab_df_add$Wald_Test_Pval), method = "BH")
    Logrank_Test_adjPval_BH <- p.adjust(as.numeric(tab_df_add$Logrank_Test_Pval), method = "BH")
    tab_df_add <- cbind(tab_df_add,Likelihood_Ratio_adjPval_BH,Wald_Test_adjPval_BH,Logrank_Test_adjPval_BH)
    tab_df_add$Likelihood_Ratio_Pval <- as.numeric(tab_df_add$Likelihood_Ratio_Pval)
    tab_df_add_ordered <- tab_df_add[order(tab_df_add$Likelihood_Ratio_Pval, decreasing = F),]
    tab_df_add_ordered$Likelihood_Ratio_Pval <- formatC(tab_df_add_ordered$Likelihood_Ratio_Pval)
    tab_df_add_ordered$Wald_Test_Pval <- formatC(as.numeric(tab_df_add_ordered$Wald_Test_Pval))
    tab_df_add_ordered$Logrank_Test_Pval <- formatC(as.numeric(tab_df_add_ordered$Logrank_Test_Pval))
    tab_df_add_ordered <- tab_df_add_ordered %>%
      relocate(variable,Hazard_Ratio,ci,High_AboveMedian_Pval,Concordance,Likelihood_Ratio_Pval,Wald_Test_Pval,
               Logrank_Test_Pval,Likelihood_Ratio_adjPval_BH,Wald_Test_adjPval_BH,Logrank_Test_adjPval_BH,
               HazardAssumption_Pval_Global,HazardAssumption_Pval_Feature,HazardAssumption_Pval_Covariate,Criteria,
               Median_Score,Mean_Score,Min_Score,Max_Score)
    new_header_add <- colnames(tab_df_add_ordered)
    colnames(tab_df_add_ordered) <- colnames(out_df_add_top)
    tab_df_add_ordered <- rbind(out_df_add_top,new_header_add,tab_df_add_ordered)
    new_file_add <- gsub(".txt","_ranked.txt",file_made_add)
    write_delim(tab_df_add_ordered, new_file_add, delim = '\t', col_names = F, na = "")
    
    
    ##--Calc Adj Pval Int--##
    out_df_int_top <- out_df_int[grep("##",out_df_int[,1]),]
    out_df_int_top[,c(49:51)] <- NA
    tab_df_int <- out_df_int[grep("##",out_df_int[,1],invert = T),]
    colnames(tab_df_int) <- header_int
    tab_df_int <- tab_df_int[-1,]
    Likelihood_Ratio_adjPval_BH <- p.adjust(as.numeric(tab_df_int$Likelihood_Ratio_Pval), method = "BH")
    Wald_Test_adjPval_BH <- p.adjust(as.numeric(tab_df_int$Wald_Test_Pval), method = "BH")
    Logrank_Test_adjPval_BH <- p.adjust(as.numeric(tab_df_int$Logrank_Test_Pval), method = "BH")
    tab_df_int <- cbind(tab_df_int,Likelihood_Ratio_adjPval_BH,Wald_Test_adjPval_BH,Logrank_Test_adjPval_BH)
    tab_df_int$Likelihood_Ratio_Pval <- as.numeric(tab_df_int$Likelihood_Ratio_Pval)
    tab_df_int_ordered <- tab_df_int[order(tab_df_int$Likelihood_Ratio_Pval, decreasing = F),]
    tab_df_int_ordered$Likelihood_Ratio_Pval <- formatC(tab_df_int_ordered$Likelihood_Ratio_Pval)
    tab_df_int_ordered$Wald_Test_Pval <- formatC(as.numeric(tab_df_int_ordered$Wald_Test_Pval))
    tab_df_int_ordered$Logrank_Test_Pval <- formatC(as.numeric(tab_df_int_ordered$Logrank_Test_Pval))
    tab_df_int_ordered <- tab_df_int_ordered %>%
      relocate(variable,Hazard_Ratio,ci,High_AboveMedian_Pval,Concordance,Likelihood_Ratio_Pval,Wald_Test_Pval,
               Logrank_Test_Pval,Likelihood_Ratio_adjPval_BH,Wald_Test_adjPval_BH,Logrank_Test_adjPval_BH,
               HazardAssumption_Pval_Global,HazardAssumption_Pval_FeatureCovariate,HazardAssumption_Pval_Feature,HazardAssumption_Pval_Covariate,
               Criteria,Median_Score,Mean_Score,Min_Score,Max_Score)
    new_header_int <- colnames(tab_df_int_ordered)
    colnames(tab_df_int_ordered) <- colnames(out_df_int_top)
    tab_df_int_ordered <- rbind(out_df_int_top,new_header_int,tab_df_int_ordered)
    new_file_int <- gsub(".txt","_ranked.txt",file_made_int)
    write_delim(tab_df_int_ordered, new_file_int, delim = '\t', col_names = F, na = "")
    print("CoxH Ranking Complete")
    
  }
  
}



ET <- Sys.time()

time <- ET-ST
print(time)








