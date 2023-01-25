

####----Input----####

ProjectName <- "Pan ICI Checkpoint iAtlas - Skin Cancer"

ExpressionMatrix_file <- "Example_Data/Expression_Data_PAN_ICI_iAtlas_SkinCancer.zip"

MetaData_file <- "Example_Data/Clinical_Data_PAN_ICI_iAtlas_SkinCancer.txt"

MetaParam_File <- "Example_Data/Clinical_Parameters_PAN_ICI_iAtlas_SkinCancer.txt"


##--Advanced Set-Up--##

## Pre-Selected Inputs
# An option from the meta, All, or NULL
PreSelect_SamplyType <- NULL
PreSelect_Feature <- "All"
# An option from the meta or NULL
PreSelect_SubFeature <- NULL
PreSelect_SecondaryFeature <- "Responder"

# DO NOT CHANGE - only adjust file path if needed
GeneSet_File <- "GeneSet_Data/GeneSet_List.RData"
GeneSetTable_File <- "GeneSet_Data/GeneSet_CatTable.zip"
About_MD_File <- "App_Markdowns/PurposeAndMethods.Rmd"



############################ Copyright 2022 Moffitt Cancer Center ############################
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##############################################################################################


####----Install and load packages----####


## Check if Immune deconvolution package is installed
immudecon <- "immunedeconv"
immudecon_check <- immudecon %in% rownames(installed.packages())
if (immudecon_check == TRUE) {
  library(immunedeconv)
}
packages <- c("shiny","shinythemes","shinyjqui","gtsummary","tidyr","RColorBrewer",
              "dplyr","DT","ggplot2","ggpubr","tibble","survival","pheatmap","stringr",
              "readr","shinycssloaders","survminer","gridExtra","viridis","plotly","ggrepel")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))
#bioconductor packages
bioCpacks <- c("GSVA")
installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
if (any(installed_packages_BIOC == FALSE)) {
  BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
}
invisible(lapply(bioCpacks, library, character.only = TRUE))




####----Read In Files----####

##--Meta--##
meta <- as.data.frame(read_delim(MetaData_file,delim = '\t', col_names = T))

##--Meta Param--##
MetaParam <- as.data.frame(read_delim(MetaParam_File,delim = '\t',col_names = F))
MetaParam[,2] <- gsub(" ","",MetaParam[,2])
## Subset meta columns by category
metacol_samplenames <- MetaParam[which(MetaParam[,2] == "SampleName"),1]
metacol_feature <- MetaParam[which(MetaParam[,2] == "Feature"),1]
metacol_description <- MetaParam[which(MetaParam[,2] == "Description"),1]
if (("SampleType" %in% MetaParam[,2]) == TRUE) {
  metacol_sampletype <- MetaParam[which(MetaParam[,2] == "SampleType"),1]
}
if (("SampleType" %in% MetaParam[,2]) == FALSE) {
  metacol_sampletype <- NULL
}
## Get surv time and ID columns - move OS to the front of the list
metacol_survtime <- MetaParam[which(MetaParam[,2] == "SurvivalTime"),1]
metacol_survid <- MetaParam[which(MetaParam[,2] == "SurvivalID"),1]

## Rename and move Sample Name column to the front
colnames(meta)[which(colnames(meta) == metacol_samplenames)] <- "SampleName"
meta <- meta %>%
  relocate(SampleName)
## Replace any special characters to make uniform with expression
meta[,1] <- gsub("[[:punct:]]","_",meta[,1])


##--Expression--##
expr <- as.data.frame(read_delim(ExpressionMatrix_file,delim = '\t', col_names = T))
rownames(expr) <- expr[,1]
expr <- expr[,-1]
## Replace any special characters to make uniform with expression
colnames(expr) <- gsub("[[:punct:]]","_",colnames(expr))

sampsame <- intersect(meta[,"SampleName"],colnames(expr))
meta <- meta[which(meta[,"SampleName"] %in% sampsame),]
expr <- expr[,sampsame]


##--Gene Set--##
# R Data list load function for naming
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
gs <- loadRData(GeneSet_File)
# Gene - "Gene Set"
exprGenes <- rownames(expr)
GeneGS_table <- data.frame(Genes = exprGenes)
# Gene Set Table
if (exists("GeneSetTable_File") == TRUE) {
  if (file.exists(GeneSetTable_File) == TRUE) {
    
    GeneSetTable_og <- as.data.frame(read_delim(GeneSetTable_File, delim = '\t', col_names = T))
    gsTab = TRUE
    
  }
  if (file.exists(GeneSetTable_File) == FALSE) {
    
    GeneSets <- names(gs)
    GeneSetTable_og <- data.frame(GeneSets)
    gsTab = FALSE
    
  }
}
if (exists("GeneSetTable_File") == FALSE) {
  
  GeneSets <- names(gs)
  GeneSetTable_og <- data.frame(GeneSets)
  gsTab = FALSE
  
}



## Pre-Selected Inputs
# An option from the meta, All, or NULL
if (is.null(PreSelect_SamplyType) == FALSE) {
  if (grepl("all",PreSelect_SamplyType, ignore.case = T) == TRUE) {
    PreSelect_SamplyType <- "All_Sample_Types"
  }
}
if (is.null(PreSelect_Feature) == FALSE) {
  if (grepl("all",PreSelect_Feature, ignore.case = T) == TRUE) {
    PreSelect_Feature <- "All_Features"
  }
}


####----Functions----####

StatCols <- c("MedianCutP","QuartileCutP","OptimalCutP","TopBottomCutP","UserCutP")

## Quartile Conversion
quartile_conversion = function(mat) {
  new_mat = mat;
  new_mat[mat <=  quantile(as.numeric(mat), na.rm = T)[2]] = "Q1_Low";
  new_mat[mat > quantile(as.numeric(mat), na.rm = T)[2] & mat <= quantile(mat)[3]] = "Q2_MedLow";
  new_mat[mat > quantile(as.numeric(mat), na.rm = T)[3] & mat <= quantile(mat)[4]] = "Q3_MedHigh";
  new_mat[mat > quantile(as.numeric(mat), na.rm = T)[4]] = "Q4_High";
  return (new_mat)
}

## High-Low
highlow = function(mat) {
  new_mat = mat;
  new_mat[mat > quantile(as.numeric(mat), na.rm = T)[3]] = "High";
  new_mat[mat <= quantile(as.numeric(mat), na.rm = T)[3]] = "Low";
  return (new_mat)
}

quantile_conversion = function(mat,cutoff) {
  new_mat = mat;
  new_mat[mat >= quantile(as.numeric(mat),1-cutoff, na.rm = T)] = "High";
  new_mat[mat <= quantile(as.numeric(mat),cutoff, na.rm = T)] = "Low";
  new_mat[mat > quantile(as.numeric(mat),cutoff, na.rm = T) & mat < quantile(mat,1-cutoff, na.rm = T)] = "BetweenCutoff";
  return (new_mat)
}

quantile_conversion2 = function(mat,cutoff) {
  new_mat = mat;
  new_mat[mat > quantile(as.numeric(mat),cutoff, na.rm = T)] = "High";
  new_mat[mat <= quantile(as.numeric(mat),cutoff, na.rm = T)] = "Low";
  return (new_mat)
}

#https://stackoverflow.com/questions/71339547/how-to-add-a-label-to-the-x-y-axis-whenever-a-vertical-horizontal-line-is-ad
add_x_intercepts <- function(p) {
  
  p2 <- ggplot_build(p)
  breaks <- p2$layout$panel_params[[1]]$x$breaks
  breaks <- breaks[!is.na(breaks)]
  
  vals <- unlist(lapply(seq_along(p$layers), function(x) {
    d <- layer_data(p, x)
    if('xintercept' %in% names(d)) d$xintercept else numeric()
  }))
  
  p + scale_x_continuous(breaks = sort(c(vals, breaks)))
}

get_lik_pval <- function(x_tab) {
  tab <- x_tab
  out <- capture.output(summary(tab))
  lik_line <- grep("^Likelihood ratio test=",out,value = T)
  lik_line_P <- as.numeric(str_split(str_split(lik_line,", ")[[1]][2],"=")[[1]][2])
  return(lik_line_P)
}

gsubCheck <- function(string) {
  string2 <- gsub("\\+","POS",string)
  string3 <- gsub("[[:punct:]]","_",string2)
  return(string3)
}


####----Immune Deconvolution----####


PreProcessed_meta_cols <- c(grep("_PreProcessedScore$",metacol_feature,value = T))
metacol_feature_noPreProcessedScore <- c(grep("_PreProcessedScore$",metacol_feature,value = T, invert = T))
#decon_bin_cols <- c()
decon_score_cols <- c()
mcp_check <- FALSE
est_check <- FALSE



## If immune deconv cols in meta
if (length(PreProcessed_meta_cols) > 0) {
  
  immDeconvCols <- c()
  ## Check in immune deconvolution scores were pre-processed
  mcp_devonv_scores <- grep("mcp_counter_PreProcessedScore$",metacol_feature,value = T)
  if (length(mcp_devonv_scores) > 0) {
    mcp_check <- TRUE
    immDeconvCols <- c(immDeconvCols,mcp_devonv_scores)
    decon_score_cols <- c(decon_score_cols,mcp_devonv_scores)
    mcp_decon_gstab <- data.frame(GeneSet_Database = "Pre-Processed Scores",
                                  GeneSet_Category = "Immune Deconvolution Cell Types",
                                  GeneSet_Sub_Category = "MCP Counter Deconvolution Method",
                                  GeneSet_Name = gsub("_PreProcessedScore","",mcp_devonv_scores))
    if (gsTab == TRUE) {
      GeneSetTable_og <- rbind(GeneSetTable_og,mcp_decon_gstab)
    }
    if (gsTab == FALSE) {
      col_name <- colnames(GeneSetTable_og)[1]
      mcp_decon_gstab <- mcp_decon_gstab[,3, drop = F]
      colnames(mcp_decon_gstab)[1] <- col_name
      GeneSetTable_og <- rbind(GeneSetTable_og,mcp_decon_gstab)
    }
  }
  
  estimate_scores <- grep("estimate_PreProcessedScore$",metacol_feature,value = T)
  if (length(estimate_scores) > 0) {
    est_check <- TRUE
    immDeconvCols <- c(immDeconvCols,estimate_scores)
    decon_score_cols <- c(decon_score_cols,estimate_scores)
    estimate_decon_gstab <- data.frame(GeneSet_Database = "Pre-Processed Scores",
                                       GeneSet_Category = "Immune Deconvolution Cell Types",
                                       GeneSet_Sub_Category = "ESTIMATE Deconvolution Method",
                                       GeneSet_Name = gsub("_PreProcessedScore","",estimate_scores))
    if (gsTab == TRUE) {
      GeneSetTable_og <- rbind(GeneSetTable_og,estimate_decon_gstab)
    }
    if (gsTab == FALSE) {
      col_name <- colnames(GeneSetTable_og)[1]
      estimate_decon_gstab <- estimate_decon_gstab[,3, drop = F]
      colnames(estimate_decon_gstab)[1] <- col_name
      GeneSetTable_og <- rbind(GeneSetTable_og,estimate_decon_gstab)
    }
  }
  
  quantiseq_scores <- grep("quantiseq_PreProcessedScore$",metacol_feature,value = T)
  if (length(quantiseq_scores) > 0) {
    immDeconvCols <- c(immDeconvCols,quantiseq_scores)
    decon_score_cols <- c(decon_score_cols,quantiseq_scores)
    quantiseq_decon_gstab <- data.frame(GeneSet_Database = "Pre-Processed Scores",
                                        GeneSet_Category = "Immune Deconvolution Cell Types",
                                        GeneSet_Sub_Category = "quanTIseq Deconvolution Method",
                                        GeneSet_Name = gsub("_PreProcessedScore","",quantiseq_scores))
    if (gsTab == TRUE) {
      GeneSetTable_og <- rbind(GeneSetTable_og,quantiseq_decon_gstab)
    }
    if (gsTab == FALSE) {
      col_name <- colnames(GeneSetTable_og)[1]
      quantiseq_decon_gstab <- quantiseq_decon_gstab[,3, drop = F]
      colnames(quantiseq_decon_gstab)[1] <- col_name
      GeneSetTable_og <- rbind(GeneSetTable_og,quantiseq_decon_gstab)
    }
  }
  
  xcell_scores <- grep("xcell_PreProcessedScore$",metacol_feature,value = T)
  if (length(xcell_scores) > 0) {
    immDeconvCols <- c(immDeconvCols,xcell_scores)
    decon_score_cols <- c(decon_score_cols,xcell_scores)
    xcell_decon_gstab <- data.frame(GeneSet_Database = "Pre-Processed Scores",
                                    GeneSet_Category = "Immune Deconvolution Cell Types",
                                    GeneSet_Sub_Category = "Xcell Deconvolution Method",
                                    GeneSet_Name = gsub("_PreProcessedScore","",xcell_scores))
    if (gsTab == TRUE) {
      GeneSetTable_og <- rbind(GeneSetTable_og,xcell_decon_gstab)
    }
    if (gsTab == FALSE) {
      col_name <- colnames(GeneSetTable_og)[1]
      xcell_decon_gstab <- xcell_decon_gstab[,3, drop = F]
      colnames(xcell_decon_gstab)[1] <- col_name
      GeneSetTable_og <- rbind(GeneSetTable_og,xcell_decon_gstab)
    }
  }
  
  epic_scores <- grep("epic_PreProcessedScore$",metacol_feature,value = T)
  if (length(epic_scores) > 0) {
    immDeconvCols <- c(immDeconvCols,epic_scores)
    decon_score_cols <- c(decon_score_cols,epic_scores)
    epic_decon_gstab <- data.frame(GeneSet_Database = "Pre-Processed Scores",
                                   GeneSet_Category = "Immune Deconvolution Cell Types",
                                   GeneSet_Sub_Category = "EPIC Deconvolution Method",
                                   GeneSet_Name = gsub("_PreProcessedScore","",epic_scores))
    if (gsTab == TRUE) {
      GeneSetTable_og <- rbind(GeneSetTable_og,epic_decon_gstab)
    }
    if (gsTab == FALSE) {
      col_name <- colnames(GeneSetTable_og)[1]
      epic_decon_gstab <- epic_decon_gstab[,3, drop = F]
      colnames(epic_decon_gstab)[1] <- col_name
      GeneSetTable_og <- rbind(GeneSetTable_og,epic_decon_gstab)
    }
  }
  
  abis_scores <- grep("abis_PreProcessedScore$",metacol_feature,value = T)
  if (length(abis_scores) > 0) {
    immDeconvCols <- c(immDeconvCols,abis_scores)
    decon_score_cols <- c(decon_score_cols,abis_scores)
    abis_decon_gstab <- data.frame(GeneSet_Database = "Pre-Processed Scores",
                                   GeneSet_Category = "Immune Deconvolution Cell Types",
                                   GeneSet_Sub_Category = "ABIS Deconvolution Method",
                                   GeneSet_Name = gsub("_PreProcessedScore","",abis_scores))
    if (gsTab == TRUE) {
      GeneSetTable_og <- rbind(GeneSetTable_og,abis_decon_gstab)
    }
    if (gsTab == FALSE) {
      col_name <- colnames(GeneSetTable_og)[1]
      abis_decon_gstab <- abis_decon_gstab[,3, drop = F]
      colnames(abis_decon_gstab)[1] <- col_name
      GeneSetTable_og <- rbind(GeneSetTable_og,abis_decon_gstab)
    }
  }
  
  cibersort_scores <- grep("cibersort_PreProcessedScore$",metacol_feature,value = T)
  if (length(cibersort_scores) > 0) {
    immDeconvCols <- c(immDeconvCols,cibersort_scores)
    decon_score_cols <- c(decon_score_cols,cibersort_scores)
    cibersort_decon_gstab <- data.frame(GeneSet_Database = "Pre-Processed Scores",
                                        GeneSet_Category = "Immune Deconvolution Cell Types",
                                        GeneSet_Sub_Category = "CIBERSORT Deconvolution Method",
                                        GeneSet_Name = gsub("_PreProcessedScore","",cibersort_scores))
    if (gsTab == TRUE) {
      GeneSetTable_og <- rbind(GeneSetTable_og,cibersort_decon_gstab)
    }
    if (gsTab == FALSE) {
      col_name <- colnames(GeneSetTable_og)[1]
      cibersort_decon_gstab <- cibersort_decon_gstab[,3, drop = F]
      colnames(cibersort_decon_gstab)[1] <- col_name
      GeneSetTable_og <- rbind(GeneSetTable_og,cibersort_decon_gstab)
    }
  }
  
  cibersort_abs_scores <- grep("cibersort_abs_PreProcessedScore$",metacol_feature,value = T)
  if (length(cibersort_abs_scores) > 0) {
    immDeconvCols <- c(immDeconvCols,cibersort_abs_scores)
    decon_score_cols <- c(decon_score_cols,cibersort_abs_scores)
    cibersort_abs_decon_gstab <- data.frame(GeneSet_Database = "Pre-Processed Scores",
                                            GeneSet_Category = "Immune Deconvolution Cell Types",
                                            GeneSet_Sub_Category = "CIBERSORT ABS Deconvolution Method",
                                            GeneSet_Name = gsub("_PreProcessedScore","",cibersort_abs_scores))
    if (gsTab == TRUE) {
      GeneSetTable_og <- rbind(GeneSetTable_og,cibersort_abs_decon_gstab)
    }
    if (gsTab == FALSE) {
      col_name <- colnames(GeneSetTable_og)[1]
      cibersort_abs_decon_gstab <- cibersort_abs_decon_gstab[,3, drop = F]
      colnames(cibersort_abs_decon_gstab)[1] <- col_name
      GeneSetTable_og <- rbind(GeneSetTable_og,cibersort_abs_decon_gstab)
    }
  }
  
  other_PreProcessedScores <- setdiff(PreProcessed_meta_cols,immDeconvCols)
  if (length(other_PreProcessedScores) > 0) {
    decon_score_cols <- c(decon_score_cols,other_PreProcessedScores)
    Other_scores_gstab <- data.frame(GeneSet_Database = "Pre-Processed Scores",
                                     GeneSet_Category = "Pre-Processed Scores",
                                     GeneSet_Sub_Category = "Pre-Processed Scores",
                                     GeneSet_Name = gsub("_PreProcessedScore","",other_PreProcessedScores))
    if (gsTab == TRUE) {
      GeneSetTable_og <- rbind(GeneSetTable_og,Other_scores_gstab)
    }
    if (gsTab == FALSE) {
      col_name <- colnames(GeneSetTable_og)[1]
      Other_scores_gstab <- Other_scores_gstab[,3, drop = F]
      colnames(Other_scores_gstab)[1] <- col_name
      GeneSetTable_og <- rbind(GeneSetTable_og,Other_scores_gstab)
    }
  }
  
  
}


## If they were not pre-processed, process mcp counter and estimate methods


if (immudecon_check == TRUE) {
  
  if (mcp_check == FALSE) {
    
    mcp_score_cols <- c()
    mcp_bin_cols <- c()
    
    mcp_counter_decon <- as.data.frame(deconvolute(expr, "mcp_counter"))
    
    rownames(mcp_counter_decon) <- mcp_counter_decon[,1]
    
    mcp_counter_decon <- mcp_counter_decon[,-1]
    
    mcp_counter_decon <- as.data.frame(t(mcp_counter_decon))
    
    colnames(mcp_counter_decon) <- paste(gsub(" ","_",colnames(mcp_counter_decon)),"mcp_counter",sep = "_")
    
    decon_score_cols <- c(decon_score_cols,colnames(mcp_counter_decon))
    
    ## Make Score label rows to add to gene set table
    mcp_decon_gstab <- data.frame(GeneSet_Database = "Pre-Processed Scores",
                                  GeneSet_Category = "Immune Deconvolution Cell Types",
                                  GeneSet_Sub_Category = "MCP Counter Deconvolution Method",
                                  GeneSet_Name = colnames(mcp_counter_decon))
    if (gsTab == TRUE) {
      GeneSetTable_og <- rbind(GeneSetTable_og,mcp_decon_gstab)
    }
    if (gsTab == FALSE) {
      col_name <- colnames(GeneSetTable_og)[1]
      mcp_decon_gstab2 <- mcp_decon_gstab[,3, drop = F]
      colnames(mcp_decon_gstab2)[1] <- col_name
      GeneSetTable_og <- rbind(GeneSetTable_og,mcp_decon_gstab2)
    }
    
    mcp_counter_decon$SampleName <- rownames(mcp_counter_decon)
    mcp_counter_decon <- mcp_counter_decon %>%
      relocate(SampleName)
    
    meta <- merge(meta,mcp_counter_decon)
    
  }
  
}


## If they were not pre-processed, process mcp counter and estimate methods
if (immudecon_check == TRUE) {
  
  if (est_check == FALSE) {
    
    estimate_decon <- as.data.frame(deconvolute(expr, "estimate"))
    
    rownames(estimate_decon) <- estimate_decon[,1]
    
    estimate_decon <- estimate_decon[,-1]
    
    estimate_decon <- as.data.frame(t(estimate_decon))
    
    colnames(estimate_decon) <- paste(gsub(" ","_",colnames(estimate_decon)),"estimate",sep = "_")
    
    decon_score_cols <- c(decon_score_cols,colnames(estimate_decon))
    
    ## Make Score label rows to add to gene set table
    est_decon_gstab <- data.frame(GeneSet_Database = "Pre-Processed Scores",
                                  GeneSet_Category = "Immune Deconvolution Cell Types",
                                  GeneSet_Sub_Category = "ESTIMATE Deconvolution Method",
                                  GeneSet_Name = colnames(estimate_decon))
    if (gsTab == TRUE) {
      GeneSetTable_og <- rbind(GeneSetTable_og,est_decon_gstab)
    }
    if (gsTab == FALSE) {
      col_name <- colnames(GeneSetTable_og)[1]
      est_decon_gstab2 <- est_decon_gstab[,3, drop = F]
      colnames(est_decon_gstab2)[1] <- col_name
      GeneSetTable_og <- rbind(GeneSetTable_og,est_decon_gstab2)
    }
    
    estimate_decon$SampleName <- rownames(estimate_decon)
    estimate_decon <- estimate_decon %>%
      relocate(SampleName)
    
    meta <- merge(meta,estimate_decon)
    
    
  }
  
}


decon_score_cols <- gsub("_PreProcessedScore","",decon_score_cols)
metacol_feature <- c(metacol_feature_noPreProcessedScore,decon_score_cols)
colnames(meta) <- gsub("_PreProcessedScore","",colnames(meta))


####----UI----####

ui <-
  navbarPage(paste("{ ",ProjectName," Survival Analysis }",sep = ""),
             
             ####----Overall Survival Tab----####
             
             tabPanel("Survival Analysis",
                      fluidPage(
                        sidebarLayout(
                          
                          ####----Sidebar Panel----####
                          
                          sidebarPanel(
                            tabsetPanel(
                              id = "survside",
                              
                              ##--Sample Parameters--##
                              
                              tabPanel("Sample Parameters",
                                       p(),
                                       uiOutput("rendSampleTypeSelection"),
                                       uiOutput("rendFeatureSelection"),
                                       uiOutput("rendSubFeatureSelection"),
                                       fluidRow(
                                         column(6,
                                                uiOutput("rendSurvivalType_time")
                                         ),
                                         column(6,
                                                uiOutput("rendSurvivalType_id")
                                         )
                                       ),
                                       uiOutput("rendScoreMethodBox"),
                                       tabsetPanel(
                                         id = "GeneSetTabs",
                                         tabPanel("Gene Sets",
                                                  p(),
                                                  uiOutput("rendGeneSetCat_Select"),
                                                  uiOutput("rendGeneSetTable"),
                                                  value = 1
                                         ),
                                         tabPanel("Single Genes",
                                                  #radioButtons("RawOrSS","Survival Analysis By:",
                                                  #             choices = c("Raw Gene Expression","Rank Normalized"),
                                                  #             selected = "Raw Gene Expression", inline = T),
                                                  uiOutput("rendGeneGeneSetTable"),
                                                  value = 2
                                         ),
                                         tabPanel("User Gene Set",
                                                  p(),
                                                  radioButtons("UserGSoption","",choices = c("Gene Set Upload","Text Box Input"), inline = T),
                                                  uiOutput("renduserGeneSet"),
                                                  uiOutput("renduserGeneSetTextName"),
                                                  uiOutput("renduserGeneSetText"),
                                                  uiOutput("rendUserGeneSetTable"),
                                                  value = 3
                                         )
                                       ),
                                       uiOutput("rendViewGeneSetGenes"),
                                       uiOutput("rendGenesInGeneSetTab")
                              ),
                              
                              ##--Survival Parameters--##
                              
                              tabPanel("Risk Strat Parameters",
                                       p(),
                                       h4("Risk Stratification Plot Parameters"),
                                       fluidRow(
                                         column(6,
                                                numericInput("cutoffTime1","High-Risk Survival Time Cutoff:", value = 364, min = 0, step = 1),
                                                selectInput("survStatus1","Survival Status Below Cutoff:", choices = c("1","0","0/1"), selected = "1")
                                         ),
                                         column(6,
                                                numericInput("cutoffTime0","Low-Risk Survival Time Cutoff:", value = 365, min = 0, step = 1),
                                                selectInput("survStatus0","Survival Status Above Cutoff:", choices = c("1","0","0/1"), selected = "0")
                                         )
                                       )
                              ),
                              
                              ##--Figure Parameters--##
                              
                              tabPanel("Figure Parameters",
                                       p(),
                                       h4("Survival Plot Parameters"),
                                       fluidRow(
                                         column(4,
                                                uiOutput("rendSurvXaxis")
                                         ),
                                         column(8,
                                                uiOutput("rendSurvPlotTitle")
                                         )
                                       ),
                                       fluidRow(
                                         column(3,
                                                selectInput("SurvLegendPos","Legend Position",choices = c("right","left","top","bottom","none"))
                                         ),
                                         column(3,
                                                checkboxInput("ShowPval","Show P.Value",value = T)
                                         ),
                                         column(3,
                                                checkboxInput("ShowConfInt","Show Confidence Interval",value = F)
                                         ),
                                         column(3,
                                                checkboxInput("ShowMedSurvLine","Show Median Survival Line",value = F)
                                         )
                                       ),
                                       hr(),
                                       h4("Boxplot Parameters"),
                                       fluidRow(
                                         column(6,
                                                numericInput("boxplotFont","Boxplot Font Size:", value = 15, step = 1),
                                                selectInput("boxoptselec","Boxplot Stat Compare Method:", choices = c("none","wilcox.test","t.test","kruskal.test","anova")) 
                                         ),
                                         column(6,
                                                numericInput("boxplotDot", "Boxplot Dot Size:", value = 0.75, step = 0.25),
                                                selectInput("boxplotTextAngle","X-Axis Text Orientation",
                                                            choices = c("Horizontal (0 degrees)" = "0","Angled (45 degrees)" = "45","Vertical (90 degrees)" = "90","Stagger"))
                                         )
                                       ),
                                       hr(),
                                       h4("Heatmap Parameters"),
                                       selectInput("ClusterMethod", "Select Cluster Method",
                                                   choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid")),
                                       fluidRow(
                                         column(6,
                                                numericInput("heatmapFontR", "Heatmap Row Font Size:", value = 9, step = 1)
                                         ),
                                         column(6,
                                                numericInput("heatmapFontC", "Heatmap Column Font Size:", value = 10, step = 1)
                                         )
                                       ),
                                       selectInput("ColorPaletteHeat", "Select Color Palette:",
                                                   choices = c("Red/Blue" = "original",
                                                               "OmniBlueRed" = "OmniBlueRed",
                                                               "LightBlue/BlackRed" = "LightBlueBlackRed",
                                                               "Green/Black/Red" = "GreenBlackRed",
                                                               "Yellow/Green/Blue" = "YlGnBu","Inferno" = "Inferno",
                                                               "Viridis" = "Viridis","Plasma" = "Plasma",
                                                               "Reds" = "OrRd","Blues" = "PuBu","Greens" = "Greens")
                                       ),
                                       hr(),
                                       h4("Forest Plot Parameters"),
                                       numericInput("ForestFontSize","Font Size",value = 1),
                                       hr(),
                                       h4("Linearity Plot Parameters"),
                                       fluidRow(
                                         column(4,
                                                numericInput("linAxisFont","X/Y Axis Font Size",
                                                             value = 14, step = 1)
                                         ),
                                         column(4,
                                                numericInput("linTickFont","Axis Tick Font Size",
                                                             value = 10, step = 1)
                                         ),
                                         column(4,
                                                numericInput("linMainFont","Title Font Size",
                                                             value = 16, step = 1)
                                         )
                                       )
                              )
                            )
                          ),
                          
                          ####----Main Panel----####
                          
                          mainPanel(
                            tabsetPanel(
                              id = "SurvPanels",
                              
                              ####----Survival Analysis Tab----####
                              
                              tabPanel("Pathway Level Survival Analysis",
                                       tabsetPanel(
                                         id = "SurvPanelsMain",
                                         
                                         ##--Median Cut Point--##
                                         
                                         tabPanel("Median Cut-Point Survival",
                                                  p(),
                                                  htmlOutput("BINSurvDescrip", style = "font-size:14px;"),
                                                  hr(),
                                                  withSpinner(jqui_resizable(plotOutput("SplotBIN", width = "850px", height = "550px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldSplotBIN_SVG","Download as SVG"),
                                                    downloadButton("dnldSplotBIN_PDF","Download as PDF")
                                                  ),
                                                  hr(),
                                                  withSpinner(jqui_resizable(plotOutput("ssgseaBINDensity", width = "650px", height = "400px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldssgseaBINDensity_SVG","Download as SVG"),
                                                    downloadButton("dnldssgseaBINDensity_PDF","Download as PDF")
                                                  ),
                                                  hr(),
                                                  h4("Cox Hazard Regression Analysis Summary"),
                                                  fluidRow(
                                                    column(6,
                                                           uiOutput("rendBINHRtab"),
                                                           style = 'border-right: 0.5px solid lightgray',
                                                    ),
                                                    column(6,
                                                           verbatimTextOutput("MedianCutPSumm")
                                                    )
                                                  ),
                                                  value = 1),
                                         
                                         ##--Quaritle Cutpoints--##
                                         
                                         tabPanel("Quartile Survival",
                                                  p(),
                                                  htmlOutput("QuartSurvDescrip", style = "font-size:14px;"),
                                                  hr(),
                                                  withSpinner(jqui_resizable(plotOutput("Splot", width = "850px", height = "550px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldSplot_SVG","Download as SVG"),
                                                    downloadButton("dnldSplot_PDF","Download as PDF")
                                                  ),
                                                  hr(),
                                                  withSpinner(jqui_resizable(plotOutput("ssgseaQuartDensity", width = "650px", height = "400px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldssgseaQuartDensity_SVG","Download as SVG"),
                                                    downloadButton("dnldssgseaQuartDensity_PDF","Download as PDF")
                                                  ),
                                                  hr(),
                                                  h4("Cox Hazard Regression Analysis Summary"),
                                                  fluidRow(
                                                    column(6,
                                                           uiOutput("rendQuartHRtab"),
                                                           style = 'border-right: 0.5px solid lightgray',
                                                    ),
                                                    column(6,
                                                           verbatimTextOutput("QuartileCutPSumm")
                                                    )
                                                  ),
                                                  value = 2),
                                         
                                         ##--Optimal Cut Point--##
                                         
                                         tabPanel("Optimal Cut-Point Survival",
                                                  p(),
                                                  htmlOutput("CutPSurvDescrip", style = "font-size:14px;"),
                                                  hr(),
                                                  withSpinner(jqui_resizable(plotOutput("ScutPointPlot", width = "850px", height = "550px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldScutPointPlot_SVG","Download as SVG"),
                                                    downloadButton("dnldScutPointPlot_PDF","Download as PDF")
                                                  ),
                                                  hr(),
                                                  withSpinner(jqui_resizable(plotOutput("ssgseaCutPDensity", width = "650px", height = "400px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldssgseaCutPDensity_SVG","Download as SVG"),
                                                    downloadButton("dnldssgseaCutPDensity_PDF","Download as PDF")
                                                  ),
                                                  hr(),
                                                  h4("Cox Hazard Regression Analysis Summary"),
                                                  fluidRow(
                                                    column(6,
                                                           uiOutput("rendCutPointHRtab"),
                                                           style = 'border-right: 0.5px solid lightgray',
                                                    ),
                                                    column(6,
                                                           verbatimTextOutput("OptimalCutPSumm")
                                                    )
                                                  ),
                                                  value = 3),
                                         
                                         ##--Quantile Cutpoints--##
                                         
                                         tabPanel("Top/Bottom Cut-Point Survival",
                                                  p(),
                                                  htmlOutput("QuantSurvDescrip", style = "font-size:14px;"),
                                                  hr(),
                                                  withSpinner(jqui_resizable(plotOutput("SquantPlot", width = "850px", height = "550px")), type = 6),
                                                  numericInput("QuantPercent","Top/Bottom Cut-Point Quantile Cutoff (%)", value = 25, min = 0, max = 100),
                                                  fluidRow(
                                                    downloadButton("dnldSquantPlot_SVG","Download as SVG"),
                                                    downloadButton("dnldSquantPlot_PDF","Download as PDF")
                                                  ),
                                                  hr(),
                                                  withSpinner(jqui_resizable(plotOutput("ssgseaQuantDensity", width = "650px", height = "400px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldssgseaQuantDensity_SVG","Download as SVG"),
                                                    downloadButton("dnldssgseaQuantDensity_PDF","Download as PDF")
                                                  ),
                                                  hr(),
                                                  h4("Cox Hazard Regression Analysis Summary"),
                                                  fluidRow(
                                                    column(6,
                                                           uiOutput("rendQuantHRtab"),
                                                           style = 'border-right: 0.5px solid lightgray',
                                                    ),
                                                    column(6,
                                                           verbatimTextOutput("QuantileCutPSumm")
                                                    )
                                                  ),
                                                  value = 4),
                                         
                                         ##--User Cut Point--##
                                         
                                         tabPanel("User Cut-Point Survival",
                                                  p(),
                                                  htmlOutput("Quant2SurvDescrip", style = "font-size:14px;"),
                                                  hr(),
                                                  withSpinner(jqui_resizable(plotOutput("SquantPlot2", width = "850px", height = "550px")), type = 6),
                                                  numericInput("QuantPercent2","Above/Below User Quantile Cut-Point (%)", value = 25, min = 0, max = 100),
                                                  fluidRow(
                                                    downloadButton("dnldSquantPlot2_SVG","Download as SVG"),
                                                    downloadButton("dnldSquantPlot2_PDF","Download as PDF")
                                                  ),
                                                  hr(),
                                                  withSpinner(jqui_resizable(plotOutput("ssgseaQuant2Density", width = "650px", height = "400px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldssgseaQuant2Density_SVG","Download as SVG"),
                                                    downloadButton("dnldssgseaQuant2Density_PDF","Download as PDF")
                                                  ),
                                                  hr(),
                                                  h4("Cox Hazard Regression Analysis Summary"),
                                                  fluidRow(
                                                    column(6,
                                                           uiOutput("rendQuantHRtab2"),
                                                           style = 'border-right: 0.5px solid lightgray',
                                                    ),
                                                    column(6,
                                                           verbatimTextOutput("UserCutPSumm")
                                                    )
                                                  ),
                                                  value = 5)
                                       ),
                                       value = 1),
                              
                              ####----Univariate Survival----####
                              
                              tabPanel("Univariate Survival Analysis",
                                       p(),
                                       fluidRow(
                                         column(5,
                                                uiOutput("rendSurvivalFeatureSingle"),
                                                ## Allows all select inputs to be wide enough to read the contents
                                                tags$head(
                                                  tags$style(HTML('
                                                                  .selectize-input {
                                                                      white-space: nowrap;
                                                                  }
                                                                  .selectize-dropdown {
                                                                      width: 500px !important;
                                                                  }'
                                                  )
                                                  )
                                                ),
                                                fluidRow(
                                                  column(3,
                                                         checkboxInput("UniVarNAcheck","Remove NA/Unknown/Inf",value = T)
                                                  ),
                                                  column(3,
                                                         checkboxInput("UniVarContCheck","Continuous Feature",value = F)
                                                  ),
                                                  column(3,
                                                         uiOutput("rendUniVarContHiLoCheck")
                                                  )
                                                ),
                                                uiOutput("rendSurvFeatVariableUni"),
                                         ),
                                         column(7,
                                                #htmlOutput("UnivarSummExpl", style = "font-size:14px;"),
                                         )
                                       ),
                                       tabsetPanel(
                                         id = "UniVarPlots",
                                         
                                         ##--Survival Plot--##
                                         
                                         tabPanel("Survival Plot",
                                                  p(),
                                                  withSpinner(jqui_resizable(plotOutput("featSplot", width = "850px", height = "550px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldfeatSplot_SVG","Download as SVG"),
                                                    downloadButton("dnldfeatSplot_PDF","Download as PDF")
                                                  )
                                         ),
                                         
                                         ##--Coxh Tables--##
                                         
                                         tabPanel("Coxh Table",
                                                  p(),
                                                  fluidRow(
                                                    column(6,
                                                           div(withSpinner(tableOutput("SSingleFeatureHRtab"), type = 7, size = 0.5), style = "font-size:12px; width:500px; overflow-X: scroll")
                                                    ),
                                                    column(6,
                                                           verbatimTextOutput("UnivarSummary")
                                                    )
                                                  )
                                         ),
                                         
                                         ##--Forest Plot--##
                                         
                                         tabPanel("Forest Plot",
                                                  p(),
                                                  withSpinner(jqui_resizable(plotOutput("SinglevarForestPlot", width = "100%", height = "800px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldUniVarForestplot_SVG","Download as SVG"),
                                                    downloadButton("dnldUniVarForestplot_PDF","Download as PDF")
                                                  )
                                         ),
                                         
                                         ##--Linearity Check--##
                                         
                                         tabPanel("Linearity Check",
                                                  p(),
                                                  fluidRow(
                                                    column(3,
                                                           selectInput("ResidualTypeUni","Select Residual Type",
                                                                       choices = c("deviance", "martingale", "score", "schoenfeld", "dfbeta", "dfbetas", "scaledsch", "partial"))
                                                    ),
                                                    column(3,
                                                           selectInput("linPredict1", "X-axis Scale:",
                                                                       choices = c("linear.predictions","observation.id","time"))
                                                    )
                                                  ),
                                                  uiOutput("timewarnmessage1"),
                                                  withSpinner(jqui_resizable(plotOutput("UnivarLinearityPlot", width = "100%", height = "500px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldUniVarLinplot_SVG","Download as SVG"),
                                                    downloadButton("dnldUniVarLinplot_PDF","Download as PDF")
                                                  )
                                         )
                                       ),
                                       value = 2),
                              
                              ####----Multivariate Survival----####
                              
                              tabPanel("Multivariate Coxh Analysis",
                                       tabsetPanel(
                                         id = "multivariate",
                                         
                                         ####----Bivariate Additive Survival----####
                                         
                                         tabPanel("Bivariate Additive Survival Analysis",
                                                  p(),
                                                  fluidRow(
                                                    column(4,
                                                           uiOutput("rendSurvivalFeatureBi1"),
                                                           fluidRow(
                                                             column(4,
                                                                    checkboxInput("BiVarAddNAcheck1","Remove NA/Unknown/Inf",value = T)
                                                             ),
                                                             column(4,
                                                                    checkboxInput("BiVarAddContCheck1","Continuous Feature",value = F)
                                                             ),
                                                             column(4,
                                                                    uiOutput("rendBiVarAddContHiLoCheck1")
                                                             )
                                                           ),
                                                           uiOutput("rendSurvFeatVariableBi1")
                                                    ),
                                                    column(4,
                                                           uiOutput("rendSurvivalFeatureBi2"),
                                                           fluidRow(
                                                             column(4,
                                                                    checkboxInput("BiVarAddNAcheck2","Remove NA/Unknown/Inf",value = T)
                                                             ),
                                                             column(4,
                                                                    checkboxInput("BiVarAddContCheck2","Continuous Feature",value = F)
                                                             ),
                                                             column(4,
                                                                    uiOutput("rendBiVarAddContHiLoCheck2")
                                                             )
                                                           ),
                                                           uiOutput("rendSurvFeatVariableBi2")
                                                    ),
                                                    column(4,
                                                           htmlOutput("BivarAddSummExpl", style = "font-size:14px;")
                                                    )
                                                  ),
                                                  tabsetPanel(
                                                    id = "BiVarPlots",
                                                    
                                                    ##--Coxh Tables--##
                                                    
                                                    tabPanel("Cox HR Table",
                                                             p(),
                                                             fluidRow(
                                                               column(6,
                                                                      div(withSpinner(tableOutput("BiFeatureHRtab"), type = 7, size = 0.5), style = "font-size:12px; width:500px; overflow-X: scroll")
                                                               ),
                                                               column(6,
                                                                      verbatimTextOutput("bivarSummary"),
                                                                      fluidRow(
                                                                        column(6,
                                                                               verbatimTextOutput("bivarAnova1")
                                                                        ),
                                                                        column(6,
                                                                               verbatimTextOutput("bivarAnova2")
                                                                        )
                                                                      )
                                                               )
                                                             )
                                                    ),
                                                    
                                                    ##--Forest Plot--##
                                                    
                                                    tabPanel("Forest Plot",
                                                             p(),
                                                             withSpinner(jqui_resizable(plotOutput("BivarForestPlot", width = "100%", height = "800px")), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldBiVarAddForest_SVG","Download as SVG"),
                                                               downloadButton("dnldBiVarAddForest_PDF","Download as PDF")
                                                             )
                                                    ),
                                                    
                                                    ##--Linearity Check--##
                                                    
                                                    tabPanel("Linearity Check",
                                                             p(),
                                                             fluidRow(
                                                               column(3,
                                                                      selectInput("ResidualTypeBi","Select Residual Type",
                                                                                  choices = c("deviance", "martingale", "score", "schoenfeld", "dfbeta", "dfbetas", "scaledsch", "partial"))
                                                               ),
                                                               column(3,
                                                                      selectInput("linPredict2", "X-axis Scale:",
                                                                                  choices = c("linear.predictions","observation.id","time"))
                                                               )
                                                             ),
                                                             uiOutput("timewarnmessage2"),
                                                             withSpinner(jqui_resizable(plotOutput("BivarLinearityPlot", width = "100%", height = "500px")), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldBiVarAddLinplot_SVG","Download as SVG"),
                                                               downloadButton("dnldBiVarAddLinplot_PDF","Download as PDF")
                                                             )
                                                    )
                                                  )
                                         ),
                                         
                                         ####----Bivariate Interaction Survival----####
                                         
                                         tabPanel("Bivariate Interaction Survival Analysis",
                                                  p(),
                                                  fluidRow(
                                                    column(4,
                                                           uiOutput("rendSurvivalFeatureBi1Inter"),
                                                           fluidRow(
                                                             column(4,
                                                                    checkboxInput("BiVarIntNAcheck1","Remove NA/Unknown/Inf",value = T)
                                                             ),
                                                             column(4,
                                                                    checkboxInput("BiVarIntContCheck1","Continuous Feature",value = F)
                                                             ),
                                                             column(4,
                                                                    uiOutput("rendBiVarIntContHiLoCheck1")
                                                             )
                                                           ),
                                                           uiOutput("rendSurvFeatVariableBi1Inter")
                                                           
                                                    ),
                                                    column(4,
                                                           uiOutput("rendSurvivalFeatureBi2Inter"),
                                                           fluidRow(
                                                             column(4,
                                                                    checkboxInput("BiVarIntNAcheck2","Remove NA/Unknown/Inf",value = T)
                                                             ),
                                                             column(4,
                                                                    checkboxInput("BiVarIntContCheck2","Continuous Feature",value = F)
                                                             ),
                                                             column(4,
                                                                    uiOutput("rendBiVarIntContHiLoCheck2")
                                                             )
                                                           ),
                                                           uiOutput("rendSurvFeatVariableBi2Inter")
                                                    ),
                                                    column(4,
                                                           verbatimTextOutput("BivarIntSummExpl")
                                                    )
                                                  ),
                                                  tabsetPanel(
                                                    id = "BiVarInterTabs",
                                                    
                                                    ##--Survival Plot--##
                                                    
                                                    tabPanel("Survival Plot",
                                                             p(),
                                                             withSpinner(jqui_resizable(plotOutput("featSplotBi", width = "850px", height = "550px")), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldfeatSplotBi_SVG","Download as SVG"),
                                                               downloadButton("dnldfeatSplotBi_PDF","Download as PDF")
                                                             )
                                                    ),
                                                    
                                                    ##--Coxh Tables--##
                                                    
                                                    tabPanel("Cox HR Table",
                                                             p(),
                                                             fluidRow(
                                                               column(6,
                                                                      div(withSpinner(tableOutput("BiFeatureHRtabInter"), type = 7, size = 0.5), style = "font-size:12px; width:500px; overflow-X: scroll")
                                                               ),
                                                               column(6,
                                                                      verbatimTextOutput("bivarSummaryInter"),
                                                                      verbatimTextOutput("bivarAnovaInter1")
                                                               )
                                                             )
                                                    ),
                                                    
                                                    ##--Forest Plot--##
                                                    
                                                    tabPanel("Forest Plot",
                                                             p(),
                                                             withSpinner(jqui_resizable(plotOutput("BivarForestPlotInter", width = "100%", height = "800px")), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldBiVarIntForest_SVG","Download as SVG"),
                                                               downloadButton("dnldBiVarIntForest_PDF","Download as PDF")
                                                             )
                                                    ),
                                                    
                                                    ##--Linearity Check--##
                                                    
                                                    tabPanel("Linearity Check",
                                                             p(),
                                                             fluidRow(
                                                               column(3,
                                                                      selectInput("ResidualTypeInter","Select Residual Type",
                                                                                  choices = c("deviance", "martingale", "score", "schoenfeld", "dfbeta", "dfbetas", "scaledsch", "partial"))
                                                               ),
                                                               column(3,
                                                                      selectInput("linPredict3", "X-axis Scale:",
                                                                                  choices = c("linear.predictions","observation.id","time"))
                                                               )
                                                             ),
                                                             uiOutput("timewarnmessage3"),
                                                             withSpinner(jqui_resizable(plotOutput("BivarLinearityPlotInter", width = "100%", height = "500px")), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldBiVarIntLinplot_SVG","Download as SVG"),
                                                               downloadButton("dnldBiVarIntLinplot_PDF","Download as PDF")
                                                             )
                                                    )
                                                    
                                                  )
                                         ),
                                         
                                         ####----Multivariate Survival----####
                                         
                                         tabPanel("Multivariate Coxh Analysis",
                                                  p(),
                                                  fluidRow(
                                                    column(4,
                                                           uiOutput("rendSurvivalFeature"),
                                                           checkboxInput("MultiVarNAcheck","Remove NA/Unknown",value = T)
                                                    )
                                                  ),
                                                  tabsetPanel(
                                                    id = "multivartabstwo",
                                                    
                                                    ##--Coxh Tables--##
                                                    
                                                    tabPanel("Coxh Tables",
                                                             fluidRow(
                                                               column(6,
                                                                      h4("Coxh Hazard Ratio (Categorical)"),
                                                                      verbatimTextOutput("multivarSummaryCat"),
                                                                      div(withSpinner(tableOutput("SFeatureHRtabCat"), type = 7, size = 0.5), style = "font-size:12px; width:500px; overflow-X: scroll")
                                                               ),
                                                               column(6,
                                                                      h4("Coxh Hazard Ratio (Continuous)"),
                                                                      verbatimTextOutput("multivarSummaryCont"),
                                                                      div(withSpinner(tableOutput("SFeatureHRtabCont"), type = 7, size = 0.5), style = "font-size:12px; width:500px; overflow-X: scroll")
                                                               )
                                                             )
                                                    ),
                                                    
                                                    ##--Forest Plot--##
                                                    
                                                    tabPanel("Forest Plot",
                                                             p(),
                                                             withSpinner(jqui_resizable(plotOutput("MultivarForestPlot", width = "100%", height = "800px")), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldMultiVarForest_SVG","Download as SVG"),
                                                               downloadButton("dnldMultiVarForest_PDF","Download as PDF")
                                                             )
                                                    )
                                                  )
                                         )
                                       ),
                                       value = 3),
                              
                              ####----Data Exploration----####
                              
                              tabPanel("Data Exploration",
                                       tabsetPanel(
                                         id = "DataExploration",
                                         tabPanel("Download Survival Data",
                                                  p(),
                                                  uiOutput("rendMetaTableCols"),
                                                  uiOutput("rendMetaTable"),
                                                  fluidRow(
                                                    column(6,
                                                           uiOutput("DnldMetaButon")
                                                    ),
                                                    column(6,
                                                           uiOutput("DnldExprButon")
                                                    )
                                                  ),
                                                  value = 5),
                                         tabPanel("Score Density",
                                                  p(),
                                                  fluidRow(
                                                    numericInput("densityPercent","User Defined Percentile (Red)",value = 15, width = "200px"),
                                                    checkboxInput("QuartileLinesCheck","Show Quartile Lines (Blue)",value = T)
                                                  ),
                                                  withSpinner(jqui_resizable(plotOutput("ssgseaDensity", width = "100%", height = "500px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldssgseaDensity_SVG","Download as SVG"),
                                                    downloadButton("dnldssgseaDensity_PDF","Download as PDF")
                                                  ),
                                                  div(DT::dataTableOutput("ssgseaDensityTable"), style = "font-size:12px"),
                                                  downloadButton("dnldssgseaDensityTable","Download Table"),
                                                  value = 6),
                                         tabPanel("Feature Comparison",
                                                  p(),
                                                  fluidRow(
                                                    column(3,
                                                           uiOutput("rendScatterFeature")
                                                    ),
                                                    column(2,
                                                           radioButtons("ColorScatterChoice","Color Plot by:",choices = c("Feature","Single Color"))
                                                    ),
                                                    column(3,
                                                           uiOutput("rendScatterColor")
                                                    ),
                                                    column(2,
                                                           checkboxGroupInput("ScatterLog","", choices = c("Log x-axis","Log y-axis"))
                                                    ),
                                                  ),
                                                  withSpinner(jqui_resizable(plotlyOutput("FeatCompScatterPlot", width = "100%", height = "500px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldFeatCompScatter_SVG","Download as SVG"),
                                                    downloadButton("dnldFeatCompScatter_PDF","Download as PDF")
                                                  ),
                                                  p(),
                                                  div(DT::dataTableOutput("FeatCompScatterTable"), style = "font-size:12px"),
                                                  downloadButton("dnldFeatCompScatterTable","Download Table"),
                                                  value = 7),
                                         tabPanel("Risk Stratification",
                                                  tabsetPanel(
                                                    tabPanel("Risk Straification Boxplot",
                                                             p("Users may adjust risk-stratification cutoff parameters under the 'Risk Strat Parameters' tab on the sidebar panel."),
                                                             checkboxInput("SBoxLog", "Log Transform Score", value = T),
                                                             withSpinner(jqui_resizable(plotOutput("Sboxplot", width = "100%", height = "500px")), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldSboxplot_SVG","Download as SVG"),
                                                               downloadButton("dnldSboxplot_PDF","Download as PDF")
                                                             ),
                                                             div(DT::dataTableOutput("SboxplotTable"), style = "font-size:12px; height:450px; overflow-Y: scroll"),
                                                             p(),
                                                             downloadButton("dnldSBoxplotTab","Download Table")
                                                    ),
                                                    tabPanel("Risk Straification Heatmap",
                                                             p("Users may adjust risk-stratification cutoff parameters under the 'Risk Strat Parameters' tab on the sidebar panel."),
                                                             uiOutput("heatmap_error_message"),
                                                             withSpinner(jqui_resizable(plotOutput("Sheatmap", width = "100%", height = "2000px")), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldSheatmap_SVG","Download as SVG"),
                                                               downloadButton("dnldSheatmap_PDF","Download as PDF"),
                                                               downloadButton("dnldSheatmapexpr","Download Expression Matrix From Heatmap")
                                                             )
                                                    )
                                                  )
                                         ),
                                         tabPanel("Feature Stratification",
                                                  tabsetPanel(
                                                    tabPanel("Feature Boxplot",
                                                             p(),
                                                             fluidRow(
                                                               column(4,
                                                                      uiOutput("rendBoxplotFeature")
                                                               ),
                                                               column(3,
                                                                      checkboxInput("BoxPRemoveNA","Remove NA/Unknowns",value = T)
                                                               ),
                                                               column(3,
                                                                      checkboxInput("FBoxLog", "Log Transform Score", value = T)
                                                               )
                                                             ),
                                                             withSpinner(jqui_resizable(plotOutput("Featureboxplot", width = "100%", height = "500px")), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldFboxplot_SVG","Download as SVG"),
                                                               downloadButton("dnldFboxplot_PDF","Download as PDF")
                                                             ),
                                                             div(DT::dataTableOutput("FeatureboxplotTable"), style = "font-size:12px; height:450px; overflow-Y: scroll"),
                                                             p(),
                                                             downloadButton("dnldFeatureboxplotTab","Download Table")
                                                    ),
                                                    tabPanel("Feature Heatmap",
                                                             p(),
                                                             uiOutput("heatmap_error_message2"),
                                                             fluidRow(
                                                               column(4,
                                                                      uiOutput("rendHeatmapFeature")
                                                               ),
                                                               column(3,
                                                                      checkboxInput("HeatRemoveNA","Remove NA/Unknowns",value = T)
                                                               )
                                                             ),
                                                             withSpinner(jqui_resizable(plotOutput("FeatureHeatmap", width = "100%", height = "2000px")), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldFheatmap_SVG","Download as SVG"),
                                                               downloadButton("dnldFheatmap_PDF","Download as PDF"),
                                                               downloadButton("dnldFheatmapexpr","Download Expression Matrix From Heatmap")
                                                             )
                                                    )
                                                  )
                                         )
                                       ),
                                       value = 4)
                            )
                          )
                        )
                      )
             ),
             tabPanel("About",
                      fluidPage(
                        mainPanel(
                          tabPanel("Purpose and Methods",
                                   uiOutput("rendPurposeAndMethodsMD")
                          )
                        )
                      )
             )
  )

####----Server----####


server <- function(input, output, session) {
  
  
  ####----Render UI----####
  
  ##--Sample Selection--##
  
  ## Select sample type to subset samples by - only render if more than one sample type
  output$rendSampleTypeSelection <- renderUI({
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      
      SampleTypeChoices <- unique(meta[,metacol_sampletype])
      SampleTypeChoices <- c("All_Sample_Types",SampleTypeChoices)
      selectInput("SampleTypeSelection",paste("Select Sample Type (",metacol_sampletype,"):",sep = ""),
                  choices = SampleTypeChoices, selected = PreSelect_SamplyType)
      
    }
    
  })
  
  ## Select primary feature to look at - All not working yet
  output$rendFeatureSelection <- renderUI({
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      
      if (input$SampleTypeSelection == "All_Sample_Types") {
        
        FeatureChoices <- c("All_Features",metacol_sampletype,metacol_feature)
        selectInput("FeatureSelection","Select Feature:", choices = FeatureChoices, selected = PreSelect_Feature)
        
      }
      else if (input$SampleTypeSelection != "All_Sample_Types") {
        
        FeatureChoices <- c("All_Features",metacol_feature)
        selectInput("FeatureSelection","Select Feature:", choices = FeatureChoices, selected = PreSelect_Feature)
        
      }
      
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      
      FeatureChoices <- c("All_Features",metacol_feature)
      selectInput("FeatureSelection","Select Feature:", choices = FeatureChoices, selected = PreSelect_Feature)
      
    }
    
  })
  
  ## Select primary features condition to look at - All not working yet
  output$rendSubFeatureSelection <- renderUI({
    
    req(input$FeatureSelection)
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      SampleType <- input$SampleTypeSelection
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      SampleType <- "All_Sample_Types"
    }
    #SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    
    if (SampleType == "All_Sample_Types") {
      meta <- meta
    }
    if (SampleType != "All_Sample_Types") {
      meta <- meta[which(meta[,metacol_sampletype] == SampleType),]
    }
    
    if (Feature != "All_Features") {
      
      SubFeatureChoices <- unique(meta[,Feature])
      # Sort options, will put 1,TRUE,yes before 0,FASLE,no, so the 'positive' value is the initial selected - puts NA last
      SubFeatureChoices <- sort(SubFeatureChoices, decreasing = T, na.last = T)
      selectInput("subFeatureSelection","Feature Condition:", choices = SubFeatureChoices, selected = PreSelect_SubFeature)
      
    }
    
  })
  
  ##--Univariate--##
  
  output$rendSurvivalFeatureSingle <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          metacol_feature <- metacol_feature[-which(metacol_feature == input$FeatureSelection)]
        }
        metacol_feature <- c(metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
        selectInput("SingleSurvivalFeature","Select Feature:",
                    choices = metacol_feature, selected = PreSelect_SecondaryFeature)
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        
        SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        }
        selectInput("SingleSurvivalFeature","Select Feature:",
                    choices = SurvFeatChoices2, selected = PreSelect_SecondaryFeature)
        
      }
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
      if (input$FeatureSelection != "All_Features") {
        SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
      }
      selectInput("SingleSurvivalFeature","Select Feature:",
                  choices = SurvFeatChoices2, selected = PreSelect_SecondaryFeature)
    }
    
  })
  
  ##--Bivariate--##
  
  output$rendSurvivalFeatureBi1 <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SurvFeatChoices <- c(metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
    SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices <- SurvFeatChoices[-which(SurvFeatChoices == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi1","Select Feature 1:",
                    choices = SurvFeatChoices, selected = "MedianCutP")
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi1","Select Feature 1:",
                    choices = SurvFeatChoices2, selected = "MedianCutP")
        
      }
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (input$FeatureSelection != "All_Features") {
        SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
      }
      selectInput("SurvivalFeatureBi1","Select Feature 1:",
                  choices = SurvFeatChoices2, selected = "MedianCutP")
    }
    
  })
  
  output$rendSurvivalFeatureBi2 <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SurvFeatChoices <- c(metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
    SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices <- SurvFeatChoices[-which(SurvFeatChoices == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi2","Select Feature 2:",
                    choices = SurvFeatChoices, selected = PreSelect_SecondaryFeature)
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi2","Select Feature 2:",
                    choices = SurvFeatChoices2, selected = PreSelect_SecondaryFeature)
        
      }
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (input$FeatureSelection != "All_Features") {
        SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
      }
      selectInput("SurvivalFeatureBi2","Select Feature 2:",
                  choices = SurvFeatChoices2, selected = PreSelect_SecondaryFeature)
    }
    
  })
  
  output$rendSurvivalFeatureBi1Inter <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SurvFeatChoices <- c(metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
    SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices <- SurvFeatChoices[-which(SurvFeatChoices == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi1Inter","Select Feature 1:",
                    choices = SurvFeatChoices, selected = "MedianCutP")
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi1Inter","Select Feature 1:",
                    choices = SurvFeatChoices2, selected = "MedianCutP")
        
      }
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (input$FeatureSelection != "All_Features") {
        SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
      }
      selectInput("SurvivalFeatureBi1Inter","Select Feature 1:",
                  choices = SurvFeatChoices2, selected = "MedianCutP")
    }
    
  })
  
  output$rendSurvivalFeatureBi2Inter <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SurvFeatChoices <- c(metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
    SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices <- SurvFeatChoices[-which(SurvFeatChoices == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi2Inter","Select Feature 2:",
                    choices = SurvFeatChoices, selected = PreSelect_SecondaryFeature)
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi2Inter","Select Feature 2:",
                    choices = SurvFeatChoices2, selected = PreSelect_SecondaryFeature)
        
      }
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (input$FeatureSelection != "All_Features") {
        SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
      }
      selectInput("SurvivalFeatureBi2Inter","Select Feature 2:",
                  choices = SurvFeatChoices2, selected = PreSelect_SecondaryFeature)
    }
    
  })
  
  ##--Multivariate--##
  
  output$rendSurvivalFeature <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SurvFeatChoices <- c(metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
    SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices <- SurvFeatChoices[-which(SurvFeatChoices == input$FeatureSelection)]
        }
        selectInput("SurvivalFeature","Select Feature(s):",
                    choices = SurvFeatChoices, multiple = T, selected = "MedianCutP")
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        }
        selectInput("SurvivalFeature","Select Feature(s):",
                    choices = SurvFeatChoices2, multiple = T, selected = "MedianCutP")
        
      }
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (input$FeatureSelection != "All_Features") {
        SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
      }
      selectInput("SurvivalFeature","Select Feature(s):",
                  choices = SurvFeatChoices2, multiple = T, selected = "MedianCutP")
    }
    
  })
  
  ##--Reference Features--##
  
  output$rendSurvFeatVariableUni <- renderUI({
    
    Feature <- input$SingleSurvivalFeature
    metaSub <- ssGSEAmeta()
    
    if (input$UniVarContCheck == FALSE) {
      
      Var_choices <- unique(metaSub[,Feature])
      if (input$UniVarNAcheck == TRUE) {
        
        Var_choices <- Var_choices[which(is.na(Var_choices) == FALSE)]
        Var_choices <- Var_choices[which(Var_choices != "Inf")]
        Var_choices <- Var_choices[grep("unknown",Var_choices,ignore.case = T, invert = T)]
        
      }
      Var_choices <- sort(Var_choices, decreasing = T, na.last = T)
      selectInput("SurvFeatVariableUni","Select Coxh Feature Reference:",
                  choices = Var_choices)
      
    }
    else if (input$UniVarContCheck == TRUE) {
      
      if (input$UniVarContHiLoCheck == TRUE) {
        
        metaSub <- metaSub[which(metaSub[,Feature] != "Inf"),]
        metaSub[,Feature] <- highlow(as.numeric(metaSub[, which(colnames(metaSub) == Feature)]))
        Var_choices <- unique(metaSub[,Feature])
        if (input$UniVarNAcheck == TRUE) {
          
          Var_choices <- Var_choices[which(is.na(Var_choices) == FALSE)]
          Var_choices <- Var_choices[which(Var_choices != "Inf")]
          Var_choices <- Var_choices[grep("unknown",Var_choices,ignore.case = T, invert = T)]
          
        }
        Var_choices <- sort(Var_choices, decreasing = T, na.last = T)
        selectInput("SurvFeatVariableUni","Select Coxh Feature Reference:",
                    choices = Var_choices)
        
      }
      
    }
    
  })
  
  output$rendSurvFeatVariableBi1 <- renderUI({
    
    Feature <- input$SurvivalFeatureBi1
    metaSub <- ssGSEAmeta()
    
    if (input$BiVarAddContCheck1 == FALSE) {
      
      Var_choices <- unique(metaSub[,Feature])
      if (input$BiVarAddNAcheck1 == TRUE) {
        
        Var_choices <- Var_choices[which(is.na(Var_choices) == FALSE)]
        Var_choices <- Var_choices[which(Var_choices != "Inf")]
        Var_choices <- Var_choices[grep("unknown",Var_choices,ignore.case = T, invert = T)]
        
      }
      Var_choices <- sort(Var_choices, decreasing = T, na.last = T)
      selectInput("SurvFeatVariableBi1","Select Coxh Feature Reference:",
                  choices = Var_choices)
      
    }
    else if (input$BiVarAddContCheck1 == TRUE) {
      
      if (input$BiVarAddContHiLoCheck1 == TRUE) {
        
        metaSub <- metaSub[which(metaSub[,Feature] != "Inf"),]
        metaSub[,Feature] <- highlow(as.numeric(metaSub[, which(colnames(metaSub) == Feature)]))
        Var_choices <- unique(metaSub[,Feature])
        if (input$BiVarAddNAcheck1 == TRUE) {
          
          Var_choices <- Var_choices[which(is.na(Var_choices) == FALSE)]
          Var_choices <- Var_choices[which(Var_choices != "Inf")]
          Var_choices <- Var_choices[grep("unknown",Var_choices,ignore.case = T, invert = T)]
          
        }
        Var_choices <- sort(Var_choices, decreasing = T, na.last = T)
        selectInput("SurvFeatVariableBi1","Select Coxh Feature Reference:",
                    choices = Var_choices)
        
      }
      
    }
    
  })
  
  output$rendSurvFeatVariableBi2 <- renderUI({
    
    Feature <- input$SurvivalFeatureBi2
    metaSub <- ssGSEAmeta()
    
    if (input$BiVarAddContCheck2 == FALSE) {
      
      Var_choices <- unique(metaSub[,Feature])
      if (input$BiVarAddNAcheck2 == TRUE) {
        
        Var_choices <- Var_choices[which(is.na(Var_choices) == FALSE)]
        Var_choices <- Var_choices[which(Var_choices != "Inf")]
        Var_choices <- Var_choices[grep("unknown",Var_choices,ignore.case = T, invert = T)]
        
      }
      Var_choices <- sort(Var_choices, decreasing = T, na.last = T)
      selectInput("SurvFeatVariableBi2","Select Coxh Feature Reference:",
                  choices = Var_choices)
      
    }
    else if (input$BiVarAddContCheck2 == TRUE) {
      
      if (input$BiVarAddContHiLoCheck2 == TRUE) {
        
        metaSub <- metaSub[which(metaSub[,Feature] != "Inf"),]
        metaSub[,Feature] <- highlow(as.numeric(metaSub[, which(colnames(metaSub) == Feature)]))
        Var_choices <- unique(metaSub[,Feature])
        if (input$BiVarAddNAcheck2 == TRUE) {
          
          Var_choices <- Var_choices[which(is.na(Var_choices) == FALSE)]
          Var_choices <- Var_choices[which(Var_choices != "Inf")]
          Var_choices <- Var_choices[grep("unknown",Var_choices,ignore.case = T, invert = T)]
          
        }
        Var_choices <- sort(Var_choices, decreasing = T, na.last = T)
        selectInput("SurvFeatVariableBi2","Select Coxh Feature Reference:",
                    choices = Var_choices)
        
      }
      
    }
    
  })
  
  output$rendSurvFeatVariableBi1Inter <- renderUI({
    
    Feature <- input$SurvivalFeatureBi1Inter
    metaSub <- ssGSEAmeta()
    
    if (input$BiVarIntContCheck1 == FALSE) {
      
      Var_choices <- unique(metaSub[,Feature])
      if (input$BiVarIntNAcheck1 == TRUE) {
        
        Var_choices <- Var_choices[which(is.na(Var_choices) == FALSE)]
        Var_choices <- Var_choices[which(Var_choices != "Inf")]
        Var_choices <- Var_choices[grep("unknown",Var_choices,ignore.case = T, invert = T)]
        
      }
      Var_choices <- sort(Var_choices, decreasing = T, na.last = T)
      selectInput("SurvFeatVariableBi1Inter","Select Coxh Feature Reference:",
                  choices = Var_choices)
      
    }
    else if (input$BiVarIntContCheck1 == TRUE) {
      
      if (input$BiVarIntContHiLoCheck1 == TRUE) {
        
        metaSub <- metaSub[which(metaSub[,Feature] != "Inf"),]
        metaSub[,Feature] <- highlow(as.numeric(metaSub[, which(colnames(metaSub) == Feature)]))
        Var_choices <- unique(metaSub[,Feature])
        if (input$BiVarIntNAcheck1 == TRUE) {
          
          Var_choices <- Var_choices[which(is.na(Var_choices) == FALSE)]
          Var_choices <- Var_choices[which(Var_choices != "Inf")]
          Var_choices <- Var_choices[grep("unknown",Var_choices,ignore.case = T, invert = T)]
          
        }
        Var_choices <- sort(Var_choices, decreasing = T, na.last = T)
        selectInput("SurvFeatVariableBi1Inter","Select Coxh Feature Reference:",
                    choices = Var_choices)
        
      }
      
    }
    
  })
  
  output$rendSurvFeatVariableBi2Inter <- renderUI({
    
    Feature <- input$SurvivalFeatureBi2Inter
    metaSub <- ssGSEAmeta()
    
    if (input$BiVarIntContCheck2 == FALSE) {
      
      Var_choices <- unique(metaSub[,Feature])
      if (input$BiVarIntNAcheck2 == TRUE) {
        
        Var_choices <- Var_choices[which(is.na(Var_choices) == FALSE)]
        Var_choices <- Var_choices[which(Var_choices != "Inf")]
        Var_choices <- Var_choices[grep("unknown",Var_choices,ignore.case = T, invert = T)]
        
      }
      Var_choices <- sort(Var_choices, decreasing = T, na.last = T)
      selectInput("SurvFeatVariableBi2Inter","Select Coxh Feature Reference:",
                  choices = Var_choices)
      
    }
    else if (input$BiVarIntContCheck2 == TRUE) {
      
      if (input$BiVarIntContHiLoCheck2 == TRUE) {
        
        metaSub <- metaSub[which(metaSub[,Feature] != "Inf"),]
        metaSub[,Feature] <- highlow(as.numeric(metaSub[, which(colnames(metaSub) == Feature)]))
        Var_choices <- unique(metaSub[,Feature])
        if (input$BiVarIntNAcheck2 == TRUE) {
          
          Var_choices <- Var_choices[which(is.na(Var_choices) == FALSE)]
          Var_choices <- Var_choices[which(Var_choices != "Inf")]
          Var_choices <- Var_choices[grep("unknown",Var_choices,ignore.case = T, invert = T)]
          
        }
        Var_choices <- sort(Var_choices, decreasing = T, na.last = T)
        selectInput("SurvFeatVariableBi2Inter","Select Coxh Feature Reference:",
                    choices = Var_choices)
        
      }
      
    }
    
  })
  
  ##--Other--##
  
  ## Survival X-axis Length
  output$rendSurvXaxis <- renderUI({
    
    meta_ssgsea <- ssGSEAmeta()
    if (is.null(input$SurvivalType_time) == TRUE & is.null(input$SurvivalType_id) == TRUE) {
      surv_time_col <- metacol_survtime[1]
      surv_id_col <- metacol_survid[1]
    }
    if (is.null(input$SurvivalType_time) == FALSE & is.null(input$SurvivalType_id) == FALSE) {
      surv_time_col <- input$SurvivalType_time
      surv_id_col <- input$SurvivalType_id
    }
    max_time <- ceiling(max(meta_ssgsea[,surv_time_col])/365.25)
    numericInput("SurvXaxis","X-Axis Limit (years)", value = max_time)
    
    
  })
  
  ## Select survival type selection based on meta columns (ex. OS vs EFS)
  output$rendSurvivalType_time <- renderUI({
    
    ## Only show if more than one option
    if (length(metacol_survtime > 1)) {
      
      selectInput("SurvivalType_time","Select Survival Time Data:", choices = metacol_survtime)
      
    }
    
  })
  
  ## Select survival type selection based on meta columns (ex. OS vs EFS)
  output$rendSurvivalType_id <- renderUI({
    
    ## Only show if more than one option
    if (length(metacol_survid > 1)) {
      
      selectInput("SurvivalType_id","Select Survival ID Data:", choices = metacol_survid)
      
    }
    
  })
  
  ## Feature selection for boxplot
  output$rendBoxplotFeature <- renderUI({
    
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        selectInput("BoxplotFeature","Select Feature:",
                    choices = metacol_feature, selected = PreSelect_SecondaryFeature)
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        selectInput("BoxplotFeature","Select Feature:",
                    choices = c(metacol_sampletype,metacol_feature), selected = PreSelect_SecondaryFeature)
        
      }
      
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      
      selectInput("BoxplotFeature","Select Feature:",
                  choices = metacol_feature, selected = PreSelect_SecondaryFeature)
      
    }
    
  })
  
  ## Feature selection for heatmap
  output$rendHeatmapFeature <- renderUI({
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        selectInput("HeatmapFeature","Select Feature:",
                    choices = metacol_feature)
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        selectInput("HeatmapFeature","Select Feature:",
                    choices = c(metacol_sampletype,metacol_feature))
        
      }
      
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      
      selectInput("HeatmapFeature","Select Feature:",
                  choices = metacol_feature)
      
    }
    
  })
  
  ## Select ssGSEA function scoring method
  output$rendScoreMethodBox <- renderUI({
    
    selectInput("ScoreMethod","Select Scoring Method",choices = c("ssgsea","gsva","zscore","plage"))
    
  })
  
  ## Survival Plot Title Input
  output$rendSurvPlotTitle <- renderUI({
    
    if (input$SurvPanels == 1) {
      if (input$SurvPanelsMain == 1) {
        textInput("SurvPlotTitleMedian","Survival Plot Title:",value = "")
      }
      else if (input$SurvPanelsMain == 2) {
        textInput("SurvPlotTitleQuartile","Survival Plot Title:",value = "")
      }
      else if (input$SurvPanelsMain == 3) {
        textInput("SurvPlotTitleOptimal","Survival Plot Title:",value = "")
      }
      else if (input$SurvPanelsMain == 4) {
        textInput("SurvPlotTitleQuantile","Survival Plot Title:",value = "")
      }
      else if (input$SurvPanelsMain == 5) {
        textInput("SurvPlotTitleUser","Survival Plot Title:",value = "")
      }
    }
    else if (input$SurvPanels == 2) {
      textInput("SurvPlotTitleUniVar","Survival Plot Title:",value = "")
    }
    else if (input$SurvPanels == 3) {
      textInput("SurvPlotTitleBiVar","Survival Plot Title:",value = "")
    }
    
    
    
  })
  
  ## Select column names of meta to view in sample box
  output$rendMetaTableCols <- renderUI({
    
    meta <- ssGSEAmeta()
    MetaColChoices <- colnames(meta)[c(2:ncol(meta))]
    selectInput("MetaTableCols","Select Meta Columns to View:", choices = MetaColChoices, selected = "", multiple = T)
    
  })
  
  ## View of meta table with choices selected
  output$rendMetaTable <- renderUI({
    
    div(DT::dataTableOutput("MetaTable"), style = "font-size:12px")
    
    
  })
  
  ## View Gene Set table
  output$rendGeneSetTable <- renderUI({
    
    div(DT::dataTableOutput("GeneSetTable"), style = "font-size:10px")
    
    
  })
  
  ## View Gene - "Gene Set" table
  output$rendGeneGeneSetTable <- renderUI({
    
    div(DT::dataTableOutput("geneGeneSetTable"), style = "font-size:10px")
    
  })
  
  ## View User Gene Set table
  output$rendUserGeneSetTable <- renderUI({
    
    req(input$userGeneSet)
    div(DT::dataTableOutput("userGeneSetTable"), style = "font-size:10px")
    
  })
  
  output$rendGenesInGeneSetTab <- renderUI({
    
    if (input$GeneSetTabs == 1) {
      
      if (input$ViewGeneSetGenes == TRUE) {
        
        div(DT::dataTableOutput("GenesInGeneSetTab"), style = "font-size:10px; height:450px; overflow-Y: scroll")
        
      }
      
    }
    
    else if (input$GeneSetTabs == 3) {
      
      req(input$userGeneSet)
      if (input$ViewGeneSetGenes == TRUE) {
        
        div(DT::dataTableOutput("GenesInGeneSetTab"), style = "font-size:10px; height:450px; overflow-Y: scroll")
        
      }
      
    }
    
  })
  
  output$rendViewGeneSetGenes <- renderUI({
    
    if (input$GeneSetTabs == 1) {
      
      checkboxInput("ViewGeneSetGenes","View Genes in Selected Gene Set", value = F)
      
    }
    
    else if (input$GeneSetTabs == 3) {
      
      req(input$userGeneSet)
      if (input$UserGSoption == "Gene Set Upload") {
        checkboxInput("ViewGeneSetGenes","View Genes in Selected Gene Set", value = F)
      }
      
    }
    
  })
  
  ## Download button for subset meta
  output$DnldMetaButon <- renderUI({
    
    downloadButton("dnldMeta", "Download Meta Subset")
    
  })
  
  ## Download button for subset expression
  output$DnldExprButon <- renderUI({
    
    downloadButton("dnldExpr", "Download Expression Subset")
    
  })
  
  output$rendQuartHRtab <- renderUI({
    
    div(withSpinner(tableOutput("SQuartileHRtab"), type = 7, size = 0.5), style = "font-size:12px")
    
  })
  output$rendBINHRtab <- renderUI({
    
    div(withSpinner(tableOutput("SBinaryHRtab"), type = 7, size = 0.5), style = "font-size:12px")
    
  })
  output$rendQuantHRtab <- renderUI({
    
    div(withSpinner(tableOutput("SQuantileHRtab"), type = 7, size = 0.5), style = "font-size:12px")
    
  })
  output$rendQuantHRtab2 <- renderUI({
    
    div(withSpinner(tableOutput("SQuantileHR2tab"), type = 7, size = 0.5), style = "font-size:12px")
    
  })
  output$rendCutPointHRtab <- renderUI({
    
    div(withSpinner(tableOutput("CutPointHRtab"), type = 7, size = 0.5), style = "font-size:12px")
    
  })
  
  output$rendPurposeAndMethodsMD <- renderUI({
    
    if (file.exists(About_MD_File)) {
      
      includeMarkdown(About_MD_File)
      
    }
    else {
      p("Page requires PurposeAndMethods.Rmd to be found")
    }
    
  })
  #output$rendTutorialMD <- renderUI({
  #  
  #  if (file.exists("App_Markdowns/Tutorial.Rmd")) {
  #    
  #    includeMarkdown("App_Markdowns/Tutorial.Rmd")
  #    
  #  }
  #  else {
  #    p("Page requires App_Markdowns/Tutorial.Rmd to be found")
  #  }
  #  
  #})
  
  output$rendGeneSetCat_Select <- renderUI({
    
    if (length(decon_score_cols) > 0) {
      selectInput("GeneSetCat_Select","Select Category",
                  choices = c("MSigDB","LINCS L1000","Cell Marker","ER Stress","Immune Signatures","Pre-Processed Scores"))
    }
    else {
      selectInput("GeneSetCat_Select","Select Category",
                  choices = c("MSigDB","LINCS L1000","Cell Marker","ER Stress","Immune Signatures"))
    }
    
  })
  
  output$renduserGeneSet <- renderUI({
    
    if (input$UserGSoption == "Gene Set Upload") {
      
      fileInput("userGeneSet","Gene Set Upload", accept = c(".gmt",".tsv",".txt",".RData"))
      
    }
    
  })
  
  output$renduserGeneSetText <- renderUI({
    
    if (input$UserGSoption == "Text Box Input") {
      
      textInput("userGeneSetText","Gene Symbols", placeholder = "Comma, space, or tab delimited")
      
    }
    
  })
  
  output$renduserGeneSetTextName <- renderUI({
    
    if (input$UserGSoption == "Text Box Input") {
      
      textInput("userGeneSetTextName","Custom Gene Set Name", value = "Custom Gene Set")
      
    }
    
  })
  
  output$rendUniVarContHiLoCheck <- renderUI({
    
    if (input$UniVarContCheck == T) {
      checkboxInput("UniVarContHiLoCheck","Continuous Feature as Median Cut-Point",value = T)
    }
    
  })
  output$rendBiVarAddContHiLoCheck1 <- renderUI({
    
    if (input$BiVarAddContCheck1 == T) {
      checkboxInput("BiVarAddContHiLoCheck1","Continuous Feature as Median Cut-Point",value = T)
    }
    
  })
  output$rendBiVarAddContHiLoCheck2 <- renderUI({
    
    if (input$BiVarAddContCheck2 == T) {
      checkboxInput("BiVarAddContHiLoCheck2","Continuous Feature as Median Cut-Point",value = T)
    }
    
  })
  output$rendBiVarIntContHiLoCheck1 <- renderUI({
    
    if (input$BiVarIntContCheck1 == T) {
      checkboxInput("BiVarIntContHiLoCheck1","Continuous Feature as Median Cut-Point",value = T)
    }
    
  })
  output$rendBiVarIntContHiLoCheck2 <- renderUI({
    
    if (input$BiVarIntContCheck2 == T) {
      checkboxInput("BiVarIntContHiLoCheck2","Continuous Feature as Median Cut-Point",value = T)
    }
    
  })
  
  output$rendScatterFeature <- renderUI({
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        selectInput("ScatterFeature","Select Feature:",
                    choices = metacol_feature, selected = PreSelect_SecondaryFeature)
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        selectInput("ScatterFeature","Select Feature:",
                    choices = c(metacol_sampletype,metacol_feature), selected = PreSelect_SecondaryFeature)
        
      }
      
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      
      selectInput("ScatterFeature","Select Feature:",
                  choices = metacol_feature, selected = PreSelect_SecondaryFeature)
      
    }
    
  })
  
  output$rendScatterColor <- renderUI({
    
    if (input$ColorScatterChoice == "Feature") {
      
      geneset <- gs_react()
      geneset_name <- names(geneset)
      
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        if (input$SampleTypeSelection != "All_Sample_Types") {
          
          if (input$FeatureSelection != "All_Features") {
            metacol_feature <- metacol_feature[-which(metacol_feature == input$FeatureSelection)]
          }
          metacol_feature <- c(metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
          metacol_feature <- c(metacol_feature,metacol_survid)
          selectInput("ScatterColor","Color By:",
                      choices = metacol_feature, selected = metacol_survid[1])
          
        }
        else if (input$SampleTypeSelection == "All_Sample_Types") {
          
          
          SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
          if (input$FeatureSelection != "All_Features") {
            SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
          }
          SurvFeatChoices2 <- c(SurvFeatChoices2,metacol_survid)
          selectInput("ScatterColor","Color By:",
                      choices = SurvFeatChoices2, selected = metacol_survid[1])
          
        }
      }
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        }
        SurvFeatChoices2 <- c(SurvFeatChoices2,metacol_survid)
        selectInput("ScatterColor","Color By:",
                    choices = SurvFeatChoices2, selected = metacol_survid[1])
      }
      
    }
    else if (input$ColorScatterChoice == "Single Color") {
      selectInput("ScatterColor","Color:",
                  choices = colors(), selected = "cadetblue")
    }
    
  })
  
  
  ####----Data Tables----####
  
  
  ## Render Meta Table
  output$MetaTable <- renderDataTable({
    
    # User selections
    SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    SubFeature <- input$subFeatureSelection
    meta <- ssGSEAmeta()
    if (is.null(input$SurvivalType_time) == TRUE & is.null(input$SurvivalType_id) == TRUE) {
      surv_time_col <- metacol_survtime[1]
      surv_id_col <- metacol_survid[1]
    }
    if (is.null(input$SurvivalType_time) == FALSE & is.null(input$SurvivalType_id) == FALSE) {
      surv_time_col <- input$SurvivalType_time
      surv_id_col <- input$SurvivalType_id
    }
    # Initial columns selected is NULL
    if (is.null(input$MetaTableCols) == TRUE) {
      
      userMetaCols <- NULL
      
    }
    # When meta tab is selected
    else if (is.null(input$MetaTableCols) == FALSE) {
      
      # If columns selected, add to list
      if (input$MetaTableCols != "") {
        userMetaCols <- input$MetaTableCols
      }
      else if (input$MetaTableCols == "") {
        userMetaCols <- NULL
      }
      
    }
    
    metaCols <- colnames(meta)[1] #select sample name column automatically
    if (Feature == "All_Features") {
      #userMetaCols <- userMetaCols[userMetaCols != Feature] #remove condition column from user selection because it is automatically added
      metaCols <- c(metaCols,surv_time_col,surv_id_col,userMetaCols) #combine column names selected
    }
    else if (Feature != "All_Features") {
      userMetaCols <- userMetaCols[userMetaCols != Feature] #remove condition column from user selection because it is automatically added
      metaCols <- c(metaCols,surv_time_col,surv_id_col,Feature,userMetaCols) #combine column names selected
    }
    meta_sub <- meta[,metaCols]
    DT::datatable(meta_sub,
                  extensions = "FixedColumns",
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 20,
                                 fixedColumns = list(leftColumns = 1),
                                 scrollX = T),
                  rownames = F)
    
  })
  
  output$GenesInGeneSetTab <- renderDataTable({
    
    if (input$GeneSetTabs == 1) {
      
      if (input$ViewGeneSetGenes == TRUE) {
        
        geneset <- gs_react()
        gs_df <- as.data.frame(geneset)
        DT::datatable(gs_df, options = list(paging = F), rownames = F)
        
      }
      
    }
    
    else if (input$GeneSetTabs == 3) {
      
      if (input$ViewGeneSetGenes == TRUE) {
        
        req(input$userGeneSet)
        geneset <- gs_react()
        gs_df <- as.data.frame(geneset)
        DT::datatable(gs_df, options = list(paging = F), rownames = F)
        
      }
      
    }
    
  })
  
  GeneSetTable_React <- reactive({
    
    GS_database <- input$GeneSetCat_Select
    sub_tab <- GeneSetTable_og[which(GeneSetTable_og[,1] == GS_database),]
    new_tab <- sub_tab[,-1]
    new_tab
    
  })
  
  ## Render Gene Set Selection Table
  output$GeneSetTable <- renderDataTable({
    
    new_GeneSetTable <- GeneSetTable_React()
    
    DT::datatable(new_GeneSetTable,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  selection = list(mode = 'single', selected = 1),
                  rownames = F)
    
  })
  
  ## Render Gene - "Gene Set" selection table
  output$geneGeneSetTable <- renderDataTable({
    
    DT::datatable(GeneGS_table,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  selection = list(mode = 'single', selected = 1),
                  rownames = F)
    
  })
  
  ## Render User Gene Set Table
  output$userGeneSetTable <- renderDataTable({
    
    req(input$userGeneSet)
    if (input$UserGSoption == "Gene Set Upload") {
      uGS_table <- userGeneSetTable_backend()
      DT::datatable(uGS_table,
                    options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                   pageLength = 10,
                                   scrollX = T),
                    selection = list(mode = 'single', selected = 1),
                    rownames = F)
    }
    
    
  })
  
  output$SboxplotTable <- renderDataTable({
    
    geneset <- gs_react()
    GeneSet <- names(geneset)
    ssGSEA_meta <- SboxplotReact()
    if (is.null(input$SurvivalType_time) == TRUE & is.null(input$SurvivalType_id) == TRUE) {
      surv_time_col <- metacol_survtime[1]
      surv_id_col <- metacol_survid[1]
    }
    if (is.null(input$SurvivalType_time) == FALSE & is.null(input$SurvivalType_id) == FALSE) {
      surv_time_col <- input$SurvivalType_time
      surv_id_col <- input$SurvivalType_id
    }
    boxTab <- ssGSEA_meta[,c("SampleName",surv_time_col,surv_id_col,GeneSet,"SurvivalCutoff")]
    DT::datatable(boxTab,
                  options = list(paging = F),
                  rownames = F)
    
  })
  
  output$ssgseaDensityTable <- renderDataTable({
    
    geneset <- gs_react()
    GeneSet <- names(geneset)
    ssgsea_meta <- ssGSEAmeta()
    table <- ssgsea_meta[,c("SampleName",GeneSet)]
    DT::datatable(table,
                  options = list(scrollY = T),
                  rownames = F)
    
  })
  
  ## Feature bloxplot table
  output$FeatureboxplotTable <- renderDataTable({
    
    geneset <- gs_react()
    GeneSet <- names(geneset)
    FeatureSelec <- input$BoxplotFeature
    meta_ssGSEA <- ssGSEAmeta()
    boxTab <- meta_ssGSEA[,c("SampleName",FeatureSelec,GeneSet)]
    if (input$BoxPRemoveNA == TRUE) {
      # Remove NA_unknown
      boxTab <- boxTab[which(is.na(boxTab[,FeatureSelec]) == FALSE),]
      boxTab <- boxTab[which(boxTab[,FeatureSelec] != "Inf"),]
      boxTab <- boxTab[grep("unknown",boxTab[,FeatureSelec],ignore.case = T, invert = T),]
    }
    DT::datatable(boxTab,
                  options = list(paging = F),
                  rownames = F)
    
  })
  
  ####----Reactives----####
  
  ## Render User Gene Set Table - Backend
  userGeneSetTable_backend <- reactive({
    
    gs.u <- input$userGeneSet
    ext <- tools::file_ext(gs.u$datapath)
    req(gs.u)
    validate(need(ext == c("gmt","tsv","txt", "RData"), "Please upload .gmt, .tsv, .txt, or .RData file"))
    
    # If user provides GMT file
    if (ext == "gmt") {
      gmt <- read.gmt(gs.u$datapath)
      uGS_table <- as.data.frame(unique(gmt[,1]))
      colnames(uGS_table)[1] <- "GeneSet"
    }
    
    # If user provides RData list file
    else if (ext == "RData") {
      gs_u <- loadRData(gs.u$datapath)
      uGS_table <- as.data.frame(names(gs_u))
      colnames(uGS_table)[1] <- "GeneSet"
    }
    
    # If user provides tab-delim two-col file
    else {
      gmt <- as.data.frame(read_delim(gs.u$datapath, delim = '\t'))
      uGS_table <- as.data.frame(unique(gmt[,1]))
      colnames(uGS_table)[1] <- "GeneSet"
    }
    uGS_table
    
  })
  
  user_gs_text <- reactive({
    
    if (input$GeneSetTabs == 3) {
      
      user_gs_text <- input$userGeneSetText
      user_gs_name <- input$userGeneSetTextName
      user_gs_name <- gsub("[[:punct:]]",".",user_gs_name)
      user_gs_name <- gsub(" ",".",user_gs_name)
      gs_text_s <- unlist(strsplit(user_gs_text, " "))
      gs_text_t <- unlist(strsplit(user_gs_text, "\t"))
      gs_text_c <- unlist(strsplit(user_gs_text, ","))
      
      gs_text <- unique(c(gs_text_s,gs_text_t,gs_text_c))
      
      gs_u_text <- list(gs_text)
      names(gs_u_text) <- user_gs_name
      gs_u_text
      
    }
    
  })
  
  ## User gene set list reactive
  user_gs <- reactive({
    
    gs.u <- input$userGeneSet
    ext <- tools::file_ext(gs.u$datapath)
    req(gs.u)
    validate(need(ext == c("gmt","tsv","txt", "RData"), "Please upload .gmt, .tsv, .txt, or .RData file"))
    
    # If user provides GMT file
    if (ext == "gmt") {
      gmt <- read.gmt(gs.u$datapath)
      colnames(gmt) <- c("term","gene")
      gs_u <- list()
      for (i in unique(gmt[,1])){
        gs_u[[i]] <- gmt[gmt[,1] == i,]$gene
      }
      gs_u
    }
    
    # If user provides RData list file
    else if (ext == "RData") {
      gs_u <- loadRData(gs.u$datapath)
    }
    
    # If user provides tab-delim two-col file
    else {
      gmt <- as.data.frame(read_delim(gs.u$datapath, delim = '\t'))
      colnames(gmt) <- c("term","gene")
      gs_u <- list()
      for (i in unique(gmt[,1])){
        gs_u[[i]] <- gmt[gmt[,1] == i,]$gene
      }
      gs_u
    }
    
    gs_u
    
  })
  
  ## Reactive to represent the chosen gene set
  gs_react <- reactive({
    
    GeneSetTable <- GeneSetTable_React()
    
    if (input$GeneSetTabs == 1) {
      if (gsTab == TRUE) {
        geneset_name <- GeneSetTable[input$GeneSetTable_rows_selected,3]
        if (geneset_name %in% decon_score_cols) {
          geneset <- list(meta[,geneset_name])
          names(geneset) <- geneset_name
        }
        else {
          geneset <- gs[geneset_name]
        }
      }
      if (gsTab == FALSE) {
        geneset_name <- GeneSetTable[input$GeneSetTable_rows_selected,1]
        if (geneset_name %in% decon_score_cols) {
          geneset <- list(meta[,geneset_name])
          names(geneset) <- geneset_name
        }
        else {
          geneset <- gs[geneset_name]
        }
      }
    }
    else if (input$GeneSetTabs == 2) {
      gene <- GeneGS_table[input$geneGeneSetTable_rows_selected,1]
      geneset <- list(gene = gene)
      names(geneset)[1] <- gene
    }
    else if (input$GeneSetTabs == 3) {
      if (input$UserGSoption == "Gene Set Upload") {
        req(input$userGeneSet)
        gs_u <- user_gs()
        uGS_table <- userGeneSetTable_backend()
        geneset <- gs_u[(uGS_table[input$userGeneSetTable_rows_selected,1])]
      }
      else if (input$UserGSoption == "Text Box Input") {
        req(input$userGeneSetText)
        geneset <- user_gs_text()
      }
      
    }
    
    geneset
    
    
  })
  
  ## Meta subset reactive - "All Primary features" not working yet
  metaSub <- reactive({
    
    req(input$FeatureSelection)
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      SampleType <- input$SampleTypeSelection
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      SampleType <- "All_Sample_Types"
    }
    Feature <- input$FeatureSelection
    SubFeature <- input$subFeatureSelection
    
    if (SampleType == "All_Sample_Types") {
      meta <- meta
    }
    if (SampleType != "All_Sample_Types") {
      meta <- meta[which(meta[,metacol_sampletype] == SampleType),]
    }
    
    if (Feature != "All_Features") {
      meta <- meta[which(meta[,Feature] == SubFeature),]
    }
    if (Feature == "All_Features") {
      meta <- meta
    }
    meta
    
  })
  
  ## Expression subset reactive
  exprSub <- reactive({
    
    meta <- metaSub()
    samples <- meta[,1]
    expr <- expr[,which(colnames(expr) %in% samples), drop = F] #Subset by sample names in subset meta table
    expr
    
  })
  
  ## Perform ssGSEA and functions on new meta table
  ssGSEAmeta <- reactive({
    
    meta <- metaSub()                     #Subset meta
    expr <- exprSub()                     #Subset expression
    geneset <- gs_react()                 #Chosen Gene Set
    geneset_name <- names(geneset)        #Name of chosen gene set
    quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
    quantCutoff2 <- input$QuantPercent2/100 #Quantile cutoff given by user
    surv_time_col <- input$SurvivalType_time
    surv_id_col <- input$SurvivalType_id
    ## Remove rows with NA in survival column
    meta <- meta[!is.na(meta[,surv_time_col]),]
    meta <- meta[!is.na(meta[,surv_id_col]),]
    
    ## Re-subset expression matrix
    samples <- meta[,1]
    expr_sub <- expr[,colnames(expr) %in% samples]
    expr_mat <- as.matrix(expr_sub)
    rownames(expr_mat) <- rownames(expr_sub)
    colnames(expr_mat) <- colnames(expr_sub)
    
    if (input$GeneSetTabs == 1) {
      if (geneset_name %in% decon_score_cols) {
        ssGSEA <- meta[,c("SampleName",geneset_name)]
        ssGSEA <- ssGSEA[!is.na(ssGSEA[,2]),]
      }
      else {
        ## Perform ssGSEA with gs and new subset data
        scoreMethod <- input$ScoreMethod      #ssGSEA score method chosen
        ssGSEA <- gsva(expr_mat,geneset,method = scoreMethod, verbose = FALSE)
        ## Transform
        ssGSEA <- as.data.frame(t(ssGSEA))
        ssGSEA$SampleName <- rownames(ssGSEA)
      }
    }
    else if (input$GeneSetTabs == 3) {
      if (input$UserGSoption == "Gene Set Upload") {
        ## Perform ssGSEA with gs and new subset data
        scoreMethod <- input$ScoreMethod      #ssGSEA score method chosen
        ssGSEA <- gsva(expr_mat,geneset,method = scoreMethod, verbose = FALSE)
        ## Transform
        ssGSEA <- as.data.frame(t(ssGSEA))
        ssGSEA$SampleName <- rownames(ssGSEA)
      }
      else if (input$UserGSoption == "Text Box Input") {
        #if (length())
        ## Perform ssGSEA with gs and new subset data
        scoreMethod <- input$ScoreMethod      #ssGSEA score method chosen
        ssGSEA <- gsva(expr_mat,geneset,method = scoreMethod, verbose = FALSE)
        ## Transform
        ssGSEA <- as.data.frame(t(ssGSEA))
        ssGSEA$SampleName <- rownames(ssGSEA)
      }
      
    }
    else if (input$GeneSetTabs == 2) {
      
      expr_sub2 <- expr_sub[geneset_name,]
      expr_sub3 <- as.data.frame(t(expr_sub2))
      expr_sub3$SampleName <- rownames(expr_sub3)
      ssGSEA <- expr_sub3
      
      #if (input$RawOrSS == "Raw Gene Expression") {
      #  
      #  expr_sub2 <- expr_sub[geneset_name,]
      #  expr_sub3 <- as.data.frame(t(expr_sub2))
      #  expr_sub3$SampleName <- rownames(expr_sub3)
      #  ssGSEA <- expr_sub3
      #  
      #}
      #else if (input$RawOrSS == "Rank Normalized") {
      #  
      #  ## Perform ssGSEA with gs and new subset data
      #  scoreMethod <- input$ScoreMethod      #ssGSEA score method chosen
      #  ssGSEA <- gsva(expr_mat,geneset,method = scoreMethod)
      #  ## Transform
      #  ssGSEA <- as.data.frame(t(ssGSEA))
      #  ssGSEA$SampleName <- rownames(ssGSEA)
      #  
      #}
    }
    
    #if (input$UserGSoption == "Text Box Input") {
    #  if (length(user_gs_text()) > 1) {
    #    
    #  }
    #}
    
    ##--Optimal cutpoint analysis--##
    
    ## Subset columns needed for plot and rename for surv function
    meta_ssgsea_sdf <- merge(meta[,c("SampleName",surv_time_col,surv_id_col)],ssGSEA[,c("SampleName",geneset_name)], by = "SampleName")
    #geneset_name <- gsub("[[:punct:]]","_",geneset_name)
    #colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
    
    if (length(meta_ssgsea_sdf[,4][meta_ssgsea_sdf[,4] > 0])/length(meta_ssgsea_sdf[,4]) > 0.01) {
      if (length(meta_ssgsea_sdf[,4]) > 1) {
        res.cut <- surv_cutpoint(meta_ssgsea_sdf,time = surv_time_col, event = surv_id_col, variable = geneset_name, minprop = 0.01)
        cutp <- res.cut$cutpoint[["cutpoint"]]
        res.cat <- surv_categorize(res.cut)
        ssGSEA$OptimalCutP <- res.cat[,3]
      }
    }
    
    
    ## Perform further functions
    ssGSEA$VAR_Q <- quartile_conversion(ssGSEA[, which(colnames(ssGSEA) == geneset_name)])
    ssGSEA$QuartileCutP <- paste("", ssGSEA$VAR_Q, sep="")
    ssGSEA$MedianCutP <- highlow(ssGSEA[, which(colnames(ssGSEA) == geneset_name)])
    ssGSEA$TopBottomCutP <- quantile_conversion(ssGSEA[, which(colnames(ssGSEA) == geneset_name)], quantCutoff)
    ssGSEA$UserCutP <- quantile_conversion2(ssGSEA[, which(colnames(ssGSEA) == geneset_name)], quantCutoff2)
    
    ## Merge with meta
    if (geneset_name %in% decon_score_cols) {
      ssGSEA <- ssGSEA[,-2] #remove score column so on merge the column does not duplicate
    }
    meta_ssGSEA <- merge(meta,ssGSEA, by = "SampleName", all = T)
    meta_ssGSEA
    
  })
  
  
  
  ####----Median Cut Point----####
  
  MedianCutP_react <- reactive({
    
    ## Assign variables
    surv_time_col <- input$SurvivalType_time
    surv_id_col <- input$SurvivalType_id
    meta_ssgsea <- ssGSEAmeta()
    
    ## Subset columns needed for plot
    meta_ssgsea_sdf <- meta_ssgsea[,c("SampleName",surv_time_col,surv_id_col,"MedianCutP")]
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "ID"
    
    meta_ssgsea_sdf
    
  })
  
  MedianCutPTab_react <- reactive({
    
    meta_ssgsea_sdf <- MedianCutP_react()
    surv_time_col <- input$SurvivalType_time
    surv_id_col <- input$SurvivalType_id
    geneset <- gs_react()
    geneset_name <- names(geneset)
    
    meta_ssgsea_sdf[,"MedianCutP"] <- as.factor(meta_ssgsea_sdf[,"MedianCutP"])
    meta_ssgsea_sdf[,"MedianCutP"] <- relevel(meta_ssgsea_sdf[,"MedianCutP"], ref = "Low")
    
    ## Survival Function
    tab <- coxph(Surv(time,ID) ~ MedianCutP, data = meta_ssgsea_sdf)
    tab
    
  })
  
  SBinaryHRtab_react <- reactive({
    
    tab <- MedianCutPTab_react()
    tab <- tab %>% 
      gtsummary::tbl_regression(exp = TRUE) %>%
      as_gt()
    
    tab_df <- as.data.frame(tab)
    
    tab_df <- tab_df %>%
      select(label,estimate,ci,p.value)
    colnames(tab_df) <- c("Characteristic","Hazard Ratio","95% Confidence Interval","P.Value")
    
    tab_df
    
  })
  
  SplotBIN_react <- reactive({
    
    ## Assign variables
    meta_ssgsea_sdf <- MedianCutP_react()
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    scoreMethod <- input$ScoreMethod
    show_pval <- input$ShowPval
    ShowConfInt <- input$ShowConfInt
    xaxlim <- input$SurvXaxis * 365.25
    surv_time_col <- input$SurvivalType_time
    showLegend <- input$SurvLegendPos
    showMedSurv <- input$ShowMedSurvLine
    if (showMedSurv == T) {
      showMedSurv <- "hv"
    }
    else if (showMedSurv == F) {
      showMedSurv <- "none"
    }
    
    ## Survival Function
    fit <- survfit(Surv(time,ID) ~ MedianCutP, data = meta_ssgsea_sdf, type="kaplan-meier")
    
    ## Determine type of survival data - OS/EFS/PFS?
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## determine Feature and Sample Type label
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (SampleType == "All_Sample_Types") {
        if (Feature == "All_Features") {
          SampleTypeLab <- "All Features in All Patients\n"
        }
        if (Feature != "All_Features") {
          SampleTypeLab <- paste(Feature," in All Patients\n")
        }
      }
      else {
        if (Feature == "All_Features") {
          SampleTypeLab <- paste("All Features (",SampleType,") Patients\n",sep = "")
        }
        if (Feature != "All_Features") {
          SampleTypeLab <- paste(Feature," (",SampleType,") Patients\n",sep = "")
        }
      }
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (Feature == "All_Features") {
        SampleTypeLab <- "All Features in All Patients\n"
      }
      if (Feature != "All_Features") {
        SampleTypeLab <- paste(Feature," in All Patients\n")
      }
    }
    
    
    ## Determine Scoring method label
    if (input$GeneSetTabs == 2) {
      scoreMethodLab <- "Gene Expression"
    }
    else if (input$GeneSetTabs == 1) {
      if (input$GeneSetCat_Select == "Pre-Processed Scores") {
        scoreMethodLab <- "Pre-Processed score"
      }
      else {
        scoreMethodLab <- paste(scoreMethod, " score", sep = "")
      }
    }
    else if (input$GeneSetTabs == 3) {
      scoreMethodLab <- paste(scoreMethod, " score", sep = "")
    }
    
    ## Determind Plot title
    if (is.null(input$SurvPlotTitleMedian)) {
      SurvPlotTitle <- paste("Survival curves of ",SampleTypeLab,
                             geneset_name," (",scoreMethodLab,")", sep = "")
    }
    else if (!is.null(input$SurvPlotTitleMedian)) {
      if (input$SurvPlotTitleMedian == "") {
        SurvPlotTitle <- paste("Survival curves of ",SampleTypeLab,
                               geneset_name," (",scoreMethodLab,")", sep = "")
      }
      else if (input$SurvPlotTitleMedian != "") {
        SurvPlotTitle <- input$SurvPlotTitleMedian
      }
    }
    
    
    
    ## Generate plot
    ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                         title = SurvPlotTitle,
                         xscale = c("d_y"),
                         break.time.by=365.25,
                         xlab = "Years", 
                         ylab = paste(SurvDateType,"Survival Probability"),
                         submain = "Based on Kaplan-Meier estimates",
                         caption = "created with survminer",
                         pval=show_pval,
                         conf.int = ShowConfInt,
                         ggtheme = theme_bw(),
                         font.title = c(16, "bold"),
                         font.submain = c(12, "italic"),
                         font.caption = c(12, "plain"),
                         font.x = c(14, "plain"),
                         font.y = c(14, "plain"),
                         font.tickslab = c(12, "plain"),
                         legend = showLegend,
                         risk.table.height = 0.20,
                         surv.median.line = showMedSurv
    )
    if (showMedSurv != "none") {
      MedSurvItem <- ggsurv[["plot"]][["layers"]][length(ggsurv[["plot"]][["layers"]])]
      MedSurvItem_df <- MedSurvItem[[1]][["data"]]
      MedSurvItem_df <- MedSurvItem_df[order(MedSurvItem_df[,1]),]
      MedSurvItem_df <- MedSurvItem_df %>%
        mutate(label = paste(round(MedSurvItem_df[,1]),"Days"))
      rownames(MedSurvItem_df) <- 1:nrow(MedSurvItem_df)
      if (nrow(MedSurvItem_df) > 1) {
        ggsurv$plot <- ggsurv$plot +
          geom_label_repel(data = MedSurvItem_df, aes(x = x1, y = y1, label = label, size = 4), label.size = NA, show.legend = FALSE)
      }
    }
    if (!is.null(input$SurvXaxis)) {
      ggsurv$plot$coordinates$limits$x <- c(0,xaxlim)
      ggsurv$table$coordinates$limits$x <- c(0,xaxlim)
    }
    
    ggsurv$table <- ggsurv$table + theme_cleantable()
    ggsurv
    
    
  })
  
  ssgseaBINDensity_react <- reactive({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    ssgsea_meta <- ssGSEAmeta()
    scoreMethod <- input$ScoreMethod
    CutPlabel <- "(Median Cut-Point)"
    
    cols_selec <- c("SampleName",geneset_name)
    ssgsea_scores <- ssgsea_meta[,cols_selec]
    
    quant_df2 <- data.frame(quantile(ssgsea_scores[,geneset_name],probs = 0.5,na.rm = T))
    #quant_df2 <- quant_df[c(2,3,4),,drop = F]
    colnames(quant_df2)[1] <- "Quantile"
    quant_df2$Quantile <- round(quant_df2$Quantile,3)
    
    #user_vline <- quantile(ssgsea_scores[,geneset_name],probs = 0.5)
    
    ## get score method for x and y labels
    if (input$GeneSetTabs == 2) {
      scoreMethodLab <- "Gene Expression Density"
      scoreMethodLab_x <- "Gene Expression (Log(exp+1))"
      ssgsea_scores[,geneset_name] <- log(ssgsea_scores[,geneset_name] + 1)
      quant_df2$Quantile <- round(log(quant_df2$Quantile + 1),3)
      
    }
    else if (input$GeneSetTabs != 2) {
      if (geneset_name %in% decon_score_cols) {
        scoreMethodLab <- "Pre-Processed Score Density"
        scoreMethodLab_x <- "Pre-Processed Score"
      }
      else {
        scoreMethodLab <- paste(scoreMethod, " Score Density", sep = "")
        scoreMethodLab_x <- paste(scoreMethod, " Score", sep = "")
      }
    }
    ## generate title based on input
    if (geneset_name %in% decon_score_cols) {
      dens_title <- paste(colnames(ssgsea_scores)[2]," Pre-Processed Score Density",sep = "")
    }
    else {
      if (input$GeneSetTabs == 2) {
        dens_title <- paste(colnames(ssgsea_scores)[2],"Log(exp+1)",scoreMethodLab)
      }
      else {
        dens_title <- paste(colnames(ssgsea_scores)[2],scoreMethodLab)
      }
    }
    
    dens_data <- density(ssgsea_scores[,geneset_name],na.rm = T)
    y_max <- max(dens_data$y)
    y_max_int <- y_max/6
    
    p <- ggplot(ssgsea_scores, aes(x=ssgsea_scores[,geneset_name])) + 
      geom_density(color="darkblue", fill="lightblue", alpha = 0.4) +
      xlab(scoreMethodLab_x) +
      ylab(scoreMethodLab) +
      ggtitle(dens_title) +
      theme(axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            plot.title = element_text(size = 20))
    p <- p + geom_vline(data = quant_df2, aes(xintercept = Quantile), linetype = "dashed", color = "darkblue", size = 1)
    p <- p + geom_text(aes(quant_df2[1,1],y_max-(y_max_int/6),label = paste(as.character(quant_df2[1,1]),CutPlabel),hjust = -0.1,vjust = -0.1),size = 6, check_overlap = T)
    p
    
  })
  
  output$BINSurvDescrip <- renderUI({
    
    ## Assign variables
    geneset <- gs_react()                      # Geneset Object
    geneset_name <- names(geneset)             # Geneset Name     
    scoreMethod <- input$ScoreMethod           # Scoring Method
    Feature <- input$FeatureSelection          # Feature Selected
    surv_time_col <- input$SurvivalType_time   # Survival Time Label
    HR_Tab <- SBinaryHRtab_react()             # Hazard Ratio Table Reactive
    ForPval <- MedianCutPTab_react()
    
    ## Survival Type
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## Determine Feature and sub feature
    if (Feature != "All_Features") {
      SubFeature <- input$subFeatureSelection
      Feature <- paste("<b>",Feature,"</b> - <b>",SubFeature,"</b></li>",sep = "")
      line2 <- paste("<li>The dataset is filtered by ",Feature,sep = "")
    }
    if (Feature == "All_Features") {
      line2 <- NULL
    }
    ## determine Sample Type
    if (is.null(metacol_sampletype) == T) {
      SampleType <- ""
    }
    else if (is.null(metacol_sampletype) == F) {
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        if (input$SampleTypeSelection == "All_Sample_Types") {
          SampleType <- " of all sample types.</li>"
        }
        else if (input$SampleTypeSelection != "All_Sample_Types") {
          SampleType <- paste(" of <b>",input$SampleTypeSelection,"</b> Patients.</li>",sep = "")
        }
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleType <- ".</li>"
      }
    }
    ## Determine Score method
    if (geneset_name %in% decon_score_cols) {
      scoreMethod <- " Pre-Processed Score"
    }
    if (input$GeneSetTabs == 2) {
      scoreMethod <- " Gene Expression"
      geneset_name <- paste(geneset_name,sep = "")
    }
    #pval <- get_lik_pval(ForPval)
    pval2 <- get_lik_pval(ForPval)
    HR <- HR_Tab[3,2]
    chacteristic <- HR_Tab[3,1]
    #if (pval == ">0.9") {
    #  pval2 <- 0.9
    #}
    #if (pval == "<0.001") {
    #  pval2 <- 0.001
    #}
    if (as.numeric(pval2) < 0.05) {
      pval_char <- "strongly associated"
    }
    if (as.numeric(pval2) >= 0.05 & as.numeric(pval2) < 0.1) {
      pval_char <- "moderately associated"
    }
    if (as.numeric(pval2) >= 0.1) {
      pval_char <- "not associated"
    }
    if (as.numeric(HR) > 1) {
      HR_char <- "high risk"
    }
    if (as.numeric(HR) <= 1) {
      HR_char <- "low risk"
    }
    
    line1 <- paste("<li><b>",SurvDateType,"</b> survival analysis ", SampleType,sep = "")
    line3 <- paste("<li>Kaplan-Meier survival curve dichotomized by <b>",geneset_name,"</b> <b>",scoreMethod,"</b> median cut-point.</li>",sep = "")
    line4 <- paste("<li>Cox hazard regression analysis finds a Likelihood Ratio P.value of <b>",pval2,"</b> and a Hazard Ratio of <b>",HR,"</b>, <b>",chacteristic,"</b> <b>",geneset_name,
                   "</b> is <b>",pval_char,"</b> with <b>",HR_char,"</b> for <b>",SurvDateType,"</b>.</li>",sep = "")
    
    if (is.null(line2)) {
      HTML(paste("<ul>",line1,line3,line4,"</ul>", sep = ""))
    }
    else if (!is.null(line2)) {
      HTML(paste("<ul>",line1,line2,line3,line4,"</ul>", sep = ""))
    }
    
  })
  
  output$MedianCutPSumm <- renderPrint({
    
    tab <- MedianCutPTab_react()
    out <- capture.output(summary(tab))
    
    con_line <- grep("^Concordance=",out,value = T)
    lik_line <- grep("^Likelihood ratio test=",out,value = T)
    wal_line <- grep("^Wald test",out,value = T)
    sco_line <- grep("^Score ",out,value = T)
    
    text <- paste("CoxH Summary:",con_line,lik_line,wal_line,sco_line,sep = "\n")
    cat(text)
    
  })
  
  output$SplotBIN <- renderPlot({
    
    plot <- SplotBIN_react()
    plot
    
  })
  
  output$SBinaryHRtab <- renderTable({
    
    tab <- SBinaryHRtab_react()
    tab
    
  })
  
  output$ssgseaBINDensity <- renderPlot({
    
    plot <- ssgseaBINDensity_react()
    plot
    
  })
  
  ####----Quartile Cut Point----####
  
  QuartileCutP_react <- reactive({
    
    ## Assign variables
    surv_time_col <- input$SurvivalType_time
    surv_id_col <- input$SurvivalType_id
    meta_ssgsea <- ssGSEAmeta()
    
    ## Subset columns needed for plot
    meta_ssgsea_sdf <- meta_ssgsea[,c("SampleName",surv_time_col,surv_id_col,"QuartileCutP")]
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "ID"
    
    meta_ssgsea_sdf
    
  })
  
  QuatileCutPTab_react <- reactive({
    
    meta_ssgsea_sdf <- QuartileCutP_react()
    surv_time_col <- input$SurvivalType_time
    surv_id_col <- input$SurvivalType_id
    geneset <- gs_react()
    geneset_name <- names(geneset)
    
    meta_ssgsea_sdf[,"QuartileCutP"] <- as.factor(meta_ssgsea_sdf[,"QuartileCutP"])
    meta_ssgsea_sdf[,"QuartileCutP"] <- relevel(meta_ssgsea_sdf[,"QuartileCutP"], ref = "Q1_Low")
    
    ## Survival Function
    tab <- coxph(Surv(time,ID) ~ QuartileCutP, data = meta_ssgsea_sdf)
    tab
    
  })
  
  SQuartileHRtab_react <- reactive({
    
    tab <- QuatileCutPTab_react()
    tab <- tab %>% 
      gtsummary::tbl_regression(exp = TRUE) %>%
      as_gt()
    
    tab_df <- as.data.frame(tab)
    
    tab_df <- tab_df %>%
      select(label,estimate,ci,p.value)
    colnames(tab_df) <- c("Characteristic","Hazard Ratio","95% Confidence Interval","P.Value")
    
    tab_df
    
  })
  
  Splot_react <- reactive({
    
    ## Assign variables
    meta_ssgsea_sdf <- QuartileCutP_react()
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    scoreMethod <- input$ScoreMethod
    show_pval <- input$ShowPval
    ShowConfInt <- input$ShowConfInt
    xaxlim <- input$SurvXaxis * 365.25
    surv_time_col <- input$SurvivalType_time
    showLegend <- input$SurvLegendPos
    showMedSurv <- input$ShowMedSurvLine
    if (showMedSurv == T) {
      showMedSurv <- "hv"
    }
    else if (showMedSurv == F) {
      showMedSurv <- "none"
    }
    
    ## Survival Function
    fit <- survfit(Surv(time,ID) ~ QuartileCutP, data = meta_ssgsea_sdf, type="kaplan-meier")
    
    ## Determine type of survival data - OS/EFS/PFS?
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## determine Feature and Sample Type label
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (SampleType == "All_Sample_Types") {
        if (Feature == "All_Features") {
          SampleTypeLab <- "All Features in All Patients\n"
        }
        if (Feature != "All_Features") {
          SampleTypeLab <- paste(Feature," in All Patients\n")
        }
      }
      else {
        if (Feature == "All_Features") {
          SampleTypeLab <- paste("All Features (",SampleType,") Patients\n",sep = "")
        }
        if (Feature != "All_Features") {
          SampleTypeLab <- paste(Feature," (",SampleType,") Patients\n",sep = "")
        }
      }
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (Feature == "All_Features") {
        SampleTypeLab <- "All Features in All Patients\n"
      }
      if (Feature != "All_Features") {
        SampleTypeLab <- paste(Feature," in All Patients\n")
      }
    }
    
    
    ## Determine Scoring method label
    if (input$GeneSetTabs == 2) {
      scoreMethodLab <- "Gene Expression"
    }
    else if (input$GeneSetTabs == 1) {
      if (input$GeneSetCat_Select == "Pre-Processed Scores") {
        scoreMethodLab <- "Pre-Processed score"
      }
      else {
        scoreMethodLab <- paste(scoreMethod, " score", sep = "")
      }
    }
    else if (input$GeneSetTabs == 3) {
      scoreMethodLab <- paste(scoreMethod, " score", sep = "")
    }
    
    ## Determind Plot title
    if (is.null(input$SurvPlotTitleQuartile)) {
      SurvPlotTitle <- paste("Survival curves of ",SampleTypeLab,
                             geneset_name," (",scoreMethodLab,")", sep = "")
    }
    else if (!is.null(input$SurvPlotTitleQuartile)) {
      if (input$SurvPlotTitleQuartile == "") {
        SurvPlotTitle <- paste("Survival curves of ",SampleTypeLab,
                               geneset_name," (",scoreMethodLab,")", sep = "")
      }
      else if (input$SurvPlotTitleQuartile != "") {
        SurvPlotTitle <- input$SurvPlotTitleQuartile
      }
    }
    
    
    ## Generate plot
    ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                         title = SurvPlotTitle,
                         xscale = c("d_y"),
                         break.time.by=365.25,
                         xlab = "Years", 
                         ylab = paste(SurvDateType,"Survival Probability"),
                         submain = "Based on Kaplan-Meier estimates",
                         caption = "created with survminer",
                         pval=show_pval,
                         conf.int = ShowConfInt,
                         ggtheme = theme_bw(),
                         font.title = c(16, "bold"),
                         font.submain = c(12, "italic"),
                         font.caption = c(12, "plain"),
                         font.x = c(14, "plain"),
                         font.y = c(14, "plain"),
                         font.tickslab = c(12, "plain"),
                         legend = showLegend,
                         risk.table.height = 0.20,
                         surv.median.line = showMedSurv
    )
    if (showMedSurv != "none") {
      MedSurvItem <- ggsurv[["plot"]][["layers"]][length(ggsurv[["plot"]][["layers"]])]
      MedSurvItem_df <- MedSurvItem[[1]][["data"]]
      MedSurvItem_df <- MedSurvItem_df[order(MedSurvItem_df[,1]),]
      MedSurvItem_df <- MedSurvItem_df %>%
        mutate(label = paste(round(MedSurvItem_df[,1]),"Days"))
      rownames(MedSurvItem_df) <- 1:nrow(MedSurvItem_df)
      if (nrow(MedSurvItem_df) > 1) {
        ggsurv$plot <- ggsurv$plot +
          geom_label_repel(data = MedSurvItem_df, aes(x = x1, y = y1, label = label, size = 4), label.size = NA, show.legend = FALSE)
      }
    }
    if (!is.null(input$SurvXaxis)) {
      ggsurv$plot$coordinates$limits$x <- c(0,xaxlim)
      ggsurv$table$coordinates$limits$x <- c(0,xaxlim)
    }
    
    ggsurv$table <- ggsurv$table + theme_cleantable()
    ggsurv
    
    
  })
  
  ssgseaQuartDensity_react <- reactive({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    ssgsea_meta <- ssGSEAmeta()
    scoreMethod <- input$ScoreMethod
    CutPlabel1 <- "(Q1 Cut-Point)"
    CutPlabel2 <- "(Q2 Cut-Point)"
    CutPlabel3 <- "(Q3 Cut-Point)"
    
    cols_selec <- c("SampleName",geneset_name)
    ssgsea_scores <- ssgsea_meta[,cols_selec]
    
    quant_df <- data.frame(quantile(ssgsea_scores[,geneset_name],na.rm = T))
    quant_df2 <- quant_df[c(2,3,4),,drop = F]
    colnames(quant_df2)[1] <- "Quantile"
    quant_df2$Quantile <- round(quant_df2$Quantile,3)
    
    #user_vline <- quantile(ssgsea_scores[,geneset_name],probs = 0.5)
    
    ## get score method for x and y labels
    if (input$GeneSetTabs == 2) {
      scoreMethodLab <- "Gene Expression Density"
      scoreMethodLab_x <- "Gene Expression (Log(exp+1))"
      ssgsea_scores[,geneset_name] <- log(ssgsea_scores[,geneset_name] + 1)
      quant_df2$Quantile <- round(log(quant_df2$Quantile + 1),3)
      
    }
    else if (input$GeneSetTabs != 2) {
      if (geneset_name %in% decon_score_cols) {
        scoreMethodLab <- "Pre-Processed Score Density"
        scoreMethodLab_x <- "Pre-Processed Score"
      }
      else {
        scoreMethodLab <- paste(scoreMethod, " Score Density", sep = "")
        scoreMethodLab_x <- paste(scoreMethod, " Score", sep = "")
      }
    }
    ## generate title based on input
    if (geneset_name %in% decon_score_cols) {
      dens_title <- paste(colnames(ssgsea_scores)[2]," Pre-Processed Score Density",sep = "")
    }
    else {
      if (input$GeneSetTabs == 2) {
        dens_title <- paste(colnames(ssgsea_scores)[2],"Log(exp+1)",scoreMethodLab)
      }
      else {
        dens_title <- paste(colnames(ssgsea_scores)[2],scoreMethodLab)
      }
    }
    
    dens_data <- density(ssgsea_scores[,geneset_name],na.rm = T)
    y_max <- max(dens_data$y)
    y_max_int <- y_max/6
    
    p <- ggplot(ssgsea_scores, aes(x=ssgsea_scores[,geneset_name])) + 
      geom_density(color="darkblue", fill="lightblue", alpha = 0.4) +
      xlab(scoreMethodLab_x) +
      ylab(scoreMethodLab) +
      ggtitle(dens_title) +
      theme(axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            plot.title = element_text(size = 20))
    p <- p + geom_vline(data = quant_df2, aes(xintercept = Quantile), linetype = "dashed", color = "darkblue", size = 1)
    p <- p + geom_text(aes(quant_df2[1,1],y_max-(y_max_int/6),label = paste(as.character(quant_df2[1,1]),CutPlabel1,sep="\n"),hjust = -0.1,vjust = 0.5),size = 6, check_overlap = T)
    p <- p + geom_text(aes(quant_df2[2,1],y_max-y_max_int,label = paste(as.character(quant_df2[2,1]),CutPlabel2,sep="\n"),hjust = -0.1,vjust = 0.5),size = 6, check_overlap = T)
    p <- p + geom_text(aes(quant_df2[3,1],y_max-(y_max_int*2),label = paste(as.character(quant_df2[3,1]),CutPlabel3,sep="\n"),hjust = -0.1,vjust = 0.5),size = 6, check_overlap = T)
    p
    
  })
  
  output$QuartSurvDescrip <- renderUI({
    
    ## Assign variables
    geneset <- gs_react()                      # Geneset Object
    geneset_name <- names(geneset)             # Geneset Name     
    scoreMethod <- input$ScoreMethod           # Scoring Method
    Feature <- input$FeatureSelection          # Feature Selected
    surv_time_col <- input$SurvivalType_time   # Survival Time Label
    HR_Tab <- SQuartileHRtab_react()             # Hazard Ratio Table Reactive
    ForPval <- QuatileCutPTab_react()
    
    ## Survival Type
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## Determine Feature and sub feature
    if (Feature != "All_Features") {
      SubFeature <- input$subFeatureSelection
      Feature <- paste("<b>",Feature,"</b> - <b>",SubFeature,"</b></li>",sep = "")
      line2 <- paste("<li>The dataset is filtered by ",Feature,sep = "")
    }
    if (Feature == "All_Features") {
      line2 <- NULL
    }
    ## determine Sample Type
    if (is.null(metacol_sampletype) == T) {
      SampleType <- ""
    }
    else if (is.null(metacol_sampletype) == F) {
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        if (input$SampleTypeSelection == "All_Sample_Types") {
          SampleType <- " of all sample types.</li>"
        }
        else if (input$SampleTypeSelection != "All_Sample_Types") {
          SampleType <- paste(" of <b>",input$SampleTypeSelection,"</b> Patients.</li>",sep = "")
        }
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleType <- ".</li>"
      }
    }
    ## Determine Score method
    if (geneset_name %in% decon_score_cols) {
      scoreMethod <- " Pre-Processed Score"
    }
    if (input$GeneSetTabs == 2) {
      scoreMethod <- " Gene Expression"
      geneset_name <- paste(geneset_name,sep = "")
    }
    #pval <- get_lik_pval(ForPval)
    pval2 <- get_lik_pval(ForPval)
    #HR <- HR_Tab[3,2]
    #chacteristic <- HR_Tab[5,1]
    #if (pval == ">0.9") {
    #  pval2 <- 0.9
    #}
    #if (pval == "<0.001") {
    #  pval2 <- 0.001
    #}
    if (as.numeric(pval2) < 0.05) {
      pval_char <- "strongly associated"
    }
    if (as.numeric(pval2) >= 0.05 & as.numeric(pval2) < 0.1) {
      pval_char <- "moderately associated"
    }
    if (as.numeric(pval2) >= 0.1) {
      pval_char <- "not associated"
    }
    #if (as.numeric(HR) > 1) {
    #  HR_char <- "high risk"
    #}
    #if (as.numeric(HR) <= 1) {
    #  HR_char <- "low risk"
    #}
    
    line1 <- paste("<li><b>",SurvDateType,"</b> survival analysis ", SampleType,sep = "")
    line3 <- paste("<li>Kaplan-Meier survival curve categorized by <b>",geneset_name,"</b> <b>",scoreMethod,"</b> quartile cut-points.</li>",sep = "")
    #line4 <- paste("<li>Analysis finds a Likelihood Ratio P.value of <b>",pval,"</b> and a Hazard Ratio of <b>",HR,"</b>, <b>",chacteristic,"</b> <b>",geneset_name,
    #               "</b> is <b>",pval_char,"</b> with <b>",HR_char,"</b> for <b>",SurvDateType,"</b>.</li>",sep = "")
    line4 <- paste("<li>Cox hazard regression analysis finds a Likelihood Ratio P.value of <b>",pval2,"</b> shows that <b>",geneset_name,
                   "</b> is <b>",pval_char,"</b> with <b>",SurvDateType,"</b>.</li>",sep = "")
    
    if (is.null(line2)) {
      HTML(paste("<ul>",line1,line3,line4,"</ul>", sep = ""))
    }
    else if (!is.null(line2)) {
      HTML(paste("<ul>",line1,line2,line3,line4,"</ul>", sep = ""))
    }
    
  })
  
  output$QuartileCutPSumm <- renderPrint({
    
    tab <- QuatileCutPTab_react()
    out <- capture.output(summary(tab))
    
    con_line <- grep("^Concordance=",out,value = T)
    lik_line <- grep("^Likelihood ratio test=",out,value = T)
    wal_line <- grep("^Wald test",out,value = T)
    sco_line <- grep("^Score ",out,value = T)
    
    text <- paste("CoxH Summary:",con_line,lik_line,wal_line,sco_line,sep = "\n")
    cat(text)
    
    
  })
  
  output$Splot <- renderPlot({
    
    plot <- Splot_react()
    plot
    
  })
  
  output$SQuartileHRtab <- renderTable({
    
    tab <- SQuartileHRtab_react()
    tab
    
  })
  
  output$ssgseaQuartDensity <- renderPlot({
    
    plot <- ssgseaQuartDensity_react()
    plot
    
  })
  
  ####----Optimal Cut Point----####
  
  OptimalCutP_react <- reactive({
    
    ## Assign variables
    surv_time_col <- input$SurvivalType_time
    surv_id_col <- input$SurvivalType_id
    meta_ssgsea <- ssGSEAmeta()
    
    ## Subset columns needed for plot
    meta_ssgsea_sdf <- meta_ssgsea[,c("SampleName",surv_time_col,surv_id_col,"OptimalCutP")]
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "ID"
    
    meta_ssgsea_sdf
    
  })
  
  OptimalCutPTab_react <- reactive({
    
    meta_ssgsea_sdf <- OptimalCutP_react()
    surv_time_col <- input$SurvivalType_time
    surv_id_col <- input$SurvivalType_id
    geneset <- gs_react()
    geneset_name <- names(geneset)
    
    meta_ssgsea_sdf[,"OptimalCutP"] <- as.factor(meta_ssgsea_sdf[,"OptimalCutP"])
    meta_ssgsea_sdf[,"OptimalCutP"] <- relevel(meta_ssgsea_sdf[,"OptimalCutP"], ref = "low")
    
    ## Survival Function
    tab <- coxph(Surv(time,ID) ~ OptimalCutP, data = meta_ssgsea_sdf)
    tab
    
  })
  
  CutPointHRtab_react <- reactive({
    
    tab <- OptimalCutPTab_react()
    tab <- tab %>% 
      gtsummary::tbl_regression(exp = TRUE) %>%
      as_gt()
    
    tab_df <- as.data.frame(tab)
    
    tab_df <- tab_df %>%
      select(label,estimate,ci,p.value)
    colnames(tab_df) <- c("Characteristic","Hazard Ratio","95% Confidence Interval","P.Value")
    
    tab_df
    
  })
  
  ScutPointPlot_react <- reactive({
    
    ## Assign variables
    meta_ssgsea_sdf <- OptimalCutP_react()
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    scoreMethod <- input$ScoreMethod
    show_pval <- input$ShowPval
    ShowConfInt <- input$ShowConfInt
    xaxlim <- input$SurvXaxis * 365.25
    surv_time_col <- input$SurvivalType_time
    showLegend <- input$SurvLegendPos
    showMedSurv <- input$ShowMedSurvLine
    if (showMedSurv == T) {
      showMedSurv <- "hv"
    }
    else if (showMedSurv == F) {
      showMedSurv <- "none"
    }
    
    
    ## Survival Function
    fit <- survfit(Surv(time,ID) ~ OptimalCutP, data = meta_ssgsea_sdf, type="kaplan-meier")
    
    ## Determine type of survival data - OS/EFS/PFS?
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## determine Feature and Sample Type label
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (SampleType == "All_Sample_Types") {
        if (Feature == "All_Features") {
          SampleTypeLab <- "All Features in All Patients\n"
        }
        if (Feature != "All_Features") {
          SampleTypeLab <- paste(Feature," in All Patients\n")
        }
      }
      else {
        if (Feature == "All_Features") {
          SampleTypeLab <- paste("All Features (",SampleType,") Patients\n",sep = "")
        }
        if (Feature != "All_Features") {
          SampleTypeLab <- paste(Feature," (",SampleType,") Patients\n",sep = "")
        }
      }
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (Feature == "All_Features") {
        SampleTypeLab <- "All Features in All Patients\n"
      }
      if (Feature != "All_Features") {
        SampleTypeLab <- paste(Feature," in All Patients\n")
      }
    }
    
    
    ## Determine Scoring method label
    if (input$GeneSetTabs == 2) {
      scoreMethodLab <- "Gene Expression"
    }
    else if (input$GeneSetTabs == 1) {
      if (input$GeneSetCat_Select == "Pre-Processed Scores") {
        scoreMethodLab <- "Pre-Processed score"
      }
      else {
        scoreMethodLab <- paste(scoreMethod, " score", sep = "")
      }
    }
    else if (input$GeneSetTabs == 3) {
      scoreMethodLab <- paste(scoreMethod, " score", sep = "")
    }
    
    ## Determind Plot title
    if (is.null(input$SurvPlotTitleOptimal)) {
      SurvPlotTitle <- paste("Survival curves of ",SampleTypeLab,
                             geneset_name," (",scoreMethodLab,")", sep = "")
    }
    else if (!is.null(input$SurvPlotTitleOptimal)) {
      if (input$SurvPlotTitleOptimal == "") {
        SurvPlotTitle <- paste("Survival curves of ",SampleTypeLab,
                               geneset_name," (",scoreMethodLab,")", sep = "")
      }
      else if (input$SurvPlotTitleOptimal != "") {
        SurvPlotTitle <- input$SurvPlotTitleOptimal
      }
    }
    
    
    ## Generate plot
    ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                         title = SurvPlotTitle,
                         xscale = c("d_y"),
                         break.time.by=365.25,
                         xlab = "Years", 
                         ylab = paste(SurvDateType,"Survival Probability"),
                         submain = "Based on Kaplan-Meier estimates",
                         caption = "created with survminer",
                         pval=show_pval,
                         conf.int = ShowConfInt,
                         ggtheme = theme_bw(),
                         font.title = c(16, "bold"),
                         font.submain = c(12, "italic"),
                         font.caption = c(12, "plain"),
                         font.x = c(14, "plain"),
                         font.y = c(14, "plain"),
                         font.tickslab = c(12, "plain"),
                         legend = showLegend,
                         risk.table.height = 0.20,
                         surv.median.line = showMedSurv
    )
    if (showMedSurv != "none") {
      MedSurvItem <- ggsurv[["plot"]][["layers"]][length(ggsurv[["plot"]][["layers"]])]
      MedSurvItem_df <- MedSurvItem[[1]][["data"]]
      MedSurvItem_df <- MedSurvItem_df[order(MedSurvItem_df[,1]),]
      MedSurvItem_df <- MedSurvItem_df %>%
        mutate(label = paste(round(MedSurvItem_df[,1]),"Days"))
      rownames(MedSurvItem_df) <- 1:nrow(MedSurvItem_df)
      if (nrow(MedSurvItem_df) > 1) {
        ggsurv$plot <- ggsurv$plot +
          geom_label_repel(data = MedSurvItem_df, aes(x = x1, y = y1, label = label, size = 4), label.size = NA, show.legend = FALSE)
      }
    }
    if (!is.null(input$SurvXaxis)) {
      ggsurv$plot$coordinates$limits$x <- c(0,xaxlim)
      ggsurv$table$coordinates$limits$x <- c(0,xaxlim)
    }
    
    ggsurv$table <- ggsurv$table + theme_cleantable()
    ggsurv
    
    
  })
  
  ssgseaCutPDensity_react <- reactive({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    ssgsea_meta <- ssGSEAmeta()
    scoreMethod <- input$ScoreMethod
    CutPlabel <- "(Optimal Cut-Point)"
    surv_time_col <- input$SurvivalType_time
    surv_id_col <- input$SurvivalType_id
    
    cols_selec <- c("SampleName",geneset_name)
    ssgsea_scores <- ssgsea_meta[,cols_selec]
    
    meta_ssgsea_sdf <- ssgsea_meta[,c("SampleName",surv_time_col,surv_id_col,geneset_name)]
    
    if (length(meta_ssgsea_sdf[,4][meta_ssgsea_sdf[,4] > 0])/length(meta_ssgsea_sdf[,4]) > 0.01) {
      res.cut <- surv_cutpoint(meta_ssgsea_sdf,time = surv_time_col, event = surv_id_col, variable = geneset_name, minprop = 0.01)
      cutp <- res.cut$cutpoint[["cutpoint"]]
      #res.cat <- surv_categorize(res.cut)
      #ssGSEA$OptimalCutP <- res.cat[,3]
      
      res.cut <- surv_cutpoint(meta_ssgsea_sdf,time = surv_time_col, event = surv_id_col, variable = geneset_name)
      cutp <- round(res.cut$cutpoint[["cutpoint"]],3)
      
      
      ## get score method for x and y labels
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "Gene Expression Density"
        scoreMethodLab_x <- "Gene Expression (Log(exp+1))"
        ssgsea_scores[,geneset_name] <- log(ssgsea_scores[,geneset_name] + 1)
        cutp <- round(log(cutp + 1),3)
        
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "Pre-Processed Score Density"
          scoreMethodLab_x <- "Pre-Processed Score"
        }
        else {
          scoreMethodLab <- paste(scoreMethod, " Score Density", sep = "")
          scoreMethodLab_x <- paste(scoreMethod, " Score", sep = "")
        }
      }
      ## generate title based on input
      if (geneset_name %in% decon_score_cols) {
        dens_title <- paste(colnames(ssgsea_scores)[2]," Pre-Processed Score Density",sep = "")
      }
      else {
        if (input$GeneSetTabs == 2) {
          dens_title <- paste(colnames(ssgsea_scores)[2],"Log(exp+1)",scoreMethodLab)
        }
        else {
          dens_title <- paste(colnames(ssgsea_scores)[2],scoreMethodLab)
        }
      }
      
      dens_data <- density(ssgsea_scores[,geneset_name],na.rm = T)
      y_max <- max(dens_data$y)
      y_max_int <- y_max/6
      
      p <- ggplot(ssgsea_scores, aes(x=ssgsea_scores[,geneset_name])) + 
        geom_density(color="darkblue", fill="lightblue", alpha = 0.4) +
        xlab(scoreMethodLab_x) +
        ylab(scoreMethodLab) +
        ggtitle(dens_title) +
        theme(axis.text = element_text(size = 14),
              axis.title = element_text(size = 16),
              plot.title = element_text(size = 20))
      p <- p + geom_vline(aes(xintercept = cutp), linetype = "dashed", color = "darkblue", size = 1)
      p <- p + geom_text(aes(cutp,y_max-(y_max_int/6),label = paste(as.character(cutp),CutPlabel),hjust = -0.1,vjust = -0.1),size = 6, check_overlap = T)
      p
      
    }
    
    
    
  })
  
  output$CutPSurvDescrip <- renderUI({
    
    ## Assign variables
    geneset <- gs_react()                      # Geneset Object
    geneset_name <- names(geneset)             # Geneset Name     
    scoreMethod <- input$ScoreMethod           # Scoring Method
    Feature <- input$FeatureSelection          # Feature Selected
    surv_time_col <- input$SurvivalType_time   # Survival Time Label
    HR_Tab <- CutPointHRtab_react()             # Hazard Ratio Table Reactive
    ForPval <- OptimalCutPTab_react()
    
    ## Survival Type
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## Determine Feature and sub feature
    if (Feature != "All_Features") {
      SubFeature <- input$subFeatureSelection
      Feature <- paste("<b>",Feature,"</b> - <b>",SubFeature,"</b></li>",sep = "")
      line2 <- paste("<li>The dataset is filtered by ",Feature,sep = "")
    }
    if (Feature == "All_Features") {
      line2 <- NULL
    }
    ## determine Sample Type
    if (is.null(metacol_sampletype) == T) {
      SampleType <- ""
    }
    else if (is.null(metacol_sampletype) == F) {
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        if (input$SampleTypeSelection == "All_Sample_Types") {
          SampleType <- " of all sample types.</li>"
        }
        else if (input$SampleTypeSelection != "All_Sample_Types") {
          SampleType <- paste(" of <b>",input$SampleTypeSelection,"</b> Patients.</li>",sep = "")
        }
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleType <- ".</li>"
      }
    }
    ## Determine Score method
    if (geneset_name %in% decon_score_cols) {
      scoreMethod <- " Pre-Processed Score"
    }
    if (input$GeneSetTabs == 2) {
      scoreMethod <- " Gene Expression"
      geneset_name <- paste(geneset_name,sep = "")
    }
    #pval <- get_lik_pval(ForPval)
    pval2 <- get_lik_pval(ForPval)
    HR <- HR_Tab[3,2]
    chacteristic <- HR_Tab[3,1]
    #if (pval == ">0.9") {
    #  pval2 <- 0.9
    #}
    #if (pval == "<0.001") {
    #  pval2 <- 0.001
    #}
    if (as.numeric(pval2) < 0.05) {
      pval_char <- "strongly associated"
    }
    if (as.numeric(pval2) >= 0.05 & as.numeric(pval2) < 0.1) {
      pval_char <- "moderately associated"
    }
    if (as.numeric(pval2) >= 0.1) {
      pval_char <- "not associated"
    }
    if (as.numeric(HR) > 1) {
      HR_char <- "high risk"
    }
    if (as.numeric(HR) <= 1) {
      HR_char <- "low risk"
    }
    
    line1 <- paste("<li><b>",SurvDateType,"</b> survival analysis ", SampleType,sep = "")
    line3 <- paste("<li>Kaplan-Meier survival curve dichotomized by <b>",geneset_name,"</b> <b>",scoreMethod,"</b> optimal cut-point.</li>",sep = "")
    line4 <- paste("<li>Cox hazard regression analysis finds a Likelihood Ratio P.value of <b>",pval2,"</b> and a Hazard Ratio of <b>",HR,"</b>, <b>",chacteristic,"</b> <b>",geneset_name,
                   "</b> is <b>",pval_char,"</b> with <b>",HR_char,"</b> for <b>",SurvDateType,"</b>.</li>",sep = "")
    
    if (is.null(line2)) {
      HTML(paste("<ul>",line1,line3,line4,"</ul>", sep = ""))
    }
    else if (!is.null(line2)) {
      HTML(paste("<ul>",line1,line2,line3,line4,"</ul>", sep = ""))
    }
    
  })
  
  output$OptimalCutPSumm <- renderPrint({
    
    tab <- OptimalCutPTab_react()
    out <- capture.output(summary(tab))
    
    con_line <- grep("^Concordance=",out,value = T)
    lik_line <- grep("^Likelihood ratio test=",out,value = T)
    wal_line <- grep("^Wald test",out,value = T)
    sco_line <- grep("^Score ",out,value = T)
    
    text <- paste("CoxH Summary:",con_line,lik_line,wal_line,sco_line,sep = "\n")
    cat(text)
    
  })
  
  output$ScutPointPlot <- renderPlot({
    
    plot <- ScutPointPlot_react()
    plot
    
  })
  
  output$CutPointHRtab <- renderTable({
    
    tab <- CutPointHRtab_react()
    tab
    
  })
  
  output$ssgseaCutPDensity <- renderPlot({
    
    plot <- ssgseaCutPDensity_react()
    plot
    
  })
  
  ####----Quantile Cut Point----####
  
  TopBottomCutP_react <- reactive({
    
    ## Assign variables
    surv_time_col <- input$SurvivalType_time
    surv_id_col <- input$SurvivalType_id
    meta_ssgsea <- ssGSEAmeta()
    
    ## Subset columns needed for plot
    meta_ssgsea_sdf <- meta_ssgsea[,c("SampleName",surv_time_col,surv_id_col,"TopBottomCutP")]
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "ID"
    
    ## Remove between cutoff samples
    meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf$TopBottomCutP != "BetweenCutoff"),]
    
    meta_ssgsea_sdf
    
  })
  
  TopBottomCutPTab_react <- reactive({
    
    meta_ssgsea_sdf <- TopBottomCutP_react()
    surv_time_col <- input$SurvivalType_time
    surv_id_col <- input$SurvivalType_id
    geneset <- gs_react()
    geneset_name <- names(geneset)
    
    meta_ssgsea_sdf[,"TopBottomCutP"] <- as.factor(meta_ssgsea_sdf[,"TopBottomCutP"])
    meta_ssgsea_sdf[,"TopBottomCutP"] <- relevel(meta_ssgsea_sdf[,"TopBottomCutP"], ref = "Low")
    
    ## Survival Function
    tab <- coxph(Surv(time,ID) ~ TopBottomCutP, data = meta_ssgsea_sdf)
    tab
    
  })
  
  SQuantileHRtab_react <- reactive({
    
    tab <- TopBottomCutPTab_react()
    tab <- tab %>% 
      gtsummary::tbl_regression(exp = TRUE) %>%
      as_gt()
    
    tab_df <- as.data.frame(tab)
    
    tab_df <- tab_df %>%
      select(label,estimate,ci,p.value)
    colnames(tab_df) <- c("Characteristic","Hazard Ratio","95% Confidence Interval","P.Value")
    
    tab_df
    
  })
  
  SquantPlot_react <- reactive({
    
    ## Assign variables
    meta_ssgsea_sdf <- TopBottomCutP_react()
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    scoreMethod <- input$ScoreMethod
    show_pval <- input$ShowPval
    ShowConfInt <- input$ShowConfInt
    xaxlim <- input$SurvXaxis * 365.25
    surv_time_col <- input$SurvivalType_time
    showLegend <- input$SurvLegendPos
    showMedSurv <- input$ShowMedSurvLine
    if (showMedSurv == T) {
      showMedSurv <- "hv"
    }
    else if (showMedSurv == F) {
      showMedSurv <- "none"
    }
    
    ## Survival Function
    fit <- survfit(Surv(time,ID) ~ TopBottomCutP, data = meta_ssgsea_sdf, type="kaplan-meier")
    
    ## Determine type of survival data - OS/EFS/PFS?
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## determine Feature and Sample Type label
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (SampleType == "All_Sample_Types") {
        if (Feature == "All_Features") {
          SampleTypeLab <- "All Features in All Patients\n"
        }
        if (Feature != "All_Features") {
          SampleTypeLab <- paste(Feature," in All Patients\n")
        }
      }
      else {
        if (Feature == "All_Features") {
          SampleTypeLab <- paste("All Features (",SampleType,") Patients\n",sep = "")
        }
        if (Feature != "All_Features") {
          SampleTypeLab <- paste(Feature," (",SampleType,") Patients\n",sep = "")
        }
      }
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (Feature == "All_Features") {
        SampleTypeLab <- "All Features in All Patients\n"
      }
      if (Feature != "All_Features") {
        SampleTypeLab <- paste(Feature," in All Patients\n")
      }
    }
    
    
    ## Determine Scoring method label
    if (input$GeneSetTabs == 2) {
      scoreMethodLab <- "Gene Expression"
    }
    else if (input$GeneSetTabs == 1) {
      if (input$GeneSetCat_Select == "Pre-Processed Scores") {
        scoreMethodLab <- "Pre-Processed score"
      }
      else {
        scoreMethodLab <- paste(scoreMethod, " score", sep = "")
      }
    }
    else if (input$GeneSetTabs == 3) {
      scoreMethodLab <- paste(scoreMethod, " score", sep = "")
    }
    
    ## Determind Plot title
    if (is.null(input$SurvPlotTitleQuantile)) {
      SurvPlotTitle <- paste("Survival curves of ",SampleTypeLab,
                             geneset_name," (",scoreMethodLab,")", sep = "")
    }
    else if (!is.null(input$SurvPlotTitleQuantile)) {
      if (input$SurvPlotTitleQuantile == "") {
        SurvPlotTitle <- paste("Survival curves of ",SampleTypeLab,
                               geneset_name," (",scoreMethodLab,")", sep = "")
      }
      else if (input$SurvPlotTitleQuantile != "") {
        SurvPlotTitle <- input$SurvPlotTitleQuantile
      }
    }
    
    
    ## Generate plot
    ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                         title = SurvPlotTitle,
                         xscale = c("d_y"),
                         break.time.by=365.25,
                         xlab = "Years", 
                         ylab = paste(SurvDateType,"Survival Probability"),
                         submain = "Based on Kaplan-Meier estimates",
                         caption = "created with survminer",
                         pval=show_pval,
                         conf.int = ShowConfInt,
                         ggtheme = theme_bw(),
                         font.title = c(16, "bold"),
                         font.submain = c(12, "italic"),
                         font.caption = c(12, "plain"),
                         font.x = c(14, "plain"),
                         font.y = c(14, "plain"),
                         font.tickslab = c(12, "plain"),
                         legend = showLegend,
                         risk.table.height = 0.20,
                         surv.median.line = showMedSurv
    )
    if (showMedSurv != "none") {
      MedSurvItem <- ggsurv[["plot"]][["layers"]][length(ggsurv[["plot"]][["layers"]])]
      MedSurvItem_df <- MedSurvItem[[1]][["data"]]
      MedSurvItem_df <- MedSurvItem_df[order(MedSurvItem_df[,1]),]
      MedSurvItem_df <- MedSurvItem_df %>%
        mutate(label = paste(round(MedSurvItem_df[,1]),"Days"))
      rownames(MedSurvItem_df) <- 1:nrow(MedSurvItem_df)
      if (nrow(MedSurvItem_df) > 1) {
        ggsurv$plot <- ggsurv$plot +
          geom_label_repel(data = MedSurvItem_df, aes(x = x1, y = y1, label = label, size = 4), label.size = NA, show.legend = FALSE)
      }
    }
    if (!is.null(input$SurvXaxis)) {
      ggsurv$plot$coordinates$limits$x <- c(0,xaxlim)
      ggsurv$table$coordinates$limits$x <- c(0,xaxlim)
    }
    
    ggsurv$table <- ggsurv$table + theme_cleantable()
    ggsurv
    
    
  })
  
  ssgseaQuantDensity_react <- reactive({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    ssgsea_meta <- ssGSEAmeta()
    scoreMethod <- input$ScoreMethod
    quantCutoff <- input$QuantPercent/100
    CutPlabel1 <- "(Bottom Cut-Point)"
    CutPlabel2 <- "(Top Cut-Point)"
    
    cols_selec <- c("SampleName",geneset_name)
    ssgsea_scores <- ssgsea_meta[,cols_selec]
    
    cutp_high <- round(quantile(ssgsea_scores[,geneset_name],1-quantCutoff,na.rm = T),3)
    cutp_low <- round(quantile(ssgsea_scores[,geneset_name],quantCutoff,na.rm = T),3)
    
    xints <- c(cutp_low,cutp_high)
    
    
    ## get score method for x and y labels
    if (input$GeneSetTabs == 2) {
      scoreMethodLab <- "Gene Expression Density"
      scoreMethodLab_x <- "Gene Expression (Log(exp+1))"
      ssgsea_scores[,geneset_name] <- log(ssgsea_scores[,geneset_name] + 1)
      cutp_high <- round(log(cutp_high + 1),3)
      cutp_low <- round(log(cutp_low + 1),3)
      xints <- c(cutp_low,cutp_high)
      
    }
    else if (input$GeneSetTabs != 2) {
      if (geneset_name %in% decon_score_cols) {
        scoreMethodLab <- "Pre-Processed Score Density"
        scoreMethodLab_x <- "Pre-Processed Score"
      }
      else {
        scoreMethodLab <- paste(scoreMethod, " Score Density", sep = "")
        scoreMethodLab_x <- paste(scoreMethod, " Score", sep = "")
      }
    }
    ## generate title based on input
    if (geneset_name %in% decon_score_cols) {
      dens_title <- paste(colnames(ssgsea_scores)[2]," Pre-Processed Score Density",sep = "")
    }
    else {
      if (input$GeneSetTabs == 2) {
        dens_title <- paste(colnames(ssgsea_scores)[2],"Log(exp+1)",scoreMethodLab)
      }
      else {
        dens_title <- paste(colnames(ssgsea_scores)[2],scoreMethodLab)
      }
    }
    
    dens_data <- density(ssgsea_scores[,geneset_name],na.rm = T)
    y_max <- max(dens_data$y)
    y_max_int <- y_max/6
    
    p <- ggplot(ssgsea_scores, aes(x=ssgsea_scores[,geneset_name])) + 
      geom_density(color="darkblue", fill="lightblue", alpha = 0.4) +
      xlab(scoreMethodLab_x) +
      ylab(scoreMethodLab) +
      ggtitle(dens_title) +
      theme(axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            plot.title = element_text(size = 20))
    p <- p + geom_vline(aes(xintercept = xints[1]), linetype = "dashed", color = "darkblue", size = 1)
    p <- p + geom_vline(aes(xintercept = xints[2]), linetype = "dashed", color = "darkblue", size = 1)
    p <- p + geom_text(aes(xints[1],y_max-(y_max_int/6),label = paste(as.character(xints[1]),CutPlabel1),hjust = -0.1,vjust = -0.1),size = 6, check_overlap = T)
    p <- p + geom_text(aes(xints[2],y_max-y_max_int,label = paste(as.character(xints[2]),CutPlabel2),hjust = -0.1,vjust = -0.1),size = 6, check_overlap = T)
    
    p
    
  })
  
  output$QuantSurvDescrip <- renderUI({
    
    ## Assign variables
    geneset <- gs_react()                      # Geneset Object
    geneset_name <- names(geneset)             # Geneset Name     
    scoreMethod <- input$ScoreMethod           # Scoring Method
    Feature <- input$FeatureSelection          # Feature Selected
    surv_time_col <- input$SurvivalType_time   # Survival Time Label
    HR_Tab <- SQuantileHRtab_react()             # Hazard Ratio Table Reactive
    ForPval <- TopBottomCutPTab_react()
    
    ## Survival Type
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## Determine Feature and sub feature
    if (Feature != "All_Features") {
      SubFeature <- input$subFeatureSelection
      Feature <- paste("<b>",Feature,"</b> - <b>",SubFeature,"</b></li>",sep = "")
      line2 <- paste("<li>The dataset is filtered by ",Feature,sep = "")
    }
    if (Feature == "All_Features") {
      line2 <- NULL
    }
    ## determine Sample Type
    if (is.null(metacol_sampletype) == T) {
      SampleType <- ""
    }
    else if (is.null(metacol_sampletype) == F) {
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        if (input$SampleTypeSelection == "All_Sample_Types") {
          SampleType <- " of all sample types.</li>"
        }
        else if (input$SampleTypeSelection != "All_Sample_Types") {
          SampleType <- paste(" of <b>",input$SampleTypeSelection,"</b> Patients.</li>",sep = "")
        }
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleType <- ".</li>"
      }
    }
    ## Determine Score method
    if (geneset_name %in% decon_score_cols) {
      scoreMethod <- " Pre-Processed Score"
    }
    if (input$GeneSetTabs == 2) {
      scoreMethod <- " Gene Expression"
      geneset_name <- paste(geneset_name,sep = "")
    }
    pval2 <- get_lik_pval(ForPval)
    HR <- HR_Tab[3,2]
    chacteristic <- HR_Tab[3,1]
    if (as.numeric(pval2) < 0.05) {
      pval_char <- "strongly associated"
    }
    if (as.numeric(pval2) >= 0.05 & as.numeric(pval2) < 0.1) {
      pval_char <- "moderately associated"
    }
    if (as.numeric(pval2) >= 0.1) {
      pval_char <- "not associated"
    }
    if (as.numeric(HR) > 1) {
      HR_char <- "high risk"
    }
    if (as.numeric(HR) <= 1) {
      HR_char <- "low risk"
    }
    
    line1 <- paste("<li><b>",SurvDateType,"</b> survival analysis ", SampleType,sep = "")
    line3 <- paste("<li>Kaplan-Meier survival curve dichotomized by <b>",geneset_name,"</b> <b>",scoreMethod,"</b> top/bottom cut-point.</li>",sep = "")
    line4 <- paste("<li>Cox hazard regression analysis finds a Likelihood Ratio P.value of <b>",pval2,"</b> and a Hazard Ratio of <b>",HR,"</b>, <b>",chacteristic,"</b> <b>",geneset_name,
                   "</b> is <b>",pval_char,"</b> with <b>",HR_char,"</b> for <b>",SurvDateType,"</b>.</li>",sep = "")
    
    if (is.null(line2)) {
      HTML(paste("<ul>",line1,line3,line4,"</ul>", sep = ""))
    }
    else if (!is.null(line2)) {
      HTML(paste("<ul>",line1,line2,line3,line4,"</ul>", sep = ""))
    }
    
  })
  
  output$QuantileCutPSumm <- renderPrint({
    
    tab <- TopBottomCutPTab_react()
    out <- capture.output(summary(tab))
    
    con_line <- grep("^Concordance=",out,value = T)
    lik_line <- grep("^Likelihood ratio test=",out,value = T)
    wal_line <- grep("^Wald test",out,value = T)
    sco_line <- grep("^Score ",out,value = T)
    
    text <- paste("CoxH Summary:",con_line,lik_line,wal_line,sco_line,sep = "\n")
    cat(text)
    
  })
  
  output$SquantPlot <- renderPlot({
    
    plot <- SquantPlot_react()
    plot
    
  })
  
  output$SQuantileHRtab <- renderTable({
    
    tab <- SQuantileHRtab_react()
    tab
    
  })
  
  output$ssgseaQuantDensity <- renderPlot({
    
    plot <- ssgseaQuantDensity_react()
    plot
    
  })
  
  ####----User Cut Point----####
  
  UserCutP_react <- reactive({
    
    ## Assign variables
    surv_time_col <- input$SurvivalType_time
    surv_id_col <- input$SurvivalType_id
    meta_ssgsea <- ssGSEAmeta()
    
    ## Subset columns needed for plot
    meta_ssgsea_sdf <- meta_ssgsea[,c("SampleName",surv_time_col,surv_id_col,"UserCutP")]
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "ID"
    
    meta_ssgsea_sdf
    
  })
  
  UserCutPTab_react <- reactive({
    
    meta_ssgsea_sdf <- UserCutP_react()
    surv_time_col <- input$SurvivalType_time
    surv_id_col <- input$SurvivalType_id
    geneset <- gs_react()
    geneset_name <- names(geneset)
    
    meta_ssgsea_sdf[,"UserCutP"] <- as.factor(meta_ssgsea_sdf[,"UserCutP"])
    meta_ssgsea_sdf[,"UserCutP"] <- relevel(meta_ssgsea_sdf[,"UserCutP"], ref = "Low")
    
    ## Survival Function
    tab <- coxph(Surv(time,ID) ~ UserCutP, data = meta_ssgsea_sdf)
    tab
    
  })
  
  SQuantileHR2tab_react <- reactive({
    
    tab <- UserCutPTab_react()
    tab <- tab %>% 
      gtsummary::tbl_regression(exp = TRUE) %>%
      as_gt()
    
    tab_df <- as.data.frame(tab)
    
    tab_df <- tab_df %>%
      select(label,estimate,ci,p.value)
    colnames(tab_df) <- c("Characteristic","Hazard Ratio","95% Confidence Interval","P.Value")
    
    tab_df
    
  })
  
  SquantPlot2_react <- reactive({
    
    ## Assign variables
    meta_ssgsea_sdf <- UserCutP_react()
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    scoreMethod <- input$ScoreMethod
    show_pval <- input$ShowPval
    ShowConfInt <- input$ShowConfInt
    xaxlim <- input$SurvXaxis * 365.25
    surv_time_col <- input$SurvivalType_time
    showLegend <- input$SurvLegendPos
    showMedSurv <- input$ShowMedSurvLine
    if (showMedSurv == T) {
      showMedSurv <- "hv"
    }
    else if (showMedSurv == F) {
      showMedSurv <- "none"
    }
    
    ## Survival Function
    fit <- survfit(Surv(time,ID) ~ UserCutP, data = meta_ssgsea_sdf, type="kaplan-meier")
    
    ## Determine type of survival data - OS/EFS/PFS?
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## determine Feature and Sample Type label
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (SampleType == "All_Sample_Types") {
        if (Feature == "All_Features") {
          SampleTypeLab <- "All Features in All Patients\n"
        }
        if (Feature != "All_Features") {
          SampleTypeLab <- paste(Feature," in All Patients\n")
        }
      }
      else {
        if (Feature == "All_Features") {
          SampleTypeLab <- paste("All Features (",SampleType,") Patients\n",sep = "")
        }
        if (Feature != "All_Features") {
          SampleTypeLab <- paste(Feature," (",SampleType,") Patients\n",sep = "")
        }
      }
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (Feature == "All_Features") {
        SampleTypeLab <- "All Features in All Patients\n"
      }
      if (Feature != "All_Features") {
        SampleTypeLab <- paste(Feature," in All Patients\n")
      }
    }
    
    
    ## Determine Scoring method label
    if (input$GeneSetTabs == 2) {
      scoreMethodLab <- "Gene Expression"
    }
    else if (input$GeneSetTabs == 1) {
      if (input$GeneSetCat_Select == "Pre-Processed Scores") {
        scoreMethodLab <- "Pre-Processed score"
      }
      else {
        scoreMethodLab <- paste(scoreMethod, " score", sep = "")
      }
    }
    else if (input$GeneSetTabs == 3) {
      scoreMethodLab <- paste(scoreMethod, " score", sep = "")
    }
    
    ## Determind Plot title
    if (is.null(input$SurvPlotTitleUser)) {
      SurvPlotTitle <- paste("Survival curves of ",SampleTypeLab,
                             geneset_name," (",scoreMethodLab,")", sep = "")
    }
    else if (!is.null(input$SurvPlotTitleUser)) {
      if (input$SurvPlotTitleUser == "") {
        SurvPlotTitle <- paste("Survival curves of ",SampleTypeLab,
                               geneset_name," (",scoreMethodLab,")", sep = "")
      }
      else if (input$SurvPlotTitleUser != "") {
        SurvPlotTitle <- input$SurvPlotTitleUser
      }
    }
    
    ## Generate plot
    ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                         title = SurvPlotTitle,
                         xscale = c("d_y"),
                         break.time.by=365.25,
                         xlab = "Years", 
                         ylab = paste(SurvDateType,"Survival Probability"),
                         submain = "Based on Kaplan-Meier estimates",
                         caption = "created with survminer",
                         pval=show_pval,
                         conf.int = ShowConfInt,
                         ggtheme = theme_bw(),
                         font.title = c(16, "bold"),
                         font.submain = c(12, "italic"),
                         font.caption = c(12, "plain"),
                         font.x = c(14, "plain"),
                         font.y = c(14, "plain"),
                         font.tickslab = c(12, "plain"),
                         legend = showLegend,
                         risk.table.height = 0.20,
                         surv.median.line = showMedSurv
    )
    if (showMedSurv != "none") {
      MedSurvItem <- ggsurv[["plot"]][["layers"]][length(ggsurv[["plot"]][["layers"]])]
      MedSurvItem_df <- MedSurvItem[[1]][["data"]]
      MedSurvItem_df <- MedSurvItem_df[order(MedSurvItem_df[,1]),]
      MedSurvItem_df <- MedSurvItem_df %>%
        mutate(label = paste(round(MedSurvItem_df[,1]),"Days"))
      rownames(MedSurvItem_df) <- 1:nrow(MedSurvItem_df)
      if (nrow(MedSurvItem_df) > 1) {
        ggsurv$plot <- ggsurv$plot +
          geom_label_repel(data = MedSurvItem_df, aes(x = x1, y = y1, label = label, size = 4), label.size = NA, show.legend = FALSE)
      }
    }
    if (!is.null(input$SurvXaxis)) {
      ggsurv$plot$coordinates$limits$x <- c(0,xaxlim)
      ggsurv$table$coordinates$limits$x <- c(0,xaxlim)
    }
    
    ggsurv$table <- ggsurv$table + theme_cleantable()
    ggsurv
    
    
  })
  
  ssgseaQuant2Density_react <- reactive({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    ssgsea_meta <- ssGSEAmeta()
    scoreMethod <- input$ScoreMethod
    quantCutoff <- input$QuantPercent2/100
    CutPlabel1 <- "(User Cut-Point)"
    
    cols_selec <- c("SampleName",geneset_name)
    ssgsea_scores <- ssgsea_meta[,cols_selec]
    
    cutp_user <- round(quantile(ssgsea_scores[,geneset_name],quantCutoff,na.rm = T),3)
    
    ## get score method for x and y labels
    if (input$GeneSetTabs == 2) {
      scoreMethodLab <- "Gene Expression Density"
      scoreMethodLab_x <- "Gene Expression (Log(exp+1))"
      ssgsea_scores[,geneset_name] <- log(ssgsea_scores[,geneset_name] + 1)
      cutp_user <- round(log(cutp_user + 1),3)
      
    }
    else if (input$GeneSetTabs != 2) {
      if (geneset_name %in% decon_score_cols) {
        scoreMethodLab <- "Pre-Processed Score Density"
        scoreMethodLab_x <- "Pre-Processed Score"
      }
      else {
        scoreMethodLab <- paste(scoreMethod, " Score Density", sep = "")
        scoreMethodLab_x <- paste(scoreMethod, " Score", sep = "")
      }
    }
    ## generate title based on input
    if (geneset_name %in% decon_score_cols) {
      dens_title <- paste(colnames(ssgsea_scores)[2]," Pre-Processed Score Density",sep = "")
    }
    else {
      if (input$GeneSetTabs == 2) {
        dens_title <- paste(colnames(ssgsea_scores)[2],"Log(exp+1)",scoreMethodLab)
      }
      else {
        dens_title <- paste(colnames(ssgsea_scores)[2],scoreMethodLab)
      }
    }
    
    dens_data <- density(ssgsea_scores[,geneset_name],na.rm = T)
    y_max <- max(dens_data$y)
    y_max_int <- y_max/6
    
    p <- ggplot(ssgsea_scores, aes(x=ssgsea_scores[,geneset_name])) + 
      geom_density(color="darkblue", fill="lightblue", alpha = 0.4) +
      xlab(scoreMethodLab_x) +
      ylab(scoreMethodLab) +
      ggtitle(dens_title) +
      theme(axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            plot.title = element_text(size = 20))
    p <- p + geom_vline(aes(xintercept = cutp_user), linetype = "dashed", color = "darkblue", size = 1)
    p <- p + geom_text(aes(cutp_user,y_max-(y_max_int/6),label = paste(as.character(cutp_user),CutPlabel1),hjust = -0.1,vjust = -0.1),size = 6, check_overlap = T)
    p
    
  })
  
  output$Quant2SurvDescrip <- renderUI({
    
    ## Assign variables
    geneset <- gs_react()                      # Geneset Object
    geneset_name <- names(geneset)             # Geneset Name     
    scoreMethod <- input$ScoreMethod           # Scoring Method
    Feature <- input$FeatureSelection          # Feature Selected
    surv_time_col <- input$SurvivalType_time   # Survival Time Label
    HR_Tab <- SQuantileHR2tab_react()             # Hazard Ratio Table Reactive
    ForPval <- UserCutPTab_react()
    
    ## Survival Type
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## Determine Feature and sub feature
    if (Feature != "All_Features") {
      SubFeature <- input$subFeatureSelection
      Feature <- paste("<b>",Feature,"</b> - <b>",SubFeature,"</b></li>",sep = "")
      line2 <- paste("<li>The dataset is filtered by ",Feature,sep = "")
    }
    if (Feature == "All_Features") {
      line2 <- NULL
    }
    ## determine Sample Type
    if (is.null(metacol_sampletype) == T) {
      SampleType <- ""
    }
    else if (is.null(metacol_sampletype) == F) {
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        if (input$SampleTypeSelection == "All_Sample_Types") {
          SampleType <- " of all sample types.</li>"
        }
        else if (input$SampleTypeSelection != "All_Sample_Types") {
          SampleType <- paste(" of <b>",input$SampleTypeSelection,"</b> Patients.</li>",sep = "")
        }
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleType <- ".</li>"
      }
    }
    ## Determine Score method
    if (geneset_name %in% decon_score_cols) {
      scoreMethod <- " Pre-Processed Score"
    }
    if (input$GeneSetTabs == 2) {
      scoreMethod <- " Gene Expression"
      geneset_name <- paste(geneset_name,sep = "")
    }
    pval2 <- get_lik_pval(ForPval)
    HR <- HR_Tab[3,2]
    chacteristic <- HR_Tab[3,1]
    if (as.numeric(pval2) < 0.05) {
      pval_char <- "strongly associated"
    }
    if (as.numeric(pval2) >= 0.05 & as.numeric(pval2) < 0.1) {
      pval_char <- "moderately associated"
    }
    if (as.numeric(pval2) >= 0.1) {
      pval_char <- "not associated"
    }
    if (as.numeric(HR) > 1) {
      HR_char <- "high risk"
    }
    if (as.numeric(HR) <= 1) {
      HR_char <- "low risk"
    }
    
    line1 <- paste("<li><b>",SurvDateType,"</b> survival analysis ", SampleType,sep = "")
    line3 <- paste("<li>Kaplan-Meier survival curve dichotomized by <b>",geneset_name,"</b> <b>",scoreMethod,"</b> above/below user cut-point.</li>",sep = "")
    line4 <- paste("<li>Cox hazard regression analysis finds a Likelihood Ratio P.value of <b>",pval2,"</b> and a Hazard Ratio of <b>",HR,"</b>, <b>",chacteristic,"</b> <b>",geneset_name,
                   "</b> is <b>",pval_char,"</b> with <b>",HR_char,"</b> for <b>",SurvDateType,"</b>.</li>",sep = "")
    
    if (is.null(line2)) {
      HTML(paste("<ul>",line1,line3,line4,"</ul>", sep = ""))
    }
    else if (!is.null(line2)) {
      HTML(paste("<ul>",line1,line2,line3,line4,"</ul>", sep = ""))
    }
    
  })
  
  output$UserCutPSumm <- renderPrint({
    
    tab <- UserCutPTab_react()
    out <- capture.output(summary(tab))
    
    con_line <- grep("^Concordance=",out,value = T)
    lik_line <- grep("^Likelihood ratio test=",out,value = T)
    wal_line <- grep("^Wald test",out,value = T)
    sco_line <- grep("^Score ",out,value = T)
    
    text <- paste("CoxH Summary:",con_line,lik_line,wal_line,sco_line,sep = "\n")
    cat(text)
    
  })
  
  output$SquantPlot2 <- renderPlot({
    
    plot <- SquantPlot2_react()
    plot
    
  })
  
  output$SQuantileHR2tab <- renderTable({
    
    tab <- SQuantileHR2tab_react()
    tab
    
  })
  
  output$ssgseaQuant2Density <- renderPlot({
    
    plot <- ssgseaQuant2Density_react()
    plot
    
  })
  
  ####----Univariate Feature Survival----####
  
  UniVarFeat_react <- reactive({
    
    if (length(input$SingleSurvivalFeature > 0)) {
      
      ## Assign variables
      Feature <- input$SingleSurvivalFeature
      surv_time_col <- input$SurvivalType_time
      surv_id_col <- input$SurvivalType_id
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      quantCutoff2 <- input$QuantPercent2/100 #Quantile cutoff given by user
      geneset_name <- names(gs_react())
      
      meta_ssgsea <- ssGSEAmeta()
      
      if (input$UniVarNAcheck == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[which(meta_ssgsea[,Feature] != "Inf"),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature],ignore.case = T, invert = T),]
        ### Re-Perform Stat functions
        #meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        #meta_ssgsea$QuartileCutP <- paste("", meta_ssgsea$VAR_Q, sep="")
        #meta_ssgsea$MedianCutP <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        #meta_ssgsea$TopBottomCutP <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        #meta_ssgsea$UserCutP <- quantile_conversion2(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff2)
        
      }
      
      if (input$UniVarContCheck == TRUE) {
        if (input$UniVarContHiLoCheck == TRUE) {
          meta_ssgsea[,Feature] <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == Feature)])
        }
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      colnames(meta_ssgsea_sdf)[4] <- gsubCheck(colnames(meta_ssgsea_sdf)[4])
      Feature <- gsubCheck(Feature)
      #colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
      #Feature <- gsub("[[:punct:]]","_",Feature)
      colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "time"
      colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "ID"
      
      ## Remove between cutoff samples
      if (Feature == "TopBottomCutP") {
        meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf$TopBottomCutP != "BetweenCutoff"),]
      }
      
      meta_ssgsea_sdf
      
    }
    
  })
  
  UniVarFeatTab_react <- reactive({
    
    meta_ssgsea_sdf <- UniVarFeat_react()
    surv_time_col <- input$SurvivalType_time
    surv_id_col <- input$SurvivalType_id
    Feature <- input$SingleSurvivalFeature
    ref_Feature <- input$SurvFeatVariableUni
    Feature <- gsubCheck(Feature)
    
    if (input$UniVarContCheck == FALSE) {
      meta_ssgsea_sdf[,Feature] <- as.factor(meta_ssgsea_sdf[,Feature])
      meta_ssgsea_sdf[,Feature] <- relevel(meta_ssgsea_sdf[,Feature], ref = ref_Feature)
    }
    if (input$UniVarContCheck == TRUE) {
      if (input$UniVarContHiLoCheck == TRUE) {
        meta_ssgsea_sdf[,Feature] <- as.factor(meta_ssgsea_sdf[,Feature])
        meta_ssgsea_sdf[,Feature] <- relevel(meta_ssgsea_sdf[,Feature], ref = ref_Feature)
      }
    }
    
    ## Survival Function
    tab <- coxph(as.formula(paste("Surv(time,ID) ~ ",Feature,sep = "")),
                 data = meta_ssgsea_sdf)
    tab
    
  })
  
  SSingleFeatureHRtab_react <- reactive({
    
    Feature <- input$SingleSurvivalFeature
    tab <- UniVarFeatTab_react()
    tab <- tab %>% 
      gtsummary::tbl_regression(exp = TRUE) %>%
      as_gt()
    
    tab_df <- as.data.frame(tab)
    
    tab_df[is.na(tab_df)] <- ""
    tab_df <- tab_df %>%
      select(label,n_obs,estimate,std.error,ci,p.value)
    tab_df[1,1] <- Feature
    colnames(tab_df) <- c("Variable","N","Hazard Ratio","Std. Error","95% Confidence Interval","P.Value")
    
    tab_df
    
  })
  
  output$SSingleFeatureHRtab <- renderTable({
    
    tab <- SSingleFeatureHRtab_react()
    tab
    
  })
  
  featSplot_react <- reactive({
    
    ## Assign variables
    meta_ssgsea_sdf <- UniVarFeat_react()
    SampleType <- input$SampleTypeSelection
    Feature_lab <- input$SingleSurvivalFeature
    show_pval <- input$ShowPval
    ShowConfInt <- input$ShowConfInt
    xaxlim <- input$SurvXaxis * 365.25
    surv_time_col <- input$SurvivalType_time
    showLegend <- input$SurvLegendPos
    showMedSurv <- input$ShowMedSurvLine
    if (showMedSurv == T) {
      showMedSurv <- "hv"
    }
    else if (showMedSurv == F) {
      showMedSurv <- "none"
    }
    
    
    Feature <- colnames(meta_ssgsea_sdf)[4]
    
    form <- paste("Surv(time,ID) ~ ",Feature,sep = "")
    form2 <- as.formula(form)
    fit <- eval(substitute(survfit(form2,data = meta_ssgsea_sdf, type="kaplan-meier")))
    
    ## Determine type of survival data - OS/EFS/PFS?
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## determine Feature and Sample Type label
    ## Adjust 'Sample Type' for label 
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      SampleTypeLab <- paste(" (",SampleType,") ",sep = "")
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      SampleTypeLab <- " "
    }
    
    ## Determind Plot title
    if (is.null(input$SurvPlotTitleUniVar)) {
      SurvPlotTitle <- paste("Survival curves of ",Feature_lab,SampleTypeLab,"Patients", sep = "")
    }
    else if (!is.null(input$SurvPlotTitleUniVar)) {
      if (input$SurvPlotTitleUniVar == "") {
        SurvPlotTitle <- paste("Survival curves of ",Feature_lab,SampleTypeLab,"Patients", sep = "")
      }
      else if (input$SurvPlotTitleUniVar != "") {
        SurvPlotTitle <- input$SurvPlotTitleUniVar
      }
    }
    
    
    ## Generate plot
    ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                         title = SurvPlotTitle,
                         xscale = c("d_y"),
                         break.time.by=365.25,
                         xlab = "Years", 
                         ylab = paste(SurvDateType,"Survival Probability"),
                         submain = "Based on Kaplan-Meier estimates",
                         caption = "created with survminer",
                         pval=show_pval,
                         conf.int = ShowConfInt,
                         ggtheme = theme_bw(),
                         font.title = c(16, "bold"),
                         font.submain = c(12, "italic"),
                         font.caption = c(12, "plain"),
                         font.x = c(14, "plain"),
                         font.y = c(14, "plain"),
                         font.tickslab = c(12, "plain"),
                         legend = showLegend,
                         risk.table.height = 0.20,
                         surv.median.line = showMedSurv
    )
    if (showMedSurv != "none") {
      MedSurvItem <- ggsurv[["plot"]][["layers"]][length(ggsurv[["plot"]][["layers"]])]
      MedSurvItem_df <- MedSurvItem[[1]][["data"]]
      MedSurvItem_df <- MedSurvItem_df[order(MedSurvItem_df[,1]),]
      MedSurvItem_df <- MedSurvItem_df %>%
        mutate(label = paste(round(MedSurvItem_df[,1]),"Days"))
      rownames(MedSurvItem_df) <- 1:nrow(MedSurvItem_df)
      if (nrow(MedSurvItem_df) > 1) {
        ggsurv$plot <- ggsurv$plot +
          geom_label_repel(data = MedSurvItem_df, aes(x = x1, y = y1, label = label, size = 4), label.size = NA, show.legend = FALSE)
      }
    }
    if (!is.null(input$SurvXaxis)) {
      ggsurv$plot$coordinates$limits$x <- c(0,xaxlim)
      ggsurv$table$coordinates$limits$x <- c(0,xaxlim)
    }
    
    ggsurv$table <- ggsurv$table + theme_cleantable()
    ggsurv
    
    
  })
  
  output$featSplot <- renderPlot({
    plot <- featSplot_react()
    plot
  })
  
  output$UnivarSummExpl <- renderUI({
    
    ## variables
    geneset <- gs_react()                      # Geneset Object
    geneset_name <- names(geneset)             # Geneset Name     
    scoreMethod <- input$ScoreMethod           # Scoring Method
    Feature <- input$FeatureSelection          # Feature Selected
    surv_time_col <- input$SurvivalType_time   # Survival Time Label
    FeatureSelec <- input$SingleSurvivalFeature
    contCheck <- input$UniVarContCheck
    hiloCheck <- input$UniVarContHiLoCheck
    tab <- UniVarFeatTab_react()
    hr_tab <- SSingleFeatureHRtab_react()
    
    ## Survival Type
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## determine Sample Type
    if (is.null(metacol_sampletype) == T) {
      SampleType <- ""
    }
    else if (is.null(metacol_sampletype) == F) {
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        if (input$SampleTypeSelection == "All_Sample_Types") {
          SampleType <- " of all sample types.</li>"
        }
        else if (input$SampleTypeSelection != "All_Sample_Types") {
          SampleType <- paste(" of <b>",input$SampleTypeSelection,"</b> Patients.</li>",sep = "")
        }
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleType <- ".</li>"
      }
    }
    
    ## Determine Feature and sub feature
    if (Feature != "All_Features") {
      SubFeature <- input$subFeatureSelection
      Feature <- paste("<b>",Feature,"</b> - <b>",SubFeature,"</b></li>",sep = "")
      line2 <- paste("<li>The dataset is filtered by ",Feature,sep = "")
    }
    else if (Feature == "All_Features") {
      line2 <- NULL
    }
    if (contCheck == TRUE) {
      if (hiloCheck == FALSE) {
        pval2 <- get_lik_pval(tab)
        HR <- hr_tab[nrow(hr_tab),3]
        referenced <- ""
      }
      if (hiloCheck == TRUE) {
        pval2 <- get_lik_pval(tab)
        HR <- hr_tab[nrow(hr_tab),3]
        referenced <- paste(" -",hr_tab[nrow(hr_tab),1])
      }
    }
    else if (contCheck == FALSE) {
      pval2 <- get_lik_pval(tab)
      HR <- hr_tab[nrow(hr_tab),3]
      referenced <- paste(" -",hr_tab[nrow(hr_tab),1])
    }
    
    
    if (as.numeric(pval2) < 0.05) {
      pval_char <- "strongly associated"
    }
    if (as.numeric(pval2) >= 0.05 & as.numeric(pval2) < 0.1) {
      pval_char <- "moderately associated"
    }
    if (as.numeric(pval2) >= 0.1) {
      pval_char <- "not associated"
    }
    if (as.numeric(HR) > 1) {
      HR_char <- "high risk"
    }
    if (as.numeric(HR) <= 1) {
      HR_char <- "low risk"
    }
    
    if (FeatureSelec %in% StatCols) {
      FeatureSelec <- paste(FeatureSelec," (",geneset_name,")",sep = "")
    }
    
    line1 <- paste("<li><b>",SurvDateType,"</b> survival analysis ", SampleType,sep = "")
    line3 <- paste("<li>Kaplan-Meier survival curve categorized by <b>",FeatureSelec,"</b>.</li>",sep = "")
    line4 <- paste("<li>Cox hazard regression analysis finds a Likelihood Ratio P.value of <b>",pval2,"</b> and a Hazard Ratio of <b>",HR,"</b>, <b>",FeatureSelec,referenced,
                   "</b> is <b>",pval_char,"</b> with <b>",HR_char,"</b> for <b>",SurvDateType,"</b>.</li>",sep = "")
    
    if (is.null(line2)) {
      HTML(paste("<ul>",line1,line3,line4,"</ul>", sep = ""))
    }
    else if (!is.null(line2)) {
      HTML(paste("<ul>",line1,line2,line3,line4,"</ul>", sep = ""))
    }
    
  })
  
  output$UnivarSummary <- renderPrint({
    
    tab <- UniVarFeatTab_react()
    out <- capture.output(summary(tab))
    
    con_line <- grep("^Concordance=",out,value = T)
    lik_line <- grep("^Likelihood ratio test=",out,value = T)
    wal_line <- grep("^Wald test",out,value = T)
    sco_line <- grep("^Score ",out,value = T)
    
    text <- paste("Coxh Summary:",con_line,lik_line,wal_line,sco_line,sep = "\n")
    cat(text)
    
  })
  
  SinglevarForestPlot_react <- reactive({
    
    if (length(input$SingleSurvivalFeature > 0)) {
      
      tab <- UniVarFeatTab_react()
      meta_ssgsea_sdf <- UniVarFeat_react()
      Feature <- input$SingleSurvivalFeature
      forextFont <- input$ForestFontSize
      forest <- ggforest(tab,
                         data = meta_ssgsea_sdf,
                         main = paste("Hazard Ratio Modeling: ",paste(Feature,collapse = ", "),sep = ""),
                         fontsize = forextFont)
      forest
    }
    
  })
  
  output$SinglevarForestPlot <- renderPlot({
    
    forest <- SinglevarForestPlot_react()
    forest
    
  })
  
  UnivarLinearityPlot_react <- reactive({
    
    if (length(input$SingleSurvivalFeature > 0)) {
      
      residType <- input$ResidualTypeUni
      AxisFont <- input$linAxisFont
      MainFont <- input$linMainFont
      TickFont <- input$linTickFont
      linpredict <- input$linPredict1
      tab <- UniVarFeatTab_react()
      Feature <- input$SingleSurvivalFeature
      
      p <- ggcoxdiagnostics(tab,
                            type = residType,
                            sline = T,
                            sline.se = T,
                            ggtheme = theme_minimal(),
                            ox.scale = linpredict)
      p <- ggpar(p,
                 font.x = AxisFont,
                 font.y = AxisFont,
                 font.main = MainFont,
                 font.tickslab = TickFont,
                 main = paste("Linearity Plot Featuring: ",Feature, sep = ""),
                 ylab = paste(str_to_title(residType)," Residuals", sep = "")
      )
      p
      
    }
    
  })
  
  
  output$UnivarLinearityPlot <- renderPlot({
    
    p <- UnivarLinearityPlot_react()
    p
    
  })
  
  ####----Bivariate Add----####
  
  BiVarAddFeature_react <- reactive({
    
    if (length(input$SurvivalFeatureBi1 > 0) & length(input$SurvivalFeatureBi2 > 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature1 <- input$SurvivalFeatureBi1
      Feature2 <- input$SurvivalFeatureBi2
      Feat1Var <- input$SurvFeatVariableBi1
      Feat2Var <- input$SurvFeatVariableBi2
      surv_time_col <- input$SurvivalType_time
      surv_id_col <- input$SurvivalType_id
      quantCutoff <- input$QuantPercent/100
      quantCutoff2 <- input$QuantPercent2/100
      meta_ssgsea <- ssGSEAmeta()
      
      ## Determine type of survival data - OS/EFS/PFS?
      SurvDateType <- sub("\\..*","",surv_time_col)
      
      ## Remove rows with NA in survival column
      meta_ssgsea <- meta_ssgsea[!is.na(meta_ssgsea[,surv_time_col]),]
      
      if (input$BiVarAddNAcheck1 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature1]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[which(meta_ssgsea[,Feature1] != "Inf"),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature1],ignore.case = T, invert = T),]
        ### Re-Perform Stat functions
        #meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        #meta_ssgsea$QuartileCutP <- paste("", meta_ssgsea$VAR_Q, sep="")
        #meta_ssgsea$MedianCutP <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        #meta_ssgsea$TopBottomCutP <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        #meta_ssgsea$UserCutP <- quantile_conversion2(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff2)
        
      }
      if (input$BiVarAddNAcheck2 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature2]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[which(meta_ssgsea[,Feature2] != "Inf"),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature2],ignore.case = T, invert = T),]
        ### Re-Perform Stat functions
        #meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        #meta_ssgsea$QuartileCutP <- paste("", meta_ssgsea$VAR_Q, sep="")
        #meta_ssgsea$MedianCutP <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        #meta_ssgsea$TopBottomCutP <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        #meta_ssgsea$UserCutP <- quantile_conversion2(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff2)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature1,Feature2)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      colnames(meta_ssgsea_sdf)[4] <- gsubCheck(colnames(meta_ssgsea_sdf)[4])
      if (!is.numeric(meta_ssgsea_sdf[,4])) {
        meta_ssgsea_sdf[,4] <- gsubCheck(meta_ssgsea_sdf[,4])
      }
      Feature1 <- gsubCheck(Feature1)
      Feat1Var <- gsubCheck(Feat1Var)
      colnames(meta_ssgsea_sdf)[5] <- gsubCheck(colnames(meta_ssgsea_sdf)[5])
      if (!is.numeric(meta_ssgsea_sdf[,5])) {
        meta_ssgsea_sdf[,5] <- gsubCheck(meta_ssgsea_sdf[,5])
      }
      Feature2 <- gsubCheck(Feature2)
      Feat2Var <- gsubCheck(Feat2Var)
      
      if (input$BiVarAddContCheck1 == FALSE) {
        meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
        meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      }
      else if (input$BiVarAddContCheck1 == TRUE) {
        if (input$BiVarAddContHiLoCheck1 == TRUE) {
          meta_ssgsea_sdf[,Feature1] <- highlow(as.numeric(meta_ssgsea_sdf[, which(colnames(meta_ssgsea_sdf) == Feature1)]))
          meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
          meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
        }
        else if (input$BiVarAddContHiLoCheck1 == FALSE) {
          meta_ssgsea_sdf[,Feature1] <- as.numeric(meta_ssgsea_sdf[,Feature1])
        }
      }
      if (input$BiVarAddContCheck2 == FALSE) {
        meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
        meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      }
      else if (input$BiVarAddContCheck2 == TRUE) {
        if (input$BiVarAddContHiLoCheck2 == TRUE) {
          meta_ssgsea_sdf[,Feature2] <- highlow(as.numeric(meta_ssgsea_sdf[, which(colnames(meta_ssgsea_sdf) == Feature2)]))
          meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
          meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
        }
        else if (input$BiVarAddContHiLoCheck2 == FALSE) {
          meta_ssgsea_sdf[,Feature2] <- as.numeric(meta_ssgsea_sdf[,Feature2])
        }
      }
    }
    
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "ID"
    
    meta_ssgsea_sdf
    
  })
  
  BiVarAddTab_react <- reactive({
    
    meta_ssgsea_sdf <- BiVarAddFeature_react()
    #Feature1 <- input$SurvivalFeatureBi1
    #Feature2 <- input$SurvivalFeatureBi2
    #Feature1 <- gsub("[[:punct:]]","_",Feature1)
    #Feature2 <- gsub("[[:punct:]]","_",Feature2)
    
    Feature1 <- colnames(meta_ssgsea_sdf)[4]
    Feature2 <- colnames(meta_ssgsea_sdf)[5]
    
    
    form <- paste("Surv(time,ID) ~ ",paste(Feature1,"+",Feature2,sep = ""),sep = "")
    form2 <- as.formula(form)
    tab <- eval(substitute(coxph(form2,data = meta_ssgsea_sdf)))
    
    tab
    
    
  })
  
  BiVarAddTabFeat1_react <- reactive({
    
    meta_ssgsea_sdf <- BiVarAddFeature_react()
    #Feature1 <- input$SurvivalFeatureBi1
    Feature1 <- colnames(meta_ssgsea_sdf)[4]
    #Feature2 <- colnames(meta_ssgsea_sdf)[5]
    
    form <- paste("Surv(time,ID) ~ ",Feature1,sep = "")
    form2 <- as.formula(form)
    tab <- eval(substitute(coxph(form2,data = meta_ssgsea_sdf)))
    
    tab
    
    
  })
  
  BiVarAddTabFeat2_react <- reactive({
    
    meta_ssgsea_sdf <- BiVarAddFeature_react()
    #Feature2 <- input$SurvivalFeatureBi2
    #Feature1 <- colnames(meta_ssgsea_sdf)[4]
    Feature2 <- colnames(meta_ssgsea_sdf)[5]
    
    form <- paste("Surv(time,ID) ~ ",Feature2,sep = "")
    form2 <- as.formula(form)
    tab <- eval(substitute(coxph(form2,data = meta_ssgsea_sdf)))
    
    tab
    
    
  })
  
  BiVarAddHRTab_react <- reactive({
    
    tab <- BiVarAddTab_react()
    tab <- tab %>% 
      gtsummary::tbl_regression(exp = TRUE) %>%
      as_gt()
    
    tab_df <- as.data.frame(tab)
    
    tab_df[is.na(tab_df)] <- ""
    tab_df <- tab_df %>%
      select(label,n_obs,estimate,std.error,ci,p.value)
    colnames(tab_df) <- c("Variable","N","Hazard Ratio","Std. Error","95% Confidence Interval","P.Value")
    
    tab_df
    
  })
  
  output$BiFeatureHRtab <- renderTable({
    
    tab <- BiVarAddHRTab_react()
    tab
    
  })
  
  output$bivarSummary <- renderPrint({
    
    tab <- BiVarAddTab_react()
    out <- capture.output(summary(tab))
    
    con_line <- grep("^Concordance=",out,value = T)
    lik_line <- grep("^Likelihood ratio test=",out,value = T)
    wal_line <- grep("^Wald test",out,value = T)
    sco_line <- grep("^Score ",out,value = T)
    
    text <- paste("Coxh Summary:",con_line,lik_line,wal_line,sco_line,sep = "\n")
    cat(text)
    
  })
  
  output$bivarAnova1 <- renderPrint({
    
    tab1 <- BiVarAddTab_react()
    tab2 <- BiVarAddTabFeat1_react()
    
    annova_res <- anova(tab1,tab2)
    
    out <- capture.output(annova_res)
    
    line1 <- out[3]
    line2 <- out[4]
    line3 <- out[5]
    line4 <- out[6]
    line5 <- out[7]
    
    text <- paste("Model Comparison:",line1,line2,line3,line4,line5,sep = "\n")
    cat(text)
    
  })
  
  output$bivarAnova2 <- renderPrint({
    
    tab1 <- BiVarAddTab_react()
    tab2 <- BiVarAddTabFeat2_react()
    
    annova_res <- anova(tab1,tab2)
    
    out <- capture.output(annova_res)
    
    line1 <- out[3]
    line2 <- out[4]
    line3 <- out[5]
    line4 <- out[6]
    line5 <- out[7]
    
    text <- paste("Model Comparison:",line1,line2,line3,line4,line5,sep = "\n")
    cat(text)
    
  })
  
  #output$BivarAddSummExpl <- renderUI({
  #  
  #  ## variables
  #  geneset <- gs_react()                      # Geneset Object
  #  geneset_name <- names(geneset)             # Geneset Name     
  #  scoreMethod <- input$ScoreMethod           # Scoring Method
  #  Feature <- input$FeatureSelection          # Feature Selected
  #  surv_time_col <- input$SurvivalType_time   # Survival Time Label
  #  Feature1 <- input$SurvivalFeatureBi1
  #  Feature2 <- input$SurvivalFeatureBi2
  #  Feat1Var <- input$SurvFeatVariableBi1
  #  Feat2Var <- input$SurvFeatVariableBi2
  #  contCheck1 <- input$BiVarAddContCheck1
  #  hiloCheck1 <- input$BiVarAddContHiLoCheck1
  #  contCheck2 <- input$BiVarAddContCheck2
  #  hiloCheck2 <- input$BiVarAddContHiLoCheck2
  #  tab <- UniVarFeatTab_react()
  #  hr_tab <- SSingleFeatureHRtab_react()
  #  
  #  ## Survival Type
  #  SurvDateType <- sub("\\..*","",surv_time_col)
  #  
  #  ## determine Sample Type
  #  if (is.null(metacol_sampletype) == T) {
  #    SampleType <- ""
  #  }
  #  else if (is.null(metacol_sampletype) == F) {
  #    if (length(unique(meta[,metacol_sampletype])) > 1) {
  #      if (input$SampleTypeSelection == "All_Sample_Types") {
  #        SampleType <- " of all sample types.</li>"
  #      }
  #      else if (input$SampleTypeSelection != "All_Sample_Types") {
  #        SampleType <- paste(" of <b>",input$SampleTypeSelection,"</b> Patients.</li>",sep = "")
  #      }
  #    }
  #    if (length(unique(meta[,metacol_sampletype])) <= 1) {
  #      SampleType <- ".</li>"
  #    }
  #  }
  #  
  #  ## Determine Feature and sub feature
  #  if (Feature != "All_Features") {
  #    SubFeature <- input$subFeatureSelection
  #    Feature <- paste("<b>",Feature,"</b> - <b>",SubFeature,"</b></li>",sep = "")
  #    line2 <- paste("<li>The dataset is filtered by ",Feature,sep = "")
  #  }
  #  if (Feature == "All_Features") {
  #    line2 <- NULL
  #  }
  #  
  #  if (contCheck == TRUE) {
  #    if (hiloCheck == FALSE) {
  #      pval2 <- get_lik_pval(tab)
  #      HR <- hr_tab[1,3]
  #      referenced <- ""
  #    }
  #    if (hiloCheck == TRUE) {
  #      pval2 <- get_lik_pval(tab)
  #      HR <- hr_tab[3,3]
  #      referenced <- paste(" -",hr_tab[3,1])
  #    }
  #  }
  #  else if (contCheck == FALSE) {
  #    pval2 <- get_lik_pval(tab)
  #    HR <- hr_tab[3,3]
  #    referenced <- paste(" -",hr_tab[3,1])
  #  }
  #  
  #  
  #  if (as.numeric(pval2) < 0.05) {
  #    pval_char <- "strongly associated"
  #  }
  #  if (as.numeric(pval2) >= 0.05 & as.numeric(pval2) < 0.1) {
  #    pval_char <- "moderately associated"
  #  }
  #  if (as.numeric(pval2) >= 0.1) {
  #    pval_char <- "not associated"
  #  }
  #  if (as.numeric(HR) > 1) {
  #    HR_char <- "high risk"
  #  }
  #  if (as.numeric(HR) <= 1) {
  #    HR_char <- "low risk"
  #  }
  #  
  #  if (FeatureSelec %in% StatCols) {
  #    FeatureSelec <- paste(FeatureSelec," (",geneset_name,")",sep = "")
  #  }
  #  
  #  line1 <- paste("<li><b>",SurvDateType,"</b> survival analysis ", SampleType,sep = "")
  #  line3 <- paste("<li>Kaplan-Meier survival curve categorized by <b>",FeatureSelec,"</b>.</li>",sep = "")
  #  line4 <- paste("<li>Cox hazard regression analysis finds a Likelihood Ratio P.value of <b>",pval2,"</b> and a Hazard Ratio of <b>",HR,"</b>, <b>",FeatureSelec,referenced,
  #                 "</b> is <b>",pval_char,"</b> with <b>",HR_char,"</b> for <b>",SurvDateType,"</b>.</li>",sep = "")
  #  
  #  if (is.null(line2)) {
  #    HTML(paste("<ul>",line1,line3,line4,"</ul>", sep = ""))
  #  }
  #  else if (!is.null(line2)) {
  #    HTML(paste("<ul>",line1,line2,line3,line4,"</ul>", sep = ""))
  #  }
  #  
  #})
  
  BivarForestPlot_react <- reactive({
    
    if (length(input$SurvivalFeatureBi1 > 0) & length(input$SurvivalFeatureBi2 > 0)) {
      
      tab <- BiVarAddTab_react()
      meta_ssgsea_sdf <- BiVarAddFeature_react()
      Feature1 <- input$SurvivalFeatureBi1
      Feature2 <- input$SurvivalFeatureBi2
      forextFont <- input$ForestFontSize
      
      
      forest <- ggforest(tab,
                         data = meta_ssgsea_sdf,
                         main = paste("Hazard Ratio Modeling: ",paste(Feature1,"+",Feature2,sep = ""),sep = ""),
                         fontsize = forextFont)
      forest
      
    }
    
  })
  
  output$BivarForestPlot <- renderPlot({
    
    forest <- BivarForestPlot_react()
    forest
    
  })
  
  BivarLinearityPlot_react <- reactive({
    
    if (length(input$SurvivalFeatureBi1 > 0) & length(input$SurvivalFeatureBi2 > 0)) {
      
      tab <- BiVarAddTab_react()
      residType <- input$ResidualTypeBi
      linpredict <- input$linPredict2
      Feature1 <- input$SurvivalFeatureBi1
      Feature2 <- input$SurvivalFeatureBi2
      AxisFont <- input$linAxisFont
      MainFont <- input$linMainFont
      TickFont <- input$linTickFont
      
      p <- ggcoxdiagnostics(tab,
                            type = residType,
                            sline = T,
                            sline.se = T,
                            ggtheme = theme_minimal(),
                            ox.scale = linpredict)
      p <- ggpar(p,
                 font.x = AxisFont,
                 font.y = AxisFont,
                 font.main = MainFont,
                 font.tickslab = TickFont,
                 main = paste("Linearity Plot Featuring: ",Feature1," + ",Feature2, sep = ""),
                 ylab = paste(str_to_title(residType)," Residuals", sep = "")
      )
      p
      
    }
    
  })
  
  output$BivarLinearityPlot <- renderPlot({
    
    p <- BivarLinearityPlot_react()
    p
    
  })
  
  ####----Bivariate Interactive----####
  
  BiVarIntFeature_react <- reactive({
    
    if (length(input$SurvivalFeatureBi1Inter > 0) & length(input$SurvivalFeatureBi2Inter> 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature1 <- input$SurvivalFeatureBi1Inter
      Feature2 <- input$SurvivalFeatureBi2Inter
      Feat1Var <- input$SurvFeatVariableBi1Inter
      Feat2Var <- input$SurvFeatVariableBi2Inter
      surv_time_col <- input$SurvivalType_time
      surv_id_col <- input$SurvivalType_id
      quantCutoff <- input$QuantPercent/100
      quantCutoff2 <- input$QuantPercent2/100
      meta_ssgsea <- ssGSEAmeta()
      
      ## Determine type of survival data - OS/EFS/PFS?
      SurvDateType <- sub("\\..*","",surv_time_col)
      
      ## Remove rows with NA in survival column
      meta_ssgsea <- meta_ssgsea[!is.na(meta_ssgsea[,surv_time_col]),]
      
      if (input$BiVarIntNAcheck1 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature1]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[which(meta_ssgsea[,Feature1] != "Inf"),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature1],ignore.case = T, invert = T),]
        ### Re-Perform Stat functions
        #meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        #meta_ssgsea$QuartileCutP <- paste("", meta_ssgsea$VAR_Q, sep="")
        #meta_ssgsea$MedianCutP <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        #meta_ssgsea$TopBottomCutP <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        #meta_ssgsea$UserCutP <- quantile_conversion2(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff2)
        
      }
      if (input$BiVarIntNAcheck2 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature2]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[which(meta_ssgsea[,Feature2] != "Inf"),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature2],ignore.case = T, invert = T),]
        ### Re-Perform Stat functions
        #meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        #meta_ssgsea$QuartileCutP <- paste("", meta_ssgsea$VAR_Q, sep="")
        #meta_ssgsea$MedianCutP <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        #meta_ssgsea$TopBottomCutP <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        #meta_ssgsea$UserCutP <- quantile_conversion2(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff2)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature1,Feature2)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      colnames(meta_ssgsea_sdf)[4] <- gsubCheck(colnames(meta_ssgsea_sdf)[4])
      if (!is.numeric(meta_ssgsea_sdf[,4])) {
        meta_ssgsea_sdf[,4] <- gsubCheck(meta_ssgsea_sdf[,4])
      }
      Feature1 <- gsubCheck(Feature1)
      Feat1Var <- gsubCheck(Feat1Var)
      colnames(meta_ssgsea_sdf)[5] <- gsubCheck(colnames(meta_ssgsea_sdf)[5])
      if (!is.numeric(meta_ssgsea_sdf[,5])) {
        meta_ssgsea_sdf[,5] <- gsubCheck(meta_ssgsea_sdf[,5])
      }
      Feature2 <- gsubCheck(Feature2)
      Feat2Var <- gsubCheck(Feat2Var)
      
      if (input$BiVarIntContCheck1 == FALSE) {
        meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
        meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      }
      else if (input$BiVarIntContCheck1 == TRUE) {
        if (input$BiVarIntContHiLoCheck1 == TRUE) {
          meta_ssgsea_sdf[,Feature1] <- highlow(as.numeric(meta_ssgsea_sdf[, which(colnames(meta_ssgsea_sdf) == Feature1)]))
          meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
          meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
        }
        else if (input$BiVarIntContHiLoCheck1 == FALSE) {
          meta_ssgsea_sdf[,Feature1] <- as.numeric(meta_ssgsea_sdf[,Feature1])
        }
      }
      if (input$BiVarIntContCheck2 == FALSE) {
        meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
        meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      }
      else if (input$BiVarIntContCheck2 == TRUE) {
        if (input$BiVarIntContHiLoCheck2 == TRUE) {
          meta_ssgsea_sdf[,Feature2] <- highlow(as.numeric(meta_ssgsea_sdf[, which(colnames(meta_ssgsea_sdf) == Feature2)]))
          meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
          meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
        }
        else if (input$BiVarIntContHiLoCheck2 == FALSE) {
          meta_ssgsea_sdf[,Feature2] <- as.numeric(meta_ssgsea_sdf[,Feature2])
        }
      }
    }
    
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "ID"
    
    meta_ssgsea_sdf
    
  })
  
  BiVarIntTab_react <- reactive({
    
    meta_ssgsea_sdf <- BiVarIntFeature_react()
    #Feature1 <- input$SurvivalFeatureBi1Inter
    #Feature2 <- input$SurvivalFeatureBi2Inter
    Feature1 <- colnames(meta_ssgsea_sdf)[4]
    Feature2 <- colnames(meta_ssgsea_sdf)[5]
    
    form <- paste("Surv(time,ID) ~ ",paste(Feature1,"*",Feature2,sep = ""),sep = "")
    form2 <- as.formula(form)
    tab <- eval(substitute(coxph(form2,data = meta_ssgsea_sdf)))
    
    tab
    
    
  })
  
  BiVarIntTab4Annova_react <- reactive({
    
    meta_ssgsea_sdf <- BiVarIntFeature_react()
    #Feature1 <- input$SurvivalFeatureBi1Inter
    #Feature2 <- input$SurvivalFeatureBi2Inter
    Feature1 <- colnames(meta_ssgsea_sdf)[4]
    Feature2 <- colnames(meta_ssgsea_sdf)[5]
    
    form <- paste("Surv(time,ID) ~ ",paste(Feature1,"+",Feature2,sep = ""),sep = "")
    form2 <- as.formula(form)
    tab <- eval(substitute(coxph(form2,data = meta_ssgsea_sdf)))
    
    tab
    
    
  })
  
  BiVarIntHRTab_react <- reactive({
    
    tab <- BiVarIntTab_react()
    tab <- tab %>% 
      gtsummary::tbl_regression(exp = TRUE) %>%
      as_gt()
    
    tab_df <- as.data.frame(tab)
    
    tab_df[is.na(tab_df)] <- ""
    tab_df <- tab_df %>%
      select(label,n_obs,estimate,std.error,ci,p.value)
    colnames(tab_df) <- c("Variable","N","Hazard Ratio","Std. Error","95% Confidence Interval","P.Value")
    
    tab_df
    
  })
  
  output$BiFeatureHRtabInter <- renderTable({
    
    tab <- BiVarIntHRTab_react()
    tab
    
  })
  
  output$bivarSummaryInter <- renderPrint({
    
    tab <- BiVarIntTab_react()
    out <- capture.output(summary(tab))
    
    con_line <- grep("^Concordance=",out,value = T)
    lik_line <- grep("^Likelihood ratio test=",out,value = T)
    wal_line <- grep("^Wald test",out,value = T)
    sco_line <- grep("^Score ",out,value = T)
    
    text <- paste("Coxh Summary:",con_line,lik_line,wal_line,sco_line,sep = "\n")
    cat(text)
    
  })
  
  output$bivarAnovaInter1 <- renderPrint({
    
    tab1 <- BiVarIntTab_react()
    tab2 <- BiVarIntTab4Annova_react()
    
    annova_res <- anova(tab1,tab2)
    
    out <- capture.output(annova_res)
    
    line1 <- out[3]
    line2 <- out[4]
    line3 <- out[5]
    line4 <- out[6]
    line5 <- out[7]
    
    text <- paste("Model Comparison:",line1,line2,line3,line4,line5,sep = "\n")
    cat(text)
    
  })
  
  ## Survival Plot - TWO FEATURE
  featSplotBi_react <- reactive({
    
    if (length(input$SurvivalFeatureBi1Inter > 0) & length(input$SurvivalFeatureBi2Inter> 0)) {
      
      ## Assign variables
      SampleType <- input$SampleTypeSelection
      Feature1_lab <- input$SurvivalFeatureBi1Inter
      Feature2_lab <- input$SurvivalFeatureBi2Inter
      show_pval <- input$ShowPval
      ShowConfInt <- input$ShowConfInt
      xaxlim <- input$SurvXaxis * 365.25
      showLegend <- input$SurvLegendPos
      surv_time_col <- input$SurvivalType_time
      surv_id_col <- input$SurvivalType_id
      meta_ssgsea_sdf <- BiVarIntFeature_react()
      Feature1 <- colnames(meta_ssgsea_sdf)[4]
      Feature2 <- colnames(meta_ssgsea_sdf)[5]
      showMedSurv <- input$ShowMedSurvLine
      if (showMedSurv == T) {
        showMedSurv <- "hv"
      }
      else if (showMedSurv == F) {
        showMedSurv <- "none"
      }
      
      ## Determine type of survival data - OS/EFS/PFS?
      SurvDateType <- sub("\\..*","",surv_time_col)
      
      form <- paste("Surv(time,ID) ~ ",paste(Feature1,"+",Feature2,sep = ""),sep = "")
      form2 <- as.formula(form)
      fit <- eval(substitute(survfit(form2,data = meta_ssgsea_sdf, type="kaplan-meier")))
      
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleTypeLab <- paste(" (",SampleType,") ",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleTypeLab <- " "
      }
      
      ## Determind Plot title
      if (is.null(input$SurvPlotTitleBiVar)) {
        SurvPlotTitle <- paste("Survival curves of ",Feature1_lab," Catagorized by\n",Feature2_lab,".", sep = "")
        #SurvPlotTitle <- paste("Survival curves of ",Feature1," and\n",Feature2," in",SampleTypeLab,"Patients", sep = "")
      }
      else if (!is.null(input$SurvPlotTitleBiVar)) {
        if (input$SurvPlotTitleBiVar == "") {
          SurvPlotTitle <- paste("Survival curves of ",Feature1_lab," Catagorized by\n",Feature2_lab,".", sep = "")
          #SurvPlotTitle <- paste("Survival curves of ",Feature1," and\n",Feature2," in",SampleTypeLab,"Patients", sep = "")
        }
        else if (input$SurvPlotTitleBiVar != "") {
          SurvPlotTitle <- input$SurvPlotTitleBiVar
        }
      }
      
      ## Generate plot
      ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                           title = SurvPlotTitle,
                           xscale = c("d_y"),
                           break.time.by=365.25,
                           xlab = "Years", 
                           ylab = paste(SurvDateType,"Survival Probability"),
                           submain = "Based on Kaplan-Meier estimates",
                           caption = "created with survminer",
                           pval=show_pval,
                           conf.int = ShowConfInt,
                           ggtheme = theme_bw(),
                           font.title = c(16, "bold"),
                           font.submain = c(12, "italic"),
                           font.caption = c(12, "plain"),
                           font.x = c(14, "plain"),
                           font.y = c(14, "plain"),
                           font.tickslab = c(12, "plain"),
                           legend = showLegend,
                           risk.table.height = 0.20,
                           surv.median.line = showMedSurv
      )
      if (showMedSurv != "none") {
        MedSurvItem <- ggsurv[["plot"]][["layers"]][length(ggsurv[["plot"]][["layers"]])]
        MedSurvItem_df <- MedSurvItem[[1]][["data"]]
        MedSurvItem_df <- MedSurvItem_df[order(MedSurvItem_df[,1]),]
        MedSurvItem_df <- MedSurvItem_df %>%
          mutate(label = paste(round(MedSurvItem_df[,1]),"Days"))
        rownames(MedSurvItem_df) <- 1:nrow(MedSurvItem_df)
        if (nrow(MedSurvItem_df) > 1) {
          ggsurv$plot <- ggsurv$plot +
            geom_label_repel(data = MedSurvItem_df, aes(x = x1, y = y1, label = label, size = 4), label.size = NA, show.legend = FALSE)
        }
      }
      if (!is.null(input$SurvXaxis)) {
        ggsurv$plot$coordinates$limits$x <- c(0,xaxlim)
        ggsurv$table$coordinates$limits$x <- c(0,xaxlim)
      }
      
      ggsurv$table <- ggsurv$table + theme_cleantable()
      ggsurv
      
    }
    
  })
  
  output$featSplotBi <- renderPlot({
    plot <- featSplotBi_react()
    plot
  })
  
  BivarForestPlotInter_react <- reactive({
    
    if (length(input$SurvivalFeatureBi1Inter > 0) & length(input$SurvivalFeatureBi2Inter > 0)) {
      
      tab <- BiVarIntTab_react()
      meta_ssgsea_sdf <- BiVarIntFeature_react()
      Feature1 <- input$SurvivalFeatureBi1Inter
      Feature2 <- input$SurvivalFeatureBi2Inter
      forextFont <- input$ForestFontSize
      
      forest <- ggforest(tab,
                         data = meta_ssgsea_sdf,
                         main = paste("Hazard Ratio Modeling: ",paste(Feature1,"*",Feature2,sep = ""),sep = ""),
                         fontsize = forextFont)
      forest
      
    }
    
  })
  
  output$BivarForestPlotInter <- renderPlot({
    
    forest <- BivarForestPlotInter_react()
    forest
    
  })
  
  
  BivarLinearityPlotInter_react <- reactive({
    
    if (length(input$SurvivalFeatureBi1Inter > 0) & length(input$SurvivalFeatureBi2Inter > 0)) {
      
      tab <- BiVarIntTab_react()
      residType <- input$ResidualTypeInter
      linpredict <- input$linPredict3
      Feature1 <- input$SurvivalFeatureBi1Inter
      Feature2 <- input$SurvivalFeatureBi2Inter
      AxisFont <- input$linAxisFont
      MainFont <- input$linMainFont
      TickFont <- input$linTickFont
      
      p <- ggcoxdiagnostics(tab,
                            type = residType,
                            sline = T,
                            sline.se = T,
                            ggtheme = theme_minimal(),
                            ox.scale = linpredict)
      p <- ggpar(p,
                 font.x = AxisFont,
                 font.y = AxisFont,
                 font.main = MainFont,
                 font.tickslab = TickFont,
                 main = paste("Linearity Plot Featuring: ",Feature1," * ",Feature2, sep = ""),
                 ylab = paste(str_to_title(residType)," Residuals", sep = "")
      )
      p
      
    }
    
  })
  
  output$BivarLinearityPlotInter <- renderPlot({
    
    p <- BivarLinearityPlotInter_react()
    p
    
  })
  
  ####----Multivariate----####
  
  MultiVarFeat_react <- reactive({
    
    if (length(input$SurvivalFeature > 0)) {
      
      Feature <- input$SurvivalFeature
      surv_time_col <- input$SurvivalType_time
      surv_id_col <- input$SurvivalType_id
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      quantCutoff2 <- input$QuantPercent2/100 #Quantile cutoff given by user
      meta_ssgsea <- ssGSEAmeta()
      geneset <- gs_react()
      geneset_name <- names(geneset)
      
      if (input$UniVarNAcheck == TRUE) {
        # Remove NA_unknown
        for (i in Feature) {
          meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,i]) == FALSE),]
          meta_ssgsea <- meta_ssgsea[which(meta_ssgsea[,i] != "Inf"),]
          meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,i],ignore.case = T, invert = T),]
        }
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      meta_ssgsea_sdf
      
    }
    
  })
  
  MultiVarFeatCat_react <- reactive({
    
    meta_ssgsea_sdf <- MultiVarFeat_react()
    Feature <- input$SurvivalFeature
    surv_time_col <- input$SurvivalType_time
    surv_id_col <- input$SurvivalType_id
    
    for (i in Feature){
      meta_ssgsea_sdf[,i] <- as.factor(meta_ssgsea_sdf[,i])
      colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == i)] <- gsubCheck(colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == i)])
    }
    
    Feature <- gsubCheck(Feature)
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "ID"
    
    meta_ssgsea_sdf
    
  })
  
  MultiVarFeatCont_react <- reactive({
    
    meta_ssgsea_sdf <- MultiVarFeat_react()
    Feature <- input$SurvivalFeature
    surv_time_col <- input$SurvivalType_time
    surv_id_col <- input$SurvivalType_id
    
    for (i in Feature){
      colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == i)] <- gsubCheck(colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == i)])
    }
    
    Feature <- gsubCheck(Feature)
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "ID"
    
    meta_ssgsea_sdf
    
  })
  
  MultiVarTabCat_react <- reactive({
    
    meta_ssgsea_sdf <- MultiVarFeatCat_react()
    Feature <- input$SurvivalFeature
    Feature <- gsubCheck(Feature)
    
    form <- paste("Surv(time,ID) ~ ",paste(Feature,collapse = "+"),sep = "")
    form2 <- as.formula(form)
    tab <- eval(substitute(coxph(form2,data = meta_ssgsea_sdf)))
    
    
    tab
    
  })
  
  MultiVarTabCont_react <- reactive({
    
    meta_ssgsea_sdf <- MultiVarFeatCont_react()
    Feature <- input$SurvivalFeature
    Feature <- gsubCheck(Feature)
    
    form <- paste("Surv(time,ID) ~ ",paste(Feature,collapse = "+"),sep = "")
    form2 <- as.formula(form)
    tab <- eval(substitute(coxph(form2,data = meta_ssgsea_sdf)))
    
    
    tab
    
  })
  
  SFeatureHRtabCat_react <- reactive({
    
    if (length(input$SurvivalFeature > 0)) {
      
      tab <- MultiVarTabCat_react()
      
      tab <- tab %>% 
        gtsummary::tbl_regression(exp = TRUE) %>%
        as_gt()
      
      tab_df <- as.data.frame(tab)
      
      tab_df[is.na(tab_df)] <- ""
      tab_df <- tab_df %>%
        select(label,n_obs,estimate,std.error,ci,p.value)
      colnames(tab_df) <- c("Variable","N","Hazard Ratio","Std. Error","95% Confidence Interval","P.Value")
      
      tab_df
      
    }
    
  })
  
  output$SFeatureHRtabCat <- renderTable({
    
    tab <- SFeatureHRtabCat_react()
    tab
    
  })
  
  SFeatureHRtabCont_react <- reactive({
    
    if (length(input$SurvivalFeature > 0)) {
      
      tab <- MultiVarTabCont_react()
      
      tab <- tab %>% 
        gtsummary::tbl_regression(exp = TRUE) %>%
        as_gt()
      
      tab_df <- as.data.frame(tab)
      
      tab_df[is.na(tab_df)] <- ""
      tab_df <- tab_df %>%
        select(label,n_obs,estimate,std.error,ci,p.value)
      colnames(tab_df) <- c("Variable","N","Hazard Ratio","Std. Error","95% Confidence Interval","P.Value")
      
      tab_df
      
    }
    
  })
  
  output$SFeatureHRtabCont <- renderTable({
    
    tab <- SFeatureHRtabCont_react()
    tab
    
  })
  
  output$multivarSummaryCat <- renderPrint({
    
    tab <- MultiVarTabCat_react()
    out <- capture.output(summary(tab))
    
    con_line <- grep("^Concordance=",out,value = T)
    lik_line <- grep("^Likelihood ratio test=",out,value = T)
    wal_line <- grep("^Wald test",out,value = T)
    sco_line <- grep("^Score ",out,value = T)
    
    text <- paste("Coxh Summary (Categorical):",con_line,lik_line,wal_line,sco_line,sep = "\n")
    cat(text)
    
  })
  
  output$multivarSummaryCont <- renderPrint({
    
    tab <- MultiVarTabCont_react()
    out <- capture.output(summary(tab))
    
    con_line <- grep("^Concordance=",out,value = T)
    lik_line <- grep("^Likelihood ratio test=",out,value = T)
    wal_line <- grep("^Wald test",out,value = T)
    sco_line <- grep("^Score ",out,value = T)
    
    text <- paste("Coxh Summary (Categorical):",con_line,lik_line,wal_line,sco_line,sep = "\n")
    cat(text)
    
  })
  
  MultivarForestPlot_react <- reactive({
    
    if (length(input$SurvivalFeature > 0)) {
      
      tab <- MultiVarTabCat_react()
      meta_ssgsea_sdf <- MultiVarFeat_react()
      Feature <- input$SurvivalFeature
      forextFont <- input$ForestFontSize
      
      forest <- ggforest(tab,
                         data = meta_ssgsea_sdf,
                         main = paste("Hazard Ratio Modeling: ",paste(Feature,collapse = ", "),sep = ""),
                         fontsize = forextFont)
      forest
      
    }
    
  })
  
  output$MultivarForestPlot <- renderPlot({
    
    forest <- MultivarForestPlot_react()
    forest
    
  })
  
  
  ####----Risk Stratification----####
  
  ## Boxplot reactive to determine high/low risk samples based on user input
  SboxplotReact <- reactive({
    
    ## Assing variables
    cutoff_1 <- input$cutoffTime1                   #High-risk cutoff time
    cutoff_0 <- input$cutoffTime0                   #Low-risk cutoff time
    OS_choice_1 <- input$survStatus1                #High-risk survival status
    OS_choice_0 <- input$survStatus0                #Low-risk survival status
    ssGSEA_meta <- ssGSEAmeta()                     #Meta with ssGSEA scores and function calculations
    surv_time_col <- input$SurvivalType_time
    surv_id_col <- input$SurvivalType_id
    
    # If both survival choices are numeric (1 or 0)
    if (is.na(as.numeric(OS_choice_1)) == F & is.na(as.numeric(OS_choice_0)) == F) {
      ssGSEA_meta <- ssGSEA_meta %>%
        mutate(SurvivalCutoff = case_when(
          ssGSEA_meta[,surv_time_col] <= cutoff_1 & ssGSEA_meta[,surv_id_col] == as.numeric(OS_choice_1) ~ "High-Risk [Below Survival Time Cutoff]",
          ssGSEA_meta[,surv_time_col] >= cutoff_0 & ssGSEA_meta[,surv_id_col] == as.numeric(OS_choice_0) ~ "Low-Risk [Above Survival Time Cutoff]"
        ))
    }
    # If either choice selects "Both" as survival outcome
    if (is.na(as.numeric(OS_choice_1)) == T | is.na(as.numeric(OS_choice_0)) == T) {
      # if both & (1 or 0)
      if (is.na(as.numeric(OS_choice_1)) == T & is.na(as.numeric(OS_choice_0)) == F) {
        ssGSEA_meta <- ssGSEA_meta %>%
          mutate(SurvivalCutoff = case_when(
            ssGSEA_meta[,surv_time_col] <= cutoff_1 ~ "High-Risk [Below Survival Time Cutoff]",
            ssGSEA_meta[,surv_time_col] >= cutoff_0 & ssGSEA_meta[,surv_id_col] == as.numeric(OS_choice_0) ~ "Low-Risk [Above Survival Time Cutoff]"
          ))
      }
      # if (1 or 0) and both
      else if (is.na(as.numeric(OS_choice_1)) == F & is.na(as.numeric(OS_choice_0)) == T) {
        ssGSEA_meta <- ssGSEA_meta %>%
          mutate(SurvivalCutoff = case_when(
            ssGSEA_meta[,surv_time_col] <= cutoff_1 & ssGSEA_meta[,surv_id_col] == as.numeric(OS_choice_1) ~ "High-Risk [Below Survival Time Cutoff]",
            ssGSEA_meta[,surv_time_col] >= cutoff_0 ~ "Low-Risk [Above Survival Time Cutoff]"
          ))
      }
      # if both and both
      else if (is.na(as.numeric(OS_choice_1)) == T & is.na(as.numeric(OS_choice_0)) == T) {
        ssGSEA_meta <- ssGSEA_meta %>%
          mutate(SurvivalCutoff = case_when(
            ssGSEA_meta[,surv_time_col] <= cutoff_1 ~ "High-Risk [Below Survival Time Cutoff]",
            ssGSEA_meta[,surv_time_col] >= cutoff_0 ~ "Low-Risk [Above Survival Time Cutoff]"
          ))
      }
    }
    
    ssGSEA_meta <- ssGSEA_meta[which(is.na(ssGSEA_meta$SurvivalCutoff) == F),]
    ssGSEA_meta$SurvivalCutoff <- factor(ssGSEA_meta$SurvivalCutoff,
                                         levels = c("High-Risk [Below Survival Time Cutoff]","Low-Risk [Above Survival Time Cutoff]"))
    ssGSEA_meta
    
  })
  
  Sboxplot_react <- reactive({
    
    ## Assign Variables
    ssGSEA_meta <- SboxplotReact()
    geneset <- gs_react()
    GeneSet <- names(geneset)
    Feature <- input$FeatureSelection
    dot <- input$boxplotDot
    font <- input$boxplotFont
    SampleType <- input$SampleTypeSelection
    scoreMethod <- input$ScoreMethod
    logchoice <- input$SBoxLog
    surv_time_col <- input$SurvivalType_time
    surv_id_col <- input$SurvivalType_id
    
    ## Determine type of survival data - OS/EFS/PFS?
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## Adjust 'Sample Type' for label 
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      SampleTypeLab <- paste(" (",SampleType,") ",sep = "")
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      SampleTypeLab <- " "
    }
    if (GeneSet %in% decon_score_cols) {
      scoreMethod <- "Pre-Processed Score"
    }
    if (input$GeneSetTabs == 2) {
      scoreMethod <- "Gene Expression Score"
      #if (input$RawOrSS == "Raw Gene Expression") {
      #  scoreMethod <- "Raw Gene Expression"
      #}
    }
    
    if (input$GeneSetTabs != 2 & !(GeneSet %in% decon_score_cols)) {
      scoreMethod <- paste(scoreMethod,"Score")
    }
    
    if (logchoice == TRUE) {
      
      ssGSEA_meta[,GeneSet] <- log(ssGSEA_meta[,GeneSet] + 1)
      scoreMethod <- paste(scoreMethod,"(log(score + 1))")
      
    }
    
    plot <- ggplot(ssGSEA_meta, aes(factor(SurvivalCutoff), ssGSEA_meta[,GeneSet], fill = SurvivalCutoff)) +
      geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
      geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = dot) +
      labs(x = "Group", y = scoreMethod,
           title = paste(GeneSet," ",scoreMethod,": ",Feature,SampleTypeLab,"Patients",sep = "")) +
      theme_bw() +
      stat_compare_means(method = input$boxoptselec) +
      theme(text = element_text(size = font))
    plot
    
  })
  
  output$Sboxplot <- renderPlot({
    
    plot <- Sboxplot_react()
    plot
    
    
  })
  
  output$heatmap_error_message <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    if (geneset_name %in% decon_score_cols | input$GeneSetTabs == 2) {
      p("Heatmap not available for Pre-Processed scores or single gene analysis")
    }
    
  })
  output$heatmap_error_message2 <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    if (geneset_name %in% decon_score_cols | input$GeneSetTabs == 2) {
      p("Heatmap not available for Pre-Processed scores or single gene analysis")
    }
    
  })
  
  Sheatmap_react  <- reactive({
    
    geneset <- gs_react()
    GeneSet <- names(geneset)
    heatgenes <- geneset[[GeneSet]]
    meta <- SboxplotReact()
    expr_start <- exprSub()
    samples <- meta$SampleName
    expr <- expr_start[which(rownames(expr_start) %in% heatgenes),colnames(expr_start) %in% samples, drop = F]
    meta <- meta[which(meta$SampleName %in% colnames(expr)),]
    clmethod <- input$ClusterMethod
    rowfont <- input$heatmapFontR
    colfont <- input$heatmapFontC
    color_choice <- input$ColorPaletteHeat
    
    dataset <- expr
    dataset <- log2(dataset + 1)
    zdataset <- apply(dataset, 1, scale)
    zdataset <- apply(zdataset, 1, rev)
    colnames(zdataset) <- names(dataset)
    dataset <- as.matrix(zdataset)
    dataset[is.na(dataset)] <- 0
    dataset = dataset[apply(dataset[,-1], 1, function(x) !all(x==0)),]
    minimum = -5;
    maximum = 5;
    if (abs(min(dataset)) > abs(max(dataset))) {
      dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
    } else {
      dataset[dataset > abs(min(dataset))] = abs(min(dataset))
    }
    meta2 <- meta[order(meta[,2]),]
    meta3 <- meta2[order(meta2$SurvivalCutoff),]
    samporder <- meta3$SampleName
    dataset2 <- dataset[,samporder]
    type <- meta3$SurvivalCutoff
    meta4 <- as.data.frame(type)
    rownames(meta4) <- meta3[,1]
    bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
    #hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
    #Heatmap color
    col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
    if (color_choice == "original") {
      HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (color_choice %in% col_sets) {
      HeatMap_Colors <- brewer.pal(n = 5, color_choice)
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (color_choice == "Inferno") {
      hmcols <- inferno(500)
    }
    else if (color_choice == "Viridis") {
      hmcols <- viridis(500)
    }
    else if (color_choice == "Plasma") {
      hmcols <- plasma(500)
    }
    else if (color_choice == "OmniBlueRed") {
      hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
    }
    else if (color_choice == "LightBlueBlackRed") {
      hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
    }
    else if (color_choice == "GreenBlackRed") {
      hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
    }
    heat <- pheatmap(dataset2,
                     cluster_col = F,
                     cluster_row = T,
                     fontsize_row = rowfont,
                     fontsize_col = colfont,
                     show_rownames = T ,
                     show_colnames = T,
                     annotation_col = meta4,
                     clustering_method = clmethod,
                     color=hmcols,
                     border_color = NA)
    heat
    
  })
  
  output$Sheatmap <- renderPlot({
    
    heat <- Sheatmap_react()
    heat
    
    
  })
  
  ####----Feature Stratification----####
  
  Featureboxplot_react <- reactive({
    
    SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    scoreMethod <- input$ScoreMethod
    dot <- input$boxplotDot
    font <- input$boxplotFont
    StatMethod <- input$boxoptselec
    geneset <- gs_react()
    GeneSet <- names(geneset)
    FeatureSelec <- input$BoxplotFeature
    boxplotang <- input$boxplotTextAngle
    meta_ssGSEA <- ssGSEAmeta()
    logchoice <- input$FBoxLog
    boxTab <- meta_ssGSEA[,c("SampleName",FeatureSelec,GeneSet)]
    if (input$BoxPRemoveNA == TRUE) {
      # Remove NA_unknown
      boxTab <- boxTab[which(is.na(boxTab[,FeatureSelec]) == FALSE),]
      boxTab <- boxTab[which(boxTab[,FeatureSelec] != "Inf"),]
      boxTab <- boxTab[grep("unknown",boxTab[,FeatureSelec],ignore.case = T, invert = T),]
    }
    boxTab[,FeatureSelec] <- as.factor(boxTab[,FeatureSelec])
    
    
    ## Adjust 'Sample Type' for label 
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      SampleTypeLab <- paste(" (",SampleType,") ",sep = "")
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      SampleTypeLab <- " "
    }
    if (GeneSet %in% decon_score_cols) {
      scoreMethod <- "Pre-Processed Score"
    }
    if (input$GeneSetTabs == 2) {
      scoreMethod <- "Gene Expression Score"
      #if (input$RawOrSS == "Raw Gene Expression") {
      #  scoreMethod <- "Raw Gene Expression"
      #}
    }
    if (input$GeneSetTabs != 2 & !(GeneSet %in% decon_score_cols)) {
      scoreMethod <- paste(scoreMethod,"Score")
    }
    
    if (logchoice == TRUE) {
      
      boxTab[,GeneSet] <- log(boxTab[,GeneSet] + 1)
      scoreMethod <- paste(scoreMethod,"(log(score + 1))")
      
    }
    
    
    if (is.na(as.numeric(boxplotang)) == TRUE) {
      ggplot(boxTab, aes(factor(boxTab[,FeatureSelec]), boxTab[,GeneSet], fill = boxTab[,FeatureSelec])) +
        geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
        geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = dot) +
        labs(x = "Group", y = scoreMethod,
             title = paste(GeneSet," ",scoreMethod,": ","\n",
                           Feature,SampleTypeLab,"Patients Featuring ",FeatureSelec,sep = ""),
             fill = FeatureSelec) +
        theme_bw() +
        stat_compare_means(method = StatMethod) +
        theme(text = element_text(size = font)) +
        scale_x_discrete(guide = guide_axis(n.dodge = 2))
    }
    
    else if (is.na(as.numeric(boxplotang)) == FALSE) {
      if (as.numeric(boxplotang) == 0) {
        ggplot(boxTab, aes(factor(boxTab[,FeatureSelec]), boxTab[,GeneSet], fill = boxTab[,FeatureSelec])) +
          geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
          geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = dot) +
          labs(x = "Group", y = scoreMethod,
               title = paste(GeneSet," ",scoreMethod,": ","\n",
                             Feature,SampleTypeLab,"Patients Featuring ",FeatureSelec,sep = ""),
               fill = FeatureSelec) +
          theme_bw() +
          stat_compare_means(method = StatMethod) +
          theme(text = element_text(size = font))
      }
      else if (as.numeric(boxplotang) == 45) {
        ggplot(boxTab, aes(factor(boxTab[,FeatureSelec]), boxTab[,GeneSet], fill = boxTab[,FeatureSelec])) +
          geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
          geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = dot) +
          labs(x = "Group", y = scoreMethod,
               title = paste(GeneSet," ",scoreMethod,": ","\n",
                             Feature,SampleTypeLab,"Patients Featuring ",FeatureSelec,sep = ""),
               fill = FeatureSelec) +
          theme_bw() +
          stat_compare_means(method = StatMethod) +
          theme(text = element_text(size = font),
                axis.text.x = element_text(angle = as.numeric(boxplotang), hjust = 1))
      }
      else if (as.numeric(boxplotang) == 90) {
        ggplot(boxTab, aes(factor(boxTab[,FeatureSelec]), boxTab[,GeneSet], fill = boxTab[,FeatureSelec])) +
          geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
          geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = dot) +
          labs(x = "Group", y = scoreMethod,
               title = paste(GeneSet," ",scoreMethod,": ","\n",
                             Feature,SampleTypeLab,"Patients Featuring ",FeatureSelec,sep = ""),
               fill = FeatureSelec) +
          theme_bw() +
          stat_compare_means(method = StatMethod) +
          theme(text = element_text(size = font),
                axis.text.x = element_text(angle = as.numeric(boxplotang)))
      }
      
    }
    
  })
  
  ## Feature Boxplot
  output$Featureboxplot <- renderPlot({
    
    plot <- Featureboxplot_react()
    plot
    
  })
  
  FeatureHeatmap_react <- reactive({
    
    geneset <- gs_react()
    GeneSet <- names(geneset)
    FeatureSelec <- input$HeatmapFeature
    meta_ssGSEA <- ssGSEAmeta()
    meta <- meta_ssGSEA[,c("SampleName",FeatureSelec,GeneSet)]
    if (input$HeatRemoveNA == TRUE) {
      # Remove NA_unknown
      meta <- meta[which(is.na(meta[,FeatureSelec]) == FALSE),]
      meta <- meta[which(meta[,FeatureSelec] != "Inf"),]
      meta <- meta[grep("unknown",meta[,FeatureSelec],ignore.case = T, invert = T),]
    }
    
    heatgenes <- geneset[[GeneSet]]
    expr_start <- exprSub()
    samples <- meta$SampleName
    expr <- expr_start[which(rownames(expr_start) %in% heatgenes),colnames(expr_start) %in% samples, drop = F]
    meta <- meta[which(meta$SampleName %in% colnames(expr)),]
    clmethod <- input$ClusterMethod
    rowfont <- input$heatmapFontR
    colfont <- input$heatmapFontC
    color_choice <- input$ColorPaletteHeat
    
    dataset <- expr
    dataset <- as.data.frame(log2(dataset + 1))
    zdataset <- as.data.frame(apply(dataset, 1, scale))
    zdataset <- as.data.frame(apply(zdataset, 1, rev))
    colnames(zdataset) <- names(dataset)
    dataset <- as.matrix(zdataset)
    dataset[is.na(dataset)] <- 0
    dataset = dataset[apply(dataset[,-1], 1, function(x) !all(x==0)),]
    minimum = -5;
    maximum = 5;
    if (abs(min(dataset)) > abs(max(dataset))) {
      dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
    } else {
      dataset[dataset > abs(min(dataset))] = abs(min(dataset))
    }
    #meta2 <- meta[order(meta[,2]),]
    meta3 <- meta[order(meta[,FeatureSelec]),]
    meta3[,FeatureSelec] <- as.character(meta3[,FeatureSelec])
    meta3[,FeatureSelec] = meta3[,FeatureSelec] %>% replace_na("NA")
    samporder <- meta3$SampleName
    dataset2 <- dataset[,samporder]
    type <- meta3[,FeatureSelec]
    meta4 <- as.data.frame(type)
    rownames(meta4) <- meta3[,1]
    bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
    #hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
    #Heatmap color
    col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
    if (color_choice == "original") {
      HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (color_choice %in% col_sets) {
      HeatMap_Colors <- brewer.pal(n = 5, color_choice)
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (color_choice == "Inferno") {
      hmcols <- inferno(500)
    }
    else if (color_choice == "Viridis") {
      hmcols <- viridis(500)
    }
    else if (color_choice == "Plasma") {
      hmcols <- plasma(500)
    }
    else if (color_choice == "OmniBlueRed") {
      hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
    }
    else if (color_choice == "LightBlueBlackRed") {
      hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
    }
    else if (color_choice == "GreenBlackRed") {
      hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
    }
    heat <- pheatmap(dataset2,
                     cluster_col = F,
                     cluster_row = T,
                     fontsize_row = rowfont,
                     fontsize_col = colfont,
                     show_rownames = T ,
                     show_colnames = T,
                     annotation_col = meta4,
                     clustering_method = clmethod,
                     color=hmcols,
                     border_color = NA)
    heat
    
  })
  
  output$FeatureHeatmap <- renderPlot({
    
    heat <- FeatureHeatmap_react()
    heat
    
  })
  
  ####----Density----####
  
  ssgseaDensity_react <- reactive({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    ssgsea_meta <- ssGSEAmeta()
    user_quant <- input$densityPercent/100
    ShowQuartile <- input$QuartileLinesCheck
    scoreMethod <- input$ScoreMethod
    
    cols_selec <- c("SampleName",geneset_name)
    ssgsea_scores <- ssgsea_meta[,cols_selec]
    
    quant_df <- data.frame(quantile(ssgsea_scores[,geneset_name],na.rm = T))
    quant_df2 <- quant_df[c(2,3,4),,drop = F]
    colnames(quant_df2)[1] <- "Quantile"
    
    user_vline <- quantile(ssgsea_scores[,geneset_name],probs = user_quant,na.rm = T)
    
    ## get score method for x and y labels
    if (input$GeneSetTabs == 2) {
      scoreMethodLab <- "Gene Expression Density"
      scoreMethodLab_x <- "Gene Expression"
      #if (input$RawOrSS == "Raw Gene Expression") {
      #  scoreMethodLab <- "Raw Gene Expression Density"
      #  scoreMethodLab_x <- "Raw Gene Expression"
      #}
      #else if (input$RawOrSS == "Rank Normalized") {
      #  scoreMethodLab <- paste(scoreMethod, " Score Density", sep = "")
      #  scoreMethodLab_x <- paste(scoreMethod, " Score", sep = "")
      #}
    }
    else if (input$GeneSetTabs != 2) {
      if (geneset_name %in% decon_score_cols) {
        scoreMethodLab <- "Pre-Processed Score Density"
        scoreMethodLab_x <- "Pre-Processed Score"
      }
      else {
        scoreMethodLab <- paste(scoreMethod, " Score Density", sep = "")
        scoreMethodLab_x <- paste(scoreMethod, " Score", sep = "")
      }
    }
    ## generate title based on input
    if (geneset_name %in% decon_score_cols) {
      dens_title <- paste(colnames(ssgsea_scores)[2]," Pre-Processed Score Density",sep = "")
    }
    else {
      dens_title <- paste(colnames(ssgsea_scores)[2],scoreMethodLab)
    }
    
    p <- ggplot(ssgsea_scores, aes(x=ssgsea_scores[,geneset_name])) + 
      geom_density(color="darkblue", fill="lightblue", alpha = 0.4) +
      xlab(scoreMethodLab_x) +
      ylab(scoreMethodLab) +
      #ggtitle(paste(colnames(ssgsea_scores)[2],scoreMethodLab)) +
      ggtitle(dens_title) +
      theme(axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            plot.title = element_text(size = 20))
    if (ShowQuartile == TRUE) {
      p <- p + geom_vline(data = quant_df2, aes(xintercept = Quantile), linetype = "dashed", color = "darkblue", size = 1)
    }
    if (user_quant != 0) {
      p <- p + geom_vline(xintercept = user_vline, linetype = "dashed", color = "darkred", size = 1)
    }
    p
    
  })
  
  output$ssgseaDensity <- renderPlot({
    
    p <- ssgseaDensity_react()
    p
    
    
  })
  
  ####----Feature Comparison----####
  
  FeatCompScatter_react <- reactive({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    Feature <- input$ScatterFeature
    ColorCol <- input$ScatterColor
    LogChoice <- input$ScatterLog
    meta <- ssGSEAmeta()
    
    if (input$ColorScatterChoice == "Feature") {
      scores <- meta[,c("SampleName",Feature,geneset_name,ColorCol)]
    }
    else if (input$ColorScatterChoice == "Single Color") {
      scores <- meta[,c("SampleName",Feature,geneset_name)]
      scores$ColorColumn <- ColorCol
    }
    
    scores[,4] <- as.factor(scores[,4])
    
    if (length(LogChoice) == 1) {
      if (LogChoice == "Log x-axis"){
        scores[,2] <- log(scores[,2] + 1)
      }
      else if (LogChoice == "Log y-axis") {
        scores[,3] <- log(scores[,3] + 1)
      }
    }
    else if (length(LogChoice) == 2) {
      scores[,2] <- log(scores[,2] + 1)
      scores[,3] <- log(scores[,3] + 1)
    }
    scores
    
  })
  
  FeatCompScatterPlot_react <- reactive({
    
    scores <- FeatCompScatter_react()
    geneset <- gs_react()
    geneset_name <- names(geneset)
    Feature <- input$ScatterFeature
    ColorCol <- input$ScatterColor
    LogChoice <- input$ScatterLog
    scoreMethod <- input$ScoreMethod
    
    if (geneset_name %in% decon_score_cols) {
      scoreMethod <- "Pre-Processed Score"
    }
    if (input$GeneSetTabs == 2) {
      scoreMethod <- "Gene Expression Score"
    }
    if (input$GeneSetTabs != 2 & !(geneset_name %in% decon_score_cols)) {
      scoreMethod <- paste(scoreMethod,"Score")
    }
    if (length(LogChoice) == 1) {
      if (LogChoice == "Log x-axis"){
        Feature <- paste(Feature,"(log(Feature + 1))")
      }
      else if (LogChoice == "Log y-axis") {
        scoreMethod <- paste(scoreMethod,"(log(score + 1))")
      }
    }
    else if (length(LogChoice) == 2) {
      Feature <- paste(Feature,"(log(Feature + 1))")
      scoreMethod <- paste(scoreMethod,"(log(score + 1))")
    }
    
    if (is.numeric(scores[,2])) {
      scores[,2] <- round(scores[,2],4)
    }
    if (is.numeric(scores[,3])) {
      scores[,3] <- round(scores[,3],4)
    }
    
    # plot
    p <- ggplot(scores, aes(x = scores[,2], y = scores[,3],
                            color = scores[,4],
                            text = paste("</br> Sample Name: ", scores[,1],
                                         "</br> ",colnames(scores)[2],": ",scores[,2],
                                         "</br> ",colnames(scores)[3],": ",scores[,3],
                                         sep =""))) +
      geom_point() +
      theme_minimal() +
      xlab(Feature)+
      ylab(scoreMethod) +
      labs(title = paste(Feature,"vs.\n",geneset_name,scoreMethod),
           color = colnames(scores)[4])
    if (input$ColorScatterChoice == "Single Color") {
      p <- p + theme(legend.position = "none") +
        geom_point(color = scores[1,4])
    }
    #p <- ggplotly(p,tooltip = "text")
    p
    
  })
  
  output$FeatCompScatterPlot <- renderPlotly({
    
    p <- FeatCompScatterPlot_react()
    ply <- ggplotly(p,tooltip = "text")
    ply
    
  })
  
  output$FeatCompScatterTable <- renderDataTable({
    
    tab <- FeatCompScatter_react()
    if (input$ColorScatterChoice == "Single Color") {
      tab <- tab[,-4]
    }
    DT::datatable(tab,
                  options = list(scrollY = T),
                  rownames = F)
    
  })
  
  
  
  ####----Text Output----####
  
  output$timewarnmessage1 <- renderUI({
    
    if (input$linPredict1 == "time") {
      p("Residual type must be schoenfeld or scaledsch.")
    }
    
  })
  
  output$timewarnmessage2 <- renderUI({
    
    if (input$linPredict2 == "time") {
      p("Residual type must be schoenfeld or scaledsch.")
    }
    
  })
  
  output$timewarnmessage3 <- renderUI({
    
    if (input$linPredict3 == "time") {
      p("Residual type must be schoenfeld or scaledsch.")
    }
    
  })
  
  
  ####----Downloaders----####
  
  ####----DNLD Path Surv----####
  
  ## quartile
  output$dnldSplot_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuartileSurvival.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuartileSurvival.svg",sep = "")
      }
    },
    content = function(file) {
      p <- Splot_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  output$dnldSplot_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuartileSurvival.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuartileSurvival.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- SplotBIN_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  
  ## binary
  output$dnldSplotBIN_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_MedianCutPSurvival.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_MedianCutPSurvival.svg",sep = "")
      }
    },
    content = function(file) {
      p <- SplotBIN_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  output$dnldSplotBIN_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_MedianCutPSurvival.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_MedianCutPSurvival.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- SplotBIN_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  
  ## cut p
  output$dnldScutPointPlot_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_OptimalCutpointSurvival.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_OptimalCutpointSurvival.svg",sep = "")
      }
    },
    content = function(file) {
      p <- ScutPointPlot_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  output$dnldScutPointPlot_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_OptimalCutpointSurvival.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_OptimalCutpointSurvival.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- ScutPointPlot_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  
  ## Quantile
  output$dnldSquantPlot_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuantileSurvival.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuantileSurvival.svg",sep = "")
      }
    },
    content = function(file) {
      p <- SquantPlot_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  output$dnldSquantPlot_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuantileSurvival.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuantileSurvival.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- SquantPlot_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  
  ## cutoff
  output$dnldSquantPlot2_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_AboveBelowCutoffSurvival.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_AboveBelowCutoffSurvival.svg",sep = "")
      }
    },
    content = function(file) {
      p <- SquantPlot2_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  output$dnldSquantPlot2_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_AboveBelowCutoffSurvival.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_AboveBelowCutoffSurvival.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- SquantPlot2_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  
  ####----DNLD PATH Density----####
  
  ## binary
  output$dnldssgseaBINDensity_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_MedianCutPDensity.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[metacol_sampletype])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_MedianCutPDensity.svg",sep = "")
      }
    },
    content = function(file) {
      p <- ssgseaBINDensity_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  output$dnldssgseaBINDensity_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_MedianCutPDensity.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_MedianCutPDensity.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- ssgseaBINDensity_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  
  ## cut p
  output$dnldssgseaCutPDensity_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_OptimalCutpointDensity.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_OptimalCutpointDensity.svg",sep = "")
      }
    },
    content = function(file) {
      p <- ssgseaCutPDensity_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  output$dnldssgseaCutPDensity_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_OptimalCutpointDensity.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_OptimalCutpointDensity.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- ssgseaCutPDensity_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  
  ## Quantile
  output$dnldssgseaQuantDensity_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuantileDensity.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuantileDensity.svg",sep = "")
      }
    },
    content = function(file) {
      p <- ssgseaQuantDensity_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  output$dnldssgseaQuantDensity_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuantileDensity.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuantileDensity.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- ssgseaQuantDensity_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  
  ## cutoff
  output$dnldssgseaQuant2Density_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_AboveBelowCutoffDensity.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_AboveBelowCutoffDensity.svg",sep = "")
      }
    },
    content = function(file) {
      p <- ssgseaQuant2Density_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  output$dnldssgseaQuant2Density_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_AboveBelowCutoffDensity.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_AboveBelowCutoffDensity.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- ssgseaQuant2Density_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  
  ####----DNLD Univar----####
  
  ##--Survival Plot--##
  
  output$dnldfeatSplot_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SingleSurvivalFeature
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateSurvival.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateSurvival.svg",sep = "")
      }
    },
    content = function(file) {
      p <- featSplot_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  output$dnldfeatSplot_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SingleSurvivalFeature
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateSurvival.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateSurvival.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- featSplot_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  
  ##--Forest Plot--##
  
  output$dnldUniVarForestplot_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SingleSurvivalFeature
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateForest.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateForest.svg",sep = "")
      }
    },
    content = function(file) {
      p <- SinglevarForestPlot_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  output$dnldUniVarForestplot_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SingleSurvivalFeature
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateForest.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateForest.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- SinglevarForestPlot_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  
  ##--Linearity Plot--##
  
  output$dnldUniVarLinplot_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SingleSurvivalFeature
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateLinearity.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateLinearity.svg",sep = "")
      }
    },
    content = function(file) {
      p <- UnivarLinearityPlot_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  output$dnldUniVarLinplot_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SingleSurvivalFeature
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateLinearity.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateLinearity.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- UnivarLinearityPlot_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  
  ####----DNLD Add Bivariate----####
  
  ##--Forest--##
  
  output$dnldBiVarAddForest_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SurvivalFeatureBi1
      Feature2 <- input$SurvivalFeatureBi2
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateForest.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateForest.svg",sep = "")
      }
    },
    content = function(file) {
      p <- BivarForestPlot_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  output$dnldBiVarAddForest_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SurvivalFeatureBi1
      Feature2 <- input$SurvivalFeatureBi2
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateForest.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateForest.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- BivarForestPlot_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  
  ##--Linearity Plot--##
  
  output$dnldBiVarAddLinplot_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SurvivalFeatureBi1
      Feature2 <- input$SurvivalFeatureBi2
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateLinearity.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateLinearity.svg",sep = "")
      }
    },
    content = function(file) {
      p <- BivarLinearityPlot_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  output$dnldBiVarAddLinplot_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SurvivalFeatureBi1
      Feature2 <- input$SurvivalFeatureBi2
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateLinearity.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateLinearity.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- BivarLinearityPlot_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  
  ####----DNLD Int Bivariate----####
  
  ##--Survival--##
  
  output$dnldfeatSplotBi_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SurvivalFeatureBi1Inter
      Feature2 <- input$SurvivalFeatureBi2Inter
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateSurvival.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateSurvival.svg",sep = "")
      }
    },
    content = function(file) {
      p <- featSplotBi_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  output$dnldfeatSplotBi_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SurvivalFeatureBi1Inter
      Feature2 <- input$SurvivalFeatureBi2Inter
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateSurvival.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateSurvival.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- featSplotBi_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  
  ##--Forest--##
  
  output$dnldBiVarIntForest_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SurvivalFeatureBi1Inter
      Feature2 <- input$SurvivalFeatureBi2Inter
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateForest.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateForest.svg",sep = "")
      }
    },
    content = function(file) {
      p <- BivarForestPlotInter_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  output$dnldBiVarIntForest_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SurvivalFeatureBi1Inter
      Feature2 <- input$SurvivalFeatureBi2Inter
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateForest.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateForest.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- BivarForestPlotInter_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  
  ##--Linearity--##
  
  output$dnldBiVarIntLinplot_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SurvivalFeatureBi1Inter
      Feature2 <- input$SurvivalFeatureBi2Inter
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateLinearity.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateLinearity.svg",sep = "")
      }
    },
    content = function(file) {
      p <- BivarLinearityPlotInter_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  output$dnldBiVarIntLinplot_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SurvivalFeatureBi1Inter
      Feature2 <- input$SurvivalFeatureBi2Inter
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateLinearity.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateLinearity.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- BivarLinearityPlotInter_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  
  ####----DNLD Multivariate----####
  
  ##--Forest--##
  
  output$dnldMultiVarForest_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SurvivalFeatureBi1Inter
      Feature2 <- input$SurvivalFeatureBi2Inter
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_MultivariateForest.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_MultivariateForest.svg",sep = "")
      }
    },
    content = function(file) {
      p <- MultivarForestPlot_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  output$dnldMultiVarForest_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SurvivalFeatureBi1Inter
      Feature2 <- input$SurvivalFeatureBi2Inter
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_MultivariateForest.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_MultivariateForest.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- MultivarForestPlot_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  
  ####----DNLD Risk Strat----####
  
  ##--Boxplot--##
  
  output$dnldSboxplot_SVG <- downloadHandler(
    filename = function() {
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      geneset <- gs_react()
      GeneSet <- names(geneset)
      scoreMethod <- input$ScoreMethod
      logchoice <- input$SBoxLog
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleTypeLab <- "_"
      }
      if (GeneSet %in% decon_score_cols) {
        scoreMethod <- "PreProcessed_Score"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethod <- "GeneExpressionScore"
      }
      if (input$GeneSetTabs != 2 & !(GeneSet %in% decon_score_cols)) {
        scoreMethod <- paste(scoreMethod,"Score",sep = "")
      }
      if (logchoice == TRUE) {
        scoreMethod <- paste(scoreMethod,"logPlus1",sep = "")
      }
      paste(gsub(" ","",ProjectName),SampleTypeLab,Feature,"_",GeneSet,"_",scoreMethod,"RiskStratBoxPlot.svg",sep = "")
    },
    content = function(file) {
      p <- Sboxplot_react()
      ggsave(file,p,width = 10, height = 8)
    }
  )
  output$dnldSboxplot_PDF <- downloadHandler(
    filename = function() {
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      geneset <- gs_react()
      GeneSet <- names(geneset)
      scoreMethod <- input$ScoreMethod
      logchoice <- input$SBoxLog
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleTypeLab <- "_"
      }
      
      if (GeneSet %in% decon_score_cols) {
        scoreMethod <- "PreProcessed_Score"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethod <- "GeneExpressionScore"
      }
      if (input$GeneSetTabs != 2 & !(GeneSet %in% decon_score_cols)) {
        scoreMethod <- paste(scoreMethod,"Score",sep = "")
      }
      if (logchoice == TRUE) {
        scoreMethod <- paste(scoreMethod,"logPlus1",sep = "")
      }
      paste(gsub(" ","",ProjectName),SampleTypeLab,Feature,"_",GeneSet,"_",scoreMethod,"RiskStratBoxPlot.pdf",sep = "")
    },
    content = function(file) {
      p <- Sboxplot_react()
      ggsave(file,p,width = 10, height = 8)
    }
  )
  
  output$dnldSBoxplotTab <- downloadHandler(
    filename = function() {
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      geneset <- gs_react()
      GeneSet <- names(geneset)
      scoreMethod <- input$ScoreMethod
      logchoice <- input$SBoxLog
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleTypeLab <- "_"
      }
      if (GeneSet %in% decon_score_cols) {
        scoreMethod <- "PreProcessed_Score"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethod <- "GeneExpressionScore"
      }
      if (input$GeneSetTabs != 2 & !(GeneSet %in% decon_score_cols)) {
        scoreMethod <- paste(scoreMethod,"Score",sep = "")
      }
      if (logchoice == TRUE) {
        scoreMethod <- paste(scoreMethod,"logPlus1",sep = "")
      }
      paste(gsub(" ","",ProjectName),SampleTypeLab,Feature,"_",GeneSet,"_",scoreMethod,".txt",sep = "")
    },
    content = function(file) {
      ssGSEA_meta <- SboxplotReact()
      geneset <- gs_react()
      GeneSet <- names(geneset)
      if (is.null(input$SurvivalType_time) == TRUE & is.null(input$SurvivalType_id) == TRUE) {
        surv_time_col <- metacol_survtime[1]
        surv_id_col <- metacol_survid[1]
      }
      if (is.null(input$SurvivalType_time) == FALSE & is.null(input$SurvivalType_id) == FALSE) {
        surv_time_col <- input$SurvivalType_time
        surv_id_col <- input$SurvivalType_id
      }
      boxTab <- as.data.frame(ssGSEA_meta[,c("SampleName",surv_time_col,surv_id_col,GeneSet,"SurvivalCutoff")])
      write_delim(boxTab,file,delim = '\t')
    }
  )
  
  ##--Heatmap--##
  
  output$dnldSheatmap_SVG <- downloadHandler(
    filename = function() {
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      geneset <- gs_react()
      GeneSet <- names(geneset)
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleTypeLab <- "_"
      }
      paste(gsub(" ","",ProjectName),SampleTypeLab,Feature,"_",GeneSet,"_","RiskStratHeatmap.svg",sep = "")
    },
    content = function(file) {
      p <- Sheatmap_react()
      ggsave(file,p,width = 10, height = 15)
    }
  )
  output$dnldSheatmap_PDF <- downloadHandler(
    filename = function() {
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      geneset <- gs_react()
      GeneSet <- names(geneset)
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleTypeLab <- "_"
      }
      paste(gsub(" ","",ProjectName),SampleTypeLab,Feature,"_",GeneSet,"_","RiskStratHeatmap.pdf",sep = "")
    },
    content = function(file) {
      p <- Sheatmap_react()
      ggsave(file,p,width = 10, height = 15)
    }
  )
  ## Download handler for expression
  output$dnldSheatmapexpr <- downloadHandler(
    filename = function() {
      # Variables
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      geneset <- gs_react()                 #Chosen Gene Set
      geneset_name <- names(geneset)        #Name of chosen gene set
      # Make file name
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"SurvivalCutoff_expr.txt",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"SurvivalCutoff_expr.txt",sep = "")
      }
    },
    content = function(file) {
      expr <- exprSub()
      geneset <- gs_react()                 #Chosen Gene Set
      geneset_name <- names(geneset)        #Name of chosen gene set
      GSgenes <- geneset[[geneset_name]]
      # Include only genes from gene set
      expr <- expr[which(rownames(expr) %in% GSgenes),]
      # Reformat to make sure genes show in file
      expr$Gene <- rownames(expr)
      expr <- expr %>%
        relocate(Gene)
      write_delim(expr,file,delim = '\t')
    }
  )
  
  ####----DNLD Feat Strat----####
  
  ##--Boxplot--##
  
  output$dnldFboxplot_SVG <- downloadHandler(
    filename = function() {
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      Feature2 <- input$BoxplotFeature
      geneset <- gs_react()
      GeneSet <- names(geneset)
      scoreMethod <- input$ScoreMethod
      logchoice <- input$FBoxLog
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleTypeLab <- "_"
      }
      if (GeneSet %in% decon_score_cols) {
        scoreMethod <- "PreProcessed_Score"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethod <- "GeneExpressionScore"
      }
      if (input$GeneSetTabs != 2 & !(GeneSet %in% decon_score_cols)) {
        scoreMethod <- paste(scoreMethod,"Score",sep = "")
      }
      if (logchoice == TRUE) {
        scoreMethod <- paste(scoreMethod,"logPlus1",sep = "")
      }
      paste(gsub(" ","",ProjectName),SampleTypeLab,Feature,"_",Feature2,"_",GeneSet,"_",scoreMethod,"FeatureStratBoxPlot.svg",sep = "")
    },
    content = function(file) {
      p <- Featureboxplot_react()
      ggsave(file,p,width = 10, height = 8)
    }
  )
  output$dnldFboxplot_PDF <- downloadHandler(
    filename = function() {
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      Feature2 <- input$BoxplotFeature
      geneset <- gs_react()
      GeneSet <- names(geneset)
      scoreMethod <- input$ScoreMethod
      logchoice <- input$FBoxLog
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleTypeLab <- "_"
      }
      
      if (GeneSet %in% decon_score_cols) {
        scoreMethod <- "PreProcessed_Score"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethod <- "GeneExpressionScore"
      }
      if (input$GeneSetTabs != 2 & !(GeneSet %in% decon_score_cols)) {
        scoreMethod <- paste(scoreMethod,"Score",sep = "")
      }
      if (logchoice == TRUE) {
        scoreMethod <- paste(scoreMethod,"logPlus1",sep = "")
      }
      paste(gsub(" ","",ProjectName),SampleTypeLab,Feature,"_",Feature2,"_",GeneSet,"_",scoreMethod,"FeatureStratBoxPlot.pdf",sep = "")
    },
    content = function(file) {
      p <- Featureboxplot_react()
      ggsave(file,p,width = 10, height = 8)
    }
  )
  
  output$dnldFeatureboxplotTab <- downloadHandler(
    filename = function() {
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      FeatureSelec <- input$BoxplotFeature
      geneset <- gs_react()
      GeneSet <- names(geneset)
      scoreMethod <- input$ScoreMethod
      logchoice <- input$FBoxLog
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleTypeLab <- "_"
      }
      if (GeneSet %in% decon_score_cols) {
        scoreMethod <- "PreProcessed_Score"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethod <- "GeneExpressionScore"
      }
      if (input$GeneSetTabs != 2 & !(GeneSet %in% decon_score_cols)) {
        scoreMethod <- paste(scoreMethod,"Score",sep = "")
      }
      if (logchoice == TRUE) {
        scoreMethod <- paste(scoreMethod,"logPlus1",sep = "")
      }
      
      paste(gsub(" ","",ProjectName),SampleTypeLab,Feature,"_",FeatureSelec,"_",GeneSet,"_",scoreMethod,".txt",sep = "")
    },
    content = function(file) {
      meta_ssGSEA <- ssGSEAmeta()
      FeatureSelec <- input$BoxplotFeature
      geneset <- gs_react()
      GeneSet <- names(geneset)
      boxTab <- meta_ssGSEA[,c("SampleName",FeatureSelec,GeneSet)]
      write_delim(boxTab,file,delim = '\t')
    }
  )
  
  ##--Heatmap--##
  
  output$dnldFheatmap_SVG <- downloadHandler(
    filename = function() {
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      geneset <- gs_react()
      GeneSet <- names(geneset)
      feature2 <- input$HeatmapFeature
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleTypeLab <- "_"
      }
      paste(gsub(" ","",ProjectName),SampleTypeLab,Feature,"_",feature2,"_",GeneSet,"_FeatureHeatmap.svg",sep = "")
    },
    content = function(file) {
      p <- FeatureHeatmap_react()
      ggsave(file,p,width = 10, height = 15)
    }
  )
  output$dnldFheatmap_PDF <- downloadHandler(
    filename = function() {
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      geneset <- gs_react()
      GeneSet <- names(geneset)
      feature2 <- input$HeatmapFeature
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleTypeLab <- "_"
      }
      paste(gsub(" ","",ProjectName),SampleTypeLab,Feature,"_",feature2,"_",GeneSet,"_FeatureHeatmap.pdf",sep = "")
    },
    content = function(file) {
      p <- FeatureHeatmap_react()
      ggsave(file,p,width = 10, height = 15)
    }
  )
  ## Download handler for expression
  output$dnldFheatmapexpr <- downloadHandler(
    filename = function() {
      # Variables
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      Feature2 <- input$HeatmapFeature
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      geneset <- gs_react()                 #Chosen Gene Set
      geneset_name <- names(geneset)        #Name of chosen gene set
      # Make file name
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"Featuring_",Feature2,"_expr.txt",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"Featuring_",Feature2,"_expr.txt",sep = "")
      }
    },
    content = function(file) {
      expr <- exprSub()
      geneset <- gs_react()                 #Chosen Gene Set
      geneset_name <- names(geneset)        #Name of chosen gene set
      GSgenes <- geneset[[geneset_name]]
      # Include only genes from gene set
      expr <- expr[which(rownames(expr) %in% GSgenes),]
      # Reformat to make sure genes show in file
      expr$Gene <- rownames(expr)
      expr <- expr %>%
        relocate(Gene)
      write_delim(expr,file,delim = '\t')
    }
  )
  
  ####----DNLD All Density----####
  
  output$dnldssgseaDensity_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_DensityPlot.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_DensityPlot.svg",sep = "")
      }
    },
    content = function(file) {
      p <- ssgseaDensity_react()
      ggsave(file,p,width = 10, height = 8)
    }
  )
  output$dnldssgseaDensity_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
        #if (input$RawOrSS == "Raw Gene Expression") {
        #  scoreMethodLab <- "RawGeneExpression"
        #}
        #else if (input$RawOrSS == "Rank Normalized") {
        #  scoreMethodLab <- scoreMethod
        #}
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_DensityPlot.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_DensityPlot.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- ssgseaDensity_react()
      ggsave(file,p,width = 10, height = 8)
    }
  )
  output$dnldssgseaDensityTable <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      score_method <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",score_method,".txt",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",score_method,".txt",sep = "")
      }
    },
    content = function(file) {
      geneset <- gs_react()
      GeneSet <- names(geneset)
      ssgsea_meta <- ssGSEAmeta()
      table <- ssgsea_meta[,c("SampleName",GeneSet)]
      write_delim(table,file,delim = '\t')
    }
  )
  
  ####---DNLD Meta----####
  
  ## Download handler for meta
  output$dnldMeta <- downloadHandler(
    filename = function() {
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      geneset <- gs_react()                 #Chosen Gene Set
      geneset_name <- names(geneset)        #Name of chosen gene set
      # Make file name
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_meta.txt",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_meta.txt",sep = "")
      }
    },
    content = function(file) {
      meta <- ssGSEAmeta()
      write_delim(meta,file,delim = '\t')
    }
  )
  
  ####---DNLD Expr----####
  
  ## Download handler for expression
  output$dnldExpr <- downloadHandler(
    filename = function() {
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      # Make file name
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"expr.txt",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"expr.txt",sep = "")
      }
    },
    content = function(file) {
      expr <- exprSub()
      # Reformat to make sure genes show in file
      expr$Gene <- rownames(expr)
      expr <- expr %>%
        relocate(Gene)
      write_delim(expr,file,delim = '\t')
    }
  )
  
  ####----DNLD Scatter Plot----####
  output$dnldFeatCompScatter_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      FeatureScatter <- input$ScatterFeature
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",FeatureScatter,"_",geneset_name,"_",scoreMethodLab,"_ScatterPlot.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",FeatureScatter,"_",geneset_name,"_",scoreMethodLab,"_ScatterPlot.svg",sep = "")
      }
    },
    content = function(file) {
      p <- FeatCompScatterPlot_react()
      ggsave(file,p,width = 10, height = 8)
    }
  )
  output$dnldFeatCompScatter_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      FeatureScatter <- input$ScatterFeature
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",FeatureScatter,"_",geneset_name,"_",scoreMethodLab,"_ScatterPlot.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",FeatureScatter,"_",geneset_name,"_",scoreMethodLab,"_ScatterPlot.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- FeatCompScatterPlot_react()
      ggsave(file,p,width = 10, height = 8)
    }
  )
  ## Download handler for expression
  output$dnldFeatCompScatterTable <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      FeatureScatter <- input$ScatterFeature
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",FeatureScatter,"_",geneset_name,"_",scoreMethodLab,"_ComparisonTable.txt",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",FeatureScatter,"_",geneset_name,"_",scoreMethodLab,"_ComparisonTable.txt",sep = "")
      }
    },
    content = function(file) {
      tab <- FeatCompScatter_react()
      if (input$ColorScatterChoice == "Single Color") {
        tab <- tab[,-4]
      }
      write_delim(tab,file,delim = '\t')
    }
  )
  
  
}





# Run the application
shinyApp(ui = ui, server = server)






