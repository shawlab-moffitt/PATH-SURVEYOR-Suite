# DRPPM - PATH SURVEIOR - Pipeline

# Introduction

The integration of patient genome expression data, phenotypye data, and clinical data can serve as an integral resource for patient prognosis. DRPPM PATH SURVEIOR: **Path**way level **Surv**ival **E**xam**i**nat**or** serves to do just that, by examining the interaction of pathway analysis with patient expression and cilinical data to discover prominent features that take part in patient outcome. This utility is comprised of 3 R Shiny apps and a pipeline script which can be employed in a cohesive manor to provide an in-depth analysis towards pathway analysis of patient survival. Gene Set pathways utilized in this workflow include the Molecular Signatures Database (MSigDB), LINCS L1000 Small-Molecule Perturbations, and Clue.io ER Stress signtatures, as well as user provided gene sets. 

Here we focus on the Pipeline mode of this workflow with the DRPPM-PATH-SURVEIOR pipeline. With the expression, phenotype, and clinical data provided by the user along with a chosen set of pathways the user may set up a parameter file which, used as input to the pipeline script, will produce a comprhensive table ranking the input pathways or genes by hazard ratio. This script can facilitate two rankings, a ranking of pathways or of genes. For pathway ranking, the single sample GSEA score is determined per pathways and the samples are binned into groups of above and below the median pathway score for each gene set. This dicotomous variable (Above or Below median) is then run through a Cox proportional regression model with the survival data to produce a hazard ratio and set of P.values associated with High-Risk patients. This process may also be run on a gene level by binning samples based on raw gene expression as above or below median for each gene. The data is then run through the same Cox proportional hazards regression model as the pathways. The analysis may also be run with a chosen covariate from the patient data to observe how the interaction of the pathway or gene affects the hazard ratio model. The pipeline results in a comprehensive table of pathways or genes ranked by P.value and the associated hazard ratio, where a hazard ratio > 1 signifies High-Risk. 

The top pathways or genes that are idetntified through this pipeline can be displayed in real-time for validation through the use of the [DRPPM-PATH-SURVEIOR R Shiny App](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR). Depending on the analysis performed the table may be further analysed through additional apps of the DRPPM-SURVEIOR family. When the analysis is run with pathways the output table can subset into top hits and be used as input in the [DRPPM-Jaccard-Pathway-Connectivity R Shiny App](https://github.com/shawlab-moffitt/DRPPM-Jaccard-Pathway-Connectivity). This app allows the user to visualize the connectivity and pathway clusters of top gene sets of interest. When the analysis is run analyzing raw gene expression, the output table of the pipeline can be subset by taking the gene symbol column and hazard ratio column and used as input for the [DRPPM-PreRanked-GSEA R Shiny App](https://github.com/shawlab-moffitt/DRPPM-PreRanked-GSEA) which uses the hazard ratio to rank the genes and generate an enriched signature table and enrichment plots based on that ranking.

## The DRPPM-PATH-SURVEIOR Family

* R Shiny Base Survival App [Interactive Mode]: https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR
* R Script for Cox Proportional Hazards Ranking [Pipeline Mode]: https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline
* R Shiny Jaccard Connectivity App: https://github.com/shawlab-moffitt/DRPPM-Jaccard-Pathway-Connectivity
* R Shiny Pre-Ranked GSEA App: https://github.com/shawlab-moffitt/DRPPM-PreRanked-GSEA

![alt text](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/blob/main/Workflow_Picture/FlowChart_PipelineMode.png?raw=true)


# Installation

## Via Download

1. Download the [Zip File](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVIOER-Pipeline/archive/refs/heads/main.zip) from this GitHub repository: https://github.com/shawlab-moffitt/DRPPM-PATH-SURVIOER-Pipeline
2. Unzip the downloaded file into the folder of your choice.

## Via Git Clone

1. Clone the [GitHub Repository](https://github.com/shawlab-moffitt/DRPPM-SURVIVE.git) into the destination of your choice.
   * Can be done in R Studio Terminal or a terminal of your choice
```bash
git clone https://github.com/shawlab-moffitt/DRPPM-PATH-SURVIOER-Pipeline.git
```

# Requirments

* `R` - https://cran.r-project.org/src/base/R-4/
* `R Studio` - https://www.rstudio.com/products/rstudio/download/

# R Dependencies

|  |  |  |  |  |
| --- | --- | --- | --- | --- |
| DT_0.23 | tibble_3.1.7 | readr_2.1.2 | dplyr_1.0.9 | tidyr_1.1.3 |
| GSVA_1.40.1 | clusterProfiler_4.0.5 | gtsummary_1.6.0 | survival_3.2-11 | survminer_0.4.9 |

# Required Files

* **Expression Matrix (.tsv/.txt):**
  * Must be tab delimited with gene names as symbols located in the first column with subsequent columns consiting of the sample name as the header and expression data down the column.

* **Meta/Clincal Data (.tsv/.txt):**
  * This should be a tab delimited file with each row depicting a sample by the same name as in the expression matrix followed by informative columns containing survival data and other features to analyze the samples by.
  * Required columns:
    * A sample name column
    * A survival time and ID column (can be more than one type of survival time)

* **Gene Set File (.gmt/.txt/.tsv/.RData):**
  * This is the file that contains the gene set names and genes for each gene set.
  * If running the analysis on a single gene set, please include a gene set name in the propper location of the parameter file described below.
  * I provide example gene sets from publicly available sources here: https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/tree/main/GeneSets
    * These include the Molecular Signatures Database (MSigDB), LINCS L1000 small molecule perturbations, ER Stress Signatures, Immune Signatures, and Cell Markers.
  * The scripts accepts different formats of gene sets
    * .gmt format ([example](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/blob/main/GeneSets/MSigDB/MSigDB_Kegg_GeneSet.gmt)) described [here](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29)
    * .txt/.tsv (two-column tab-delimited) ([example](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/blob/main/GeneSets/MSigDB/MSigDB_Hallmark_GeneSets.txt)) with the first column being the gene set name repeating for every gene symbol that would be placed in the second column.
    * .RData list ([example readable with R](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/blob/main/GeneSets/MSigDB/MSigDB_HS_GeneSet.RData)) which is a named list of gene sets and genes. A script to generate this list is provided here: [GeneSetRDataListGen.R](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/blob/main/GeneSets/GeneSetRDataListGen.R)
  * If the user chooses to run the analysis on the genes only, no gene set file is needed
    
* **Gene Set .lst File (.lst) (OPTIONAL):**
  * In the case the user would like to run the pipeline with multiple gene set files, the user can provide a two column file containing the gene set name in the first column and the path to the gene set file in the second column.
  * [Example here](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/blob/main/GeneSets/Example_GeneSets.lst)
  * The gene set files listed should follow the format described above.
    
* **Parameter File (.txt/.tsv):**
  * This contains the file paths of the files described above, as well as run parameters.
  * This file is described in detail in the Set-Up section below.

# Set-Up

## Parameter File

* General Notes
  * This file will contain all the necessary information for the pipeline to run.
  * Some parameters are optional (Gene_Set_File,Gene_Set_Name,Covariate_Column_Label,Covariate_Reference).
  * Do not include a header in this file.
  * The order does not matter, but the parameter names do matter
  * Blank example parameter file [here](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/blob/main/Example_Run_Files/BLANK_ParamFile.txt)
  * Filled in example files [here](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/tree/main/Example_Run_Files)
* Project Name
  * Please be descriptive and add the survival data type in the name (OS/EFS/PFI/ect)
  * The script will add the gene set name, if given, to the output file
* Gene Set File and Name
  * This can be just a singular gene set file or a .lst file described [above](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline#required-files)
  * If a singular file please include a gene set name, though not required
* Output File Path
  * If no output file path is given, current working directory will be used
* Survival Data
  * The Survival Time and ID Label should be the column names of the respective survival data you want to analyze.
* Raw Gene Expression Analysis Option
  * The "Rank_Genes" parameter is an optional function (if answered TRUE) that will run the analysis based on each gene in the expression data.
  * This will produce a table of genes ranked by hazard ratio as well as the gene set(s) being run.
  * The user may choose to only rank genes and leave the Gene_Set_File parameter blank if they choose.
  * **This analysis may take an extended period of time**
* Covariate Analysis
  * The user may perform an interactive Coxh survival analysis between the pathway High/Low score and a specified categorical column from the meta data.
  * The user input should be the column name of the covariate.
  * If performing this analysis, it is requested the user provide a reference variable for the Coxh analysis, this would go in the "Covariate_Reference" parameter below.
    * If this is not specified the function will use the first level from the factored covariate column.
  * **This analysis may take an extended period of time** 

| Parameter | User Input |
| --- | --- |
| Project_Name | [User Defined Project Name] |
| Expression_Matrix_file | [Expression Matrix File] |
| Meta_Data_File | [Clinical Meta Information File] |
| Gene_Set_File | [Gene Set File] |
| Gene_Set_Name | [Name of Gene Set]
| Jaccard_Top_Hits | [Number of Top Pathways for Jaccard Connectivity] |
| Output_File_Path | [Path to Output Directory] |
| Survival_Time_Label | [Survival/Event Time Column Name] |
| Survival_ID_Label | [Survival/Event ID Column Name] |
| Rank_Genes | [TRUE or FALSE] |
| Covariate_Column_Label | [Covariate Column Name (Optional)] |
| Covariate_Reference | [Covariate Reference Feature (Optional)] |

## Running the script

**Please keep in mind, depending on how large the gene set files are and if the user chooses to perform the individual gene analysis, this script can take multiple hours to run.**

### R Studio

* The user can input the path to the parameter file at the top of the [RStudio version](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/blob/main/RStudio_Scripts/ssGSEA_Coxh_Ranking.R) of the script
  * The initial script is set up with the PAN ICI iAtlas Skin Cancer data with the MSigDB Hallmark gene set
* It is recommended to run the script as a local job in R Studio (See pictures below)
  * On the bottom console of R Studio select the "Jobs" tab and select "Start a Local Job" then choose the edited script where you saved it and "DRPPM-PATH-SURVEIOR-Pipeline" as your working directory 
  * On the top right of the script section, selcting the dropdown from "Source" will also allow you to start a local job, select the script file where is it located and "DRPPM-PATH-SURVEIOR-Pipeline" as your working directory 
* The user can select all of the contents and run the script or run section by section
  * Running the script this way will keep the R studio session busy and you will be able to run anything else while the script is running.

<p align="center">
  <img src="https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/blob/main/Workflow_Picture/RStudio_LocalJob1.PNG?raw=true"/>
  <img src="https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/blob/main/Workflow_Picture/RStudio_LocalJob2.PNG?raw=true"/>
</p>

### Command Line
* The user can run the [command line verison](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/blob/main/CommandLine_Scripts/ssGSEA_Coxh_Ranking.R) of the script in a command line environment as long as the requirments are met in the environment you are using.
* When the script and parameter file are in the desired directory run the command below:
```{linux}
Rscript ssGSEA_Coxh_Ranking.R [parameter_file]
```

## Pre-Set Example

* Multiple example cases with their own parameter files are set up in the [Example_Run_Files](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/tree/main/Example_Run_Files) folder.
* The may be run by following the [Intallation Steps](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline#installation) and copying the path to the selected parameter file to be read by the script
* The given example uses PAN ICI iAtlas checkpoint data from skin cancer patients. The expression and meta data is provided.
  * This data is the skin cancer subset from the PAN ICI iAtlas example data used in the [DRPPM-PATH-SURVEIOR R Shiny App](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR)
* [PAN_ICI_iAtlas_Skin_OS_MSigDBHallmark_Params](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/blob/main/Example_Run_Files/PAN_ICI_iAtlas_Skin_OS_MSigDBHallmark_Params.txt)
  * This parameter file runs the example data with only the 50 MSigDB Hallmark gene sets provided in the [GeneSets](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/tree/main/GeneSets) folder
    * Output example for this run is shown [here](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/tree/main/Example_Run_Files/PAN_ICI_iAtlas_OS_MSigDBHallmark_ExampleRun_Outputs)
  * Example provided with covariate analysis [here](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/blob/main/Example_Run_Files/PAN_ICI_iAtlas_Skin_OS_MSigDBHallmark_ResponderCovariate_Params.txt)
    * Output example for this run is shown [here](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/tree/main/Example_Run_Files/PAN_ICI_iAtlas_OS_MSigDBHallmark_ResponderCovariate_ExampleRun_Outputs)
* [PAN_ICI_iAtlas_Skin_OS_MSigDBHallmark_withGenes_Params](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/blob/main/Example_Run_Files/PAN_ICI_iAtlas_Skin_OS_MSigDBHallmark_withGenes_Params.txt)
  * This parameter file runs the 50 Hallmark Gene sets, as well as the raw gene expression analysis based on the genes in the expression data
  * Output example for this run is shown [here[(https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/tree/main/Example_Run_Files/PAN_ICI_iAtlas_OS_MSigDBHallmartk_withGenes_ExampleRun_Outputs)
    * Please note that that [Pan_ICI_iAtlas_Skin_OS_Genes_BIN.zip](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/blob/main/Example_Run_Files/PAN_ICI_iAtlas_OS_MSigDBHallmartk_withGenes_ExampleRun_Outputs/Pan_ICI_iAtlas_Skin_OS_Genes_BIN.zip) file was zipped after the run to upload to GitHub due to its 100+mb size
  * Example provided with covariate analysis [here](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/blob/main/Example_Run_Files/PAN_ICI_iAtlas_Skin_OS_MSigDBHallmark_withGenes_ResponderCovariate_Params.txt)
  * **This analysis will take an extended period of time**
    * For reference, running this locally (non-covariate analysis) on an Intel i7 10610U machine with 16GB of memory took 15+ hours, performing covariate analysis could take almost double this time 
* [PAN_ICI_iAtlas_Skin_OS_GeneSetLst_withGenes_Params](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/blob/main/Example_Run_Files/PAN_ICI_iAtlas_Skin_OS_GeneSetLst_withGenes_Params.txt)
  * This parameter file runs all of the gene sets provided as listed in the [gene set .lst file](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/blob/main/GeneSets/Example_GeneSets.lst) as well as the raw gene expression analysis
  * Example provided with covariate analysis [here](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline/blob/main/Example_Run_Files/PAN_ICI_iAtlas_Skin_OS_GeneSetLst_withGenes_ResponderCovariate_Params.txt)
  * **This analysis will take an extended period of time and is recommended to be ran on a command line prompt if available**

## Advanced Use Case



# Pipeline Steps and Output

## Prepping Data

* The script is set up to install any packages that are not found. The user may pre-load these is they choose
* The script begins by reading in the parameter file and files within the parameters file to set up the environment
* If survival data or covariate data is missing from any samples they will be removed, as they could skew the median and binning of samples by the ssGSEA scrore 

## ssGSEA Scoring

* If performing an analysis with gene sets, the ssGSEA score is determined for each gene sets within each gene set file provided
  * This file is output with the project name, gene set name, and "ssGSEA_score.txt" as the suffix for the file name
* This step is skipped if only analyzing raw gene expression

## Above and Below Median Calculation

* If raw gene expression is being analyzed the median raw gene expression for each gene is determined and the samples are binned into a group of "above" or "below" median for each gene
* For gene set analysis, each the median ssGSEA score for each gene set is determined and the samples are binned into a group of "above" or "below" median for each gene set
* These files are output as project name, gene set name (or "Genes"), and "BIN.txt" as the suffix for each file name

## Cox Propotional Hazard Analysis

* This step takes the longest to process
* Each binned gene set or gene is ran through the Cox proportional hazard function with the chosen survival data to generate the hazard ratio, concordance index, and p.values for each gene set
* If a covariate is included in the analysis, this variable will be taken into account in the CoxPH formula.
* As this formula runs it will continuously write out a line in the output file which is labeled with project name, gene set name (or "Genes"), and "coxh.txt" as the suffix for each file name
  * **Please DO NOT open this file while the script is running,** this will likely cause an error is R tries to write to a file that is currently open. 
    * If you want to check the process it is recommended that you copy the file to another location and open to view it there

## Ranking

* When the CoxPH formula is finished running though the gene sets or genes, the script with rank the table by the overall P.value
* The table will be complete with paramater lines at the top which are preceeded by "##"
* The first column will consist of the variable name which will be the gene set name or the gene depending on what was analyzed.
* The adjusted P.values ("BH" method) of Likelihood Ratio, Wald Test, and Logrank Test will be calculated and included


# Further Applications

The comprehensive table that is output from the analysis will elucidate some of the top pathways or genes that contribute to high-risk patients. If the [DRPPM-PATH-SURVIOER Shiny App](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR) was generated, these top pathways can be visualized within the app for validation. 

Additionally, to investigate the top pathways the user may take the top user specified number of pathways or a grouping of pathways the user wishes to explore further and input this list of pathways to the DRPPM-Jaccard-Pathway-Connectivity Shiny App and visualize the connectedness of the different pathways based on Jaccard Distance. The app set-up and function is explained in detail on the [GitHub page](https://github.com/shawlab-moffitt/DRPPM-Jaccard-Pathway-Connectivity).

Furthermore, the hazard ratio ranked gene list output from the script can be used as input for the Pre-Ranked GSEA Shiny App which allows users to perform Gene Set Enrichment Analysis with a pre-ranked set of genes on various different gene sets and provides enriched signature tables and enrichment plots for visualization. More information on set-up and funciton on out [GitHub page](https://github.com/shawlab-moffitt/DRPPM-Pre-Rank-GSEA).


# Quesions and Comments

Please email Alyssa Obermayer at alyssa.obermayer@moffitt.org if you have any further comments or questions.
