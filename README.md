# PATH-SURVEYOR-Suite

# Introduction

The integration of patient genome expression data, phenotype data, and clinical data can serve as an integral resource for patient prognosis and treatment guidance. The PATH SURVEYOR Suite: **PATH**way level **SURV**ival **E**nquir**Y** for Immuno-**O**ncology and Drug **R**epurposing serves to do just that, by examining the interaction of gene or gene set pathway expression with clinical data to discover prominent features that play a role in patient outcome. This utility is comprised of R Shiny apps and a pipeline script which can be employed in a cohesive manor to provide an in-depth analysis towards underlying features affective survival. Through a systematic Cox proportional regression pipeline the expression of individual genes or gene set pathways is linked with a hazard ratio and p.values to determine signifigance and patient risk, which can be followed by interactive visualization in the suite of R Shiny applications.

Users can download and employ all of the tools locally to perform analyses and gain an in-depth perspective of their own data. To start, we provide example skin cancer data from the PAN ICI (Immune Checkpoint Inhibition) iAtlas database. Additionally, we provied a comprehensive set of pathways which havew been derieved from a variety of databases, such as MSigDB, LINCS L1000, Cell Marker, as well as derived gene sets focusing on ER stress and immune signatures. The example data and provided gene sets will be referenced and utilized throughout the guided overview of the PATH-SURVEYOR Suite of tools. 

## The PATH-SURVEYOR Family

* [PATH-SURVEYOR File Prep App](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/1-Getting_Started/2-FilePrep)
* [PATH-SURVEYOR App [Interactive Mode]](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/2-PATH-SURVEYOR-Interactive-App)
* [PATH-SURVEYOR Cox Proportional Hazards Ranking [Pipeline Mode]](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/3-PATH-SURVEYOR-Pipeline)
* [Pathway Connectivity App](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/4-Pathway-Connectivity-App)
* [Pre-Ranked Hazard Ratio GSEA App](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/5-PreRanked-HazardRatio-GSEA-App)
* [PATH-SURVEYOR App with User File Input](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/6-PATH-SURVEYOR-UserInput-App)
* [PATH-SURVEYOR Docker Suite](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/7-PATH-SURVEYOR-Docker)

![alt text](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/blob/main/2-PATH-SURVEYOR-Interactive-App/App_Demo_Pictures/PATH_SURVEYOR_Main_schematic.PNG?raw=true)

* Users can now explore various use cases of the PATH-SURVEYOR Suite of tools which are currently pending further review for publication, but can be seen here: https://github.com/shawlab-moffitt/PATH-SURVEYOR_Manuscript_Supplementary

# Installation

* Install PATH-SURVEYOR Suite GitHub repository
  * `git clone https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite.git`
  * Download and unzip repository https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/archive/refs/heads/main.zip
  * Install with Docker, instructions here: https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/7-PATH-SURVEYOR-Docker
* Install required R packages
  * Suite of tools was built on R version 4.1
  * R package installation script here: [1-Getting_Started/1-R_Package_Installation](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/1-Getting_Started/1-R_Package_Installation)

# Required Files

![alt text](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/blob/main/2-PATH-SURVEYOR-Interactive-App/App_Demo_Pictures/ExampleData.png?raw=true)

* Expression Matrix
  * Tab delimited matrix with **HGNC gene symbols** in the first column and sample names as the first-row header
  * Remove duplicate gene symbols is preferred, will check on app startup
  * Depending on the size, users may want to remove lowly expressed genes to reduce load time
* Clinical Meta Information
  * Tab-delimited file which should contain a first column of sample names matching the names in the expression matrix, followed by event and time-to-event information, as well as additional covariates or pre-processed scores
*	Clinical Feature Parameter
    * Tab-delimited two-column file where the first column consists of column names of the “Clinical Meta Information File” and the second column defining the column type
      * **SampleName** (mandatory): Contains sample names matching the expression data 
      * **SurvivalTime** (mandatory): Contains the overall survival time in days for the samples (can be other types of survival)
      * **SurvivalID** (mandatory): Contains the survival ID for the samples, should be in a 0/1 format, 0 for alive/no event or 1 for dead/event (can be other types of survival)
      * **SampleType** (optional): Higher level grouping of patient samples 
      * **Feature** (mandatory): One of more clinical or non-clinical features that can be included in the Cox-hazard analysis model.

* **NEW FEATURE** added to clean, properly format, and generate a Clinical Feature Parameter File for you here: https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/1-Getting_Started/2-FilePrep

# Immune Deconvolution Pre-Processing (Optional)

* Script
  * [1-Getting_Started/2-Immune_Deconvolution/Immune_Deconvolution.R](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/blob/main/1-Getting_Started/3-Immune_Deconvolution/Immune_Deconvolution.R)
  * Only available for R version 4.1 or greater
* Input
  * ProjectName: A descriptive name for your project/data
  * File inputs: Supply the path and file name or your expression matrix, clinical meta information, and clinical meta feature parameter files
  * Output_Path: Provide a path to write the output files to
  * Immune Deconvolution Methods: Indicate with TRUE or FALSE which methods to run
* Output
  * Updated clinical meta information and clinical meta feature parameter file which can be used as input to the interactive Shiny app


More Information Here: https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/1-Getting_Started/3-Immune_Deconvolution

# PATH-SURVEYOR: Interactive Mode

* Script
  * [2-PATH-SURVEYOR-Interactive-App/app.R](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/blob/main/2-PATH-SURVEYOR-Interactive-App/app.R)
* Input
  * Project Name: A descriptive name for your project/data
  * File inputs: Supply the path and file name or your expression matrix, clinical meta information, and clinical meta feature parameter files
  * Advanced User Input
    * Pre-set UI input options to be chosen upon startup, further described in GitHub README
  * Gene Set and Markdown files
    * These are provided files in the GitHub repository. Please ensure the path to these files are correct

More Information Here: https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/2-PATH-SURVEYOR-Interactive-App

# PATH-SURVEYOR: Pipeline Mode

* Script
  * Scripts to run in R Studio and in a command line interface found here: [3-PATH-SURVEYOR-Pipeline/](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/3-PATH-SURVEYOR-Pipeline)
* Input
  * Parameter File: Tab-delimited two-column file containing input file paths and run parameters described in the GitHub [3-PATH-SURVEYOR-Pipeline#parameter-file](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/3-PATH-SURVEYOR-Pipeline#parameter-file)
    * Pathway Level: When ranking gene set pathways according to Cox proportional hazards, a gene set file and name is required. Users can include a number of top pathways ranked on significance, and a Jaccard connectivity matrix will be included in the output
    * Gene Level: When ranking individual genes, no gene set file is required, though, if one is included, users can select to perform GSEA with a hazard ratio ranked list of genes upon Cox regression completion.
* Output
  * ssGSEA score table (gene set file required)
  * Median Cut-Point table
  * Cox Proportional Hazard Regression output for all pathways or genes

More Information Here: https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/3-PATH-SURVEYOR-Pipeline

# PATH-SURVEYOR: Pathway Connectivity

* Script
  * [4-Pathway-Connectivity-App/app.R](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/blob/main/4-Pathway-Connectivity-App/app.R)
* Input
  * Gene Set File: 
    * GeneSet_Data/Comprehensive_GeneSet.Rdata
  * User-Derived (Input upon app start-up)
    * CoxPH output file from PATH-SURVEYOR Pipeline pathway analysis
    * CoxPH output file from PATH-SURVEYOR Pipeline gene analysis 
* For gene cluster annotation (Optional) 
    * GMT file of gene set pathways is also accepted

More Information Here: https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/4-Pathway-Connectivity-App

# PATH-SURVEYOR: Pre-Ranked Hazard Ratio Ranked GSEA

* Script
    * [5-PreRanked-HazardRatio-GSEA-App/app.R](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/blob/main/5-PreRanked-HazardRatio-GSEA-App/app.R)
* Input
    * Gene Set File:
      * GeneSet_Data/GeneSets.zip
    * User-Derived (Input upon app start-up)
      * CoxPH output file from PATH-SURVEYOR Pipeline gene analysis

More Information Here: https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/5-PreRanked-HazardRatio-GSEA-App



# Disclamer

Copyright 2022 Moffitt Cancer Center
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
