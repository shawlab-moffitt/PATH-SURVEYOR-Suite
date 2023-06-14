# Pathway-Connectivity

# Introduction

The integration of patient genome expression data, phenotypye data, and clinical data can serve as an integral resource for patient prognosis. PATH SURVEYOR: **PATH**way level **SURV**ival **E**nquir**Y** for Immuno-**O**ncology and Drug **R**epurposing serves to do just that, by examining the interaction of pathway analysis with patient expression and cilinical data to discover prominent features that take part in patient outcome. This utility is comprised of 3 R Shiny apps and a pipeline script which can be employed in a cohesive manor to provide an in-depth analysis towards pathway analysis of patient survival. Gene Set pathways utilized in this workflow include the Molecular Signatures Database (MSigDB), LINCS L1000 Small-Molecule Perturbations, and Clue.io ER Stress signtatures, as well as user provided gene sets. 

Here we focus on the Pathway Connectivity portion of this workflow with the DRPPM-Jaccard-Pathway-Connectivity R Shiny App. This app takes a list of gene sets as input and performs a Jaccard distance calculation to determine the proximity on the gene sets to one another. Working in tandem with the [PATH-SURVEYOR pipeline](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/3-PATH-SURVEYOR-Pipeline), the user may subset a top portion of the gene sets that were output from the comprehensive Cox Proportional Hazard table and use that as input to the Jaccard Connectivity app. This allows the user to gain another perspective on the top gene sets identified and how they cluster together by utilizing visualiztions of heatmaps, dendrograms, and phylogeny-type branched outputs.

An example Jaccard Connectivity R Shiny App is hosted here: https://shawlab-moffitt.shinyapps.io/pathway_connectivity/ where you are welcome to use the example inputs provided in the GitHub or your own to explore.

## The PATH-SURVEYOR Family

* [PATH-SURVEYOR File Prep App](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/1-Getting_Started/2-FilePrep)
* [PATH-SURVEYOR App [Interactive Mode]](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/2-PATH-SURVEYOR-Interactive-App)
* [PATH-SURVEYOR Cox Proportional Hazards Ranking [Pipeline Mode]](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/3-PATH-SURVEYOR-Pipeline)
* [Pathway Connectivity App](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/4-Pathway-Connectivity-App)
* [Pre-Ranked Hazard Ratio GSEA App](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/5-PreRanked-HazardRatio-GSEA-App)
* [PATH-SURVEYOR App with User File Input](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/6-PATH-SURVEYOR-UserInput-App)
* [PATH-SURVEYOR Docker Suite](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/7-PATH-SURVEYOR-Docker)

![alt text](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/blob/main/2-PATH-SURVEYOR-Interactive-App/App_Demo_Pictures/PATH_SURVEYOR_Main_schematic.PNG?raw=true)

# Installation

## Via Download

1. Download the [Zip File](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/archive/refs/heads/main.zip) from this GitHub repository: https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main
2. Unzip the downloaded file into the folder of your choice.

## Via Git Clone

1. Clone the [GitHub Repository](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite.git) into the destination of your choice.
   * Can be done in R Studio Terminal or a terminal of your choice
```bash
git clone https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite.git
```
3. Set your working directory in R to the cloned repository
   * This can be done through the "More" settings in the bottom-right box in R Stuido
   * You may also use the `setwd()` function in R Console.

# Requirments

* `R` - https://cran.r-project.org/src/base/R-4/
* `R Studio` - https://www.rstudio.com/products/rstudio/download/

# R Dependencies

|  |  |  |  |
| --- | --- | --- | --- |
| shiny_1.7.1 | shinythemes_1.2.0 | shinyjqui_0.4.1 | shinycssloaders_1.0.0 |
| DT_0.23 | pheatmap_1.0.12 | readr_2.1.2 | dplyr_1.0.9 |
| plotly_4.10.0 | clusterProfiler_4.0.5 | ggdendro_0.1.23 | factoextra_1.0.7 |
| reshape2_1.4.4 | stringr_1.4.0 | viridis_0.6.2 | RColorBrewer_1.1-3 |


# Required Files

* **Comprehensive Gene Set File:**
  * This is a provided file [GeneSet_Data/GeneSet_List_HS.RData](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/blob/main/4-Pathway-Connectivity-App/GeneSet_Data/GeneSet_List_HS.RData)
  * This is for use in the back end and provides the genes for the gene sets that are input.
  * The gene set names should match the ones provdide when running the PATH-SURVEYOR-Pipeline
    * If you ran the pipeline with a user provided gene set the genes for those gene sets will unlikely be found to compare distance between gene sets.

* **User Provided Pathways of Interest (.txt/.tsv):**
  * The only requirements for the file is that it it tab delimited and thr first column is the list of gene sets. 
    * The file will work if there is only one column or multiple, the app will only use the first column
  * This input should be a subset of the table that was output from the [PATH-SURVEYOR-Pipeline](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/3-PATH-SURVEYOR-Pipeline)
  * Example Input files are provided [here](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/4-Pathway-Connectivity-App/Example_Input_Files)

 * **Gene Annotation File (Optional):**
  * A tab delimited file with gene symbols as the first column followed by annotation columns 
  * It is recommended to use the output from the raw gene expression run of the [PATH-SURVEYOR pipeline](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Pipeline)
    * Though any file where the first column is gene symbols will do
    * An example input file was derived from a raw gene expression run of the [PATH-SURVEYOR-Pipeline](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/3-PATH-SURVEYOR-Pipeline) with the PAN ICI iAtlas Skin Cancer data to use for input can be found here: [Example_Input_Files/Gene_Annoation.txt](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/blob/main/4-Pathway-Connectivity-App/Example_Input_Files/Gene_Annoation.txt)

# App Set-Up

* It is important to ensure that the comprehensive gene set file that is provided is in the proper location for the app to locate it when running.
  * If the [Installation Section](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/4-Pathway-Connectivity-App#installation) is followed properly there should be no issue.
* To run the app:
  * The user can select the "Run App" button at the top write of the script in R Studio
  * Or the user can user the runApp() function in R Console
* When the app is running the user can select to input a file in the user interface and proceed with analysis

# App Features

## Sidebar Panel

###   Pathway and Clustering Parameters

![alt text](https://github.com/shawlab-moffitt/Pathway-Connectivity/blob/main/App_Pictures/Jaccard_Conn_Sidebar1.png?raw=true)

1. The user may pload their pathways of interest here
2. The user has the ability to choose which clustering method they want to use for the hclust() function, as well as the number of clusters they want to form with their data
3. A distance cutoff can be used input to generate a SIF file that can be downloaded (4) into a further application by the user. All gene set pair below the designated cutoff will be included in the file.

### Figure Parameters

![alt text](https://github.com/shawlab-moffitt/Pathway-Connectivity/blob/main/App_Pictures/Jaccard_Conn_Sidebar2.png?raw=true)

1. Figure parameter for theheatmap may be adjusted, such as color palette, column and row names, and dendrogram height.
2. The connectivity visualization can be customized to be viewed as a rectangular or circular denrogram or a phylogeny figure. 
   * If phylogeny is chosen a veriety of options to view the phylogeny figure as is provided.

## Main Panel

### Jaccard Pathway Connectivity Table

![alt text](https://github.com/shawlab-moffitt/Pathway-Connectivity/blob/main/App_Pictures/Jaccard_Conn_MainTable.png?raw=true)

1. When a file is uploaded to the app, after a few moments a Jaccard Connectivity able will appear showing the jaccard distance, 0-1, (similarity) between gene sets
   * The smaller the number, the more similar the gene set
2. The table can be downloaded for further use

### Connectivity Heatmap

![alt text](https://github.com/shawlab-moffitt/Pathway-Connectivity/blob/main/App_Pictures/Jaccard_Conn_heatmap.png?raw=true)

1. The heatmap give a global picture of similarit between gene sets

### Clustering

![alt text](https://github.com/shawlab-moffitt/Pathway-Connectivity/blob/main/App_Pictures/Jaccard_Conn_clustering.png?raw=true)

1. Clustering can be shown as a phylogenetic object, with or without names displyed. The names can be displayed also by hovering the points
   * The visualization is a plotly object, so the user may zoom in to interact with the plot
2. A dendrogram is another form of visualization available to see the clustering. This is also made with plotly, so the user may zoom in to examine the branches
3. The clusters can also be viewed as a circular dendrogram. This is not a plotly object and can not be interacted with.

### Clustering Annotation

![alt text](https://github.com/shawlab-moffitt/Pathway-Connectivity/blob/main/App_Pictures/Jaccard_Conn_annotation.png?raw=true)

* A data frame is displayed on the last tab starting with the gene set, cluster, and gene for each row. This allows users to see what genes are in the gene sets and cluster
1. A table provided by the user can be uploaded to annotate the genes. The uploaded table must list the gene symbol in the firsst column, the corresponding column can be any annotation the user chooses.
   * It is recommended to use the input from the PATH-SURVEYOR Pipeline raw gene expresison ranking output.

# Quesions and Comments

Please email Alyssa Obermayer at alyssa.obermayer@moffitt.org if you have any further comments or questions.


# Disclamer

Copyright 2022 Moffitt Cancer Center
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
