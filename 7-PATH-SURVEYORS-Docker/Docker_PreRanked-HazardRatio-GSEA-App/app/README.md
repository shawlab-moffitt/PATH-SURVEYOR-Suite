# DRPPM-Pre-Rank-GSEA

# Introduction

The integration of patient genome expression data, phenotypye data, and clinical data can serve as an integral resource for patient prognosis. DRPPM PATH SURVEIOR: **Path**way level **Surv**ival **E**xam**i**nat**or** serves to do just that, by examining the interaction of pathway analysis with patient expression and cilinical data to discover prominent features that take part in patient outcome. This utility is comprised of 3 R Shiny apps and a pipeline script which can be employed in a cohesive manor to provide an in-depth analysis towards pathway analysis of patient survival. Gene Set pathways utilized in this workflow include the Molecular Signatures Database (MSigDB), LINCS L1000 Small-Molecule Perturbations, and Clue.io ER Stress signtatures, as well as user provided gene sets. 

Here we focus on Gene Set Enrichment Analysis (GSEA) with a pre-ranked list of genes through the use of the DRPPM-PreRanked-GSEA R Shiny App. This app utilizes the list of hazard ratio ranked genes that would be output from the [DRPPM-PATH-SURVEIOR pipeline](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline) when the user chooses to analyze raw gene expression, in general this app can accept any file of pre-ranked genes if desired. This allows the user to control what the method for ranking and alaysis to determine which gene sets may be significantly enriched from the expression data. The user can perform GSEA with many available gene sets within the app, such as MSigDB, LINCS L1000, Cell Marker, ER Stress, Immune Signatures, as well as upload their own gene set file to use. Furthermore, the user can visualize each gene set as an enrichment plot and the leading edge genes within that gene set.

An example Pre-Ranked GSEA R Shiny App is hosted here: http://shawlab.science/shiny/DRPPM_PATH_SURVEIOR_PreRanked_GSEA_App/ where you are welcome to use the example inputs provided in the GitHub or your own to explore.

## The DRPPM-PATH-SURVEIOR Family

* R Shiny Base Survival App [Interactive Mode]: https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR
* R Script for Cox Proportional Hazards Ranking [Pipeline Mode]: https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline
* R Shiny Jaccard Connectivity App: https://github.com/shawlab-moffitt/DRPPM-Jaccard-Pathway-Connectivity
* R Shiny Pre-Ranked GSEA App: https://github.com/shawlab-moffitt/DRPPM-PreRanked-GSEA

![alt text](https://github.com/shawlab-moffitt/DRPPM-PreRanked-GSEA/blob/main/App_Pictures/FlowChart_PreRankedGSEA.png?raw=true)

# Installation

## Via Download

1. Download the [Zip File](https://github.com/shawlab-moffitt/DRPPM-PreRanked-GSEA/archive/refs/heads/main.zip) from this GitHub repository: https://github.com/shawlab-moffitt/DRPPM-PreRanked-GSEA
2. Unzip the downloaded file into the folder of your choice.
4. Set your working directory in R to the local version of the repository
   * This can be done through the "More" settings in the bottom-right box in R Stuido
   * You may also use the `setwd()` function in R Console.

## Via Git Clone

1. Clone the [GitHub Repository](https://github.com/shawlab-moffitt/DRPPM-PreRanked-GSEA.git) into the destination of your choice.
   * Can be done in R Studio Terminal or a terminal of your choice
```bash
git clone https://github.com/shawlab-moffitt/DRPPM-PreRanked-GSEA.git
```
3. Set your working directory in R to the cloned repository
   * This can be done through the "More" settings in the bottom-right box in R Stuido
   * You may also use the `setwd()` function in R Console.

# Requirments

* `R` - https://cran.r-project.org/src/base/R-4/
* `R Studio` - https://www.rstudio.com/products/rstudio/download/

# R Dependencies

|  |  |  |  |  |
| --- | --- | --- | --- | --- |
| shiny_1.7.1 | shinythemes_1.2.0 | shinyjqui_0.4.1 | shinycssloaders_1.0.0 | enrichplot_1.12.3 |
| DT_0.23 | ggplot2_3.3.6 | readr_2.1.2 | dplyr_1.0.9 | clusterProfiler_4.0.5 |


# Required Files

* **Comprehensive Gene Set File:**
  * This is a provided file that contains various different gene sets to run the pre-ranked GSEA on
    * MSigDB, LINCS L1000, and Cell Marker
    * The app may also take user uploaded gene set files
      * Gene Set file must be in .gmt/.txt/.tsv format described [here](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline#required-files)

* **User Provided Pre Ranked Gene List:**
  * In the app, the user must upload a list of gene symbols, which can be formatted in different ways
  * Tab-delimited txt or tsv file
    * Preferred Format:
      * A two column file of gene symbols with their corresponding hazard ratio is preferred. This can be obtained from the output of the [DRPPM-PATH-SURVIOER-Pipeline](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline)
      * To extract this information, subset the gene symbol and hazard ratio column (the first two columns) to their own file
    * Other Formats
      * A single column of genes
        * The app will assume these are in ranked order and give them a ranking from 0-n in decreasing order
      * A table of more than two columns, the first of which being gene symbols
        * The app with take the first column and assume these are in ranked order and give them a ranking from 0-n in decreasing order

# App Set-Up

* It is important to ensure that the comprehensive gene set file that is provided is in the proper location for the app to locate it when running.
  * If the [Installation Section](https://github.com/shawlab-moffitt/DRPPM-PreRanked-GSEA#installation) is followed properly there should be no issue.
* To run the app:
  * The user can select the "Run App" button at the top write of the script in R Studio
  * Or the user can user the runApp() function in R Console
* When the app is running the user can select to input a file in the user interface and proceed with analysis

# App Features

## Enriched Signatures Table

![alt text](https://github.com/shawlab-moffitt/DRPPM-PreRanked-GSEA/blob/main/App_Pictures/PreRankGSEA_FirstTab.png?raw=true)

1. The user may upload their file here
   * Please check off if there is a header in your file
2. The GSEA function allows for an adjusted P.value cutoff to cutoff gene sets included in the enrichment table
3. The user may upload a ranked gene list of gene symbols corresponding to Human or Mouse
   * There are currently limited gene set databases available in the app for the mouse genome
   * The available gene set databases in the drop-down box will update depending on what specimen is chosen
   * When seleting a gene set, please note that depending on the size of the gene set, the enriched signatures table may take an extended period of time to generate
4. The user may upload their own gene set to perform GSEA
5. The user may preview all of the gene sets within the selected or uploaded gene set list
6. The Enriched signatures table will display when finished
7. When completed the table can be downloaded for furture use

## Enrichment Plots

![alt text](https://github.com/shawlab-moffitt/DRPPM-PreRanked-GSEA/blob/main/App_Pictures/PreRankGSEA_SecondTab.png?raw=true)

1. Based on the file uploaded and the enrichment table generated there gene sets that are enriched below the adjusted P.value cutoff will display to view as an enrichment plot
   * The will show with their corresponding Normalized Enrichment Score (NES) and P.Value, ranked by P.value (lowest to highest)
2. Above the enrichment plor the gene sets NES and P.value will be shown
3. The Enrichment plot will be displayed
4. The enrichment plot is available for download
5. The leading edge genes list can be downloaded
6. The leading edge gene list is also viewable, where you can see their rank within the gene set and their rank provided when the user file was input.

# Quesions and Comments

Please email Alyssa Obermayer at alyssa.obermayer@moffitt.org if you have any further comments or questions.
