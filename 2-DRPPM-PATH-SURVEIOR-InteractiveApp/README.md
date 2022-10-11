# DRPPM - PATH SURVEIOR Shiny App

# Introduction

The integration of patient genome expression data, phenotype data, and clinical data can serve as an integral resource for patient prognosis. DRPPM PATH SURVEIOR: **Path**way level **Surv**ival **E**xam**i**nat**or** serves to do just that, by examining the interaction of pathway analysis with patient expression and cilinical data to discover prominent features that take part in patient outcome. This utility is comprised of 3 R Shiny apps and a pipeline script which can be employed in a cohesive manor to provide an in-depth analysis towards pathway analysis of patient survival. Gene Set pathways utilized in this workflow include the Molecular Signatures Database (MSigDB), LINCS L1000 Small-Molecule Perturbations, and Clue.io ER Stress signtatures, as well as user provided gene sets. 

Here we focus on the Interactive mode of this workflow with the DRPPM-PATH-SURVEIOR R Shiny App. With the expression, phenotype, and clincial data provided by the user we can integrate single sample GSEA (ssGSEA) pathway analysis with the comprehensive list of gene set pathways provided (or user provided). Additionally, upon app start-up, if the "immunedeconv" package is installed, immune deconvolution using ESTIMATE and MCP Counter methods will be performed. The score data is able to be partitioned and viewed in the survival plots or as a feature using median cut-point. The user can view a variety of survival plots based on binning the score data into quartile, quantile, above/below median, and optimal cut-point, or look through the lense of univariate, bivariate, and multivariate analysis with the integration of additional phenotype and clincal patient data. Further data exploration is available within the app to observe ssGSEA score density across the cohort as well as box plots and heatmaps to examine risk and feature stratification. The Shiny app comes complete with the ability to subset your cohort of patients, upload your own gene set data, along with customization and download of plots and tables throughout the app.

An example app using the PAN ICI iAtlas Checkpoint data can be see here: http://shawlab.science/shiny/DRPPM_PATH_SURVEIOR_PAN_ICI_iAtlas_Survival_App/ and is free to explore. This is the app that would be set up with the test data provided in this GitHub.

To facilitate identifying significant genes and pathways for further analysis, we have developed a Cox Proportional Hazard ranking script which ranks pathways or genes based on ssGSEA score or raw gene expression, respectively, above and below the median which returns a comprehensive table of pathways or genes ranked by Hazard Ratio which allows the user to find high-risk features with ease. When these are identified, the user can return to the interactive R Shiny App and visualize these features in real-time and perform additional bivariate or multivariate analyses to observe how the pathway survival interacts with covariates. More information on this pipeline can be found in our GitHub repository here: https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline.

## The DRPPM-PATH-SURVEIOR Family

* R Shiny Base Survival App [Interactive Mode]: https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR
* R Script for Cox Proportional Hazards Ranking [Pipeline Mode]: https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Pipeline
* R Shiny Jaccard Connectivity App: https://github.com/shawlab-moffitt/DRPPM-Jaccard-Pathway-Connectivity
* R Shiny Pre-Ranked GSEA App: https://github.com/shawlab-moffitt/DRPPM-PreRanked-GSEA

![alt text](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/App_Demo_Pictures/FlowChart_InteractiveMode.png?raw=true)

# Installation

## R Requirments

* `R` - https://cran.r-project.org/src/base/R-4/
   * R version >= 4.1 is recommended for use of the immune deconvolution R package
* `R Studio` - https://www.rstudio.com/products/rstudio/download/

## R Dependencies

|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| shiny_1.7.1 | shinythemes_1.2.0 | shinyjqui_0.4.1 | shinycssloaders_1.0.0 | dplyr_1.0.9 | tidyr_1.1.3 |
| readr_2.1.2 | gtsummary_1.6.0 | ggplot2_3.3.6 | pheatmap_1.0.12 | viridis_0.6.2 | plotly_4.10.0 |
| GSVA_1.40.1 | clusterProfiler_4.0.5 | ggpubr_0.4.0 | RColorBrewer_1.1-3 | tibble_3.1.7 | DT_0.23 |
| gridExtra_2.3 | survival_3.2-11 | survminer_0.4.9 | immunedeconv_2.1.0 |  |  |

* Above are the required packages for the R Shiny Application
* Users are encouraged to pre-intall the packages for quicker initial start-up of the application
* A package installation script is provided [R_Package_Installation.R](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/R_Package_Installation.R)
  * This is also in the app.R script in case the user does not pre-install

### Immune Deconvolution

* **Please note the immunedeconv package requires R version >= 4.1 to install**
* More information on this package can be found here: https://github.com/omnideconv/immunedeconv
  * This package requires additional packages to install completely, once all is installed you only need to load the immundeconv package for the app
* The user may use the app without the package installed and it will work, but without the immune deconvolution data

* The application processes MCP-counter and ESTIMATE deconvolution methods if the package is installed
* The immune deconvolution scores may be pre-process with the help of the R script we developed, this will supply you with the application update file inputs if performed according to instructions.
* More information found in the [Immune_Deconvolution folder](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/tree/main/Immune_Deconvolution) of this repository.

## Via Download

1. Download the [Zip File](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/archive/refs/heads/main.zip) from this GitHub repository: https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR
2. Unzip the downloaded file into the folder of your choice.
4. Set your working directory in R to the local version of the repository
   * This can be done through the "More" settings in the bottom-right box in R Stuido
   * You may also use the `setwd()` function in R Console.

## Via Git Clone

1. Clone the [GitHub Repository](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR.git) into the destination of your choice.
   * Can be done in R Studio Terminal or a terminal of your choice
```bash
git clone https://github.com/shawlab-moffitt/DRPPM-SURVIVE.git
```
3. Set your working directory in R to the cloned repository
   * This can be done through the "More" settings in the bottom-right box in R Stuido
   * You may also use the `setwd()` function in R Console.

# Required Files - User Provided

* **Expression Matrix (.tsv/.txt):**
  * Must be tab delimited with gene names as symbols located in the first column with subsequent columns consiting of the sample name as the header and expression data down the column.
  *  The App expects lowly expressed genes filtered out and normalized data either to FPKM or TMM.
     * Larger files might inflict memory issues for you local computer.
  * An example file with data from the PAN ICI iAtlas study on Skin and Kidney cancer tissue is loacted here [Pan_ICI_Example_Data/Pan_ICI_iAtlas_Skin_Kidney_Expression.zip](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/Pan_ICI_Example_Data/Pan_ICI_iAtlas_Skin_Kidney_Expression.zip)

* **Meta Data (.tsv/.txt):**
  * This should be a tab delimited file with each row depicting a sample by the same name as in the expression matrix followed by informative columns containing survival data and other features to analyze the samples by.
  * Required columns:
    * A sample name column
    * A survival time (in days) and ID column (can be more than one type of survival time)
    * Feature column(s) that allow for grouping of samples for analysis
  * Optional Column
    * A sample type column that allows for an initial subsetting of samples followed by grouping by feature (ex. Tissue or Disease Type)
    * Description column(s) that give additional information on the samples
  * An example file with data from the PAN ICI iAtlas study on Skin and Kidney cancer tissue is loacted here [Pan_ICI_Example_Data/Pan_ICI_iAtlas_Skin_Kidney_Meta.txt](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/Pan_ICI_Example_Data/Pan_ICI_iAtlas_Skin_Kidney_Meta.txt)
  * If you chose to perform the immune deconvolution pre-processing, the output meta data from that script should be used
    * Example here: [Pan_ICI_Example_Data/PAN_ICI_Skin_Kidney_Meta_ImmDeconv.txt](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/Pan_ICI_Example_Data/PAN_ICI_Skin_Kidney_Meta_ImmDeconv.txt)

* **Meta Data Parameters (.tsv/.txt):**
  * This should be a two-column tab-delimited file with the first column containing the column names of the meta file and the second column containing the column type of that meta column
  * The first column containing the users meta column names can be named however the user chooses, the **second column must have names matching the format provided below**.
    * For example, the survival ID column containing the 0/1 for events in the meta data could be nammed "OS" or "ID" or "OS.ID", but the column type in the second column meta parameter file must say "SurvivalID" for that specific column.
  * Column types and exmplinations:
    * **SampleName:** Contains sample names matching exprssion data (ONLY ONE ALLOWED)
    * **SampleType:** Contains a way to group and subset samples for further analysis (ONLY ONE ALLOWED and OPTIONAL)
    * **SurvivalTime:** Contains the overall survival time in days for the samples (can be other types of survival)
    * **SurvivalID:** Contains the survival ID for the samples, should be in a 0/1 format, 0 for alive/no event or 1 for dead/event (can be other types of survival)
    * **Feature:** Contains a feature that allows the samples to be grouped for analysis (More than one feature column allowed)
    * **Description:** Contains descriptions for the samples that may be viewed in the app (OPTIONAL)
  * Below is an example of these catagories and a full example file can be found here: [Pan_ICI_Example_Data/Pan_ICI_iAtlas_MetaData_Params.txt](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/Pan_ICI_Example_Data/Pan_ICI_iAtlas_MetaData_Params.txt)
  * If you chose to perform the immune deconvolution pre-processing, the output meta data from that script should be used
    * Example here: [Pan_ICI_Example_Data/PAN_ICI_Skin_Kidney_Meta_Params_ImmDeconv.txt](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/Pan_ICI_Example_Data/PAN_ICI_Skin_Kidney_Meta_Params_ImmDeconv.txt)

|  |  |
| --- | --- |
| sample | SampleName |
| Cancer_Tissue | SampleType |
| OS.time | SurvivalTime |
| OS.ID | SurivalID |
| EFS.time | SurvivalTime |
| EFS.ID | SurivalID |
| Responder | Feature |
| Race | Feature |
| Age | Description |
| Center | Description |

# Required Files - Provided

**Please keep in mind the paths to these documents, as they are relative to the path on the Github page. If the Github folder is your working directory the script should find these files**

* **Gene Set File:**
  * This is the file that contains the gene set names and genes for each gene set.
  * It is provided in RData list format
  * This file ([GeneSet_Data/GeneSet_List.RData](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/GeneSet_Data/GeneSet_List.RData)) contains gene sets from the following resources:
    * [The Molecular Signatures Database (MSigDB)](http://www.gsea-msigdb.org/gsea/msigdb/index.jsp)
    * [LINCS L1000 Small Molecule Perturbations](https://lincsproject.org/LINCS/)
    * [Cell Marker](https://academic.oup.com/nar/article/47/D1/D721/5115823)
    * [ER Stress Signatures from Clue.io](https://clue.io/)
    * Immunne Signatures
 
* **Gene Set Master Table:**
  * This is a three-column tab-delimited table the catagorizes and subcatagorizes the gene sets provided
  * It allows for organization of the large gene set list in the UI of the gene set selection for the Shiny App.
  * This file can be found here: [GeneSet_Data/GeneSet_Table.zip](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/GeneSet_Data/GeneSet_Table.zip)

# App Set-Up

* When using the DRPPM_SURVIVE App with the Pan ICI Checkpoint example data, you may follow the [Installation Section](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR#installation) of the README and may press the 'Run App' button in R studio or use the `runApp()` function in your console with the path to the directory of the app.R file.
* **When using your own files, please ensure all the required files above are in the proper directory**
  * The user must provide an expression matrix, meta data, and mata data parameter file.
* At the top of the app.R script you will find the user data input section where you input the path of your desired files (previewed below)
* When these files are entered the user may run the App by pressing the 'Run App' button in R studio or use the `runApp()` function in your terminal with the path to the app.R file
```{r}
####----User Input----####

ProjectName <- "Pan ICI Checkpoint Atlas"

ExpressionMatrix_file <- "Pan_ICI_Example_Data/Pan_ICI_iAtlas_Skin_Kidney_Expression.zip"

MetaData_file <- "Pan_ICI_Example_Data/Pan_ICI_iAtlas_Skin_Kidney_Meta.txt"

MetaParam_File <- "Pan_ICI_Example_Data/Pan_ICI_iAtlas_MetaData_Params.txt"
```


## Advanced Set-Up

### Pre-Set Selection Inputs

* There is a section at the top of the app.R script to allow users to pre-select their sample selection to be loaded on app start up.
* Users can select:
  * Sample Type choice with **PreSelect_SamplyType**
    * If there is a SampleType column in the meta, the user may select a variable from that column and write it in quotations
    * If NULL no choice will be pre-selected, the app will show the first option
    * If "all" (case ignored) the pre-selected option will be "All_Sample_Types"
  * Feature choice with **PreSelect_Feature**
    * Feature options come from feature columns in the meta data, the user would write the column name in quotations
    * If NULL no choice will be pre-selected, the app will show the first option
    * If "all" (case ignored) the pre-selected option will be "All_Features"
  * Sub-Feature Choice with **PreSelect_SubFeature**
    * This is only used when a Feature is pre-selected
    * SubFeature options come from unique values of the feature column that is pre-selected from the meta data, this would be writen in quotations
    * If NULL no choice will be pre-selected, the app will show the first option
  * Secondary Feature with **PreSelect_SecondaryFeature**
    * This is the feature that is used with univariate and multivariate analyses
    * Feature options come from feature columns in the meta data, the user would write the column name in quotations
    * If NULL no choice will be pre-selected, the app will show the first option
* Below is an example of this secion take from the app. It is located just below the user file input.

```{r}
## Pre-Selected Inputs
# An option from the meta, All, or NULL
PreSelect_SamplyType <- NULL
PreSelect_Feature <- "All"
# An option from the meta or NULL
PreSelect_SubFeature <- NULL
PreSelect_SecondaryFeature <- NULL
```

### Immune Deconvolution

Due to the variety of immune deconvolution methods available, we have developed a pre-processing script which allows the user to process through an array of different methods, add those features to the meta data and meta parameters, and view those results in correlation with survival data within the application.

More information on the pre-processing of this data can be found in the [Immune_Deconvolution folder](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/tree/main/Immune_Deconvolution) of this repository

# App Features

## Sidebar Panel

### Sample Selection and Parameters

![alt text](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/App_Demo_Pictures/SideBar_SampleParameters.png?raw=true)

1. Sample Type selection is an optional parameter that will appear if the user has a SampleType column to subset their data by. 
   * The user can select a single sample type to analyze or select all sample types
2. Feature selection allows the user to observe a specified feature from the feature columns annoated in the meta data parameter file.
   * The user has the option to select all features
3. Feature condition selection is related to Feature selection, where the Feature Condition options are updated based on the unique values of the Feature chosen.
   * If user selects all features, this eliminates the Feature Condition selection option
4. The user has a scoring method option based on the `gsva()` function perfomed
   * ssGSEA, GSVA, plage, or zscore
5. The gene set of interest is selected through the selection table.
   * The user may select a specific gene of interest or upload their own gene set file
6. The genes within the gene set that is chosen can be viewed by checking the box

### Survival Parameters

![alt text](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/App_Demo_Pictures/SideBar_SurvivalParameters.png?raw=true)

1. The user may select to view a specific type of survival analysis based on the available survival types in the meta data provided. 
   * For exmaple, OS, EFS, or PFS amoung others
2. One of the survival plots shown is a Quantile Survival Plot, this numeric input allows the user to choose their top and bottom quantile cutoff
3. The user may also specify the survival time and status cutoff when viewing the Risk Stratification Box Plot and Heatmap
   * The time is in days and the status signifies 0 for 'no-event' and 1 for an 'event'

### Figure Parameters

![alt text](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/App_Demo_Pictures/SideBar_FigureParamaters.png?raw=true)

1. The limit on years for the survival plots may be adjusted
2. The user may adjust the font and dot size, as well as the text orientation and stat comparison method for the boxplots within the app
   * The selections show "Wilcox.text" and "t.test" for 2 group boxplots and show "Wilcox.text", "t.test", "Kruskal.test" and "anova" for 3+ group boxplots
3. Heatmap clustering methods for the rows, font sizes, and color palettes
   * The Cluster options are "complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", and "centroid"
4. The font size for the forest plots my be adjusted
5. The font size for the linearity plots may be adjusted

## Main Panel

### Survival Analysis of Pathway Activity

![alt text](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/App_Demo_Pictures/MainPanel_SurvivalPlot.png?raw=true)

1. The Quartile Survival Plot shows at the top with a descriptive title indicating the feature, gene set, and score method
2. Each plot on the screen allows for the display of a hazard ratio table by selecting the checkbox
3. The Binary Survival Plot is displyed second with a title describing the feature, gene set, and score method
4. The Quantile Survival Plot is displayed third with a title describing the feature, gene set, score method, and user quantile input as indicated in the Survival Parameters side panel tab

### Univariate Survival Analysis

![alt text](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/App_Demo_Pictures/MainPanel_Univar_Survival.png?raw=true)

1. The user may view survival outcome based on a selected feature from the meta data
2. If the feature chosen is continuous, please check the box so the proper Cox Proportional Hazard analysis is performed
3. The reference feature for the Coxh analysis can be specified
4. The survival plot is show below, along with options to view the Coxh table, a forest plot, and a linearity check for continuous variables

### Bivariate Additive Survival Analysis

![alt text](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/App_Demo_Pictures/MainPanel_BivarAdd_Survival.png?raw=true)

1. Two features may be selected to view an additive Coxh survival analyis
2. If either feature is continuous the box should be checked
3. A reference variable for both features my be selected as well
4. The Cox Hazard Ratio table and summary is viewable, along with the Forest plot and a linearity check for continuous variables
5. An model comparison between the two features is performed through ANOVA

### Bivariate Interaction Survival Analysis

![alt text](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/App_Demo_Pictures/MainPanel_Bivar_Inter_Survival.png?raw=true)

1. Two features may be selected to view an additive Coxh survival analyis
2. If either feature is continuous the box should be checked
3. A reference variable for both features my be selected as well
4. The survival plot displaying the interaction of the two features is shown, along with the Cox HR table, forest plot, and a linearity plot

### Multivariate Survival Analysis

![alt text](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/App_Demo_Pictures/MainPanel_Multivar_Survival.png?raw=true)

1. The user may select multiple features to perform a Cox Proportion Hazard regression analysis on
   * If too many features are added the model may become convoluted
2. The Coxh tables and summaries are shown, along with options to view the forest and linearity plot
3. The Continous Coxh table is currently shown along side in case there is a mix of categorical and continuous features

### Meta Data Exploration

![alt text](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/App_Demo_Pictures/MainPanel_MetaTable.png?raw=true)

1. The user may select columns from the cumulative meta data to view in the UI table.
   * The table appears standard with the survival time, status, and current feature of interest
   * Additional columns are added to the selection options, such as the Quartile, Binary, Quantile, and ssGSEA calculations
2. The meta and expression data are available for download based on the subset critiria from the Sample Parameters

### ssGSEA Score Density

![alt text](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/App_Demo_Pictures/MainPanel_ssGSEA_Density.png?raw=true)

1. The user may define a percentile to view as a red vertical line on the ssGSEA score density plot
   * The plot originates with three blue dashed lines representing the 25, 50, and 75 percentile. This can be turned off through the check box below
2. The plot can be downloaded as an svg or pdf
3. The table below will show the samples along with their ssGSEA score for the gene set selected
4. This table may be downloaded

### Risk Stratification Box Plot

![alt text](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/App_Demo_Pictures/MainPanel_RiskStrat_BoxPlot.png?raw=true)

1. Survival parameters may be set in the side panel to bin the samples into two groups based on a user specified survival time and event status
2. A Survival Boxplot is generated based on the sample and feature parameters selected with the title indicating the gene set, score method, and feature
3. A table below displays the sample names being view along with their survival time, survival status, gene set score, and cutoff indicator column
4. This table can be downloaded for further anaylysis here

### Risk Stratification Heatmap

![alt text](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/App_Demo_Pictures/MainPanel_RiskStrat_Heatmap.png?raw=true)

1. Survival parameters may be set in the side panel to bin the samples into two groups based on a user specified survival time and event status
   * The Survival Heatmap is generated with the same Survival Parameters as the Survival Boxplot with the cutoff indication in the annotation at the top of the heatmap
   * The genes shown in the heatmap are the genes of the gene set selected
2. An expression matrix based on the genes and samples of the heatmap can be downloaded by the button at the bottom.
3. The heatmap size can be adjusted with the small triangle at the bottom-right
   * The row and column font size can also be adjusted in the Figure Parameters side panel tab

### Feature Box Plot

![alt text](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/App_Demo_Pictures/MainPanel_Feature_Boxplot.png?raw=true)

1. A feature may be selected to view amoung the samples that have been already subset by sample type, feature and feature condition
   * The boxplot title will indicate the gene set, score method, feature, and the additional feature being observed
2. A table will appear below listing the subset samples, selected additional feature, and the gene set score.
3. This table can be downloaded for further analysis

### Feature Heatmap

![alt text](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/App_Demo_Pictures/MainPanel_Feature_Heatmap.png?raw=true)

1. The Feature Heatmap is similar to the Feature Boxplot described above, where you may select an additional feature to view amoung you previously subset samples
2. The feature grouping is indicated and annotated at the top of the heatmap
   * The genes shown in the heatmap are the genes from the selected gene set
3. An expression matrix based on the genes and samples of the heatmap can be downloaded by the button at the bottom.
4. The heatmap size can be adjusted with the small triangle at the bottom-right
   * The row and column font size can also be adjusted in the Figure Parameters side panel tab

# Quesions and Comments

Please email Alyssa Obermayer at alyssa.obermayer@moffitt.org if you have any further comments or questions in regards to the R Shiny Application.
