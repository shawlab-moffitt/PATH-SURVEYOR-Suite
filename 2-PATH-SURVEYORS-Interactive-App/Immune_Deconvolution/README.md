# Guide to Immune Deconvolution Pre-Processing


# Introduction

Immune deconvolution allows the user the elucidate gene expression data and identify key cell populations identified across samples. The understanding of cell type composition allows for further clarification and additional infrances to be made on the data being examined which could be leveraged to find specific celltypes that may be used as therapeutic targets. When evaluating the cell type population scores, we can pair it with patient survival data to take a deeper look at how patients outcomes may be affected by the population levels of immune cell types. 

Many computational methodologies have been developed which produce a variety of cell type categories for the user to delineate their samples. To consolidate this effort of manually running through each method separately, we have developed a simple plug and run script which allows users to run through an array of prominent methods, including marker gene based methods (MCP-counter, Xcell, ESTIMATE) and formula based deconvolution methoeds (quanTIseq,EPIC,CIBERSORT,ABIS). The output of this script can be examined on its own as just the samples and their deconvolution score, but these scores can further be appended onto the users existing meta data to be used as input to our [PATH-SURVEYOR App](https://github.com/shawlab-moffitt/PATH-SURVEYOR). While the app can perform immune deconvolution using MCP-counter and ESTIMATE methods in the R version is up to date, but is helpful if users want to perform a deconvolution methods that may take an extended amount of time, or the machine they run the application on does not have R v4.1 or greater.

# R Requirments

* `R` - https://cran.r-project.org/src/base/R-4/
   * R version >= 4.1 is required for this script
* `R Studio` - https://www.rstudio.com/products/rstudio/download/

|  |  |  |  |
| --- | --- | --- | --- |
| dplyr_1.0.9 | readr_2.1.2 | stringr_1.4.0 | immunedeconv_2.1.0 |

## Immune Deconvolution Package

* Please go to this link for further information on package installation: https://github.com/omnideconv/immunedeconv
* Installation of this package could take up to 30 minutes depending on how many dependency packages need to be compiled
* The immune deconvolution methods CIBERSORT and CIBERSORT_abs require additional files which can be obtained from https://cibersortx.stanford.edu/
  * CIBERSORT.R and LM22.txt files are required
  * A license is required to download these files and may be obtained for free, upon approval
  * The script may still run without these files, but the CIBERSORT deconvolution method will not process
  
# Required Files

* **Expression Martix (.txt/.tsv):**
  * Must be tab delimited with gene names as symbols located in the first column with subsequent columns consiting of the sample name as the header and expression data down the column.
  * Example file here: [Pan_ICI_Example_Data/Pan_ICI_iAtlas_Skin_Kidney_Expression.zip](https://github.com/shawlab-moffitt/PATH-SURVEYOR/blob/main/Pan_ICI_Example_Data/Pan_ICI_iAtlas_Skin_Kidney_Expression.zip)

* **Meta Data (.txt/.tsv):**
  * This file is optional, but recommended if the user plans to view the immune deconvolution results in the [PATH-SURVEYOR App](https://github.com/shawlab-moffitt/PATH-SURVEYORS)
    * Including this file allows the script to output an updated meta file with the immune deconvolution result columns added
  * This should be a tab delimited file with each row depicting a sample by the same name as in the expression matrix followed by informative columns containing survival data and other features to analyze the samples by.
  * An example file here [Pan_ICI_Example_Data/Pan_ICI_iAtlas_Skin_Kidney_Meta.txt](https://github.com/shawlab-moffitt/PATH-SURVEYOR/blob/main/Pan_ICI_Example_Data/Pan_ICI_iAtlas_Skin_Kidney_Meta.txt)

* **Meta Data Parameters (.txt/.tsv):**
  * This file is optional, but recommended if the user plans to view the immune deconvolution results in the [PATH-SURVEYOR App](https://github.com/shawlab-moffitt/PATH-SURVEYOR)
    * Including this file allows the script to output an updated meta parameter file with additional immune deconvolution feature rows appended to the bottom
    * If this is not included, the script will output the feature rows and the user could append them to the parameter file manually
  * This should be a two-column tab-delimited file with the first column containing the column names of the meta file and the second column containing the column type of that meta column
    * Described in further detail [here](https://github.com/shawlab-moffitt/PATH-SURVEYOR#required-files---user-provided)
  * And example file here: [Pan_ICI_Example_Data/Pan_ICI_iAtlas_MetaData_Params.txt](https://github.com/shawlab-moffitt/PATH-SURVEYOR/blob/main/Pan_ICI_Example_Data/Pan_ICI_iAtlas_MetaData_Params.txt)

# Set-Up

* Input desired Project name, file names and paths, and an output path for the processed data to be held
```{r}
####----User Input----####
ProjectName <- "PAN_ICI_Skin_Kidney"
Expression_Matrix_File <- "Pan_ICI_Example_Data/Pan_ICI_iAtlas_Skin_Kidney_Expression.zip"
Meta_Data_File <- "Pan_ICI_Example_Data/Pan_ICI_iAtlas_Skin_Kidney_Meta.txt"
Meta_Data_Param_File <- "Pan_ICI_Example_Data/Pan_ICI_iAtlas_MetaData_Params.txt"
Output_Path <- "Pan_ICI_Example_Data/"
```
* Denote `TRUE` or `FALSE` for which methods you want to run through in the script
  * `mcp_counter` and `estimate` are ran in the [PATH-SURVEYOR App](https://github.com/shawlab-moffitt/PATH-SURVEYOR) upon start up if it is not included in the app sart-up files
```{r}
quantiseq <- TRUE
mcp_counter <- TRUE
xcell <- TRUE
epic <- TRUE
abis <- TRUE
estimate <- TRUE
```
* If the user would like to run the CIBERSORT methods please follow the instructions to obtain the proper files [in the installation step](https://github.com/shawlab-moffitt/PATH-SURVEYOR/blob/main/Immune_Deconvolution/README.md#immune-deconvolution-package)
* The script must be able to find the required files in the paths provided
  * If the CIBERSORT methods denote `TRUE` but the files are not found, the script will not run those methods
```{r}
cibersort <- TRUE
cibersort_abs <- TRUE
# CIBERSORT.R path and file name
CIBERSORT_Script <- "Path/To/CIBERSORT.R"
# LM22 path and file name
LM22_File <- "Path/To/LM22.txt"
```

# Running the Script

* When the set-up is complete the script can be run in sections if the user prefers, but should not be tampered with
* It is recommended to run the script as a local job in R Studio
  * On the bottom console of R Studio select the "Jobs" tab and select "Start a Local Job" then choose the edited script where you saved it and "PATH-SURVEYORS" as your working directory
    * The working directory may be adjusted of not touched if the user includes absolute paths to the files
* Please keep in mind, depending on the methods chosen and the size of the expression matrix the script may take several minutes to run
<p align="center">
  <img src="https://github.com/shawlab-moffitt/PATH-SURVEYOR-Pipeline/blob/main/Workflow_Picture/RStudio_LocalJob1.PNG?raw=true"/>
  <img src="https://github.com/shawlab-moffitt/PATH-SURVEYOR-Pipeline/blob/main/Workflow_Picture/RStudio_LocalJob2.PNG?raw=true"/>
</p>

# Script Output

Files will be output the the output path designated upon app set-up

* **Immune Deconvolution Score File:**
  * This will be a tab delimited file with the first column consiting of the sample names and the subsequent columns immune deconvolution scores name by cell type and method used to obtain the score

* **Updated Meta File:**
  * This will be output if a meta file is included during set-up
  * It will be composed of the input meta file with added columns of immune deconvolution scores and the scores partitions into high/low based on median cut-point
  * The added immune deconvolution columns are named in a way to be found by the [PATH-SURVEYOR App](https://github.com/shawlab-moffitt/PATH-SURVEYOR) and should not be tampered with
  * This updated file can be used with the updated meta parameter file as input to the [PATH-SURVEYOR App](https://github.com/shawlab-moffitt/PATH-SURVEYOR)

* **Updated Meta Parameter File:**
  * This will be output if a meta parameter file is included during set-up
  * If will be composed of the input meta parameter file with added rows of the column names of the immune deconvolution columns added to the meta and "Feature"
  * The added immune deconvolution columns in each row are named in a way to be found by the [PATH-SURVEYOR App](https://github.com/shawlab-moffitt/PATH-SURVEYOR) and should not be tampered with
  * This updated file can be used with the updated meta file as input to the [PATH-SURVEYOR App](https://github.com/shawlab-moffitt/PATH-SURVEYOR)


