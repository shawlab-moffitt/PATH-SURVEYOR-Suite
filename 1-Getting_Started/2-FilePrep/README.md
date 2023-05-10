# PATH-SURVEYOR FilePrep R Shiny App

To help ease the set-up process of the PATH-SURVEYOR suite of tools, we have developed two methods for the user to prepare their files to be compatible with the applications and various scripts. Users can cole the Github repo and run the application or script locally (both perform the same task), or the user can visit the following link to a server hosting the app: http://shawlab.science/shiny/Dev/PATH_SURVEYOR_FilePrep/

# Features

- [x] Remove troublesome special characters from sample names
- [x] Matching sample names between expression and meta data
- [x] Remove duplicate gene symbols from expression by summarizing to the gene with highest average expression
- [x] Filtering for genes with an average expression of zero
- [x] Convert survival time columns to days, if needed
- [x] Generate meta parameter file for use in the [PATH-SURVEYOR Shiny App](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/2-PATH-SURVEYOR-Interactive-App)

# User Setup

Users can run the app in R Studio, run the script in R, or run the [docker](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/7-PATH-SURVEYOR-Docker/Docker_PATH-SURVEYOR-FilePrep-App) image which will open a local R Shiny application

## R Requirments

* `R` - https://cran.r-project.org/src/base/R-4/
* `R Studio` - https://www.rstudio.com/products/rstudio/download/

## R Dependencies

|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| shiny_1.7.1 | stringr_1.5.0 | tools_4.2.2 | readr_2.1.2 | dplyr_1.0.9 | DT_0.23 |

# Required Files

* **Expression Matrix:**
  * Please have the first column as the gene names and the subsequent columns consiting of the sample name as the header and expression data down the column.
* **Meta Data:**
  * Please have the first column consist of the sample names that match the headers of the expression matrix.


# User Interface

![alt text](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/1-Getting_Started/2-FilePrep/Example_UI_Screenshots/PATH_SURVEYOR_FilePrepApp1.PNG?raw=true)

1. User may input their expression and meta data here and denmote how their file is delimited, as to be read in properly.
2. The app will try to detect survival time and event columns and will fill in any it thinks it finds. The user can delete any that are not true or search for ones that the app did not detect. 
   * Additionally, be sure to select the unit of time your survival time columns are in.
3. The app will run through the checks that are noted [above](https://github.com/shawlab-moffitt/PATH-SURVEYOR-FilePrep_ShinyApp#features) and will write out any findings at the top of the main panel.
4. Below that, the user can visualize the original expression data input compared to the new expression file they can download.
5. Lastly, the user can download the cleaned and properly formatted expression, meta, and meta parameter files for use in further analysis.

![alt text](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/1-Getting_Started/2-FilePrep/Example_UI_Screenshots/PATH_SURVEYOR_FilePrepApp2.PNG?raw=true)

* Along with cleaning and formatting the expression and meta data, the app with also generate a meta parameter file that the user can use as input to the [PATH-SURVEYOR Shiny App](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/2-PATH-SURVEYOR-Interactive-App)
* This file give each column in the meta data a small definition that allows the survival app to find the survival time and event columns with little effort. 
* As highlighted in the red, the file prep app with try to detect the survival time and event columns and generate a parameter file that defines these in the second column.


