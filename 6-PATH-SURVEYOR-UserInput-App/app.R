###########################
##                       ##
##  DRPPM PATH SURVEYOR  ##
##      User Input       ##
##                       ##
###########################


##--Provided Input--##
## User make sure paths are correct
GeneSet_File <- "GeneSet_Data/GeneSet_List_HS.RData"
GeneSetTable_File <- "GeneSet_Data/GeneSet_CatTable.zip"
About_MD_File <- "App_Markdowns/PurposeAndMethods.Rmd"


ExampleExpr_File <- "Example_Data/Expression_Data.txt"
ExampleClin_File <- "Example_Data/Clinical_Data.txt"
ExampleParam_File <- "Example_Data/Clinical_Parameters.txt"




####----Install and load packages----####

## Check if Immune deconvolution package is installed
immudecon <- "immunedeconv"
immudecon_check <- immudecon %in% rownames(installed.packages())
if (immudecon_check == TRUE) {
  library(immunedeconv)
}
packages <- c("shiny","shinyjqui","gtsummary","tidyr","RColorBrewer","shinyjs",
              "dplyr","DT","ggplot2","ggpubr","tibble","survival","pheatmap","stringr",
              "readr","shinycssloaders","survminer","gridExtra","viridis","plotly","ggrepel")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))
#bioconductor packages
bioCpacks <- c("GSVA","clusterProfiler","AnnotationDbi")
installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
if (any(installed_packages_BIOC == FALSE)) {
   BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
}
invisible(lapply(bioCpacks, library, character.only = TRUE))



#increase file upload size
options(shiny.maxRequestSize=5000*1024^2)
##--Advanced Setup--##
SurvPlot_Height <- "850px"
SurvPlot_Width <- "550px"



##--Load in Gene Sets--##
# R Data list load function for naming
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
gs <- loadRData(GeneSet_File)
# Gene Set Table
GeneSetTable_og <- as.data.frame(read_delim(GeneSetTable_File, delim = '\t', col_names = T))
GeneSetTable_noCat <- GeneSetTable_og[,-1]
gsTab = TRUE
geneset_name_first <- GeneSetTable_og[1,4]


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

####----UI----####

ui <-
  navbarPage("{ DRPPM-PATH-SURVEYOR }",
             shinyjs::useShinyjs(),
             
             ####----Data Input----####
             
             tabPanel("Data Input",
                      fluidPage(
                        sidebarPanel(
                          width = 3,
                          id = "DataInputPanel",
                          p(),
                          #uiOutput("rendUserProjectName"),
                          textInput("UserProjectName","Project Name:", value = "Survival Analysis"),
                          uiOutput("rendExprFileInput"),
                          #fileInput("ExprFileInput","Expression Matrix",accept = c(".tsv",".txt",".csv",".zip")),
                          fluidRow(
                            column(5, style = 'padding-right:2px;margin-top:-30px;',
                                   #uiOutput("rendLogExprFile")
                                   checkboxInput("LogExprFile","Log2 Expression", value = F)
                            ),
                            column(7, style = 'padding-left:2px;margin-top:-30px;',
                                   #uiOutput("rendScaleNormExprFile")
                                   checkboxInput("ScaleNormExprFile","Scale Normalize Expression", value = F)
                            )
                          ),
                          uiOutput("rendClinFileInput"),
                          #fileInput("ClinFileInput","Clinical Data",accept = c(".tsv",".txt",".csv",".zip")),
                          fluidRow(
                            column(12, style = 'margin-top:-15px;',
                                   actionButton("UseExpData","Load Example Data"),
                                   tags$a(href="http://shawlab.science/shiny/PATH_SURVEYOR_ExampleData/PATH_SURVEYOR_App/", "Download example data", target='_blank'),
                                   )
                          ),
                          shiny::hr(),
                          uiOutput("rendClinParamHeader"),
                          h4("Clinical Parameters"),
                          #uiOutput("rendParamChoice"),
                          fluidRow(
                            column(12, style = 'margin-top:-20px;',
                                   radioButtons("ParamChoice","", choices = c("Define Parameters","Upload Parameter File"), inline = T)
                                   )
                            ),
                          fluidRow(
                            column(8,
                                   uiOutput("rendClinParamFileInput"),
                                   uiOutput("rendSurvTimeColSelect"),
                                   uiOutput("rendSurvIDColSelect")
                            ),
                            column(4,
                                   uiOutput("rendSurvTimeUnits")
                            )
                          ),
                        ),
                        mainPanel(
                          #uiOutput("rendAlertMessage"),
                          #h3("Expression File Preview"),
                          uiOutput("rendExprFilePrevHeader"),
                          div(DT::dataTableOutput("ExprFile_Preview"), style = "font-size:10px"),
                          #h3("Clinical File Preview"),
                          uiOutput("rendClinFilePrevHeader"),
                          div(DT::dataTableOutput("ClinFile_Preview"), style = "font-size:10px"),
                          #h3("Clinical Parameter File Preview"),
                          uiOutput("rendParamFilePrevHeader"),
                          div(DT::dataTableOutput("ClinParamFile_Preview"), style = "font-size:10px")
                        )
                      )
             ),
             
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
                                         column(4,
                                                uiOutput("rendSurvivalType_time")
                                         ),
                                         column(4,
                                                uiOutput("rendSurvivalType_id")
                                         ),
                                         column(4,
                                                uiOutput("rendScoreMethodBox")
                                         )
                                       ),
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
                                       shiny::hr(),
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
                                       shiny::hr(),
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
                                       shiny::hr(),
                                       h4("Forest Plot Parameters"),
                                       numericInput("ForestFontSize","Font Size",value = 1),
                                       shiny::hr(),
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
                                                  shiny::hr(),
                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("SplotBIN", width = SurvPlot_Height, height = SurvPlot_Width)), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldSplotBIN_SVG","Download as SVG"),
                                                    downloadButton("dnldSplotBIN_PDF","Download as PDF")
                                                  ),
                                                  shiny::hr(),
                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("ssgseaBINDensity", width = "650px", height = "400px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldssgseaBINDensity_SVG","Download as SVG"),
                                                    downloadButton("dnldssgseaBINDensity_PDF","Download as PDF")
                                                  ),
                                                  shiny::hr(),
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
                                                  shiny::hr(),
                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("Splot", width = SurvPlot_Height, height = SurvPlot_Width)), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldSplot_SVG","Download as SVG"),
                                                    downloadButton("dnldSplot_PDF","Download as PDF")
                                                  ),
                                                  shiny::hr(),
                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("ssgseaQuartDensity", width = "650px", height = "400px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldssgseaQuartDensity_SVG","Download as SVG"),
                                                    downloadButton("dnldssgseaQuartDensity_PDF","Download as PDF")
                                                  ),
                                                  shiny::hr(),
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
                                                  shiny::hr(),
                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("ScutPointPlot", width = SurvPlot_Height, height = SurvPlot_Width)), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldScutPointPlot_SVG","Download as SVG"),
                                                    downloadButton("dnldScutPointPlot_PDF","Download as PDF")
                                                  ),
                                                  shiny::hr(),
                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("ssgseaCutPDensity", width = "650px", height = "400px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldssgseaCutPDensity_SVG","Download as SVG"),
                                                    downloadButton("dnldssgseaCutPDensity_PDF","Download as PDF")
                                                  ),
                                                  shiny::hr(),
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
                                                  shiny::hr(),
                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("SquantPlot", width = SurvPlot_Height, height = SurvPlot_Width)), type = 6),
                                                  numericInput("QuantPercent","Top/Bottom Cut-Point Quantile Cutoff (%)", value = 25, min = 0, max = 100),
                                                  fluidRow(
                                                    downloadButton("dnldSquantPlot_SVG","Download as SVG"),
                                                    downloadButton("dnldSquantPlot_PDF","Download as PDF")
                                                  ),
                                                  shiny::hr(),
                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("ssgseaQuantDensity", width = "650px", height = "400px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldssgseaQuantDensity_SVG","Download as SVG"),
                                                    downloadButton("dnldssgseaQuantDensity_PDF","Download as PDF")
                                                  ),
                                                  shiny::hr(),
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
                                                  shiny::hr(),
                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("SquantPlot2", width = SurvPlot_Height, height = SurvPlot_Width)), type = 6),
                                                  numericInput("QuantPercent2","Above/Below User Quantile Cut-Point (%)", value = 25, min = 0, max = 100),
                                                  fluidRow(
                                                    downloadButton("dnldSquantPlot2_SVG","Download as SVG"),
                                                    downloadButton("dnldSquantPlot2_PDF","Download as PDF")
                                                  ),
                                                  shiny::hr(),
                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("ssgseaQuant2Density", width = "650px", height = "400px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldssgseaQuant2Density_SVG","Download as SVG"),
                                                    downloadButton("dnldssgseaQuant2Density_PDF","Download as PDF")
                                                  ),
                                                  shiny::hr(),
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
                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("featSplot", width = SurvPlot_Height, height = SurvPlot_Width)), type = 6),
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
                                                           div(shinycssloaders::withSpinner(tableOutput("SSingleFeatureHRtab"), type = 7, size = 0.5), style = "font-size:12px; width:500px; overflow-X: scroll")
                                                    ),
                                                    column(6,
                                                           verbatimTextOutput("UnivarSummary")
                                                    )
                                                  )
                                         ),
                                         
                                         ##--Forest Plot--##
                                         
                                         tabPanel("Forest Plot",
                                                  p(),
                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("SinglevarForestPlot", width = "100%", height = "800px")), type = 6),
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
                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("UnivarLinearityPlot", width = "100%", height = "500px")), type = 6),
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
                                                                      div(shinycssloaders::withSpinner(tableOutput("BiFeatureHRtab"), type = 7, size = 0.5), style = "font-size:12px; width:500px; overflow-X: scroll")
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
                                                             shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("BivarForestPlot", width = "100%", height = "800px")), type = 6),
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
                                                             shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("BivarLinearityPlot", width = "100%", height = "500px")), type = 6),
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
                                                             shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("featSplotBi", width = SurvPlot_Height, height = SurvPlot_Width)), type = 6),
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
                                                                      div(shinycssloaders::withSpinner(tableOutput("BiFeatureHRtabInter"), type = 7, size = 0.5), style = "font-size:12px; width:500px; overflow-X: scroll")
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
                                                             shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("BivarForestPlotInter", width = "100%", height = "800px")), type = 6),
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
                                                             shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("BivarLinearityPlotInter", width = "100%", height = "500px")), type = 6),
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
                                                                      div(shinycssloaders::withSpinner(tableOutput("SFeatureHRtabCat"), type = 7, size = 0.5), style = "font-size:12px; width:500px; overflow-X: scroll")
                                                               ),
                                                               column(6,
                                                                      h4("Coxh Hazard Ratio (Continuous)"),
                                                                      verbatimTextOutput("multivarSummaryCont"),
                                                                      div(shinycssloaders::withSpinner(tableOutput("SFeatureHRtabCont"), type = 7, size = 0.5), style = "font-size:12px; width:500px; overflow-X: scroll")
                                                               )
                                                             )
                                                    ),
                                                    
                                                    ##--Forest Plot--##
                                                    
                                                    tabPanel("Forest Plot",
                                                             p(),
                                                             shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("MultivarForestPlot", width = "100%", height = "800px")), type = 6),
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
                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("ssgseaDensity", width = "100%", height = "500px")), type = 6),
                                                  p(),
                                                  fluidRow(
                                                    downloadButton("dnldssgseaDensity_SVG","Download as SVG"),
                                                    downloadButton("dnldssgseaDensity_PDF","Download as PDF")
                                                  ),
                                                  p(),
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
                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotlyOutput("FeatCompScatterPlot", width = "100%", height = "500px")), type = 6),
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
                                                             shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("Sboxplot", width = "100%", height = "500px")), type = 6),
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
                                                             shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("Sheatmap", width = "100%", height = "2000px")), type = 6),
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
                                                             shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("Featureboxplot", width = "100%", height = "500px")), type = 6),
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
                                                             shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("FeatureHeatmap", width = "100%", height = "2000px")), type = 6),
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
             
             ####----About Panel----####
             
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




# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  observeEvent(input$ResetButton, {
    shinyjs::refresh()
  })
  
  
  ####----Data Input Panel----####
  
  ####----Render UI----####
  
  output$rendExprFileInput <- renderUI({
    
    refresh <- input$UseExpData
    fileInput("ExprFileInput","Expression Matrix",accept = c(".tsv",".txt",".csv",".zip"))
    
  })
  output$rendClinFileInput <- renderUI({
    
    refresh <- input$UseExpData
    fileInput("ClinFileInput","Clinical Data",accept = c(".tsv",".txt",".csv",".zip"))
    
  })
  
  output$rendExprFilePrevHeader <- renderUI({
    req((isTruthy(input$ExprFileInput) && isTruthy(input$ClinFileInput)) | (isTruthy(input$UseExpData) & (!isTruthy(input$ExprFileInput) && !isTruthy(input$ClinFileInput))))
    h3("Expression File Preview")
    
  })
  
  output$rendClinFilePrevHeader <- renderUI({
    req((isTruthy(input$ExprFileInput) && isTruthy(input$ClinFileInput)) | (isTruthy(input$UseExpData) & (!isTruthy(input$ExprFileInput) && !isTruthy(input$ClinFileInput))))
    h3("Clinical File Preview")
    
  })
  
  output$rendParamFilePrevHeader <- renderUI({
    req((isTruthy(input$ExprFileInput) && isTruthy(input$ClinFileInput)) | (isTruthy(input$UseExpData) & (!isTruthy(input$ExprFileInput) && !isTruthy(input$ClinFileInput))))
    h3("Clinical Parameter File Preview")
    
  })
  
  output$rendClinParamFileInput <- renderUI({
    
    if (!is.null(input$ParamChoice)) {
      if (input$ParamChoice == "Upload Parameter File") {
        refresh <- input$UseExpData
        fileInput("ClinParamFileInput","Clinical Parameter File", accept = c(".tsv",".txt",".csv",".zip"))
      }
    }
    
  })
  
  output$rendSurvTimeColSelect <- renderUI({
    req((isTruthy(input$ExprFileInput) && isTruthy(input$ClinFileInput)) | (isTruthy(input$UseExpData) & (!isTruthy(input$ExprFileInput) && !isTruthy(input$ClinFileInput))))
    
    if (!is.null(input$ParamChoice)) {
      if (input$ParamChoice == "Define Parameters") {
        clin <- clin_PostID_react()
        ColumnCheck <- grep(paste(c("OS_time","EFS_time","PFI_time","PFS_time","RFS_time","DSS_time","DFI_time"),collapse = "|"),
                            colnames(clin), ignore.case = T, value = T)
        ColCheckFound <- colnames(clin)[which(colnames(clin) %in% ColumnCheck)]
        if (length(ColCheckFound) == 0) {
          ColCheckFound <- ""
        }
        selectInput("SurvTimeColSelect","Survival Time Column(s)",
                    choices = colnames(clin), selected = ColCheckFound, multiple = T)
      }
    }
    
  })
  
  output$rendSurvTimeUnits <- renderUI({
    req((isTruthy(input$ExprFileInput) && isTruthy(input$ClinFileInput)) | (isTruthy(input$UseExpData) & (!isTruthy(input$ExprFileInput) && !isTruthy(input$ClinFileInput))))
    
    #req((isTruthy(input$ExprFileInput) && isTruthy(input$ClinFileInput)) | isTruthy(input$UseExpData))
    #req(input$ClinFileInput)
    #if (input$ParamChoice == "Define Parameters") {
    
    selectizeInput("SurvTimeUnits","Time Units:", choices = c("Days","Months","Years"))
    
    #}
    
  })
  
  output$rendSurvIDColSelect <- renderUI({
    req((isTruthy(input$ExprFileInput) && isTruthy(input$ClinFileInput)) | (isTruthy(input$UseExpData) & (!isTruthy(input$ExprFileInput) && !isTruthy(input$ClinFileInput))))
    
    #req((isTruthy(input$ExprFileInput) && isTruthy(input$ClinFileInput)) | isTruthy(input$UseExpData))
    if (!is.null(input$ParamChoice)) {
      if (input$ParamChoice == "Define Parameters") {
        clin <- clin_PostID_react()
        ColumnCheckTime <- grep(paste(c("OS_time","EFS_time","PFI_time","PFS_time","RFS_time","DSS_time","DFI_time"),collapse = "|"),
                                colnames(clin), ignore.case = T, value = T)
        ColumnCheck <- grep(paste(c("^OS","^EFS","^PFI","^PFS","^RFS","^DSS","^DFI",
                                    "OS_ID","EFS_ID","PFI_ID","PFS_ID","RFS_ID","DSS_ID","DFI_ID"),collapse = "|"),
                            colnames(clin), ignore.case = T, value = T)
        ColumnCheck <- ColumnCheck[which(!ColumnCheck %in% ColumnCheckTime)]
        ColCheckFound <- colnames(clin)[which(colnames(clin) %in% ColumnCheck)]
        if (length(ColCheckFound) == 0) {
          ColCheckFound <- ""
        }
        selectInput("SurvIDColSelect","Survival ID Column(s)",
                    choices = colnames(clin), selected = ColCheckFound, multiple = T)
      }
    }
    
  })
  
  ####----Reactives----####
  
  exprIn_nonT_react <- reactiveVal()
  clinIn_react <- reactiveVal()
  NAmetaCols_val <- reactiveVal()
  clinP_react <- reactiveVal()
  
  exprFileName <- reactiveVal()
  metaFileName <- reactiveVal()
  metaPFileName <- reactiveVal()
  
  observeEvent(input$UseExpData, {
    
    ##--Load example expression--##
    expr <- as.data.frame(readr::read_delim(ExampleExpr_File, delim = '\t', col_names = T))
    
    num_test <- apply(expr[,-1],2, is.numeric)
    if (all(num_test) == FALSE) {
      expr[,-1] <- mutate_all(expr[,-1], function(x) as.numeric(as.character(x)))
    }
    ## Replace any special characters to make uniform with expression
    colnames(expr) <- gsub(" ","_",colnames(expr))
    colnames(expr) <- gsub("[[:punct:]]","_",colnames(expr))
    colnames(expr)[1] <- "Symbol"
    if (TRUE %in% duplicated(expr[,1])) {
      expr <- expr %>%
        group_by(Symbol) %>%
        summarise_all(max)
    }
    expr <- as.data.frame(expr)
    rownames(expr) <- expr[,1]
    expr <- expr[,-1]
    exprIn_nonT_react(expr)
    
    ext <- tools::file_ext(ExampleExpr_File)
    ExampleExpr_File <- gsub(paste0(".",ext), "", basename(ExampleExpr_File))
    exprFileName(ExampleExpr_File)
    
    ##--Load example clinical--##
    clin <- as.data.frame(readr::read_delim(ExampleClin_File, delim = '\t', col_names = T))
    
    colnames(clin)[1] <- "SampleName"
    ## Replace any special characters to make uniform with expression
    clin[,1] <- gsub(" ","_",clin[,1])
    clin[,1] <- gsub("[[:punct:]]","_",clin[,1])
    
    colnames(clin) <- gsub(" ","_",colnames(clin))
    colnames(clin) <- gsub("[[:punct:]]","_",colnames(clin))
    
    clinIn_react(clin)
    
    ext <- tools::file_ext(ExampleClin_File)
    ExampleClin_File <- gsub(paste0(".",ext), "", basename(ExampleClin_File))
    metaFileName(ExampleClin_File)
    
    ##--Load example clinical parameters--##
    clinP_react(as.data.frame(readr::read_delim(ExampleParam_File, delim = '\t', col_names = T)))
    ext <- tools::file_ext(ExampleParam_File)
    ExampleParam_File <- gsub(paste0(".",ext), "", basename(ExampleParam_File))
    metaPFileName(ExampleParam_File)
  })
  
  #observeEvent(input$UseExpData, {
  #  clin <- as.data.frame(readr::read_delim(ExampleClin_File, delim = '\t', col_names = T))
  #  
  #  colnames(clin)[1] <- "SampleName"
  #  ## Replace any special characters to make uniform with expression
  #  clin[,1] <- gsub(" ","_",clin[,1])
  #  clin[,1] <- gsub("[[:punct:]]","_",clin[,1])
  #  
  #  colnames(clin) <- gsub(" ","_",colnames(clin))
  #  colnames(clin) <- gsub("[[:punct:]]","_",colnames(clin))
  #  
  #  clinIn_react(clin)
  #  
  #  ext <- tools::file_ext(ExampleClin_File)
  #  ExampleClin_File <- gsub(paste0(".",ext), "", basename(ExampleClin_File))
  #  metaFileName(ExampleClin_File)
  #})
  #
  #observeEvent(input$UseExpData, {
  #  clinP_react(as.data.frame(readr::read_delim(ExampleParam_File, delim = '\t', col_names = T)))
  #  ext <- tools::file_ext(ExampleParam_File)
  #  ExampleParam_File <- gsub(paste0(".",ext), "", basename(ExampleParam_File))
  #  metaPFileName(ExampleParam_File)
  #})
  
  ##--Expression--##
  
  ## Read in Expression Matrix - Pre user transform
  observe({
    
    ##--Load user expression--##
    gs.u <- input$ExprFileInput
    ext <- tools::file_ext(gs.u$datapath)
    req(gs.u)
    validate(need(ext == c("tsv","txt", "csv", "zip"), "Please upload .tsv, .txt, or .csv file"))
    
    if (ext == "csv") {
      matrix.u <- as.data.frame(read_delim(gs.u$datapath, delim = ',', col_names = T))
    }
    else {
      matrix.u <- as.data.frame(read_delim(gs.u$datapath, delim = '\t', col_names = T))
    }
    num_test <- apply(matrix.u[,-1],2, is.numeric)
    if (all(num_test) == FALSE) {
      matrix.u[,-1] <- mutate_all(matrix.u[,-1], function(x) as.numeric(as.character(x)))
    }
    ## Replace any special characters to make uniform with expression
    colnames(matrix.u) <- gsub(" ","_",colnames(matrix.u))
    colnames(matrix.u) <- gsub("[[:punct:]]","_",colnames(matrix.u))
    colnames(matrix.u)[1] <- "Symbol"
    if (TRUE %in% duplicated(matrix.u[,1])) {
      matrix.u <- matrix.u %>%
        group_by(Symbol) %>%
        summarise_all(max)
    }
    matrix.u <- as.data.frame(matrix.u)
    rownames(matrix.u) <- matrix.u[,1]
    matrix.u <- matrix.u[,-1]
    exprIn_nonT_react(matrix.u)
    
  })
  #exprIn_nonT_react <- reactive({
  #  
  #  gs.u <- input$ExprFileInput
  #  ext <- tools::file_ext(gs.u$datapath)
  #  req(gs.u)
  #  validate(need(ext == c("tsv","txt", "csv", "zip"), "Please upload .tsv, .txt, or .csv file"))
  #  
  #  if (ext == "csv") {
  #    matrix.u <- as.data.frame(read_delim(gs.u$datapath, delim = ',', col_names = T))
  #  }
  #  else {
  #    matrix.u <- as.data.frame(read_delim(gs.u$datapath, delim = '\t', col_names = T))
  #  }
  #  num_test <- apply(matrix.u[,-1],2, is.numeric)
  #  if (all(num_test) == FALSE) {
  #    matrix.u[,-1] <- mutate_all(matrix.u[,-1], function(x) as.numeric(as.character(x)))
  #    
  #  }
  #  ## Replace any special characters to make uniform with expression
  #  colnames(matrix.u) <- gsub(" ","_",colnames(matrix.u))
  #  colnames(matrix.u) <- gsub("[[:punct:]]","_",colnames(matrix.u))
  #  colnames(matrix.u)[1] <- "Symbol"
  #  if (TRUE %in% duplicated(matrix.u[,1])) {
  #    matrix.u <- matrix.u %>%
  #      group_by(Symbol) %>%
  #      summarise_all(max)
  #  }
  #  rownames(matrix.u) <- matrix.u[,1]
  #  matrix.u <- matrix.u[,-1]
  #  
  #  matrix.u
  #  
  #})
  
  GeneGS_table_react <- reactive({
    
    expr <- exprIn_nonT_react()
    exprGenes <- rownames(expr)
    GeneGS_table <- data.frame(Genes = exprGenes)
    GeneGS_table
    
  })
  
  ## Transform expression matrix if user choose
  exprIn_react <- reactive({
    
    if (!is.null(exprIn_nonT_react())) {
      matrix.u <- exprIn_nonT_react()
      matrix.u2 <- matrix.u
      if (input$LogExprFile == T) {
        matrix.u <- log2(matrix.u + 0.00001)
      }
      if (input$ScaleNormExprFile == T) {
        matrix.u = apply(matrix.u, 1, scale)
        matrix.u = apply(matrix.u, 1, rev)
        colnames(matrix.u) <- colnames(matrix.u2)
      }
      matrix.u <- as.matrix(matrix.u[sort(rownames(matrix.u)),])
      matrix.u
    }
    
  })
  
  ##--Clinical--##
  
  ## User Meta Upload
  observe({
    
    gs.u <- input$ClinFileInput
    ext <- tools::file_ext(gs.u$datapath)
    req(gs.u)
    validate(need(ext == c("tsv","txt","csv", "zip"), "Please upload .tsv, .txt, or .csv file"))
    
    if (ext == "csv") {
      meta.u <- as.data.frame(read_delim(gs.u$datapath, delim = ',', col_names = T))
    }
    else {
      meta.u <- as.data.frame(read_delim(gs.u$datapath, delim = '\t', col_names = T))
    }
    NAmetaCols <- colnames(meta.u[,colSums(is.na(meta.u))==nrow(meta.u)])
    NAmetaCols_val(NAmetaCols)
    meta.u <- meta.u[,colSums(is.na(meta.u))<nrow(meta.u)]
    colnames(meta.u)[1] <- "SampleName"
    ## Replace any special characters to make uniform with expression
    meta.u[,1] <- gsub(" ","_",meta.u[,1])
    meta.u[,1] <- gsub("[[:punct:]]","_",meta.u[,1])
    
    colnames(meta.u) <- gsub(" ","_",colnames(meta.u))
    colnames(meta.u) <- gsub("[[:punct:]]","_",colnames(meta.u))
    
    
    clinIn_react(meta.u)
    
  })
  #clinIn_react <- reactive({
  #  
  #  gs.u <- input$ClinFileInput
  #  ext <- tools::file_ext(gs.u$datapath)
  #  req(gs.u)
  #  validate(need(ext == c("tsv","txt","csv", "zip"), "Please upload .tsv, .txt, or .csv file"))
  #  
  #  if (ext == "csv") {
  #    meta.u <- as.data.frame(read_delim(gs.u$datapath, delim = ',', col_names = T))
  #  }
  #  else {
  #    meta.u <- as.data.frame(read_delim(gs.u$datapath, delim = '\t', col_names = T))
  #  }
  #  colnames(meta.u)[1] <- "SampleName"
  #  ## Replace any special characters to make uniform with expression
  #  meta.u[,1] <- gsub(" ","_",meta.u[,1])
  #  meta.u[,1] <- gsub("[[:punct:]]","_",meta.u[,1])
  #  
  #  colnames(meta.u) <- gsub(" ","_",colnames(meta.u))
  #  colnames(meta.u) <- gsub("[[:punct:]]","_",colnames(meta.u))
  #  
  #  
  #  meta.u
  #  
  #})
  
  ##--Compare Sample Names--##
  
  ## Make sure sample names are same in both expression and clinical
  expr_react <- reactive({
    
    req((isTruthy(input$ExprFileInput) && isTruthy(input$ClinFileInput)) | (isTruthy(input$UseExpData) & (!isTruthy(input$ExprFileInput) && !isTruthy(input$ClinFileInput))))
    #if (!is.null(exprIn_nonT_react()) && !is.null(clinIn_react())) {
      expr <- exprIn_react()
      clin <- clinIn_react()
      
      sampsame <- intersect(clin[,"SampleName"],colnames(expr))
      expr <- expr[,which(colnames(expr) %in% sampsame)]
      expr
    #}
    
  })
  
  ## Make sure sample names are same in both expression and clinical
  clin_PreID_react <- reactive({
    #req((isTruthy(input$ExprFileInput) && isTruthy(input$ClinFileInput)) | isTruthy(input$UseExpData))
    req((isTruthy(input$ExprFileInput) && isTruthy(input$ClinFileInput)) | (isTruthy(input$UseExpData) & (!isTruthy(input$ExprFileInput) && !isTruthy(input$ClinFileInput))))
    #if (!is.null(exprIn_nonT_react()) && !is.null(clinIn_react())) {
    #if (!is.null(exprIn_nonT_react()) && !is.null(clinIn_react())) {
      expr <- exprIn_react()
      clin <- clinIn_react()
      
      sampsame <- intersect(clin[,"SampleName"],colnames(expr))
      clin <- clin[which(clin[,"SampleName"] %in% sampsame),]
      
      clin
    #}
    
  })
  
  clin_PostID_react <- reactive({
    #req((isTruthy(input$ExprFileInput) && isTruthy(input$ClinFileInput)) | isTruthy(input$UseExpData))
    req((isTruthy(input$ExprFileInput) && isTruthy(input$ClinFileInput)) | (isTruthy(input$UseExpData) & (!isTruthy(input$ExprFileInput) && !isTruthy(input$ClinFileInput))))
    #if (!is.null(exprIn_nonT_react()) && !is.null(clinIn_react())) {
      expr <- expr_react()
      clin <- clin_PreID_react()
      
      if (nrow(clin) > 0 & ncol(expr) > 0) {
        
        PreProcessed_meta_cols <- c(grep("_PreProcessedScore$",colnames(clin),value = T))
        if (length(PreProcessed_meta_cols) == 0) {
          if (immudecon_check == TRUE) {
            mcp_counter_decon <- as.data.frame(deconvolute(expr, "mcp_counter"))
            rownames(mcp_counter_decon) <- mcp_counter_decon[,1]
            mcp_counter_decon <- mcp_counter_decon[,-1]
            mcp_counter_decon <- as.data.frame(t(mcp_counter_decon))
            colnames(mcp_counter_decon) <- paste(gsub(" ","_",colnames(mcp_counter_decon)),"mcp_counter_PreProcessedScore",sep = "_")
            mcp_counter_decon$SampleName <- rownames(mcp_counter_decon)
            clin <- merge(clin,mcp_counter_decon)
            
            estimate_decon <- as.data.frame(deconvolute(expr, "estimate"))
            rownames(estimate_decon) <- estimate_decon[,1]
            estimate_decon <- estimate_decon[,-1]
            estimate_decon <- as.data.frame(t(estimate_decon))
            colnames(estimate_decon) <- paste(gsub(" ","_",colnames(estimate_decon)),"estimate_PreProcessedScore",sep = "_")
            estimate_decon$SampleName <- rownames(estimate_decon)
            clin <- merge(clin,estimate_decon)
          }
        }
        clin
        
      }
      
    #}
    
  })
  
  clin_react <- reactive({
    
    #req((isTruthy(input$ExprFileInput) && isTruthy(input$ClinFileInput)) | isTruthy(input$UseExpData))
    req((isTruthy(input$ExprFileInput) && isTruthy(input$ClinFileInput)) | (isTruthy(input$UseExpData) & (!isTruthy(input$ExprFileInput) && !isTruthy(input$ClinFileInput))))
    #if (!is.null(exprIn_nonT_react()) && !is.null(clinIn_react())) {
      
      clin <- clin_PostID_react()
      if (!is.null(clin)) {
        
        clinp <- clinP_react()
        colnames(clin) <- gsub("_PreProcessedScore","",colnames(clin))
        TimeUnits <- input$SurvTimeUnits
        TimCols <- clinp[which(clinp[,2] == "SurvivalTime"),1]
        #metacol_survtime <- input$SurvTimeColSelect
        
        if (!is.null(TimeUnits)) {
          if (TimeUnits == "Months") {
            for (i in TimCols) {
              clin[,i] <- clin[,i] * 30.4375
            }
          }
          if (TimeUnits == "Years") {
            for (i in TimCols) {
              clin[,i] <- clin[,i] * 365.25
            }
          }
        }
        clin
      }
      
      
      #clin
    #}
    
  })
  
  ##--Clinical Parameters--##
  
  ## User Clinical Parameter Upload
  observe({
    
    req((isTruthy(input$ExprFileInput) && isTruthy(input$ClinFileInput)) | (isTruthy(input$UseExpData) & (!isTruthy(input$ExprFileInput) && !isTruthy(input$ClinFileInput))))
    if (!is.null(input$ParamChoice)) {
      if (input$ParamChoice == "Upload Parameter File") {
        
        gs.u <- input$ClinParamFileInput
        ext <- tools::file_ext(gs.u$datapath)
        req(gs.u)
        validate(need(ext == c("tsv","txt","csv", "zip"), "Please upload .tsv, .txt, or .csv file"))
        
        if (ext == "csv") {
          MetaParam <- as.data.frame(read_delim(gs.u$datapath, delim = ',', col_names = F))
        }
        else {
          MetaParam <- as.data.frame(read_delim(gs.u$datapath, delim = '\t', col_names = F))
        }
        
        NAmetaCols <- NAmetaCols_val()
        MetaParam <- MetaParam[which(!MetaParam[,1] %in% NAmetaCols),]
        MetaParam <- as.data.frame(apply(MetaParam,2,trimws))
        MetaParam <- MetaParam[which(!MetaParam[,1] %in% NAmetaCols),]
        MetaParam[,1] <- gsub(" ","_",MetaParam[,1])
        MetaParam[,1] <- gsub("[[:punct:]]","_",MetaParam[,1])
        MetaParam[which(MetaParam[,2] == "Description"),2] <- "Feature"
        MetaParam[which(MetaParam[,2] == "SampleName"),1] <- "SampleName"
        
        colnames(MetaParam) <- c("Clinical_Column_Name","Clinical_Column_Type")
        
        clin <- clin_PostID_react()
        if (nrow(clin) > 0) {
          
          missingFeat <- colnames(clin)[which(!colnames(clin) %in% MetaParam[,1])]
          if (length(missingFeat) > 0) {
            MetaParam2 <- data.frame(Clinical_Column_Name = missingFeat,
                                     Clinical_Column_Type = "Feature")
            MetaParam <- rbind(MetaParam,MetaParam2)
          }
          clinP_react(MetaParam)
          
        }
        
      }
      else if (input$ParamChoice == "Define Parameters") {
        
        clin <- clin_PostID_react()
        if (!is.null(clin)) {
          
          if (nrow(clin) > 0) {
            
            clin <- clin[,colSums(is.na(clin))<nrow(clin)]
            metacol_survid <- input$SurvIDColSelect
            metacol_survtime <- input$SurvTimeColSelect
            if (length(metacol_survtime) > 0 && length(metacol_survid) > 0) {
              metacol_feature <- colnames(clin)[-1]
              metacol_feature <- metacol_feature[!metacol_feature %in% c(metacol_survtime,metacol_survid)]
            }
            
            if (length(metacol_survtime) > 0 && length(metacol_survid) > 0) {
              MetaParam1 <- data.frame(Clinical_Column_Name = "SampleName",
                                       Clinical_Column_Type = "SampleName")
              MetaParam2 <- data.frame(Clinical_Column_Name = metacol_survtime,
                                       Clinical_Column_Type = "SurvivalTime")
              MetaParam3 <- data.frame(Clinical_Column_Name = metacol_survid,
                                       Clinical_Column_Type = "SurvivalID")
              MetaParam4 <- data.frame(Clinical_Column_Name = metacol_feature,
                                       Clinical_Column_Type = "Feature")
              MetaParam <- rbind(MetaParam1,MetaParam2,MetaParam3,MetaParam4)
              clinP_react(MetaParam)
            }
            
          }
          
        }
        
        
      }
    }
    
  })
  #clinP_react <- reactive({
  #  
  #  if (!is.null(input$ParamChoice)) {
  #    if (input$ParamChoice == "Upload Parameter File") {
  #      
  #      gs.u <- input$ClinParamFileInput
  #      ext <- tools::file_ext(gs.u$datapath)
  #      req(gs.u)
  #      validate(need(ext == c("tsv","txt","csv", "zip"), "Please upload .tsv, .txt, or .csv file"))
  #      
  #      if (ext == "csv") {
  #        MetaParam <- as.data.frame(read_delim(gs.u$datapath, delim = ',', col_names = F))
  #      }
  #      else {
  #        MetaParam <- as.data.frame(read_delim(gs.u$datapath, delim = '\t', col_names = F))
  #      }
  #      
  #      MetaParam[,2] <- gsub(" ","",MetaParam[,2])
  #      MetaParam[,1] <- gsub(" ","_",MetaParam[,1])
  #      MetaParam[,1] <- gsub("[[:punct:]]","_",MetaParam[,1])
  #      MetaParam[which(MetaParam[,2] == "Description"),2] <- "Feature"
  #      MetaParam[which(MetaParam[,2] == "SampleName"),1] <- "SampleName"
  #      
  #      colnames(MetaParam) <- c("Clinical_Column_Name","Clinical_Column_Type")
  #      
  #      clin <- clin_PostID_react()
  #      missingFeat <- colnames(clin)[which(!colnames(clin) %in% MetaParam[,1])]
  #      if (length(missingFeat) > 0) {
  #        MetaParam2 <- data.frame(Clinical_Column_Name = missingFeat,
  #                                 Clinical_Column_Type = "Feature")
  #        MetaParam <- rbind(MetaParam,MetaParam2)
  #      }
  #      MetaParam
  #    }
  #    else if (input$ParamChoice == "Define Parameters") {
  #      
  #      clin <- clin_PostID_react()
  #      metacol_survid <- input$SurvIDColSelect
  #      metacol_survtime <- input$SurvTimeColSelect
  #      if (length(metacol_survtime) > 0 && length(metacol_survid) > 0) {
  #        metacol_feature <- colnames(clin)[-1]
  #        metacol_feature <- metacol_feature[!metacol_feature %in% c(metacol_survtime,metacol_survid)]
  #      }
  #      
  #      if (length(metacol_survtime) > 0 && length(metacol_survid) > 0) {
  #        MetaParam1 <- data.frame(Clinical_Column_Name = "SampleName",
  #                                 Clinical_Column_Type = "SampleName")
  #        MetaParam2 <- data.frame(Clinical_Column_Name = metacol_survtime,
  #                                 Clinical_Column_Type = "SurvivalTime")
  #        MetaParam3 <- data.frame(Clinical_Column_Name = metacol_survid,
  #                                 Clinical_Column_Type = "SurvivalID")
  #        MetaParam4 <- data.frame(Clinical_Column_Name = metacol_feature,
  #                                 Clinical_Column_Type = "Feature")
  #        MetaParam <- rbind(MetaParam1,MetaParam2,MetaParam3,MetaParam4)
  #        MetaParam
  #      }
  #    }
  #  }
  #  
  #})
  
  # "SampleType" category depreciated
  metacol_sampletype <- reactive({
    
    metacol_sampletype <- NULL
    if (!is.null(input$ParamChoice)) {
      if (input$ParamChoice == "Upload Parameter File") {
        MetaParam <- clinP_react()
        if (("SampleType" %in% MetaParam[,2]) == TRUE) {
          metacol_sampletype <- MetaParam[which(MetaParam[,2] == "SampleType"),1]
        }
        if (("SampleType" %in% MetaParam[,2]) == FALSE) {
          metacol_sampletype <- NULL
        }
      }
    }
    metacol_sampletype
    
  })
  
  metacol_survtime <- reactive({
    
    if (!is.null(input$ParamChoice)) {
      if (input$ParamChoice == "Upload Parameter File") {
        MetaParam <- clinP_react()
        metacol_survtime <- MetaParam[which(MetaParam[,2] == "SurvivalTime"),1]
        metacol_survtime
      }
      else if (input$ParamChoice == "Define Parameters") {
        metacol_survtime <- input$SurvTimeColSelect
        metacol_survtime
      }
    }
    
  })
  
  metacol_survid <- reactive({
    
    if (!is.null(input$ParamChoice)) {
      if (input$ParamChoice == "Upload Parameter File") {
        MetaParam <- clinP_react()
        metacol_survid <- MetaParam[which(MetaParam[,2] == "SurvivalID"),1]
        metacol_survid
      }
      else if (input$ParamChoice == "Define Parameters") {
        metacol_survid <- input$SurvIDColSelect
        metacol_survid
      }
    }
    
  })
  
  ## feature with pre process disignation
  metacol_feature_og <- reactive({
    
    MetaParam <- clinP_react()
    metacol_feature <- MetaParam$Clinical_Column_Name
    metacol_feature
    
  })
  
  decon_score_cols <- reactive({
    
    metacol_feature <- metacol_feature_og()
    PreProcessed_meta_cols <- c(grep("_PreProcessedScore$",metacol_feature,value = T))
    metacol_feature_noPreProcessedScore <- c(grep("_PreProcessedScore$",metacol_feature,value = T, invert = T))
    decon_score_cols <- c()
    
    if (length(PreProcessed_meta_cols) > 0) {
      
      immDeconvCols <- c()
      mcp_devonv_scores <- grep("mcp_counter_PreProcessedScore$",metacol_feature,value = T)
      if (length(mcp_devonv_scores) > 0) {
        immDeconvCols <- c(immDeconvCols,mcp_devonv_scores)
        decon_score_cols <- c(decon_score_cols,mcp_devonv_scores)
      }
      estimate_scores <- grep("estimate_PreProcessedScore$",metacol_feature,value = T)
      if (length(estimate_scores) > 0) {
        immDeconvCols <- c(immDeconvCols,estimate_scores)
        decon_score_cols <- c(decon_score_cols,estimate_scores)
      }
      quantiseq_scores <- grep("quantiseq_PreProcessedScore$",metacol_feature,value = T)
      if (length(quantiseq_scores) > 0) {
        immDeconvCols <- c(immDeconvCols,quantiseq_scores)
        decon_score_cols <- c(decon_score_cols,quantiseq_scores)
      }
      xcell_scores <- grep("xcell_PreProcessedScore$",metacol_feature,value = T)
      if (length(xcell_scores) > 0) {
        immDeconvCols <- c(immDeconvCols,xcell_scores)
        decon_score_cols <- c(decon_score_cols,xcell_scores)
      }
      epic_scores <- grep("epic_PreProcessedScore$",metacol_feature,value = T)
      if (length(epic_scores) > 0) {
        immDeconvCols <- c(immDeconvCols,epic_scores)
        decon_score_cols <- c(decon_score_cols,epic_scores)
      }
      abis_scores <- grep("abis_PreProcessedScore$",metacol_feature,value = T)
      if (length(abis_scores) > 0) {
        immDeconvCols <- c(immDeconvCols,abis_scores)
        decon_score_cols <- c(decon_score_cols,abis_scores)
      }
      cibersort_scores <- grep("cibersort_PreProcessedScore$",metacol_feature,value = T)
      if (length(cibersort_scores) > 0) {
        immDeconvCols <- c(immDeconvCols,cibersort_scores)
        decon_score_cols <- c(decon_score_cols,cibersort_scores)
      }
      cibersort_abs_scores <- grep("cibersort_abs_PreProcessedScore$",metacol_feature,value = T)
      if (length(cibersort_abs_scores) > 0) {
        immDeconvCols <- c(immDeconvCols,cibersort_abs_scores)
        decon_score_cols <- c(decon_score_cols,cibersort_abs_scores)
      }
      other_PreProcessedScores <- setdiff(PreProcessed_meta_cols,immDeconvCols)
      if (length(other_PreProcessedScores) > 0) {
        decon_score_cols <- c(decon_score_cols,other_PreProcessedScores)
      }
      
      decon_score_cols <- gsub("_PreProcessedScore","",decon_score_cols)
      decon_score_cols
      
    }
    
  })
  
  ## features with pre process designation removed
  metacol_feature <- reactive({
    
    decon_score_cols <- decon_score_cols()
    metacol_feature_og <- metacol_feature_og()
    metacol_feature_noPreProcessedScore <- c(grep("_PreProcessedScore$",metacol_feature_og,value = T, invert = T))
    metacol_feature <- c(metacol_feature_noPreProcessedScore,decon_score_cols)
    metacol_feature
    
  })
  
  GeneSetTable_new <- reactive({
    
    metacol_feature <- metacol_feature_og()
    PreProcessed_meta_cols <- c(grep("_PreProcessedScore$",metacol_feature,value = T))
    immDeconvCols <- c()
    GeneSetTable_og_new <- GeneSetTable_og
    if (length(PreProcessed_meta_cols) > 0) {
      
      ## Check in immune deconvolution scores were pre-processed
      mcp_devonv_scores <- grep("mcp_counter_PreProcessedScore$",metacol_feature,value = T)
      if (length(mcp_devonv_scores) > 0) {
        immDeconvCols <- c(immDeconvCols,mcp_devonv_scores)
        mcp_decon_gstab <- data.frame(GeneSet_Database = "Pre-Processed Scores",
                                      GeneSet_Category = "Immune Deconvolution Cell Types",
                                      GeneSet_Sub_Category = "MCP Counter Deconvolution Method",
                                      GeneSet_Name = gsub("_PreProcessedScore","",mcp_devonv_scores))
        GeneSetTable_og_new <- rbind(GeneSetTable_og_new,mcp_decon_gstab)
      }
      
      estimate_scores <- grep("estimate_PreProcessedScore$",metacol_feature,value = T)
      if (length(estimate_scores) > 0) {
        immDeconvCols <- c(immDeconvCols,estimate_scores)
        estimate_decon_gstab <- data.frame(GeneSet_Database = "Pre-Processed Scores",
                                           GeneSet_Category = "Immune Deconvolution Cell Types",
                                           GeneSet_Sub_Category = "ESTIMATE Deconvolution Method",
                                           GeneSet_Name = gsub("_PreProcessedScore","",estimate_scores))
        GeneSetTable_og_new <- rbind(GeneSetTable_og_new,estimate_decon_gstab)
      }
      
      quantiseq_scores <- grep("quantiseq_PreProcessedScore$",metacol_feature,value = T)
      if (length(quantiseq_scores) > 0) {
        immDeconvCols <- c(immDeconvCols,quantiseq_scores)
        quantiseq_decon_gstab <- data.frame(GeneSet_Database = "Pre-Processed Scores",
                                            GeneSet_Category = "Immune Deconvolution Cell Types",
                                            GeneSet_Sub_Category = "quanTIseq Deconvolution Method",
                                            GeneSet_Name = gsub("_PreProcessedScore","",quantiseq_scores))
        GeneSetTable_og_new <- rbind(GeneSetTable_og_new,quantiseq_decon_gstab)
      }
      
      xcell_scores <- grep("xcell_PreProcessedScore$",metacol_feature,value = T)
      if (length(xcell_scores) > 0) {
        immDeconvCols <- c(immDeconvCols,xcell_scores)
        xcell_decon_gstab <- data.frame(GeneSet_Database = "Pre-Processed Scores",
                                        GeneSet_Category = "Immune Deconvolution Cell Types",
                                        GeneSet_Sub_Category = "Xcell Deconvolution Method",
                                        GeneSet_Name = gsub("_PreProcessedScore","",xcell_scores))
        GeneSetTable_og_new <- rbind(GeneSetTable_og_new,xcell_decon_gstab)
      }
      
      epic_scores <- grep("epic_PreProcessedScore$",metacol_feature,value = T)
      if (length(epic_scores) > 0) {
        immDeconvCols <- c(immDeconvCols,epic_scores)
        epic_decon_gstab <- data.frame(GeneSet_Database = "Pre-Processed Scores",
                                       GeneSet_Category = "Immune Deconvolution Cell Types",
                                       GeneSet_Sub_Category = "EPIC Deconvolution Method",
                                       GeneSet_Name = gsub("_PreProcessedScore","",epic_scores))
        GeneSetTable_og_new <- rbind(GeneSetTable_og_new,epic_decon_gstab)
      }
      
      abis_scores <- grep("abis_PreProcessedScore$",metacol_feature,value = T)
      if (length(abis_scores) > 0) {
        immDeconvCols <- c(immDeconvCols,abis_scores)
        abis_decon_gstab <- data.frame(GeneSet_Database = "Pre-Processed Scores",
                                       GeneSet_Category = "Immune Deconvolution Cell Types",
                                       GeneSet_Sub_Category = "ABIS Deconvolution Method",
                                       GeneSet_Name = gsub("_PreProcessedScore","",abis_scores))
        GeneSetTable_og_new <- rbind(GeneSetTable_og_new,abis_decon_gstab)
      }
      
      cibersort_scores <- grep("cibersort_PreProcessedScore$",metacol_feature,value = T)
      if (length(cibersort_scores) > 0) {
        immDeconvCols <- c(immDeconvCols,cibersort_scores)
        cibersort_decon_gstab <- data.frame(GeneSet_Database = "Pre-Processed Scores",
                                            GeneSet_Category = "Immune Deconvolution Cell Types",
                                            GeneSet_Sub_Category = "CIBERSORT Deconvolution Method",
                                            GeneSet_Name = gsub("_PreProcessedScore","",cibersort_scores))
        GeneSetTable_og_new <- rbind(GeneSetTable_og_new,cibersort_decon_gstab)
      }
      
      cibersort_abs_scores <- grep("cibersort_abs_PreProcessedScore$",metacol_feature,value = T)
      if (length(cibersort_abs_scores) > 0) {
        immDeconvCols <- c(immDeconvCols,cibersort_abs_scores)
        cibersort_abs_decon_gstab <- data.frame(GeneSet_Database = "Pre-Processed Scores",
                                                GeneSet_Category = "Immune Deconvolution Cell Types",
                                                GeneSet_Sub_Category = "CIBERSORT ABS Deconvolution Method",
                                                GeneSet_Name = gsub("_PreProcessedScore","",cibersort_abs_scores))
        GeneSetTable_og_new <- rbind(GeneSetTable_og_new,cibersort_abs_decon_gstab)
      }
      
      other_PreProcessedScores <- setdiff(PreProcessed_meta_cols,immDeconvCols)
      if (length(other_PreProcessedScores) > 0) {
        Other_scores_gstab <- data.frame(GeneSet_Database = "Pre-Processed Scores",
                                         GeneSet_Category = "Pre-Processed Scores",
                                         GeneSet_Sub_Category = "Pre-Processed Scores",
                                         GeneSet_Name = gsub("_PreProcessedScore","",other_PreProcessedScores))
        GeneSetTable_og_new <- rbind(GeneSetTable_og_new,Other_scores_gstab)
      }
    }
    GeneSetTable_og_new
    
    
  })
  
  
  
  ####----Data Tables----####
  
  output$ExprFile_Preview <- DT::renderDataTable({
    
    expr <- expr_react()
    DT::datatable(expr,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 5,
                                 scrollX = T),
                  rownames = T)
    
  })
  
  output$ClinFile_Preview <- DT::renderDataTable({
    
    clin <- clin_react()
    DT::datatable(clin,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 5,
                                 scrollX = T),
                  rownames = F)
    
  })
  
  output$ClinParamFile_Preview <- DT::renderDataTable({
    
    clinP <- clinP_react()
    DT::datatable(clinP,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 5,
                                 scrollX = T),
                  rownames = F)
    
  })
  
  #output$rendAlertMessage <- renderUI({
  #  
  #  expr <- expr_react()
  #  clin <- clin_react()
  #  clinP <- clinP_react()
  #  
  #  if (any(sort(colnames(expr)) != sort(clin[,"SampleName"]))) {
  #    SampNameError <- "Sample Names are missing or do not match between expression and clinical data."
  #  }
  #  
  #  if (ncol(clin) = nrow(clinP)) {
  #    
  #  }
  #  
  #})
  
  ####----Survival Analysis Panel----####
  
  ####----Render UI----####
  
  ####----Sample Selection----####
  
  ## Select sample type to subset samples by - only render if more than one sample type
  output$rendSampleTypeSelection <- renderUI({
    
    metacol_sampletype <- metacol_sampletype()
    clin <- clin_react()
    
    if (length(unique(clin[,metacol_sampletype])) > 1) {
      
      SampleTypeChoices <- unique(clin[,metacol_sampletype])
      SampleTypeChoices <- c(SampleTypeChoices,"All_Sample_Types")
      selectInput("SampleTypeSelection",paste("Select Sample Type (",metacol_sampletype,"):",sep = ""),
                  choices = SampleTypeChoices)
      
    }
    
  })
  
  ## Select primary feature to look at
  output$rendFeatureSelection <- renderUI({
    
    metacol_sampletype <- metacol_sampletype()
    meta <- clin_react()
    metacol_feature <- metacol_feature()
    
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      
      if (input$SampleTypeSelection == "All_Sample_Types") {
        
        FeatureChoices <- c(metacol_sampletype,metacol_feature,"Show All Samples")
        MetaClass <- sapply(as_tibble(meta), typeof)
        MetaClass_num <- names(MetaClass[which(MetaClass %in% c("integer","double"))])
        IntMetaCols <- apply(meta[,MetaClass_num],2,function(x) any(round(x) != x))
        MaybeShow <- names(IntMetaCols[which(IntMetaCols == FALSE)])
        ## Checks if integer columns might be categorical or continuous
        MaybeShowCols <- apply(meta[,MaybeShow],2,function(x) any(length(levels(as.factor(x)))<(nrow(meta)*0.75)))
        MaybeNotShow <- names(MaybeShowCols[which(MaybeShowCols == FALSE)])
        DoNotShow <- c(MaybeNotShow,names(IntMetaCols[which(IntMetaCols == TRUE)]))
        FeatureChoices <- FeatureChoices[!FeatureChoices %in% DoNotShow]
        selectInput("FeatureSelection","Select Feature to Subset Samples:", choices = FeatureChoices, selected = "Show All Samples")
        
      }
      else if (input$SampleTypeSelection != "All_Sample_Types") {
        
        FeatureChoices <- c(metacol_feature,"Show All Samples")
        MetaClass <- sapply(as_tibble(meta), typeof)
        MetaClass_num <- names(MetaClass[which(MetaClass %in% c("integer","double"))])
        IntMetaCols <- apply(meta[,MetaClass_num],2,function(x) any(round(x) != x))
        MaybeShow <- names(IntMetaCols[which(IntMetaCols == FALSE)])
        ## Checks if integer columns might be categorical or continuous
        MaybeShowCols <- apply(meta[,MaybeShow],2,function(x) any(length(levels(as.factor(x)))<(nrow(meta)*0.75)))
        MaybeNotShow <- names(MaybeShowCols[which(MaybeShowCols == FALSE)])
        DoNotShow <- c(MaybeNotShow,names(IntMetaCols[which(IntMetaCols == TRUE)]))
        FeatureChoices <- FeatureChoices[!FeatureChoices %in% DoNotShow]
        selectInput("FeatureSelection","Select Feature to Subset Samples:", choices = FeatureChoices, selected = "Show All Samples")
        
      }
      
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      
      FeatureChoices <- c(metacol_feature,"Show All Samples")
      MetaClass <- sapply(as_tibble(meta), typeof) #Get column type
      MetaClass_num <- names(MetaClass[which(MetaClass %in% c("integer","double"))]) #subset column names that are types of numeric
      IntMetaCols <- apply(meta[,MetaClass_num],2,function(x) any(round(x) != x)) #Find which column have decimals
      MaybeShow <- names(IntMetaCols[which(IntMetaCols == FALSE)]) #subset column names that are not decimals
      ## Checks if integer columns might be categorical or continuous
      MaybeShowCols <- apply(meta[,MaybeShow],2,function(x) any(length(levels(as.factor(x)))<(nrow(meta)*0.75)))
      MaybeNotShow <- names(MaybeShowCols[which(MaybeShowCols == FALSE)])
      DoNotShow <- c(MaybeNotShow,names(IntMetaCols[which(IntMetaCols == TRUE)]))
      FeatureChoices <- FeatureChoices[!FeatureChoices %in% DoNotShow]
      selectInput("FeatureSelection","Select Feature to Subset Samples:", choices = FeatureChoices, selected = "Show All Samples")
      
    }
    
  })
  
  ## Select primary features condition to look at - All not working yet
  output$rendSubFeatureSelection <- renderUI({
    
    req(input$FeatureSelection)
    metacol_sampletype <- metacol_sampletype()
    HasSampleTypes <- FALSE
    clin <- clin_react()
    if (!is.null(metacol_sampletype)) {
      if (length(unique(clin[,metacol_sampletype])) > 1) {
        HasSampleTypes <- TRUE
      }
      else if (length(unique(clin[,metacol_sampletype])) <= 1) {
        HasSampleTypes <- FALSE
      }
    }
    
    if (HasSampleTypes == TRUE) {
      SampleType <- input$SampleTypeSelection
    }
    if (HasSampleTypes == FALSE) {
      SampleType <- "All_Sample_Types"
    }
    
    Feature <- input$FeatureSelection
    
    if (SampleType == "All_Sample_Types") {
      clin <- clin
    }
    if (SampleType != "All_Sample_Types") {
      clin <- clin[which(clin[,metacol_sampletype] == SampleType),]
    }
    
    if (Feature != "Show All Samples") {
      
      SubFeatureChoices <- unique(clin[,Feature])
      # Sort options, will put 1,TRUE,yes before 0,FASLE,no, so the 'positive' value is the initial selected - puts NA last
      SubFeatureChoices <- sort(SubFeatureChoices, decreasing = T, na.last = T)
      selectInput("subFeatureSelection","Feature Condition:", choices = SubFeatureChoices)
      
    }
    
  })
  
  ####----Univariate----####
  
  output$rendSurvivalFeatureSingle <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    metacol_sampletype <- metacol_sampletype()
    meta <- clin_react()
    metacol_feature <- metacol_feature()
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        if (input$FeatureSelection != "Show All Samples") {
          metacol_feature <- metacol_feature[-which(metacol_feature == input$FeatureSelection)]
        }
        metacol_feature <- c(metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
        selectInput("SingleSurvivalFeature","Select Feature:",
                    choices = metacol_feature, selected = "MedianCutP")
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        
        SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
        if (input$FeatureSelection != "Show All Samples") {
          SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        }
        selectInput("SingleSurvivalFeature","Select Feature:",
                    choices = SurvFeatChoices2, selected = "MedianCutP")
        
      }
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
      if (input$FeatureSelection != "Show All Samples") {
        SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
      }
      selectInput("SingleSurvivalFeature","Select Feature:",
                  choices = SurvFeatChoices2, selected = "MedianCutP")
    }
    
  })
  
  output$rendUniVarContHiLoCheck <- renderUI({
    
    if (input$UniVarContCheck == T) {
      checkboxInput("UniVarContHiLoCheck","Continuous Feature as Median Cut-Point",value = T)
    }
    
  })
  
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
  
  ####----BiVar Add----####
  
  output$rendSurvivalFeatureBi1 <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    metacol_feature <- metacol_feature()
    meta <- clin_react()
    metacol_sampletype <- metacol_sampletype()
    SurvFeatChoices <- c(metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
    SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        if (input$FeatureSelection != "Show All Samples") {
          SurvFeatChoices <- SurvFeatChoices[-which(SurvFeatChoices == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi1","Select Feature 1:",
                    choices = SurvFeatChoices, selected = "MedianCutP")
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        if (input$FeatureSelection != "Show All Samples") {
          SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi1","Select Feature 1:",
                    choices = SurvFeatChoices2, selected = "MedianCutP")
        
      }
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (input$FeatureSelection != "Show All Samples") {
        SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
      }
      selectInput("SurvivalFeatureBi1","Select Feature 1:",
                  choices = SurvFeatChoices2, selected = "MedianCutP")
    }
    
  })
  
  output$rendSurvivalFeatureBi2 <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    meta <- clin_react()
    metacol_feature <- metacol_feature()
    metacol_sampletype <- metacol_sampletype()
    SurvFeatChoices <- c(metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
    SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        if (input$FeatureSelection != "Show All Samples") {
          SurvFeatChoices <- SurvFeatChoices[-which(SurvFeatChoices == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi2","Select Feature 2:",
                    choices = SurvFeatChoices, selected = 1)
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        if (input$FeatureSelection != "Show All Samples") {
          SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi2","Select Feature 2:",
                    choices = SurvFeatChoices2, selected = 1)
        
      }
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (input$FeatureSelection != "Show All Samples") {
        SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
      }
      selectInput("SurvivalFeatureBi2","Select Feature 2:",
                  choices = SurvFeatChoices2, selected = 1)
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
  
  ####----BiVar Int----####
  
  output$rendSurvivalFeatureBi1Inter <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    metacol_feature <- metacol_feature()
    meta <- clin_react()
    metacol_sampletype <- metacol_sampletype()
    SurvFeatChoices <- c(metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
    SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        if (input$FeatureSelection != "Show All Samples") {
          SurvFeatChoices <- SurvFeatChoices[-which(SurvFeatChoices == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi1Inter","Select Feature 1:",
                    choices = SurvFeatChoices, selected = "MedianCutP")
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        if (input$FeatureSelection != "Show All Samples") {
          SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi1Inter","Select Feature 1:",
                    choices = SurvFeatChoices2, selected = "MedianCutP")
        
      }
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (input$FeatureSelection != "Show All Samples") {
        SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
      }
      selectInput("SurvivalFeatureBi1Inter","Select Feature 1:",
                  choices = SurvFeatChoices2, selected = "MedianCutP")
    }
    
  })
  
  output$rendSurvivalFeatureBi2Inter <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    metacol_feature <- metacol_feature()
    meta <- clin_react()
    metacol_sampletype <- metacol_sampletype()
    SurvFeatChoices <- c(metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
    SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        if (input$FeatureSelection != "Show All Samples") {
          SurvFeatChoices <- SurvFeatChoices[-which(SurvFeatChoices == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi2Inter","Select Feature 2:",
                    choices = SurvFeatChoices, selected = 1)
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        if (input$FeatureSelection != "Show All Samples") {
          SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi2Inter","Select Feature 2:",
                    choices = SurvFeatChoices2, selected = 1)
        
      }
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (input$FeatureSelection != "Show All Samples") {
        SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
      }
      selectInput("SurvivalFeatureBi2Inter","Select Feature 2:",
                  choices = SurvFeatChoices2, selected = 1)
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
  
  ####---MultiVar----####
  
  output$rendSurvivalFeature <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    metacol_feature <- metacol_feature()
    meta <- clin_react()
    metacol_sampletype <- metacol_sampletype()
    SurvFeatChoices <- c(metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
    SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        if (input$FeatureSelection != "Show All Samples") {
          SurvFeatChoices <- SurvFeatChoices[-which(SurvFeatChoices == input$FeatureSelection)]
        }
        selectInput("SurvivalFeature","Select Feature(s):",
                    choices = SurvFeatChoices, multiple = T, selected = "MedianCutP")
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        if (input$FeatureSelection != "Show All Samples") {
          SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        }
        selectInput("SurvivalFeature","Select Feature(s):",
                    choices = SurvFeatChoices2, multiple = T, selected = "MedianCutP")
        
      }
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (input$FeatureSelection != "Show All Samples") {
        SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
      }
      selectInput("SurvivalFeature","Select Feature(s):",
                  choices = SurvFeatChoices2, multiple = T, selected = "MedianCutP")
    }
    
  })
  
  ## Select survival type selection based on meta columns (ex. OS vs EFS)
  output$rendSurvivalType_time <- renderUI({
    
    ## Only show if more than one option
    metacol_survtime <- metacol_survtime()
    #if (length(metacol_survtime) > 1) {
      selectInput("SurvivalType_time","Survival Time Data:", choices = metacol_survtime)
    #}
    
  })
  
  survTime <- reactive({
    
    metacol_survtime <- metacol_survtime()
    if (length(metacol_survtime) > 1) {
      metacol_survtime <- input$SurvivalType_time
    }
    metacol_survtime
    
  })
  
  ## Select survival type selection based on meta columns (ex. OS vs EFS)
  output$rendSurvivalType_id <- renderUI({
    
    ## Only show if more than one option
    metacol_survid <- metacol_survid()
    #if (length(metacol_survid) > 1) {
      selectInput("SurvivalType_id","Survival ID Data:", choices = metacol_survid)
    #}
    
  })
  
  survID <- reactive({
    
    metacol_survid <- metacol_survid()
    if (length(metacol_survid) > 1) {
      metacol_survid <- input$SurvivalType_id
    }
    metacol_survid
    
  })
  
  ## Select ssGSEA function scoring method
  output$rendScoreMethodBox <- renderUI({
    
      selectInput("ScoreMethod","Scoring Method",choices = c("ssgsea","gsva","zscore","plage"))
    
  })
  
  output$rendGeneSetCat_Select <- renderUI({
    
    if (length(decon_score_cols()) > 0) {
      selectInput("GeneSetCat_Select","Select Category",
                  choices = c("MSigDB","LINCS L1000","Cell Marker","ER Stress","Immune Signatures","TCGA","Pre-Processed Scores"))
    }
    else {
      selectInput("GeneSetCat_Select","Select Category",
                  choices = c("MSigDB","LINCS L1000","Cell Marker","ER Stress","TCGA","Immune Signatures"))
    }
    
  })
  
  ## View Gene Set table
  output$rendGeneSetTable <- renderUI({
    
    div(DT::dataTableOutput("GeneSetTable"), style = "font-size:10px")
    
    
  })
  
  ## View Gene - "Gene Set" table
  output$rendGeneGeneSetTable <- renderUI({
    
    div(DT::dataTableOutput("geneGeneSetTable"), style = "font-size:10px")
    
  })
  
  output$renduserGeneSet <- renderUI({
    
    if (input$UserGSoption == "Gene Set Upload") {
      
      fileInput("userGeneSet","Gene Set Upload", accept = c(".gmt",".tsv",".txt",".RData"))
      
    }
    
  })
  
  output$renduserGeneSetTextName <- renderUI({
    
    if (input$UserGSoption == "Text Box Input") {
      
      textInput("userGeneSetTextName","Custom Gene Set Name", value = "Custom_Geneset")
      
    }
    
  })
  
  output$renduserGeneSetText <- renderUI({
    
    if (input$UserGSoption == "Text Box Input") {
      
      textInput("userGeneSetText","Gene Symbols", placeholder = "Comma, space, or tab delimited")
      
    }
    
  })
  
  ## View User Gene Set table
  output$rendUserGeneSetTable <- renderUI({
    
    req(input$userGeneSet)
    div(DT::dataTableOutput("userGeneSetTable"), style = "font-size:10px")
    
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
  
  output$rendGenesInGeneSetTab <- renderUI({
    
    if (input$GeneSetTabs == 1) {
      
      if (!is.null(input$ViewGeneSetGenes)) {
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
    }
    
  })
  
  ## Survival X-axis Length
  output$rendSurvXaxis <- renderUI({
    
    meta_ssgsea <- ssGSEAmeta()
    #if (is.null(input$SurvivalType_time) == TRUE & is.null(input$SurvivalType_id) == TRUE) {
    #  surv_time_col <- metacol_survtime[1]
    #  surv_id_col <- metacol_survid[1]
    #}
    #if (is.null(input$SurvivalType_time) == FALSE & is.null(input$SurvivalType_id) == FALSE) {
    #  surv_time_col <- input$SurvivalType_time
    #  surv_id_col <- input$SurvivalType_id
    #}
    surv_time_col <- input$SurvivalType_time
    surv_id_col <- input$SurvivalType_id
    max_time <- ceiling(max(meta_ssgsea[,surv_time_col])/365.25)
    numericInput("SurvXaxis","X-Axis Limit (years)", value = max_time)
    
    
  })
  
  output$rendQuartHRtab <- renderUI({
    
    div(shinycssloaders::withSpinner(tableOutput("SQuartileHRtab"), type = 7, size = 0.5), style = "font-size:12px")
    
  })
  output$rendBINHRtab <- renderUI({
    
    div(shinycssloaders::withSpinner(tableOutput("SBinaryHRtab"), type = 7, size = 0.5), style = "font-size:12px")
    
  })
  output$rendQuantHRtab <- renderUI({
    
    div(shinycssloaders::withSpinner(tableOutput("SQuantileHRtab"), type = 7, size = 0.5), style = "font-size:12px")
    
  })
  output$rendQuantHRtab2 <- renderUI({
    
    div(shinycssloaders::withSpinner(tableOutput("SQuantileHR2tab"), type = 7, size = 0.5), style = "font-size:12px")
    
  })
  output$rendCutPointHRtab <- renderUI({
    
    div(shinycssloaders::withSpinner(tableOutput("CutPointHRtab"), type = 7, size = 0.5), style = "font-size:12px")
    
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
  
  ####----Data Exploration----####
  
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
  
  ## Download button for subset meta
  output$DnldMetaButon <- renderUI({
    
    downloadButton("dnldMeta", "Download Meta Subset")
    
  })
  
  ## Download button for subset expression
  output$DnldExprButon <- renderUI({
    
    downloadButton("dnldExpr", "Download Expression Subset")
    
  })
  
  output$rendScatterFeature <- renderUI({
    
    meta <- clin_react()
    metacol_sampletype <- metacol_sampletype()
    metacol_feature <- metacol_feature()
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        selectInput("ScatterFeature","Select Feature:",
                    choices = metacol_feature, selected = 1)
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        selectInput("ScatterFeature","Select Feature:",
                    choices = c(metacol_sampletype,metacol_feature), selected = 1)
        
      }
      
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      
      selectInput("ScatterFeature","Select Feature:",
                  choices = metacol_feature, selected = 1)
      
    }
    
  })
  
  output$rendScatterColor <- renderUI({
    
    if (input$ColorScatterChoice == "Feature") {
      
      geneset <- gs_react()
      geneset_name <- names(geneset)
      metacol_sampletype <- metacol_sampletype()
      meta <- clin_react()
      metacol_feature <- metacol_feature()
      metacol_survid <- metacol_survid()
      survID <- survID()
      
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        if (input$SampleTypeSelection != "All_Sample_Types") {
          
          if (input$FeatureSelection != "Show All Samples") {
            metacol_feature <- metacol_feature[-which(metacol_feature == input$FeatureSelection)]
          }
          metacol_feature <- c(metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
          metacol_feature <- c(metacol_survid,metacol_feature)
          selectInput("ScatterColor","Color By:",
                      choices = metacol_feature, selected = 1)
          
        }
        else if (input$SampleTypeSelection == "All_Sample_Types") {
          
          
          SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
          if (input$FeatureSelection != "Show All Samples") {
            SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
          }
          SurvFeatChoices2 <- c(metacol_survid,SurvFeatChoices2)
          selectInput("ScatterColor","Color By:",
                      choices = SurvFeatChoices2, selected = 1)
          
        }
      }
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
        if (input$FeatureSelection != "Show All Samples") {
          SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        }
        SurvFeatChoices2 <- c(metacol_survid,SurvFeatChoices2)
        selectInput("ScatterColor","Color By:",
                    choices = SurvFeatChoices2, selected = 1)
      }
      
    }
    else if (input$ColorScatterChoice == "Single Color") {
      selectInput("ScatterColor","Color:",
                  choices = colors(), selected = "cadetblue")
    }
    
  })
  
  ## Feature selection for boxplot
  output$rendBoxplotFeature <- renderUI({
    
    metacol_sampletype <- metacol_sampletype()
    metacol_feature <- metacol_feature()
    meta <- clin_react()
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        selectInput("BoxplotFeature","Select Feature:",
                    choices = metacol_feature, selected = 1)
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        selectInput("BoxplotFeature","Select Feature:",
                    choices = c(metacol_sampletype,metacol_feature), selected = 1)
        
      }
      
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      
      selectInput("BoxplotFeature","Select Feature:",
                  choices = metacol_feature, selected = 1)
      
    }
    
  })
  
  ## Feature selection for heatmap
  output$rendHeatmapFeature <- renderUI({
    
    metacol_sampletype <- metacol_sampletype()
    metacol_feature <- metacol_feature()
    meta <- clin_react()
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
  
  output$heatmap_error_message <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    decon_score_cols <- decon_score_cols()
    if (geneset_name %in% decon_score_cols | input$GeneSetTabs == 2) {
      p("Heatmap not available for Pre-Processed scores or single gene analysis")
    }
    
  })
  output$heatmap_error_message2 <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    decon_score_cols <- decon_score_cols()
    if (geneset_name %in% decon_score_cols | input$GeneSetTabs == 2) {
      p("Heatmap not available for Pre-Processed scores or single gene analysis")
    }
    
  })
  
  ####----Data Tables----####
  
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
    
    GeneGS_table <- GeneGS_table_react()
    DT::datatable(GeneGS_table,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  selection = list(mode = 'single', selected = 1),
                  rownames = F)
    
  })
  
  output$GenesInGeneSetTab <- DT::renderDataTable({
    
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
  
  ## Render Meta Table
  output$MetaTable <- renderDataTable({
    
    # User selections
    SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    SubFeature <- input$subFeatureSelection
    meta <- ssGSEAmeta()
    surv_time_col <- input$SurvivalType_time
    surv_id_col <- input$SurvivalType_id
    # Initial columns selected is NULL
    if (is.null(input$MetaTableCols) == TRUE) {
      
      userMetaCols <- NULL
      
    }
    # When meta tab is selected
    else if (is.null(input$MetaTableCols) == FALSE) {
      
      # If columns selected, add to list
      if (length(input$MetaTableCols) > 0) {
        userMetaCols <- input$MetaTableCols
      }
      else if (input$MetaTableCols == "") {
        userMetaCols <- NULL
      }
      
    }
    
    metaCols <- colnames(meta)[1] #select sample name column automatically
    if (Feature == "Show All Samples") {
      #userMetaCols <- userMetaCols[userMetaCols != Feature] #remove condition column from user selection because it is automatically added
      metaCols <- c(metaCols,surv_time_col,surv_id_col,userMetaCols) #combine column names selected
    }
    else if (Feature != "Show All Samples") {
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
  
  output$ssgseaDensityTable <- renderDataTable({
    
    geneset <- gs_react()
    GeneSet <- names(geneset)
    ssgsea_meta <- ssGSEAmeta()
    table <- ssgsea_meta[,c("SampleName",GeneSet)]
    DT::datatable(table,
                  options = list(scrollY = T),
                  rownames = F)
    
  })
  
  output$SboxplotTable <- renderDataTable({
    
    geneset <- gs_react()
    GeneSet <- names(geneset)
    ssGSEA_meta <- SboxplotReact()
    surv_time_col <- input$SurvivalType_time
    surv_id_col <- input$SurvivalType_id
    boxTab <- ssGSEA_meta[,c("SampleName",surv_time_col,surv_id_col,GeneSet,"SurvivalCutoff")]
    DT::datatable(boxTab,
                  options = list(paging = F),
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
  
  ####----Reactives----####
  
  GeneSetTable_React <- reactive({
    
    GeneSetTable_new <- GeneSetTable_new()
    GS_database <- "MSigDB"
    if (!is.null(input$GeneSetCat_Select)) {
      GS_database <- input$GeneSetCat_Select
    }
    sub_tab <- GeneSetTable_new[which(GeneSetTable_new[,1] == GS_database),]
    new_tab <- sub_tab[,-1]
    new_tab
    
  })
  
  ## Render User Gene Set Table - Backend
  userGeneSetTable_backend <- reactive({
    
    gs.u <- input$userGeneSet
    ext <- tools::file_ext(gs.u$datapath)
    req(gs.u)
    validate(need(ext == c("gmt","tsv","txt", "RData"), "Please upload .gmt, .tsv, .txt, or .RData file"))
    
    # If user provides GMT file
    if (ext == "gmt") {
      gmt <- clusterProfiler::read.gmt(gs.u$datapath)
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
      gmt <- clusterProfiler::read.gmt(gs.u$datapath)
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
    
    GeneGS_table <- GeneGS_table_react()
    GeneSetTable <- GeneSetTable_React()
    clin <- clin_react()
    RowSelected <- 1
    RowSelected2 <- 1
    RowSelected3 <- 1
    if (!is.null(input$GeneSetTable_rows_selected)) {
      RowSelected <- input$GeneSetTable_rows_selected
    }
    if (!is.null(input$geneGeneSetTable_rows_selected)) {
      RowSelected2 <- input$geneGeneSetTable_rows_selected
    }
    if (!is.null(input$userGeneSetTable_rows_selected)) {
      RowSelected3 <- input$userGeneSetTable_rows_selected
    }
    
    if (input$GeneSetTabs == 1) {
      geneset_name <- GeneSetTable[RowSelected,3]
      if (geneset_name %in% names(gs)) {
        geneset <- gs[geneset_name]
      }
      else if (geneset_name %in% decon_score_cols()) {
        geneset <- list(clin[,geneset_name])
        names(geneset) <- geneset_name
      }
      geneset 
    }
    else if (input$GeneSetTabs == 2) {
      gene <- GeneGS_table[RowSelected2,1]
      geneset <- list(gene = gene)
      names(geneset)[1] <- gene
      geneset
    }
    else if (input$GeneSetTabs == 3) {
      if (input$UserGSoption == "Gene Set Upload") {
        req(input$userGeneSet)
        gs_u <- user_gs()
        uGS_table <- userGeneSetTable_backend()
        geneset <- gs_u[(uGS_table[RowSelected3,1])]
      }
      else if (input$UserGSoption == "Text Box Input") {
        req(input$userGeneSetText)
        geneset <- user_gs_text()
      }
      geneset
    }
    
    
    
    
    
  })
  
  ## Meta subset reactive - "All Primary features" not working yet
  metaSub <- reactive({
    
    req(input$FeatureSelection)
    
    clin <- clin_react()
    
    metacol_sampletype <- metacol_sampletype()
    if (length(unique(clin[,metacol_sampletype])) > 1) {
      SampleType <- input$SampleTypeSelection
    }
    if (length(unique(clin[,metacol_sampletype])) <= 1) {
      SampleType <- "All_Sample_Types"
    }
    Feature <- input$FeatureSelection
    SubFeature <- input$subFeatureSelection
    
    if (SampleType == "All_Sample_Types") {
      clin <- clin
    }
    if (SampleType != "All_Sample_Types") {
      clin <- clin[which(clin[,metacol_sampletype] == SampleType),]
    }
    
    if (Feature != "Show All Samples") {
      clin <- clin[which(clin[,Feature] == SubFeature),]
    }
    if (Feature == "Show All Samples") {
      clin <- clin
    }
    clin
    
  })
  
  ## Expression subset reactive
  exprSub <- reactive({
    
    meta <- metaSub()
    expr <- expr_react()
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
    surv_time_col <- survTime()
    surv_id_col <- survID()
    ## Remove rows with NA in survival column
    meta <- meta[!is.na(meta[,surv_time_col]),]
    meta <- meta[!is.na(meta[,surv_id_col]),]
    
    ## Re-subset expression matrix
    samples <- meta[,1]
    expr_sub <- expr[,colnames(expr) %in% samples, drop = F]
    expr_mat <- as.matrix(expr_sub)
    rownames(expr_mat) <- rownames(expr_sub)
    colnames(expr_mat) <- colnames(expr_sub)
    
    if (input$GeneSetTabs == 1) {
      if (geneset_name %in% decon_score_cols()) {
        ssGSEA <- meta[,c("SampleName",geneset_name)]
        ssGSEA <- ssGSEA[!is.na(ssGSEA[,2]),]
        ssGSEA <- ssGSEA[which(ssGSEA[,2] != "Inf"),]
        ssGSEA[,2] <- as.numeric(ssGSEA[,2])
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
      
      expr_sub <- as.data.frame(expr_sub)
      expr_sub2 <- expr_sub[geneset_name,]
      expr_sub3 <- as.data.frame(t(expr_sub2))
      rownames(expr_sub3)[1]
      expr_sub3$SampleName <- rownames(expr_sub3)
      ssGSEA <- expr_sub3
    }
    
    ##--Optimal cutpoint analysis--##
    
    ## Subset columns needed for plot and rename for surv function
    meta_ssgsea_sdf <- merge(meta[,c("SampleName",surv_time_col,surv_id_col)],ssGSEA[,c("SampleName",geneset_name)], by = "SampleName")
    
    if (length(meta_ssgsea_sdf[,4][meta_ssgsea_sdf[,4] > 0])/length(meta_ssgsea_sdf[,4]) > 0.01) {
      if (length(meta_ssgsea_sdf[,4]) > 1) {
        res.cut <- survminer::surv_cutpoint(meta_ssgsea_sdf,time = surv_time_col, event = surv_id_col, variable = geneset_name, minprop = 0.01)
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
    if (geneset_name %in% decon_score_cols()) {
      ssGSEA <- ssGSEA[,-2] #remove score column so on merge the column does not duplicate
    }
    meta_ssGSEA <- merge(meta,ssGSEA, by = "SampleName", all = T)
    meta_ssGSEA
    
  })
  
  
  ####----Median Cut Point----####
  
  MedianCutP_react <- reactive({
    
    ## Assign variables
    surv_time_col <- survTime()
    surv_id_col <- survID()
    meta_ssgsea <- ssGSEAmeta()
    
    ## Subset columns needed for plot
    meta_ssgsea_sdf <- meta_ssgsea[,c("SampleName",surv_time_col,surv_id_col,"MedianCutP")]
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "ID"
    
    meta_ssgsea_sdf
    
  })
  
  MedianCutPTab_react <- reactive({
    
    meta_ssgsea_sdf <- MedianCutP_react()
    
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
      dplyr::select(label,estimate,ci,p.value)
    colnames(tab_df) <- c("Characteristic","Hazard Ratio","95% Confidence Interval","P.Value")
    
    tab_df
    
  })
  
  SplotBIN_react <- reactive({
    
    ## Assign variables
    meta_ssgsea_sdf <- MedianCutP_react()
    geneset <- gs_react()
    geneset_name <- names(geneset)
    Feature <- input$FeatureSelection
    scoreMethod <- input$ScoreMethod
    show_pval <- input$ShowPval
    ShowConfInt <- input$ShowConfInt
    xaxlim <- input$SurvXaxis * 365.25
    surv_time_col <- survTime()
    showLegend <- input$SurvLegendPos
    metacol_sampletype <- metacol_sampletype()
    clin <- clin_react()
    SampleType <- input$SampleTypeSelection
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
    if (length(unique(clin[,metacol_sampletype])) > 1) {
      if (SampleType == "All_Sample_Types") {
        if (Feature == "Show All Samples") {
          SampleTypeLab <- "All Features in All Patients\n"
        }
        if (Feature != "Show All Samples") {
          SampleTypeLab <- paste(Feature," in All Patients\n")
        }
      }
      else {
        if (Feature == "Show All Samples") {
          SampleTypeLab <- paste("All Features (",SampleType,") Patients\n",sep = "")
        }
        if (Feature != "Show All Samples") {
          SampleTypeLab <- paste(Feature," (",SampleType,") Patients\n",sep = "")
        }
      }
    }
    if (length(unique(clin[,metacol_sampletype])) <= 1) {
      if (Feature == "Show All Samples") {
        SampleTypeLab <- "All Features in All Patients\n"
      }
      if (Feature != "Show All Samples") {
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
    
    ## Determine Plot title
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
    ggsurv <- survminer::ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
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
    decon_score_cols <- decon_score_cols()
    CutPlabel <- "(Median Cut-Point)"
    
    cols_selec <- c("SampleName",geneset_name)
    ssgsea_scores <- ssgsea_meta[,cols_selec]
    
    quant_df2 <- data.frame(quantile(as.numeric(ssgsea_scores[,geneset_name]),probs = 0.5,na.rm = T))
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
    clin <- clin_react()
    metacol_sampletype <- metacol_sampletype()
    decon_score_cols <- decon_score_cols()
    
    ## Survival Type
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## Determine Feature and sub feature
    if (Feature != "Show All Samples") {
      SubFeature <- input$subFeatureSelection
      Feature <- paste("<b>",Feature,"</b> - <b>",SubFeature,"</b></li>",sep = "")
      line2 <- paste("<li>The dataset is filtered by ",Feature,sep = "")
    }
    if (Feature == "Show All Samples") {
      line2 <- NULL
    }
    ## determine Sample Type
    if (is.null(metacol_sampletype) == T) {
      SampleType <- ""
    }
    else if (is.null(metacol_sampletype) == F) {
      if (length(unique(clin[,metacol_sampletype])) > 1) {
        if (input$SampleTypeSelection == "All_Sample_Types") {
          SampleType <- " of all sample types.</li>"
        }
        else if (input$SampleTypeSelection != "All_Sample_Types") {
          SampleType <- paste(" of <b>",input$SampleTypeSelection,"</b> Patients.</li>",sep = "")
        }
      }
      if (length(unique(clin[,metacol_sampletype])) <= 1) {
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
      dplyr::select(label,estimate,ci,p.value)
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
    metacol_sampletype <- metacol_sampletype()
    clin <- clin_react()
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
    if (length(unique(clin[,metacol_sampletype])) > 1) {
      if (SampleType == "All_Sample_Types") {
        if (Feature == "Show All Samples") {
          SampleTypeLab <- "All Features in All Patients\n"
        }
        if (Feature != "Show All Samples") {
          SampleTypeLab <- paste(Feature," in All Patients\n")
        }
      }
      else {
        if (Feature == "Show All Samples") {
          SampleTypeLab <- paste("All Features (",SampleType,") Patients\n",sep = "")
        }
        if (Feature != "Show All Samples") {
          SampleTypeLab <- paste(Feature," (",SampleType,") Patients\n",sep = "")
        }
      }
    }
    if (length(unique(clin[,metacol_sampletype])) <= 1) {
      if (Feature == "Show All Samples") {
        SampleTypeLab <- "All Features in All Patients\n"
      }
      if (Feature != "Show All Samples") {
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
    ggsurv <- survminer::ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
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
    decon_score_cols <- decon_score_cols()
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
    metacol_sampletype <- metacol_sampletype()
    meta <- clin_react()
    decon_score_cols <- decon_score_cols()
    
    ## Survival Type
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## Determine Feature and sub feature
    if (Feature != "Show All Samples") {
      SubFeature <- input$subFeatureSelection
      Feature <- paste("<b>",Feature,"</b> - <b>",SubFeature,"</b></li>",sep = "")
      line2 <- paste("<li>The dataset is filtered by ",Feature,sep = "")
    }
    if (Feature == "Show All Samples") {
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
      dplyr::select(label,estimate,ci,p.value)
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
    metacol_sampletype <- metacol_sampletype()
    meta <- clin_react()
    decon_score_cols <- decon_score_cols()
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
        if (Feature == "Show All Samples") {
          SampleTypeLab <- "All Features in All Patients\n"
        }
        if (Feature != "Show All Samples") {
          SampleTypeLab <- paste(Feature," in All Patients\n")
        }
      }
      else {
        if (Feature == "Show All Samples") {
          SampleTypeLab <- paste("All Features (",SampleType,") Patients\n",sep = "")
        }
        if (Feature != "Show All Samples") {
          SampleTypeLab <- paste(Feature," (",SampleType,") Patients\n",sep = "")
        }
      }
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (Feature == "Show All Samples") {
        SampleTypeLab <- "All Features in All Patients\n"
      }
      if (Feature != "Show All Samples") {
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
    ggsurv <- survminer::ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
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
    metacol_sampletype <- metacol_sampletype()
    meta <- clin_react()
    decon_score_cols <- decon_score_cols()
    
    cols_selec <- c("SampleName",geneset_name)
    ssgsea_scores <- ssgsea_meta[,cols_selec]
    
    meta_ssgsea_sdf <- ssgsea_meta[,c("SampleName",surv_time_col,surv_id_col,geneset_name)]
    
    if (length(meta_ssgsea_sdf[,4][meta_ssgsea_sdf[,4] > 0])/length(meta_ssgsea_sdf[,4]) > 0.01) {
      res.cut <- survminer::surv_cutpoint(meta_ssgsea_sdf,time = surv_time_col, event = surv_id_col, variable = geneset_name, minprop = 0.01)
      cutp <- res.cut$cutpoint[["cutpoint"]]
      #res.cat <- surv_categorize(res.cut)
      #ssGSEA$OptimalCutP <- res.cat[,3]
      
      res.cut <- survminer::surv_cutpoint(meta_ssgsea_sdf,time = surv_time_col, event = surv_id_col, variable = geneset_name)
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
    metacol_sampletype <- metacol_sampletype()
    meta <- clin_react()
    decon_score_cols <- decon_score_cols()
    
    ## Survival Type
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## Determine Feature and sub feature
    if (Feature != "Show All Samples") {
      SubFeature <- input$subFeatureSelection
      Feature <- paste("<b>",Feature,"</b> - <b>",SubFeature,"</b></li>",sep = "")
      line2 <- paste("<li>The dataset is filtered by ",Feature,sep = "")
    }
    if (Feature == "Show All Samples") {
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
      dplyr::select(label,estimate,ci,p.value)
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
    metacol_sampletype <- metacol_sampletype()
    meta <- clin_react()
    decon_score_cols <- decon_score_cols()
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
        if (Feature == "Show All Samples") {
          SampleTypeLab <- "All Features in All Patients\n"
        }
        if (Feature != "Show All Samples") {
          SampleTypeLab <- paste(Feature," in All Patients\n")
        }
      }
      else {
        if (Feature == "Show All Samples") {
          SampleTypeLab <- paste("All Features (",SampleType,") Patients\n",sep = "")
        }
        if (Feature != "Show All Samples") {
          SampleTypeLab <- paste(Feature," (",SampleType,") Patients\n",sep = "")
        }
      }
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (Feature == "Show All Samples") {
        SampleTypeLab <- "All Features in All Patients\n"
      }
      if (Feature != "Show All Samples") {
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
    ggsurv <- survminer::ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
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
    metacol_sampletype <- metacol_sampletype()
    meta <- clin_react()
    decon_score_cols <- decon_score_cols()
    
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
    metacol_sampletype <- metacol_sampletype()
    meta <- clin_react()
    decon_score_cols <- decon_score_cols()
    
    ## Survival Type
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## Determine Feature and sub feature
    if (Feature != "Show All Samples") {
      SubFeature <- input$subFeatureSelection
      Feature <- paste("<b>",Feature,"</b> - <b>",SubFeature,"</b></li>",sep = "")
      line2 <- paste("<li>The dataset is filtered by ",Feature,sep = "")
    }
    if (Feature == "Show All Samples") {
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
      dplyr::select(label,estimate,ci,p.value)
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
    metacol_sampletype <- metacol_sampletype()
    meta <- clin_react()
    decon_score_cols <- decon_score_cols()
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
        if (Feature == "Show All Samples") {
          SampleTypeLab <- "All Features in All Patients\n"
        }
        if (Feature != "Show All Samples") {
          SampleTypeLab <- paste(Feature," in All Patients\n")
        }
      }
      else {
        if (Feature == "Show All Samples") {
          SampleTypeLab <- paste("All Features (",SampleType,") Patients\n",sep = "")
        }
        if (Feature != "Show All Samples") {
          SampleTypeLab <- paste(Feature," (",SampleType,") Patients\n",sep = "")
        }
      }
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (Feature == "Show All Samples") {
        SampleTypeLab <- "All Features in All Patients\n"
      }
      if (Feature != "Show All Samples") {
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
    ggsurv <- survminer::ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
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
    metacol_sampletype <- metacol_sampletype()
    meta <- clin_react()
    decon_score_cols <- decon_score_cols()
    
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
    metacol_sampletype <- metacol_sampletype()
    meta <- clin_react()
    decon_score_cols <- decon_score_cols()
    
    ## Survival Type
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## Determine Feature and sub feature
    if (Feature != "Show All Samples") {
      SubFeature <- input$subFeatureSelection
      Feature <- paste("<b>",Feature,"</b> - <b>",SubFeature,"</b></li>",sep = "")
      line2 <- paste("<li>The dataset is filtered by ",Feature,sep = "")
    }
    if (Feature == "Show All Samples") {
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
      dplyr::select(label,n_obs,estimate,std.error,ci,p.value)
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
    metacol_sampletype <- metacol_sampletype()
    meta <- clin_react()
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
    ggsurv <- survminer::ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
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
    metacol_sampletype <- metacol_sampletype()
    meta <- clin_react()
    
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
    if (Feature != "Show All Samples") {
      SubFeature <- input$subFeatureSelection
      Feature <- paste("<b>",Feature,"</b> - <b>",SubFeature,"</b></li>",sep = "")
      line2 <- paste("<li>The dataset is filtered by ",Feature,sep = "")
    }
    else if (Feature == "Show All Samples") {
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
      forest <- survminer::ggforest(tab,
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
      
      p <- survminer::ggcoxdiagnostics(tab,
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
      }
      if (input$BiVarAddNAcheck2 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature2]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[which(meta_ssgsea[,Feature2] != "Inf"),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature2],ignore.case = T, invert = T),]
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
      dplyr::select(label,n_obs,estimate,std.error,ci,p.value)
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
  #  if (Feature != "Show All Samples") {
  #    SubFeature <- input$subFeatureSelection
  #    Feature <- paste("<b>",Feature,"</b> - <b>",SubFeature,"</b></li>",sep = "")
  #    line2 <- paste("<li>The dataset is filtered by ",Feature,sep = "")
  #  }
  #  if (Feature == "Show All Samples") {
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
      
      
      forest <- survminer::ggforest(tab,
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
      
      p <- survminer::ggcoxdiagnostics(tab,
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
      }
      if (input$BiVarIntNAcheck2 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature2]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[which(meta_ssgsea[,Feature2] != "Inf"),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature2],ignore.case = T, invert = T),]
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
      dplyr::select(label,n_obs,estimate,std.error,ci,p.value)
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
      metacol_sampletype <- metacol_sampletype()
      meta <- clin_react()
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
      ggsurv <- survminer::ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
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
      
      forest <- survminer::ggforest(tab,
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
      
      p <- survminer::ggcoxdiagnostics(tab,
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
        #meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature]) == FALSE),]
        #meta_ssgsea <- meta_ssgsea[which(meta_ssgsea[,Feature] != "Inf"),]
        #meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature],ignore.case = T, invert = T),]hich(colnames(meta_ssgsea) == geneset_name)], quantCutoff2)
        
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
    
    #Feature <- gsubCheck(Feature)
    if (ncol(meta_ssgsea_sdf) >= 3) {
      colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "time"
      colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "ID"
    }
    
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
    
    if (ncol(meta_ssgsea_sdf) >= 3) {
      colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "time"
      colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "ID"
    }
    
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
        dplyr::select(label,n_obs,estimate,std.error,ci,p.value)
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
        dplyr::select(label,n_obs,estimate,std.error,ci,p.value)
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
      
      forest <- survminer::ggforest(tab,
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
  
  ####----Data Exploration----####
  
  ####----Density----####
  
  ssgseaDensity_react <- reactive({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    ssgsea_meta <- ssGSEAmeta()
    user_quant <- input$densityPercent/100
    ShowQuartile <- input$QuartileLinesCheck
    scoreMethod <- input$ScoreMethod
    decon_score_cols <- decon_score_cols()
    
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
  
  ####----Feature Comparison---####
  
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
    decon_score_cols <- decon_score_cols()
    
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
  
  output$FeatCompScatterTable <- DT::renderDataTable({
    
    tab <- FeatCompScatter_react()
    if (input$ColorScatterChoice == "Single Color") {
      tab <- tab[,-4]
    }
    DT::datatable(tab,
                  options = list(scrollY = T),
                  rownames = F)
    
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
    decon_score_cols <- decon_score_cols()
    metacol_sampletype <- metacol_sampletype()
    meta <- clin_react()
    
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
      ggpubr::stat_compare_means(method = input$boxoptselec) +
      theme(text = element_text(size = font),
            legend.position = "none")
    plot
    
  })
  
  output$Sboxplot <- renderPlot({
    
    plot <- Sboxplot_react()
    plot
    
    
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
    
    dataset_og <- expr
    dataset <- expr
    dataset <- log2(dataset + 1)
    zdataset <- apply(dataset, 1, scale)
    zdataset <- apply(zdataset, 1, rev)
    colnames(zdataset) <- colnames(dataset_og)
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
      hmcols <- viridis::viridis(500)
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
    heat <- pheatmap::pheatmap(dataset2,
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
    metacol_sampletype <- metacol_sampletype()
    decon_score_cols <- decon_score_cols()
    meta <- clin_react()
    
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
        ggpubr::stat_compare_means(method = StatMethod) +
        theme(text = element_text(size = font),
                legend.position = "none") +
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
          ggpubr::stat_compare_means(method = StatMethod) +
          theme(text = element_text(size = font),
                legend.position = "none")
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
          ggpubr::stat_compare_means(method = StatMethod) +
          theme(text = element_text(size = font),
                axis.text.x = element_text(angle = as.numeric(boxplotang), hjust = 1),
                legend.position = "none")
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
          ggpubr::stat_compare_means(method = StatMethod) +
          theme(text = element_text(size = font),
                axis.text.x = element_text(angle = as.numeric(boxplotang)),
                legend.position = "none")
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
    
    dataset_og <- expr
    dataset <- expr
    dataset <- as.data.frame(log2(dataset + 1))
    zdataset <- as.data.frame(apply(dataset, 1, scale))
    zdataset <- as.data.frame(apply(zdataset, 1, rev))
    colnames(zdataset) <- colnames(dataset_og)
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
      hmcols <- viridis::viridis(500)
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
    heat <- pheatmap::pheatmap(dataset2,
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
  
  output$rendPurposeAndMethodsMD <- renderUI({
    
    if (file.exists(About_MD_File)) {
      
      includeMarkdown(About_MD_File)
      
    }
    else {
      p("Page requires PurposeAndMethods.Rmd to be found")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuartileSurvival.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuartileSurvival.svg",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuartileSurvival.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuartileSurvival.pdf",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_MedianCutPSurvival.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_MedianCutPSurvival.svg",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_MedianCutPSurvival.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_MedianCutPSurvival.pdf",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_OptimalCutpointSurvival.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_OptimalCutpointSurvival.svg",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_OptimalCutpointSurvival.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_OptimalCutpointSurvival.pdf",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuantileSurvival.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuantileSurvival.svg",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuantileSurvival.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuantileSurvival.pdf",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_AboveBelowCutoffSurvival.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_AboveBelowCutoffSurvival.svg",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_AboveBelowCutoffSurvival.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_AboveBelowCutoffSurvival.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- SquantPlot2_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  
  ####----DNLD Path Density----####
  
  ## quartile
  output$dnldssgseaQuartDensity_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuartileDensity.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuartileDensity.svg",sep = "")
      }
    },
    content = function(file) {
      p <- ssgseaQuartDensity_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  output$dnldssgseaQuartDensity_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuartileDensity.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuartileDensity.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- ssgseaQuartDensity_react()
      ggsave(file,p,width = 10, height = 8)
      
    }
  )
  
  ## binary
  output$dnldssgseaBINDensity_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_MedianCutPDensity.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_MedianCutPDensity.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_OptimalCutpointDensity.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_OptimalCutpointDensity.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuantileDensity.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuantileDensity.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_AboveBelowCutoffDensity.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_AboveBelowCutoffDensity.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateSurvival.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateSurvival.svg",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateSurvival.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateSurvival.pdf",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateForest.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateForest.svg",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateForest.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateForest.pdf",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateLinearity.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateLinearity.svg",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateLinearity.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateLinearity.pdf",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateForest.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateForest.svg",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateForest.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateForest.pdf",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateLinearity.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateLinearity.svg",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateLinearity.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateLinearity.pdf",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateSurvival.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateSurvival.svg",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateSurvival.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateSurvival.pdf",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateForest.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateForest.svg",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateForest.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateForest.pdf",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateLinearity.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateLinearity.svg",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateLinearity.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateLinearity.pdf",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_MultivariateForest.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_MultivariateForest.svg",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_MultivariateForest.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_MultivariateForest.pdf",sep = "")
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
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        SampleTypeLab <- "_"
      }
      if (GeneSet %in% decon_score_cols()) {
        scoreMethod <- "PreProcessed_Score"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethod <- "GeneExpressionScore"
      }
      if (input$GeneSetTabs != 2 & !(GeneSet %in% decon_score_cols())) {
        scoreMethod <- paste(scoreMethod,"Score",sep = "")
      }
      if (logchoice == TRUE) {
        scoreMethod <- paste(scoreMethod,"logPlus1",sep = "")
      }
      paste(gsub(" ","",input$UserProjectName),SampleTypeLab,Feature,"_",GeneSet,"_",scoreMethod,"RiskStratBoxPlot.svg",sep = "")
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
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        SampleTypeLab <- "_"
      }
      
      if (GeneSet %in% decon_score_cols()) {
        scoreMethod <- "PreProcessed_Score"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethod <- "GeneExpressionScore"
      }
      if (input$GeneSetTabs != 2 & !(GeneSet %in% decon_score_cols())) {
        scoreMethod <- paste(scoreMethod,"Score",sep = "")
      }
      if (logchoice == TRUE) {
        scoreMethod <- paste(scoreMethod,"logPlus1",sep = "")
      }
      paste(gsub(" ","",input$UserProjectName),SampleTypeLab,Feature,"_",GeneSet,"_",scoreMethod,"RiskStratBoxPlot.pdf",sep = "")
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
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        SampleTypeLab <- "_"
      }
      if (GeneSet %in% decon_score_cols()) {
        scoreMethod <- "PreProcessed_Score"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethod <- "GeneExpressionScore"
      }
      if (input$GeneSetTabs != 2 & !(GeneSet %in% decon_score_cols())) {
        scoreMethod <- paste(scoreMethod,"Score",sep = "")
      }
      if (logchoice == TRUE) {
        scoreMethod <- paste(scoreMethod,"logPlus1",sep = "")
      }
      paste(gsub(" ","",input$UserProjectName),SampleTypeLab,Feature,"_",GeneSet,"_",scoreMethod,".txt",sep = "")
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
      meta <- clin_react()
      if (length(unique(meta[,metacol_sampletype()])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype()])) <= 1) {
        SampleTypeLab <- "_"
      }
      paste(gsub(" ","",input$UserProjectName),SampleTypeLab,Feature,"_",GeneSet,"_","RiskStratHeatmap.svg",sep = "")
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
      meta <- clin_react()
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype()])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype()])) <= 1) {
        SampleTypeLab <- "_"
      }
      paste(gsub(" ","",input$UserProjectName),SampleTypeLab,Feature,"_",GeneSet,"_","RiskStratHeatmap.pdf",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
        SubFeature <- "_"
      }
      geneset <- gs_react()                 #Chosen Gene Set
      geneset_name <- names(geneset)        #Name of chosen gene set
      # Make file name
      meta <- clin_react()
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype()])) > 1) {
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"SurvivalCutoff_expr.txt",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"SurvivalCutoff_expr.txt",sep = "")
      }
    },
    content = function(file) {
      expr <- exprSub()
      expr <- as.data.frame(expr)
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
      meta <- clin_react()
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype()])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype()])) <= 1) {
        SampleTypeLab <- "_"
      }
      if (GeneSet %in% decon_score_cols()) {
        scoreMethod <- "PreProcessed_Score"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethod <- "GeneExpressionScore"
      }
      if (input$GeneSetTabs != 2 & !(GeneSet %in% decon_score_cols())) {
        scoreMethod <- paste(scoreMethod,"Score",sep = "")
      }
      if (logchoice == TRUE) {
        scoreMethod <- paste(scoreMethod,"logPlus1",sep = "")
      }
      paste(gsub(" ","",input$UserProjectName),SampleTypeLab,Feature,"_",Feature2,"_",GeneSet,"_",scoreMethod,"FeatureStratBoxPlot.svg",sep = "")
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
      meta <- clin_react()
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype()])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype()])) <= 1) {
        SampleTypeLab <- "_"
      }
      
      if (GeneSet %in% decon_score_cols()) {
        scoreMethod <- "PreProcessed_Score"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethod <- "GeneExpressionScore"
      }
      if (input$GeneSetTabs != 2 & !(GeneSet %in% decon_score_cols())) {
        scoreMethod <- paste(scoreMethod,"Score",sep = "")
      }
      if (logchoice == TRUE) {
        scoreMethod <- paste(scoreMethod,"logPlus1",sep = "")
      }
      paste(gsub(" ","",input$UserProjectName),SampleTypeLab,Feature,"_",Feature2,"_",GeneSet,"_",scoreMethod,"FeatureStratBoxPlot.pdf",sep = "")
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
      meta <- clin_react()
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype()])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype()])) <= 1) {
        SampleTypeLab <- "_"
      }
      if (GeneSet %in% decon_score_cols()) {
        scoreMethod <- "PreProcessed_Score"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethod <- "GeneExpressionScore"
      }
      if (input$GeneSetTabs != 2 & !(GeneSet %in% decon_score_cols())) {
        scoreMethod <- paste(scoreMethod,"Score",sep = "")
      }
      if (logchoice == TRUE) {
        scoreMethod <- paste(scoreMethod,"logPlus1",sep = "")
      }
      
      paste(gsub(" ","",input$UserProjectName),SampleTypeLab,Feature,"_",FeatureSelec,"_",GeneSet,"_",scoreMethod,".txt",sep = "")
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
      meta <- clin_react()
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype()])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype()])) <= 1) {
        SampleTypeLab <- "_"
      }
      paste(gsub(" ","",input$UserProjectName),SampleTypeLab,Feature,"_",feature2,"_",GeneSet,"_FeatureHeatmap.svg",sep = "")
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
      meta <- clin_react()
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype()])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype()])) <= 1) {
        SampleTypeLab <- "_"
      }
      paste(gsub(" ","",input$UserProjectName),SampleTypeLab,Feature,"_",feature2,"_",GeneSet,"_FeatureHeatmap.pdf",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
        SubFeature <- "_"
      }
      geneset <- gs_react()                 #Chosen Gene Set
      geneset_name <- names(geneset)        #Name of chosen gene set
      # Make file name
      meta <- clin_react()
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype()])) > 1) {
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"Featuring_",Feature2,"_expr.txt",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"Featuring_",Feature2,"_expr.txt",sep = "")
      }
    },
    content = function(file) {
      expr <- exprSub()
      expr <- as.data.frame(expr)
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      meta <- clin_react()
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype()])) > 1) {
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_DensityPlot.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_DensityPlot.svg",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
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
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      meta <- clin_react()
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype()])) > 1) {
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_DensityPlot.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_DensityPlot.pdf",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
        SubFeature <- "_"
      }
      meta <- clin_react()
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype()])) > 1) {
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",score_method,".txt",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",score_method,".txt",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
        SubFeature <- "_"
      }
      geneset <- gs_react()                 #Chosen Gene Set
      geneset_name <- names(geneset)        #Name of chosen gene set
      # Make file name
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_meta.txt",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",geneset_name,"_meta.txt",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
        SubFeature <- "_"
      }
      # Make file name
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"expr.txt",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"expr.txt",sep = "")
      }
    },
    content = function(file) {
      expr <- exprSub()
      expr <- as.data.frame(expr)
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",FeatureScatter,"_",geneset_name,"_",scoreMethodLab,"_ScatterPlot.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",FeatureScatter,"_",geneset_name,"_",scoreMethodLab,"_ScatterPlot.svg",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",FeatureScatter,"_",geneset_name,"_",scoreMethodLab,"_ScatterPlot.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",FeatureScatter,"_",geneset_name,"_",scoreMethodLab,"_ScatterPlot.pdf",sep = "")
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
      if (Feature != "Show All Samples") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "Show All Samples") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        scoreMethodLab <- "RawGeneExpression"
      }
      else if (input$GeneSetTabs != 2) {
        if (geneset_name %in% decon_score_cols()) {
          scoreMethodLab <- "PreProcessed"
        }
        else {
          scoreMethodLab <- scoreMethod
        }
      }
      # If more than one sample type
      if (length(unique(clin_react()[,metacol_sampletype()])) > 1) {
        paste(gsub(" ","",input$UserProjectName),"_",SampleType,"_",Feature,SubFeature,"_",FeatureScatter,"_",geneset_name,"_",scoreMethodLab,"_ComparisonTable.txt",sep = "")
      }
      # If only one sample type
      else if (length(unique(clin_react()[,metacol_sampletype()])) <= 1) {
        paste(gsub(" ","",input$UserProjectName),"_",Feature,SubFeature,"_",FeatureScatter,"_",geneset_name,"_",scoreMethodLab,"_ComparisonTable.txt",sep = "")
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



shinyApp(ui,server)





















