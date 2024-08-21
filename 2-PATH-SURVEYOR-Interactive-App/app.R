version_id <- paste0("v2.0.20240821")

# User Input -------------------------------------------------------------------

ProjectName <- ''

ExpressionMatrix_file <- ''

MetaData_file <- ''

MetaParam_File <- ''


## Password Settings -----------------------------------------------------------
Password_Protected <- FALSE
PasswordSet <- "teamscience"

## Pre-Selected Inputs ---------------------------------------------------------
# An option from the meta, All, or NULL
PreSelect_SamplyType <- "All"
PreSelect_Feature <- "All"
# An option from the meta or NULL
PreSelect_SubFeature <- NULL
PreSelect_SecondaryFeature <- NULL

## Provided Input --------------------------------------------------------------
## User make sure paths are correct
GeneSet_File <- "Genesets/GeneSet_List_HS_v6.RData"
GeneSetTable_File <- "Genesets/GeneSet_CatTable_v6.txt"
About_MD_File <- "PurposeAndMethods.Rmd"
ExampleExpr_File <- "Example_Data/TCGA_CHOL_Expression_PatientID.txt"
ExampleClin_File <- "Example_Data/TCGA_CHOL_Clinical_PatientID.txt"
ExampleParam_File <- "Example_Data/TCGA_SurvivalApp_Global_Parameters.txt"

## Check if Immune deconvolution package is installed
immudecon <- "immunedeconv"
immudecon_check <- immudecon %in% rownames(installed.packages())


# Input Check ------------------------------------------------------------------

if (!file.exists(ExpressionMatrix_file)) {
  FileProvided <- FALSE
} else { FileProvided <- TRUE }

if (!file.exists(MetaParam_File)) {
  ParamFileProvided <- FALSE
} else { ParamFileProvided <- TRUE }

if (!isTruthy(ProjectName)) {
  ProjectName <- paste("{ Survival Analysis }")
} else {
  ProjectName <- paste("{",ProjectName,"Survival Analysis }")
}
if (isTruthy(PreSelect_Feature)) {
  if (toupper(PreSelect_Feature) == "ALL") {
    PreSelect_Feature <- "Show all Samples"
  }
}

gs <- loadRData(GeneSet_File)
# Gene Set Table
GeneSetTable <- as.data.frame(fread(GeneSetTable_File, sep = '\t', header = T))
GeneSetCats <- unique(GeneSetTable[,1])

# Functions --------------------------------------------------------------------

#increase file upload size
options(shiny.maxRequestSize=5000*1024^2)

SurvPlot_Height <- "550px"
SurvPlot_Width <- "850px"

StatCols <- c("MedianCutP","QuartileCutP","OptimalCutP","TopBottomCutP","UserCutP")


# Password Table ---------------------------------------------------------------
# user database for login
if (Password_Protected) {
  user_base <- tibble::tibble(
    user = "user",
    password = PasswordSet,
    permissions = "admin",
    name = "User"
  )
}

# UI Tabs ----------------------------------------------------------------------
## Login Tab -------------------------------------------------------------------
login_tab <- tabPanel(
  title = icon("lock"),
  value = "login",
  loginUI("login")
)

## Data Input Tab --------------------------------------------------------------
#if (!FileProvided) {
DataInput_tab <- tabPanel("Data Input",
                          fluidPage(
                            sidebarPanel(
                              width = 3,
                              id = "DataInputPanel",
                              p(),
                              textInput("UserProjectName","Project Name:", value = "Survival Analysis"),
                              uiOutput("rendExprFileInput"),
                              fluidRow(
                                column(5, style = 'padding-right:2px;margin-top:-30px;',
                                       checkboxInput("LogExprFile","Log2 Expression", value = F)
                                ),
                                column(7, style = 'padding-left:2px;margin-top:-30px;',
                                       checkboxInput("ScaleNormExprFile","Scale Normalize Expression", value = F)
                                )
                              ),
                              uiOutput("rendClinFileInput"),
                              fluidRow(
                                column(12, style = 'margin-top:-15px;',
                                       actionButton("UseExpData","Load Example Data"),
                                       tags$a(href="http://shawlab.science/shiny/PATH_SURVEYOR_ExampleData/PATH_SURVEYOR_App/", "Download example data", target='_blank'),
                                )
                              ),
                              shiny::hr(),
                              uiOutput("rendClinParamHeader"),
                              h4("Clinical Parameters"),
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
                              )
                            ),
                            mainPanel(
                              uiOutput("rendExprFilePrevHeader"),
                              div(DT::dataTableOutput("ExprFile_Preview"), style = "font-size:10px"),
                              uiOutput("rendClinFilePrevHeader"),
                              div(DT::dataTableOutput("ClinFile_Preview"), style = "font-size:10px"),
                              uiOutput("rendParamFilePrevHeader"),
                              div(DT::dataTableOutput("ClinParamFile_Preview"), style = "font-size:10px")
                            ),
                            tagList(
                              tags$head(
                                tags$style(
                                  HTML("
                                     .info_box {
                                     width: auto;
                                     height: auto;
                                     color: #000000;
                                     background-color: #f5f5f5;
                                     padding: 3px 8px;
                                     font-size: 12px;
                                     z-index : 9999;
                                     }",
                                     glue::glue("#{'AppVersion'} {{
                                                position: {'fixed'};
                                                top: 0;
                                                right: 0;
                                                }}")
                                  )
                                )
                              ),
                              div(id = "AppVersion", class = "info_box", version_id)
                            )
                          )
)
#}

## Survival Tab ----------------------------------------------------------------

Survival_tab <- tabPanel("Survival Analysis",
                         fluidPage(
                           sidebarLayout(
                             
                             ### Sidebar Panel ---------------------------------
                             
                             sidebarPanel(
                               conditionalPanel(condition = "input.SurvPanels == '1' | input.SurvPanels == '2' | input.SurvPanels == '3' | input.SurvPanels == '4'",
                                                tabsetPanel(
                                                  id = "survside",
                                                  tabPanel("Sample Parameters",
                                                           p(),
                                                           uiOutput("rendSampleTypeSelection"),
                                                           selectizeInput("FeatureSelection","Select Feature to Subset Samples:",choices = NULL, selected = 1),
                                                           uiOutput("rendSubFeatureSelection"),
                                                           fluidRow(
                                                             column(4,
                                                                    selectizeInput("SurvivalType_time","Survival Time Data:",choices = NULL, selected = 1),
                                                             ),
                                                             column(4,
                                                                    selectizeInput("SurvivalType_id","Survival ID Data:",choices = NULL, selected = 1),
                                                             ),
                                                             column(4,
                                                                    selectInput("ScoreMethod","Scoring Method",choices = c("ssgsea","gsva","zscore","plage"))
                                                             )
                                                           ),
                                                           tabsetPanel(
                                                             id = "GeneSetTabs",
                                                             tabPanel("Gene Sets",
                                                                      p(),
                                                                      selectInput("GeneSetCat_Select","Select Geneset Category:", choices = GeneSetCats),
                                                                      div(DT::dataTableOutput("GeneSetTable"), style = "font-size:10px"),
                                                                      value = 1
                                                             ),
                                                             tabPanel("Single Genes",
                                                                      p(),
                                                                      div(DT::dataTableOutput("GeneGeneSetTable"), style = "font-size:10px"),
                                                                      value = 2
                                                             ),
                                                             tabPanel("User Gene Set",
                                                                      p(),
                                                                      radioButtons("UserGSoption","",choices = c("Gene Set Upload","Text Box Input"), inline = T),
                                                                      conditionalPanel(condition = "input.UserGSoption == 'Gene Set Upload'",
                                                                                       fileInput("userGeneSet","Gene Set Upload", accept = c(".gmt",".tsv",".txt",".RData")),
                                                                                       div(DT::dataTableOutput("userGeneSetTable"), style = "font-size:10px")
                                                                      ),
                                                                      conditionalPanel(condition = "input.UserGSoption == 'Text Box Input'",
                                                                                       textInput("userGeneSetTextName","Custom Gene Set Name", value = "Custom_Geneset"),
                                                                                       textInput("userGeneSetText","Gene Symbols", placeholder = "Comma, space, or tab delimited")
                                                                      ),
                                                                      value = 3
                                                             )
                                                           ),
                                                           checkboxInput("ViewGeneSetGenes","View Genes in Selected Gene Set", value = F),
                                                           uiOutput("rendGenesInGeneSetTab")
                                                  ),
                                                  tabPanel("Figure Parameters",
                                                           p(),
                                                           conditionalPanel(condition = "input.SurvPanels == '1' | input.SurvPanels == '2' | input.SurvPanels == '3'",
                                                                            h3("Survival Plot Parameters"),
                                                                            fluidRow(
                                                                              column(6,
                                                                                     numericInput("SurvXaxis","Survival Age Limit (years)", value = NA),
                                                                                     #uiOutput("rendSurvXaxis"),
                                                                                     numericInput("SurvXaxisBreaks","Survival X-Axis Breaks (Years):",value = 1, min = 0, step = 0.25),
                                                                                     selectInput("SurvLegendPos","Legend Position",choices = c("right","left","top","bottom","none"))
                                                                              ),
                                                                              column(6, style = "margin-top:15px",
                                                                                     radioButtons("SurvYearOrMonth","Survival X-Axis Units:",choices = c("Years","Months"), inline = T),
                                                                                     checkboxInput("ShowPval","Show P.Value",value = T),
                                                                                     checkboxInput("ShowConfInt","Show Confidence Interval",value = F),
                                                                                     checkboxInput("ShowMedSurvLine","Show Median Survival Line",value = F)
                                                                              )
                                                                            )
                                                           ),
                                                           conditionalPanel(condition = "input.SurvPanels == '2' | input.SurvPanels == '3'",
                                                                            shiny::hr(),
                                                                            h3("Forest Plot Parameters"),
                                                                            numericInput("ForestFontSize","Font Size",value = 1),
                                                                            shiny::hr(),
                                                                            h3("Linearity Plot Parameters"),
                                                                            fluidRow(
                                                                              column(4,
                                                                                     numericInput("linAxisFont","X/Y Axis Font Size",value = 14, step = 1)
                                                                              ),
                                                                              column(4,
                                                                                     numericInput("linTickFont","Axis Tick Font Size",value = 10, step = 1)
                                                                              ),
                                                                              column(4,
                                                                                     numericInput("linMainFont","Title Font Size",value = 16, step = 1)
                                                                              )
                                                                            )
                                                           ),
                                                           conditionalPanel(condition = "input.SurvPanels == '4'",
                                                                            p(),
                                                                            h3("Boxplot Parameters"),
                                                                            fluidRow(
                                                                              column(4,
                                                                                     numericInput("boxplotFont","Font Size:", value = 15, step = 1)
                                                                              ),
                                                                              column(4,
                                                                                     numericInput("boxplotDot", "Dot Size:", value = 0.75, step = 0.25)
                                                                              ),
                                                                              column(4,
                                                                                     selectInput("boxplotTextAngle","X-Axis Text:",
                                                                                                 choices = c("Horizontal (0 degrees)" = 0,"Angled (45 degrees)" = 45,"Vertical (90 degrees)" = 90))
                                                                              )
                                                                            ),
                                                                            hr(),
                                                                            h3("Heatmap Parameters"),
                                                                            fluidRow(
                                                                              column(6,
                                                                                     numericInput("heatmapFontR", "Row Font Size:", value = 9, step = 1)
                                                                              ),
                                                                              column(6,
                                                                                     numericInput("heatmapFontC", "Column Font Size:", value = 10, step = 1)
                                                                              )
                                                                            )
                                                           ),
                                                           hr(),
                                                           h3("Plot Download Parameters"),
                                                           fluidRow(
                                                             column(4,
                                                                    numericInput("PlotDnldHight","Plot Height",value = 8, min = 0, step = 1)
                                                             ),
                                                             column(4,
                                                                    numericInput("PlotDnldWidth","Plot Width",value = 8, min = 0, step = 1)
                                                             ),
                                                             column(4,
                                                                    selectInput("PlotDnldUnits","Units",choices = c("in","cm","mm","px"))
                                                             )
                                                           )
                                                           
                                                  )
                                                )
                               )
                               #conditionalPanel(condition = "input.SurvPanels == '5'",
                               #                 tabsetPanel(
                               #                   tabPanel("Input Parameters",
                               #                            p(),
                               #                            h4("Sample Selection"),
                               #                            uiOutput("rendSampleTypeSelection_lasso"),
                               #                            selectizeInput("FeatureSelection_lasso","Select Feature:", choices = NULL, selected = 1),
                               #                            uiOutput("rendSubFeatureSelection_lasso"),
                               #                            h4("Feature Selection"),
                               #                            textInput("LassoModelName","Lasso Model Name:",value = "Custom_Lasso_Model"),
                               #                            selectizeInput("LassoFeatureSelection_lasso","Select or Paste Features to Generate Lasso Model:", choices = NULL, selected = 1, multiple = T),
                               #                            h4("Lasso Parameters"),
                               #                            fluidRow(
                               #                              column(6,
                               #                                     selectizeInput("SurvTimeSelec_lasso","Survival Time Data:", choices = NULL, selected = 1),
                               #                              ),
                               #                              column(6,
                               #                                     selectizeInput("SurvIDSelect_lasso","Survival ID Data:", choices = NULL, selected = 1),
                               #                              )
                               #                            ),
                               #                            fluidRow(
                               #                              column(5,
                               #                                     numericInput("LassoTrainProp","Training Sample Proportion:",value = 50, step = 1, max = 100,min = 1)
                               #                              ),
                               #                              column(4,
                               #                                     numericInput("LassoAlpha","Set Lasso Alpha:", value = 1, min = 0, step = 0.1)
                               #                              ),
                               #                              column(3,
                               #                                     numericInput("LassoSeedSelection","Set Seed:", value = 100, min = 1, step = 1)
                               #                              )
                               #                            ),
                               #                            p(),
                               #                            fluidRow(
                               #                              column(8,
                               #                                     radioButtons("viewLassoMinOrSE","Select Lambda to Generate Risk Score:", inline = T,choices = c("Lambda Min","Lambda SE","Custom")),
                               #                              ),
                               #                              column(4,
                               #                                     uiOutput("rendCustomLambda")
                               #                              )
                               #                            ),
                               #                            uiOutput("rednLassoCoefTable"),
                               #                            #,
                               #                            fluidRow(
                               #                              column(6,
                               #                                     dnld_ui("dnldLassoModel","Download Lasso Model")
                               #                              ),
                               #                              column(6,
                               #                                     dnld_ui("dnldLassoRunData","Download Lasso Run Data")
                               #                              )
                               #                            )
                               #                   ),
                               #                   tabPanel("Figure Parameters",
                               #                            p(),
                               #                            h4("Survival Plot Parameters"),
                               #                            radioButtons("LassoPlotCutP","Plot Cut-Point", choices = c("Median","Quartile","Optimal","Quantile","User Specified"),
                               #                                         inline = T),
                               #                            uiOutput("rendCutPinput"),
                               #                            fluidRow(
                               #                              column(4,
                               #                                     uiOutput("rendSurvXaxis_lasso")
                               #                              ),
                               #                              column(8,
                               #                                     textInput("SurvPlotTitle_lasso","Survival Plot Title:",value = "")
                               #                              )
                               #                            ),
                               #                            fluidRow(
                               #                              column(3,
                               #                                     selectInput("SurvLegendPos_lasso","Legend Position",choices = c("top","right","left","bottom","none"))
                               #                              ),
                               #                              column(3,
                               #                                     checkboxInput("ShowPval_lasso","Show P.Value",value = T)
                               #                              ),
                               #                              column(3,
                               #                                     checkboxInput("ShowConfInt_lasso","Show Confidence Interval",value = F)
                               #                              ),
                               #                              column(3,
                               #                                     checkboxInput("ShowMedSurvLine_lasso","Show Median Survival Line",value = F)
                               #                              )
                               #                            )
                               #                   )
                               #                 )
                               #)
                             ),
                             
                             ### Main Panel ------------------------------------
                             
                             mainPanel(
                               tabsetPanel(
                                 id = "SurvPanels",
                                 
                                 #### Main Surv --------------------------------
                                 
                                 tabPanel("Pathway Level Survival Analysis",
                                          tabsetPanel(
                                            id = "SurvPanelsMain",
                                            tabPanel("Median Cut-Point Survival",
                                                     p(),
                                                     tableOutput("ssgseaCheck"),
                                                     htmlOutput("BINSurvDescrip", style = "font-size:14px;"),
                                                     shiny::hr(),
                                                     shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("SplotBIN", height = SurvPlot_Height, width = SurvPlot_Width)), type = 6),
                                                     dnld_ui("dnldSplotBIN_SVG","SVG"),
                                                     shiny::hr(),
                                                     shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("ssgseaBINDensity", width = "650px", height = "400px")), type = 6),
                                                     dnld_ui("dnldssgseaBINDensity_SVG","SVG"),
                                                     shiny::hr(),
                                                     h4("Cox Hazard Regression Analysis Summary"),
                                                     fluidRow(
                                                       column(6,
                                                              div(shinycssloaders::withSpinner(tableOutput("SBinaryHRtab"), type = 7, size = 0.5), style = "font-size:12px"),
                                                              style = 'border-right: 0.5px solid lightgray',
                                                       ),
                                                       column(6,
                                                              verbatimTextOutput("MedianCutPSumm")
                                                       )
                                                     ),
                                                     value = 1),
                                            tabPanel("Quartile Survival",
                                                     p(),
                                                     htmlOutput("QuartSurvDescrip", style = "font-size:14px;"),
                                                     shiny::hr(),
                                                     shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("Splot", height = SurvPlot_Height, width = SurvPlot_Width)), type = 6),
                                                     dnld_ui("dnldSplot_SVG","SVG"),
                                                     shiny::hr(),
                                                     shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("ssgseaQuartDensity", width = "650px", height = "400px")), type = 6),
                                                     dnld_ui("dnldssgseaQuartDensity_SVG","SVG"),
                                                     shiny::hr(),
                                                     h4("Cox Hazard Regression Analysis Summary"),
                                                     fluidRow(
                                                       column(6,
                                                              div(shinycssloaders::withSpinner(tableOutput("SQuartileHRtab"), type = 7, size = 0.5), style = "font-size:12px"),
                                                              style = 'border-right: 0.5px solid lightgray',
                                                       ),
                                                       column(6,
                                                              verbatimTextOutput("QuartileCutPSumm")
                                                       )
                                                     ),
                                                     value = 2),
                                            tabPanel("Optimal Cut-Point Survival",
                                                     p(),
                                                     htmlOutput("CutPSurvDescrip", style = "font-size:14px;"),
                                                     shiny::hr(),
                                                     shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("ScutPointPlot", height = SurvPlot_Height, width = SurvPlot_Width)), type = 6),
                                                     dnld_ui("dnldScutPointPlot_SVG","SVG"),
                                                     shiny::hr(),
                                                     shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("ssgseaCutPDensity", width = "650px", height = "400px")), type = 6),
                                                     dnld_ui("dnldssgseaCutPDensity_SVG","SVG"),
                                                     shiny::hr(),
                                                     h4("Cox Hazard Regression Analysis Summary"),
                                                     fluidRow(
                                                       column(6,
                                                              div(shinycssloaders::withSpinner(tableOutput("CutPointHRtab"), type = 7, size = 0.5), style = "font-size:12px"),
                                                              #uiOutput("rendCutPointHRtab"),
                                                              style = 'border-right: 0.5px solid lightgray',
                                                       ),
                                                       column(6,
                                                              verbatimTextOutput("OptimalCutPSumm")
                                                       )
                                                     ),
                                                     value = 3),
                                            tabPanel("Top/Bottom Cut-Point Survival",
                                                     p(),
                                                     htmlOutput("QuantSurvDescrip", style = "font-size:14px;"),
                                                     shiny::hr(),
                                                     shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("SquantPlot", height = SurvPlot_Height, width = SurvPlot_Width)), type = 6),
                                                     numericInput("QuantPercent","Top/Bottom Cut-Point Quantile Cutoff (%)", value = 25, min = 0, max = 100),
                                                     dnld_ui("dnldSquantPlot_SVG","SVG"),
                                                     shiny::hr(),
                                                     shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("ssgseaQuantDensity", width = "650px", height = "400px")), type = 6),
                                                     dnld_ui("dnldssgseaQuantDensity_SVG","SVG"),
                                                     shiny::hr(),
                                                     h4("Cox Hazard Regression Analysis Summary"),
                                                     fluidRow(
                                                       column(6,
                                                              div(shinycssloaders::withSpinner(tableOutput("SQuantileHRtab"), type = 7, size = 0.5), style = "font-size:12px"),
                                                              #uiOutput("rendQuantHRtab"),
                                                              style = 'border-right: 0.5px solid lightgray',
                                                       ),
                                                       column(6,
                                                              verbatimTextOutput("QuantileCutPSumm")
                                                       )
                                                     ),
                                                     value = 4),
                                            tabPanel("User Cut-Point Survival",
                                                     p(),
                                                     htmlOutput("Quant2SurvDescrip", style = "font-size:14px;"),
                                                     shiny::hr(),
                                                     shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("SquantPlot2", height = SurvPlot_Height, width = SurvPlot_Width)), type = 6),
                                                     numericInput("QuantPercent2","Above/Below User Quantile Cut-Point (%)", value = 25, min = 0, max = 100),
                                                     dnld_ui("dnldSquantPlot2_SVG","SVG"),
                                                     shiny::hr(),
                                                     shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("ssgseaQuant2Density", width = "650px", height = "400px")), type = 6),
                                                     dnld_ui("dnldssgseaQuant2Density_SVG","SVG"),
                                                     shiny::hr(),
                                                     h4("Cox Hazard Regression Analysis Summary"),
                                                     fluidRow(
                                                       column(6,
                                                              div(shinycssloaders::withSpinner(tableOutput("SQuantileHR2tab"), type = 7, size = 0.5), style = "font-size:12px"),
                                                              style = 'border-right: 0.5px solid lightgray',
                                                       ),
                                                       column(6,
                                                              verbatimTextOutput("UserCutPSumm")
                                                       )
                                                     ),
                                                     value = 5)
                                          ),
                                          value = 1),
                                 
                                 #### Univar Surv ------------------------------
                                 
                                 tabPanel("Univariate Survival Analysis",
                                          p(),
                                          fluidRow(
                                            column(6,
                                                   selectizeInput("SingleSurvivalFeature","Select Feature:",choices = NULL, selected = 1, width = "80%"),
                                                   # Allows all select inputs to be wide enough to read the contents
                                                   tags$head(
                                                     tags$style(HTML('
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
                                                            conditionalPanel(condition = "input.UniVarContCheck == true",
                                                                             checkboxInput("UniVarContHiLoCheck","Continuous Feature as Median Cut-Point",value = T)
                                                            )
                                                            
                                                     )
                                                   ),
                                                   uiOutput("rendSurvFeatVariableUni"),
                                            )
                                          ),
                                          tabsetPanel(
                                            id = "UniVarPlots",
                                            tabPanel("Survival Plot",
                                                     p(),
                                                     shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("featSplot", width = SurvPlot_Width, height = SurvPlot_Height)), type = 6),
                                                     dnld_ui("dnldfeatSplot_SVG","SVG")
                                            ),
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
                                            tabPanel("Forest Plot",
                                                     p(),
                                                     shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("SinglevarForestPlot", width = "100%", height = "800px")), type = 6),
                                                     downloadButton("dnldUniVarForestplot_SVG","SVG")
                                            ),
                                            tabPanel("Multi-Feature Forest Plot",
                                                     p(),
                                                     fluidRow(
                                                       column(6,
                                                              selectizeInput("MultiFeatUnivarSelect","Select Features for Forest Plot",
                                                                             choices = NULL, multiple = T, selected = 1, width = "100%")
                                                       ),
                                                       column(2,
                                                              textInput("UnivarForestPlotXlim","X-Limits (Low,High)",value = "0,5",
                                                                        placeholder = "Low,High", width = "150px")
                                                       ),
                                                       column(2,
                                                              selectInput("UnivarForestPlotXtrans","X-Axis Transform",
                                                                          choices = c("none", "log", "log2", "log10"))
                                                       )
                                                     ),
                                                     uiOutput("rendMultiFeatUnivarForestPlot"),
                                                     fluidRow(
                                                       dnld_ui("dnldMultiFeatUnivarForestPlot_SVG","SVG"),
                                                       dnld_ui("dnldMultiFeatUnivarForestPlot_table","Download as Table")
                                                     ),
                                                     p(),
                                                     div(DT::dataTableOutput("univarForestPlotTable"), style = "font-size:12px"),
                                                     dnld_ui("dnldunivarForestPlotTable","Download Table")
                                                     
                                            ),
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
                                                     dnld_ui("dnldUniVarLinplot_SVG","SVG")
                                            )
                                          ),
                                          value = 2),
                                 
                                 #### Multivar Surv ----------------------------
                                 
                                 tabPanel("Multivariate Coxh Analysis",
                                          tabsetPanel(
                                            id = "multivariate",
                                            
                                            ##### Multivar Add Surv ----------------------------
                                            
                                            tabPanel("Bivariate Additive Survival Analysis",
                                                     p(),
                                                     fluidRow(
                                                       column(6,
                                                              selectizeInput("SurvivalFeatureBi1","Select Feature 1:",choices = NULL, selected = 1, width = "80%"),
                                                              #uiOutput("rendSurvivalFeatureBi1"),
                                                              fluidRow(
                                                                column(4,
                                                                       checkboxInput("BiVarAddNAcheck1","Remove NA/Unknown/Inf",value = T)
                                                                ),
                                                                column(4,
                                                                       checkboxInput("BiVarAddContCheck1","Continuous Feature",value = F)
                                                                ),
                                                                column(4,
                                                                       conditionalPanel(condition = "input.BiVarAddContCheck1 == true",
                                                                                        checkboxInput("BiVarAddContHiLoCheck1","Continuous Feature as Median Cut-Point",value = T)
                                                                       )
                                                                       #uiOutput("rendBiVarAddContHiLoCheck1")
                                                                )
                                                              ),
                                                              uiOutput("rendSurvFeatVariableBi1")
                                                       ),
                                                       column(6,
                                                              selectizeInput("SurvivalFeatureBi2","Select Feature 2:",choices = NULL, selected = 1, width = "80%"),
                                                              #uiOutput("rendSurvivalFeatureBi2"),
                                                              fluidRow(
                                                                column(4,
                                                                       checkboxInput("BiVarAddNAcheck2","Remove NA/Unknown/Inf",value = T)
                                                                ),
                                                                column(4,
                                                                       checkboxInput("BiVarAddContCheck2","Continuous Feature",value = F)
                                                                ),
                                                                column(4,
                                                                       conditionalPanel(condition = "input.BiVarAddContCheck2 == true",
                                                                                        checkboxInput("BiVarAddContHiLoCheck2","Continuous Feature as Median Cut-Point",value = T)
                                                                       )
                                                                       #uiOutput("rendBiVarAddContHiLoCheck2")
                                                                )
                                                              ),
                                                              uiOutput("rendSurvFeatVariableBi2")
                                                       ),
                                                       #column(4,
                                                       #       htmlOutput("BivarAddSummExpl", style = "font-size:14px;")
                                                       #)
                                                     ),
                                                     tabsetPanel(
                                                       id = "BiVarPlots",
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
                                                       tabPanel("Forest Plot",
                                                                p(),
                                                                shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("BivarForestPlot", width = "100%", height = "800px")), type = 6),
                                                                downloadButton("dnldBiVarAddForest_SVG","SVG")
                                                                #dnld_ui("dnldBiVarAddForest_SVG","SVG")
                                                       ),
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
                                                                dnld_ui("dnldBiVarAddLinplot_SVG","SVG")
                                                       )
                                                     )
                                            ),
                                            
                                            ##### Multivar Int Surv ----------------------------
                                            
                                            tabPanel("Bivariate Interaction Survival Analysis",
                                                     p(),
                                                     fluidRow(
                                                       column(6,
                                                              selectizeInput("SurvivalFeatureBi1Inter","Select Feature 1:",choices = NULL, selected = 1, width = "80%"),
                                                              fluidRow(
                                                                column(4,
                                                                       checkboxInput("BiVarIntNAcheck1","Remove NA/Unknown/Inf",value = T)
                                                                ),
                                                                column(4,
                                                                       checkboxInput("BiVarIntContCheck1","Continuous Feature",value = F)
                                                                ),
                                                                column(4,
                                                                       conditionalPanel(condition = "input.BiVarIntContCheck1 == true",
                                                                                        checkboxInput("BiVarIntContHiLoCheck1","Continuous Feature as Median Cut-Point",value = T)
                                                                       )
                                                                )
                                                              ),
                                                              uiOutput("rendSurvFeatVariableBi1Inter")
                                                              
                                                       ),
                                                       column(6,
                                                              selectizeInput("SurvivalFeatureBi2Inter","Select Feature 1:",choices = NULL, selected = 1, width = "80%"),
                                                              fluidRow(
                                                                column(4,
                                                                       checkboxInput("BiVarIntNAcheck2","Remove NA/Unknown/Inf",value = T)
                                                                ),
                                                                column(4,
                                                                       checkboxInput("BiVarIntContCheck2","Continuous Feature",value = F)
                                                                ),
                                                                column(4,
                                                                       conditionalPanel(condition = "input.BiVarIntContCheck2 == true",
                                                                                        checkboxInput("BiVarIntContHiLoCheck2","Continuous Feature as Median Cut-Point",value = T)
                                                                       )
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
                                                       tabPanel("Survival Plot",
                                                                p(),
                                                                shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("featSplotBi", width = SurvPlot_Width, height = SurvPlot_Height)), type = 6),
                                                                dnld_ui("dnldfeatSplotBi_SVG","SVG")
                                                       ),
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
                                                                dnld_ui("dnldBiVarIntLinplot_SVG","SVG")
                                                       )
                                                       
                                                     )
                                            ),
                                            
                                            ##### Multivar Surv ----------------------------
                                            
                                            tabPanel("Multivariate Coxh Analysis",
                                                     p(),
                                                     fluidRow(
                                                       column(6,
                                                              selectizeInput("SurvivalFeature","Select Feature(s):",choices = NULL, selected = 1, multiple = T, width = "100%"),
                                                              #uiOutput("rendSurvivalFeature"),
                                                              checkboxInput("MultiVarNAcheck","Remove NA/Unknown",value = T)
                                                       )
                                                     ),
                                                     tabsetPanel(
                                                       id = "multivartabstwo",
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
                                                       tabPanel("Forest Plot",
                                                                p(),
                                                                shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("MultivarForestPlot", width = "100%", height = "800px")), type = 6),
                                                                downloadButton("dnldMultiVarForest_SVG","SVG")
                                                       )
                                                     )
                                            ),
                                            ##### Multivar Forest -----------------
                                            tabPanel("Multivariate Forest Plot",
                                                     p(),
                                                     fluidRow(
                                                       column(4,
                                                              selectizeInput("MultiFeatMultivarSelect","Select Main Forest Feature",
                                                                             choices = NULL, multiple = F, selected = 1, width = "100%"),
                                                              #uiOutput("rendMultiFeatMultivarSelect"),
                                                              uiOutput("rendMultiFeatMultivarRefSelect")
                                                       ),
                                                       column(3,
                                                              selectizeInput("MultiFeatMultivarSubSelect","Select and Add Features to Forest Plot:",
                                                                             choices = NULL, multiple = T, selected = 1, width = "100%")
                                                       ),
                                                       column(2,
                                                              textInput("MultivarForestPlotXlim","X-Limits (Low,High)",value = "0,5",
                                                                        placeholder = "Low,High")
                                                       ),
                                                       column(2,
                                                              selectInput("MultivarForestPlotXtrans","X-Axis Transform",
                                                                          choices = c("none", "log", "log2", "log10"))
                                                       )
                                                     ),
                                                     uiOutput("rendMultiFeatMultivarForestPlot"),
                                                     fluidRow(
                                                       dnld_ui("dnldMultiFeatMultivarForestPlot_SVG","SVG"),
                                                       dnld_ui("dnldMultiFeatMultivarForestPlot_table","Download as Table")
                                                     ),
                                                     p(),
                                                     div(DT::dataTableOutput("multiForestPlotTable"), style = "font-size:12px"),
                                                     dnld_ui("dnldmultiForestPlotTable","Download Table")
                                                     
                                            )
                                          ),
                                          value = 3),
                                 
                                 #### Lasso
                                 
                                 #tabPanel("Lasso Cox Model",
                                 #         p(),
                                 #         fluidRow(
                                 #           column(2,
                                 #                  br(),
                                 #                  actionButton("RunLassoModelGen","Generate Lasso Model")
                                 #                  ),
                                 #           column(3,
                                 #                  uiOutput("rendSurvTimeSelecView_lasso")
                                 #                  ),
                                 #           column(3,
                                 #                  uiOutput("rendSurvIDSelecView_lasso")
                                 #                  )
                                 #         ),
                                 #         fluidRow(
                                 #           column(6,
                                 #                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("Lasso_Train_Splot", width = "100%", height = "500px")),type = 6),
                                 #                  dnld_ui("dnldSplotLassoTrain_SVG","SVG"),
                                 #                  shiny::hr(),
                                 #                  h4("Path of Coefficients"),
                                 #                  shinyjqui::jqui_resizable(plotOutput("Lasso_CoeffPlot", width = "100%", height = "400px")),
                                 #                  shiny::hr(),
                                 #                  h4("Training Cox Hazard Regression Analysis Summary"),
                                 #                  uiOutput("rendLassoTrainHRtab"),
                                 #                  verbatimTextOutput("LassoTrainCoxSumm")
                                 #           ),
                                 #           column(6,
                                 #                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("Lasso_Test_Splot", width = "100%", height = "500px")),type = 6),
                                 #                  dnld_ui("dnldSplotLassoTest_SVG","SVG"),
                                 #                  shiny::hr(),
                                 #                  h4("Lambda Cross-Validation"),
                                 #                  shinyjqui::jqui_resizable(plotOutput("Lasso_LambdaPlot", width = "100%", height = "400px")),
                                 #                  shiny::hr(),
                                 #                  h4("Testing Cox Hazard Regression Analysis Summary"),
                                 #                  uiOutput("rendLassoTestHRtab"),
                                 #                  verbatimTextOutput("LassoTestCoxSumm")
                                 #           )
                                 #           ),
                                 #         value = 5),
                                 
                                 #### Data Explor ------------------------------
                                 
                                 tabPanel("Data Exploration",
                                          tabsetPanel(
                                            id = "DataExploration",
                                            ##### Surv Data ---------------------------------
                                            tabPanel("Download Survival Data",
                                                     p(),
                                                     selectizeInput("MetaTableCols","Select Meta Columns to View:",multiple = T, choices = NULL, selected = ""),
                                                     div(DT::dataTableOutput("MetaTable"), style = "font-size:12px"),
                                                     fluidRow(
                                                       column(2,
                                                              dnld_ui("DnldClin","Download Clinical Data")
                                                       ),
                                                       column(2,
                                                              dnld_ui("DnldExpr","Download Expression Data")
                                                       )
                                                     ),
                                                     value = 5),
                                            ##### Density ---------------------------------
                                            tabPanel("Score Density",
                                                     p(),
                                                     fluidRow(
                                                       numericInput("densityPercent","User Defined Percentile (Red)",value = 15, width = "200px"),
                                                       checkboxInput("QuartileLinesCheck","Show Quartile Lines (Blue)",value = T)
                                                     ),
                                                     shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("ssgseaDensity", width = "100%", height = "500px")), type = 6),
                                                     dnld_ui("dnldssgseaDensity_SVG","SVG"),
                                                     p(),
                                                     div(DT::dataTableOutput("ssgseaDensityTable"), style = "font-size:12px"),
                                                     dnld_ui("dnldssgseaDensityTable","Download Table"),
                                                     value = 6),
                                            ##### Feat Comp ---------------------------------
                                            tabPanel("Feature Comparison",
                                                     p(),
                                                     fluidRow(
                                                       column(3,
                                                              selectizeInput("ScatterFeature","Select Feature:",multiple = F, choices = NULL, selected = 1)
                                                       ),
                                                       column(2,
                                                              radioButtons("ColorScatterChoice","Color Plot by:",choices = c("Feature","Single Color"))
                                                       ),
                                                       column(3,
                                                              conditionalPanel(condition = "input.ColorScatterChoice == 'Feature'",
                                                                               selectizeInput("ScatterColor","Color By:",multiple = F, choices = NULL, selected = 1)
                                                              ),
                                                              conditionalPanel(condition = "input.ColorScatterChoice == 'Single Color'",
                                                                               selectInput("ScatterColor2","Color:",choices = colors(), selected = "cadetblue"))
                                                       ),
                                                       column(1, style = "padding-right:2px",
                                                              checkboxGroupInput("ScatterLog","", choices = c("Log x-axis","Log y-axis"))
                                                       ),
                                                       column(2, style = "margin-top:12px",
                                                              checkboxInput("RegressionLine","Add Regression Line", value = T)
                                                       )
                                                     ),
                                                     shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotlyOutput("FeatCompScatterPlot", width = "100%", height = "500px")), type = 6),
                                                     p(),
                                                     div(DT::dataTableOutput("FeatCompScatterTable"), style = "font-size:12px"),
                                                     dnld_ui("dnldFeatCompScatterTable","Download Table"),
                                                     value = 7),
                                            ##### Risk Strat ---------------------------------
                                            tabPanel("Risk Stratification",
                                                     p(),
                                                     fluidRow(
                                                       column(4,
                                                              numericInput("cutoffTime1","High-Risk Survival Time Cutoff:", value = 364, min = 0, step = 1),
                                                              selectInput("survStatus1","Survival Status Below Cutoff:", choices = c("1","0","0/1"), selected = "1")
                                                       ),
                                                       column(4,
                                                              numericInput("cutoffTime0","Low-Risk Survival Time Cutoff:", value = 365, min = 0, step = 1),
                                                              selectInput("survStatus0","Survival Status Above Cutoff:", choices = c("1","0","0/1"), selected = "0")
                                                       ),
                                                       column(4,
                                                              conditionalPanel(condition = "input.riskstrat == 'box'",
                                                                               selectInput("boxoptselecRisk","Stat Compare Method:", choices = c("none","wilcox.test","t.test")),
                                                                               checkboxInput("SBoxLog", "Log Transform Score", value = T)
                                                              ),
                                                              conditionalPanel(condition = "input.riskstrat == 'heat'",
                                                                               selectizeInput("riskHeatAnno","Add Annotation:", choices = NULL, selected = 1, multiple = T),
                                                                               selectInput("ClusterMethod", "Select Cluster Method:",
                                                                                           choices = c("ward.D", "ward.D2", "complete", "single", "average", "mcquitty", "median", "centroid"))
                                                              )
                                                              
                                                       )
                                                     ),
                                                     tabsetPanel(id = "riskstrat",
                                                                 tabPanel("Risk Straification Boxplot",
                                                                          p(),
                                                                          shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("Sboxplot", width = "100%", height = "400px")), type = 6),
                                                                          dnld_ui("dnldSboxplot_SVG","SVG"),
                                                                          div(DT::dataTableOutput("SboxplotTable"), style = "font-size:12px; height:450px; overflow-Y: scroll"),
                                                                          p(),
                                                                          dnld_ui("dnldSBoxplotTab","Download Table"),
                                                                          value = "box"
                                                                 ),
                                                                 tabPanel("Risk Straification Heatmap",
                                                                          p(),
                                                                          uiOutput("heatmap_error_message"),
                                                                          shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("Sheatmap", width = "100%", height = "2000px")), type = 6),
                                                                          fluidRow(
                                                                            dnld_ui("dnldSheatmap_SVG","SVG"),
                                                                            dnld_ui("dnldSheatmapexpr","Download Expression Matrix From Heatmap")
                                                                          ),
                                                                          value = "heat"
                                                                 )
                                                     )
                                            ),
                                            ##### Feat Strat ---------------------------------
                                            tabPanel("Feature Stratification",
                                                     p(),
                                                     fluidRow(
                                                       column(4,
                                                              selectizeInput("BoxplotFeature","Select Feature:",choices = NULL, selected = 1, width = "100%")
                                                       ),
                                                       column(2, style = "margin-top:12px",
                                                              checkboxGroupInput("FeatBoxOpts",NULL,choices = c("Remove NA/Unknowns","Log Transform Score"))
                                                       ),
                                                       column(3,
                                                              conditionalPanel(condition = "input.featstrat == 'box'",
                                                                               selectInput("boxoptselec","Boxplot Stat Compare Method:", choices = c("none","wilcox.test","t.test","kruskal.test","anova"))
                                                              ),
                                                              conditionalPanel(condition = "input.featstrat == 'heat'",
                                                                               selectizeInput("stratHeatAnno","Add Annotation:", choices = NULL, selected = 1, multiple = T)
                                                              )
                                                              
                                                       ),
                                                       column(3,
                                                              conditionalPanel(condition = "input.featstrat == 'heat'",
                                                                               selectInput("ClusterMethod2", "Select Cluster Method:",
                                                                                           choices = c("ward.D", "ward.D2", "complete", "single", "average", "mcquitty", "median", "centroid"))
                                                              )
                                                       )
                                                     ),
                                                     tabsetPanel(id = "featstrat",
                                                                 tabPanel("Feature Boxplot",
                                                                          p(),
                                                                          shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("Featureboxplot", width = "100%", height = "500px")), type = 6),
                                                                          dnld_ui("dnldFboxplot_SVG","SVG"),
                                                                          div(DT::dataTableOutput("FeatureboxplotTable"), style = "font-size:12px; height:450px; overflow-Y: scroll"),
                                                                          p(),
                                                                          dnld_ui("dnldFeatureboxplotTab","Download Table"),
                                                                          value = "box"
                                                                 ),
                                                                 tabPanel("Feature Heatmap",
                                                                          p(),
                                                                          uiOutput("heatmap_error_message2"),
                                                                          shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("FeatureHeatmap", width = "100%", height = "2000px")), type = 6),
                                                                          fluidRow(
                                                                            dnld_ui("dnldFheatmap_SVG","SVG"),
                                                                            dnld_ui("dnldFheatmapexpr","Download Expression Matrix From Heatmap")
                                                                          ),
                                                                          value = "heat"
                                                                 )
                                                     )
                                            )
                                          ),
                                          value = 4)
                               )
                             )
                           )
                         )
)

About_tab <- tabPanel("About",
                      fluidPage(
                        mainPanel(
                          tabPanel("Purpose and Methods",
                                   uiOutput("rendPurposeAndMethodsMD")
                          )
                        )
                      )
)

# Render UI --------------------------------------------------------------------

if (Password_Protected) {
  ui <- navbarPage(paste(ProjectName),
                   id = "tabs",
                   collapsible = TRUE,
                   login_tab)
} else {
  #if (FileProvided) {
  #  ui <- navbarPage(paste(ProjectName),
  #                   id = "tabs",
  #                   collapsible = TRUE,
  #                   Survival_tab,
  #                   About_tab)
  #} else {
  ui <- navbarPage(paste(ProjectName),
                   id = "tabs",
                   collapsible = TRUE,
                   DataInput_tab,
                   Survival_tab,
                   About_tab)
  #}
}



# Server --------------------------------------------------------------------

server <- function(input, output, session) {
  
  
  
  # hack to add the logout button to the navbar on app launch
  if (Password_Protected) {
    insertUI(
      selector = ".navbar .container-fluid .navbar-collapse",
      ui = tags$ul(
        class="nav navbar-nav navbar-right",
        tags$li(
          div(
            style = "padding: 10px; padding-top: 8px; padding-bottom: 0;",
            logoutUI("logout")
          )
        )
      )
    )
  }
  
  # call the shinyauthr login and logout server modules
  if (Password_Protected) {
    credentials <- loginServer(
      id = "login",
      data = user_base,
      user_col = "user",
      pwd_col = "password",
      sodium_hashed = FALSE,
      reload_on_logout = TRUE,
      log_out = reactive(logout_init())
    )
  } else {
    credentials <- reactive({
      list(user_auth = TRUE)
    })
  }
  
  if (Password_Protected) {
    logout_init <- logoutServer(
      id = "logout",
      active = reactive(credentials()$user_auth)
    )
  }
  
  
  observeEvent(credentials()$user_auth, {
    
    #browser()
    
    if (Password_Protected) {
      if (credentials()$user_auth) {
        # remove the login tab
        removeTab("tabs", "login")
        # add home tab
        if (!FileProvided) {
          appendTab("tabs", DataInput_tab, select = TRUE)
          appendTab("tabs", Survival_tab, select = FALSE)
          appendTab("tabs", About_tab, select = FALSE)
        } else {
          appendTab("tabs", DataInput_tab, select = TRUE)
          appendTab("tabs", Survival_tab, select = FALSE)
          appendTab("tabs", About_tab, select = FALSE)
        }
      }
    }
    
    if (credentials()$user_auth) {
      
      ## Set User Inputs
      ProjectName_react <- reactiveVal(ProjectName)
      ExpressionMatrix_file_react <- reactiveVal(ExpressionMatrix_file)
      MetaData_file_react <- reactiveVal(MetaData_file)
      MetaParam_File_react <- reactiveVal(MetaParam_File)
      
      PreSelect_SamplyType_react <- reactiveVal(PreSelect_SamplyType)
      PreSelect_Feature_react <- reactiveVal(PreSelect_Feature)
      PreSelect_SubFeature_react <- reactiveVal(PreSelect_SubFeature)
      PreSelect_SecondaryFeature_react <- reactiveVal(PreSelect_SecondaryFeature)
      
      ## Set file reactives
      expr_raw <- reactiveVal()
      expr_react <- reactiveVal()
      meta_react <- reactiveVal()
      metaP_react <- reactiveVal()
      metacol_samplename_react <- reactiveVal()
      
      # Data Input Tab ---------------------------------------------------------
      
      output$rendExprFileInput <- renderUI({
        
        refresh <- input$UseExpData
        fileInput("ExprFileInput","Expression Matrix")
        
      })
      output$rendClinFileInput <- renderUI({
        
        refresh <- input$UseExpData
        fileInput("ClinFileInput","Clinical Data")
        
      })
      
      output$rendClinParamFileInput <- renderUI({
        
        req(input$ParamChoice)
        if (input$ParamChoice == "Upload Parameter File") {
          refresh <- input$UseExpData
          fileInput("ClinParamFileInput","Clinical Parameter File")
        }
        
      })
      
      output$rendSurvTimeColSelect <- renderUI({
        req((isTruthy(ExpressionMatrix_file_react()) && isTruthy(MetaData_file_react())))
        req(input$ParamChoice)
        if (input$ParamChoice == "Define Parameters") {
          clin <- meta_react()
          ColumnCheck <- grep(paste(c("OS[[:punct:]]time","EFS[[:punct:]]time","PFI[[:punct:]]time","PFS[[:punct:]]time",
                                      "RFS[[:punct:]]time","DSS[[:punct:]]time","DFI[[:punct:]]time"),collapse = "|"),
                              colnames(clin), ignore.case = T, value = T)
          ColCheckFound <- colnames(clin)[which(colnames(clin) %in% ColumnCheck)]
          if (length(ColCheckFound) == 0) {
            ColCheckFound <- ""
          }
          selectInput("SurvTimeColSelect","Survival Time Column(s)",
                      choices = colnames(clin), selected = ColCheckFound, multiple = T)
        }
        
      })
      
      output$rendSurvTimeUnits <- renderUI({
        req((isTruthy(ExpressionMatrix_file_react()) && isTruthy(MetaData_file_react())))
        selectInput("SurvTimeUnits","Time Units:", choices = c("Days","Months","Years"))
        
      })
      
      output$rendSurvIDColSelect <- renderUI({
        req((isTruthy(ExpressionMatrix_file_react()) && isTruthy(MetaData_file_react())))
        req(input$ParamChoice)
        if (input$ParamChoice == "Define Parameters") {
          clin <- meta_react()
          ColumnCheckTime <- grep(paste(c("OS[[:punct:]]time","EFS[[:punct:]]time","PFI[[:punct:]]time","PFS[[:punct:]]time",
                                          "RFS[[:punct:]]time","DSS[[:punct:]]time","DFI[[:punct:]]time"),collapse = "|"),
                                  colnames(clin), ignore.case = T, value = T)
          ColumnCheck <- grep(paste(c("^OS","^EFS","^PFI","^PFS","^RFS","^DSS","^DFI",
                                      "OS[[:punct:]]ID","EFS[[:punct:]]ID","PFI[[:punct:]]ID","PFS[[:punct:]]ID",
                                      "RFS[[:punct:]]ID","DSS[[:punct:]]ID","DFI[[:punct:]]ID"),collapse = "|"),
                              colnames(clin), ignore.case = T, value = T)
          ColumnCheck <- ColumnCheck[which(!ColumnCheck %in% ColumnCheckTime)]
          ColCheckFound <- colnames(clin)[which(colnames(clin) %in% ColumnCheck)]
          if (length(ColCheckFound) == 0) {
            ColCheckFound <- ""
          }
          selectInput("SurvIDColSelect","Survival ID Column(s)",
                      choices = colnames(clin), selected = ColCheckFound, multiple = T)
        }
        
      })
      
      ## Update Inputs -------------------------------------------------------------
      
      output$rendExprFilePrevHeader <- renderUI({
        req(ExpressionMatrix_file_react())
        h3("Expression File Preview")
        
      })
      
      output$rendClinFilePrevHeader <- renderUI({
        req(MetaData_file_react())
        h3("Clinical File Preview")
        
      })
      
      output$rendParamFilePrevHeader <- renderUI({
        req(MetaParam_File_react())
        h3("Clinical Parameter File Preview")
        
      })
      
      observe({
        if (ParamFileProvided) {
          updateRadioButtons(session,"ParamChoice",selected = "Upload Parameter File")
        }
      })
      
      observe({
        query <- parseQueryString(session$clientData$url_search)
        print(query)
        if (isTruthy(query[['expr']]) && isTruthy(query[['meta']])) {
          ExpressionMatrix_file_react(query[['expr']])
          MetaData_file_react(query[['meta']])
        }
        if (isTruthy(query[['proj']])) {
          ProjectName_react(query[['proj']])
        }
        if (isTruthy(query[['param']])) {
          MetaParam_File_react(query[['param']])
          updateRadioButtons(session,"ParamChoice",selected = "Upload Parameter File")
        }
        if (isTruthy(query[['type']])) {
          PreSelect_SamplyType_react(query[['type']])
        }
        if (isTruthy(query[['feature']])) {
          PreSelect_Feature_react(query[['feature']])
        }
        if (isTruthy(query[['subfeature']])) {
          PreSelect_SubFeature_react(query[['subfeature']])
        }
        if (isTruthy(query[['secfeature']])) {
          PreSelect_SecondaryFeature_react(query[['secfeature']])
        }
      })
      
      # If load example data selected
      observeEvent(input$UseExpData, {
        
        ProjectName_react("TCGA_CHOL")
        ExpressionMatrix_file_react(ExampleExpr_File)
        MetaData_file_react(ExampleClin_File)
        MetaParam_File_react(ExampleParam_File)
        PreSelect_SamplyType_react("All")
        PreSelect_Feature_react("all")
        PreSelect_SubFeature_react(NULL)
        PreSelect_SecondaryFeature_react("ajcc_pathologic_tumor_stage")
        updateRadioButtons(session,"ParamChoice",selected = "Upload Parameter File")
        updateTextInput(session,"UserProjectName",value = "TCGA_CHOL")
      })
      
      # If user uploads data
      observe({
        if (isTruthy(input$ExprFileInput$datapath) & isTruthy(input$ClinFileInput$datapath)) {
          ProjectName_react(input$UserProjectName)
          ExpressionMatrix_file_react(input$ExprFileInput$datapath)
          MetaData_file_react(input$ClinFileInput$datapath)
          updateRadioButtons(session,"ParamChoice",selected = "Define Parameters")
        }
        if (isTruthy(input$ClinParamFileInput$datapath)) {
          MetaParam_File_react(input$ClinParamFileInput$datapath)
        }
        
      })
      
      ## Load Data --------------------------------------------------------------
      ## User Clinical Parameter Upload
      observe({
        
        if (isTruthy(input$ParamChoice)) {
          if (input$ParamChoice == "Define Parameters") {
            req(meta_react())
            clin <- meta_react()
            if (nrow(clin) > 0) {
              metacol_survid <- input$SurvIDColSelect
              metacol_survtime <- input$SurvTimeColSelect
              if (length(metacol_survtime) > 0 && length(metacol_survid) > 0) {
                metacol_feature <- colnames(clin)[-1]
                metacol_feature <- metacol_feature[!metacol_feature %in% c(metacol_survtime,metacol_survid)]
              }
              metacol_samplename <- colnames(clin)[1]
              if (length(metacol_survtime) > 0 && length(metacol_survid) > 0) {
                MetaParam <- rbind(data.frame(Clinical_Column_Name = metacol_samplename,
                                              Clinical_Column_Type = "SampleName"),
                                   data.frame(Clinical_Column_Name = metacol_survtime,
                                              Clinical_Column_Type = "SurvivalTime"),
                                   data.frame(Clinical_Column_Name = metacol_survid,
                                              Clinical_Column_Type = "SurvivalID"),
                                   data.frame(Clinical_Column_Name = metacol_feature,
                                              Clinical_Column_Type = "Feature"))
                metaP_react(MetaParam)
              }
            }
          } else {
            req(MetaParam_File_react())
            if (file.exists(MetaParam_File_react())) {
              if (tools::file_ext(MetaParam_File_react()) %in% c("txt","tsv","zip","gz","TXT","TSV","ZIP","GZ")) {
                metaP <- as.data.frame(fread(MetaParam_File_react(), sep = '\t', header = F))
              } else if (tools::file_ext(MetaParam_File_react()) %in% c("csv","CSV")) {
                metaP <- as.data.frame(fread(MetaParam_File_react(), sep = ',', header = F))
              } else if (tools::file_ext(MetaParam_File_react()) %in% c("RData","rdata")) {
                metaP <- loadRData(MetaParam_File_react())
              } else if (tools::file_ext(MetaParam_File_react()) %in% c("rds","RDS")) {
                metaP <- readRDS(MetaParam_File_react())
              }
              colnames(metaP) <- c("Clinical_Column_Name","Clinical_Column_Type")
              SampleNameRow <- which(metaP[,2] == "SampleName")
              metaP <- rbind(metaP[SampleNameRow,],metaP[-SampleNameRow,])
              metaP_react(metaP)
              # If expression is URL
            } else {
              if (tools::file_ext(MetaParam_File_react()) %in% c("txt","tsv","TXT","TSV")) {
                metaP <- as.data.frame(read_delim(url(MetaParam_File_react()), delim = '\t', col_names = F))
              } else if (tools::file_ext(MetaParam_File_react()) %in% c("zip","gz","ZIP","GZ")) {
                metaP <- as.data.frame(read_delim(getZip(MetaParam_File_react()), delim = '\t', col_names = F))
              } else if (tools::file_ext(MetaParam_File_react()) %in% c("csv","CSV")) {
                metaP <- as.data.frame(read_delim(url(MetaParam_File_react()), delim = ',', col_names = F))
              } else if (tools::file_ext(MetaParam_File_react()) %in% c("RData","rdata")) {
                metaP <- loadRData(url(MetaParam_File_react()))
              } else if (tools::file_ext(MetaParam_File_react()) %in% c("rds","RDS")) {
                metaP <- readRDS(url(MetaParam_File_react()))
              }
              colnames(metaP) <- c("Clinical_Column_Name","Clinical_Column_Type")
              SampleNameRow <- which(metaP[,2] == "SampleName")
              metaP <- rbind(metaP[SampleNameRow,],metaP[-SampleNameRow,])
              metaP_react(metaP)
            }
          }
        } else {
          req(MetaParam_File_react())
          if (file.exists(MetaParam_File_react())) {
            if (tools::file_ext(MetaParam_File_react()) %in% c("txt","tsv","zip","gz","TXT","TSV","ZIP","GZ")) {
              metaP <- as.data.frame(fread(MetaParam_File_react(), sep = '\t', header = F))
            } else if (tools::file_ext(MetaParam_File_react()) %in% c("csv","CSV")) {
              metaP <- as.data.frame(fread(MetaParam_File_react(), sep = ',', header = F))
            } else if (tools::file_ext(MetaParam_File_react()) %in% c("RData","rdata")) {
              metaP <- loadRData(MetaParam_File_react())
            } else if (tools::file_ext(MetaParam_File_react()) %in% c("rds","RDS")) {
              metaP <- readRDS(MetaParam_File_react())
            }
            colnames(metaP) <- c("Clinical_Column_Name","Clinical_Column_Type")
            SampleNameRow <- which(metaP[,2] == "SampleName")
            metaP <- rbind(metaP[SampleNameRow,],metaP[-SampleNameRow,])
            metaP_react(metaP)
            # If expression is URL
          }
        }
        
      })
      
      
      observe({
        
        if (isTruthy(ExpressionMatrix_file_react()) & isTruthy(MetaData_file_react())) {
          # Expression file
          # If expression is file
          if (file.exists(ExpressionMatrix_file_react())) {
            print(paste0("Loading in local file: ",ExpressionMatrix_file_react()))
            if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("txt","tsv","zip","gz","TXT","TSV","ZIP","GZ")) {
              expr <- as.data.frame(fread(ExpressionMatrix_file_react(), sep = '\t', header = T))
            } else if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("csv","CSV")) {
              expr <- as.data.frame(fread(ExpressionMatrix_file_react(), sep = ',', header = T))
            } else if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("RData","rdata")) {
              expr <- loadRData(ExpressionMatrix_file_react())
            } else if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("rds","RDS")) {
              expr <- readRDS(ExpressionMatrix_file_react())
            }
            # If expression is URL
          } else {
            print(paste0("Loading in url file: ",ExpressionMatrix_file_react()))
            if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("txt","tsv","TXT","TSV")) {
              expr <- as.data.frame(read_delim(url(ExpressionMatrix_file_react()), delim = '\t', col_names = T))
            } else if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("zip","gz","ZIP","GZ")) {
              expr <- as.data.frame(read_delim(getZip(ExpressionMatrix_file_react()), delim = '\t', col_names = T))
            } else if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("csv","CSV")) {
              expr <- as.data.frame(read_delim(url(ExpressionMatrix_file_react()), delim = ',', col_names = T))
            } else if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("RData","rdata")) {
              expr <- loadRData(url(ExpressionMatrix_file_react()))
            } else if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("rds","RDS")) {
              expr <- readRDS(url(ExpressionMatrix_file_react()))
            }
          }
          # Remove Expression with NA
          expr[,1] <- date_to_gene(expr[,1])
          expr <- expr %>%
            drop_na() %>%
            as.data.frame()
          # Check that expression data is numeric
          isChar <- unname(which(sapply(expr, function(x) is.character(x))))
          isChar <-  isChar[-1]
          if (length(isChar) > 0) {
            expr[isChar] <- sapply(expr[isChar],as.numeric)
          }
          colnames(expr)[1] <- "Gene"
          # Remove Duplicate genes
          expr_dup <- expr[which(expr[,1] %in% expr[,1][duplicated(expr[,1])]),]
          expr_nondup <- expr[which(!expr[,1] %in% expr[,1][duplicated(expr[,1])]),]
          if (nrow(expr_dup) > 0) {
            expr_dup <- expr_dup %>%
              group_by(Gene) %>%
              summarise_all(max) %>%
              as.data.frame()
          }
          expr <- rbind(expr_dup,expr_nondup)
          row.names(expr) <- expr[,1]
          expr <- expr[,-1]
          
          expr <- as.matrix(expr[order(rownames(expr)),])
          
          # Meta file
          # If Meta is file
          if (file.exists(MetaData_file_react())) {
            print(paste0("Loading in local file: ",MetaData_file_react()))
            if (tools::file_ext(MetaData_file_react()) %in% c("txt","tsv","zip","gz","TXT","TSV","ZIP","GZ")) {
              meta <- as.data.frame(fread(MetaData_file_react(), sep = '\t', header = T))
            } else if (tools::file_ext(MetaData_file_react()) %in% c("csv","CSV")) {
              meta <- as.data.frame(fread(MetaData_file_react(), sep = ',', header = T))
            } else if (tools::file_ext(MetaData_file_react()) %in% c("RData","rdata")) {
              meta <- loadRData(MetaData_file_react())
            } else if (tools::file_ext(MetaData_file_react()) %in% c("rds","RDS")) {
              meta <- readRDS(MetaData_file_react())
            }
            # If Meta is URL
          } else {
            print(paste0("Loading in url file: ",MetaData_file_react()))
            if (tools::file_ext(MetaData_file_react()) %in% c("txt","tsv","TXT","TSV")) {
              meta <- as.data.frame(read_delim(url(MetaData_file_react()), delim = '\t', col_names = T))
            } else if (tools::file_ext(MetaData_file_react()) %in% c("zip","gz","ZIP","GZ")) {
              meta <- as.data.frame(read_delim(getZip(MetaData_file_react()), delim = '\t', col_names = T))
            } else if (tools::file_ext(MetaData_file_react()) %in% c("csv","CSV")) {
              meta <- as.data.frame(read_delim(url(MetaData_file_react()), delim = ',', col_names = T))
            } else if (tools::file_ext(MetaData_file_react()) %in% c("RData","rdata")) {
              meta <- loadRData(url(MetaData_file_react()))
            } else if (tools::file_ext(MetaData_file_react()) %in% c("rds","RDS")) {
              meta <- readRDS(url(MetaData_file_react()))
            }
          }
          
          expr_raw(expr)
          expr_react(expr)
          meta_react(meta)
          
        }
        
      })
      
      observe({
        
        if (isTruthy(expr_react()) & isTruthy(meta_react()) & isTruthy(metaP_react())) {
          meta <- meta_react()
          metaP <- metaP_react()
          expr <- expr_react()
          
          metacol_samplename <- metaP[which(metaP[,2] == "SampleName"),1]
          metacol_samplename_react(metacol_samplename)
          if (metacol_samplename %in% colnames(meta)) {
            sampsames <- intersect(colnames(expr),meta[,metacol_samplename])
            #ensure expression samples and meta are exact
            expr <- expr[,sampsames]
            meta <- meta[which(meta[,metacol_samplename] %in% sampsames),]
            meta <- meta %>% relocate(any_of(metacol_samplename))
            
            PreProcessed_meta_cols <- c(grep("_PreProcessedScore$",colnames(meta),value = T))
            if (length(PreProcessed_meta_cols) == 0) {
              if (immudecon_check == TRUE) {
                mcp_counter_decon <- as.data.frame(deconvolute(expr, "mcp_counter"))
                rownames(mcp_counter_decon) <- mcp_counter_decon[,1]
                mcp_counter_decon <- mcp_counter_decon[,-1]
                mcp_counter_decon <- as.data.frame(t(mcp_counter_decon))
                colnames(mcp_counter_decon) <- paste(gsub(" ","_",colnames(mcp_counter_decon)),"mcp_counter_InApp_PreProcessedScore",sep = "_")
                metaP <- rbind(metaP,
                               data.frame(Clinical_Column_Name = colnames(mcp_counter_decon),
                                          Clinical_Column_Type = "Feature"))
                mcp_counter_decon[,metacol_samplename] <- rownames(mcp_counter_decon)
                meta <- merge(meta,mcp_counter_decon)
                estimate_decon <- as.data.frame(deconvolute(expr, "estimate"))
                rownames(estimate_decon) <- estimate_decon[,1]
                estimate_decon <- estimate_decon[,-1]
                estimate_decon <- as.data.frame(t(estimate_decon))
                colnames(estimate_decon) <- paste(gsub(" ","_",colnames(estimate_decon)),"estimate_InApp_PreProcessedScore",sep = "_")
                metaP <- rbind(metaP,
                               data.frame(Clinical_Column_Name = colnames(estimate_decon),
                                          Clinical_Column_Type = "Feature"))
                estimate_decon[,metacol_samplename] <- rownames(estimate_decon)
                meta <- merge(meta,estimate_decon)
                metaP_react(metaP)
                meta_react(meta)
                expr_react(expr)
              } else {
                expr_react(expr)
                meta_react(meta)
              }
            } else {
              expr_react(expr)
              meta_react(meta)
            }
            expr_react(expr)
            meta_react(meta)
          }
        }
        
      })
      
      decon_score_cols <- reactive({
        req(meta_react())
        meta <- meta_react()
        decon_score_cols <- grep("_PreProcessedScore$", colnames(meta), value = T, ignore.case = T)
        decon_score_cols
      })
      
      observe({
        req(decon_score_cols())
        if (isTruthy(decon_score_cols())) {
          updateSelectInput(session,"GeneSetCat_Select",choices = c(GeneSetCats,"Pre-Processed Scores"))
        }
      })
      
      observeEvent(input$SurvTimeUnits,{
        if (isTruthy(input$SurvTimeUnits)) {
          meta <- meta_react()
          metaP <- metaP_react()
          TimCols <- metaP[which(metaP[,2] == "SurvivalTime"),1]
          if (input$SurvTimeUnits == "Months") {
            for (i in TimCols) {
              meta[,i] <- meta[,i] * 30.4375
            }
          } else if (input$SurvTimeUnits == "Years") {
            for (i in TimCols) {
              meta[,i] <- meta[,i] * 365.25
            }
          }
          meta_react(meta)
        }
        
      })
      
      observeEvent(input$LogExprFile | input$ScaleNormExprFile,{
        req(input$LogExprFile)
        req(input$ScaleNormExprFile)
        req(meta_react())
        req(metaP_react())
        req(expr_raw())
        meta <- meta_react()
        metaP <- metaP_react()
        expr <- expr_raw()
        
        if (input$LogExprFile == T) {
          expr <- log2(expr + 0.00001)
          meta <- meta[,grep("_InApp_PreProcessedScore$",colnames(meta), invert = T)]
          metaP <- metaP[which(metaP[,1] %in% grep("_InApp_PreProcessedScore$",metaP[,1], invert = T, value = T)),]
          meta_react(meta)
          metaP_react(metaP)
          expr <- as.matrix(expr[sort(rownames(expr)),])
          expr_react(expr)
        } else {
          meta_react(meta)
          metaP_react(metaP)
          expr_react(expr)
        }
        if (input$ScaleNormExprFile == T) {
          expr_col <- colnames(expr)
          expr = apply(expr, 1, scale)
          expr = apply(expr, 1, rev)
          colnames(expr) <- expr_col
          meta <- meta[,grep("_InApp_PreProcessedScore$",colnames(meta), invert = T)]
          metaP <- metaP[which(metaP[,1] %in% grep("_InApp_PreProcessedScore$",metaP[,1], invert = T, value = T)),]
          meta_react(meta)
          metaP_react(metaP)
          expr <- as.matrix(expr[sort(rownames(expr)),])
          expr_react(expr)
        } else {
          meta_react(meta)
          metaP_react(metaP)
          expr_react(expr)
        }
      })
      
      metacol_sampletype <- reactive({
        req(metaP_react())
        metaP <- metaP_react()
        if ("SampleType" %in% metaP[,2]) {
          metacol_sampletype <- metaP[which(metaP[,2] == "SampleType"),1]
        }
        if (!"SampleType" %in% metaP[,2]) {
          metacol_sampletype <- NULL
        }
        metacol_sampletype
        
      })
      
      metacol_survtime <- reactive({
        
        req(metaP_react())
        metaP <- metaP_react()
        metacol_survtime <- metaP[which(metaP[,2] == "SurvivalTime"),1]
        metacol_survtime
        
      })
      
      metacol_survid <- reactive({
        
        req(metaP_react())
        metaP <- metaP_react()
        metacol_survid <- metaP[which(metaP[,2] == "SurvivalID"),1]
        metacol_survid
        
      })
      
      metacol_feature <- reactive({
        req(metaP_react)
        metaP <- metaP_react()
        geneset <- gs_react()
        geneset_name <- names(geneset)
        metacol_feature <- metaP[which(metaP[,2] == "Feature"),1]
        if ("SampleType" %in% metaP[,2]) {
          metacol_feature <- ifelse(is.null(metacol_sampletype()),metacol_feature,c(metacol_sampletype(),metacol_feature))
          metacol_feature <- metacol_feature[-which(metacol_feature == input$FeatureSelection)]
          metacol_feature <- ifelse(!is.null(input$SampleTypeSelection),metacol_feature[-which(metacol_feature == input$SampleTypeSelection)],metacol_feature)
        }
        metacol_feature <- c(metacol_feature,geneset_name,paste0(geneset_name,"_QuartileCutP"),paste0(geneset_name,"_MedianCutP"),
                             paste0(geneset_name,"_OptimalCutP"),paste0(geneset_name,"_TopBottomCutP"),paste0(geneset_name,"_UserCutP"))
        metacol_feature
        
      })
      
      ## Data Preview -----------------------------------------------------------
      output$ExprFile_Preview <- DT::renderDataTable({
        
        expr <- expr_react()
        DT::datatable(expr,
                      options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                     pageLength = 5,
                                     scrollX = T),
                      rownames = T)
        
      })
      
      output$ClinFile_Preview <- DT::renderDataTable({
        
        clin <- meta_react()
        DT::datatable(clin,
                      options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                     pageLength = 5,
                                     scrollX = T),
                      rownames = F)
        
      })
      
      output$ClinParamFile_Preview <- DT::renderDataTable({
        
        clinP <- metaP_react()
        DT::datatable(clinP,
                      options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                     pageLength = 5,
                                     scrollX = T),
                      rownames = F)
        
      })
      
      # Surv Tab -----------------------------------------------------------
      
      ## Sidebar ---------------------------------------------------------------
      ### Choices --------------------------------------------------------------
      ## Select sample type to subset samples by - only render if more than one sample type
      output$rendSampleTypeSelection <- renderUI({
        metacol_sampletype <- metacol_sampletype()
        meta <- meta_react()
        if (length(unique(meta[,metacol_sampletype])) > 1) {
          SampleTypeChoices <- unique(meta[,metacol_sampletype])
          SampleTypeChoices <- c("All Sample Types",SampleTypeChoices)
          if (toupper(PreSelect_SamplyType_react()) == "ALL") {
            preselect <- "All Sample Types"
          } else {
            preselect <- PreSelect_SamplyType_react()
          }
          selectInput("SampleTypeSelection",paste("Select Sample Type (",metacol_sampletype,"):",sep = ""),
                      choices = SampleTypeChoices, selected = preselect)
        }
      })
      
      observe({
        req(metaP_react())
        MetaParam <- metaP_react()
        metacol_feature <- c("Show all Samples",MetaParam[,1])
        featSelec <- 1
        if (!PreSelect_Feature_react() %in% MetaParam[,1] & toupper(PreSelect_Feature_react()) != "ALL") {
          featSelec <- "Show all Samples"
        }
        if (toupper(PreSelect_Feature_react()) == "ALL") {
          featSelec <- "Show all Samples"
        }
        if (is.null(PreSelect_Feature_react())) {
          featSelec <- "Show all Samples"
        }
        if (PreSelect_Feature_react() %in% MetaParam[,1]) {
          featSelec <- PreSelect_Feature_react()
        }
        updateSelectizeInput(session = session, inputId = "FeatureSelection",
                             choices = metacol_feature, selected = featSelec, server = T)
      })
      
      output$rendSubFeatureSelection <- renderUI({
        req(input$FeatureSelection)
        meta <- meta_react()
        if (isTruthy(input$SampleTypeSelection)) {
          if (input$SampleTypeSelection != "All Sample Types") {
            metacol_sampletype <- metacol_sampletype()
            meta <- meta[which(meta[,metacol_sampletype] == input$SampleTypeChoices),]
          } else {
            meta <- meta
          }
        }
        if (input$FeatureSelection != "Show all Samples") {
          SubFeatureChoices <- unique(meta[,input$FeatureSelection])
          SubFeatureChoices <- sort(SubFeatureChoices, decreasing = T, na.last = T)
          selectInput("subFeatureSelection","Feature Condition:",choices = SubFeatureChoices, selected = PreSelect_SubFeature_react())
        }
      })
      
      observe({
        req(metacol_survtime())
        SurTimeChoices <- metacol_survtime()
        updateSelectizeInput(session = session, inputId = "SurvivalType_time",choices = SurTimeChoices, server = T)
      })
      observe({
        req(metacol_survid())
        SurIDChoices <- metacol_survid()
        updateSelectizeInput(session = session, inputId = "SurvivalType_id",choices = SurIDChoices, server = T)
      })
      
      #observe({
      #  req(meta_react())
      #  req(input$SurvivalType_time)
      #  req(input$SurvivalType_id)
      #  req(input$SurvYearOrMonth)
      #  meta_ssgsea <- meta_react()
      #  yOm <- input$SurvYearOrMonth
      #  surv_time_col <- input$SurvivalType_time
      #  max_time <- ceiling(max(meta_ssgsea[,surv_time_col], na.rm = T)/365.25)
      #  if (yOm == "Years") {
      #    updateNumericInput(session,"SurvXaxis","Survival Age Limit (years)", value = max_time)
      #  } else if (yOm == "Months") {
      #    updateNumericInput(session,"SurvXaxis","Survival Age Limit (months)", value = max_time*12)
      #  }
      #})
      
      #output$rendSurvXaxis <- renderUI({
      #  req(ssGSEAmeta())
      #  req(input$SurvivalType_time)
      #  req(input$SurvivalType_id)
      #  req(input$SurvYearOrMonth)
      #  meta_ssgsea <- ssGSEAmeta()
      #  yOm <- input$SurvYearOrMonth
      #  surv_time_col <- input$SurvivalType_time
      #  surv_id_col <- input$SurvivalType_id
      #  max_time <- ceiling(max(meta_ssgsea[,surv_time_col])/365.25)
      #  if (yOm == "Years") {
      #    numericInput("SurvXaxis","Survival Age Limit (years)", value = max_time)
      #  } else if (yOm == "Months") {
      #    numericInput("SurvXaxis","Survival Age Limit (months)", value = max_time*12)
      #  }
      #})
      
      observe({
        req(input$SurvYearOrMonth)
        if (input$SurvYearOrMonth == "Months") {
          updateNumericInput(session,"SurvXaxisBreaks",label = "Survival X-Axis Breaks (Months):", value = 12, step = 1)
        } else if (input$SurvYearOrMonth == "Years") {
          updateNumericInput(session,"SurvXaxisBreaks",label = "Survival X-Axis Breaks (Years):", value = 1, step = 1)
        }
      })
      
      ### Genesets -------------------------------------------------------------
      GeneSetTable_React <- reactive({
        GS_database <- input$GeneSetCat_Select
        if (GS_database %in% GeneSetCats) {
          sub_tab <- GeneSetTable[which(GeneSetTable[,1] == GS_database),]
          new_tab <- sub_tab[,-1]
        } else if (GS_database == "Pre-Processed Scores") {
          new_tab <- data.frame(scores = decon_score_cols())
          colnames(new_tab)[1] <- "Pre-Processed Scores"
        }
        new_tab
      })
      
      ## Render Gene Set Selection Table
      output$GeneSetTable <- DT::renderDataTable({
        DT::datatable(GeneSetTable_React(),
                      options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                     pageLength = 10,
                                     scrollX = T),
                      selection = list(mode = 'single', selected = 1),
                      rownames = F)
      })
      
      GeneGS_table <- reactive({
        req(expr_react())
        expr <- expr_react()
        GenesDF <- data.frame(Genes = rownames(expr))
        GenesDF
      })
      
      output$GeneGeneSetTable <- DT::renderDataTable({
        GenesDF <- GeneGS_table()
        DT::datatable(GenesDF,
                      options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                     pageLength = 10,
                                     scrollX = T),
                      selection = list(mode = 'single', selected = 1),
                      rownames = F)
      })
      
      ## Render User Gene Set Table - Backend
      userGeneSet_loaded <- reactive({
        gs.u <- input$userGeneSet
        ext <- tools::file_ext(gs.u$datapath)
        req(gs.u)
        validate(need(ext == c("gmt","tsv","txt", "RData"), "Please upload .gmt, .tsv, .txt, or .RData file"))
        # If user provides GMT file
        if (ext == "gmt") {
          gmt <- clusterProfiler::read.gmt(gs.u$datapath)
        } else if (ext == "RData") {
          gmt <- loadRData(gs.u$datapath)
        } else {
          #gmt <- as.data.frame(read_delim(gs.u$datapath, delim = '\t'))
          gmt <- as.data.frame(fread(gs.u$datapath, sep = '\t'))
        }
        gmt
      })
      
      userGeneSet_table <- reactive({
        gs.u <- input$userGeneSet
        ext <- tools::file_ext(gs.u$datapath)
        gmt <- userGeneSet_loaded()
        # If user provides GMT file
        if (ext == "gmt") {
          uGS_table <- as.data.frame(unique(gmt[,1]))
          colnames(uGS_table)[1] <- "GeneSet"
        } else if (ext == "RData") {
          uGS_table <- as.data.frame(names(gmt))
          colnames(uGS_table)[1] <- "GeneSet"
        } else {
          uGS_table <- as.data.frame(unique(gmt[,1]))
          colnames(uGS_table)[1] <- "GeneSet"
        }
        uGS_table
      })
      
      output$userGeneSetTable <- DT::renderDataTable({
        req(input$userGeneSet)
        uGS_table <- userGeneSet_table()
        DT::datatable(uGS_table,
                      options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                     pageLength = 10,
                                     scrollX = T),
                      selection = list(mode = 'single', selected = 1),
                      rownames = F)
      })
      
      gs_react <- reactive({
        req(input$GeneSetTabs)
        if (input$GeneSetTabs == 1) {
          req(GeneSetTable_React())
          if (input$GeneSetCat_Select == "Pre-Processed Scores") {
            geneset_name <- GeneSetTable_React()[input$GeneSetTable_rows_selected,1]
            meta <- meta_react()
            geneset <- list(meta[,geneset_name])
            names(geneset) <- geneset_name
            geneset
          } else {
            geneset_name <- GeneSetTable_React()[input$GeneSetTable_rows_selected,3]
            geneset <- gs[geneset_name]
            geneset
          }
        } else if (input$GeneSetTabs == 2) {
          req(GeneGS_table())
          gene <- GeneGS_table()[input$GeneGeneSetTable_rows_selected,1]
          if (isTruthy(gene)) {
            geneset <- list(gene = gene)
            names(geneset)[1] <- gene
            geneset
          }
        } else if (input$GeneSetTabs == 3) {
          if (input$UserGSoption == "Gene Set Upload") {
            req(input$userGeneSet)
            gs.u <- input$userGeneSet
            ext <- tools::file_ext(gs.u$datapath)
            gmt <- userGeneSet_loaded()
            uGS_table <- userGeneSet_table()
            geneset_name <- uGS_table[input$userGeneSetTable_rows_selected,1]
            if (isTruthy(geneset_name)) {
              if (ext == "gmt") {
                geneset <- list(geneset_name = gmt[which(gmt[,1] == geneset_name),2])
                names(geneset) <- geneset_name
              } else if (ext == "RData") {
                geneset <- gmt[geneset_name]
              } else {
                geneset <- list(geneset_name = gmt[which(gmt[,1] == geneset_name),2])
                names(geneset) <- geneset_name
              }
              geneset
            }
          }
          else if (input$UserGSoption == "Text Box Input") {
            user_gs_text <- input$userGeneSetText
            user_gs_name <- input$userGeneSetTextName
            user_gs_name <- gsub(" ",".",user_gs_name)
            gs_text_s <- unlist(strsplit(user_gs_text, " "))
            gs_text_t <- unlist(strsplit(user_gs_text, "\t"))
            gs_text_c <- unlist(strsplit(user_gs_text, ","))
            gs_text <- unique(c(gs_text_s,gs_text_t,gs_text_c))
            geneset <- list(gs_text)
            names(geneset) <- user_gs_name
            geneset
          }
        }
      })
      
      output$rendGenesInGeneSetTab <- renderUI({
        if (input$ViewGeneSetGenes) {
          div(DT::dataTableOutput("GenesInGeneSetTab"), style = "font-size:10px")
        }
      })
      
      output$GenesInGeneSetTab <- DT::renderDataTable({
        geneset <- gs_react()
        geneset_df <- as.data.frame(geneset)
        DT::datatable(geneset_df,
                      options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                     pageLength = 10,
                                     scrollX = T),
                      rownames = F)
      })
      
      metaSub <- reactive({
        req(meta_react())
        req(input$FeatureSelection)
        meta <- meta_react()
        metacol_sampletype <- metacol_sampletype()
        SampleTypeSelected <- input$SampleTypeChoices
        Feature <- input$FeatureSelection
        subFeature <- input$subFeatureSelection
        if (isTruthy(SampleTypeSelected)) {
          if (SampleTypeSelected != "All Sample Types") {
            meta <- meta[which(meta[,metacol_sampletype] == SampleTypeSelected),]
          } else {
            meta <- meta
          }
        } else {
          meta <- meta
        }
        if (Feature != "Show all Samples") {
          meta <- meta[which(meta[,Feature] == subFeature),]
        } else {
          meta <- meta
        }
        meta
      })
      
      exprSub <- reactive({
        req(metaSub())
        req(expr_react())
        expr <- expr_react()
        meta <- metaSub()
        expr <- expr[,meta[,1], drop = F]
        expr
      })
      
      ### ssGSEA Reactive ------------------------------------------------------
      ## Perform ssGSEA and functions on new meta table
      ssGSEAmeta <- reactive({
        req(metaSub())
        req(exprSub())
        req(gs_react())
        #req(input$SurvXaxis)
        meta <- metaSub()
        expr <- exprSub()
        geneset <- gs_react()
        geneset_name <- names(geneset)
        if (length(unname(unlist(geneset))) > 0) {
          if (any(unname(unlist(geneset)) %in% rownames(expr))) {
            if (isTruthy(input$SurvivalType_time) & isTruthy(input$SurvivalType_id) & isTruthy(geneset_name)) {
              SampleNameCol <- metacol_samplename_react()
              quantCutoff <- input$QuantPercent/100
              quantCutoff2 <- input$QuantPercent2/100
              surv_time_col <- input$SurvivalType_time
              surv_id_col <- input$SurvivalType_id
              scoreMethod <- input$ScoreMethod
              ## Remove rows with NA in survival column
              meta <- meta[!is.na(meta[,surv_time_col]),]
              meta <- meta[!is.na(meta[,surv_id_col]),]
              meta[,surv_id_col] <- as.numeric(meta[,surv_id_col])
              meta[,surv_time_col] <- as.numeric(meta[,surv_time_col])
              
              ## Re-subset expression matrix
              if (length(!is.na(meta[,surv_id_col])) > 0) {
                if (any(colnames(expr) %in% meta[,1])) {
                  expr_sub <- expr[,colnames(expr) %in% meta[,1], drop = F]
                  expr_mat <- as.matrix(expr_sub)
                  rownames(expr_mat) <- rownames(expr_sub)
                  colnames(expr_mat) <- colnames(expr_sub)
                  if (input$GeneSetTabs == 1 | input$GeneSetTabs == 3) {
                    if (geneset_name %in% decon_score_cols()) {
                      ssGSEA <- meta[,c(SampleNameCol,geneset_name)]
                      ssGSEA <- ssGSEA[!is.na(ssGSEA[,2]),]
                      ssGSEA <- ssGSEA[which(ssGSEA[,2] != "Inf"  & ssGSEA[,2] != "N/A" & ssGSEA[,2] != "n/a"),]
                      ssGSEA[,2] <- as.numeric(ssGSEA[,2])
                    }
                    else {
                      if (as.numeric(tools::file_path_sans_ext(packageVersion("GSVA"))) >= 1.5) {
                        if (scoreMethod == "ssgsea") {
                          ssGSEA_param <- GSVA::ssgseaParam(expr_mat,geneset)
                        } else if (scoreMethod == "gsva") {
                          ssGSEA_param <- GSVA::gsvaParam(expr_mat,geneset)
                        } else if (scoreMethod == "plage") {
                          ssGSEA_param <- GSVA::plageParam(expr_mat,geneset)
                        } else if (scoreMethod == "zscore") {
                          ssGSEA_param <- GSVA::zscoreParam(expr_mat,geneset)
                        }
                        ssGSEA <- GSVA::gsva(ssGSEA_param)
                        ssGSEA <- as.data.frame(t(ssGSEA))
                        ssGSEA[,SampleNameCol] <- rownames(ssGSEA)
                      } else {
                        ssGSEA <- gsva(expr_mat,geneset,method = scoreMethod, verbose = FALSE)
                        ssGSEA <- as.data.frame(t(ssGSEA))
                        ssGSEA[,SampleNameCol] <- rownames(ssGSEA)
                      }
                    }
                  } else if (input$GeneSetTabs == 2) {
                    expr_sub <- as.data.frame(expr_sub[geneset_name,])
                    colnames(expr_sub)[1] <- geneset_name
                    expr_sub[,SampleNameCol] <- rownames(expr_sub)
                    ssGSEA <- expr_sub
                  }
                  
                  #ssGSEA
                  
                  # Subset columns needed for plot and rename for surv function
                  #print(head(meta[,c(SampleNameCol,surv_time_col,surv_id_col)]))
                  #print(head(ssGSEA[,c(SampleNameCol,geneset_name)]))
                  meta_ssgsea_sdf <- merge(meta[,c(SampleNameCol,surv_time_col,surv_id_col)],ssGSEA[,c(SampleNameCol,geneset_name)])
                  
                  #print(input$SurvXaxis)
                  #print(input$SurvYearOrMonth)
                  if (isTruthy(input$SurvXaxis)) {
                    if (input$SurvYearOrMonth == "Months") {
                      xaxlim <- input$SurvXaxis * 30.4375
                      meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,surv_time_col] <= xaxlim),]
                      ssGSEA <- ssGSEA[meta_ssgsea_sdf[,SampleNameCol],]
                    } else {
                      xaxlim <- input$SurvXaxis * 365.25
                      meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,surv_time_col] <= xaxlim),]
                      ssGSEA <- ssGSEA[meta_ssgsea_sdf[,SampleNameCol],]
                    }
                  } else {
                    xaxlim <- max(meta_ssgsea_sdf[,surv_time_col],na.rm = T)
                    meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,surv_time_col] <= xaxlim),]
                    ssGSEA <- ssGSEA[meta_ssgsea_sdf[,SampleNameCol],]
                  }
                  
                  #print(head(meta_ssgsea_sdf))
                  #save(list = ls(), file = "ssgsea_env.RData", envir = environment())
                  
                  if (length(meta_ssgsea_sdf[,4][meta_ssgsea_sdf[,4] > 0])/length(meta_ssgsea_sdf[,4]) > 0.01) {
                    if (length(meta_ssgsea_sdf[,4]) > 1) {
                      res.cut <- survminer::surv_cutpoint(meta_ssgsea_sdf,time = surv_time_col, event = surv_id_col, variable = geneset_name, minprop = 0.01)
                      cutp <- res.cut$cutpoint[["cutpoint"]]
                      res.cat <- surv_categorize(res.cut)
                      ssGSEA[,paste0(geneset_name,"_OptimalCutP")] <- res.cat[,3]
                    }
                  }
                  
                  ssGSEA[,paste0(geneset_name,"_VAR_Q")] <- quartile_conversion(ssGSEA[, which(colnames(ssGSEA) == geneset_name)])
                  ssGSEA[,paste0(geneset_name,"_QuartileCutP")] <- paste("", ssGSEA[,paste0(geneset_name,"_VAR_Q")], sep="")
                  ssGSEA[,paste0(geneset_name,"_MedianCutP")] <- highlow(ssGSEA[, which(colnames(ssGSEA) == geneset_name)])
                  ssGSEA[,paste0(geneset_name,"_TopBottomCutP")] <- quantile_conversion(ssGSEA[, which(colnames(ssGSEA) == geneset_name)], quantCutoff)
                  ssGSEA[which(ssGSEA[,paste0(geneset_name,"_TopBottomCutP")] == "BetweenCutoff"),paste0(geneset_name,"_TopBottomCutP")] <- NA
                  ssGSEA[,paste0(geneset_name,"_UserCutP")] <- quantile_conversion2(ssGSEA[, which(colnames(ssGSEA) == geneset_name)], quantCutoff2)
                  
                  meta_ssGSEA <- merge(meta,ssGSEA)
                  meta_ssGSEA
                }
              }
              
              
            }
          }
        }
        
      })
      
      #ssGSEAmeta <- reactive({
      #  
      #  req(ssGSEAmeta_pre())
      #  req(metaSub())
      #  ssGSEA <- ssGSEAmeta_pre()
      #  meta <- metaSub()
      #  geneset_name <- colnames(ssGSEA)[1]
      #  quantCutoff <- input$QuantPercent/100
      #  quantCutoff2 <- input$QuantPercent2/100
      #  surv_time_col <- input$SurvivalType_time
      #  surv_id_col <- input$SurvivalType_id
      #  SampleNameCol <- metacol_samplename_react()
      #  
      #  ## Remove rows with NA in survival column
      #  meta <- meta[!is.na(meta[,surv_time_col]),]
      #  meta <- meta[!is.na(meta[,surv_id_col]),]
      #  meta[,surv_id_col] <- as.numeric(meta[,surv_id_col])
      #  meta[,surv_time_col] <- as.numeric(meta[,surv_time_col])
      #  meta_ssgsea_sdf <- merge(meta[,c(SampleNameCol,surv_time_col,surv_id_col)],ssGSEA[,c(SampleNameCol,geneset_name)])
      #  
      #  save(list = ls(), file = "ssGSEA_React.RData",envir = environment())
      #
      #  if (isTruthy(input$SurvXaxis)) {
      #    if (input$SurvYearOrMonth == "Months") {
      #      xaxlim <- input$SurvXaxis * 30.4375
      #      meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,2] <= xaxlim),]
      #      ssGSEA <- ssGSEA[meta_ssgsea_sdf[,1],]
      #    } else {
      #      xaxlim <- input$SurvXaxis * 365.25
      #      meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,2] <= xaxlim),]
      #      ssGSEA <- ssGSEA[meta_ssgsea_sdf[,1],]
      #    }
      #  }
      #  
      #  if (isTruthy(meta_ssgsea_sdf)) {
      #    if (length(meta_ssgsea_sdf[,4][meta_ssgsea_sdf[,4] > 0])/length(meta_ssgsea_sdf[,4]) > 0.01) {
      #      if (length(meta_ssgsea_sdf[,4]) > 1) {
      #        res.cut <- survminer::surv_cutpoint(meta_ssgsea_sdf,time = surv_time_col, event = surv_id_col, variable = geneset_name, minprop = 0.01)
      #        cutp <- res.cut$cutpoint[["cutpoint"]]
      #        res.cat <- surv_categorize(res.cut)
      #        ssGSEA[,paste0(geneset_name,"_OptimalCutP")] <- res.cat[,3]
      #      }
      #    }
      #  }
      #  
      #  
      #  ssGSEA[,paste0(geneset_name,"_VAR_Q")] <- quartile_conversion(ssGSEA[,1])
      #  ssGSEA[,paste0(geneset_name,"_QuartileCutP")] <- paste("", ssGSEA[,paste0(geneset_name,"_VAR_Q")], sep="")
      #  ssGSEA[,paste0(geneset_name,"_MedianCutP")] <- highlow(ssGSEA[,1])
      #  ssGSEA[,paste0(geneset_name,"_TopBottomCutP")] <- quantile_conversion(ssGSEA[,1], quantCutoff)
      #  ssGSEA[which(ssGSEA[,paste0(geneset_name,"_TopBottomCutP")] == "BetweenCutoff"),paste0(geneset_name,"_TopBottomCutP")] <- NA
      #  ssGSEA[,paste0(geneset_name,"_UserCutP")] <- quantile_conversion2(ssGSEA[,1], quantCutoff2)
      #  
      #  meta_ssGSEA <- merge(meta,ssGSEA)
      #  meta_ssGSEA
      #  
      #})
      
      
      ## Median CutP -----------------------------------------------------------
      
      MedianCutP_react <- reactive({
        req(ssGSEAmeta())
        geneset <- gs_react()
        geneset_name <- names(geneset)
        SubsetSurvData(ssGSEAmeta(),input$SurvivalType_time,input$SurvivalType_id,paste0(geneset_name,"_MedianCutP"))
      })
      
      MedianCutPTab_react <- reactive({
        req(MedianCutP_react())
        geneset <- gs_react()
        geneset_name <- names(geneset)
        MedianCutP_react <- MedianCutP_react()
        #save(list = ls(), file = "Surv_env.RData", envir = environment())
        CoxPHobj(MedianCutP_react(),paste0(geneset_name,"_MedianCutP"),"Low")
      })
      
      SBinaryHRtab_react <- reactive({
        req(MedianCutPTab_react())
        CoxPHtabUni(MedianCutPTab_react())
      })
      
      output$BINSurvDescrip <- renderUI({
        HR_Tab <- SBinaryHRtab_react()
        Pval_Tab <- MedianCutPTab_react()
        metacol_sampletype <- metacol_sampletype()
        SampleTypeSelected <- input$SampleTypeChoices
        Feature <- input$FeatureSelection
        subFeature <- input$subFeatureSelection
        geneset <- gs_react()
        geneset_name <- names(geneset)
        decon_score_cols <- decon_score_cols()
        surv_time_col <- input$SurvivalType_time
        scoreMethod <- input$ScoreMethod
        
        if (geneset_name %in% decon_score_cols) {
          scoreMethod <- "Pre-Processed Score"
        }
        if (input$GeneSetTabs == 2) {
          scoreMethod <- "Gene Expression"
        }
        SurvPlotExpl("median cut-point",surv_time_col,geneset_name,scoreMethod,metacol_sampletype,SampleTypeSelected,Feature,subFeature,Pval_Tab,HR_Tab)
      })
      
      SplotBIN_react <- reactive({
        req(MedianCutP_react())
        meta_ssgsea_sdf <- MedianCutP_react()
        geneset <- gs_react()
        geneset_name <- names(geneset)
        SurvFeature <- paste0(geneset_name,"_MedianCutP")
        Feature <- input$FeatureSelection
        subFeature <- input$subFeatureSelection
        SampleTypeSelected <- input$SampleTypeChoices
        scoreMethod <- input$ScoreMethod
        show_pval <- input$ShowPval
        ShowConfInt <- input$ShowConfInt
        if (input$SurvYearOrMonth == "Years") {
          xBreaks <- input$SurvXaxisBreaks * 365.25
        } else {
          xBreaks <- input$SurvXaxisBreaks * 30.4375
        }
        if (isTruthy(input$SurvXaxis)) {
          if (input$SurvYearOrMonth == "Years") {
            xaxlim <- input$SurvXaxis * 365.25
            #meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,2] <= xaxlim),]
            #xBreaks <- input$SurvXaxisBreaks * 365.25
          } else {
            #xBreaks <- input$SurvXaxisBreaks * 30.4375
            xaxlim <- input$SurvXaxis * 30.4375
            #meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,2] <= xaxlim),]
          }
        } else {
          xaxlim <- NULL
          #xBreaks <- 365.25
        }
        surv_time_col <- input$SurvivalType_time
        showLegend <- input$SurvLegendPos
        showMedSurv <- input$ShowMedSurvLine
        if (showMedSurv == T) {
          showMedSurv <- "hv"
        }
        else if (showMedSurv == F) {
          showMedSurv <- "none"
        }
        SurvDateType <- sub("\\..*","",surv_time_col)
        if (input$GeneSetTabs == 2) {
          scoreMethodLab <- "Gene Expression"
        } else if (input$GeneSetTabs == 1 | input$GeneSetTabs == 3) {
          if (input$GeneSetCat_Select == "Pre-Processed Scores") {
            scoreMethodLab <- "Pre-Processed score"
          } else {
            scoreMethodLab <- paste(scoreMethod, " score", sep = "")
          }
        }
        
        PlotTitle <- SurvPlotTitle(SampleTypeSelected,geneset_name,scoreMethodLab,Feature,subFeature,"Median Cut-Point")
        #SurvFeature <- sprintf(ifelse((grepl(" ", SurvFeature) | !is.na(suppressWarnings(as.numeric(substring(SurvFeature, 1, 1))))), "`%s`", "%s"), SurvFeature)
        #save(list = ls(), file = "SurvPlotEnv.RData", envir = environment())
        #form <- as.formula(paste0("Surv(time,ID) ~ ",SurvFeature))
        #colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == SurvFeature)] <- sprintf(ifelse((grepl(" ", SurvFeature) | !is.na(suppressWarnings(as.numeric(substring(SurvFeature, 1, 1))))), "`%s`", "%s"), SurvFeature)
        #fit <- eval(substitute(survfit(form,data = meta_ssgsea_sdf, type="kaplan-meier")))
        
        
        colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == SurvFeature)] <- "Feature"
        form <- paste0("Surv(time,ID) ~ Feature")
        fit <- eval(substitute(survfit(as.formula(form),data = meta_ssgsea_sdf, type="kaplan-meier")))
        names(fit[["strata"]]) <- gsub("^Feature=",paste0(SurvFeature,"="),names(fit[["strata"]]))
        
        
        
        
        SurvPlot(fit,meta_ssgsea_sdf,PlotTitle,ylab = paste(SurvDateType,"Survival Probability"),
                 pval = show_pval,conf = ShowConfInt,legend = showLegend,median = showMedSurv,xlim = xaxlim,
                 xScale = input$SurvYearOrMonth, xBreaks = xBreaks)
      })
      output$SplotBIN <- renderPlot({
        plot <- SplotBIN_react()
        plot
      })
      
      ssgseaBINDensity_react <- reactive({
        
        req(ssGSEAmeta())
        geneset <- gs_react()
        geneset_name <- names(geneset)
        ssGSEAmeta <- ssGSEAmeta()
        scoreMethod <- input$ScoreMethod
        decon_score_cols <- decon_score_cols()
        CutPlabel <- "(Median Cut-Point)"
        
        SampleNameCol <- colnames(ssGSEAmeta)[1]
        ssgsea_scores <- ssGSEAmeta[,c(SampleNameCol,geneset_name)]
        quant_df <- data.frame(quantile(as.numeric(ssgsea_scores[,geneset_name]),probs = 0.5,na.rm = T))
        colnames(quant_df)[1] <- "Quantile"
        quant_df$Quantile <- round(quant_df$Quantile,3)
        
        ## get score method for x and y labels
        if (input$GeneSetTabs == 2) {
          ylab <- "Gene Expression Density"
          xlab <- "Gene Expression (Log(exp+1))"
          ssgsea_scores[,geneset_name] <- log(ssgsea_scores[,geneset_name] + 1)
          quant_df$Quantile <- round(log(quant_df$Quantile + 1),3)
        } else if (input$GeneSetTabs != 2) {
          if (geneset_name %in% decon_score_cols) {
            ylab <- "Pre-Processed Score Density"
            xlab <- "Pre-Processed Score"
          } else {
            ylab <- paste(scoreMethod, " Score Density", sep = "")
            xlab <- paste(scoreMethod, " Score", sep = "")
          }
        }
        ## generate title based on input
        if (geneset_name %in% decon_score_cols) {
          dens_title <- paste(geneset_name," Pre-Processed Score Density",sep = "")
        } else {
          if (input$GeneSetTabs == 2) {
            dens_title <- paste(geneset_name,"Log(exp+1)",ylab)
          } else {
            dens_title <- paste(geneset_name,ylab)
          }
        }
        
        densPlot(ssgsea_scores[,geneset_name,drop = F],quant_df,xlab,ylab,dens_title,CutPlabel)
        
      })
      output$ssgseaBINDensity <- renderPlot({
        
        plot <- ssgseaBINDensity_react()
        plot
        
      })
      
      output$MedianCutPSumm <- renderPrint({
        req(MedianCutPTab_react())
        CoxPHsumm(MedianCutPTab_react())
      })
      
      output$SBinaryHRtab <- renderTable({
        SBinaryHRtab_react()
      })
      
      ## Quartile CutP ---------------------------------------------------------
      
      QuartileCutP_react <- reactive({
        req(ssGSEAmeta())
        geneset <- gs_react()
        geneset_name <- names(geneset)
        SubsetSurvData(ssGSEAmeta(),input$SurvivalType_time,input$SurvivalType_id,paste0(geneset_name,"_QuartileCutP"))
      })
      
      QuartileCutPTab_react <- reactive({
        req(QuartileCutP_react())
        geneset <- gs_react()
        geneset_name <- names(geneset)
        CoxPHobj(QuartileCutP_react(),paste0(geneset_name,"_QuartileCutP"),"Q1_Low")
      })
      
      SQuartileHRtab_react <- reactive({
        req(QuartileCutPTab_react())
        CoxPHtabUni(QuartileCutPTab_react())
      })
      
      output$QuartSurvDescrip <- renderUI({
        HR_Tab <- SQuartileHRtab_react()
        Pval_Tab <- QuartileCutPTab_react()
        metacol_sampletype <- metacol_sampletype()
        SampleTypeSelected <- input$SampleTypeChoices
        Feature <- input$FeatureSelection
        subFeature <- input$subFeatureSelection
        geneset <- gs_react()
        geneset_name <- names(geneset)
        decon_score_cols <- decon_score_cols()
        surv_time_col <- input$SurvivalType_time
        scoreMethod <- input$ScoreMethod
        
        if (geneset_name %in% decon_score_cols) {
          scoreMethod <- "Pre-Processed Score"
        }
        if (input$GeneSetTabs == 2) {
          scoreMethod <- "Gene Expression"
        }
        SurvPlotExpl("quartile cut-point",surv_time_col,geneset_name,scoreMethod,metacol_sampletype,SampleTypeSelected,Feature,subFeature,Pval_Tab,HR_Tab)
      })
      
      Splot_react <- reactive({
        req(QuartileCutP_react())
        meta_ssgsea_sdf <- QuartileCutP_react()
        geneset <- gs_react()
        geneset_name <- names(geneset)
        SurvFeature <- paste0(geneset_name,"_QuartileCutP")
        Feature <- input$FeatureSelection
        subFeature <- input$subFeatureSelection
        SampleTypeSelected <- input$SampleTypeChoices
        scoreMethod <- input$ScoreMethod
        show_pval <- input$ShowPval
        ShowConfInt <- input$ShowConfInt
        if (input$SurvYearOrMonth == "Years") {
          xBreaks <- input$SurvXaxisBreaks * 365.25
        } else {
          xBreaks <- input$SurvXaxisBreaks * 30.4375
        }
        if (isTruthy(input$SurvXaxis)) {
          if (input$SurvYearOrMonth == "Years") {
            xaxlim <- input$SurvXaxis * 365.25
            #meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,2] <= xaxlim),]
            #xBreaks <- input$SurvXaxisBreaks * 365.25
          } else {
            #xBreaks <- input$SurvXaxisBreaks * 30.4375
            xaxlim <- input$SurvXaxis * 30.4375
            #meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,2] <= xaxlim),]
          }
        } else {
          xaxlim <- NULL
          #xBreaks <- 365.25
        }
        surv_time_col <- input$SurvivalType_time
        showLegend <- input$SurvLegendPos
        showMedSurv <- input$ShowMedSurvLine
        if (showMedSurv == T) {
          showMedSurv <- "hv"
        }
        else if (showMedSurv == F) {
          showMedSurv <- "none"
        }
        SurvDateType <- sub("\\..*","",surv_time_col)
        if (input$GeneSetTabs == 2) {
          scoreMethodLab <- "Gene Expression"
        } else if (input$GeneSetTabs == 1 | input$GeneSetTabs == 3) {
          if (input$GeneSetCat_Select == "Pre-Processed Scores") {
            scoreMethodLab <- "Pre-Processed score"
          } else {
            scoreMethodLab <- paste(scoreMethod, " score", sep = "")
          }
        }
        
        PlotTitle <- SurvPlotTitle(SampleTypeSelected,geneset_name,scoreMethodLab,Feature,subFeature,"Quartile Cut-Point")
        
        SurvFeature <- sprintf(ifelse((grepl(" ", SurvFeature) | !is.na(suppressWarnings(as.numeric(substring(SurvFeature, 1, 1))))), "`%s`", "%s"), SurvFeature)
        form <- as.formula(paste0("Surv(time,ID) ~ ",SurvFeature))
        fit <- eval(substitute(survfit(form,data = meta_ssgsea_sdf, type="kaplan-meier")))
        
        SurvPlot(fit,meta_ssgsea_sdf,PlotTitle,ylab = paste(SurvDateType,"Survival Probability"),
                 pval = show_pval,conf = ShowConfInt,legend = showLegend,median = showMedSurv,xlim = xaxlim,
                 xScale = input$SurvYearOrMonth, xBreaks = xBreaks)
      })
      output$Splot <- renderPlot({
        Splot_react()
      })
      ssgseaQuartDensity_react <- reactive({
        req(decon_score_cols())
        req(input$ScoreMethod)
        req(ssGSEAmeta())
        geneset <- gs_react()
        geneset_name <- names(geneset)
        ssGSEAmeta <- ssGSEAmeta()
        scoreMethod <- input$ScoreMethod
        decon_score_cols <- decon_score_cols()
        CutPlabel <- c("(Q1 Cut-Point)","(Q2 Cut-Point)","(Q3 Cut-Point)")
        
        SampleNameCol <- colnames(ssGSEAmeta)[1]
        ssgsea_scores <- ssGSEAmeta[,c(SampleNameCol,geneset_name)]
        
        quant_df <- data.frame(quantile(as.numeric(ssgsea_scores[,geneset_name]),na.rm = T))
        quant_df$quantiles <- as.numeric(gsub("%","",rownames(quant_df)))
        quant_df <- quant_df[which(quant_df$quantiles > 0 & quant_df$quantiles < 100),]
        quant_df <- quant_df[,which(colnames(quant_df) != "quantiles"),drop = F]
        colnames(quant_df)[1] <- "Quantile"
        quant_df$Quantile <- round(quant_df$Quantile,3)
        
        ## get score method for x and y labels
        if (input$GeneSetTabs == 2) {
          ylab <- "Gene Expression Density"
          xlab <- "Gene Expression (Log(exp+1))"
          ssgsea_scores[,geneset_name] <- log(ssgsea_scores[,geneset_name] + 1)
          quant_df$Quantile <- round(log(quant_df$Quantile + 1),3)
        } else if (input$GeneSetTabs != 2) {
          if (geneset_name %in% decon_score_cols) {
            ylab <- "Pre-Processed Score Density"
            xlab <- "Pre-Processed Score"
          } else {
            ylab <- paste(scoreMethod, " Score Density", sep = "")
            xlab <- paste(scoreMethod, " Score", sep = "")
          }
        }
        ## generate title based on input
        if (geneset_name %in% decon_score_cols) {
          dens_title <- paste(geneset_name," Pre-Processed Score Density",sep = "")
        } else {
          if (input$GeneSetTabs == 2) {
            dens_title <- paste(geneset_name,"Log(exp+1)",ylab)
          } else {
            dens_title <- paste(geneset_name,ylab)
          }
        }
        
        densPlot(ssgsea_scores[,geneset_name,drop = F],quant_df,xlab,ylab,dens_title,CutPlabel)
        
      })
      output$ssgseaQuartDensity <- renderPlot({
        ssgseaQuartDensity_react()
      })
      
      output$QuartileCutPSumm <- renderPrint({
        req(QuartileCutPTab_react())
        CoxPHsumm(QuartileCutPTab_react())
      })
      
      output$SQuartileHRtab <- renderTable({
        SQuartileHRtab_react()
      })
      
      ## Optimal CutP ----------------------------------------------------------
      
      OptimalCutP_react <- reactive({
        req(ssGSEAmeta())
        geneset <- gs_react()
        geneset_name <- names(geneset)
        SubsetSurvData(ssGSEAmeta(),input$SurvivalType_time,input$SurvivalType_id,paste0(geneset_name,"_OptimalCutP"))
      })
      
      OptimalCutPTab_react <- reactive({
        req(OptimalCutP_react())
        geneset <- gs_react()
        geneset_name <- names(geneset)
        CoxPHobj(OptimalCutP_react(),paste0(geneset_name,"_OptimalCutP"),"low")
      })
      
      CutPointHRtab_react <- reactive({
        req(OptimalCutPTab_react())
        CoxPHtabUni(OptimalCutPTab_react())
      })
      
      output$CutPSurvDescrip <- renderUI({
        HR_Tab <- CutPointHRtab_react()
        Pval_Tab <- OptimalCutPTab_react()
        metacol_sampletype <- metacol_sampletype()
        SampleTypeSelected <- input$SampleTypeChoices
        Feature <- input$FeatureSelection
        subFeature <- input$subFeatureSelection
        geneset <- gs_react()
        geneset_name <- names(geneset)
        decon_score_cols <- decon_score_cols()
        surv_time_col <- input$SurvivalType_time
        scoreMethod <- input$ScoreMethod
        
        if (geneset_name %in% decon_score_cols) {
          scoreMethod <- "Pre-Processed Score"
        }
        if (input$GeneSetTabs == 2) {
          scoreMethod <- "Gene Expression"
        }
        SurvPlotExpl("optimal cut-point",surv_time_col,geneset_name,scoreMethod,metacol_sampletype,SampleTypeSelected,Feature,subFeature,Pval_Tab,HR_Tab)
      })
      
      ScutPointPlot_react <- reactive({
        req(OptimalCutP_react())
        meta_ssgsea_sdf <- OptimalCutP_react()
        geneset <- gs_react()
        geneset_name <- names(geneset)
        SurvFeature <- paste0(geneset_name,"_OptimalCutP")
        Feature <- input$FeatureSelection
        subFeature <- input$subFeatureSelection
        SampleTypeSelected <- input$SampleTypeChoices
        scoreMethod <- input$ScoreMethod
        show_pval <- input$ShowPval
        ShowConfInt <- input$ShowConfInt
        if (input$SurvYearOrMonth == "Years") {
          xBreaks <- input$SurvXaxisBreaks * 365.25
        } else {
          xBreaks <- input$SurvXaxisBreaks * 30.4375
        }
        if (isTruthy(input$SurvXaxis)) {
          if (input$SurvYearOrMonth == "Years") {
            xaxlim <- input$SurvXaxis * 365.25
            #meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,2] <= xaxlim),]
            #xBreaks <- input$SurvXaxisBreaks * 365.25
          } else {
            #xBreaks <- input$SurvXaxisBreaks * 30.4375
            xaxlim <- input$SurvXaxis * 30.4375
            #meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,2] <= xaxlim),]
          }
        } else {
          xaxlim <- NULL
          #xBreaks <- 365.25
        }
        surv_time_col <- input$SurvivalType_time
        showLegend <- input$SurvLegendPos
        showMedSurv <- input$ShowMedSurvLine
        if (showMedSurv == T) {
          showMedSurv <- "hv"
        }
        else if (showMedSurv == F) {
          showMedSurv <- "none"
        }
        SurvDateType <- sub("\\..*","",surv_time_col)
        if (input$GeneSetTabs == 2) {
          scoreMethodLab <- "Gene Expression"
        } else if (input$GeneSetTabs == 1 | input$GeneSetTabs == 3) {
          if (input$GeneSetCat_Select == "Pre-Processed Scores") {
            scoreMethodLab <- "Pre-Processed score"
          } else {
            scoreMethodLab <- paste(scoreMethod, " score", sep = "")
          }
        }
        
        PlotTitle <- SurvPlotTitle(SampleTypeSelected,geneset_name,scoreMethodLab,Feature,subFeature,"Optimal Cut-Point")
        
        SurvFeature <- sprintf(ifelse((grepl(" ", SurvFeature) | !is.na(suppressWarnings(as.numeric(substring(SurvFeature, 1, 1))))), "`%s`", "%s"), SurvFeature)
        form <- as.formula(paste0("Surv(time,ID) ~ ",SurvFeature))
        fit <- eval(substitute(survfit(form,data = meta_ssgsea_sdf, type="kaplan-meier")))
        
        SurvPlot(fit,meta_ssgsea_sdf,PlotTitle,ylab = paste(SurvDateType,"Survival Probability"),
                 pval = show_pval,conf = ShowConfInt,legend = showLegend,median = showMedSurv,xlim = xaxlim,
                 xScale = input$SurvYearOrMonth, xBreaks = xBreaks)
      })
      output$ScutPointPlot <- renderPlot({
        ScutPointPlot_react()
      })
      ssgseaCutPDensity_react <- reactive({
        
        req(decon_score_cols())
        req(input$ScoreMethod)
        req(ssGSEAmeta())
        geneset <- gs_react()
        geneset_name <- names(geneset)
        ssGSEAmeta <- ssGSEAmeta()
        scoreMethod <- input$ScoreMethod
        decon_score_cols <- decon_score_cols()
        CutPlabel <- "(Optimal Cut-Point)"
        
        SampleNameCol <- colnames(ssGSEAmeta)[1]
        ssgsea_scores <- ssGSEAmeta[,c(SampleNameCol,geneset_name)]
        meta_ssgsea_sdf <- merge(OptimalCutP_react(),ssgsea_scores)
        
        if (length(meta_ssgsea_sdf[,4][meta_ssgsea_sdf[,4] > 0])/length(meta_ssgsea_sdf[,4]) > 0.01) {
          res.cut <- survminer::surv_cutpoint(meta_ssgsea_sdf,time = "time", event = "ID", variable = geneset_name, minprop = 0.01)
          cutp <- res.cut$cutpoint[["cutpoint"]]
        }
        
        quant_df <- data.frame(Quantile = cutp)
        quant_df$Quantile <- round(quant_df$Quantile,3)
        
        ## get score method for x and y labels
        if (input$GeneSetTabs == 2) {
          ylab <- "Gene Expression Density"
          xlab <- "Gene Expression (Log(exp+1))"
          ssgsea_scores[,geneset_name] <- log(ssgsea_scores[,geneset_name] + 1)
          quant_df$Quantile <- round(log(quant_df$Quantile + 1),3)
        } else if (input$GeneSetTabs != 2) {
          if (geneset_name %in% decon_score_cols) {
            ylab <- "Pre-Processed Score Density"
            xlab <- "Pre-Processed Score"
          } else {
            ylab <- paste(scoreMethod, " Score Density", sep = "")
            xlab <- paste(scoreMethod, " Score", sep = "")
          }
        }
        ## generate title based on input
        if (geneset_name %in% decon_score_cols) {
          dens_title <- paste(geneset_name," Pre-Processed Score Density",sep = "")
        } else {
          if (input$GeneSetTabs == 2) {
            dens_title <- paste(geneset_name,"Log(exp+1)",ylab)
          } else {
            dens_title <- paste(geneset_name,ylab)
          }
        }
        
        densPlot(ssgsea_scores[,geneset_name,drop = F],quant_df,xlab,ylab,dens_title,CutPlabel)
        
      })
      output$ssgseaCutPDensity <- renderPlot({
        ssgseaCutPDensity_react()
      })
      
      output$OptimalCutPSumm <- renderPrint({
        req(OptimalCutPTab_react())
        CoxPHsumm(OptimalCutPTab_react())
      })
      
      output$CutPointHRtab <- renderTable({
        CutPointHRtab_react()
      })
      
      ## Top/Bottom CutP -------------------------------------------------------
      
      TopBottomCutP_react <- reactive({
        req(ssGSEAmeta())
        geneset <- gs_react()
        geneset_name <- names(geneset)
        meta_ssgsea_sdf <- SubsetSurvData(ssGSEAmeta(),input$SurvivalType_time,input$SurvivalType_id,paste0(geneset_name,"_TopBottomCutP"))
        meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,paste0(geneset_name,"_TopBottomCutP")] != "BetweenCutoff"),]
        meta_ssgsea_sdf
      })
      
      TopBottomCutPTab_react <- reactive({
        req(TopBottomCutP_react())
        geneset <- gs_react()
        geneset_name <- names(geneset)
        CoxPHobj(TopBottomCutP_react(),paste0(geneset_name,"_TopBottomCutP"),"Low")
      })
      
      SQuantileHRtab_react <- reactive({
        req(TopBottomCutPTab_react())
        CoxPHtabUni(TopBottomCutPTab_react())
      })
      
      output$QuantSurvDescrip <- renderUI({
        HR_Tab <- SQuantileHRtab_react()
        Pval_Tab <- TopBottomCutPTab_react()
        metacol_sampletype <- metacol_sampletype()
        SampleTypeSelected <- input$SampleTypeChoices
        Feature <- input$FeatureSelection
        subFeature <- input$subFeatureSelection
        geneset <- gs_react()
        geneset_name <- names(geneset)
        decon_score_cols <- decon_score_cols()
        surv_time_col <- input$SurvivalType_time
        scoreMethod <- input$ScoreMethod
        
        if (geneset_name %in% decon_score_cols) {
          scoreMethod <- "Pre-Processed Score"
        }
        if (input$GeneSetTabs == 2) {
          scoreMethod <- "Gene Expression"
        }
        SurvPlotExpl("top/bottom quantile cut-point",surv_time_col,geneset_name,scoreMethod,metacol_sampletype,SampleTypeSelected,Feature,subFeature,Pval_Tab,HR_Tab)
      })
      
      SquantPlot_react <- reactive({
        req(TopBottomCutP_react())
        meta_ssgsea_sdf <- TopBottomCutP_react()
        geneset <- gs_react()
        geneset_name <- names(geneset)
        SurvFeature <- paste0(geneset_name,"_TopBottomCutP")
        Feature <- input$FeatureSelection
        subFeature <- input$subFeatureSelection
        SampleTypeSelected <- input$SampleTypeChoices
        scoreMethod <- input$ScoreMethod
        show_pval <- input$ShowPval
        ShowConfInt <- input$ShowConfInt
        if (input$SurvYearOrMonth == "Years") {
          xBreaks <- input$SurvXaxisBreaks * 365.25
        } else {
          xBreaks <- input$SurvXaxisBreaks * 30.4375
        }
        if (isTruthy(input$SurvXaxis)) {
          if (input$SurvYearOrMonth == "Years") {
            xaxlim <- input$SurvXaxis * 365.25
            #meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,2] <= xaxlim),]
            #xBreaks <- input$SurvXaxisBreaks * 365.25
          } else {
            #xBreaks <- input$SurvXaxisBreaks * 30.4375
            xaxlim <- input$SurvXaxis * 30.4375
            #meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,2] <= xaxlim),]
          }
        } else {
          xaxlim <- NULL
          #xBreaks <- 365.25
        }
        surv_time_col <- input$SurvivalType_time
        showLegend <- input$SurvLegendPos
        showMedSurv <- input$ShowMedSurvLine
        if (showMedSurv == T) {
          showMedSurv <- "hv"
        }
        else if (showMedSurv == F) {
          showMedSurv <- "none"
        }
        SurvDateType <- sub("\\..*","",surv_time_col)
        if (input$GeneSetTabs == 2) {
          scoreMethodLab <- "Gene Expression"
        } else if (input$GeneSetTabs == 1 | input$GeneSetTabs == 3) {
          if (input$GeneSetCat_Select == "Pre-Processed Scores") {
            scoreMethodLab <- "Pre-Processed score"
          } else {
            scoreMethodLab <- paste(scoreMethod, " score", sep = "")
          }
        }
        
        PlotTitle <- SurvPlotTitle(SampleTypeSelected,geneset_name,scoreMethodLab,Feature,subFeature,"top/bottom quantile cut-points")
        
        SurvFeature <- sprintf(ifelse((grepl(" ", SurvFeature) | !is.na(suppressWarnings(as.numeric(substring(SurvFeature, 1, 1))))), "`%s`", "%s"), SurvFeature)
        form <- as.formula(paste0("Surv(time,ID) ~ ",SurvFeature))
        fit <- eval(substitute(survfit(form,data = meta_ssgsea_sdf, type="kaplan-meier")))
        
        SurvPlot(fit,meta_ssgsea_sdf,PlotTitle,ylab = paste(SurvDateType,"Survival Probability"),
                 pval = show_pval,conf = ShowConfInt,legend = showLegend,median = showMedSurv,xlim = xaxlim,
                 xScale = input$SurvYearOrMonth, xBreaks = xBreaks)
      })
      output$SquantPlot <- renderPlot({
        SquantPlot_react()
      })
      ssgseaQuantDensity_react <- reactive({
        
        req(decon_score_cols())
        req(input$ScoreMethod)
        req(ssGSEAmeta())
        geneset <- gs_react()
        geneset_name <- names(geneset)
        ssGSEAmeta <- ssGSEAmeta()
        scoreMethod <- input$ScoreMethod
        decon_score_cols <- decon_score_cols()
        quantCutoff <- input$QuantPercent/100
        CutPlabel1 <- "(Bottom Cut-Point)"
        CutPlabel2 <- "(Top Cut-Point)"
        CutPlabel <- c(CutPlabel1,CutPlabel2)
        
        SampleNameCol <- colnames(ssGSEAmeta)[1]
        ssgsea_scores <- ssGSEAmeta[,c(SampleNameCol,geneset_name)]
        
        cutp_high <- round(quantile(ssgsea_scores[,geneset_name],1-quantCutoff,na.rm = T),3)
        cutp_low <- round(quantile(ssgsea_scores[,geneset_name],quantCutoff,na.rm = T),3)
        xints <- c(cutp_low,cutp_high)
        
        quant_df <- data.frame(Quantile = xints)
        quant_df$Quantile <- round(quant_df$Quantile,3)
        
        ## get score method for x and y labels
        if (input$GeneSetTabs == 2) {
          ylab <- "Gene Expression Density"
          xlab <- "Gene Expression (Log(exp+1))"
          ssgsea_scores[,geneset_name] <- log(ssgsea_scores[,geneset_name] + 1)
          quant_df$Quantile <- round(log(quant_df$Quantile + 1),3)
        } else if (input$GeneSetTabs != 2) {
          if (geneset_name %in% decon_score_cols) {
            ylab <- "Pre-Processed Score Density"
            xlab <- "Pre-Processed Score"
          } else {
            ylab <- paste(scoreMethod, " Score Density", sep = "")
            xlab <- paste(scoreMethod, " Score", sep = "")
          }
        }
        ## generate title based on input
        if (geneset_name %in% decon_score_cols) {
          dens_title <- paste(geneset_name," Pre-Processed Score Density",sep = "")
        } else {
          if (input$GeneSetTabs == 2) {
            dens_title <- paste(geneset_name,"Log(exp+1)",ylab)
          } else {
            dens_title <- paste(geneset_name,ylab)
          }
        }
        
        densPlot(ssgsea_scores[,geneset_name,drop = F],quant_df,xlab,ylab,dens_title,CutPlabel)
        
      })
      output$ssgseaQuantDensity <- renderPlot({
        ssgseaQuantDensity_react()
      })
      
      output$QuantileCutPSumm <- renderPrint({
        req(TopBottomCutPTab_react())
        CoxPHsumm(TopBottomCutPTab_react())
      })
      
      output$SQuantileHRtab <- renderTable({
        SQuantileHRtab_react()
      })
      
      ## User CutP -------------------------------------------------------
      
      UserCutP_react <- reactive({
        req(ssGSEAmeta())
        geneset <- gs_react()
        geneset_name <- names(geneset)
        meta_ssgsea_sdf <- SubsetSurvData(ssGSEAmeta(),input$SurvivalType_time,input$SurvivalType_id,paste0(geneset_name,"_UserCutP"))
        meta_ssgsea_sdf
      })
      
      UserCutPTab_react <- reactive({
        req(UserCutP_react())
        geneset <- gs_react()
        geneset_name <- names(geneset)
        CoxPHobj(UserCutP_react(),paste0(geneset_name,"_UserCutP"),"Low")
      })
      
      SQuantileHR2tab_react <- reactive({
        req(UserCutPTab_react())
        CoxPHtabUni(UserCutPTab_react())
      })
      
      output$Quant2SurvDescrip <- renderUI({
        HR_Tab <- SQuantileHR2tab_react()
        Pval_Tab <- UserCutPTab_react()
        metacol_sampletype <- metacol_sampletype()
        SampleTypeSelected <- input$SampleTypeChoices
        Feature <- input$FeatureSelection
        subFeature <- input$subFeatureSelection
        geneset <- gs_react()
        geneset_name <- names(geneset)
        decon_score_cols <- decon_score_cols()
        surv_time_col <- input$SurvivalType_time
        scoreMethod <- input$ScoreMethod
        
        if (geneset_name %in% decon_score_cols) {
          scoreMethod <- "Pre-Processed Score"
        }
        if (input$GeneSetTabs == 2) {
          scoreMethod <- "Gene Expression"
        }
        SurvPlotExpl("user cut-point",surv_time_col,geneset_name,scoreMethod,metacol_sampletype,SampleTypeSelected,Feature,subFeature,Pval_Tab,HR_Tab)
      })
      
      SquantPlot2_react <- reactive({
        req(UserCutP_react())
        meta_ssgsea_sdf <- UserCutP_react()
        geneset <- gs_react()
        geneset_name <- names(geneset)
        SurvFeature <- paste0(geneset_name,"_UserCutP")
        Feature <- input$FeatureSelection
        subFeature <- input$subFeatureSelection
        SampleTypeSelected <- input$SampleTypeChoices
        scoreMethod <- input$ScoreMethod
        show_pval <- input$ShowPval
        ShowConfInt <- input$ShowConfInt
        if (input$SurvYearOrMonth == "Years") {
          xBreaks <- input$SurvXaxisBreaks * 365.25
        } else {
          xBreaks <- input$SurvXaxisBreaks * 30.4375
        }
        if (isTruthy(input$SurvXaxis)) {
          if (input$SurvYearOrMonth == "Years") {
            xaxlim <- input$SurvXaxis * 365.25
            #meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,2] <= xaxlim),]
            #xBreaks <- input$SurvXaxisBreaks * 365.25
          } else {
            #xBreaks <- input$SurvXaxisBreaks * 30.4375
            xaxlim <- input$SurvXaxis * 30.4375
            #meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,2] <= xaxlim),]
          }
        } else {
          xaxlim <- NULL
          #xBreaks <- 365.25
        }
        surv_time_col <- input$SurvivalType_time
        showLegend <- input$SurvLegendPos
        showMedSurv <- input$ShowMedSurvLine
        if (showMedSurv == T) {
          showMedSurv <- "hv"
        }
        else if (showMedSurv == F) {
          showMedSurv <- "none"
        }
        SurvDateType <- sub("\\..*","",surv_time_col)
        if (input$GeneSetTabs == 2) {
          scoreMethodLab <- "Gene Expression"
        } else if (input$GeneSetTabs == 1 | input$GeneSetTabs == 3) {
          if (input$GeneSetCat_Select == "Pre-Processed Scores") {
            scoreMethodLab <- "Pre-Processed score"
          } else {
            scoreMethodLab <- paste(scoreMethod, " score", sep = "")
          }
        }
        
        PlotTitle <- SurvPlotTitle(SampleTypeSelected,geneset_name,scoreMethodLab,Feature,subFeature,"user cut-point")
        
        SurvFeature <- sprintf(ifelse((grepl(" ", SurvFeature) | !is.na(suppressWarnings(as.numeric(substring(SurvFeature, 1, 1))))), "`%s`", "%s"), SurvFeature)
        form <- as.formula(paste0("Surv(time,ID) ~ ",SurvFeature))
        fit <- eval(substitute(survfit(form,data = meta_ssgsea_sdf, type="kaplan-meier")))
        
        SurvPlot(fit,meta_ssgsea_sdf,PlotTitle,ylab = paste(SurvDateType,"Survival Probability"),
                 pval = show_pval,conf = ShowConfInt,legend = showLegend,median = showMedSurv,xlim = xaxlim,
                 xScale = input$SurvYearOrMonth, xBreaks = xBreaks)
      })
      output$SquantPlot2 <- renderPlot({
        SquantPlot2_react()
      })
      ssgseaQuant2Density_react <- reactive({
        
        req(decon_score_cols())
        req(input$ScoreMethod)
        req(ssGSEAmeta())
        geneset <- gs_react()
        geneset_name <- names(geneset)
        ssGSEAmeta <- ssGSEAmeta()
        scoreMethod <- input$ScoreMethod
        decon_score_cols <- decon_score_cols()
        quantCutoff <- input$QuantPercent2/100
        CutPlabel <- "(User Cut-Point)"
        
        SampleNameCol <- colnames(ssGSEAmeta)[1]
        ssgsea_scores <- ssGSEAmeta[,c(SampleNameCol,geneset_name)]
        
        cutp_user <- round(quantile(ssgsea_scores[,geneset_name],quantCutoff,na.rm = T),3)
        
        quant_df <- data.frame(Quantile = cutp_user)
        quant_df$Quantile <- round(quant_df$Quantile,3)
        
        ## get score method for x and y labels
        if (input$GeneSetTabs == 2) {
          ylab <- "Gene Expression Density"
          xlab <- "Gene Expression (Log(exp+1))"
          ssgsea_scores[,geneset_name] <- log(ssgsea_scores[,geneset_name] + 1)
          quant_df$Quantile <- round(log(quant_df$Quantile + 1),3)
        } else if (input$GeneSetTabs != 2) {
          if (geneset_name %in% decon_score_cols) {
            ylab <- "Pre-Processed Score Density"
            xlab <- "Pre-Processed Score"
          } else {
            ylab <- paste(scoreMethod, " Score Density", sep = "")
            xlab <- paste(scoreMethod, " Score", sep = "")
          }
        }
        ## generate title based on input
        if (geneset_name %in% decon_score_cols) {
          dens_title <- paste(geneset_name," Pre-Processed Score Density",sep = "")
        } else {
          if (input$GeneSetTabs == 2) {
            dens_title <- paste(geneset_name,"Log(exp+1)",ylab)
          } else {
            dens_title <- paste(geneset_name,ylab)
          }
        }
        
        densPlot(ssgsea_scores[,geneset_name,drop = F],quant_df,xlab,ylab,dens_title,CutPlabel)
        
      })
      output$ssgseaQuant2Density <- renderPlot({
        ssgseaQuant2Density_react()
      })
      
      output$UserCutPSumm <- renderPrint({
        req(UserCutPTab_react())
        CoxPHsumm(UserCutPTab_react())
      })
      
      output$SQuantileHR2tab <- renderTable({
        SQuantileHR2tab_react()
      })
      
      ## Univariate ----------------------------------------------------------------
      
      observe({
        geneset <- gs_react()
        geneset_name <- names(geneset)
        updateSelectizeInput(session = session,inputId = "SingleSurvivalFeature", choices = metacol_feature(), selected = paste0(geneset_name,"_MedianCutP"), server = T)
      })
      
      output$rendSurvFeatVariableUni <- renderUI({
        req(input$SingleSurvivalFeature)
        req(ssGSEAmeta())
        Feature <- input$SingleSurvivalFeature
        meta <- ssGSEAmeta()
        Var_choices <- survFeatRefSelect(meta,Feature,input$UniVarNAcheck,input$UniVarContCheck,input$UniVarContCheck)
        selectInput("SurvFeatVariableUni","Select Coxh Feature Reference:", choices = Var_choices, width = "80%")
      })
      
      UniVarFeat_react <- reactive({
        req(ssGSEAmeta())
        req(input$SurvivalType_time)
        req(input$SurvivalType_id)
        req(input$SingleSurvivalFeature)
        Feature <- input$SingleSurvivalFeature
        meta_ssgsea <- ssGSEAmeta()
        if (Feature %in% colnames(meta_ssgsea)) {
          if (input$UniVarNAcheck == TRUE) {
            meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature]) == FALSE),]
            meta_ssgsea <- meta_ssgsea[which(meta_ssgsea[,Feature] != "Inf" & meta_ssgsea[,Feature] != "N/A" & meta_ssgsea[,Feature] != "n/a"),]
            meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature],ignore.case = T, invert = T),]
          }
          if (input$UniVarContCheck == TRUE) {
            if (input$UniVarContHiLoCheck == TRUE) {
              meta_ssgsea[,Feature] <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == Feature)])
            }
          }
          meta_ssgsea_sdf <- SubsetSurvData(meta_ssgsea,input$SurvivalType_time,input$SurvivalType_id,Feature)
          meta_ssgsea_sdf
        }
        
      })
      
      UniVarFeatTab_react <- reactive({
        
        req(UniVarFeat_react())
        req(input$SingleSurvivalFeature)
        req(input$SurvFeatVariableUni)
        meta_ssgsea_sdf <- UniVarFeat_react()
        Feature <- input$SingleSurvivalFeature
        ref_Feature <- input$SurvFeatVariableUni
        #save(list = ls(), file = "UniVarFeatTab_react_env.RData", envir = environment())
        if ((input$UniVarContCheck == FALSE) | (input$UniVarContCheck == TRUE & input$UniVarContHiLoCheck == TRUE)) {
          if (ref_Feature %in% levels(factor(meta_ssgsea_sdf[,Feature]))) {
            tab <- CoxPHobj(meta_ssgsea_sdf,Feature,ref_Feature)
            tab
          }
        } else {
          colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == Feature)] <- gsub("-","_",Feature)
          Feature <- gsub("-","_",Feature)
          Feature <- sprintf(ifelse((grepl(" ", Feature) | !is.na(suppressWarnings(as.numeric(substring(Feature, 1, 1))))), "`%s`", "%s"), Feature)
          tab <- coxph(as.formula(paste0("Surv(time,ID) ~ ",Feature)),data = meta_ssgsea_sdf)
          tab
        }
        
      })
      
      SSingleFeatureHRtab_react <- reactive({
        
        #Feature <- input$SingleSurvivalFeature
        #meta_ssgsea_sdf <- UniVarFeat_react()
        Feature <- input$SingleSurvivalFeature
        #ref_Feature <- input$SurvFeatVariableUni
        
        tab <- UniVarFeatTab_react()
        #names(tab[["assign"]]) <- gsub("^Feature",Feature,names(tab[["assign"]]))
        tabOut <- CoxPHtabUni(tab,Feature)
        tabOut
      })
      
      featSplot_react <- reactive({
        
        ## Assign variables
        req(UniVarFeat_react())
        meta_ssgsea_sdf <- UniVarFeat_react()
        SampleType <- input$SampleTypeSelection
        Feature <- input$SingleSurvivalFeature
        Feature_sub <- input$FeatureSelection
        subFeature <- input$subFeatureSelection
        show_pval <- input$ShowPval
        ShowConfInt <- input$ShowConfInt
        if (input$SurvYearOrMonth == "Years") {
          xBreaks <- input$SurvXaxisBreaks * 365.25
        } else {
          xBreaks <- input$SurvXaxisBreaks * 30.4375
        }
        if (isTruthy(input$SurvXaxis)) {
          if (input$SurvYearOrMonth == "Years") {
            xaxlim <- input$SurvXaxis * 365.25
            #meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,2] <= xaxlim),]
            #xBreaks <- input$SurvXaxisBreaks * 365.25
          } else {
            #xBreaks <- input$SurvXaxisBreaks * 30.4375
            #meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,2] <= xaxlim),]
            xaxlim <- input$SurvXaxis * 30.4375
          }
        } else {
          xaxlim <- NULL
          #xBreaks <- 365.25
        }
        surv_time_col <- input$SurvivalType_time
        showLegend <- input$SurvLegendPos
        showMedSurv <- input$ShowMedSurvLine
        if (showMedSurv == T) {
          showMedSurv <- "hv"
        }
        else if (showMedSurv == F) {
          showMedSurv <- "none"
        }
        SurvDateType <- sub("\\..*","",surv_time_col)
        
        ## determine Feature and Sample Type label
        ## Adjust 'Sample Type' for label
        if (!is.null(SampleType)) {
          SampleTypeLab <- paste(" (",SampleType,") ",sep = "")
        } else { SampleTypeLab <- " " }
        
        ## Determine Plot title
        if (is.null(input$SurvPlotTitleUniVar)) {
          SurvPlotTitle <- paste("Survival curves of ",Feature,SampleTypeLab,"Patients", sep = "")
        }
        else if (!is.null(input$SurvPlotTitleUniVar)) {
          if (input$SurvPlotTitleUniVar == "") {
            SurvPlotTitle <- paste("Survival curves of ",Feature,SampleTypeLab,"Patients", sep = "")
          }
          else if (input$SurvPlotTitleUniVar != "") {
            SurvPlotTitle <- input$SurvPlotTitleUniVar
          }
        }
        
        
        PlotTitle <- SurvPlotTitle(SampleTypeSelected = SampleType,Feature = Feature_sub, subFeature = subFeature, univar = Feature)
        colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == Feature)] <- "Feature"
        form <- paste0("Surv(time,ID) ~ Feature")
        fit <- eval(substitute(survfit(as.formula(form),data = meta_ssgsea_sdf, type="kaplan-meier")))
        names(fit[["strata"]]) <- gsub("^Feature=",paste0(Feature,"="),names(fit[["strata"]]))
        
        SurvPlot(fit,meta_ssgsea_sdf,PlotTitle,ylab = paste(SurvDateType,"Survival Probability"),
                 pval = show_pval,conf = ShowConfInt,legend = showLegend,median = showMedSurv,xlim = xaxlim,
                 xScale = input$SurvYearOrMonth, xBreaks = xBreaks)
        
      })
      
      
      output$featSplot <- renderPlot({
        plot <- featSplot_react()
        plot
      })
      
      output$SSingleFeatureHRtab <- renderTable({
        
        tab <- SSingleFeatureHRtab_react()
        tab
        
      })
      
      output$UnivarSummary <- renderPrint({
        
        tab <- UniVarFeatTab_react()
        Feature <- input$SingleSurvivalFeature
        CoxPHsumm(tab)
        
      })
      
      SinglevarForestPlot_react <- reactive({
        req(UniVarFeatTab_react())
        req(UniVarFeat_react())
        req(input$SingleSurvivalFeature)
        req(input$ForestFontSize)
        obj <- UniVarFeatTab_react()
        Feature <- input$SingleSurvivalFeature
        meta_ssgsea_sdf <- UniVarFeat_react()
        forest <- forestPlot_Simple(obj,meta_ssgsea_sdf,input$SingleSurvivalFeature,input$ForestFontSize)
        forest
      })
      output$SinglevarForestPlot <- renderPlot({
        forest <- SinglevarForestPlot_react()
        forest
      })
      
      observe({
        req(ssGSEAmeta())
        meta <- ssGSEAmeta()
        geneset <- gs_react()
        geneset_name <- names(geneset)
        ColLevels <- apply(meta,2,function(x) length(levels(as.factor(x))))
        DicotCols <- names(ColLevels[ColLevels==2])
        ContCols <- GetColsOfType(meta,"continuous")
        FeatureChoices <- unique(c(DicotCols,ContCols))
        SurvDataCols <- c(metacol_survid,metacol_survtime)
        FeatureChoices <- FeatureChoices[which(!FeatureChoices %in% SurvDataCols)]
        updateSelectizeInput(session,"MultiFeatUnivarSelect",choices = FeatureChoices,selected = paste0(geneset_name,"_MedianCutP"), server = T)
      })
      
      
      MultiFeat_ForestMeta <- reactive({
        
        req(ssGSEAmeta())
        req(input$SurvivalType_time)
        req(input$SurvivalType_id)
        meta <- ssGSEAmeta()
        geneset <- gs_react()
        geneset_name <- names(geneset)
        surv_time_col <- input$SurvivalType_time
        surv_id_col <- input$SurvivalType_id
        colnames(meta)[1] <- "SampleName"
        FeatCols <- c("SampleName",surv_time_col,surv_id_col,input$MultiFeatUnivarSelect)
        metaSub <- meta[,FeatCols]
        continuousCols <- GetColsOfType(metaSub[,-c(1:3), drop = F],"continuous")
        metaSub[,continuousCols] <- apply(metaSub[,continuousCols, drop = F],2,function(x) highlow2(x))
        if (paste0(geneset_name,"_MedianCutP") %in% colnames(metaSub)) {
          metaSub[,paste0(geneset_name,"_MedianCutP")] <- as.factor(metaSub[,paste0(geneset_name,"_MedianCutP")])
          metaSub[,paste0(geneset_name,"_MedianCutP")] <- relevel(metaSub[,paste0(geneset_name,"_MedianCutP")], ref = "Low")
        }
        metaSub
        
      })
      output$univarForestPlotTable <- DT::renderDataTable({
        df <- MultiFeat_ForestMeta()
        DT::datatable(df,
                      extensions = "FixedColumns",
                      options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                     pageLength = 20,
                                     fixedColumns = list(leftColumns = 1),
                                     scrollX = T),
                      rownames = F)
        
      })
      plotht <- reactiveVal(300)
      
      MultiFeatUnivarForestPlotTab_react <- reactive({
        req(MultiFeat_ForestMeta())
        req(SBinaryHRtab_react())
        req(input$SurvivalType_time)
        req(input$SurvivalType_id)
        req(input$MultiFeatUnivarSelect)
        metaSub <- MultiFeat_ForestMeta()
        geneset <- gs_react()                 #Chosen Gene Set
        geneset_name <- names(geneset)        #Name of chosen gene set
        surv_time_col <- input$SurvivalType_time
        surv_id_col <- input$SurvivalType_id
        FeatCols <- c(input$MultiFeatUnivarSelect)
        HR_Tab <- SBinaryHRtab_react()             # Hazard Ratio Table Reactive
        HR <- HR_Tab[3,2]
        if (as.numeric(HR) > 1) {
          HighRisk <- TRUE
        } else {HighRisk <- FALSE}
        FeatColsList <- list()
        for (feat in FeatCols) {
          metaSubSet <- metaSub[,c("SampleName",surv_time_col,surv_id_col,feat)]
          metaSubSet <- metaSubSet[which(is.na(metaSubSet[,feat]) == FALSE),]
          metaSubSet <- metaSubSet[which(metaSubSet[,feat] != "Inf"),]
          metaSubSet <- metaSubSet[grep("unknown",metaSubSet[,feat],ignore.case = T, invert = T),]
          colnames(metaSubSet)[c(2:4)] <- c("time","status","feature")
          tab <- coxph(Surv(time,status) ~ feature, data = metaSubSet)
          tabData <- get_tabData(tab)
          tabData[1] <- feat
          colnames(metaSubSet)[4] <- feat
          if (HighRisk) {
            if (as.numeric(tabData[3]) > 1) {
              levelLabel <- paste0("(",levels(as.factor(metaSubSet[,feat]))[2],")")
              tabData[1] <- feat
              tabData[1] <- paste(tabData[1],levelLabel)
              FeatColsList[[feat]] <- tabData
            } else {
              metaSubSet[,feat]  <- metaSubSet[,feat] %>% as.factor() %>% forcats::fct_rev()
              colnames(metaSubSet)[c(2:4)] <- c("time","status","feature")
              tab <- coxph(Surv(time,status) ~ feature, data = metaSubSet)
              colnames(metaSubSet)[4] <- feat
              tabData <- get_tabData(tab)
              levelLabel <- paste0("(",levels(as.factor(metaSubSet[,feat]))[2],")")
              tabData[1] <- feat
              tabData[1] <- paste(tabData[1],levelLabel)
              FeatColsList[[feat]] <- tabData
            }
          } else {
            if (as.numeric(tabData[3]) < 1) {
              levelLabel <- paste0("(",levels(as.factor(metaSubSet[,feat]))[2],")")
              tabData[1] <- feat
              tabData[1] <- paste(tabData[1],levelLabel)
              FeatColsList[[feat]] <- tabData
            } else {
              metaSubSet[,feat]  <- metaSubSet[,feat] %>% as.factor() %>% forcats::fct_rev()
              colnames(metaSubSet)[c(2:4)] <- c("time","status","feature")
              tab <- coxph(Surv(time,status) ~ feature, data = metaSubSet)
              colnames(metaSubSet)[4] <- feat
              tabData <- get_tabData(tab)
              levelLabel <- paste0("(",levels(as.factor(metaSubSet[,feat]))[2],")")
              tabData[1] <- feat
              tabData[1] <- paste(tabData[1],levelLabel)
              FeatColsList[[feat]] <- tabData
            }
          }
        }
        FeatColsDF <- as.data.frame(do.call(rbind, FeatColsList))
        rownames(FeatColsDF) <- NULL
        FeatColsDF$Low <- as.numeric(FeatColsDF$Low)
        FeatColsDF$High <- as.numeric(FeatColsDF$High)
        FeatColsDF$`Hazard Ratio` <- as.numeric(FeatColsDF$`Hazard Ratio`)
        FeatColsDF$P.Value <- as.numeric(FeatColsDF$P.Value)
        FeatColsDF$` ` <- paste(rep(" ", 20), collapse = " ")
        FeatColsDF$`HR (95% CI)` <- ifelse(is.na(FeatColsDF$`Standard Error`), "",
                                           sprintf("%.2f (%.2f to %.2f)",
                                                   FeatColsDF$`Hazard Ratio`, FeatColsDF$Low, FeatColsDF$High))
        FeatColsDF <- FeatColsDF %>%
          relocate(P.Value, .after = `HR (95% CI)`)
        FeatColsDF$P.Value <- ifelse(FeatColsDF$P.Value < 0.1,formatC(FeatColsDF$P.Value, format = "e", digits = 1),FeatColsDF$P.Value)
        
        if (any(grepl(paste(StatCols,collapse = "|"),FeatColsDF$Variable))) {
          found <- grep(paste(StatCols,collapse = "|"),FeatColsDF$Variable)
          for (i in found) {
            coln <- strsplit(FeatColsDF[i,1],"\\s+")[[1]][1]
            ref <- strsplit(FeatColsDF[i,1],"\\s+")[[1]][2]
            FeatColsDF[i,1] <- paste0(geneset_name,"_",coln," ",ref)
          }
        }
        FeatColsDF$Variable <- paste(FeatColsDF$Variable,"     ")
        FeatColsDF
      })
      
      MultiFeatUnivarForestPlot_react <- reactive({
        FeatColsDF <- MultiFeatUnivarForestPlotTab_react()
        Xlims <- input$UnivarForestPlotXlim
        Xtrans <- input$UnivarForestPlotXtrans
        Xlims <- gsub("\\s+","",Xlims)
        Xlims_low <- as.numeric(strsplit(Xlims,",")[[1]][1])
        Xlims_high <- as.numeric(strsplit(Xlims,",")[[1]][2])
        Xtrans <- input$UnivarForestPlotXtrans
        if (Xtrans != "none") {
          XLabel <- paste0("HR (95% CI ",Xtrans,")")
        } else {
          XLabel <- "HR (95% CI)"
        }
        tm <- forest_theme(base_size = 16,
                           # Confidence interval point shape, line type/color/width
                           ci_pch = 15,ci_col = "black",ci_fill = "black",ci_alpha = 0.8,
                           ci_lty = 1,ci_lwd = 1.5,ci_Theight = 0.2, # Set an T end at the end of CI
                           # Reference line width/type/color
                           refline_lwd = 1,refline_lty = "dashed",refline_col = "grey20",
                           # Vertical line width/type/color
                           vertline_lwd = 1,vertline_lty = "dashed",vertline_col = "grey20",
                           # Change summary color for filling and borders
                           summary_fill = "#4575b4",summary_col = "#4575b4")
        p_OS <- forest(FeatColsDF[,c(1:2,7:9)],
                       #title = "PROMOTE OS",
                       est = FeatColsDF$`Hazard Ratio`,
                       lower = FeatColsDF$Low,
                       upper = FeatColsDF$High,
                       #sizes = coxOS$`Standard Error`,
                       x_trans = Xtrans,
                       xlab = XLabel,
                       ci_column = 3,
                       ref_line = 1,
                       xaxis_cex = 0.5,
                       arrow_lab = c("Low Risk", "High Risk"),
                       xlim = c(Xlims_low, Xlims_high),
                       #ticks_at = c(0.5, 1, 2, 3, 4, 5),
                       theme = tm
                       
        )
        p_OS$heights[length(p_OS$heights)-1] <- unit((as.numeric(p_OS$heights[length(p_OS$heights)-1])+8),"mm")
        p_OS2 <- edit_plot(p_OS, row = 1, gp = gpar(col = "red", fontface = "italic"))
        p_wh<-get_wh(p_OS2)
        plotht(round(unname(p_wh[2])*100))
        p_OS2
        
      })
      output$rendMultiFeatUnivarForestPlot <- renderUI({
        shinycssloaders::withSpinner(plotOutput("MultiFeatUnivarForestPlot", width = "100%", height = plotht()), type = 6)
      })
      
      output$MultiFeatUnivarForestPlot <- renderPlot({
        forest <- MultiFeatUnivarForestPlot_react()
        forest
      })
      
      UnivarLinearityPlot_react <- reactive({
        req(UniVarFeatTab_react())
        req(input$SingleSurvivalFeature)
        req(input$ResidualTypeUni)
        req(input$linPredict1)
        req(input$linAxisFont)
        req(input$linMainFont)
        req(input$linTickFont)
        linearityPlot(UniVarFeatTab_react(),input$SingleSurvivalFeature,input$ResidualTypeUni,input$linPredict1,
                      input$linAxisFont,input$linMainFont,input$linTickFont)
      })
      output$UnivarLinearityPlot <- renderPlot({
        UnivarLinearityPlot_react()
      })
      
      ## Bivar Add -------------------------------------------------------------
      
      observe({
        req(ssGSEAmeta())
        meta <- ssGSEAmeta()
        geneset <- gs_react()
        geneset_name <- names(geneset)
        FeatureChoices <- meta %>%
          dplyr::select(where(~ n_distinct(.x[nzchar(.x)], na.rm = TRUE) > 1)) %>%
          names
        #metacol_feature()
        updateSelectizeInput(session = session,inputId = "SurvivalFeatureBi1", choices = FeatureChoices, selected = paste0(geneset_name,"_MedianCutP"), server = T)
      })
      observe({
        req(ssGSEAmeta())
        meta <- ssGSEAmeta()
        FeatureChoices <- meta %>%
          dplyr::select(where(~ n_distinct(.x[nzchar(.x)], na.rm = TRUE) > 1)) %>%
          names
        if (isTruthy(PreSelect_SecondaryFeature_react())) {
          if (!PreSelect_SecondaryFeature_react() %in% FeatureChoices) {
            selec <- FeatureChoices[1]
          } else { selec <- PreSelect_SecondaryFeature_react() }
        } else { selec <- FeatureChoices[1] }
        updateSelectizeInput(session = session,inputId = "SurvivalFeatureBi2", choices = FeatureChoices,selected = selec, server = T)
      })
      output$rendSurvFeatVariableBi1 <- renderUI({
        req(ssGSEAmeta())
        req(input$SurvivalFeatureBi1)
        Feature <- input$SurvivalFeatureBi1
        meta <- ssGSEAmeta()
        Var_choices <- survFeatRefSelect(meta,Feature,input$BiVarAddNAcheck1,input$BiVarAddContCheck1,input$BiVarAddContHiLoCheck1)
        selectInput("SurvFeatVariableBi1","Select Coxh Feature Reference:", choices = Var_choices, width = "80%")
      })
      output$rendSurvFeatVariableBi2 <- renderUI({
        req(ssGSEAmeta())
        req(input$SurvivalFeatureBi2)
        Feature <- input$SurvivalFeatureBi2
        meta <- ssGSEAmeta()
        Var_choices <- survFeatRefSelect(meta,Feature,input$BiVarAddNAcheck2,input$BiVarAddContCheck2,input$BiVarAddContHiLoCheck2)
        selectInput("SurvFeatVariableBi2","Select Coxh Feature Reference:", choices = Var_choices, width = "80%")
      })
      
      BiVarAddFeature_react <- reactive({
        req(input$SurvivalFeatureBi1)
        req(input$SurvFeatVariableBi1)
        req(input$SurvivalFeatureBi2)
        req(input$SurvFeatVariableBi2)
        req(ssGSEAmeta())
        meta_ssgsea <- ssGSEAmeta()
        Feature1 <- input$SurvivalFeatureBi1
        Feature1ref <- input$SurvFeatVariableBi1
        Feature2 <- input$SurvivalFeatureBi2
        Feature2ref <- input$SurvFeatVariableBi2
        
        if ((Feature1ref %in% levels(factor(meta_ssgsea[,Feature1]))) & (Feature2ref %in% levels(factor(meta_ssgsea[,Feature2])))) {
          if (input$BiVarAddNAcheck1 == TRUE) {
            meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature1]) == FALSE),]
            meta_ssgsea <- meta_ssgsea[which(meta_ssgsea[,Feature1] != "Inf" & meta_ssgsea[,Feature1] != "N/A" & meta_ssgsea[,Feature1] != "n/a"),]
            meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature1],ignore.case = T, invert = T),]
          }
          if (input$BiVarAddContCheck1 == TRUE) {
            if (input$BiVarAddContHiLoCheck1 == TRUE) {
              meta_ssgsea[,Feature1] <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == Feature1)])
            }
          }
          if (input$BiVarAddContCheck1 == FALSE) {
            meta_ssgsea[,Feature1] <- factor(meta_ssgsea[,Feature1])
            meta_ssgsea[,Feature1] <- relevel(meta_ssgsea[,Feature1], ref = Feature1ref)
          }
          
          if (input$BiVarAddNAcheck2 == TRUE) {
            meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature2]) == FALSE),]
            meta_ssgsea <- meta_ssgsea[which(meta_ssgsea[,Feature2] != "Inf" & meta_ssgsea[,Feature2] != "N/A" & meta_ssgsea[,Feature2] != "n/a"),]
            meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature2],ignore.case = T, invert = T),]
          }
          if (input$BiVarAddContCheck2 == TRUE) {
            if (input$BiVarAddContHiLoCheck2 == TRUE) {
              meta_ssgsea[,Feature2] <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == Feature2)])
            }
          }
          if (input$BiVarAddContCheck2 == FALSE) {
            meta_ssgsea[,Feature2] <- factor(meta_ssgsea[,Feature2])
            meta_ssgsea[,Feature2] <- relevel(meta_ssgsea[,Feature2], ref = Feature2ref)
          }
          
          SubsetSurvData(meta_ssgsea,input$SurvivalType_time,input$SurvivalType_id,Feature1,Feature2)
        }
        
      })
      BiVarAddTab_react <- reactive({
        req(input$SurvivalFeatureBi1)
        req(input$SurvivalFeatureBi2)
        req(BiVarAddFeature_react())
        df <- BiVarAddFeature_react()
        Feature1 <- input$SurvivalFeatureBi1
        Feature2 <- input$SurvivalFeatureBi2
        colnames(df)[which(colnames(df) == Feature1)] <- gsub("-","_",Feature1)
        Feature1 <- gsub("-","_",Feature1)
        colnames(df)[which(colnames(df) == Feature2)] <- gsub("-","_",Feature2)
        Feature2 <- gsub("-","_",Feature2)
        Feature1 <- sprintf(ifelse((grepl(" ", Feature1) | !is.na(suppressWarnings(as.numeric(substring(Feature1, 1, 1))))), "`%s`", "%s"), Feature1)
        Feature2 <- sprintf(ifelse((grepl(" ", Feature2) | !is.na(suppressWarnings(as.numeric(substring(Feature2, 1, 1))))), "`%s`", "%s"), Feature2)
        tab <- coxph(as.formula(paste0("Surv(time,ID) ~ ",paste0(Feature1,"+",Feature2))),data = df)
        tab
      })
      BiVarAddTabFeat1_react <- reactive({
        Feature1 <- input$SurvivalFeatureBi1
        df <- BiVarAddFeature_react()
        colnames(df)[which(colnames(df) == Feature1)] <- gsub("-","_",Feature1)
        Feature1 <- gsub("-","_",Feature1)
        Feature1 <- sprintf(ifelse((grepl(" ", Feature1) | !is.na(suppressWarnings(as.numeric(substring(Feature1, 1, 1))))), "`%s`", "%s"), Feature1)
        tab <- coxph(as.formula(paste0("Surv(time,ID) ~ ",Feature1)),data = df)
        tab
      })
      BiVarAddTabFeat2_react <- reactive({
        Feature2 <- input$SurvivalFeatureBi2
        df <- BiVarAddFeature_react()
        colnames(df)[which(colnames(df) == Feature2)] <- gsub("-","_",Feature2)
        Feature2 <- gsub("-","_",Feature2)
        Feature2 <- sprintf(ifelse((grepl(" ", Feature2) | !is.na(suppressWarnings(as.numeric(substring(Feature2, 1, 1))))), "`%s`", "%s"), Feature2)
        tab <- coxph(as.formula(paste0("Surv(time,ID) ~ ",Feature2)),data = df)
        tab
      })
      BiVarAddHRTab_react <- reactive({
        if (isTruthy(BiVarAddTab_react())) {
          tab_df <- CoxPHtabUni(BiVarAddTab_react())
          tab_df
        }
      })
      output$BiFeatureHRtab <- renderTable({
        tab <- BiVarAddHRTab_react()
        tab
      })
      output$bivarSummary <- renderPrint({
        req(BiVarAddTab_react())
        CoxPHsumm(BiVarAddTab_react(),bivarAdd = TRUE)
      })
      
      output$bivarAnova1 <- renderPrint({
        biVarAnova(BiVarAddTab_react(),BiVarAddTabFeat1_react())
      })
      
      output$bivarAnova2 <- renderPrint({
        biVarAnova(BiVarAddTab_react(),BiVarAddTabFeat2_react())
      })
      
      BivarForestPlot_react <- reactive({
        req(BiVarAddTab_react())
        req(input$SurvivalFeatureBi1)
        req(input$SurvivalFeatureBi2)
        req(BiVarAddFeature_react())
        req(input$ForestFontSize)
        df <- BiVarAddFeature_react()
        obj <- BiVarAddTab_react()
        fontSize <- input$ForestFontSize
        Feature1 <- input$SurvivalFeatureBi1
        Feature2 <- input$SurvivalFeatureBi2
        #save(list = ls(), file = "Forest_env.RData", envir = environment())
        forest <- forestPlot_Simple(obj,df,paste0(Feature1,"+",Feature2),fontSize)
      })
      output$BivarForestPlot <- renderPlot({
        forest <- BivarForestPlot_react()
        forest
      })
      
      BivarLinearityPlot_react <- reactive({
        
        if (length(input$SurvivalFeatureBi1 > 0) & length(input$SurvivalFeatureBi2 > 0)) {
          P <- linearityPlot(BiVarAddTab_react(),paste0(input$SurvivalFeatureBi1,"+",input$SurvivalFeatureBi2),
                             input$ResidualTypeBi,input$linPredict2,input$linAxisFont,input$linMainFont,input$linTickFont)
        }
      })
      
      output$BivarLinearityPlot <- renderPlot({
        p <- BivarLinearityPlot_react()
        p
      })
      
      ## Bivar Int -------------------------------------------------------------
      
      observe({
        req(ssGSEAmeta())
        meta <- ssGSEAmeta()
        geneset <- gs_react()
        geneset_name <- names(geneset)
        FeatureChoices <- meta %>%
          dplyr::select(where(~ n_distinct(.x[nzchar(.x)], na.rm = TRUE) > 1)) %>%
          names
        updateSelectizeInput(session = session,inputId = "SurvivalFeatureBi1Inter", choices = FeatureChoices, selected = paste0(geneset_name,"_MedianCutP"), server = T)
      })
      observe({
        req(ssGSEAmeta())
        meta <- ssGSEAmeta()
        FeatureChoices <- meta %>%
          dplyr::select(where(~ n_distinct(.x[nzchar(.x)], na.rm = TRUE) > 1)) %>%
          names
        if (isTruthy(PreSelect_SecondaryFeature_react())) {
          if (!PreSelect_SecondaryFeature_react() %in% FeatureChoices) {
            selec <- FeatureChoices[1]
          } else { selec <- PreSelect_SecondaryFeature_react() }
        } else { selec <- FeatureChoices[1] }
        updateSelectizeInput(session = session,inputId = "SurvivalFeatureBi2Inter", choices = FeatureChoices,selected = selec, server = T)
      })
      output$rendSurvFeatVariableBi1Inter <- renderUI({
        Feature <- input$SurvivalFeatureBi1Inter
        req(ssGSEAmeta())
        meta <- ssGSEAmeta()
        Var_choices <- survFeatRefSelect(meta,Feature,input$BiVarIntNAcheck1,input$BiVarIntContCheck1,input$BiVarIntContHiLoCheck1)
        selectInput("SurvFeatVariableBi1Inter","Select Coxh Feature Reference:", choices = Var_choices, width = "80%")
      })
      output$rendSurvFeatVariableBi2Inter <- renderUI({
        Feature <- input$SurvivalFeatureBi2Inter
        req(ssGSEAmeta())
        meta <- ssGSEAmeta()
        Var_choices <- survFeatRefSelect(meta,Feature,input$BiVarIntNAcheck2,input$BiVarIntContCheck2,input$BiVarIntContHiLoCheck2)
        selectInput("SurvFeatVariableBi2Inter","Select Coxh Feature Reference:", choices = Var_choices, width = "80%")
      })
      
      BiVarIntFeature_react <- reactive({
        req(input$SurvivalFeatureBi1Inter)
        req(input$SurvivalFeatureBi2Inter)
        req(input$SurvFeatVariableBi1Inter)
        req(input$SurvFeatVariableBi2Inter)
        req(ssGSEAmeta())
        meta_ssgsea <- ssGSEAmeta()
        Feature1 <- input$SurvivalFeatureBi1Inter
        Feature1ref <- input$SurvFeatVariableBi1Inter
        if (Feature1 %in% colnames(meta_ssgsea)) {
          if (input$BiVarIntNAcheck1 == TRUE) {
            meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature1]) == FALSE),]
            meta_ssgsea <- meta_ssgsea[which(meta_ssgsea[,Feature1] != "Inf" & meta_ssgsea[,Feature1] != "N/A" & meta_ssgsea[,Feature1] != "n/a"),]
            meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature1],ignore.case = T, invert = T),]
          }
          if (input$BiVarIntContCheck1 == TRUE) {
            if (input$BiVarIntContHiLoCheck1 == TRUE) {
              meta_ssgsea[,Feature1] <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == Feature1)])
            }
          }
          if (input$BiVarIntContCheck1 == FALSE) {
            meta_ssgsea[,Feature1] <- factor(meta_ssgsea[,Feature1])
            meta_ssgsea[,Feature1] <- relevel(meta_ssgsea[,Feature1], ref = Feature1ref)
          }
          
          Feature2 <- input$SurvivalFeatureBi2Inter
          Feature2ref <- input$SurvFeatVariableBi2Inter
          if (input$BiVarIntNAcheck2 == TRUE) {
            meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature2]) == FALSE),]
            meta_ssgsea <- meta_ssgsea[which(meta_ssgsea[,Feature2] != "Inf" & meta_ssgsea[,Feature2] != "N/A" & meta_ssgsea[,Feature2] != "n/a"),]
            meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature2],ignore.case = T, invert = T),]
          }
          if (input$BiVarIntContCheck2 == TRUE) {
            if (input$BiVarIntContHiLoCheck2 == TRUE) {
              meta_ssgsea[,Feature2] <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == Feature2)])
            }
          }
          if (input$BiVarIntContCheck2 == FALSE) {
            meta_ssgsea[,Feature2] <- factor(meta_ssgsea[,Feature2])
            meta_ssgsea[,Feature2] <- relevel(meta_ssgsea[,Feature2], ref = Feature2ref)
          }
          SubsetSurvData(meta_ssgsea,input$SurvivalType_time,input$SurvivalType_id,Feature1,Feature2)
        }
        
      })
      
      BiVarIntTab_react <- reactive({
        req(input$SurvivalFeatureBi1Inter)
        req(input$SurvivalFeatureBi2Inter)
        df <- BiVarIntFeature_react()
        Feature1 <- input$SurvivalFeatureBi1Inter
        Feature2 <- input$SurvivalFeatureBi2Inter
        colnames(df)[which(colnames(df) == Feature1)] <- gsub("-","_",Feature1)
        Feature1 <- gsub("-","_",Feature1)
        colnames(df)[which(colnames(df) == Feature2)] <- gsub("-","_",Feature2)
        Feature2 <- gsub("-","_",Feature2)
        Feature1 <- sprintf(ifelse((grepl(" ", Feature1) | !is.na(suppressWarnings(as.numeric(substring(Feature1, 1, 1))))), "`%s`", "%s"), Feature1)
        Feature2 <- sprintf(ifelse((grepl(" ", Feature2) | !is.na(suppressWarnings(as.numeric(substring(Feature2, 1, 1))))), "`%s`", "%s"), Feature2)
        form <- as.formula(paste0("Surv(time,ID) ~ ",paste0(Feature1,"*",Feature2)))
        tab <- eval(substitute(coxph(form,data = df)))
        tab
      })
      BiVarIntTab4Annova_react <- reactive({
        req(input$SurvivalFeatureBi1Inter)
        req(input$SurvivalFeatureBi2Inter)
        df <- BiVarIntFeature_react()
        Feature1 <- input$SurvivalFeatureBi1Inter
        Feature2 <- input$SurvivalFeatureBi2Inter
        colnames(df)[which(colnames(df) == Feature1)] <- gsub("-","_",Feature1)
        Feature1 <- gsub("-","_",Feature1)
        colnames(df)[which(colnames(df) == Feature2)] <- gsub("-","_",Feature2)
        Feature2 <- gsub("-","_",Feature2)
        Feature1 <- sprintf(ifelse((grepl(" ", Feature1) | !is.na(suppressWarnings(as.numeric(substring(Feature1, 1, 1))))), "`%s`", "%s"), Feature1)
        Feature2 <- sprintf(ifelse((grepl(" ", Feature2) | !is.na(suppressWarnings(as.numeric(substring(Feature2, 1, 1))))), "`%s`", "%s"), Feature2)
        form <- as.formula(paste0("Surv(time,ID) ~ ",paste0(Feature1,"+",Feature2)))
        tab <- eval(substitute(coxph(form,data = df)))
        tab
      })
      
      BiVarIntHRTab_react <- reactive({
        if (isTruthy(BiVarIntTab_react())) {
          tab_df <- CoxPHtabUni(BiVarIntTab_react())
          tab_df
        }
      })
      output$BiFeatureHRtabInter <- renderTable({
        tab <- BiVarIntHRTab_react()
        tab
      })
      output$bivarSummaryInter <- renderPrint({
        req(BiVarIntTab_react())
        CoxPHsumm(BiVarIntTab_react(),bivarInt = TRUE)
      })
      
      output$bivarAnovaInter1 <- renderPrint({
        biVarAnova(BiVarIntTab_react(),BiVarIntTab4Annova_react())
      })
      
      featSplotBi_react <- reactive({
        
        ## Assign variables
        req(BiVarIntFeature_react())
        meta_ssgsea_sdf <- BiVarIntFeature_react()
        SampleType <- input$SampleTypeSelection
        Feature1 <- input$SurvivalFeatureBi1Inter
        Feature2 <- input$SurvivalFeatureBi2Inter
        Feature_sub <- input$FeatureSelection
        subFeature <- input$subFeatureSelection
        show_pval <- input$ShowPval
        ShowConfInt <- input$ShowConfInt
        if (input$SurvYearOrMonth == "Years") {
          xBreaks <- input$SurvXaxisBreaks * 365.25
        } else {
          xBreaks <- input$SurvXaxisBreaks * 30.4375
        }
        if (isTruthy(input$SurvXaxis)) {
          if (input$SurvYearOrMonth == "Years") {
            xaxlim <- input$SurvXaxis * 365.25
            #meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,2] <= xaxlim),]
            #xBreaks <- input$SurvXaxisBreaks * 365.25
          } else {
            #xBreaks <- input$SurvXaxisBreaks * 30.4375
            #meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,2] <= xaxlim),]
            xaxlim <- input$SurvXaxis * 30.4375
          }
        } else {
          xaxlim <- NULL
          #xBreaks <- 365.25
        }
        surv_time_col <- input$SurvivalType_time
        showLegend <- input$SurvLegendPos
        showMedSurv <- input$ShowMedSurvLine
        if (showMedSurv == T) {
          showMedSurv <- "hv"
        }
        else if (showMedSurv == F) {
          showMedSurv <- "none"
        }
        SurvDateType <- sub("\\..*","",surv_time_col)
        
        PlotTitle <- SurvPlotTitle(SampleTypeSelected = SampleType,Feature = Feature_sub, subFeature = subFeature,
                                   multivar = paste0(Feature1," + ",Feature2))
        
        
        #Feature1 <- sprintf(ifelse(grepl(" ", Feature1), "`%s`", "%s"), Feature1)
        #Feature2 <- sprintf(ifelse(grepl(" ", Feature2), "`%s`", "%s"), Feature2)
        
        colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == Feature1)] <- "Feature1"
        colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == Feature2)] <- "Feature2"
        form <- paste0("Surv(time,ID) ~ Feature1 + Feature2")
        
        #save(list = ls(), file = "shiny_env_mutliInt.Rdata", envir = environment())
        fit <- eval(substitute(survfit(as.formula(form),data = meta_ssgsea_sdf, type="kaplan-meier")))
        names(fit[["strata"]]) <- gsub("^Feature1=",paste0(Feature1,"="),names(fit[["strata"]]))
        names(fit[["strata"]]) <- gsub(", Feature2=",paste0(Feature2,"="),names(fit[["strata"]]))
        
        #form <- as.formula(paste0("Surv(time,ID) ~ ",paste0(Feature1,"+",Feature2)))
        #fit <- eval(substitute(survfit(form,data = meta_ssgsea_sdf, type="kaplan-meier")))
        
        SurvPlot(fit,meta_ssgsea_sdf,PlotTitle,ylab = paste(SurvDateType,"Survival Probability"),
                 pval = show_pval,conf = ShowConfInt,legend = showLegend,median = showMedSurv,xlim = xaxlim,
                 xScale = input$SurvYearOrMonth, xBreaks = xBreaks)
      })
      
      output$featSplotBi <- renderPlot({
        plot <- featSplotBi_react()
        plot
      })
      
      BivarLinearityPlotInter_react <- reactive({
        
        if (length(input$SurvivalFeatureBi1Inter > 0) & length(input$SurvivalFeatureBi2Inter > 0)) {
          P <- linearityPlot(BiVarIntTab_react(),paste0(input$SurvivalFeatureBi1Inter,"*",input$SurvivalFeatureBi2Inter),
                             input$ResidualTypeInter,input$linPredict3,input$linAxisFont,input$linMainFont,input$linTickFont)
        }
      })
      
      output$BivarLinearityPlotInter <- renderPlot({
        p <- BivarLinearityPlotInter_react()
        p
      })
      
      ## Multivar --------------------------------------------------------------
      
      observe({
        req(ssGSEAmeta())
        meta <- ssGSEAmeta()
        geneset <- gs_react()
        geneset_name <- names(geneset)
        FeatureChoices <- meta %>%
          dplyr::select(where(~ n_distinct(.x[nzchar(.x)], na.rm = TRUE) > 1)) %>%
          names
        updateSelectizeInput(session = session,inputId = "SurvivalFeature", choices = FeatureChoices, selected = paste0(geneset_name,"_MedianCutP"), server = T)
      })
      
      MultiVarFeat_react <- reactive({
        if (length(input$SurvivalFeature > 0)) {
          Feature <- input$SurvivalFeature
          surv_time_col <- input$SurvivalType_time
          surv_id_col <- input$SurvivalType_id
          meta_ssgsea <- ssGSEAmeta()
          SampleNameCol <- colnames(meta_ssgsea)[1]
          if (input$UniVarNAcheck == TRUE) {
            for (i in Feature) {
              meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,i]) == FALSE),]
              meta_ssgsea <- meta_ssgsea[which(meta_ssgsea[,i] != "Inf"),]
              meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,i],ignore.case = T, invert = T),]
            }
          }
          select_cols <- c(SampleNameCol,surv_time_col,surv_id_col,Feature)
          meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
          colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "time"
          colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "ID"
          meta_ssgsea_sdf
        }
      })
      
      MultiVarFeatCat_react <- reactive({
        meta_ssgsea_sdf <- MultiVarFeat_react()
        Feature <- input$SurvivalFeature
        for (i in Feature){
          meta_ssgsea_sdf[,i] <- as.factor(meta_ssgsea_sdf[,i])
        }
        meta_ssgsea_sdf
      })
      
      MultiVarFeatCont_react <- reactive({
        meta_ssgsea_sdf <- MultiVarFeat_react()
        meta_ssgsea_sdf
      })
      
      MultiVarTabCat_react <- reactive({
        req(MultiVarFeatCat_react())
        req(input$SurvivalFeature)
        meta_ssgsea_sdf <- MultiVarFeatCat_react()
        Feature <- input$SurvivalFeature
        colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == Feature)] <- gsub("-","_",Feature)
        Feature <- gsub("-","_",Feature)
        Feature <- sprintf(ifelse((grepl(" ", Feature) | !is.na(suppressWarnings(as.numeric(substring(Feature, 1, 1))))), "`%s`", "%s"), Feature)
        form <- as.formula(paste0("Surv(time,ID) ~ ",paste(Feature,collapse = "+")))
        tab <- eval(substitute(coxph(form,data = meta_ssgsea_sdf)))
        tab
      })
      
      MultiVarTabCont_react <- reactive({
        meta_ssgsea_sdf <- MultiVarFeatCont_react()
        Feature <- input$SurvivalFeature
        colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == Feature)] <- gsub("-","_",Feature)
        Feature <- gsub("-","_",Feature)
        Feature <- sprintf(ifelse((grepl(" ", Feature) | !is.na(suppressWarnings(as.numeric(substring(Feature, 1, 1))))), "`%s`", "%s"), Feature)
        form <- as.formula(paste0("Surv(time,ID) ~ ",paste(Feature,collapse = "+")))
        tab <- eval(substitute(coxph(form,data = meta_ssgsea_sdf)))
        tab
      })
      
      SFeatureHRtabCat_react <- reactive({
        if (length(input$SurvivalFeature > 0)) {
          tab_df <- CoxPHtabUni(MultiVarTabCat_react())
          tab_df
        }
      })
      output$SFeatureHRtabCat <- renderTable({
        tab <- SFeatureHRtabCat_react()
        tab
      })
      
      SFeatureHRtabCont_react <- reactive({
        if (length(input$SurvivalFeature > 0)) {
          tab_df <- CoxPHtabUni(MultiVarTabCont_react())
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
        xph <- capture.output(cox.zph(tab))
        con_line <- grep("^Concordance=",out,value = T)
        lik_line <- grep("^Likelihood ratio test=",out,value = T)
        wal_line <- grep("^Wald test",out,value = T)
        sco_line <- grep("^Score ",out,value = T)
        text <- paste("CoxH Summary:",con_line,lik_line,wal_line,sco_line,"","Proportional Hazards assumption:",paste(xph[1:length(xph)],collapse = "\n"),sep = "\n")
        cat(text)
      })
      
      output$multivarSummaryCont <- renderPrint({
        tab <- MultiVarTabCont_react()
        out <- capture.output(summary(tab))
        xph <- capture.output(cox.zph(tab))
        con_line <- grep("^Concordance=",out,value = T)
        lik_line <- grep("^Likelihood ratio test=",out,value = T)
        wal_line <- grep("^Wald test",out,value = T)
        sco_line <- grep("^Score ",out,value = T)
        text <- paste("CoxH Summary:",con_line,lik_line,wal_line,sco_line,"","Proportional Hazards assumption:",paste(xph[1:length(xph)],collapse = "\n"),sep = "\n")
        cat(text)
      })
      
      MultivarForestPlot_react <- reactive({
        req(MultiVarTabCat_react())
        req(MultiVarFeat_react())
        req(input$SurvivalFeature)
        req(input$ForestFontSize)
        tab <- MultiVarTabCat_react()
        meta_ssgsea_sdf <- MultiVarFeat_react()
        Feature <- input$SurvivalFeature
        forest <- forestPlot_Simple(tab,meta_ssgsea_sdf,paste0(Feature,collapse = ", "),input$ForestFontSize)
      })
      output$MultivarForestPlot <- renderPlot({
        forest <- MultivarForestPlot_react()
        forest
      })
      
      ## Multivar Forest -------------------------------------------------------
      
      observe({
        req(ssGSEAmeta())
        meta <- ssGSEAmeta()
        geneset <- gs_react()
        geneset_name <- names(geneset)
        ColLevels <- apply(meta,2,function(x) length(levels(as.factor(x))))
        DicotCols <- names(ColLevels[ColLevels==2])
        ContCols <- GetColsOfType(meta,"continuous")
        FeatureChoices <- unique(c(DicotCols,ContCols))
        SurvDataCols <- c(metacol_survid,metacol_survtime)
        FeatureChoices <- FeatureChoices[which(!FeatureChoices %in% SurvDataCols)]
        #save(list = ls(), file = "Shiny_env.RData", envir = environment())
        updateSelectizeInput(session,"MultiFeatMultivarSelect",choices = FeatureChoices, selected = paste0(geneset_name,"_MedianCutP"))
        #selectInput("MultiFeatMultivarSelect", "Select Main Forest Feature", choices = FeatureChoices, selected = "MedianCutP", multiple = F)
      })
      
      output$rendMultiFeatMultivarRefSelect <- renderUI({
        req(ssGSEAmeta())
        meta <- ssGSEAmeta()
        Feature <- input$MultiFeatMultivarSelect
        IsItCont <- GetColsOfType(meta[,Feature, drop = F],"continuous")
        if (isTruthy(IsItCont)) {
          if (IsItCont == Feature) {
            meta[,Feature] <- apply(meta[,Feature, drop = F],2,function(x) highlow2(x))
          }
        }
        opts <- unique(meta[,Feature])
        opts <- sort(opts, decreasing = T, na.last = T)
        selectInput("MultiFeatMultivarRefSelect", "Select Reference:", choices = opts)
      })
      
      observe({
        req(ssGSEAmeta())
        meta <- ssGSEAmeta()
        FeatureChoices <- colnames(meta)
        SurvDataCols <- c(metacol_survid,metacol_survtime)
        FeatureChoices <- FeatureChoices[which(!FeatureChoices %in% SurvDataCols)]
        ColLevels <- apply(meta,2,function(x) length(levels(as.factor(x))))
        NonNameCols <- names(ColLevels[ColLevels<nrow(meta)])
        FeatureChoices <- intersect(FeatureChoices,NonNameCols)
        if (isTruthy(PreSelect_SecondaryFeature_react())) {
          if (!PreSelect_SecondaryFeature_react() %in% FeatureChoices) {
            selec <- FeatureChoices[1]
          } else { selec <- PreSelect_SecondaryFeature_react() }
        } else { selec <- FeatureChoices[1] }
        updateSelectizeInput(session,"MultiFeatMultivarSubSelect",choices = FeatureChoices, selected = selec)
      })
      
      plothtMulti <- reactiveVal(300)
      output$rendMultiFeatMultivarForestPlot <- renderUI({
        shinycssloaders::withSpinner(plotOutput("MultiFeatMultivarForestPlot", width = "100%", height = plothtMulti()), type = 6)
      })
      
      MultiFeat_InterForestMeta <- reactive({
        
        req(ssGSEAmeta())
        req(input$SurvivalType_time)
        req(input$SurvivalType_id)
        req(input$MultiFeatMultivarSelect)
        req(input$MultiFeatMultivarRefSelect)
        meta <- ssGSEAmeta()
        surv_time_col <- input$SurvivalType_time
        surv_id_col <- input$SurvivalType_id
        mainFeat <- input$MultiFeatMultivarSelect
        mainFeatRef <- input$MultiFeatMultivarRefSelect
        subFeats <- input$MultiFeatMultivarSubSelect
        FeatCols <- c(colnames(meta)[1],surv_time_col,surv_id_col,mainFeat,subFeats)
        metaSub <- meta[,FeatCols]
        continuousCols <- GetColsOfType(metaSub[,-c(1:3), drop = F],"continuous")
        metaSub[,continuousCols] <- apply(metaSub[,continuousCols, drop = F],2,function(x) highlow2(x))
        if (mainFeatRef %in% metaSub[,mainFeat]) {
          metaSub[,mainFeat] <- as.factor(metaSub[,mainFeat])
          metaSub[,mainFeat] <- relevel(metaSub[,mainFeat], ref = mainFeatRef)
          metaSub
        }
      })
      
      output$multiForestPlotTable <- DT::renderDataTable({
        
        df <- MultiFeat_InterForestMeta()
        
        DT::datatable(df,
                      extensions = "FixedColumns",
                      options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                     pageLength = 20,
                                     fixedColumns = list(leftColumns = 1),
                                     scrollX = T),
                      rownames = F)
        
      })
      
      
      plothtMulti <- reactiveVal(300)
      
      MultiFeatMultivarForestPlotTab_react <- reactive({
        
        metaSub <- MultiFeat_InterForestMeta()
        geneset <- gs_react()                 #Chosen Gene Set
        geneset_name <- names(geneset)        #Name of chosen gene set
        surv_time_col <- input$SurvivalType_time
        surv_id_col <- input$SurvivalType_id
        mainFeat <- input$MultiFeatMultivarSelect
        subFeats <- input$MultiFeatMultivarSubSelect
        tabTemp <- c(N = NA,`Hazard Ratio` = NA,`Standard Error` = NA,Low = NA,High = NA,P.Value = NA)
        
        FeatColsList <- list()
        j <- 1
        for (feat in subFeats) {
          featVals <- unique(metaSub[,feat])
          FeatColsList[[feat]] <- c(Variable = feat,tabTemp)
          for (feat2 in featVals) {
            metaSubSet <- metaSub[which(metaSub[,feat] == feat2),c(colnames(metaSub)[1],surv_time_col,surv_id_col,mainFeat)]
            metaSubSet <- metaSubSet[which(is.na(metaSubSet[,mainFeat]) == FALSE),]
            metaSubSet <- metaSubSet[which(metaSubSet[,mainFeat] != "Inf"),]
            metaSubSet <- metaSubSet[grep("unknown",metaSubSet[,mainFeat],ignore.case = T, invert = T),]
            if (nrow(metaSubSet) >= 6){
              #if (length(levels(metaSubSet[,mainFeat])) >= 2){
              colnames(metaSubSet)[2:4] <- c("time","status","feature")
              tab <- coxph(Surv(time,status) ~ feature, dat = metaSubSet)
              #tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",mainFeat,sep = "")), data = metaSubSet)
              tabData <- get_tabData(tab)
              tabData[1] <- feat2
              FeatColsList[[as.character(j)]] <- tabData
              j <- j+1
            }
            
          }
          
        }
        
        FeatColsDF <- as.data.frame(do.call(rbind, FeatColsList))
        rownames(FeatColsDF) <- NULL
        
        FeatColsDF$Variable <- ifelse(is.na(FeatColsDF$N),
                                      FeatColsDF$Variable,
                                      paste0("   ", FeatColsDF$Variable))
        
        FeatColsDF$Low <- as.numeric(FeatColsDF$Low)
        FeatColsDF$High <- as.numeric(FeatColsDF$High)
        FeatColsDF$`Hazard Ratio` <- as.numeric(FeatColsDF$`Hazard Ratio`)
        FeatColsDF$P.Value <- as.numeric(FeatColsDF$P.Value)
        
        FeatColsDF$` ` <- paste(rep(" ", 20), collapse = " ")
        
        FeatColsDF$`HR (95% CI)` <- ifelse(is.na(FeatColsDF$`Standard Error`), "",
                                           sprintf("%.2f (%.2f to %.2f)",
                                                   FeatColsDF$`Hazard Ratio`, FeatColsDF$Low, FeatColsDF$High))
        FeatColsDF <- FeatColsDF %>%
          relocate(P.Value, .after = `HR (95% CI)`)
        FeatColsDF$P.Value <- ifelse(FeatColsDF$P.Value < 0.1,formatC(FeatColsDF$P.Value, format = "e", digits = 1),FeatColsDF$P.Value)
        
        FeatColsDF$P.Value <- ifelse(is.na(FeatColsDF$P.Value), "", FeatColsDF$P.Value)
        FeatColsDF$N <- ifelse(is.na(FeatColsDF$N), "", FeatColsDF$N)
        
        if (any(grepl(paste(StatCols,collapse = "|"),FeatColsDF$Variable))) {
          found <- grep(paste(StatCols,collapse = "|"),FeatColsDF$Variable)
          for (i in found) {
            coln <- strsplit(FeatColsDF[i,1],"\\s+")[[1]][1]
            ref <- strsplit(FeatColsDF[i,1],"\\s+")[[1]][2]
            FeatColsDF[i,1] <- paste0(geneset_name,"_",coln," ",ref)
          }
        }
        FeatColsDF$Variable <- paste(FeatColsDF$Variable,"     ")
        
        FeatColsDF
        
        
      })
      
      MultiFeatMultivarForestPlot_react <- reactive({
        
        FeatColsDF <- MultiFeatMultivarForestPlotTab_react()
        metaSub <- MultiFeat_InterForestMeta()
        Xtrans <- input$MultivarForestPlotXtrans
        
        Xlims <- input$MultivarForestPlotXlim
        Xlims <- gsub("\\s+","",Xlims)
        Xlims_low <- as.numeric(strsplit(Xlims,",")[[1]][1])
        Xlims_high <- as.numeric(strsplit(Xlims,",")[[1]][2])
        
        if (Xtrans != "none") {
          XLabel <- paste0("HR (95% CI ",Xtrans,")")
        } else {
          XLabel <- "HR (95% CI)"
        }
        
        surv_id_col <- input$SurvivalType_id
        mainFeat <- input$MultiFeatMultivarSelect
        mainFeatRef <- input$MultiFeatMultivarRefSelect
        metaSub[,mainFeat] <- as.factor(metaSub[,mainFeat])
        metaSub[,mainFeat] <- relevel(metaSub[,mainFeat], ref = mainFeatRef)
        ref <- levels(metaSub[,mainFeat])[2]
        PlotTitle <- paste0("Forest plot featuring ",mainFeat," - ",ref," amoung selected features\n")
        
        tm <- forest_theme(base_size = 16,
                           title_cex = 1.5,
                           # Confidence interval point shape, line type/color/width
                           ci_pch = 15,
                           ci_col = "black",
                           ci_fill = "black",
                           ci_alpha = 0.8,
                           ci_lty = 1,
                           ci_lwd = 1.5,
                           ci_Theight = 0.2, # Set an T end at the end of CI
                           # Reference line width/type/color
                           refline_lwd = 1,
                           refline_lty = "dashed",
                           refline_col = "grey20",
                           # Vertical line width/type/color
                           vertline_lwd = 1,
                           vertline_lty = "dashed",
                           vertline_col = "grey20",
                           # Change summary color for filling and borders
                           summary_fill = "#4575b4",
                           summary_col = "#4575b4")
        
        
        p_OS <- forest(FeatColsDF[,c(1:2,7:9)],
                       title = PlotTitle,
                       est = FeatColsDF$`Hazard Ratio`,
                       lower = FeatColsDF$Low,
                       upper = FeatColsDF$High,
                       #sizes = coxOS$`Standard Error`,
                       x_trans = Xtrans,
                       xlab = XLabel,
                       ci_column = 3,
                       ref_line = 1,
                       arrow_lab = c("Low Risk", "High Risk"),
                       xlim = c(Xlims_low, Xlims_high),
                       #ticks_at = c(0.5, 1, 2, 3, 4, 5),
                       theme = tm
                       
        )
        
        p_OS$heights[length(p_OS$heights)-1] <- unit((as.numeric(p_OS$heights[length(p_OS$heights)-1])+8),"mm")
        p_wh<-get_wh(p_OS)
        plothtMulti(round(unname(p_wh[2])*100)-100)
        p_OS
        
      })
      
      output$MultiFeatMultivarForestPlot <- renderPlot({
        MultiFeatMultivarForestPlot_react()
      })
      
      ## Lasso
      #
      #SampleType_Selec <- reactiveValues(SampleType = PreSelect_SamplyType_react())
      #observeEvent(input$SampleTypeSelection, {
      #  SampleType_Selec$SampleType <- input$SampleTypeSelection
      #})
      ### Select sample type to subset samples by - only render if more than one sample type
      #output$rendSampleTypeSelection_lasso <- renderUI({
      #  metacol_sampletype <- metacol_sampletype()
      #  meta <- meta_react()
      #  if (length(unique(meta[,metacol_sampletype])) > 1) {
      #    SampleTypeChoices <- unique(meta[,metacol_sampletype])
      #    SampleTypeChoices <- c("All Sample Types",SampleTypeChoices)
      #    selectInput("SampleTypeSelection_lasso",paste("Select Sample Type (",metacol_sampletype,"):",sep = ""),
      #                choices = SampleTypeChoices, selected = SampleType_Selec$SampleType)
      #  }
      #})
      #
      #Feature_Selec <- reactiveValues(Feature = PreSelect_Feature_react())
      #observeEvent(input$FeatureSelection, {
      #  Feature_Selec$Feature <- input$FeatureSelection
      #})
      #observe({
      #  req(metaP_react())
      #  MetaParam <- metaP_react()
      #  metacol_feature <- c("Show all Samples",MetaParam[,1])
      #  updateSelectizeInput(session = session, inputId = "FeatureSelection_lasso",
      #                       choices = metacol_feature, selected = Feature_Selec$Feature, server = T)
      #})
      #
      #SubFeature_Selec <- reactiveValues(SubFeature = PreSelect_SubFeature_react())
      #observeEvent(input$SubFeatureSelection, {
      #  SubFeature_Selec$SubFeature <- input$SubFeatureSelection
      #})
      #output$rendSubFeatureSelection_lasso <- renderUI({
      #  req(input$FeatureSelection_lasso)
      #  meta <- meta_react()
      #  if (isTruthy(input$SampleTypeSelection_lasso)) {
      #    if (input$SampleTypeSelection_lasso != "All Sample Types") {
      #      metacol_sampletype <- metacol_sampletype()
      #      meta <- meta[which(meta[,metacol_sampletype] == input$SampleTypeSelection_lasso),]
      #    } else {
      #      meta <- meta
      #    }
      #  }
      #  if (input$FeatureSelection_lasso != "Show all Samples") {
      #    SubFeatureChoices <- unique(meta[,input$FeatureSelection_lasso])
      #    SubFeatureChoices <- sort(SubFeatureChoices, decreasing = T, na.last = T)
      #    selectInput("SubFeatureSelection_lasso","Feature Condition:",choices = SubFeatureChoices, selected = PreSelect_SubFeature_react())
      #  }
      #})
      #
      #observe({
      #  req(metacol_survtime())
      #  SurTimeChoices <- metacol_survtime()
      #  updateSelectizeInput(session = session, inputId = "SurvTimeSelec_lasso",choices = SurTimeChoices, selected = input$SurvivalType_time, server = T)
      #})
      #
      #observe({
      #  req(metacol_survid())
      #  SurIDChoices <- metacol_survid()
      #  updateSelectizeInput(session = session, inputId = "SurvIDSelect_lasso",choices = SurIDChoices, selected = input$SurvivalType_id, server = T)
      #})
      #
      #observe({
      #  exprGenes <- rownames(expr_react())
      #  LassoFeatures <- c(exprGenes,decon_score_cols())
      #  updateSelectizeInput(session, "LassoFeatureSelection_lasso", choices = LassoFeatures,
      #                       selected = "", options = list(delimiter = " ", create = T), server = T)
      #})
      #
      #output$rendCutPinput <- renderUI({
      #  if (input$LassoPlotCutP == "Quantile") {
      #    numericInput("CutPinput","Top/Bottom Cut-Point Quantile Cutoff (%)", value = 25, min = 0, max = 100, width = "200px")
      #  }
      #  else if (input$LassoPlotCutP == "User Specified") {
      #    numericInput("CutPinput","Above/Below User Quantile Cut-Point (%)", value = 25, min = 0, max = 100, width = "200px")
      #  }
      #})
      #
      #output$rendCustomLambda <- renderUI({
      #
      #  if (input$viewLassoMinOrSE == "Custom") {
      #    numericInput("CustomLambda","Custom Lambda:",min = 0, value = "")
      #  }
      #  else if(input$viewLassoMinOrSE == "Lambda Min") {
      #    model <- LassoRun_train_model()
      #    l_min <- model$lambda.min
      #    p(paste("Lambda Min:",l_min))
      #  }
      #  else if(input$viewLassoMinOrSE == "Lambda SE") {
      #    model <- LassoRun_train_model()
      #    l_se <- model$lambda.1se
      #    p(paste("Lambda SE:",l_se))
      #  }
      #})
      #
      #output$rednLassoCoefTable <- renderUI({
      #
      #  if (input$viewLassoMinOrSE == "Lambda Min" || input$viewLassoMinOrSE == "Lambda SE") {
      #    div(DT::dataTableOutput("LassoCoefTable"), style = "font-size:12px")
      #  }
      #
      #})
      #
      #output$rendLassoTrainHRtab <- renderUI({
      #  div(shinycssloaders::withSpinner(tableOutput("LassoTrainHRtab"), type = 7, size = 0.5), style = "font-size:12px")
      #})
      #output$rendLassoTestHRtab <- renderUI({
      #  div(shinycssloaders::withSpinner(tableOutput("LassoTestHRtab"), type = 7, size = 0.5), style = "font-size:12px")
      #})
      #
      #
      #lasso_runData <- eventReactive(input$RunLassoModelGen, {
      #
      #  if (length(input$LassoFeatureSelection_lasso) > 1 ) {
      #
      #    req(metaSub())
      #    req(exprSub())
      #
      #    SeedSelected <- input$LassoSeedSelection
      #    SampleType <- input$SampleTypeSelection_lasso
      #    Feature <- input$FeatureSelection_lasso
      #    SubFeature <- input$SubFeatureSelection_lasso
      #    LassoFeatures <- input$LassoFeatureSelection_lasso
      #    LassoSurvTimeCol <- input$SurvTimeSelec_lasso
      #    LassoSurvIDCol <- input$SurvIDSelect_lasso
      #    LassoTrainProportion <- input$LassoTrainProp
      #
      #    set.seed(SeedSelected)
      #    meta <- meta_react()
      #    expr <- expr_react()
      #    metacol_sampletype <- metacol_sampletype()
      #
      #    if (isTruthy(SampleType)) {
      #      if (SampleType != "All Sample Types") {
      #        meta <- meta[which(meta[,metacol_sampletype] == SampleType),]
      #      } else {
      #        meta <- meta
      #      }
      #    } else {
      #      meta <- meta
      #    }
      #    if (Feature != "Show all Samples") {
      #      meta <- meta[which(meta[,Feature] == SubFeature),]
      #    } else {
      #      meta <- meta
      #    }
      #    metaSub <- meta
      #    exprSub <- expr[,metaSub[,1]]
      #
      #    if (any(LassoFeatures %in% rownames(expr))) {
      #      expr_feats <- exprSub[which(rownames(exprSub) %in% LassoFeatures),]
      #      expr_feats <- as.data.frame(t(expr_feats))
      #      expr_feats$SampleName <- rownames(expr_feats)
      #    } else {
      #      expr_feats <- data.frame(SampleName = metaSub[,1])
      #    }
      #
      #    if (any(LassoFeatures %in% colnames(metaSub))) {
      #      meta_feats <- metaSub[,which(colnames(metaSub) %in% LassoFeatures), drop = F]
      #    } else {
      #      meta_feats <- data.frame(SampleName = metaSub[,1])
      #    }
      #
      #    Lasso_Score_df <- merge(meta_feats,expr_feats, by.x = colnames(meta_feats)[1], by.y = "SampleName")
      #    Lasso_Score_df <- Lasso_Score_df[complete.cases(Lasso_Score_df),]
      #    colnames(Lasso_Score_df)[1] <- "SampleName"
      #
      #    survData <- metaSub[metaSub[,1] %in% Lasso_Score_df[,1],c(colnames(metaSub)[1],LassoSurvIDCol,LassoSurvTimeCol)]
      #    rownames(survData) <- survData[,1]
      #    survData <- survData[,-1]
      #    colnames(survData) <- c("status","time")
      #    survData <- survData[which(survData$time > 0),]
      #
      #    Lasso_Score_df <- Lasso_Score_df[which(Lasso_Score_df$SampleName %in% rownames(survData)),]
      #
      #    train_num <- round(length(Lasso_Score_df$SampleName) * (LassoTrainProportion/100))
      #    train_samp <- sample(Lasso_Score_df$SampleName,train_num)
      #    test_samp <- setdiff(Lasso_Score_df$SampleName,train_samp)
      #
      #    rownames(Lasso_Score_df) <- Lasso_Score_df$SampleName
      #    Lasso_Score_df <- Lasso_Score_df[,-1, drop = F]
      #
      #    score_train <- as.matrix(Lasso_Score_df[train_samp,])
      #    score_test <- as.matrix(Lasso_Score_df[test_samp,])
      #
      #    survData_test <- as.matrix(survData[test_samp,])
      #    survData_train <- as.matrix(survData[train_samp,])
      #
      #    runData <- list(train_samp = train_samp,
      #                    test_samp = test_samp,
      #                    score_train = score_train,
      #                    score_test = score_test,
      #                    survData_test = survData_test,
      #                    survData_train = survData_train)
      #    runData
      #  }
      #
      #})
      #
      #LassoRun_train_model <- eventReactive(input$RunLassoModelGen, {
      #  #LassoRun_train_model <- reactive({
      #
      #  runData <- lasso_runData()
      #  expr <- expr_react()
      #  score_train <- as.matrix(runData$score_train)
      #  survData_train <- as.matrix(runData$survData_train)
      #  AlphaIn <- input$LassoAlpha
      #
      #  save(list = ls(), file = "~/R/ShinyAppsIO/PATH_SURVEYOR/shiny_env.Rdata", envir = environment())
      #
      #  print(summary(runData))
      #  print(head(score_train,c(5,5)))
      #  print(dim(score_train))
      #  print(head(survData_train,c(5,5)))
      #  print(dim(survData_train))
      #  print(AlphaIn)
      #
      #  model <- cv.glmnet(score_train, survData_train, family = "cox", type.measure = "C", alpha=AlphaIn)
      #
      #  score_train_save <- score_train
      #
      #  rownames(score_train) <- NULL
      #  survData_train_save <- survData_train
      #
      #  survData_train <- survData_train %>% as.data.frame() %>% select(time,status) %>% as.matrix()
      #  survData_train_save <- survData_train_save %>% as.data.frame() %>% select(status,time) %>% as.matrix()
      #  survData_train_save <- apply(survData_train_save,2,as.numeric)
      #
      #  survData_train <- with(survData_train_save, Surv(time, status))
      #
      #
      #  samples <- intersect(rownames(survData_train_save),rownames(score_train_save))
      #  score_train_save <- as.data.frame(score_train_save)
      #  score_train_save <- mutate_all(score_train_save, function(x) as.numeric(as.character(x)))
      #  score_train_save <- as.matrix(score_train_save[samples,])
      #  survData_train_save <- survData_train_save[samples,]
      #  survData_train_save_surv <- with(as.data.frame(survData_train_save), Surv(time, status) )
      #  model <- cv.glmnet(score_train_save, survData_train_save_surv, family = "cox", type.measure = "C", alpha=1)
      #
      #
      #  survData_train <- Surv(time = unname(survData_train_save[,1]),event = unname(survData_train_save[,2]))
      #
      #
      #  survData_train[,1] <- survData_train[,1]/365.25
      #
      #  sparsematrix <- as(score_train_save, "sparseMatrix")
      #
      #  model <- cv.glmnet(score_train, survData_train, family = "cox", type.measure = "C", alpha=1)
      #  model
      #  x <- CoxExample$x
      #  y <- CoxExample$y
      #  fit <- cv.glmnet(x, y, family = "cox")
      #})
      #
      #LassoLambdaPlot_react <- reactive({
      #
      #  model <- LassoRun_train_model()
      #  plot(model)
      #
      #})
      #output$Lasso_LambdaPlot <- renderPlot({
      #
      #  p <- LassoLambdaPlot_react()
      #  p
      #
      #})
      #
      #
      #LassoRun_train_model2 <- eventReactive(input$RunLassoModelGen, {
      #  #LassoRun_train_model <- reactive({
      #
      #  runData <- lasso_runData()
      #  score_train <- as.matrix(runData$score_train)
      #  survData_train <- as.matrix(runData$survData_train)
      #  AlphaIn <- input$LassoAlpha
      #
      #  model2 <- glmnet(score_train, survData_train, family = "cox", type.measure = "C", alpha=AlphaIn)
      #
      #  model2
      #
      #})
      #
      #LassoCoeffPlot_react <- reactive({
      #
      #  model2 <- LassoRun_train_model2()
      #  plot(model2)
      #
      #})
      #output$Lasso_CoeffPlot <- renderPlot({
      #
      #  p <- LassoCoeffPlot_react()
      #  p
      #
      #})
      #
      #LassoRun_Lmin_coef_table <- reactive({
      #
      #  model <- LassoRun_train_model()
      #  l_min <- model$lambda.min
      #  model_coef_min <- coef(model ,s=l_min)
      #  model_coef_min_df <- data.frame(Feature = names(sort(model_coef_min[,1])), Coefficient = unname(sort(model_coef_min[,1])))
      #  model_coef_min_df
      #
      #})
      #LassoRun_Lse_coef_table <- reactive({
      #
      #  model <- LassoRun_train_model()
      #  l_se <- model$lambda.1se
      #  model_coef_se <- coef(model ,s=l_se)
      #  model_coef_se_df <- data.frame(Feature = names(sort(model_coef_se[,1])), Coefficient = unname(sort(model_coef_se[,1])))
      #  model_coef_se_df
      #
      #})
      #
      #output$LassoCoefTable <- DT::renderDataTable({
      #
      #  if (input$viewLassoMinOrSE == "Lambda Min") {
      #    df <- LassoRun_Lmin_coef_table()
      #  }
      #  else if (input$viewLassoMinOrSE == "Lambda SE") {
      #    df <- LassoRun_Lse_coef_table()
      #  }
      #  DT::datatable(df, options = list(paging = F,searching = FALSE), rownames = F)
      #
      #})
      #LassoRun_Lmin_Pred_train <- reactive({
      #
      #  ModelName <- input$LassoModelName
      #  runData <- lasso_runData()
      #  model <- LassoRun_train_model()
      #  l_min <- model$lambda.min
      #  score_test <- runData$score_train
      #  Pred_Lmin_train <- as.data.frame(predict(model, s=l_min, newx=score_test, type = "response"))
      #  Pred_Lmin_train_ScoreName <- paste("LassoLmin_Train_",ModelName,"_RiskScore", sep = "")
      #  colnames(Pred_Lmin_train)[1] <- Pred_Lmin_train_ScoreName
      #  Pred_Lmin_train$SampleName <- rownames(Pred_Lmin_train)
      #  Pred_Lmin_train <- Pred_Lmin_train %>%
      #    relocate(SampleName)
      #  Pred_Lmin_train
      #
      #})
      #LassoRun_Lse_Pred_train <- reactive({
      #
      #  ModelName <- input$LassoModelName
      #  runData <- lasso_runData()
      #  model <- LassoRun_train_model()
      #  l_se <- model$lambda.1se
      #  score_test <- runData$score_train
      #  Pred_Lse_train <- as.data.frame(predict(model, s=l_se, newx=score_test, type = "response"))
      #  Pred_Lse_train_ScoreName <- paste("LassoLse_Train_",ModelName,"_RiskScore", sep = "")
      #  colnames(Pred_Lse_train)[1] <- Pred_Lse_train_ScoreName
      #  Pred_Lse_train$SampleName <- rownames(Pred_Lse_train)
      #  Pred_Lse_train<- Pred_Lse_train %>%
      #    relocate(SampleName)
      #  Pred_Lse_train
      #
      #})
      #LassoRun_LCustom_Pred_train <- reactive({
      #
      #  if (input$viewLassoMinOrSE == "Custom") {
      #    if (!is.na(input$CustomLambda)) {
      #      ModelName <- input$LassoModelName
      #      runData <- lasso_runData()
      #      model <- LassoRun_train_model()
      #      l_Custom <- input$CustomLambda
      #      score_test <- runData$score_train
      #      Pred_LCustom_train <- as.data.frame(predict(model, s=l_Custom, newx=score_test, type = "response"))
      #      Pred_LCustom_train_ScoreName <- paste("LassoCustomLambda_Train_",ModelName,"_RiskScore", sep = "")
      #      colnames(Pred_LCustom_train)[1] <- Pred_LCustom_train_ScoreName
      #      Pred_LCustom_train$SampleName <- rownames(Pred_LCustom_train)
      #      Pred_LCustom_train<- Pred_LCustom_train %>%
      #        relocate(SampleName)
      #      Pred_LCustom_train
      #    }
      #  }
      #
      #})
      #LassoRun_Lmin_Pred_test <- reactive({
      #
      #  ModelName <- input$LassoModelName
      #  runData <- lasso_runData()
      #  model <- LassoRun_train_model()
      #  l_min <- model$lambda.min
      #  score_test <- runData$score_test
      #  Pred_Lmin_test <- as.data.frame(predict(model, s=l_min, newx=score_test, type = "response"))
      #  Pred_Lmin_test_ScoreName <- paste("LassoLmin_Test_",ModelName,"_RiskScore", sep = "")
      #  colnames(Pred_Lmin_test)[1] <- Pred_Lmin_test_ScoreName
      #  Pred_Lmin_test$SampleName <- rownames(Pred_Lmin_test)
      #  Pred_Lmin_test <- Pred_Lmin_test %>%
      #    relocate(SampleName)
      #  Pred_Lmin_test
      #
      #})
      #LassoRun_Lse_Pred_test <- reactive({
      #
      #  ModelName <- input$LassoModelName
      #  runData <- lasso_runData()
      #  model <- LassoRun_train_model()
      #  l_se <- model$lambda.1se
      #  score_test <- runData$score_test
      #  Pred_Lse_test <- as.data.frame(predict(model, s=l_se, newx=score_test, type = "response"))
      #  Pred_Lse_test_ScoreName <- paste("LassoLse_Test_",ModelName,"_RiskScore", sep = "")
      #  colnames(Pred_Lse_test)[1] <- Pred_Lse_test_ScoreName
      #  Pred_Lse_test$SampleName <- rownames(Pred_Lse_test)
      #  Pred_Lse_test <- Pred_Lse_test %>%
      #    relocate(SampleName)
      #  Pred_Lse_test
      #
      #})
      #LassoRun_LCustom_Pred_test <- reactive({
      #
      #  if (input$viewLassoMinOrSE == "Custom") {
      #    if (!is.na(input$CustomLambda)) {
      #      ModelName <- input$LassoModelName
      #      runData <- lasso_runData()
      #      model <- LassoRun_train_model()
      #      l_Custom <- input$CustomLambda
      #      score_test <- runData$score_test
      #      Pred_LCustom_test <- as.data.frame(predict(model, s=l_Custom, newx=score_test, type = "response"))
      #      Pred_LCustom_test_ScoreName <- paste("LassoCustomLambda_Test_",ModelName,"_RiskScore", sep = "")
      #      colnames(Pred_LCustom_test)[1] <- Pred_LCustom_test_ScoreName
      #      Pred_LCustom_test$SampleName <- rownames(Pred_LCustom_test)
      #      Pred_LCustom_test <- Pred_LCustom_test %>%
      #        relocate(SampleName)
      #      Pred_LCustom_test
      #    }
      #  }
      #
      #})
      #
      #LassoDnldTable <- reactive({
      #
      #  runData <- lasso_runData()
      #  RiskScore_List <- list()
      #  RiskScore_LminTrain <- LassoRun_Lmin_Pred_train()
      #  RiskScore_List[["RiskScore_LminTrain"]] <- RiskScore_LminTrain
      #  RiskScore_LseTrain <- LassoRun_Lse_Pred_train()
      #  RiskScore_List[["RiskScore_LseTrain"]] <- RiskScore_LseTrain
      #  RiskScore_LminTest <- LassoRun_Lmin_Pred_test()
      #  RiskScore_List[["RiskScore_LminTest"]] <- RiskScore_LminTest
      #  RiskScore_LseTest <- LassoRun_Lse_Pred_test()
      #  RiskScore_List[["RiskScore_LseTest"]] <- RiskScore_LseTest
      #  if (input$viewLassoMinOrSE == "Custom") {
      #    RiskScore_LCustomTrain <- LassoRun_LCustom_Pred_train()
      #    RiskScore_List[["RiskScore_LCustomTrain"]] <- RiskScore_LCustomTrain
      #    RiskScore_LCustomTest <- LassoRun_LCustom_Pred_test()
      #    RiskScore_List[["RiskScore_LCustomTest"]] <- RiskScore_LCustomTest
      #  }
      #  RiskScore_df <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "SampleName", all = TRUE),
      #                         RiskScore_List)
      #
      #  RiskScore_df$TrainOrTest <- ifelse(RiskScore_df$SampleName %in% runData$train_samp, "Training", "Testing")
      #  RiskScore_df <- RiskScore_df %>%
      #    relocate(SampleName,TrainOrTest)
      #  RiskScore_df
      #
      #})
      #
      #Lasso_Train_Surv_df <- reactive({
      #
      #  ## Assign variables
      #  CutPoption <- input$LassoPlotCutP
      #  surv_time_col <- input$SurvTimeSelecView_lasso
      #  surv_id_col <- input$SurvIDSelecView_lasso
      #  userCutP <- input$CutPinput
      #  if (input$viewLassoMinOrSE == "Lambda Min") {
      #    Pred_train_df <- LassoRun_Lmin_Pred_train()
      #  }
      #  else if (input$viewLassoMinOrSE == "Lambda SE") {
      #    Pred_train_df <- LassoRun_Lse_Pred_train()
      #  }
      #  else if (input$viewLassoMinOrSE == "Custom") {
      #    req(input$CustomLambda)
      #    Pred_train_df <- LassoRun_LCustom_Pred_train()
      #  }
      #
      #  meta_surv <- meta[,c("SampleName",surv_time_col,surv_id_col)]
      #
      #  KP_df_train <- merge(Pred_train_df,meta_surv, by = "SampleName")
      #
      #  ## Subset columns needed for plot
      #  colnames(KP_df_train)[which(colnames(KP_df_train) == surv_time_col)] <- "time"
      #  colnames(KP_df_train)[which(colnames(KP_df_train) == surv_id_col)] <- "ID"
      #
      #  if (CutPoption == "Median") {
      #    KP_df_train$MedianCutP <- highlow(KP_df_train[,2])
      #  }
      #  else if (CutPoption == "Quartile") {
      #    KP_df_train$VAR_Q <- quartile_conversion(KP_df_train[,2])
      #    KP_df_train$QuartileCutP <- paste("", KP_df_train$VAR_Q, sep="")
      #  }
      #  else if (CutPoption == "Optimal") {
      #    if (length(KP_df_train[,2][KP_df_train[,2] > 0])/length(KP_df_train[,2]) > 0.01) {
      #      if (length(KP_df_train[,4]) > 1) {
      #        res.cut <- survminer::surv_cutpoint(KP_df_train,time = "time", event = "ID", variable = colnames(KP_df_train)[2], minprop = 0.01)
      #        cutp <- res.cut$cutpoint[["cutpoint"]]
      #        res.cat <- surv_categorize(res.cut)
      #        KP_df_train$OptimalCutP <- res.cat[,3]
      #      }
      #    }
      #  }
      #  else if (CutPoption == "Quantile") {
      #    userCutP <- userCutP/100
      #    KP_df_train$TopBottomCutP <- quantile_conversion(KP_df_train[,2], userCutP)
      #    KP_df_train <- KP_df_train[which(KP_df_train$TopBottomCutP != "BetweenCutoff"),]
      #  }
      #  else if (CutPoption == "User Specified") {
      #    userCutP <- userCutP/100
      #    KP_df_train$UserCutP <- quantile_conversion2(KP_df_train[,2], userCutP)
      #  }
      #
      #  KP_df_train
      #
      #})
      #Lasso_Test_Surv_df <- reactive({
      #
      #  req(metaSub())
      #  req(exprSub())
      #  metaSub <- metaSub()
      #  expr <- exprSub()
      #
      #  ## Assign variables
      #  CutPoption <- input$LassoPlotCutP
      #  surv_time_col <- input$SurvTimeSelecView_lasso
      #  surv_id_col <- input$SurvIDSelecView_lasso
      #  userCutP <- input$CutPinput
      #  if (input$viewLassoMinOrSE == "Lambda Min") {
      #    Pred_test_df <- LassoRun_Lmin_Pred_test()
      #  }
      #  else if (input$viewLassoMinOrSE == "Lambda SE") {
      #    Pred_test_df <- LassoRun_Lse_Pred_test()
      #  }
      #  else if (input$viewLassoMinOrSE == "Custom") {
      #    req(input$CustomLambda)
      #    Pred_test_df <- LassoRun_LCustom_Pred_test()
      #  }
      #
      #
      #
      #  meta_surv <- meta[,c("SampleName",surv_time_col,surv_id_col)]
      #
      #  KP_df_test <- merge(Pred_test_df,meta_surv, by = "SampleName")
      #
      #  ## Subset columns needed for plot
      #  colnames(KP_df_test)[which(colnames(KP_df_test) == surv_time_col)] <- "time"
      #  colnames(KP_df_test)[which(colnames(KP_df_test) == surv_id_col)] <- "ID"
      #
      #  if (CutPoption == "Median") {
      #    KP_df_test$MedianCutP <- highlow(KP_df_test[,2])
      #  }
      #  else if (CutPoption == "Quartile") {
      #    KP_df_test$VAR_Q <- quartile_conversion(KP_df_test[,2])
      #    KP_df_test$QuartileCutP <- paste("", KP_df_test$VAR_Q, sep="")
      #  }
      #  else if (CutPoption == "Optimal") {
      #    if (length(KP_df_test[,2][KP_df_test[,2] > 0])/length(KP_df_test[,2]) > 0.01) {
      #      if (length(KP_df_test[,4]) > 1) {
      #        res.cut <- survminer::surv_cutpoint(KP_df_test,time = "time", event = "ID", variable = colnames(KP_df_test)[2], minprop = 0.01)
      #        cutp <- res.cut$cutpoint[["cutpoint"]]
      #        res.cat <- surv_categorize(res.cut)
      #        KP_df_test$OptimalCutP <- res.cat[,3]
      #      }
      #    }
      #  }
      #  else if (CutPoption == "Quantile") {
      #    userCutP <- userCutP/100
      #    KP_df_test$TopBottomCutP <- quantile_conversion(KP_df_test[,2], userCutP)
      #    KP_df_test <- KP_df_test[which(KP_df_test$TopBottomCutP != "BetweenCutoff"),]
      #  }
      #  else if (CutPoption == "User Specified") {
      #    userCutP <- userCutP/100
      #    KP_df_test$UserCutP <- quantile_conversion2(KP_df_test[,2], userCutP)
      #  }
      #
      #  KP_df_test
      #
      #})
      #
      #Lasso_Train_Surv_Tab_react <- reactive({
      #
      #  KP_df_train <- Lasso_Train_Surv_df()
      #  CutPoption <- input$LassoPlotCutP
      #  KP_df_train[,5] <- as.factor(KP_df_train[,5])
      #  KP_df_train[,5] <- relevel(KP_df_train[,5],
      #                             ref = grep("low",unique(KP_df_train[,5]),
      #                                        value = T, ignore.case = T))
      #
      #  ## Survival Function
      #  ## Survival Function
      #  tab_train <- coxph(as.formula(paste("Surv(time,ID) ~ ",colnames(KP_df_train)[5],sep = "")),
      #                     data = KP_df_train)
      #  tab_train
      #
      #})
      #Lasso_Test_Surv_Tab_react <- reactive({
      #
      #  KP_df_test <- Lasso_Test_Surv_df()
      #  CutPoption <- input$LassoPlotCutP
      #  KP_df_test[,5] <- as.factor(KP_df_test[,5])
      #  KP_df_test[,5] <- relevel(KP_df_test[,5],
      #                            ref = grep("low",unique(KP_df_test[,5]),
      #                                       value = T, ignore.case = T))
      #
      #  ## Survival Function
      #  tab_test <- coxph(as.formula(paste("Surv(time,ID) ~ ",colnames(KP_df_test)[5],sep = "")),
      #                    data = KP_df_test)
      #  tab_test
      #
      #})
      #
      #output$LassoTrainCoxSumm <- renderPrint({
      #
      #  tab <- Lasso_Train_Surv_Tab_react()
      #  out <- capture.output(summary(tab))
      #
      #  con_line <- grep("^Concordance=",out,value = T)
      #  lik_line <- grep("^Likelihood ratio test=",out,value = T)
      #  wal_line <- grep("^Wald test",out,value = T)
      #  sco_line <- grep("^Score ",out,value = T)
      #
      #  text <- paste("CoxH Summary:",con_line,lik_line,wal_line,sco_line,sep = "\n")
      #  cat(text)
      #
      #})
      #output$LassoTestCoxSumm <- renderPrint({
      #
      #  tab <- Lasso_Test_Surv_Tab_react()
      #  out <- capture.output(summary(tab))
      #
      #  con_line <- grep("^Concordance=",out,value = T)
      #  lik_line <- grep("^Likelihood ratio test=",out,value = T)
      #  wal_line <- grep("^Wald test",out,value = T)
      #  sco_line <- grep("^Score ",out,value = T)
      #
      #  text <- paste("CoxH Summary:",con_line,lik_line,wal_line,sco_line,sep = "\n")
      #  cat(text)
      #
      #})
      #
      #
      #Lasso_Train_Surv_HRTab_react <- reactive({
      #
      #  tab_train <- Lasso_Train_Surv_Tab_react()
      #  tab_train <- tab_train %>%
      #    gtsummary::tbl_regression(exp = TRUE) %>%
      #    as_gt()
      #
      #  tab_train_df <- as.data.frame(tab_train)
      #
      #  tab_train_df <- tab_train_df %>%
      #    dplyr::select(label,estimate,ci,p.value)
      #  colnames(tab_train_df) <- c("Characteristic","Hazard Ratio","95% Confidence Interval","P.Value")
      #
      #  tab_train_df
      #
      #})
      #Lasso_Test_Surv_HRTab_react <- reactive({
      #
      #  tab_test <- Lasso_Test_Surv_Tab_react()
      #  tab_test <- tab_test %>%
      #    gtsummary::tbl_regression(exp = TRUE) %>%
      #    as_gt()
      #
      #  tab_test_df <- as.data.frame(tab_test)
      #
      #  tab_test_df <- tab_test_df %>%
      #    dplyr::select(label,estimate,ci,p.value)
      #  colnames(tab_test_df) <- c("Characteristic","Hazard Ratio","95% Confidence Interval","P.Value")
      #
      #  tab_test_df
      #
      #})
      #
      #output$LassoTrainHRtab <- renderTable({
      #
      #  tab <- Lasso_Train_Surv_HRTab_react()
      #  tab
      #
      #})
      #output$LassoTestHRtab <- renderTable({
      #
      #  tab <- Lasso_Test_Surv_HRTab_react()
      #  tab
      #
      #})
      #
      #Lasso_Train_Splot_react <- reactive({
      #
      #  ## Assign variables
      #  KP_df_train <- Lasso_Train_Surv_df()
      #  #LassoModelName <- colnames(KP_df_train)[2]
      #  LassoModelName <- input$LassoModelName
      #  LambdaChoice <- input$viewLassoMinOrSE
      #  CutPoption <- input$LassoPlotCutP
      #  SampleType <- input$SampleTypeSelection_lasso
      #  Feature <- input$FeatureSelection_lasso
      #  show_pval <- input$ShowPval_lasso
      #  ShowConfInt <- input$ShowConfInt_lasso
      #  xaxlim <- input$SurvXaxis_lasso * 365.25
      #  surv_time_col <- input$SurvTimeSelecView_lasso
      #  showLegend <- input$SurvLegendPos_lasso
      #  showMedSurv <- input$ShowMedSurvLine_lasso
      #  LambdaSelect <- input$viewLassoMinOrSE
      #  if (showMedSurv == T) {
      #    showMedSurv <- "hv"
      #  }
      #  else if (showMedSurv == F) {
      #    showMedSurv <- "none"
      #  }
      #
      #  Feature <- colnames(KP_df_train)[5]
      #
      #  form <- paste("Surv(time,ID) ~ ",Feature,sep = "")
      #  form2 <- as.formula(form)
      #  fit_train <- eval(substitute(survfit(form2,data = KP_df_train, type="kaplan-meier")))
      #
      #  ## Survival Function
      #  #fit_train <- survfit(Surv(time,ID) ~ FeatureCut, data = KP_df_train, type="kaplan-meier")
      #
      #  ## Determine type of survival data - OS/EFS/PFS?
      #  SurvDateType <- sub("\\..*","",surv_time_col)
      #
      #
      #
      #  ### determine Feature and Sample Type label
      #  #if (length(unique(meta[,metacol_sampletype])) > 1) {
      #  #  if (SampleType == "Show All Sample Types") {
      #  #    if (Feature == "Show All Samples") {
      #  #      SampleTypeLab <- "All Features in All Patients\n"
      #  #    }
      #  #    if (Feature != "Show All Samples") {
      #  #      SampleTypeLab <- paste(Feature," in All Patients\n")
      #  #    }
      #  #  }
      #  #  else {
      #  #    if (Feature == "Show All Samples") {
      #  #      SampleTypeLab <- paste("All Features (",SampleType,") Patients\n",sep = "")
      #  #    }
      #  #    if (Feature != "Show All Samples") {
      #  #      SampleTypeLab <- paste(Feature," (",SampleType,") Patients\n",sep = "")
      #  #    }
      #  #  }
      #  #}
      #  #if (length(unique(meta[,metacol_sampletype])) <= 1) {
      #  #  if (Feature == "Show All Samples") {
      #  #    SampleTypeLab <- "All Features in All Patients\n"
      #  #  }
      #  #  if (Feature != "Show All Samples") {
      #  #    SampleTypeLab <- paste(Feature," in All Patients\n")
      #  #  }
      #  #}
      #
      #  CutPMethodLab <- paste(CutPoption, " Cut-Point", sep = "")
      #  if (LambdaChoice == "Custom") {
      #    LambdaChoice <- "Custom Lambda"
      #  }
      #
      #  ## Determine Plot title
      #  if (is.null(input$SurvPlotTitleLasso)) {
      #    SurvPlotTitle <- paste("Survival curves of Training Data (",LambdaChoice,")\n",
      #                           LassoModelName," (",CutPMethodLab,")", sep = "")
      #  }
      #  else if (!is.null(input$SurvPlotTitleLasso)) {
      #    if (input$SurvPlotTitleMedian == "") {
      #      SurvPlotTitle <- paste("Survival curves of Training Data (",LambdaChoice,")\n",
      #                             LassoModelName," (",CutPMethodLab,")", sep = "")
      #    }
      #    else if (input$SurvPlotTitleLasso != "") {
      #      SurvPlotTitle <- input$SurvPlotTitleLasso
      #    }
      #  }
      #
      #  breakTime <- 365.25
      #  if (max(KP_df_train[,"time"]) < 365.25) {
      #    breakTime <- NULL
      #  }
      #
      #  ## Generate plot
      #  ggsurv <- survminer::ggsurvplot(fit_train, data = KP_df_train, risk.table = TRUE,
      #                                  title = SurvPlotTitle,
      #                                  xscale = c("d_y"),
      #                                  break.time.by=breakTime,
      #                                  xlab = "Years",
      #                                  ylab = paste(SurvDateType,"Survival Probability"),
      #                                  submain = "Based on Kaplan-Meier estimates",
      #                                  caption = "created with survminer",
      #                                  pval=show_pval,
      #                                  conf.int = ShowConfInt,
      #                                  ggtheme = theme_bw(),
      #                                  font.title = c(16, "bold"),
      #                                  font.submain = c(12, "italic"),
      #                                  font.caption = c(12, "plain"),
      #                                  font.x = c(14, "plain"),
      #                                  font.y = c(14, "plain"),
      #                                  font.tickslab = c(12, "plain"),
      #                                  legend = showLegend,
      #                                  risk.table.height = 0.20,
      #                                  surv.median.line = showMedSurv
      #  )
      #  if (showMedSurv != "none") {
      #    MedSurvItem <- ggsurv[["plot"]][["layers"]][length(ggsurv[["plot"]][["layers"]])]
      #    MedSurvItem_df <- MedSurvItem[[1]][["data"]]
      #    MedSurvItem_df <- MedSurvItem_df[order(MedSurvItem_df[,1]),]
      #    MedSurvItem_df <- MedSurvItem_df %>%
      #      mutate(label = paste(round(MedSurvItem_df[,1]),"Days"))
      #    rownames(MedSurvItem_df) <- 1:nrow(MedSurvItem_df)
      #    if (nrow(MedSurvItem_df) > 1) {
      #      ggsurv$plot <- ggsurv$plot +
      #        geom_label_repel(data = MedSurvItem_df, aes(x = x1, y = y1, label = label, size = 4), label.size = NA, show.legend = FALSE)
      #    }
      #  }
      #  if (isTruthy(input$SurvXaxis)) {
      #    ggsurv$plot$coordinates$limits$x <- c(0,xaxlim)
      #    ggsurv$table$coordinates$limits$x <- c(0,xaxlim)
      #  }
      #
      #  ggsurv$table <- ggsurv$table + theme_cleantable()
      #  ggsurv
      #
      #
      #})
      #
      #Lasso_Test_Splot_react <- reactive({
      #
      #  ## Assign variables
      #  KP_df_test <- Lasso_Test_Surv_df()
      #  #LassoModelName <- colnames(KP_df_test)[2]
      #  LassoModelName <- input$LassoModelName
      #  CutPoption <- input$LassoPlotCutP
      #  LambdaChoice <- input$viewLassoMinOrSE
      #  SampleType <- input$SampleTypeSelection_lasso
      #  Feature <- input$FeatureSelection_lasso
      #  show_pval <- input$ShowPval_lasso
      #  ShowConfInt <- input$ShowConfInt_lasso
      #  xaxlim <- input$SurvXaxis_lasso * 365.25
      #  surv_time_col <- input$SurvTimeSelecView_lasso
      #  showLegend <- input$SurvLegendPos_lasso
      #  showMedSurv <- input$ShowMedSurvLine_lasso
      #  LambdaSelect <- input$viewLassoMinOrSE
      #  if (showMedSurv == T) {
      #    showMedSurv <- "hv"
      #  }
      #  else if (showMedSurv == F) {
      #    showMedSurv <- "none"
      #  }
      #
      #  Feature <- colnames(KP_df_test)[5]
      #
      #  form <- paste("Surv(time,ID) ~ ",Feature,sep = "")
      #  form2 <- as.formula(form)
      #  fit_test <- eval(substitute(survfit(form2,data = KP_df_test, type="kaplan-meier")))
      #
      #  ## Survival Function
      #  #fit_train <- survfit(Surv(time,ID) ~ FeatureCut, data = KP_df_train, type="kaplan-meier")
      #
      #  ## Determine type of survival data - OS/EFS/PFS?
      #  SurvDateType <- sub("\\..*","",surv_time_col)
      #
      #  ### determine Feature and Sample Type label
      #  #if (length(unique(meta[,metacol_sampletype])) > 1) {
      #  #  if (SampleType == "Show All Sample Types") {
      #  #    if (Feature == "Show All Samples") {
      #  #      SampleTypeLab <- "All Features in All Patients\n"
      #  #    }
      #  #    if (Feature != "Show All Samples") {
      #  #      SampleTypeLab <- paste(Feature," in All Patients\n")
      #  #    }
      #  #  }
      #  #  else {
      #  #    if (Feature == "Show All Samples") {
      #  #      SampleTypeLab <- paste("All Features (",SampleType,") Patients\n",sep = "")
      #  #    }
      #  #    if (Feature != "Show All Samples") {
      #  #      SampleTypeLab <- paste(Feature," (",SampleType,") Patients\n",sep = "")
      #  #    }
      #  #  }
      #  #}
      #  #if (length(unique(meta[,metacol_sampletype])) <= 1) {
      #  #  if (Feature == "Show All Samples") {
      #  #    SampleTypeLab <- "All Features in All Patients\n"
      #  #  }
      #  #  if (Feature != "Show All Samples") {
      #  #    SampleTypeLab <- paste(Feature," in All Patients\n")
      #  #  }
      #  #}
      #
      #  CutPMethodLab <- paste(CutPoption, " Cut-Point", sep = "")
      #
      #
      #  ## Determine Plot title
      #  if (is.null(input$SurvPlotTitleLasso)) {
      #    SurvPlotTitle <- paste("Survival curves of Testing Data (",LambdaChoice,")\n",
      #                           LassoModelName," (",CutPMethodLab,")", sep = "")
      #  }
      #  else if (!is.null(input$SurvPlotTitleLasso)) {
      #    if (input$SurvPlotTitleMedian == "") {
      #      SurvPlotTitle <- paste("Survival curves of Testing Data (",LambdaChoice,")\n",
      #                             LassoModelName," (",CutPMethodLab,")", sep = "")
      #    }
      #    else if (input$SurvPlotTitleLasso != "") {
      #      SurvPlotTitle <- input$SurvPlotTitleLasso
      #    }
      #  }
      #
      #  breakTime <- 365.25
      #  if (max(KP_df_test[,"time"]) < 365.25) {
      #    breakTime <- NULL
      #  }
      #
      #  ## Generate plot
      #  ggsurv <- survminer::ggsurvplot(fit_test, data = KP_df_test, risk.table = TRUE,
      #                                  title = SurvPlotTitle,
      #                                  xscale = c("d_y"),
      #                                  break.time.by=breakTime,
      #                                  xlab = "Years",
      #                                  ylab = paste(SurvDateType,"Survival Probability"),
      #                                  submain = "Based on Kaplan-Meier estimates",
      #                                  caption = "created with survminer",
      #                                  pval=show_pval,
      #                                  conf.int = ShowConfInt,
      #                                  ggtheme = theme_bw(),
      #                                  font.title = c(16, "bold"),
      #                                  font.submain = c(12, "italic"),
      #                                  font.caption = c(12, "plain"),
      #                                  font.x = c(14, "plain"),
      #                                  font.y = c(14, "plain"),
      #                                  font.tickslab = c(12, "plain"),
      #                                  legend = showLegend,
      #                                  risk.table.height = 0.20,
      #                                  surv.median.line = showMedSurv
      #  )
      #  if (showMedSurv != "none") {
      #    MedSurvItem <- ggsurv[["plot"]][["layers"]][length(ggsurv[["plot"]][["layers"]])]
      #    MedSurvItem_df <- MedSurvItem[[1]][["data"]]
      #    MedSurvItem_df <- MedSurvItem_df[order(MedSurvItem_df[,1]),]
      #    MedSurvItem_df <- MedSurvItem_df %>%
      #      mutate(label = paste(round(MedSurvItem_df[,1]),"Days"))
      #    rownames(MedSurvItem_df) <- 1:nrow(MedSurvItem_df)
      #    if (nrow(MedSurvItem_df) > 1) {
      #      ggsurv$plot <- ggsurv$plot +
      #        geom_label_repel(data = MedSurvItem_df, aes(x = x1, y = y1, label = label, size = 4), label.size = NA, show.legend = FALSE)
      #    }
      #  }
      #  if (isTruthy(input$SurvXaxis)) {
      #    ggsurv$plot$coordinates$limits$x <- c(0,xaxlim)
      #    ggsurv$table$coordinates$limits$x <- c(0,xaxlim)
      #  }
      #
      #  ggsurv$table <- ggsurv$table + theme_cleantable()
      #  ggsurv
      #
      #
      #})
      #
      #output$Lasso_Train_Splot <- renderPlot({
      #  plot <- Lasso_Train_Splot_react()
      #  plot
      #})
      #output$Lasso_Test_Splot <- renderPlot({
      #  plot <- Lasso_Test_Splot_react()
      #  plot
      #})
      
      # Data Exploration -------------------------------------------------------
      
      ## Meta Exploration ------------------------------------------------------
      
      observe({
        req(ssGSEAmeta())
        meta <- ssGSEAmeta()
        MetaColChoices <- colnames(meta)[c(2:ncol(meta))]
        updateSelectizeInput(session,"MetaTableCols",choices = MetaColChoices, selected = "",server = T)
      })
      
      output$MetaTable <- DT::renderDataTable({
        req(ssGSEAmeta())
        meta <- ssGSEAmeta()
        surv_time_col <- input$SurvivalType_time
        surv_id_col <- input$SurvivalType_id
        Feature <- input$FeatureSelection
        userMetaCols <- input$MetaTableCols
        metaCols <- colnames(meta)[1] #select sample name column automatically
        if (Feature == "Show all Samples") {
          metaCols <- c(metaCols,surv_time_col,surv_id_col,userMetaCols) #combine column names selected
        }
        else if (Feature != "Show all Samples") {
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
      
      ## Density ---------------------------------------------------------------
      
      ssgseaDensity_react <- reactive({
        req(ssGSEAmeta())
        geneset <- gs_react()
        geneset_name <- names(geneset)
        ssgsea_meta <- ssGSEAmeta()
        user_quant <- input$densityPercent/100
        ShowQuartile <- input$QuartileLinesCheck
        scoreMethod <- input$ScoreMethod
        decon_score_cols <- decon_score_cols()
        
        cols_selec <- c(colnames(ssgsea_meta)[1],geneset_name)
        ssgsea_scores <- ssgsea_meta[,cols_selec]
        
        quant_df <- data.frame(quantile(ssgsea_scores[,geneset_name],na.rm = T))
        quant_df <- quant_df[c(2,3,4),,drop = F]
        colnames(quant_df)[1] <- "Quantile"
        
        user_vline <- quantile(ssgsea_scores[,geneset_name],probs = user_quant,na.rm = T)
        
        ## get score method for x and y labels
        if (input$GeneSetTabs == 2) {
          ylab <- "Gene Expression Density"
          xlab <- "Gene Expression (Log(exp+1))"
          ssgsea_scores[,geneset_name] <- log(ssgsea_scores[,geneset_name] + 1)
          quant_df$Quantile <- round(log(quant_df$Quantile + 1),3)
        } else if (input$GeneSetTabs != 2) {
          if (geneset_name %in% decon_score_cols) {
            ylab <- "Pre-Processed Score Density"
            xlab <- "Pre-Processed Score"
          } else {
            ylab <- paste(scoreMethod, " Score Density", sep = "")
            xlab <- paste(scoreMethod, " Score", sep = "")
          }
        }
        ## generate title based on input
        if (geneset_name %in% decon_score_cols) {
          dens_title <- paste(geneset_name," Pre-Processed Score Density",sep = "")
        } else {
          if (input$GeneSetTabs == 2) {
            dens_title <- paste(geneset_name,"Log(exp+1)",ylab)
          } else {
            dens_title <- paste(geneset_name,ylab)
          }
        }
        
        p <- densPlot(ssgsea_scores[,geneset_name,drop = F],quant_df,xlab,ylab,dens_title,CutPlabel = NULL,ShowQuartile = ShowQuartile,user_vline = user_vline)
        p
        
      })
      output$ssgseaDensity <- renderPlot({
        p <- ssgseaDensity_react()
        p
      })
      output$ssgseaDensityTable <- renderDataTable({
        geneset <- gs_react()
        GeneSet <- names(geneset)
        req(ssGSEAmeta())
        ssgsea_meta <- ssGSEAmeta()
        table <- ssgsea_meta[,c(colnames(ssgsea_meta)[1],GeneSet)]
        DT::datatable(table,
                      options = list(scrollY = T),
                      rownames = F)
      })
      
      ## Feature Comp ----------------------------------------------------------
      
      observe({
        Features <- c(metacol_feature(),metacol_survtime(),metacol_survid())
        updateSelectizeInput(session,"ScatterFeature",choices = Features,selected = input$SurvivalType_time, server = T)
      })
      observe({
        req(input$ColorScatterChoice)
        if (input$ColorScatterChoice == "Feature") {
          Features <- c(metacol_feature(),metacol_survtime(),metacol_survid())
          updateSelectizeInput(session,"ScatterColor",choices = Features,selected = input$SurvivalType_id, server = T)
        }
      })
      
      FeatCompScatter_react <- reactive({
        req(ssGSEAmeta())
        req(gs_react())
        req(input$ScatterFeature)
        geneset <- gs_react()
        geneset_name <- names(geneset)
        Feature <- input$ScatterFeature
        ColorCol <- input$ScatterColor
        ColorCol2 <- input$ScatterColor2
        LogChoice <- input$ScatterLog
        meta <- ssGSEAmeta()
        scores <- NULL
        if (input$ColorScatterChoice == "Feature") {
          if (ColorCol %in% colnames(meta)) {
            scores <- meta[,c(colnames(meta)[1],Feature,geneset_name,ColorCol)]
          }
        } else if (input$ColorScatterChoice == "Single Color") {
          scores <- meta[,c(colnames(meta)[1],Feature,geneset_name)]
          scores$ColorColumn <- ColorCol2
        }
        if (isTruthy(scores)) {
          if (ncol(scores) == 4) {
            scores[,4] <- as.factor(scores[,4])
            
            if ("Log x-axis" %in% LogChoice) {
              scores[,2] <- log(scores[,2] + 1)
            }
            if ("Log y-axis" %in% LogChoice) {
              scores[,3] <- log(scores[,3] + 1)
            }
            
            scores
          }
        }
      })
      
      FeatCompScatterPlot_react <- reactive({
        
        req(FeatCompScatter_react())
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
        if ("Log x-axis" %in% LogChoice) {
          Feature <- paste(Feature,"(log(Feature + 1))")
        }
        if ("Log y-axis" %in% LogChoice) {
          scoreMethod <- paste(scoreMethod,"(log(score + 1))")
        }
        
        if (is.numeric(scores[,2])) {
          scores[,2] <- round(scores[,2],4)
        }
        if (is.numeric(scores[,3])) {
          scores[,3] <- round(scores[,3],4)
        }
        
        if (input$ColorScatterChoice == "Feature") {
          p <- ggplot(scores, aes(x = scores[,2], y = scores[,3],
                                  color = scores[,4],
                                  text = paste("</br> Sample Name: ", scores[,1],
                                               "</br> ",colnames(scores)[2],": ",scores[,2],
                                               "</br> ",colnames(scores)[3],": ",scores[,3],
                                               sep =""))) +
            geom_point() +
            theme(legend.position = "none") +
            labs(color = colnames(scores)[4])
        }
        if (input$ColorScatterChoice == "Single Color") {
          p <- ggplot(scores, aes(x = scores[,2], y = scores[,3],
                                  text = paste("</br> Sample Name: ", scores[,1],
                                               "</br> ",colnames(scores)[2],": ",scores[,2],
                                               "</br> ",colnames(scores)[3],": ",scores[,3],
                                               sep =""))) +
            geom_point(color = scores[1,4]) +
            theme(legend.position = "none")
        }
        p <- p +
          theme_minimal() +
          xlab(Feature) +
          ylab(scoreMethod) +
          theme(plot.margin = margin(2, 3, 0, 0, "cm"))
        p
        
      })
      
      output$FeatCompScatterPlot <- renderPlotly({
        
        geneset_name <- names(gs_react())
        Feature <- input$ScatterFeature
        scoreMethod <- input$ScoreMethod
        scores <- FeatCompScatter_react()
        regLine <- input$RegressionLine
        p <- FeatCompScatterPlot_react()
        p <- ggplotly(p,tooltip = "text")
        p <- p %>%
          config(
            toImageButtonOptions = list(
              format = "svg",
              height = input$PlotDnldHight,
              width = input$PlotDnldWidth,
              filename = paste0(ProjectName_react(),"_",Feature,"_vs_",geneset_name)
            )
          )
        plot_df_sub2 <- scores[,c(1,2,3)]
        colnames(plot_df_sub2) <- c("SampleName","xVar","yVar")
        ScatterTitle_in <- paste(Feature,"vs.",geneset_name,scoreMethod)
        
        ## Generate Regression Line
        if (length(plot_df_sub2$yVar) > 0 & length(plot_df_sub2$xVar) > 0) {
          reg = lm(plot_df_sub2$yVar ~ plot_df_sub2$xVar)
          R2 = summary(reg)$r.squared
          if (regLine == T) {
            xdf <- data.frame(plot_df_sub2$xVar)
            colnames(xdf) <- c('xVar')
            ydf <- reg %>% predict(xdf) %>% data.frame()
            colnames(ydf) <- c('yVar')
            xy <- data.frame(xdf, ydf)
            p <- p %>%
              add_lines(data = xy, x = ~xVar, y = ~yVar, name = "Regression Fit",
                        line = list(color = "black", width=2, dash="dash"))
          }
          # Add title and subtitle
          coef <- paste0("Coefficients: y = ",round(coefficients(reg)[2],3),"x"," + ",round(coefficients(reg)[1],3))
          rSqu <- paste0("R-Squared: ",R2)
          p <- p %>% layout(title = list(text = paste0(ScatterTitle_in,
                                                       '<br>',
                                                       '<sup>',
                                                       coef,
                                                       '<br>',
                                                       rSqu,
                                                       '</sup>'),
                                         x = 0,
                                         xref='paper',
                                         yref='paper',
                                         align = "left"
          )
          )
          p
        }
        
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
      
      ## Risk Strat ------------------------------------------------------------
      ## Boxplot reactive to determine high/low risk samples based on user input
      SboxplotReact <- reactive({
        
        ## Assing variables
        req(ssGSEAmeta())
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
        ssGSEA_meta <- ssGSEA_meta[order(ssGSEA_meta$SurvivalCutoff),]
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
        
        ## Determine type of survival data - OS/EFS/PFS?
        SurvDateType <- sub("\\..*","",surv_time_col)
        
        ## Adjust 'Sample Type' for label
        if (isTruthy(SampleType)) {
          SampleTypeLab <- paste(" (",SampleType,") ",sep = "")
        } else {
          SampleTypeLab <- " "
        }
        if (GeneSet %in% decon_score_cols) {
          scoreMethod <- "Pre-Processed Score"
        }
        if (input$GeneSetTabs == 2) {
          scoreMethod <- "Gene Expression Score"
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
          ggpubr::stat_compare_means(method = input$boxoptselecRisk,label.x = 0.5) +
          theme(text = element_text(size = font),
                legend.position = "none")
        plot
        
      })
      
      output$Sboxplot <- renderPlot({
        plot <- Sboxplot_react()
        plot
      })
      output$SboxplotTable <- renderDataTable({
        boxTab <- SboxplotReact()
        DT::datatable(boxTab,
                      options = list(paging = F),
                      rownames = F)
      })
      
      output$heatmap_error_message <- renderUI({
        geneset <- gs_react()
        geneset_name <- names(geneset)
        if (geneset_name %in% decon_score_cols() | input$GeneSetTabs == 2) {
          p("Heatmap not available for Pre-Processed scores or single gene analysis", style = "color:red")
        }
      })
      output$heatmap_error_message2 <- renderUI({
        geneset <- gs_react()
        geneset_name <- names(geneset)
        if (geneset_name %in% decon_score_cols() | input$GeneSetTabs == 2) {
          p("Heatmap not available for Pre-Processed scores or single gene analysis", style = "color:red")
        }
      })
      
      observe({
        Features <- c("SurvivalCutoff",metacol_feature(),metacol_survid(),metacol_survtime())
        updateSelectizeInput(session,"riskHeatAnno", choices = Features, selected = "SurvivalCutoff", server = T)
      })
      
      riskheat_colAnn <- reactive({
        if (length(input$riskHeatAnno) > 0) {
          meta <- SboxplotReact()
          rownames(meta) <- meta[,1]
          meta_sub <- meta[,input$riskHeatAnno, drop = F]
          meta_sub[,input$riskHeatAnno] <- as.data.frame(lapply(meta_sub[,input$riskHeatAnno, drop = F], factor))
          meta_sub <- meta_sub %>% relocate(any_of("SurvivalCutoff"), .after = last_col())
          rownames(meta_sub) <- rownames(meta)
          colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
                                                      name = input$riskHeatAnno,
                                                      which = 'col')
          colAnn
        }  else {
          colAnn <- NULL
          colAnn
        }
      })
      Sheatmap_react  <- reactive({
        
        req(gs_react())
        if (length(gs_react()) > 0) {
          geneset <- gs_react()
          GeneSet <- names(geneset)
          if (GeneSet %in% names(geneset)) {
            heatgenes <- geneset[[GeneSet]]
            meta <- SboxplotReact()
            expr_start <- exprSub()
            samples <- meta[,1]
            expr <- expr_start[which(rownames(expr_start) %in% heatgenes),colnames(expr_start) %in% samples, drop = F]
            clmethod <- input$ClusterMethod
            rowfont <- input$heatmapFontR
            colfont <- input$heatmapFontC
            colAnno <- riskheat_colAnn()
            
            dataset <- expr
            dataset <- log2(dataset + 1)
            zdataset <- apply(dataset, 1, scale)
            zdataset <- apply(zdataset, 1, rev)
            colnames(zdataset) <- names(dataset)
            dataset <- as.matrix(zdataset)
            dataset[is.na(dataset)] <- 0
            dataset = dataset[apply(dataset[,-1], 1, function(x) !all(x==0)),]
            
            #col_fun = colorRamp2(c(min(dataset), 0, max(dataset)), c("blue", "white", "red"))
            #lgd = Legend(col_fun = col_fun, title = "Expression")
            
            p <- suppressMessages(ComplexHeatmap::Heatmap(dataset,
                                                          top_annotation = colAnno,
                                                          clustering_method_rows = clmethod,
                                                          show_row_names = T, show_column_names = T,
                                                          #cluster_rows = clust_rows_opt,
                                                          cluster_columns = FALSE,
                                                          row_names_gp = gpar(fontsize = rowfont), column_names_gp = gpar(fontsize = colfont),
                                                          heatmap_legend_param = list(title = "Expression"),
                                                          border = F))
            draw(p, padding = unit(c(50, 50, 2, 2), "mm")) # unit(c(bottom,left,right,top))
            #draw(lgd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))
          }
        }
        
      })
      output$Sheatmap <- renderPlot({
        heat <- Sheatmap_react()
        heat
      })
      
      ## Feat Strat ------------------------------------------------------------
      
      observe({
        FeatureChoices <- c(metacol_feature(),metacol_survid())
        if (isTruthy(PreSelect_SecondaryFeature_react())) {
          if (!PreSelect_SecondaryFeature_react() %in% FeatureChoices) {
            selec <- FeatureChoices[1]
          } else { selec <- PreSelect_SecondaryFeature_react() }
        } else { selec <- FeatureChoices[1] }
        updateSelectizeInput(session,"BoxplotFeature",choices = c(metacol_feature(),metacol_survid()), selected = selec, server = T)
      })
      
      FeatureStrat_df <- reactive({
        req(input$BoxplotFeature)
        req(input$ScoreMethod)
        req(ssGSEAmeta())
        geneset <- gs_react()
        GeneSet <- names(geneset)
        FeatureSelec <- input$BoxplotFeature
        scoreMethod <- input$ScoreMethod
        meta_ssGSEA <- ssGSEAmeta()
        boxTab <- meta_ssGSEA[,c(colnames(meta_ssGSEA)[1],FeatureSelec,GeneSet)]
        if ("Remove NA/Unknowns" %in% input$FeatBoxOpts) {
          boxTab <- boxTab[which(is.na(boxTab[,FeatureSelec]) == FALSE),]
          boxTab <- boxTab[which(boxTab[,FeatureSelec] != "Inf"),]
          boxTab <- boxTab[grep("unknown",boxTab[,FeatureSelec],ignore.case = T, invert = T),]
        }
        boxTab[,FeatureSelec] <- as.factor(boxTab[,FeatureSelec])
        if ("Log Transform Score" %in% input$FeatBoxOpts) {
          boxTab[,GeneSet] <- log(boxTab[,GeneSet] + 1)
          scoreMethod <- paste(scoreMethod,"(log(score + 1))")
        }
        boxTab <- boxTab[order(boxTab[,FeatureSelec]),]
        boxTab
      })
      
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
        decon_score_cols <- decon_score_cols()
        boxTab <- FeatureStrat_df()
        
        ## Adjust 'Sample Type' for label
        if (isTruthy(SampleType)) {
          SampleTypeLab <- paste(" (",SampleType,") ",sep = "")
        } else {
          SampleTypeLab <- " "
        }
        if (GeneSet %in% decon_score_cols) {
          scoreMethod <- "Pre-Processed Score"
        }
        if (input$GeneSetTabs == 2) {
          scoreMethod <- "Gene Expression Score"
        }
        if (input$GeneSetTabs != 2 & !(GeneSet %in% decon_score_cols)) {
          scoreMethod <- paste(scoreMethod,"Score")
        }
        
        p <- ggplot(boxTab, aes(factor(boxTab[,FeatureSelec]), boxTab[,GeneSet], fill = boxTab[,FeatureSelec])) +
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
        
        if (boxplotang == 45) {
          p <- p + theme(text = element_text(size = font),
                         axis.text.x = element_text(angle = as.numeric(boxplotang), hjust = 1),
                         legend.position = "none")
        }
        if (boxplotang == 90) {
          p <- p + theme(text = element_text(size = font),
                         axis.text.x = element_text(angle = as.numeric(boxplotang)),
                         legend.position = "none")
        }
        p
        
      })
      ## Feature Boxplot
      output$Featureboxplot <- renderPlot({
        plot <- Featureboxplot_react()
        plot
      })
      
      ## Feature bloxplot table
      output$FeatureboxplotTable <- renderDataTable({
        boxTab <- FeatureStrat_df()
        DT::datatable(boxTab,
                      options = list(paging = F),
                      rownames = F)
      })
      
      
      observe({
        Features <- c(metacol_feature(),metacol_survid(),metacol_survtime())
        updateSelectizeInput(session,"stratHeatAnno", choices = Features, selected = input$BoxplotFeature, server = T)
      })
      
      stratheat_colAnn <- reactive({
        req(input$stratHeatAnno)
        req(input$BoxplotFeature)
        if (length(input$stratHeatAnno) > 0) {
          req(ssGSEAmeta())
          meta <- ssGSEAmeta()
          meta <- meta[order(meta[,input$BoxplotFeature]),]
          rownames(meta) <- meta[,1]
          meta_sub <- meta[,input$stratHeatAnno, drop = F]
          meta_sub[,input$stratHeatAnno] <- as.data.frame(lapply(meta_sub[,input$stratHeatAnno, drop = F], factor))
          meta_sub <- meta_sub %>% relocate(any_of(input$BoxplotFeature), .after = last_col())
          rownames(meta_sub) <- rownames(meta)
          colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
                                                      name = input$stratHeatAnno,
                                                      which = 'col')
          colAnn
        }  else {
          colAnn <- NULL
          colAnn
        }
      })
      FeatureHeatmap_react  <- reactive({
        
        req(gs_react())
        req(stratheat_colAnn())
        req(exprSub())
        req(FeatureStrat_df())
        req(input$ClusterMethod2)
        req(input$heatmapFontR)
        req(input$heatmapFontC)
        if (length(gs_react()) > 0) {
          geneset <- gs_react()
          GeneSet <- names(geneset)
          if (GeneSet %in% names(geneset)) {
            heatgenes <- geneset[[GeneSet]]
            meta <- FeatureStrat_df()
            expr_start <- exprSub()
            samples <- meta[,1]
            expr <- expr_start[which(rownames(expr_start) %in% heatgenes),colnames(expr_start) %in% samples, drop = F]
            clmethod <- input$ClusterMethod2
            rowfont <- input$heatmapFontR
            colfont <- input$heatmapFontC
            colAnno <- stratheat_colAnn()
            
            dataset <- expr
            dataset <- log2(dataset + 1)
            zdataset <- apply(dataset, 1, scale)
            zdataset <- apply(zdataset, 1, rev)
            colnames(zdataset) <- names(dataset)
            dataset <- as.matrix(zdataset)
            dataset[is.na(dataset)] <- 0
            dataset = dataset[apply(dataset[,-1], 1, function(x) !all(x==0)),]
            
            #col_fun = colorRamp2(c(min(dataset), 0, max(dataset)), c("blue", "white", "red"))
            #lgd = Legend(col_fun = col_fun, title = "Expression")
            
            p <- suppressMessages(ComplexHeatmap::Heatmap(dataset,
                                                          top_annotation = colAnno,
                                                          clustering_method_rows = clmethod,
                                                          show_row_names = T, show_column_names = T,
                                                          #cluster_rows = clust_rows_opt,
                                                          cluster_columns = FALSE,
                                                          row_names_gp = gpar(fontsize = rowfont), column_names_gp = gpar(fontsize = colfont),
                                                          heatmap_legend_param = list(title = "Expression"),
                                                          border = F))
            draw(p, padding = unit(c(50, 50, 2, 2), "mm")) # unit(c(bottom,left,right,top))
            #draw(lgd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))
          }
        }
        
      })
      output$FeatureHeatmap <- renderPlot({
        req(FeatureHeatmap_react())
        heat <- FeatureHeatmap_react()
        heat
      })
      
      
      ####----Text Output----####
      
      output$timewarnmessage1 <- renderUI({
        req(input$linPredict1)
        if (input$linPredict1 == "time") {
          p("Residual type must be schoenfeld or scaledsch.")
        }
        
      })
      output$timewarnmessage1S <- renderUI({
        req(input$linPredict1S)
        if (input$linPredict1S == "time") {
          p("Residual type must be schoenfeld or scaledsch.")
        }
        
      })
      
      output$timewarnmessage2 <- renderUI({
        req(input$linPredict2)
        if (input$linPredict2 == "time") {
          p("Residual type must be schoenfeld or scaledsch.")
        }
        
      })
      output$timewarnmessage2S <- renderUI({
        req(input$linPredict2S)
        if (input$linPredict2S == "time") {
          p("Residual type must be schoenfeld or scaledsch.")
        }
        
      })
      
      output$timewarnmessage3 <- renderUI({
        req(input$linPredict3)
        if (input$linPredict3 == "time") {
          p("Residual type must be schoenfeld or scaledsch.")
        }
        
      })
      output$timewarnmessage3S <- renderUI({
        req(input$linPredict3S)
        if (input$linPredict3S == "time") {
          p("Residual type must be schoenfeld or scaledsch.")
        }
        
      })
      
      
      # Downloaders ------------------------------------------------------------
      
      ## Path Surv Plots -------------------------------------------------------
      #observe({
      #  req(SplotBIN_react())
      #  plot <- SplotBIN_react()$plot
      #  dnldPlot_server("dnldSplot_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_QuaretileCutPoint_Survival_",Sys.Date(),".svg")),
      #                  input$PlotDnldHight,input$PlotDnldWidth,input$PlotDnldUnits)
      #})
      #observe({
      #  req(ScutPointPlot_react())
      #  plot <- ScutPointPlot_react()$plot
      #  dnldPlot_server("dnldScutPointPlot_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_OptimalCutPoint_Survival_",Sys.Date(),".svg")),
      #                  input$PlotDnldHight,input$PlotDnldWidth,input$PlotDnldUnits)
      #})
      #observe({
      #  req(SquantPlot_react())
      #  plot <- SquantPlot_react()$plot
      #  dnldPlot_server("dnldSquantPlot_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_TopBottomCutPoint_Survival_",Sys.Date(),".svg")),
      #                  input$PlotDnldHight,input$PlotDnldWidth,input$PlotDnldUnits)
      #})
      #observe({
      #  req(SquantPlot2_react())
      #  plot <- SquantPlot2_react()$plot
      #  dnldPlot_server("dnldSquantPlot2_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_UserCutPoint_Survival_",Sys.Date(),".svg")),
      #                  input$PlotDnldHight,input$PlotDnldWidth,input$PlotDnldUnits)
      #})
      #
      ### Path Density Plots ----------------------------------------------------
      #
      #observe({
      #  req(ssgseaBINDensity_react())
      #  plot <- ssgseaBINDensity_react()
      #  dnldPlot_server("dnldssgseaBINDensity_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_MedianCutPoint_Density_",Sys.Date(),".svg")),
      #                  input$PlotDnldHight,input$PlotDnldWidth,input$PlotDnldUnits)
      #})
      #observe({
      #  req(ssgseaQuartDensity_react())
      #  plot <- ssgseaQuartDensity_react()
      #  dnldPlot_server("dnldssgseaQuartDensity_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_QuaretileCutPoint_Density_",Sys.Date(),".svg")),
      #                  input$PlotDnldHight,input$PlotDnldWidth,input$PlotDnldUnits)
      #})
      #observe({
      #  req(ssgseaCutPDensity_react())
      #  plot <- ssgseaCutPDensity_react()
      #  dnldPlot_server("dnldssgseaCutPDensity_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_OptimalCutPoint_Density_",Sys.Date(),".svg")),
      #                  input$PlotDnldHight,input$PlotDnldWidth,input$PlotDnldUnits)
      #})
      #observe({
      #  req(ssgseaQuantDensity_react())
      #  plot <- ssgseaQuantDensity_react()
      #  dnldPlot_server("dnldssgseaQuantDensity_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_TopBottomCutPoint_Density_",Sys.Date(),".svg")),
      #                  input$PlotDnldHight,input$PlotDnldWidth,input$PlotDnldUnits)
      #})
      #observe({
      #  req(ssgseaQuant2Density_react())
      #  plot <- ssgseaQuant2Density_react()
      #  dnldPlot_server("dnldssgseaQuant2Density_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_UserCutPoint_Density_",Sys.Date(),".svg")),
      #                  input$PlotDnldHight,input$PlotDnldWidth,input$PlotDnldUnits)
      #})
      #
      ### Univariate Plots ------------------------------------------------------
      #
      #observe({
      #  req(featSplot_react())
      #  plot <- featSplot_react()$plot
      #  dnldPlot_server("dnldfeatSplot_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_Univariate_Survival_",Sys.Date(),".svg")),
      #                  input$PlotDnldHight,input$PlotDnldWidth,input$PlotDnldUnits)
      #})
      ##observe({
      ##  req(SinglevarForestPlot_react())
      ##  plot <- SinglevarForestPlot_react()
      ##  dnldPlot_server("dnldUniVarForestplot_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_Univariate_Forest_",Sys.Date(),".svg")),
      ##                  input$PlotDnldHight,input$PlotDnldWidth,input$PlotDnldUnits)
      ##})
      #output$dnldUniVarForestplot_SVG <- shiny::downloadHandler(
      #  filename = function() {
      #    gsub("[[:space:]]","",paste0(ProjectName_react(),"_Univariate_Forest_",Sys.Date(),".svg"))
      #  },
      #  content = function(file) {
      #    ggplot2::ggsave(filename = file, plot = SinglevarForestPlot_react(),
      #                    width = input$PlotDnldWidth, height = input$PlotDnldHight,units = input$PlotDnldUnits)
      #  })
      #observe({
      #  req(MultiFeatUnivarForestPlot_react())
      #  plot <- MultiFeatUnivarForestPlot_react()
      #  dnldPlot_server("dnldMultiFeatUnivarForestPlot_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_Univariate_MultiFeatureForest_",Sys.Date(),".svg")),
      #                  type = "forest")
      #})
      #observe({
      #  req(MultiFeatUnivarForestPlotTab_react())
      #  plot <- MultiFeatUnivarForestPlotTab_react()[,-7]
      #  dnldDF_server("dnldMultiFeatUnivarForestPlot_table",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_Univariate_MultiFeatureForest_",Sys.Date(),".txt")))
      #})
      #observe({
      #  req(MultiFeat_ForestMeta())
      #  plot <- MultiFeat_ForestMeta()
      #  dnldDF_server("dnldunivarForestPlotTable",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_Univariate_MultiFeatureForestData_",Sys.Date(),".txt")))
      #})
      #observe({
      #  req(UnivarLinearityPlot_react())
      #  plot <- UnivarLinearityPlot_react()
      #  dnldPlot_server("dnldUniVarLinplot_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_Univariate_Linearity_",Sys.Date(),".svg")),
      #                  input$PlotDnldHight,input$PlotDnldWidth,input$PlotDnldUnits)
      #})
      #
      ### Bivariate Add Plots ------------------------------------------------------
      #
      ##observe({
      ##  req(BivarForestPlot_react())
      ##  plot <- BivarForestPlot_react()
      ##  dnldPlot_server("dnldBiVarAddForest_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_BivariateAdditive_Forest_",Sys.Date(),".svg")),
      ##                  input$PlotDnldHight,input$PlotDnldWidth,input$PlotDnldUnits)
      ##})
      #output$dnldBiVarAddForest_SVG <- shiny::downloadHandler(
      #  filename = function() {
      #    gsub("[[:space:]]","",paste0(ProjectName_react(),"_BivariateAdditive_Forest_",Sys.Date(),".svg"))
      #  },
      #  content = function(file) {
      #    ggplot2::ggsave(filename = file, plot = BivarForestPlot_react(),
      #                    width = input$PlotDnldWidth, height = input$PlotDnldHight,units = input$PlotDnldUnits)
      #  })
      #observe({
      #  req(BivarLinearityPlot_react())
      #  plot <- BivarLinearityPlot_react()
      #  dnldPlot_server("dnldBiVarAddLinplot_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_BivariateAdditive_Linearity_",Sys.Date(),".svg")),
      #                  input$PlotDnldHight,input$PlotDnldWidth,input$PlotDnldUnits)
      #})
      #
      ### Bivariate Int Plots ------------------------------------------------------
      #
      #observe({
      #  req(featSplotBi_react())
      #  plot <- featSplotBi_react()$plot
      #  dnldPlot_server("dnldfeatSplotBi_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_BivariateInteractive_Survival_",Sys.Date(),".svg")),
      #                  input$PlotDnldHight,input$PlotDnldWidth,input$PlotDnldUnits)
      #})
      #observe({
      #  req(BivarLinearityPlotInter_react())
      #  plot <- BivarLinearityPlotInter_react()
      #  dnldPlot_server("dnldBiVarIntLinplot_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_BivariateInteractive_Linearity_",Sys.Date(),".svg")),
      #                  input$PlotDnldHight,input$PlotDnldWidth,input$PlotDnldUnits)
      #})
      #
      ### Multivariate Plots ------------------------------------------------------
      #
      ##observe({
      ##  req(MultivarForestPlot_react())
      ##  plot <- MultivarForestPlot_react()
      ##  dnldPlot_server("dnldMultiVarForest_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_Multivariate_Forest_",Sys.Date(),".svg")),
      ##                  input$PlotDnldHight,input$PlotDnldWidth,input$PlotDnldUnits)
      ##})
      #
      #output$dnldMultiVarForest_SVG <- shiny::downloadHandler(
      #  filename = function() {
      #    gsub("[[:space:]]","",paste0(ProjectName_react(),"_Multivariate_Forest_",Sys.Date(),".svg"))
      #  },
      #  content = function(file) {
      #    ggplot2::ggsave(filename = file, plot = MultivarForestPlot_react(),
      #                    width = input$PlotDnldWidth, height = input$PlotDnldHight,units = input$PlotDnldUnits)
      #  })
      #observe({
      #  req(MultiFeatMultivarForestPlot_react())
      #  plot <- MultiFeatMultivarForestPlot_react()
      #  dnldPlot_server("dnldMultiFeatMultivarForestPlot_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_Multivariate_MultiFeatureForest_",Sys.Date(),".svg")),
      #                  type = "forest")
      #})
      #observe({
      #  req(MultiFeatMultivarForestPlotTab_react())
      #  plot <- MultiFeatMultivarForestPlotTab_react()[,-7]
      #  dnldDF_server("dnldMultiFeatMultivarForestPlot_table",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_Multivariate_MultiFeatureForest_",Sys.Date(),".txt")))
      #})
      #observe({
      #  req(MultiFeat_InterForestMeta())
      #  plot <- MultiFeat_InterForestMeta()
      #  dnldDF_server("dnldmultiForestPlotTable",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_Multivariate_MultiFeatureForestData_",Sys.Date(),".txt")))
      #})
      #
      ### Data Exploration Plots ------------------------------------------------------
      #
      #observe({
      #  req(ssGSEAmeta())
      #  plot <- ssGSEAmeta()
      #  dnldDF_server("DnldClin",plot,
      #                gsub("[[:space:]]","",paste0(ProjectName_react(),"_Clinical_",Sys.Date(),".txt")))
      #})
      #exprDnld <- reactive({
      #  req(exprSub())
      #  expr <- as.data.frame(exprSub())
      #  expr$Gene <- rownames(expr)
      #  expr <- expr %>%
      #    relocate(Gene)
      #  expr
      #})
      #observe({
      #  req(exprDnld())
      #  plot <- exprDnld()
      #  dnldDF_server("DnldExpr",plot,
      #                gsub("[[:space:]]","",paste0(ProjectName_react(),"_Expression_",Sys.Date(),".txt")))
      #})
      #observe({
      #  req(ssgseaDensity_react())
      #  plot <- ssgseaDensity_react()
      #  dnldPlot_server("dnldssgseaDensity_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_",names(gs_react()),"_Density_",Sys.Date(),".svg")),
      #                  input$PlotDnldHight,input$PlotDnldWidth,input$PlotDnldUnits)
      #})
      #observe({
      #  req(ssGSEAmeta())
      #  plot <- ssGSEAmeta()[,c(colnames(ssGSEAmeta())[1],names(gs_react()))]
      #  dnldDF_server("dnldssgseaDensityTable",plot,
      #                gsub("[[:space:]]","",paste0(ProjectName_react(),"_",names(gs_react()),"_Density_",Sys.Date(),".txt")))
      #})
      #observe({
      #  req(FeatCompScatter_react())
      #  plot <- FeatCompScatter_react()
      #  dnldDF_server("dnldFeatCompScatterTable",plot,
      #                gsub("[[:space:]]","",paste0(ProjectName_react(),"_ScatterComparisonTable_",Sys.Date(),".txt")))
      #})
      #observe({
      #  req(Sboxplot_react())
      #  plot <- Sboxplot_react()
      #  dnldPlot_server("dnldSboxplot_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_RiskStrat_Boxplot_",Sys.Date(),".svg")),
      #                  input$PlotDnldHight,input$PlotDnldWidth,input$PlotDnldUnits)
      #})
      #observe({
      #  req(SboxplotReact())
      #  req(input$SurvivalType_id)
      #  req(input$SurvivalType_time)
      #  plot <- as.data.frame(SboxplotReact()[,c(colnames(SboxplotReact())[1],input$SurvivalType_time,input$SurvivalType_id,names(gs_react()),"SurvivalCutoff")])
      #  dnldDF_server("dnldSBoxplotTab",plot,
      #                gsub("[[:space:]]","",paste0(ProjectName_react(),"_RiskStrat_BoxplotTable_",Sys.Date(),".txt")))
      #})
      #
      #dnldHeatExpr <- reactive({
      #  req(exprSub())
      #  req(gs_react())
      #  expr <- exprSub()
      #  GSgenes <- unname(unlist(gs_react()))
      #  if (length(GSgenes) > 0) {
      #    expr <- as.data.frame(expr[which(rownames(expr) %in% GSgenes),])
      #    expr$Gene <- rownames(expr)
      #    expr <- expr %>%
      #      relocate(Gene)
      #    expr
      #  }
      #})
      #
      #observe({
      #  req(Sheatmap_react())
      #  plot <- Sheatmap_react()
      #  dnldPlot_server("dnldSheatmap_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_RiskStrat_Heatmap_",Sys.Date(),".svg")),
      #                  input$PlotDnldHight,input$PlotDnldWidth,type = "complex")
      #})
      #observe({
      #  req(dnldHeatExpr())
      #  plot <- dnldHeatExpr()
      #  dnldDF_server("dnldSheatmapexpr",plot,
      #                gsub("[[:space:]]","",paste0(ProjectName_react(),"_",names(gs_react()),"_RiskStrat_Heatmap_Expression",Sys.Date(),".txt")))
      #})
      #observe({
      #  req(Featureboxplot_react())
      #  plot <- Featureboxplot_react()
      #  dnldPlot_server("dnldFboxplot_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_FeatureStrat_Boxplot_",Sys.Date(),".svg")),
      #                  input$PlotDnldHight,input$PlotDnldWidth,input$PlotDnldUnits)
      #})
      #observe({
      #  req(ssGSEAmeta())
      #  req(input$BoxplotFeature)
      #  plot <- as.data.frame(ssGSEAmeta()[,c(colnames(ssGSEAmeta())[1],input$BoxplotFeature,names(gs_react()))])
      #  dnldDF_server("dnldFeatureboxplotTab",plot,
      #                gsub("[[:space:]]","",paste0(ProjectName_react(),"_FeatureStrat_BoxplotTable_",Sys.Date(),".txt")))
      #})
      #observe({
      #  req(FeatureHeatmap_react())
      #  plot <- FeatureHeatmap_react()
      #  dnldPlot_server("dnldFheatmap_SVG",plot,gsub("[[:space:]]","",paste0(ProjectName_react(),"_FeatureStrat_Heatmap_",Sys.Date(),".svg")),
      #                  input$PlotDnldHight,input$PlotDnldWidth,type = "complex")
      #})
      #observe({
      #  req(dnldHeatExpr())
      #  plot <- dnldHeatExpr()
      #  dnldDF_server("dnldFheatmapexpr",plot,
      #                gsub("[[:space:]]","",paste0(ProjectName_react(),"_",names(gs_react()),"_FeatureStrat_Heatmap_Expression",Sys.Date(),".txt")))
      #})
      
      
      
    }
    
  }
  
  )
  
}


shinyApp(ui = ui, server = server)

