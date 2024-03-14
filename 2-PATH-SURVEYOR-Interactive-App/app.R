
# User Input -------------------------------------------------------------------

ProjectName <- 'PAN ICI Melanoma - Van Allen anti-CTLA4'

ExpressionMatrix_file <- 'Example_Data/Expression_Data.zip'

MetaData_file <- 'Example_Data/Clinical_Data.txt'

MetaParam_File <- 'Example_Data/Clinical_Parameters.txt'


## Password Settings -----------------------------------------------------------
Password_Protected <- FALSE
PasswordSet <- ""

## Pre-Selected Inputs ---------------------------------------------------------
# An option from the meta, All, or NULL
PreSelect_SamplyType <- "All"
PreSelect_Feature <- "All"
# An option from the meta or NULL
PreSelect_SubFeature <- NULL
PreSelect_SecondaryFeature <- "Responder"

## Provided Input --------------------------------------------------------------
## User make sure paths are correct
GeneSet_File <- "GeneSet_Data/GeneSet_List.RData"
GeneSetTable_File <- "GeneSet_Data/GeneSet_CatTable.zip"
About_MD_File <- "App_Markdowns/PurposeAndMethods.Rmd"
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

# Load Geneset Data
# R Data list load function for naming
#loadRData <- function(fileName){
#  #loads an RData file, and returns it
#  load(fileName)
#  get(ls()[ls() != "fileName"])
#}
gs <- loadRData(GeneSet_File)
# Gene Set Table
#GeneSetTable <- as.data.frame(read_delim(GeneSetTable_File, delim = '\t', col_names = T))
GeneSetTable <- as.data.frame(fread(GeneSetTable_File, sep = '\t', header = T))
GeneSetCats <- unique(GeneSetTable[,1])

# Functions --------------------------------------------------------------------

#increase file upload size
options(shiny.maxRequestSize=5000*1024^2)

SurvPlot_Height <- "550px"
SurvPlot_Width <- "850px"

StatCols <- c("MedianCutP","QuartileCutP","OptimalCutP","TopBottomCutP","UserCutP")

## Quartile Conversion
#quartile_conversion = function(mat) {
#  new_mat = mat;
#  new_mat[mat <=  quantile(as.numeric(mat), na.rm = T)[2]] = "Q1_Low";
#  new_mat[mat > quantile(as.numeric(mat), na.rm = T)[2] & mat <= quantile(mat)[3]] = "Q2_MedLow";
#  new_mat[mat > quantile(as.numeric(mat), na.rm = T)[3] & mat <= quantile(mat)[4]] = "Q3_MedHigh";
#  new_mat[mat > quantile(as.numeric(mat), na.rm = T)[4]] = "Q4_High";
#  return (new_mat)
#}
#
### High-Low
#highlow = function(mat) {
#  new_mat = mat;
#  new_mat[mat > quantile(as.numeric(mat), na.rm = T)[3]] = "High";
#  new_mat[mat <= quantile(as.numeric(mat), na.rm = T)[3]] = "Low";
#  return (new_mat)
#}
#
#highlow2 = function(mat) {
#  new_mat = mat;
#  new_mat[mat > quantile(as.numeric(mat), na.rm = T)[3]] = "Above_Median";
#  new_mat[mat <= quantile(as.numeric(mat), na.rm = T)[3]] = "Below_Median";
#  return (new_mat)
#}
#
#quantile_conversion = function(mat,cutoff) {
#  new_mat = mat;
#  new_mat[mat >= quantile(as.numeric(mat),1-cutoff, na.rm = T)] = "High";
#  new_mat[mat <= quantile(as.numeric(mat),cutoff, na.rm = T)] = "Low";
#  new_mat[mat > quantile(as.numeric(mat),cutoff, na.rm = T) & mat < quantile(mat,1-cutoff, na.rm = T)] = "BetweenCutoff";
#  return (new_mat)
#}
#
#quantile_conversion2 = function(mat,cutoff) {
#  new_mat = mat;
#  new_mat[mat > quantile(as.numeric(mat),cutoff, na.rm = T)] = "High";
#  new_mat[mat <= quantile(as.numeric(mat),cutoff, na.rm = T)] = "Low";
#  return (new_mat)
#}
#
#add_x_intercepts <- function(p) {
#  
#  p2 <- ggplot_build(p)
#  breaks <- p2$layout$panel_params[[1]]$x$breaks
#  breaks <- breaks[!is.na(breaks)]
#  
#  vals <- unlist(lapply(seq_along(p$layers), function(x) {
#    d <- layer_data(p, x)
#    if('xintercept' %in% names(d)) d$xintercept else numeric()
#  }))
#  
#  p + scale_x_continuous(breaks = sort(c(vals, breaks)))
#}
#
#get_lik_pval <- function(x_tab) {
#  tab <- x_tab
#  out <- capture.output(summary(tab))
#  lik_line <- grep("^Likelihood ratio test=",out,value = T)
#  lik_line_P <- as.numeric(str_split(str_split(lik_line,", ")[[1]][2],"=")[[1]][2])
#  return(lik_line_P)
#}
#
#gsubCheck <- function(string) {
#  string2 <- gsub("\\+","POS",string)
#  string3 <- gsub("[[:punct:]]","_",string2)
#  string3 <- gsub(" ","_",string3)
#  return(string3)
#}
#
#lm_eqn <- function(df){
#  m <- lm(y ~ x, df);
#  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
#                   list(a = format(unname(coef(m)[1]), digits = 2),
#                        b = format(unname(coef(m)[2]), digits = 2),
#                        r2 = format(summary(m)$r.squared, digits = 3)))
#  as.character(as.expression(eq));
#}
#
#get_tabData <- function(tab) {
#  ## Developed on package "suvival v3.4-0 and survminder v0.4.9
#  Variable <- as.character(tab[["terms"]][[3]])
#  N <- tab[["n"]]
#  out <- capture.output(summary(tab))
#  coef_line1 <- strsplit(grep(paste0("^",Variable),out, value = T)[1], "\\s+")[[1]]
#  coef_line2 <- strsplit(grep(paste0("^",Variable),out, value = T)[2], "\\s+")[[1]]
#  lik_line <- grep("^Likelihood ratio test=",out,value = T)
#  lik_line_P <- as.numeric(str_split(str_split(lik_line,", ")[[1]][2],"=")[[1]][2])
#  tabData <- c(Variable = Variable, N = as.numeric(N), `Hazard Ratio` = as.numeric(coef_line1[3]), `Standard Error` = as.numeric(coef_line1[4]),
#               Low = as.numeric(coef_line2[4]), High = as.numeric(coef_line2[5]), `P.Value` = as.numeric(lik_line_P))
#  return(tabData)
#}
#
#GetColsOfType <- function(meta,type = c("discrete","continuous"), threshold = 0.75) {
#  require(dplyr)
#  MetaClass_num <- meta %>%
#    type.convert(as.is = TRUE) %>% 
#    dplyr::select(where(is.numeric)) %>%
#    names()
#  IntMetaCols <- apply(meta[,MetaClass_num, drop = F],2,function(x) any(round(as.numeric(x)) != as.numeric(x)))
#  IntMetaCols <- names(IntMetaCols)[which(IntMetaCols == T)]
#  MetaClass_NonNum <- colnames(meta)[which(!colnames(meta) %in% MetaClass_num)]
#  ShowOrNot <- apply(meta[,MetaClass_num,drop = F],2,function(x) any(length(levels(as.factor(x)))<(nrow(meta)*threshold)))
#  discreteCols <- c(MetaClass_NonNum,names(ShowOrNot)[which(ShowOrNot == T)])
#  discreteCols <- discreteCols[which(!discreteCols %in% IntMetaCols)]
#  continuousCols <- colnames(meta)[which(!colnames(meta) %in% discreteCols)]
#  if (toupper(type) == "DISCRETE") {
#    return(discreteCols)
#  } else if (toupper(type) == "CONTINUOUS") {
#    return(continuousCols)
#  } else {
#    print("ERROR: Argument 'type' must be 'discrete' or 'continuous'")
#  }
#}



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
if (!FileProvided) {
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
                                ),
                              ),
                              mainPanel(
                                uiOutput("rendExprFilePrevHeader"),
                                div(DT::dataTableOutput("ExprFile_Preview"), style = "font-size:10px"),
                                uiOutput("rendClinFilePrevHeader"),
                                div(DT::dataTableOutput("ClinFile_Preview"), style = "font-size:10px"),
                                uiOutput("rendParamFilePrevHeader"),
                                div(DT::dataTableOutput("ClinParamFile_Preview"), style = "font-size:10px")
                              )
                            )
  )
}

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
                                                                    numericInput("boxplotFont","Boxplot Font Size:", value = 15, step = 1)
                                                                     
                                                             ),
                                                             column(6,
                                                                    numericInput("boxplotDot", "Boxplot Dot Size:", value = 0.75, step = 0.25),
                                                                    selectInput("boxplotTextAngle","X-Axis Text Orientation",
                                                                                choices = c("Horizontal (0 degrees)" = 0,"Angled (45 degrees)" = 45,"Vertical (90 degrees)" = 90))
                                                             )
                                                           ),
                                                           shiny::hr(),
                                                           h4("Heatmap Parameters"),
                                                           fluidRow(
                                                             column(6,
                                                                    numericInput("heatmapFontR", "Heatmap Row Font Size:", value = 9, step = 1)
                                                             ),
                                                             column(6,
                                                                    numericInput("heatmapFontC", "Heatmap Column Font Size:", value = 10, step = 1)
                                                             )
                                                           ),
                                                           #selectInput("ColorPaletteHeat", "Select Color Palette:",
                                                           #            choices = c("Red/Blue" = "original",
                                                           #                        "OmniBlueRed" = "OmniBlueRed",
                                                           #                        "LightBlue/BlackRed" = "LightBlueBlackRed",
                                                           #                        "Green/Black/Red" = "GreenBlackRed",
                                                           #                        "Yellow/Green/Blue" = "YlGnBu","Inferno" = "Inferno",
                                                           #                        "Viridis" = "Viridis","Plasma" = "Plasma",
                                                           #                        "Reds" = "OrRd","Blues" = "PuBu","Greens" = "Greens")
                                                           #),
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
                               )
                               #conditionalPanel(condition = "input.SurvPanels == '5'",
                               #                 tabsetPanel(
                               #                   tabPanel("Input Parameters",
                               #                            p(),
                               #                            h4("Sample Selection"),
                               #                            uiOutput("rendSampleTypeSelection_lasso"),
                               #                            uiOutput("rendFeatureSelection_lasso"),
                               #                            uiOutput("rendSubFeatureSelection_lasso"),
                               #                            h4("Feature Selection"),
                               #                            textInput("LassoModelName","Lasso Model Name:",value = "Custom_Lasso_Model"),
                               #                            uiOutput("rendLassoFeatureSelection_lasso"),
                               #                            h4("Lasso Parameters"),
                               #                            fluidRow(
                               #                              column(6,
                               #                                     uiOutput("rendSurvTimeSelec_lasso")
                               #                              ),
                               #                              column(6,
                               #                                     uiOutput("rendSurvIDSelect_lasso")
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
                               #                                     downloadButton("dnldLassoModel","Download Lasso Model")
                               #                              ),
                               #                              column(6,
                               #                                     downloadButton("dnldLassoRunData","Download Lasso Run Data")
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
                                                              div(shinycssloaders::withSpinner(tableOutput("SQuantileHR2tab"), type = 7, size = 0.5), style = "font-size:12px"),
                                                              #uiOutput("rendQuantHRtab2"),
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
                                                   selectizeInput("SingleSurvivalFeature","Select Feature:",choices = NULL, selected = 1),
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
                                                     fluidRow(
                                                       downloadButton("dnldfeatSplot_SVG","Download as SVG"),
                                                       downloadButton("dnldfeatSplot_PDF","Download as PDF")
                                                     )
                                                     
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
                                                     fluidRow(
                                                       downloadButton("dnldUniVarForestplot_SVG","Download as SVG"),
                                                       downloadButton("dnldUniVarForestplot_PDF","Download as PDF")
                                                     )
                                            ),
                                            tabPanel("Multi-Feature Forest Plot",
                                                     p(),
                                                     fluidRow(
                                                       column(6,
                                                              selectizeInput("MultiFeatUnivarSelect","Select Features for Forest Plot",
                                                                             choices = NULL, multiple = T, selected = 1, width = "100%")
                                                              #uiOutput("rendMultiFeatUnivarSelect")
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
                                                       downloadButton("dnldMultiFeatUnivarForestPlot_table","Download as Table"),
                                                       downloadButton("dnldMultiFeatUnivarForestPlot_SVG","Download as SVG")
                                                     ),
                                                     p(),
                                                     div(DT::dataTableOutput("univarForestPlotTable"), style = "font-size:12px"),
                                                     downloadButton("dnldunivarForestPlotTable","Download Table")
                                                     
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
                                                     fluidRow(
                                                       downloadButton("dnldUniVarLinplot_SVG","Download as SVG"),
                                                       downloadButton("dnldUniVarLinplot_PDF`","Download as PDF")
                                                     )
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
                                                              selectizeInput("SurvivalFeatureBi1","Select Feature 1:",choices = NULL, selected = 1),
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
                                                              selectizeInput("SurvivalFeatureBi2","Select Feature 2:",choices = NULL, selected = 1),
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
                                                                fluidRow(
                                                                  downloadButton("dnldBiVarAddForest_SVG","Download as SVG"),
                                                                  downloadButton("dnldBiVarAddForest_PDF","Download as PDF")
                                                                )
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
                                                                fluidRow(
                                                                  downloadButton("dnldBiVarAddLinplot_SVG","Download as SVG"),
                                                                  downloadButton("dnldBiVarAddLinplot_PDF","Download as PDF")
                                                                )
                                                       )
                                                     )
                                            ),
                                            
                                            ##### Multivar Int Surv ----------------------------
                                            
                                            tabPanel("Bivariate Interaction Survival Analysis",
                                                     p(),
                                                     fluidRow(
                                                       column(6,
                                                              selectizeInput("SurvivalFeatureBi1Inter","Select Feature 1:",choices = NULL, selected = 1),
                                                              #uiOutput("rendSurvivalFeatureBi1Inter"),
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
                                                                       #uiOutput("rendBiVarIntContHiLoCheck1")
                                                                )
                                                              ),
                                                              uiOutput("rendSurvFeatVariableBi1Inter")
                                                              
                                                       ),
                                                       column(6,
                                                              selectizeInput("SurvivalFeatureBi2Inter","Select Feature 1:",choices = NULL, selected = 1),
                                                              #uiOutput("rendSurvivalFeatureBi2Inter"),
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
                                                                       #uiOutput("rendBiVarIntContHiLoCheck2")
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
                                                                fluidRow(
                                                                  downloadButton("dnldfeatSplotBi_SVG","Download as SVG"),
                                                                  downloadButton("dnldfeatSplotBi_PDF","Download as PDF")
                                                                )
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
                                                                fluidRow(
                                                                  downloadButton("dnldBiVarIntLinplot_SVG","Download as SVG"),
                                                                  downloadButton("dnldBiVarIntLinplot_PDF","Download as PDF")
                                                                )
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
                                                                fluidRow(
                                                                  downloadButton("dnldMultiVarForest_SVG","Download as SVG"),
                                                                  downloadButton("dnldMultiVarForest_PDF","Download as PDF")
                                                                )
                                                       )
                                                     )
                                            ),
                                            ##### Multivar Forest -----------------
                                            tabPanel("Multivariate Forest Plot",
                                                     p(),
                                                     fluidRow(
                                                       column(3,
                                                              selectizeInput("MultiFeatMultivarSelect","Select Main Forest Feature",
                                                                             choices = NULL, multiple = F, selected = 1, width = "100%"),
                                                              #uiOutput("rendMultiFeatMultivarSelect"),
                                                              uiOutput("rendMultiFeatMultivarRefSelect")
                                                       ),
                                                       column(4,
                                                              selectizeInput("MultiFeatMultivarSubSelect","Select and Add Features to Forest Plot:",
                                                                             choices = NULL, multiple = T, selected = 1, width = "100%")
                                                              #uiOutput("rendMultiFeatMultivarSubSelect")
                                                       ),
                                                       column(2,
                                                              textInput("MultivarForestPlotXlim","X-Limits (Low,High)",value = "0,5",
                                                                        placeholder = "Low,High") #, width = "150px"
                                                       ),
                                                       column(2,
                                                              selectInput("MultivarForestPlotXtrans","X-Axis Transform",
                                                                          choices = c("none", "log", "log2", "log10"))
                                                       )
                                                     ),
                                                     uiOutput("rendMultiFeatMultivarForestPlot"),
                                                     fluidRow(
                                                       downloadButton("dnldMultiFeatMultivarForestPlot_table","Download as Table"),
                                                       downloadButton("dnldMultiFeatMultivarForestPlot_SVG","Download as SVG")
                                                     ),
                                                     p(),
                                                     div(DT::dataTableOutput("multiForestPlotTable"), style = "font-size:12px"),
                                                     downloadButton("dnldmultiForestPlotTable","Download Table")
                                                     
                                            )
                                          ),
                                          value = 3),
                                 
                                 #### Lasso ----------------------------       
                                 
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
                                 #                  downloadButton("dnldSplotLassoTrain_SVG","Download as SVG"),
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
                                 #                  downloadButton("dnldSplotLassoTest_SVG","Download as SVG"),
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
                                                       column(6,
                                                              uiOutput("DnldMetaButon")
                                                       ),
                                                       column(6,
                                                              uiOutput("DnldExprButon")
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
                                                     fluidRow(
                                                       downloadButton("dnldssgseaDensity_SVG","Download as SVG"),
                                                       downloadButton("dnldssgseaDensity_PDF","Download as PDF")
                                                     ),
                                                     div(DT::dataTableOutput("ssgseaDensityTable"), style = "font-size:12px"),
                                                     downloadButton("dnldssgseaDensityTable","Download Table"),
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
                                                     fluidRow(
                                                       downloadButton("dnldFeatCompScatter_SVG","Download as SVG"),
                                                       downloadButton("dnldFeatCompScatter_PDF","Download as PDF")
                                                     ),
                                                     p(),
                                                     div(DT::dataTableOutput("FeatCompScatterTable"), style = "font-size:12px"),
                                                     downloadButton("dnldFeatCompScatterTable","Download Table"),
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
                                                                fluidRow(
                                                                  downloadButton("dnldSboxplot_SVG","Download as SVG"),
                                                                  downloadButton("dnldSboxplot_PDF","Download as PDF")
                                                                ),
                                                                div(DT::dataTableOutput("SboxplotTable"), style = "font-size:12px; height:450px; overflow-Y: scroll"),
                                                                p(),
                                                                downloadButton("dnldSBoxplotTab","Download Table"),
                                                                value = "box"
                                                       ),
                                                       tabPanel("Risk Straification Heatmap",
                                                                p(),
                                                                uiOutput("heatmap_error_message"),
                                                                shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("Sheatmap", width = "100%", height = "2000px")), type = 6),
                                                                fluidRow(
                                                                  downloadButton("dnldSheatmap_SVG","Download as SVG"),
                                                                  downloadButton("dnldSheatmap_PDF","Download as PDF"),
                                                                  downloadButton("dnldSheatmapexpr","Download Expression Matrix From Heatmap")
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
                                                                fluidRow(
                                                                  downloadButton("dnldFboxplot_SVG","Download as SVG"),
                                                                  downloadButton("dnldFboxplot_PDF","Download as PDF")
                                                                ),
                                                                div(DT::dataTableOutput("FeatureboxplotTable"), style = "font-size:12px; height:450px; overflow-Y: scroll"),
                                                                p(),
                                                                downloadButton("dnldFeatureboxplotTab","Download Table"),
                                                                value = "box"
                                                       ),
                                                       tabPanel("Feature Heatmap",
                                                                p(),
                                                                uiOutput("heatmap_error_message2"),
                                                                shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput("FeatureHeatmap", width = "100%", height = "2000px")), type = 6),
                                                                fluidRow(
                                                                  downloadButton("dnldFheatmap_SVG","Download as SVG"),
                                                                  downloadButton("dnldFheatmap_PDF","Download as PDF"),
                                                                  downloadButton("dnldFheatmapexpr","Download Expression Matrix From Heatmap")
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
  if (FileProvided) {
    ui <- navbarPage(paste(ProjectName),
                     id = "tabs",
                     collapsible = TRUE,
                     Survival_tab,
                     About_tab)
  } else {
    ui <- navbarPage(paste(ProjectName),
                     id = "tabs",
                     collapsible = TRUE,
                     DataInput_tab,
                     Survival_tab,
                     About_tab)
  }
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
          appendTab("tabs", Survival_tab, select = TRUE)
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
                #metaP <- as.data.frame(read_delim(MetaParam_File_react(), delim = '\t', col_names = F))
                metaP <- as.data.frame(fread(MetaParam_File_react(), sep = '\t', header = F))
              } else if (tools::file_ext(MetaParam_File_react()) %in% c("csv","CSV")) {
                #metaP <- as.data.frame(read_delim(MetaParam_File_react(), delim = ',', col_names = F))
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
                #metaP <- as.data.frame(read_delim(url(MetaParam_File_react()), delim = '\t', col_names = F))
                metaP <- as.data.frame(fread(url(MetaParam_File_react()), sep = '\t', header = F))
              } else if (tools::file_ext(MetaParam_File_react()) %in% c("zip","gz","ZIP","GZ")) {
                #metaP <- as.data.frame(read_delim(getZip(MetaParam_File_react()), delim = '\t', col_names = F))
                metaP <- as.data.frame(fread(getZip(MetaParam_File_react()), sep = '\t', header = F))
              } else if (tools::file_ext(MetaParam_File_react()) %in% c("csv","CSV")) {
                #metaP <- as.data.frame(read_delim(url(MetaParam_File_react()), delim = ',', col_names = F))
                metaP <- as.data.frame(fread(url(MetaParam_File_react()), sep = ',', header = F))
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
              #metaP <- as.data.frame(read_delim(MetaParam_File_react(), delim = '\t', col_names = F))
              metaP <- as.data.frame(fread(MetaParam_File_react(), sep = '\t', header = F))
            } else if (tools::file_ext(MetaParam_File_react()) %in% c("csv","CSV")) {
              #metaP <- as.data.frame(read_delim(MetaParam_File_react(), delim = ',', col_names = F))
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
              #expr <- as.data.frame(read_delim(ExpressionMatrix_file_react(), delim = '\t', col_names = T))
              expr <- as.data.frame(fread(ExpressionMatrix_file_react(), sep = '\t', header = T))
            } else if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("csv","CSV")) {
              #expr <- as.data.frame(read_delim(ExpressionMatrix_file_react(), delim = ',', col_names = T))
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
              #expr <- as.data.frame(fread(url(ExpressionMatrix_file_react()), sep = '\t', header = T))
            } else if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("zip","gz","ZIP","GZ")) {
              expr <- as.data.frame(read_delim(getZip(ExpressionMatrix_file_react()), delim = '\t', col_names = T))
              #expr <- as.data.frame(fread(getZip(ExpressionMatrix_file_react()), sep = '\t', header = T))
            } else if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("csv","CSV")) {
              expr <- as.data.frame(read_delim(url(ExpressionMatrix_file_react()), delim = ',', col_names = T))
              #expr <- as.data.frame(fread(url(ExpressionMatrix_file_react()), sep = ',', header = T))
            } else if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("RData","rdata")) {
              expr <- loadRData(url(ExpressionMatrix_file_react()))
            } else if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("rds","RDS")) {
              expr <- readRDS(url(ExpressionMatrix_file_react()))
            }
          }
          # Remove Expression with NA
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
          if (TRUE %in% duplicated(expr[,1])) {
            expr <- expr %>%
              group_by(Gene) %>%
              summarise_all(max) %>%
              as.data.frame()
          }
          row.names(expr) <- expr[,1]
          expr <- expr[,-1]
          
          expr <- as.matrix(expr[order(rownames(expr)),])
          
          # Meta file
          # If Meta is file
          if (file.exists(MetaData_file_react())) {
            print(paste0("Loading in local file: ",MetaData_file_react()))
            if (tools::file_ext(MetaData_file_react()) %in% c("txt","tsv","zip","gz","TXT","TSV","ZIP","GZ")) {
              #meta <- as.data.frame(read_delim(MetaData_file_react(), delim = '\t', col_names = T))
              meta <- as.data.frame(fread(MetaData_file_react(), sep = '\t', header = T))
            } else if (tools::file_ext(MetaData_file_react()) %in% c("csv","CSV")) {
              #meta <- as.data.frame(read_delim(MetaData_file_react(), delim = ',', col_names = T))
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
              #meta <- as.data.frame(fread(url(MetaData_file_react()), sep = '\t', header = T))
            } else if (tools::file_ext(MetaData_file_react()) %in% c("zip","gz","ZIP","GZ")) {
              meta <- as.data.frame(read_delim(getZip(MetaData_file_react()), delim = '\t', col_names = T))
              #meta <- as.data.frame(fread(getZip(MetaData_file_react()), sep = '\t', header = T))
            } else if (tools::file_ext(MetaData_file_react()) %in% c("csv","CSV")) {
              meta <- as.data.frame(read_delim(url(MetaData_file_react()), delim = ',', col_names = T))
              #meta <- as.data.frame(fread(url(MetaData_file_react()), sep = ',', header = T))
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
            
            #metacol_samplename <- metaP[which(metaP[,2] == "SampleName"),1]
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
        meta <- meta_react()
        decon_score_cols <- grep("_PreProcessedScore$", colnames(meta), value = T, ignore.case = T)
        decon_score_cols
      })
      
      observe({
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
          }
          if (input$SurvTimeUnits == "Years") {
            for (i in TimCols) {
              meta[,i] <- meta[,i] * 365.25
            }
          }
          meta_react(meta)
        }
        
      })
      
      observeEvent(input$LogExprFile | input$ScaleNormExprFile,{
        if (isTruthy(input$LogExprFile) | isTruthy(input$ScaleNormExprFile)) {
          if (input$LogExprFile & input$ScaleNormExprFile) {
            meta <- meta_react()
            metaP <- metaP_react()
            meta <- meta[,grep("_InApp_PreProcessedScore$",colnames(meta), invert = T)]
            metaP <- metaP[which(metaP[,1] %in% grep("_InApp_PreProcessedScore$",metaP[,1], invert = T, value = T)),]
            meta_react(meta)
            metaP_react(metaP)
            expr <- expr_raw()
            expr <- log2(expr+1)
            expr_col <- colnames(expr)
            expr <- apply(expr, 1, scale)
            expr <- apply(expr, 1, rev)
            colnames(expr) <- expr_col
            expr <- as.matrix(expr[order(rownames(expr)),])
            expr_react(expr)
          } else if (input$LogExprFile & !input$ScaleNormExprFile) {
            meta <- meta_react()
            metaP <- metaP_react()
            meta <- meta[,grep("_InApp_PreProcessedScore$",colnames(meta), invert = T)]
            metaP <- metaP[which(metaP[,1] %in% grep("_InApp_PreProcessedScore$",metaP[,1], invert = T, value = T)),]
            meta_react(meta)
            metaP_react(metaP)
            expr <- expr_raw()
            expr <- log2(expr+1)
            expr <- as.matrix(expr[order(rownames(expr)),])
            expr_react(expr)
          } else if (!input$LogExprFile & input$ScaleNormExprFile) {
            meta <- meta_react()
            metaP <- metaP_react()
            meta <- meta[,grep("_InApp_PreProcessedScore$",colnames(meta), invert = T)]
            metaP <- metaP[which(metaP[,1] %in% grep("_InApp_PreProcessedScore$",metaP[,1], invert = T, value = T)),]
            meta_react(meta)
            metaP_react(metaP)
            expr <- expr_raw()
            expr_col <- colnames(expr)
            expr <- apply(expr, 1, scale)
            expr <- apply(expr, 1, rev)
            colnames(expr) <- expr_col
            expr <- as.matrix(expr[order(rownames(expr)),])
            expr_react(expr)
          } else {
            meta <- meta_react()
            metaP <- metaP_react()
            meta <- meta[,grep("_InApp_PreProcessedScore$",colnames(meta), invert = T)]
            metaP <- metaP[which(metaP[,1] %in% grep("_InApp_PreProcessedScore$",metaP[,1], invert = T, value = T)),]
            meta_react(meta)
            metaP_react(metaP)
            expr <- expr_raw()
            expr_react(expr)
          }
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
        metacol_feature <- c(metacol_feature,geneset_name,"QuartileCutP","MedianCutP","OptimalCutP","TopBottomCutP","UserCutP")
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
          SampleTypeChoices <- c(SampleTypeChoices,"All Sample Types")
          selectInput("SampleTypeSelection",paste("Select Sample Type (",metacol_sampletype,"):",sep = ""),
                      choices = SampleTypeChoices, selected = PreSelect_SamplyType_react())
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
        meta <- metaSub()
        expr <- exprSub()
        geneset <- gs_react()
        geneset_name <- names(geneset)
        if (isTruthy(input$SurvivalType_time) & isTruthy(input$SurvivalType_id) & isTruthy(geneset_name)) {
          SampleNameCol <- metacol_samplename_react()
          geneset <- gs_react()
          geneset_name <- names(geneset)
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
            
            ## Subset columns needed for plot and rename for surv function
            meta_ssgsea_sdf <- merge(meta[,c(SampleNameCol,surv_time_col,surv_id_col)],ssGSEA[,c(SampleNameCol,geneset_name)])
            
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
            
            meta_ssGSEA <- merge(meta,ssGSEA)
            meta_ssGSEA
          }
          
        }
        
        
      })
      
      ## Median CutP -----------------------------------------------------------
      
      #SubsetSurvData <- function(df,time,id,feat,feat2 = NULL) {
      #  SampleNameCol <- colnames(df)[1]
      #  if (is.null(feat2)) {
      #    df <- df[,c(SampleNameCol,time,id,feat)]
      #    if ("TopBottomCutP" %in% c(feat,feat2)) {
      #      df <- df[which(df$TopBottomCutP != "BetweenCutoff"),]
      #    }
      #    colnames(df)[which(colnames(df) == time)] <- "time"
      #    colnames(df)[which(colnames(df) == id)] <- "ID"
      #  } else {
      #    df <- df[,c(SampleNameCol,time,id,feat,feat2)]
      #    if ("TopBottomCutP" %in% c(feat,feat2)) {
      #      df <- df[which(df$TopBottomCutP != "BetweenCutoff"),]
      #    }
      #    colnames(df)[which(colnames(df) == time)] <- "time"
      #    colnames(df)[which(colnames(df) == id)] <- "ID"
      #  }
      #  return(df)
      #}
      
      MedianCutP_react <- reactive({
        req(ssGSEAmeta())
        SubsetSurvData(ssGSEAmeta(),input$SurvivalType_time,input$SurvivalType_id,"MedianCutP")
      })
      
      #CoxPHobj <- function(df,feat,ref) {
      #  df[,feat] <- as.factor(df[,feat])
      #  df[,feat] <- relevel(df[,feat], ref = ref)
      #  tab <- coxph(as.formula(paste0("Surv(time,ID) ~ ",feat)),data = df)
      #  #tab <- coxph(Surv(time,ID) ~ Feature, data = df)
      #  return(tab)
      #}
      
      MedianCutPTab_react <- reactive({
        req(MedianCutP_react())
        CoxPHobj(MedianCutP_react(),"MedianCutP","Low")
      })
      
      #CoxPHtabUni <- function(obj) {
      #  obj <- obj %>% 
      #    gtsummary::tbl_regression(exp = TRUE) %>%
      #    as_gt()
      #  tab_df <- as.data.frame(obj)
      #  tab_df <- tab_df %>%
      #    dplyr::select(label,estimate,ci,p.value)
      #  colnames(tab_df) <- c("Characteristic","Hazard Ratio","95% Confidence Interval","P.Value")
      #  tab_df <- sapply(tab_df,function(x) { gsub("<br />", "", x) })
      #  return(tab_df)
      #}
      
      SBinaryHRtab_react <- reactive({
        req(MedianCutPTab_react())
        #CoxPHtabUni(MedianCutPTab_react(),"MedianCutP")
        CoxPHtabUni(MedianCutPTab_react())
      })
      
      #SurvPlotExpl <- function(CutPlabel,surv_time_col,geneset_name,scoreMethod,metacol_sampletype,SampleTypeSelected,Feature,subFeature,Pval_Tab,HR_Tab) {
      #  SurvDateType <- sub("\\..*","",surv_time_col)
      #  pval <- get_lik_pval(Pval_Tab)
      #  if (as.numeric(pval) < 0.05) {
      #    pval_char <- "strongly associated"
      #  } else if (as.numeric(pval) >= 0.05 & as.numeric(pval) < 0.1) {
      #    pval_char <- "moderately associated"
      #  } else if (as.numeric(pval) >= 0.1) {
      #    pval_char <- "not associated"
      #  }
      #  if (isTruthy(SampleTypeSelected)) {
      #    if (SampleTypeSelected != "All Sample Types") {
      #      if (Feature != "Show all Samples") {
      #        line1 <- paste0("<li><b>",SurvDateType,"</b> survival analysis ", metacol_sampletype," of <b>",SampleTypeSelected,"</b> Patients.</li>")
      #        line2 <- paste0("<li>The dataset is filtered by <b>",Feature,"</b> - <b>",subFeature,"</b>.</li>")
      #      } else {
      #        line1 <- paste0("<li><b>",SurvDateType,"</b> survival analysis ", metacol_sampletype," of <b>",SampleTypeSelected,"</b> Patients.</li>")
      #        line2 <- NULL
      #      }
      #    } else {
      #      if (Feature != "Show all Samples") {
      #        line1 <- paste0("<li><b>",SurvDateType,"</b> survival analysis of all sample types.</li>")
      #        line2 <- paste0("<li>The dataset is filtered by <b>",Feature,"</b> - <b>",subFeature,"</b>.</li>")
      #      } else {
      #        line1 <- paste0("<li><b>",SurvDateType,"</b> survival analysis of all sample types.</li>")
      #        line2 <- NULL
      #      }
      #    }
      #  } else {
      #    if (Feature != "Show all Samples") {
      #      line1 <- paste0("<li><b>",SurvDateType,"</b> survival analysis.</li>")
      #      line2 <- paste0("<li>The dataset is filtered by <b>",Feature,"</b> - <b>",subFeature,"</b>.</li>")
      #    } else {
      #      line1 <- paste0("<li><b>",SurvDateType,"</b> survival analysis.</li>")
      #      line2 <- NULL
      #    }
      #  }
      #  line3 <- paste0("<li>Kaplan-Meier survival curve categorized by <b>",geneset_name,"</b> <b>",scoreMethod,"</b> ",CutPlabel,".</li>")
      #  if (nrow(HR_Tab) == 3) {
      #    HR <- HR_Tab[3,2]
      #    HR_char <- ifelse(as.numeric(gsub(",","",(HR))) > 1,"high risk","low risk")
      #    chacteristic <- str_squish(HR_Tab[3,1])
      #    line4 <- paste("<li>Cox hazard regression analysis finds a Likelihood Ratio P.value of <b>",pval,"</b> and a Hazard Ratio of <b>",HR,"</b>, <b>",chacteristic,"</b> <b>",geneset_name,
      #                   "</b> is <b>",pval_char,"</b> with <b>",HR_char,"</b> for <b>",SurvDateType,"</b>.</li>",sep = "")
      #  } else {
      #    line4 <- paste("<li>Cox hazard regression analysis finds a Likelihood Ratio P.value of <b>",pval,"</b> shows that <b>",geneset_name,
      #                   "</b> is <b>",pval_char,"</b> with <b>",SurvDateType,"</b>.</li>",sep = "")
      #  }
      #  if (is.null(line2)) {
      #    return(HTML(paste0("<ul>",line1,line3,line4,"</ul>")))
      #  } else if (!is.null(line2)) {
      #    return(HTML(paste0("<ul>",line1,line2,line3,line4,"</ul>")))
      #  }
      #  
      #}
      
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
      
      #SurvPlot <- function(fit,df,title,ylab,pval,conf,legend,median,xlim) {
      #  breakTime <- ifelse(max(df[,"time"]) < 365.25,NULL,365.25)
      #  ggsurv <- survminer::ggsurvplot(fit, data = df, risk.table = TRUE,
      #                                  title = title,
      #                                  xscale = c("d_y"),
      #                                  break.time.by=breakTime,
      #                                  xlab = "Years", 
      #                                  ylab = ylab,
      #                                  submain = "Based on Kaplan-Meier estimates",
      #                                  caption = "created with survminer",
      #                                  pval = pval,
      #                                  conf.int = conf,
      #                                  ggtheme = theme_bw(),
      #                                  font.title = c(16, "bold"),
      #                                  font.submain = c(12, "italic"),
      #                                  font.caption = c(12, "plain"),
      #                                  font.x = c(14, "plain"),
      #                                  font.y = c(14, "plain"),
      #                                  font.tickslab = c(12, "plain"),
      #                                  legend = legend,
      #                                  risk.table.height = 0.20,
      #                                  surv.median.line = median
      #  )
      #  if (median != "none") {
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
      #  if (!is.null(xlim)) {
      #    ggsurv$plot$coordinates$limits$x <- c(0,xlim)
      #    ggsurv$table$coordinates$limits$x <- c(0,xlim)
      #  }
      #  ggsurv$table <- ggsurv$table + theme_cleantable()
      #  ggsurv
      #}
      
      #SurvPlotTitle <- function(SampleTypeSelected,geneset_name = NULL,scoreMethodLab,Feature,subFeature,CutPLabel,univar = NULL,multivar = NULL) {
      #  if (isTruthy(geneset_name)) {
      #    primeFeatLab <- paste0(geneset_name," (",scoreMethodLab,") at ",CutPLabel)
      #  }
      #  if (isTruthy(univar)) {
      #    primeFeatLab <- paste0(univar)
      #  }
      #  if (isTruthy(multivar)) {
      #    primeFeatLab <- paste0(multivar)
      #  }
      #  if (isTruthy(SampleTypeSelected)) {
      #    if (SampleTypeSelected != "All Sample Types") {
      #      if (Feature != "Show all Samples") {
      #        PlotTitle <- paste0("Survival Curve of ",primeFeatLab,"\nAcross Patients Classified as ",SampleTypeSelected," - ",Feature," (",subFeature,")")
      #      } else {
      #        PlotTitle <- paste0("Survival Curve of ",primeFeatLab,"\nAcross Patients Classified as ",SampleTypeSelected)
      #      }
      #    } else {
      #      if (Feature != "Show all Samples") {
      #        PlotTitle <- paste0("Survival Curve of ",primeFeatLab,"\nAcross Patients Classified as ",Feature," (",subFeature,")")
      #      } else {
      #        PlotTitle <- paste0("Survival Curve of ",primeFeatLab,"\nAcross All Patients")
      #      }
      #    }
      #  } else {
      #    if (Feature != "Show all Samples") {
      #      PlotTitle <- paste0("Survival Curve of ",primeFeatLab,"\nAcross Patients Classified as ",Feature," (",subFeature,")")
      #    } else {
      #      PlotTitle <- paste0("Survival Curve of ",primeFeatLab,"\nAcross All Patients")
      #    }
      #  }
      #  return(PlotTitle)
      #}
      
      SplotBIN_react <- reactive({
        req(MedianCutP_react())
        meta_ssgsea_sdf <- MedianCutP_react()
        geneset <- gs_react()
        geneset_name <- names(geneset)
        SurvFeature <- "MedianCutP"
        Feature <- input$FeatureSelection
        subFeature <- input$subFeatureSelection
        SampleTypeSelected <- input$SampleTypeChoices
        scoreMethod <- input$ScoreMethod
        show_pval <- input$ShowPval
        ShowConfInt <- input$ShowConfInt
        if (!is.null(input$SurvXaxis)) {
          xaxlim <- input$SurvXaxis * 365.25
        } else {
          xaxlim <- NULL
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
        
        form <- as.formula(paste0("Surv(time,ID) ~ ",SurvFeature))
        fit <- eval(substitute(survfit(form,data = meta_ssgsea_sdf, type="kaplan-meier")))
        
        #fit <- survfit(Surv(time,ID) ~ Feature, data = meta_ssgsea_sdf, type="kaplan-meier")
        #attr(fit$strata,"names") <- gsub("Feature=",paste0(geneset_name,"="),attr(fit$strata,"names"))
        SurvPlot(fit,meta_ssgsea_sdf,PlotTitle,ylab = paste(SurvDateType,"Survival Probability"),
                 pval = show_pval,conf = ShowConfInt,legend = showLegend,median = showMedSurv,xlim = xaxlim)
      })
      output$SplotBIN <- renderPlot({
        plot <- SplotBIN_react()
        plot
      })
      
      #densPlot <- function(score,quant,xlab,ylab,title,CutPlabel,ShowQuartile = TRUE,user_vline = 0){
      #  dens_data <- density(score[,1],na.rm = T)
      #  y_max <- max(dens_data$y)
      #  y_max_int <- y_max/6
      #  p <- ggplot(score, aes(x=score[,1])) + 
      #    geom_density(color="darkblue", fill="lightblue", alpha = 0.4) +
      #    xlab(xlab) +
      #    ylab(ylab) +
      #    ggtitle(title) +
      #    theme(axis.text = element_text(size = 14),
      #          axis.title = element_text(size = 16),
      #          plot.title = element_text(size = 20))
      #  if (!is.null(CutPlabel)) {
      #    p <- p + geom_vline(data = quant, aes(xintercept = Quantile), linetype = "dashed", color = "darkblue", linewidth = 1)
      #    if (nrow(quant) == 1) {
      #      p <- p + geom_text(aes(quant[1,1],y_max-(y_max_int/6),label = paste(as.character(quant[1,1]),CutPlabel),hjust = -0.1,vjust = -0.1),size = 6, check_overlap = T)
      #    } else if (nrow(quant) == 2) {
      #      p <- p + geom_text(aes(quant[1,1],y_max-(y_max_int/6),label = paste(as.character(quant[1,1]),CutPlabel[1]),hjust = -0.1,vjust = -0.1),size = 6, check_overlap = T)
      #      p <- p + geom_text(aes(quant[2,1],y_max-y_max_int,label = paste(as.character(quant[2,1]),CutPlabel[2]),hjust = -0.1,vjust = -0.1),size = 6, check_overlap = T)
      #    } else if (nrow(quant) == 3) {
      #      p <- p + geom_text(aes(quant[1,1],y_max-(y_max_int/6),label = paste(as.character(quant[1,1]),CutPlabel[1],sep="\n"),hjust = -0.1,vjust = 0.5),size = 6, check_overlap = T)
      #      p <- p + geom_text(aes(quant[2,1],y_max-y_max_int,label = paste(as.character(quant[2,1]),CutPlabel[2],sep="\n"),hjust = -0.1,vjust = 0.5),size = 6, check_overlap = T)
      #      p <- p + geom_text(aes(quant[3,1],y_max-(y_max_int*2),label = paste(as.character(quant[3,1]),CutPlabel[3],sep="\n"),hjust = -0.1,vjust = 0.5),size = 6, check_overlap = T)
      #    }
      #  } else {
      #    if (ShowQuartile == TRUE) {
      #      p <- p + geom_vline(data = quant, aes(xintercept = Quantile), linetype = "dashed", color = "darkblue", linewidth = 1)
      #    }
      #    if (user_vline != 0) {
      #      p <- p + geom_vline(xintercept = user_vline, linetype = "dashed", color = "darkred", linewidth = 1)
      #    }
      #  }
      #  
      #  return(p)
      #}
      
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
      
      #CoxPHsumm <- function(CoxPHobj,bivarAdd = FALSE,bivarInt = FALSE) {
      #  out <- capture.output(summary(CoxPHobj))
      #  xph <- capture.output(cox.zph(CoxPHobj))
      #  con_line <- grep("^Concordance=",out,value = T)
      #  lik_line <- grep("^Likelihood ratio test=",out,value = T)
      #  wal_line <- grep("^Wald test",out,value = T)
      #  sco_line <- grep("^Score ",out,value = T)
      #  text <- paste("CoxH Summary:",con_line,lik_line,wal_line,sco_line,"","Proportional Hazards assumption:",xph[1],xph[2],xph[3],sep = "\n")
      #  if (bivarAdd) {
      #    text <- paste("CoxH Summary:",con_line,lik_line,wal_line,sco_line,"","Proportional Hazards assumption:",xph[1],xph[2],xph[3],xph[4],sep = "\n")
      #  }
      #  if (bivarInt) {
      #    text <- paste("CoxH Summary:",con_line,lik_line,wal_line,sco_line,"","Proportional Hazards assumption:",xph[1],xph[2],xph[3],xph[4],xph[5], sep = "\n")
      #  }
      #  return(cat(text))
      #}
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
        SubsetSurvData(ssGSEAmeta(),input$SurvivalType_time,input$SurvivalType_id,"QuartileCutP")
      })
      
      QuartileCutPTab_react <- reactive({
        req(QuartileCutP_react())
        CoxPHobj(QuartileCutP_react(),"QuartileCutP","Q1_Low")
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
        SurvFeature <- "QuartileCutP"
        Feature <- input$FeatureSelection
        subFeature <- input$subFeatureSelection
        SampleTypeSelected <- input$SampleTypeChoices
        scoreMethod <- input$ScoreMethod
        show_pval <- input$ShowPval
        ShowConfInt <- input$ShowConfInt
        if (!is.null(input$SurvXaxis)) {
          xaxlim <- input$SurvXaxis * 365.25
        } else {
          xaxlim <- NULL
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
        
        form <- as.formula(paste0("Surv(time,ID) ~ ",SurvFeature))
        fit <- eval(substitute(survfit(form,data = meta_ssgsea_sdf, type="kaplan-meier")))
        
        #fit <- survfit(Surv(time,ID) ~ Feature, data = meta_ssgsea_sdf, type="kaplan-meier")
        #attr(fit$strata,"names") <- gsub("Feature=",paste0(geneset_name,"="),attr(fit$strata,"names"))
        SurvPlot(fit,meta_ssgsea_sdf,PlotTitle,ylab = paste(SurvDateType,"Survival Probability"),
                 pval = show_pval,conf = ShowConfInt,legend = showLegend,median = showMedSurv,xlim = xaxlim)
      })
      output$Splot <- renderPlot({
        Splot_react()
      })
      ssgseaQuartDensity_react <- reactive({
        
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
        SubsetSurvData(ssGSEAmeta(),input$SurvivalType_time,input$SurvivalType_id,"OptimalCutP")
      })
      
      OptimalCutPTab_react <- reactive({
        req(OptimalCutP_react())
        CoxPHobj(OptimalCutP_react(),"OptimalCutP","low")
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
        SurvFeature <- "OptimalCutP"
        Feature <- input$FeatureSelection
        subFeature <- input$subFeatureSelection
        SampleTypeSelected <- input$SampleTypeChoices
        scoreMethod <- input$ScoreMethod
        show_pval <- input$ShowPval
        ShowConfInt <- input$ShowConfInt
        if (!is.null(input$SurvXaxis)) {
          xaxlim <- input$SurvXaxis * 365.25
        } else {
          xaxlim <- NULL
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
        form <- as.formula(paste0("Surv(time,ID) ~ ",SurvFeature))
        fit <- eval(substitute(survfit(form,data = meta_ssgsea_sdf, type="kaplan-meier")))
        #fit <- survfit(Surv(time,ID) ~ Feature, data = meta_ssgsea_sdf, type="kaplan-meier")
        #attr(fit$strata,"names") <- gsub("Feature=",paste0(geneset_name,"="),attr(fit$strata,"names"))
        SurvPlot(fit,meta_ssgsea_sdf,PlotTitle,ylab = paste(SurvDateType,"Survival Probability"),
                 pval = show_pval,conf = ShowConfInt,legend = showLegend,median = showMedSurv,xlim = xaxlim)
      })
      output$ScutPointPlot <- renderPlot({
        ScutPointPlot_react()
      })
      ssgseaCutPDensity_react <- reactive({
        
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
          res.cut <- survminer::surv_cutpoint(meta_ssgsea_sdf,time = "time", event = "ID", variable = geneset_name)
          cutp <- round(res.cut$cutpoint[["cutpoint"]],3) 
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
        meta_ssgsea_sdf <- SubsetSurvData(ssGSEAmeta(),input$SurvivalType_time,input$SurvivalType_id,"TopBottomCutP")
        meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf[,"TopBottomCutP"] != "BetweenCutoff"),]
        meta_ssgsea_sdf
      })
      
      TopBottomCutPTab_react <- reactive({
        req(TopBottomCutP_react())
        CoxPHobj(TopBottomCutP_react(),"TopBottomCutP","Low")
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
        SurvFeature <- "TopBottomCutP"
        Feature <- input$FeatureSelection
        subFeature <- input$subFeatureSelection
        SampleTypeSelected <- input$SampleTypeChoices
        scoreMethod <- input$ScoreMethod
        show_pval <- input$ShowPval
        ShowConfInt <- input$ShowConfInt
        if (!is.null(input$SurvXaxis)) {
          xaxlim <- input$SurvXaxis * 365.25
        } else {
          xaxlim <- NULL
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
        
        form <- as.formula(paste0("Surv(time,ID) ~ ",SurvFeature))
        fit <- eval(substitute(survfit(form,data = meta_ssgsea_sdf, type="kaplan-meier")))
        #fit <- survfit(Surv(time,ID) ~ Feature, data = meta_ssgsea_sdf, type="kaplan-meier")
        #attr(fit$strata,"names") <- gsub("Feature=",paste0(geneset_name,"="),attr(fit$strata,"names"))
        SurvPlot(fit,meta_ssgsea_sdf,PlotTitle,ylab = paste(SurvDateType,"Survival Probability"),
                 pval = show_pval,conf = ShowConfInt,legend = showLegend,median = showMedSurv,xlim = xaxlim)
      })
      output$SquantPlot <- renderPlot({
        SquantPlot_react()
      })
      ssgseaQuantDensity_react <- reactive({
        
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
        meta_ssgsea_sdf <- SubsetSurvData(ssGSEAmeta(),input$SurvivalType_time,input$SurvivalType_id,"UserCutP")
        meta_ssgsea_sdf
      })
      
      UserCutPTab_react <- reactive({
        req(UserCutP_react())
        CoxPHobj(UserCutP_react(),"UserCutP","Low")
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
        SurvFeature <- "UserCutP"
        Feature <- input$FeatureSelection
        subFeature <- input$subFeatureSelection
        SampleTypeSelected <- input$SampleTypeChoices
        scoreMethod <- input$ScoreMethod
        show_pval <- input$ShowPval
        ShowConfInt <- input$ShowConfInt
        if (!is.null(input$SurvXaxis)) {
          xaxlim <- input$SurvXaxis * 365.25
        } else {
          xaxlim <- NULL
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
        form <- as.formula(paste0("Surv(time,ID) ~ ",SurvFeature))
        fit <- eval(substitute(survfit(form,data = meta_ssgsea_sdf, type="kaplan-meier")))
        #fit <- survfit(Surv(time,ID) ~ Feature, data = meta_ssgsea_sdf, type="kaplan-meier")
        #attr(fit$strata,"names") <- gsub("Feature=",paste0(geneset_name,"="),attr(fit$strata,"names"))
        SurvPlot(fit,meta_ssgsea_sdf,PlotTitle,ylab = paste(SurvDateType,"Survival Probability"),
                 pval = show_pval,conf = ShowConfInt,legend = showLegend,median = showMedSurv,xlim = xaxlim)
      })
      output$SquantPlot2 <- renderPlot({
        SquantPlot2_react()
      })
      ssgseaQuant2Density_react <- reactive({
        
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
        updateSelectizeInput(session = session,inputId = "SingleSurvivalFeature", choices = metacol_feature(), selected = "MedianCutP", server = T)
      })
      
      #survFeatRefSelect <- function(meta,Feature,na.rm = TRUE,cont = FALSE,hilo = TRUE) {
      #  Var_choices <- meta[,Feature]
      #  if (na.rm == TRUE) {
      #    Var_choices <- Var_choices[which(is.na(Var_choices) == FALSE)]
      #    Var_choices <- Var_choices[which(Var_choices != "Inf" & Var_choices != "N/A" & Var_choices != "n/a")]
      #    Var_choices <- Var_choices[grep("unknown",Var_choices,ignore.case = T, invert = T)]
      #  }
      #  if (cont == FALSE) {
      #    Var_choices <- unique(meta[,Feature])
      #    Var_choices <- sort(Var_choices, decreasing = T, na.last = T)
      #  } else if (cont == TRUE) {
      #    if (hilo == TRUE) {
      #      Var_choices <- c("Low","High")
      #    }
      #  }
      #  return(Var_choices)
      #}
      
      output$rendSurvFeatVariableUni <- renderUI({
        req(input$SingleSurvivalFeature)
        req(ssGSEAmeta())
        Feature <- input$SingleSurvivalFeature
        meta <- ssGSEAmeta()
        Var_choices <- survFeatRefSelect(meta,Feature,input$UniVarNAcheck,input$UniVarContCheck,input$UniVarContCheck)
        selectInput("SurvFeatVariableUni","Select Coxh Feature Reference:", choices = Var_choices)
      })
      
      #observe({
      #  Feature <- input$SingleSurvivalFeature
      #  meta_ssgsea <- ssGSEAmeta()
      #  if (isTruthy(meta_ssgsea[,Feature])) {
      #    if (is.numeric(meta_ssgsea[,Feature])) {
      #      updateSelectInput(session,"UniVarContCheck",selected = TRUE)
      #    }
      #  }
      #})
      
      UniVarFeat_react <- reactive({
        Feature <- input$SingleSurvivalFeature
        meta_ssgsea <- ssGSEAmeta()
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
        SubsetSurvData(meta_ssgsea,input$SurvivalType_time,input$SurvivalType_id,Feature)

      })
      
      UniVarFeatTab_react <- reactive({
        
        meta_ssgsea_sdf <- UniVarFeat_react()
        Feature <- input$SingleSurvivalFeature
        ref_Feature <- input$SurvFeatVariableUni
        if ((input$UniVarContCheck == FALSE) | (input$UniVarContCheck == TRUE & input$UniVarContHiLoCheck == TRUE)) {
          tab <- CoxPHobj(meta_ssgsea_sdf,Feature,ref_Feature)
        } else {
          tab <- coxph(as.formula(paste0("Surv(time,ID) ~ ",Feature)),data = meta_ssgsea_sdf)
          #tab <- coxph(Surv(time,ID) ~ Feature ,data = meta_ssgsea_sdf)
        }
        tab
        
      })
      
      SSingleFeatureHRtab_react <- reactive({
        
        #Feature <- input$SingleSurvivalFeature
        tab <- UniVarFeatTab_react()
        CoxPHtabUni(tab)
        
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
        if (!is.null(input$SurvXaxis)) {
          xaxlim <- input$SurvXaxis * 365.25
        } else {
          xaxlim <- NULL
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
        
        form <- as.formula(paste0("Surv(time,ID) ~ ",Feature))
        fit <- eval(substitute(survfit(form,data = meta_ssgsea_sdf, type="kaplan-meier")))
        
        #fit <- survfit(Surv(time,ID) ~ Feature, data = meta_ssgsea_sdf, type="kaplan-meier")
        #attr(fit$strata,"names") <- gsub("Feature=",paste0(Feature,"="),attr(fit$strata,"names"))
        SurvPlot(fit,meta_ssgsea_sdf,PlotTitle,ylab = paste(SurvDateType,"Survival Probability"),
                 pval = show_pval,conf = ShowConfInt,legend = showLegend,median = showMedSurv,xlim = xaxlim)
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
        CoxPHsumm(tab)
        
      })
      
      #forestPlot_Simple <- function(obj,df,Feature,Font) {
      #  forest <- survminer::ggforest(obj,
      #                                data = df,
      #                                main = paste("Hazard Ratio Modeling: ",Feature,sep = ""),
      #                                fontsize = Font)
      #  return(forest)
      #}
      
      SinglevarForestPlot_react <- reactive({
        forestPlot_Simple(UniVarFeatTab_react(),UniVarFeat_react(),input$SingleSurvivalFeature,input$ForestFontSize)
      })
      output$SinglevarForestPlot <- renderPlot({
        forest <- SinglevarForestPlot_react()
        forest
      })
      
      observe({
        req(ssGSEAmeta())
        meta <- ssGSEAmeta()
        ColLevels <- apply(meta,2,function(x) length(levels(as.factor(x))))
        DicotCols <- names(ColLevels[ColLevels==2])
        ContCols <- GetColsOfType(meta,"continuous")
        FeatureChoices <- unique(c(DicotCols,ContCols))
        SurvDataCols <- c(metacol_survid,metacol_survtime)
        FeatureChoices <- FeatureChoices[which(!FeatureChoices %in% SurvDataCols)]
        updateSelectizeInput(session,"MultiFeatUnivarSelect",choices = FeatureChoices,selected = "MedianCutP", server = T)
      })
      
      
      MultiFeat_ForestMeta <- reactive({
        
        meta <- ssGSEAmeta()
        surv_time_col <- input$SurvivalType_time
        surv_id_col <- input$SurvivalType_id
        colnames(meta)[1] <- "SampleName"
        FeatCols <- c("SampleName",surv_time_col,surv_id_col,input$MultiFeatUnivarSelect)
        metaSub <- meta[,FeatCols]
        continuousCols <- GetColsOfType(metaSub[,-c(1:3), drop = F],"continuous")
        metaSub[,continuousCols] <- apply(metaSub[,continuousCols, drop = F],2,function(x) highlow2(x))
        if ("MedianCutP" %in% colnames(metaSub)) {
          metaSub[,"MedianCutP"] <- as.factor(metaSub[,"MedianCutP"])
          metaSub[,"MedianCutP"] <- relevel(metaSub[,"MedianCutP"], ref = "Low")
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
      #linearityPlot <- function(obj,Feature,resid,pred,axisFont,mainFont,tickFont) {
      #  p <- survminer::ggcoxdiagnostics(obj,
      #                                   type = resid,
      #                                   sline = T,
      #                                   sline.se = T,
      #                                   ggtheme = theme_minimal(),
      #                                   ox.scale = pred)
      #  p <- ggpar(p,
      #             font.x = axisFont,
      #             font.y = axisFont,
      #             font.main = mainFont,
      #             font.tickslab = tickFont,
      #             main = paste("Linearity Plot Featuring: ",Feature, sep = ""),
      #             ylab = paste(str_to_title(resid)," Residuals", sep = "")
      #  )
      #  p
      #}
      UnivarLinearityPlot_react <- reactive({
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
        FeatureChoices <- meta %>%
          dplyr::select(where(~ n_distinct(.x[nzchar(.x)], na.rm = TRUE) > 1)) %>%
          names
        #metacol_feature()
        updateSelectizeInput(session = session,inputId = "SurvivalFeatureBi1", choices = FeatureChoices, selected = "MedianCutP", server = T)
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
        Feature <- input$SurvivalFeatureBi1
        meta <- ssGSEAmeta()
        Var_choices <- survFeatRefSelect(meta,Feature,input$BiVarAddNAcheck1,input$BiVarAddContCheck1,input$BiVarAddContHiLoCheck1)
        selectInput("SurvFeatVariableBi1","Select Coxh Feature Reference:", choices = Var_choices)
      })
      output$rendSurvFeatVariableBi2 <- renderUI({
        Feature <- input$SurvivalFeatureBi2
        meta <- ssGSEAmeta()
        Var_choices <- survFeatRefSelect(meta,Feature,input$BiVarAddNAcheck2,input$BiVarAddContCheck2,input$BiVarAddContHiLoCheck2)
        selectInput("SurvFeatVariableBi2","Select Coxh Feature Reference:", choices = Var_choices)
      })
      
      BiVarAddFeature_react <- reactive({
        req(input$SurvFeatVariableBi1)
        req(input$SurvFeatVariableBi2)
        meta_ssgsea <- ssGSEAmeta()
        Feature1 <- input$SurvivalFeatureBi1
        Feature1ref <- input$SurvFeatVariableBi1
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
        
        Feature2 <- input$SurvivalFeatureBi2
        Feature2ref <- input$SurvFeatVariableBi2
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
      })
      BiVarAddTab_react <- reactive({
        Feature1 <- input$SurvivalFeatureBi1
        Feature2 <- input$SurvivalFeatureBi2
        tab <- coxph(as.formula(paste0("Surv(time,ID) ~ ",paste0(Feature1,"+",Feature2))),data = BiVarAddFeature_react())
        tab
      })
      BiVarAddTabFeat1_react <- reactive({
        Feature1 <- input$SurvivalFeatureBi1
        tab <- coxph(as.formula(paste0("Surv(time,ID) ~ ",Feature1)),data = BiVarAddFeature_react())
        tab
      })
      BiVarAddTabFeat2_react <- reactive({
        Feature2 <- input$SurvivalFeatureBi2
        tab <- coxph(as.formula(paste0("Surv(time,ID) ~ ",Feature2)),data = BiVarAddFeature_react())
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
      
      #biVarAnova <- function(obj,obj2) {
      #  annova_res <- anova(obj,obj2)
      #  out <- capture.output(annova_res)
      #  line1 <- out[3]
      #  line2 <- out[4]
      #  line3 <- out[5]
      #  line4 <- out[6]
      #  line5 <- out[7]
      #  text <- paste("Model Comparison:",line1,line2,line3,line4,line5,sep = "\n")
      #  return(cat(text))
      #}
      output$bivarAnova1 <- renderPrint({
        biVarAnova(BiVarAddTab_react(),BiVarAddTabFeat1_react())
      })
      
      output$bivarAnova2 <- renderPrint({
        biVarAnova(BiVarAddTab_react(),BiVarAddTabFeat2_react())
      })
      
      BivarForestPlot_react <- reactive({
        df <- BiVarAddFeature_react()
        Feature1 <- input$SurvivalFeatureBi1
        Feature2 <- input$SurvivalFeatureBi2
        forest <- forestPlot_Simple(BiVarAddTab_react(),df,paste0(Feature1,"+",Feature2),input$ForestFontSize)
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
        FeatureChoices <- meta %>%
          dplyr::select(where(~ n_distinct(.x[nzchar(.x)], na.rm = TRUE) > 1)) %>%
          names
        updateSelectizeInput(session = session,inputId = "SurvivalFeatureBi1Inter", choices = FeatureChoices, selected = "MedianCutP", server = T)
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
        meta <- ssGSEAmeta()
        Var_choices <- survFeatRefSelect(meta,Feature,input$BiVarIntNAcheck1,input$BiVarIntContCheck1,input$BiVarIntContHiLoCheck1)
        selectInput("SurvFeatVariableBi1Inter","Select Coxh Feature Reference:", choices = Var_choices)
      })
      output$rendSurvFeatVariableBi2Inter <- renderUI({
        Feature <- input$SurvivalFeatureBi2Inter
        meta <- ssGSEAmeta()
        Var_choices <- survFeatRefSelect(meta,Feature,input$BiVarIntNAcheck2,input$BiVarIntContCheck2,input$BiVarIntContHiLoCheck2)
        selectInput("SurvFeatVariableBi2Inter","Select Coxh Feature Reference:", choices = Var_choices)
      })
      
      BiVarIntFeature_react <- reactive({
        req(input$SurvivalFeatureBi1Inter)
        req(input$SurvivalFeatureBi2Inter)
        meta_ssgsea <- ssGSEAmeta()
        Feature1 <- input$SurvivalFeatureBi1Inter
        Feature1ref <- input$SurvFeatVariableBi1Inter
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
      })
      
      BiVarIntTab_react <- reactive({
        Feature1 <- input$SurvivalFeatureBi1Inter
        Feature2 <- input$SurvivalFeatureBi2Inter
        form <- as.formula(paste0("Surv(time,ID) ~ ",paste0(Feature1,"*",Feature2)))
        tab <- eval(substitute(coxph(form,data = BiVarIntFeature_react())))
        tab
      })
      BiVarIntTab4Annova_react <- reactive({
        Feature1 <- input$SurvivalFeatureBi1Inter
        Feature2 <- input$SurvivalFeatureBi2Inter
        form <- as.formula(paste0("Surv(time,ID) ~ ",paste0(Feature1,"+",Feature2)))
        tab <- eval(substitute(coxph(form,data = BiVarIntFeature_react())))
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
        if (!is.null(input$SurvXaxis)) {
          xaxlim <- input$SurvXaxis * 365.25
        } else {
          xaxlim <- NULL
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
        
        form <- as.formula(paste0("Surv(time,ID) ~ ",paste0(Feature1,"+",Feature2)))
        fit <- eval(substitute(survfit(form,data = meta_ssgsea_sdf, type="kaplan-meier")))
        
        SurvPlot(fit,meta_ssgsea_sdf,PlotTitle,ylab = paste(SurvDateType,"Survival Probability"),
                 pval = show_pval,conf = ShowConfInt,legend = showLegend,median = showMedSurv,xlim = xaxlim)
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
        FeatureChoices <- meta %>%
          dplyr::select(where(~ n_distinct(.x[nzchar(.x)], na.rm = TRUE) > 1)) %>%
          names
        updateSelectizeInput(session = session,inputId = "SurvivalFeature", choices = FeatureChoices, selected = "MedianCutP", server = T)
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
        meta_ssgsea_sdf <- MultiVarFeatCat_react()
        Feature <- input$SurvivalFeature
        form <- as.formula(paste0("Surv(time,ID) ~ ",paste(Feature,collapse = "+")))
        tab <- eval(substitute(coxph(form,data = meta_ssgsea_sdf)))
        tab
      })
      
      MultiVarTabCont_react <- reactive({
        meta_ssgsea_sdf <- MultiVarFeatCont_react()
        Feature <- input$SurvivalFeature
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
        ColLevels <- apply(meta,2,function(x) length(levels(as.factor(x))))
        DicotCols <- names(ColLevels[ColLevels==2])
        ContCols <- GetColsOfType(meta,"continuous")
        FeatureChoices <- unique(c(DicotCols,ContCols))
        SurvDataCols <- c(metacol_survid,metacol_survtime)
        FeatureChoices <- FeatureChoices[which(!FeatureChoices %in% SurvDataCols)]
        updateSelectizeInput(session,"MultiFeatMultivarSelect",choices = FeatureChoices, selected = "MedianCutP")
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
        if (input$ColorScatterChoice == "Feature") {
          Features <- c(metacol_feature(),metacol_survtime(),metacol_survid())
          updateSelectizeInput(session,"ScatterColor",choices = Features,selected = input$SurvivalType_id, server = T)
        }
      })
      
      FeatCompScatter_react <- reactive({
        req(ssGSEAmeta())
        req(gs_react())
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
          theme(text = element_text(size = font))
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
        
        geneset <- gs_react()
        GeneSet <- names(geneset)
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
                                                      #show_row_names = row_names_choice, show_column_names = col_names_choice,
                                                      #cluster_rows = clust_rows_opt,
                                                      cluster_columns = FALSE,
                                                      row_names_gp = gpar(fontsize = rowfont), column_names_gp = gpar(fontsize = colfont),
                                                      heatmap_legend_param = list(title = "Expression"),
                                                      border = F))
        draw(p, padding = unit(c(50, 50, 2, 2), "mm")) # unit(c(bottom,left,right,top))
        #draw(lgd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))
        
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
        if (length(input$stratHeatAnno) > 0) {
          meta <- ssGSEAmeta()
          meta <- meta[order(meta[,input$BoxplotFeature]),]
          rownames(meta) <- meta[,1]
          meta_sub <- meta[,input$stratHeatAnno, drop = F]
          meta_sub[,input$stratHeatAnno] <- as.data.frame(lapply(meta_sub[,input$stratHeatAnno, drop = F], factor))
          meta_sub <- meta_sub %>% relocate(any_of(input$BoxplotFeature), .after = last_col())
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
        
        geneset <- gs_react()
        GeneSet <- names(geneset)
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
                                                      #show_row_names = row_names_choice, show_column_names = col_names_choice,
                                                      #cluster_rows = clust_rows_opt,
                                                      cluster_columns = FALSE,
                                                      row_names_gp = gpar(fontsize = rowfont), column_names_gp = gpar(fontsize = colfont),
                                                      heatmap_legend_param = list(title = "Expression"),
                                                      border = F))
        draw(p, padding = unit(c(50, 50, 2, 2), "mm")) # unit(c(bottom,left,right,top))
        #draw(lgd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))
        
      })
      output$FeatureHeatmap <- renderPlot({
        heat <- FeatureHeatmap_react()
        heat
      })
      
      
      ####----Text Output----####
      
      output$timewarnmessage1 <- renderUI({
        
        if (input$linPredict1 == "time") {
          p("Residual type must be schoenfeld or scaledsch.")
        }
        
      })
      output$timewarnmessage1S <- renderUI({
        
        if (input$linPredict1S == "time") {
          p("Residual type must be schoenfeld or scaledsch.")
        }
        
      })
      
      output$timewarnmessage2 <- renderUI({
        
        if (input$linPredict2 == "time") {
          p("Residual type must be schoenfeld or scaledsch.")
        }
        
      })
      output$timewarnmessage2S <- renderUI({
        
        if (input$linPredict2S == "time") {
          p("Residual type must be schoenfeld or scaledsch.")
        }
        
      })
      
      output$timewarnmessage3 <- renderUI({
        
        if (input$linPredict3 == "time") {
          p("Residual type must be schoenfeld or scaledsch.")
        }
        
      })
      output$timewarnmessage3S <- renderUI({
        
        if (input$linPredict3S == "time") {
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
      
      output$dnldMultiFeatUnivarForestPlot_table <- downloadHandler(
        filename = function() {
          ProjectName <- gsub(" ","",ProjectName)
          ProjectName <- gsub("[[:punct:]]","",ProjectName)
          Todaydate <- gsub("-","",Sys.Date())
          paste0(ProjectName,"_UnivarForestPlotTable_",Todaydate,".txt")
        },
        content = function(file) {
          df <- MultiFeatUnivarForestPlotTab_react()
          df <- df[,-7]
          write_delim(df,file,delim = '\t', col_names = T)
        }
      )
      
      output$dnldMultiFeatUnivarForestPlot_SVG <- downloadHandler(
        filename = function() {
          ProjectName <- gsub(" ","",ProjectName)
          ProjectName <- gsub("[[:punct:]]","",ProjectName)
          Todaydate <- gsub("-","",Sys.Date())
          paste0(ProjectName,"_UnivarForestPlot_",Todaydate,".svg")
        },
        content = function(file) {
          #p_OS <- MultiFeatUnivarForestPlot_react()
          #dims <- get_wh(MultiFeatUnivarForestPlot_react())
          pdf(NULL)
          #p_OS <- MultiFeatUnivarForestPlot_react()
          dims <- get_wh(MultiFeatUnivarForestPlot_react())
          #pdf(file="RPlots.pdf")
          #svg(NULL)
          #dev.off()
          #dims <- get_wh(MultiFeatUnivarForestPlot_react())
          #pdf(file, width = unname(dims[1]), height = unname(dims[2]),family = "sans",pointsize = 12)
          svg(file, width = unname(dims[1]+1), height = unname(dims[2]+1))
          plot(MultiFeatUnivarForestPlot_react())
          dev.off()
        }
      )
      
      output$dnldunivarForestPlotTable <- downloadHandler(
        filename = function() {
          ProjectName <- gsub(" ","",ProjectName)
          ProjectName <- gsub("[[:punct:]]","",ProjectName)
          Todaydate <- gsub("-","",Sys.Date())
          paste0(ProjectName,"_UnivarForestPlot_MetaData_",Todaydate,".txt")
        },
        content = function(file) {
          df <- MultiFeat_ForestMeta()
          write_delim(df,file,delim = '\t', col_names = T)
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
      
      
      output$dnldMultiFeatMultivarForestPlot_table <- downloadHandler(
        filename = function() {
          ProjectName <- gsub(" ","",ProjectName)
          ProjectName <- gsub("[[:punct:]]","",ProjectName)
          Todaydate <- gsub("-","",Sys.Date())
          paste0(ProjectName,"_MultivarForestPlotTable_",Todaydate,".txt")
        },
        content = function(file) {
          df <- MultiFeatMultivarForestPlotTab_react()
          df <- df[,-7]
          df[is.na(df)] <- ""
          write_delim(df,file,delim = '\t', col_names = T)
        }
      )
      
      output$dnldMultiFeatMultivarForestPlot_SVG <- downloadHandler(
        filename = function() {
          ProjectName <- gsub(" ","",ProjectName)
          ProjectName <- gsub("[[:punct:]]","",ProjectName)
          Todaydate <- gsub("-","",Sys.Date())
          paste0(ProjectName,"_MultivarForestPlot_",Todaydate,".svg")
        },
        content = function(file) {
          #p_OS <- MultiFeatUnivarForestPlot_react()
          #dims <- get_wh(MultiFeatUnivarForestPlot_react())
          pdf(NULL)
          #p_OS <- MultiFeatUnivarForestPlot_react()
          dims <- get_wh(MultiFeatMultivarForestPlot_react())
          #pdf(file="RPlots.pdf")
          #svg(NULL)
          #dev.off()
          #dims <- get_wh(MultiFeatUnivarForestPlot_react())
          #pdf(file, width = unname(dims[1]), height = unname(dims[2]),family = "sans",pointsize = 12)
          svg(file, width = unname(dims[1]+1), height = unname(dims[2]+1))
          plot(MultiFeatMultivarForestPlot_react())
          dev.off()
        }
      )
      
      output$dnldmultiForestPlotTable <- downloadHandler(
        filename = function() {
          ProjectName <- gsub(" ","",ProjectName)
          ProjectName <- gsub("[[:punct:]]","",ProjectName)
          Todaydate <- gsub("-","",Sys.Date())
          paste0(ProjectName,"_MultivarForestPlot_MetaData_",Todaydate,".txt")
        },
        content = function(file) {
          df <- MultiFeat_InterForestMeta()
          write_delim(df,file,delim = '\t', col_names = T)
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
          if (Feature != "Show All Samples") {
            SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
          }
          if (Feature == "Show All Samples") {
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
          if (Feature != "Show All Samples") {
            SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
          }
          if (Feature == "Show All Samples") {
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
    
  }
  
  )
  
}


shinyApp(ui = ui, server = server)

