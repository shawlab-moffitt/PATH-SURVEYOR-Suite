

library(shiny)
library(readr)
library(dplyr)
library(tools)
library(DT)

#increase file upload size
options(shiny.maxRequestSize=5000*1024^2)

ui <-
  fluidPage(
    titlePanel("{ File Prep }"),
    sidebarPanel(
      h4("File Input"),
      fluidRow(
        column(8,
               fileInput("ExprFileInput","Expression File"),
        ),
        column(4,
               selectInput("ExprDelim","Delimeter", choices = c("Tab","Comma","Sapce"))
        )
      ),
      fluidRow(
        column(8, style = 'margin-top:-20px;',
               fileInput("MetaFileInput","Meta File"),
        ),
        column(4, style = 'margin-top:-20px;',
               selectInput("MetaDelim","Delimeter", choices = c("Tab","Comma","Sapce"))
        )
      ),
      #uiOutput("rendSampNameCol")
      h4("Parameter File Formatting"),
      fluidRow(
        column(8,
               uiOutput("rendSurvTimeColSelect"),
               uiOutput("rendSurvIDColSelect")
        ),
        column(4,
               uiOutput("rendSurvTimeUnits")
        )
      ),
      h4("Download Prepped Files"),
      downloadButton("dnldExprOutFile","Expression File"),
      downloadButton("dnldMetaOutFile","Meta File"),
      downloadButton("dnldParamOutFile","Meta Parameter File")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Expression Data Preview",
                 p(),
                 verbatimTextOutput("FileCheckAlerts"),
                 h4("Original Input Expression File"),
                 div(DT::dataTableOutput("InExprPreview"), style = "font-size:12px"),
                 h4("Reformatted Expression File"),
                 div(DT::dataTableOutput("OutExprPreview"), style = "font-size:12px")
        ),
        tabPanel("Meta Data Preview",
                 p(),
                 h4("Original Input Meta File"),
                 div(DT::dataTableOutput("InMetaPreview"), style = "font-size:12px"),
                 h4("Reformatted Meta File"),
                 div(DT::dataTableOutput("OutMetaPreview"), style = "font-size:12px")
        ),
        tabPanel("Meta Parameters Preview",
                 h4("Meta Parameter File"),
                 div(DT::dataTableOutput("MetaParamPreview"), style = "font-size:12px")
        )
      )
    )
  )


server <- function(input, output, session) {
  
  ####----Reactives----####
  
  exprIn <- reactive({
    
    exprFile <- input$ExprFileInput
    delim <- input$ExprDelim
    ext <- tools::file_ext(exprFile$datapath)
    req(exprFile)
    
    if (delim == "Tab") {
      expr <- as.data.frame(readr::read_delim(exprFile$datapath, delim = '\t', col_names = T))
    }
    else if (delim == "Comma") {
      expr <- as.data.frame(readr::read_delim(exprFile$datapath, delim = ',', col_names = T))
    }
    else if (delim == "Space") {
      expr <- as.data.frame(readr::read_delim(exprFile$datapath, delim = ' ', col_names = T))
    }
    expr
    
  })
  
  metaIn <- reactive({
    
    metaFile <- input$MetaFileInput
    delim <- input$MetaDelim
    ext <- tools::file_ext(metaFile$datapath)
    req(metaFile)
    
    if (delim == "Tab") {
      meta <- as.data.frame(readr::read_delim(metaFile$datapath, delim = '\t', col_names = T))
    }
    else if (delim == "Comma") {
      meta <- as.data.frame(readr::read_delim(metaFile$datapath, delim = ',', col_names = T))
    }
    else if (delim == "Space") {
      meta <- as.data.frame(readr::read_delim(metaFile$datapath, delim = ' ', col_names = T))
    }
    meta
    
  })
  
  exprReForm <- reactive({
    
    expr <- exprIn()
    meta <- metaIn()
    #metaSampCol <- input$SampNameCol
    
    ## Clean Names
    colnames(expr) <- gsub("[[:punct:]]",".",colnames(expr))
    #meta[,metaSampCol] <- gsub("[[:punct:]]",".",meta[,metaSampCol])
    meta[,1] <- gsub("[[:punct:]]",".",meta[,1])
    
    ## Check sample names match
    exprSamp <- colnames(expr)[-1]
    metaSamp <- meta[,1]
    SampSame <- intersect(exprSamp,metaSamp)
    GeneCol <- colnames(expr)[1]
    ColsToSelec <- c(GeneCol,SampSame)
    if (length(ColsToSelec) > 1) {
      expr2 <- expr[,ColsToSelec]
      
      ## Check that expression columns are numeric
      expr3 <- expr2
      isChar <- unname(which(sapply(expr3, function(x) is.character(x))))
      isChar <-  isChar[-1]
      if (length(isChar) > 0) {
        expr3[isChar] <- sapply(expr3[isChar],as.numeric)
      }
      
      ## Check for duplicated genes
      expr4 <- expr3
      colnames(expr4)[1] <- "Symbol"
      if (TRUE %in% duplicated(expr4[,1])) {
        expr4 <- expr4 %>%
          dplyr::group_by(Symbol) %>%
          dplyr::summarise_all(max)
      }
      colnames(expr4)[1] <- GeneCol
      
      ## Filter out rows with all 0
      expr5 <- expr4[which(rowSums(expr4[,c(2:ncol(expr4))], na.rm = T) > 0),]
      
      expr5
      
    }
    
  })
  
  metaReForm <- reactive({
    
    expr <- exprIn()
    meta <- metaIn()
    TimeUnits <- input$SurvTimeUnits
    TimCols <- input$SurvTimeColSelect
    
    ## Clean Names
    colnames(expr) <- gsub("[[:punct:]]",".",colnames(expr))
    meta[,1] <- gsub("[[:punct:]]",".",meta[,1])
    
    ## Check sample names match
    exprSamp <- colnames(expr)[-1]
    metaSamp <- meta[,1]
    SampSame <- intersect(exprSamp,metaSamp)
    meta <- meta[which(meta[,1] %in% SampSame),]
    
    ## Check survival time units
    if (!is.null(TimeUnits)) {
      if (TimeUnits != "Days") {
        if (TimeUnits == "Months") {
          for (i in TimCols) {
            meta[,i] <- meta[,i] * 30.4375
          }
        }
        if (TimeUnits == "Years") {
          for (i in TimCols) {
            meta[,i] <- meta[,i] * 365.25
          }
        }
      }
    }
    
    ## Remove column white space
    colnames(meta) <- gsub(" ","_",colnames(meta))
    
    meta
    
  })
  
  metaParam <- reactive({
    
    meta <- metaReForm()
    metacol_survtime <- input$SurvTimeColSelect
    metacol_survid <- input$SurvIDColSelect
    
    metaCols <- colnames(meta)
    SampNameCol <- colnames(meta)[1]
    
    if (length(metacol_survtime) > 0 && length(metacol_survid) > 0) {
      metacol_feature <- metaCols[-1]
      metacol_feature <- metacol_feature[!metacol_feature %in% c(metacol_survtime,metacol_survid)]
      
      MetaParam1 <- data.frame(Clinical_Column_Name = SampNameCol,
                               Clinical_Column_Type = "SampleName")
      MetaParam2 <- data.frame(Clinical_Column_Name = metacol_survtime,
                               Clinical_Column_Type = "SurvivalTime")
      MetaParam3 <- data.frame(Clinical_Column_Name = metacol_survid,
                               Clinical_Column_Type = "SurvivalID")
      MetaParam4 <- data.frame(Clinical_Column_Name = metacol_feature,
                               Clinical_Column_Type = "Feature")
      MetaParam <- rbind(MetaParam1,MetaParam2,MetaParam3,MetaParam4)
      MetaParam
    }
    
  })
  
  ####----Render UI----####
  
  #output$rendSampNameCol <- renderUI({
  #  
  #  meta <- metaIn()
  #  ColOpt <- colnames(meta)
  #  selectInput("SampNameCol","Select Sample Name Column", choices = ColOpt, selected = 1)
  #  
  #})
  
  output$rendSurvTimeColSelect <- renderUI({
    
    req(input$MetaFileInput)
    meta <- metaIn()
    ColumnCheck <- grep(paste(c("OS_time","EFS_time","PFI_time","PFS_time","RFS_time","DSS_time","DFI_time",
                                "OS.time","EFS.time","PFI.time","PFS.time","RFS.time","DSS.time","DFI.time"),collapse = "|"),
                        colnames(meta), ignore.case = T, value = T)
    ColCheckFound <- colnames(meta)[which(colnames(meta) %in% ColumnCheck)]
    if (length(ColCheckFound) == 0) {
      ColCheckFound <- ""
    }
    selectInput("SurvTimeColSelect","Select Survival Time Column(s)",
                choices = colnames(meta), selected = ColCheckFound, multiple = T)
    
  })
  
  output$rendSurvTimeUnits <- renderUI({
    
    req(input$MetaFileInput)
    selectizeInput("SurvTimeUnits","Time Units:", choices = c("Days","Months","Years"))
    
  })
  
  output$rendSurvIDColSelect <- renderUI({
    
    req(input$MetaFileInput)
    meta <- metaIn()
    ColumnCheckTime <- grep(paste(c("OS_time","EFS_time","PFI_time","PFS_time","RFS_time","DSS_time","DFI_time",
                                    "OS.time","EFS.time","PFI.time","PFS.time","RFS.time","DSS.time","DFI.time"),collapse = "|"),
                            colnames(meta), ignore.case = T, value = T)
    ColumnCheck <- grep(paste(c("^OS","^EFS","^PFI","^PFS","^RFS","^DSS","^DFI",
                                "OS_ID","EFS_ID","PFI_ID","PFS_ID","RFS_ID","DSS_ID","DFI_ID",
                                "OS.ID","EFS.ID","PFI.ID","PFS.ID","RFS.ID","DSS.ID","DFI.ID"),collapse = "|"),
                        colnames(meta), ignore.case = T, value = T)
    ColumnCheck <- ColumnCheck[which(!ColumnCheck %in% ColumnCheckTime)]
    ColCheckFound <- colnames(meta)[which(colnames(meta) %in% ColumnCheck)]
    if (length(ColCheckFound) == 0) {
      ColCheckFound <- ""
    }
    selectInput("SurvIDColSelect","Select Survival ID Column(s)",
                choices = colnames(meta), selected = ColCheckFound, multiple = T)
    
  })
  
  ####----Text Output----####
  
  output$FileCheckAlerts <- renderPrint({
    
    expr <- exprIn()
    meta <- metaIn()
    TimeUnits <- input$SurvTimeUnits
    
    ## Clean Names
    colnames(expr) <- gsub("[[:punct:]]",".",colnames(expr))
    meta[,1] <- gsub("[[:punct:]]",".",meta[,1])
    
    ## Check that samples are the same in Expr and Meta
    SampSame <- intersect(colnames(expr)[-1],meta[,1])
    if (length(SampSame) != length(colnames(expr)[-1]) | length(SampSame) != length(rownames(meta))) {
      SampCheck_line1 <- print("Uneven number of samples or mismatching names between expression and meta file.")
      SampCheck_line2 <- print(paste("Number of matching sample names:",length(SampSame)))
      SampCheck_line2 <- print("New files based on similar samples")
      SampCheck_lines <- paste(SampCheck_line1,SampCheck_line2,SampCheck_line3,sep="\n")
    } else {
      SampCheck_lines <- "All sample names match between expression and meta data."
    }
    
    ## Check for duplicated genes
    if (TRUE %in% duplicated(expr[,1])) {
      DupGene_Line <- "Duplicated genes found. Summarizing to gene with highest average expression."
      colnames(expr)[1] <- "Symbol"
      if (TRUE %in% duplicated(expr[,1])) {
        expr <- expr %>%
          dplyr::group_by(Symbol) %>%
          dplyr::summarise_all(max)
      }
    } else {
      DupGene_Line <- "No duplicated genes found."
    }
    
    ## Lowly expressed countisChar <- unname(which(sapply(expr3, function(x) is.character(x))))
    isChar <- unname(which(sapply(expr, function(x) is.character(x))))
    isChar <-  isChar[-1]
    if (length(isChar) > 0) {
      expr[isChar] <- sapply(expr[isChar],as.numeric)
    }
    expr_0 <- length(which(rowSums(expr[,c(2:ncol(expr))], na.rm = T) == 0))
    if (expr_0 > 0) {
      Expr0_Line <- paste0(expr_0," genes found expression level of zero across all samples.")
    } else {
      Expr0_Line <- "All genes have an average expression level above zero."
    }
    
    ## Check survival time column conversion
    TimeUnit_Line <- NULL
    if (!is.null(TimeUnits)) {
      if (TimeUnits != "Days") {
        TimeUnit_Line <- "Survival time column(s) converted to days."
      }
      else {
        TimeUnit_Line <- NULL
      }
    }
    
    text <- paste(SampCheck_lines,DupGene_Line,Expr0_Line,TimeUnit_Line, sep = "\n")
    cat(text)
    
  })
  
  ####----Data Tables----####
  
  
  output$InExprPreview <- DT::renderDataTable({
    
    df <- exprIn()
    
    DT::datatable(df,
                  extensions = "FixedColumns",
                  options = list(lengthMenu = c(5,10,20,50,100,1000,5000,10000),
                                 pageLength = 10,
                                 scrollX = TRUE,
                                 autoWidth = TRUE,
                                 fixedColumns = list(leftColumns = 1)),
                  rownames = F,
                  selection=list(mode = "multiple")) %>%
      formatRound(columns = c(2:ncol(df)), digits = 4)
    
  })
  
  output$InMetaPreview <- DT::renderDataTable({
    
    df <- metaIn()
    
    DT::datatable(df,
                  extensions = "FixedColumns",
                  options = list(lengthMenu = c(5,10,20,50,100,1000,5000,10000),
                                 pageLength = 10,
                                 scrollX = TRUE,
                                 autoWidth = TRUE,
                                 fixedColumns = list(leftColumns = 1)),
                  rownames = F,
                  selection=list(mode = "multiple"))
    
  })
  
  output$OutExprPreview <- DT::renderDataTable({
    
    df <- exprReForm()
    
    DT::datatable(df,
                  extensions = "FixedColumns",
                  options = list(lengthMenu = c(5,10,20,50,100,1000,5000,10000),
                                 pageLength = 10,
                                 scrollX = TRUE,
                                 autoWidth = TRUE,
                                 fixedColumns = list(leftColumns = 1)),
                  rownames = F,
                  selection=list(mode = "multiple")) %>%
      formatRound(columns = c(2:ncol(df)), digits = 4)
    
  })
  
  output$OutMetaPreview <- DT::renderDataTable({
    
    df <- metaReForm()
    
    DT::datatable(df,
                  extensions = "FixedColumns",
                  options = list(lengthMenu = c(5,10,20,50,100,1000,5000,10000),
                                 pageLength = 10,
                                 scrollX = TRUE,
                                 autoWidth = TRUE,
                                 fixedColumns = list(leftColumns = 1)),
                  rownames = F,
                  selection=list(mode = "multiple"))
    
  })
  
  output$MetaParamPreview <- DT::renderDataTable({
    
    df <- metaParam()
    
    DT::datatable(df,
                  options = list(lengthMenu = c(5,10,20,50,100,1000,5000,10000),
                                 pageLength = 20,
                                 scrollX = TRUE),
                  rownames = F,
                  selection=list(mode = "multiple"))
    
  })
  
  ####----Download Handelers----####
  
  output$dnldExprOutFile <- downloadHandler(
    filename = function() {
      exprFile <- input$ExprFileInput
      ext <- tools::file_ext(exprFile$datapath)
      fileName <- gsub(paste0(".",ext), "", basename(exprFile$name))
      paste0(fileName,"_Cleaned.txt")
    },
    content = function(file) {
      expr <- exprReForm()
      write_delim(expr,file,delim = '\t')
    }
  )
  
  output$dnldMetaOutFile <- downloadHandler(
    filename = function() {
      metaFile <- input$MetaFileInput
      ext <- tools::file_ext(metaFile$datapath)
      fileName <- gsub(paste0(".",ext), "", basename(metaFile$name))
      paste0(fileName,"_Cleaned.txt")
    },
    content = function(file) {
      meta <- metaReForm()
      write_delim(meta,file,delim = '\t')
    }
  )
  
  output$dnldParamOutFile <- downloadHandler(
    filename = function() {
      metaFile <- input$MetaFileInput
      ext <- tools::file_ext(metaFile$datapath)
      fileName <- gsub(paste0(".",ext), "", basename(metaFile$name))
      paste0(fileName,"_Parameters.txt")
    },
    content = function(file) {
      expr <- metaParam()
      write_delim(expr,file,delim = '\t',col_names = F)
    }
  )
  
  
  
  
  
}
# Run the application 
shinyApp(ui = ui, server = server)
