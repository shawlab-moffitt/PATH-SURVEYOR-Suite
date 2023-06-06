


####----File Input----####

GeneSetCat_File <- "GeneSet_Data/GeneSet_CatCountTable.zip"

GeneSet_File <- "GeneSet_Data/GeneSet_gsNsym.zip"

ExampleGeneList_File <- "Example_Input_File/PreRanked_GeneList.txt"

##--Advanced Set Up--##
# Optional preset input file
PreSet_RankedGeneList_File <- ""





####----Install and load packages----####

packages <- c("shiny","shinycssloaders","DT","readr")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))
#bioconductor packages
bioCpacks <- c("clusterProfiler","enrichplot")
installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
if (any(installed_packages_BIOC == FALSE)) {
  BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
}
invisible(lapply(bioCpacks, library, character.only = TRUE))




####----Read In and Sort File----####

GeneSetCats <- as.data.frame(readr::read_delim(GeneSetCat_File, delim = '\t', col_names = T))
GeneSetCats_tab <- GeneSetCats[,-4]
GeneSetCats_tab <- unique(GeneSetCats_tab)

GeneSet <- as.data.frame(readr::read_delim(GeneSet_File, delim = '\t', col_names = T))

if (is.character(PreSet_RankedGeneList_File)) {
  if (file.exists(PreSet_RankedGeneList_File)) {
    PreSet_RankedGeneList <- as.data.frame(readr::read_delim(PreSet_RankedGeneList_File, delim = '\t', col_names = T, comment = "#"))
  } else { PreSet_RankedGeneList <- NULL }
} else { PreSet_RankedGeneList <- NULL }


#increase file upload size
options(shiny.maxRequestSize=50*1024^2)


ui <- 
  
  #tagList(
    
    #shinyjs::useShinyjs(),
    #div(
    #  id = "loading_page",
    #  h1("Loading...")
    #),
    #hidden(
      #div(
      #  id = "main_content",
        
        navbarPage("{ Hazard Ratio Ranked GSEA }",
                   
                   ####----Enrichment Table Tab----####
                   
                   tabPanel("Enrichment Table",
                            fluidPage(
                              title = "Enrichment Table",
                              sidebarPanel(
                                width = 3,
                                fluidRow(
                                  column(9,
                                         fileInput("UserGeneList","Upload Ranked Gene List")
                                  ),
                                  column(3,
                                         checkboxInput("UserGeneListheaderCheck","Header",value = T)
                                  )
                                ),
                                fluidRow(
                                  column(6, style = 'margin-top:-10px;',
                                         uiOutput("rendUseExpData")
                                         #actionButton("UseExpData","Load Example Data")
                                         ),
                                  column(6, style = 'margin-top:-10px;',
                                         uiOutput("rendDownloadLink")
                                         #tags$a(href="http://shawlab.science/shiny/PATH_SURVEYOR_ExampleData/HazardRatio_PreRanked_GSEA_App/", "Download example data", target='_blank')
                                         )
                                ),
                                p(),
                                fluidRow(
                                  column(7,
                                         numericInput("gseaPval","GSEA Adj Pvalue Cutoff",value = 1)
                                  ),
                                  column(5,
                                         numericInput("UserSeedSet","Set Seed:",value = 101)
                                  )
                                ),
                                h4("Select a Gene Set:"),
                                div(DT::dataTableOutput("GenesetCatTable"), style = "font-size:10px"),
                                checkboxInput("UserUploadCheck","Upload Gene Set"),
                                fluidRow(
                                  column(9,
                                         uiOutput("rendUserGSupload")
                                  ),
                                  column(3,
                                         uiOutput("rendUserGSheaderCheck")
                                  )
                                ),
                                checkboxInput("PreviewGeneSet","Preview Selected Gene Set"),
                                uiOutput("rendGeneSetPreview"),
                              ),
                              mainPanel(
                                p("Please note this may take excess time depending on the number of gene sets being analyzed."),
                                uiOutput("rendEnrichTable"),
                                uiOutput("renddnldEnrichTable"),
                              )
                            )
                   ),
                   
                   ####----Enrichment Plot Tab----####
                   
                   tabPanel("Enrichment Plot",
                            fluidPage(
                              title = "Enrichment Plot",
                              sidebarPanel(
                                h4("Select Gene Set:"),
                                div(DT::dataTableOutput("GeneSetTable"), style = "font-size:10px")
                              ),
                              mainPanel(
                                h3("GSEA Enrichment Plot"),
                                verbatimTextOutput("NESandPval"),
                                shinycssloaders::withSpinner(plotOutput("enrichplot", width = "500px", height = "450px"), type = 6),
                                fluidRow(
                                  downloadButton("dnldPlotSVG_gsea","Download as SVG"),
                                  downloadButton("dnldPlotPDF_gsea","Download as PDF")
                                ),
                                h3("Leading Edge Genes"),
                                downloadButton("LEGdownload", "Download Leading Edge Gene Table"),
                                div(DT::dataTableOutput("LeadingEdgeGenes"), style = "font-size:12px; height:500px; width:400px")
                              )
                            )
                   )
        )
      #)
    #)
  #)


server <- function(input, output, session) {
  
  #load_data()
  
  ####----Render UI----####
  
  
  output$rendUseExpData <- renderUI({
    
    if (is.null(PreSet_RankedGeneList)) {
      actionButton("UseExpData","Load Example Data")
    }
    
  })
  
  output$rendDownloadLink <- renderUI({
    
    if (is.null(PreSet_RankedGeneList)) {
      tags$a(href="http://shawlab.science/shiny/PATH_SURVEYOR_ExampleData/HazardRatio_PreRanked_GSEA_App/", "Download example data", target='_blank')
    }
    
  })
  
  #output$rendPreSetDataUI <- renderUI({
  #  
  #  if (is.null(PreSet_RankedGeneList)) {
  #    fluidRow(
  #      column(6, style = 'margin-top:-10px;',
  #             #uiOutput("rendUseExpData")
  #             actionButton("UseExpData","Load Example Data")
  #      ),
  #      column(6, style = 'margin-top:-10px;',
  #             #uiOutput("rendDownloadLink")
  #             tags$a(href="http://shawlab.science/shiny/PATH_SURVEYOR_ExampleData/HazardRatio_PreRanked_GSEA_App/", "Download example data", target='_blank')
  #      )
  #    )
  #    p()
  #  }
  #  
  #})
  
  
  output$rendEnrichTable <- renderUI({
    
    #req(input$UserGeneList)
    req(isTruthy(input$UserGeneList) | isTruthy(input$UseExpData) | isTruthy(PreSet_RankedGeneList))
    withSpinner(DT::dataTableOutput("EnrichTable", width = "100%"), type = 6)
    
  })
  
  output$renddnldEnrichTable <- renderUI({
    
    #req(input$UserGeneList)
    req(isTruthy(input$UserGeneList) | isTruthy(input$UseExpData) | isTruthy(PreSet_RankedGeneList))
    downloadButton("dnldEnrichTable","Download Enriched Signatures Table")
    
  })
  
  output$rendUserGSupload <- renderUI({
    
    #if (is.null(input$GeneSetChosen) == FALSE) {
      
      if (input$UserUploadCheck == TRUE) {
        
        fileInput("UserGeneSet","Upload Gene Set File")
        
      }
      
    #}
    
  })
  
  output$rendUserGSheaderCheck <- renderUI({
    
    #if (is.null(input$GeneSetChosen) == FALSE) {
      
      if (input$UserUploadCheck == TRUE) {
        
        checkboxInput("UserGSheaderCheck","Header",value = T)
        
      }
      
    #}
    
  })
  
  ####----Reactives----####
  
  user_gl_df <- reactiveVal()
  user_gl_df_fileName <- reactiveVal()
  
  ## User Upload gene set reactive
  user_gs <- reactive({
    
    header_check <- input$UserGSheaderCheck
    gs.u <- input$UserGeneSet
    ext <- tools::file_ext(gs.u$datapath)
    req(gs.u)
    validate(need(ext == c("gmt","tsv","txt"), "Please upload .gmt, .tsv, or .txt file"))
    
    # If user provides GMT file
    if (ext == "gmt") {
      gmt <- clusterProfiler::read.gmt(gs.u$datapath)
      colnames(gmt) <- c("term","gene")
      gs_u <- gmt
    }
    
    # If user provides tab-delim file
    else if (ext == "txt" | ext == "tsv") {
      
      gmt <- as.data.frame(readr::read_delim(gs.u$datapath, delim = '\t', col_names = header_check))
      
      if (ncol(gmt) == 1) {
        paths <- as.vector(gmt[,1])
        paths <- gsub("[[:punct:]]",".",paths)
        gs_u <- GeneSet[which(GeneSet$term %in% paths),]
      }
      else if (ncol(gmt) == 2) {
        colnames(gmt) <- c("term","gene")
        gs_u <- gmt
      }
      else if (ncol(gmt) > 2) {
        paths <- as.vector(gmt[,1])
        paths <- gsub("[[:punct:]]",".",paths)
        gs_u <- GeneSet[which(GeneSet$term %in% paths),]
      }
      
    }
    
    gs_u
    
  })
  
  #gmt reactive based on user choice for GSEA
  gmt_react <- reactive({
    
    if (input$UserUploadCheck == TRUE) {
      
      #req(input$UserGeneList)
      gmt <- user_gs()
      
    } else {
      
      gs_chosen <- as.character(GeneSetCats_tab[input$GenesetCatTable_rows_selected,3])
      
      GeneSets <- GeneSetCats[which(GeneSetCats$GeneSet_Sub_Category == gs_chosen),"GeneSet_Name"]
      gmt <- GeneSet[which(GeneSet[,1] %in% GeneSets),]
      
    }
    
  })
  
  output$rendGeneSetPreview <- renderUI({
    
    if (input$PreviewGeneSet == TRUE) {
      
      div(DT::dataTableOutput("GeneSetPreview"), style = "font-size:10px; height:450px; overflow-X: scroll")
      
    }
    
  })
  
  output$GeneSetPreview <- DT::renderDataTable({
    
    #if (input$PreviewGeneSet == TRUE) {
      
      gmt <- gmt_react()
      geneset_names <- data.frame("Gene_Set_Name" = unique(gmt[,1]))
      DT::datatable(geneset_names,
                    options = list(keys = TRUE,
                                   searchHighlight = TRUE,
                                   pageLength = 10,
                                   lengthMenu = c("10", "20", "50", "100")),
                    rownames = F)
      
    #}
    
    
  })
  
  observeEvent(input$UseExpData, {
    
    df <- as.data.frame(readr::read_delim(ExampleGeneList_File, delim = '\t', col_names = T, comment = "#"))
    user_gl_df(df)
    user_gl_df_fileName(ExampleGeneList_File)
    
  })
  
  observe({
    
    if (!is.null(PreSet_RankedGeneList)) {
      df <- PreSet_RankedGeneList
      user_gl_df(df)
    } else {
      header_check <- input$UserGeneListheaderCheck
      gs.u <- input$UserGeneList
      ext <- tools::file_ext(gs.u$datapath)
      req(gs.u)
      validate(need(ext == c("tsv","txt"), "Please upload .tsv or .txt file"))
      
      df <- as.data.frame(readr::read_delim(gs.u$datapath, delim = '\t', col_names = header_check, comment = "#"))
      user_gl_df(df)
    }
    
  })
  
  #generate ranked list of genes
  user_gl <- reactive({
    
    req(isTruthy(input$UserGeneList) | isTruthy(input$UseExpData) | isTruthy(PreSet_RankedGeneList))
    df <- user_gl_df()
    
    # assume just list of gene symbols
    if (ncol(df) == 1) {
      
      # generate rank from high to low
      ranking <- seq(from = nrow(df),to = 1)
      names(ranking) <- as.character(df[,1])
      
    }
    # assume gene symbols and rank
    else if (ncol(df) == 2) {
      
      ranking <- df[,2]
      names(ranking) <- as.character(df[,1])
      ranking <- sort(ranking, decreasing = T)
      
    }
    # assume gene symbols and general annotation columns
    else if (ncol(df) > 2) {
      
      # assume Coxh Ranking table is output
      if (colnames(df)[2] == "Hazard_Ratio" | colnames(df)[2] == "estimate") {
        
        ranking <- df[,2]
        names(ranking) <- as.character(df[,1])
        ranking <- sort(ranking, decreasing = T)
        
      }
      else {
        
        # generate rank from high to low
        ranking <- seq(from = nrow(df),to = 1)
        names(ranking) <- as.character(df[,1])
        
      }
      
    }
    
    ranking
    
  })
  
  EnrichTabReact <- reactive({
    
    req(isTruthy(input$UserGeneList) | isTruthy(input$UseExpData) | isTruthy(PreSet_RankedGeneList))
    set.seed(input$UserSeedSet)
    gmt <- gmt_react()
    ranked_genes <- user_gl()
    head(ranked_genes)
    pval <- input$gseaPval
    gsea.res <- clusterProfiler::GSEA(ranked_genes,
                     TERM2GENE = gmt,
                     eps = NA, pvalueCutoff = pval,
                     minGSSize = 5, maxGSSize = 10000,
                     nPermSimple = 5000,
                     seed = T, verbose = T)
    gsea.res
    
  })
  
  
  LeadingEdgeReact <- reactive({
    
    req(isTruthy(input$UserGeneList) | isTruthy(input$UseExpData) | isTruthy(PreSet_RankedGeneList))
    res <- EnrichTabReact()
    res_df <- res@result
    geneset_name <- as.character(res_df[input$GeneSetTable_rows_selected,1])
    NES = res_df$NES[which(res_df[,'ID']==geneset_name)]
    geneList <- user_gl()
    ## Subset core enriched genes
    genes1 <- as.matrix(res_df[which(res_df$ID==geneset_name),"core_enrichment"])
    genes2 <- strsplit(genes1,"/")
    if (length(genes2) > 0) {
      GeneSymbol <- as.data.frame(genes2)
      colnames(GeneSymbol)[1] <- "Gene_Symbol"
      GeneSymbol$Leading_Edge_Rank <- as.numeric(rownames(GeneSymbol))
      GeneSymbol <- GeneSymbol[,c("Leading_Edge_Rank","Gene_Symbol")]
      
      ranking_df <- stack(geneList)
      uploaded_df <- user_gl_df()
      if (input$UserGeneListheaderCheck == TRUE) {
        if (ncol(uploaded_df) == 2) {
          rank_col_name <- colnames(uploaded_df)[2]
        }
        if (ncol(uploaded_df) == 1) {
          rank_col_name <- "Overall_Rank"
        }
        if (ncol(uploaded_df) > 2) {
          if (colnames(uploaded_df)[2] == "Hazard_Ratio" | colnames(uploaded_df)[2] == "estimate") {
            rank_col_name <- "Hazard_Ratio"
          }
          else {
            rank_col_name <- "Overall_Rank"
          }
        }
      }
      if (input$UserGeneListheaderCheck == FALSE) {
        if (ncol(uploaded_df) == 2) {
          rank_col_name <- "User_Uploaded_Ranking"
        }
        if (ncol(uploaded_df) != 2) {
          rank_col_name <- "Overall_Rank"
        }
      }
      
      colnames(ranking_df) <- c(rank_col_name,"Gene_Symbol")
      
      Leading_Merge <- merge(GeneSymbol,ranking_df,by = "Gene_Symbol",all.x = T)
      if (NES < 0) {
        Leading_Merge <- Leading_Merge[order(Leading_Merge[,1]),]
      }
      else if (NES > 0) {
        Leading_Merge <- Leading_Merge[order(Leading_Merge[,1], decreasing = F),]
      }
      Leading_Merge <- Leading_Merge[,c("Leading_Edge_Rank","Gene_Symbol",rank_col_name)]
      Leading_Merge
    }
    
  })
  
  ####----Data Tables----####
  
  output$GenesetCatTable <- DT::renderDataTable({
    
    colnames(GeneSetCats_tab) <- gsub("_"," ",colnames(GeneSetCats_tab))
    DT::datatable(GeneSetCats_tab,
                  extensions = c("KeyTable"),
                  selection = list(mode = 'single', selected = 1),
                  options = list(keys = T,
                                 searchHighlight = T,
                                 pageLength = 5,
                                 lengthMenu = c("5","10", "25", "50", "100"),
                                 scrollX = T),
                  rownames = F)
    
  })
  
  output$EnrichTable <- DT::renderDataTable({
    
    req(isTruthy(input$UserGeneList) | isTruthy(input$UseExpData) | isTruthy(PreSet_RankedGeneList))
    res <- EnrichTabReact()
    res_df <- res@result
    
    if (nrow(res@result) == 0) {
      
      res_df <- data.frame("No_term_enriched_under_specific_Pvalue_Cutoff")
      
    }
    
    DT::datatable(res_df,
                  extensions = c("KeyTable", "FixedHeader","FixedColumns"),
                  options = list(keys = T,
                                 searchHighlight = T,
                                 pageLength = 25,
                                 lengthMenu = c("10", "25", "50", "100"),
                                 scrollX = T,
                                 autoWidth = TRUE,
                                 fixedColumns = list(leftColumns = 1)),
                  rownames = F) %>%
      formatRound(columns = c(2:10), digits = 4)
    
  })
  
  output$GeneSetTable <- DT::renderDataTable({
    
    res <- EnrichTabReact()
    res_df <- res@result
    GeneSets <- res_df[,c(1,5,6)]
    
    DT::datatable(GeneSets,
                  selection = list(mode = 'single', selected = 1),
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 20,
                                 lengthMenu = c("10", "20", "50", "100")),
                  rownames = F) %>%
      formatRound(columns = c(2,3), digits = 4)
    
  })
  
  output$LeadingEdgeGenes <- DT::renderDataTable({
    
    req(isTruthy(input$UserGeneList) | isTruthy(input$UseExpData) | isTruthy(PreSet_RankedGeneList))
    LEGs <- LeadingEdgeReact()
    DT::datatable(LEGs,
                  options = list(paging = F,
                                 autoWidth = TRUE,
                                 columnDefs = list(list(className = 'dt-center', targets = 0:2))),
                  rownames = F)
    
  })
  
  ####----Enrichment Plot----####
  
  output$enrichplot <- renderPlot({
    
    req(isTruthy(input$UserGeneList) | isTruthy(input$UseExpData) | isTruthy(PreSet_RankedGeneList))
    if (length(input$GeneSetTable_rows_selected) > 0) {
      
      res <- EnrichTabReact()
      res_df <- res@result
      geneset_name <- as.character(res_df[input$GeneSetTable_rows_selected,1])
      enrichplot::gseaplot2(res,
                geneset_name,
                geneset_name,
                pvalue_table = F)
      
    }
    
  })
  
  output$NESandPval <- renderText({
    
    req(isTruthy(input$UserGeneList) | isTruthy(input$UseExpData) | isTruthy(PreSet_RankedGeneList))
    if (length(input$GeneSetTable_rows_selected) > 0){
      res <- EnrichTabReact()
      gsea.df <- as.data.frame(res@result)
      GS = as.character(gsea.df[input$GeneSetTable_rows_selected,1])
      NES = gsea.df$NES[which(gsea.df[,'ID']==GS)]
      Pval = gsea.df$pvalue[which(gsea.df[,'ID']==GS)]
      NES.o <- paste0("NES: ", NES)
      Pval.o <- paste0("Pvalue: ", Pval)
      paste(NES.o, Pval.o, sep = '\n')
    }
    else if (length(input$msigdbTable_rows_selected) == 0){
      paste("Please select gene set from side panel table to begin.", sep = '')
    }
    
  })
  
  ####----Downloads----####
  
  output$dnldEnrichTable <- downloadHandler(
    filename = function() {
      paste("EnrichedSignaturesTable.txt",sep = "")
    },
    content = function(file) {
      res <- EnrichTabReact()
      res_df <- res@result
      write_delim(res_df,file,delim = '\t')
    }
  )
  
  output$dnldPlotSVG_gsea <- downloadHandler(
    filename = function() {
      
      res <- EnrichTabReact()
      gsea.df <- as.data.frame(res@result)
      GS = as.character(gsea.df[input$GeneSetTable_rows_selected,1])
      paste(GS,"_EnrichmentPlot.svg",sep = "")
      
    },
    content = function(file) {
      
      res <- EnrichTabReact()
      res_df <- res@result
      geneset_name <- as.character(res_df[input$GeneSetTable_rows_selected,1])
      p <- gseaplot2(res,
                     geneset_name,
                     geneset_name,
                     pvalue_table = F)
      ggsave(file,p, width = 10, height = 8)
      
    }
  )
  output$dnldPlotPDF_gsea <- downloadHandler(
    filename = function() {
      
      res <- EnrichTabReact()
      gsea.df <- as.data.frame(res@result)
      GS = as.character(gsea.df[input$GeneSetTable_rows_selected,1])
      paste(GS,"_EnrichmentPlot.pdf",sep = "")
      
    },
    content = function(file) {
      
      res <- EnrichTabReact()
      res_df <- res@result
      geneset_name <- as.character(res_df[input$GeneSetTable_rows_selected,1])
      p <- gseaplot2(res,
                     geneset_name,
                     geneset_name,
                     pvalue_table = F)
      ggsave(file,p, width = 10, height = 8)
      
    }
  )
  
  output$LEGdownload <- downloadHandler(
    
    filename = function() {
      
      res <- EnrichTabReact()
      res_df <- res@result
      geneset_name <- as.character(res_df[input$GeneSetTable_rows_selected,1])
      paste(geneset_name,"_LeadingEdgeGenes.txt",sep ="")
      
    },
    content = function(file) {
      
      LEGs <- LeadingEdgeReact()
      write_delim(LEGs,file,delim = '\t')
      
    }
    
  )
  
}



shinyApp(ui,server)












