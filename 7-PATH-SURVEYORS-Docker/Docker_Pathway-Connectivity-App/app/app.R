####----Install and load packages----####

packages <- c("shiny","shinythemes","shinyjqui","pheatmap","RColorBrewer",
              "ggdendro","factoextra","dplyr","DT","viridis","readr",
              "shinycssloaders","stringr","tools","plotly","reshape2")

#installed_packages <- packages %in% rownames(installed.packages())
#if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages])
#}
invisible(lapply(packages, library, character.only = TRUE))
#bioconductor packages
bioCpacks <- c("clusterProfiler")
#installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
#if (any(installed_packages_BIOC == FALSE)) {
#  BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
#}
invisible(lapply(bioCpacks, library, character.only = TRUE))


####---- Read in Files----####

GeneSet_File <- "GeneSet_Data/Comprehensive_GeneSet.RData"






# R Data list load function for naming
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
gs <- loadRData(GeneSet_File)

#learn Jaccard index function
jaccard <- function(a, b) {  
  a = unique(a)
  b = unique(b)
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (1-(intersection/union))
}

#increase file upload size
options(shiny.maxRequestSize=50*1024^2)

shinytheme("sandstone")

ui <- 
  navbarPage("{ Jaccard Pathway Connectivity Index }",
             
             ####----Intra-Pathway Connectivity----####
             
             tabPanel("Pathway Connectivity",
                      fluidPage(
                        title = "Pathway Connectivity",
                        sidebarLayout(
                          sidebarPanel(
                            tabsetPanel(
                              
                              ####----Pathway Input----####
                              
                              tabPanel("Pathway Parameters",
                                       p(),
                                       h4("Upload Pathways of Interest"),
                                       fluidRow(
                                         column(9,
                                                 fileInput("UserPathwayFile","Upload File (.gmt/.tsv/.txt/.RData)",
                                                           accept = c(".gmt",".tsv",".txt"))
                                                ),
                                         column(3,
                                                checkboxInput("HeaderCheckIntra","Header",value = T))
                                       ),
                                       h4("Clustering Parameters"),
                                       fluidRow(
                                         column(4,
                                                selectInput("ClustMethodIntra","Clustering Method:",
                                                            choices = c("ward.D", "complete", "ward.D2", "single", "average", "mcquitty", "median", "centroid"))
                                                ),
                                         column(8,
                                                numericInput("NumClusters", step = 1, label = "Number of Clusters (Cut Tree with ~k)", value = 10)
                                                )
                                         ),
                                       checkboxInput("ViewClustTabIntra","View Cluster Results Table"),
                                       uiOutput("rendClustTabIntra"),
                                       downloadButton("dnldClustTabIntra","Download Cluster Result"),
                                       h4("SIF Download"),
                                       numericInput("JaccDistCutoff","Jaccard Distance Cutoff",
                                                    min = 0,max = 1, step = 0.1, value = 0.9, width = "200px"),
                                       checkboxInput("PrevSIF","Preview SIF File",value = F),
                                       uiOutput("rendSIFPreview"),
                                       downloadButton("dnldSIFTabIntra","Download SIF File")
                                       ),
                              
                              ####----Figure Parameters----####
                              
                              tabPanel("Figure Parameters",
                                       h4("Heatmap Parameters"),
                                       selectInput("ColorPaletteIntra", "Select Color Palette:",
                                                   choices = c("Red/Blue" = "original",
                                                               "OmniBlueRed" = "OmniBlueRed",
                                                               "LightBlue/BlackRed" = "LightBlueBlackRed",
                                                               "Green/Black/Red" = "GreenBlackRed",
                                                               "Yellow/Green/Blue" = "YlGnBu","Inferno" = "Inferno",
                                                               "Viridis" = "Viridis","Plasma" = "Plasma",
                                                               "Reds" = "OrRd","Blues" = "PuBu","Greens" = "Greens")),
                                       fluidRow(
                                         column(6,
                                                checkboxInput("HeatColNamesIntra","Show Heatmap Column Names", value = F),
                                                numericInput("HeatColFontIntra", "Heatmap Column Font Size:",
                                                             min = 5, max = 75,
                                                             value = 12, step = 1),
                                                numericInput("HeatColDendHeight","Column Dendrogram Height:",
                                                             value = 50, step = 1)
                                                ),
                                         column(6,
                                                checkboxInput("HeatRowNamesIntra","Show Heatmap Row Names", value = F),
                                                numericInput("HeatRowFontIntra", "Heatmap Row Font Size:",
                                                             min = 5, max = 75,
                                                             value = 9, step = 1),
                                                numericInput("HeatRowDendHeight","Row Dendrogram Height:",
                                                             value = 50, step = 1)
                                                )
                                         ),
                                       h4("Connectivity Visualization Parameters"),
                                       selectInput("ConnView","View Connectivity as:",
                                                   choices = c("Phylogeny" = "phylogenic","Dendrogram" = "rectangle","Circular" = "circular")),
                                       uiOutput("rendPhyloLayout"),
                                       fluidRow(
                                         column(6,
                                                checkboxInput("ShowConnLabels","Show Labels:",value = T)
                                                ),
                                         column(6,
                                                numericInput("ConnFontSize","Font Size:",
                                                             value = 0.4, step = 0.1)
                                                )
                                         )
                                       )
                              )
                            ),
                          mainPanel(
                            tabsetPanel(
                              
                              ####----Jaccard Table----####
                              
                              tabPanel("Jaccard Pathway Connectivity Table",
                                       p(),
                                       div(DT::dataTableOutput("JaccTableIntra"), style = "font-size:12px"),
                                       uiOutput("renddnldJaccTabIntra")
                                       ),
                              
                              ####----Jaccard Heatmap----####
                              
                              tabPanel("Heatmap",
                                       withSpinner(jqui_resizable(plotOutput('JaccHeatmapIntra', width = "100%", height = "800px")), type = 6)
                                       ),
                              
                              ####----Jaccard Clustering----####
                              
                              tabPanel("Clustering",
                                       uiOutput("rendJaccDendo")
                                       ),
                              
                              ####----Jaccard Cluster Annotation----####
                              
                              tabPanel("Gene Clusters and Annoation",
                                       fileInput("UserAnnotationFile","Upload Annotation File",
                                                 accept = c(".tsv",".txt",".csv")),
                                       div(DT::dataTableOutput("ClusterTabAnno"), style = "font-size:12px"),
                                       uiOutput("renddnldClusterTabAnno")
                                       )
                              )
                            )
                          )
                        )
             )
  )



server <- function(input, output, session) {
  
  ####----Render UI----####
  
  output$rendClustTabIntra <- renderUI({
    
    if (input$ViewClustTabIntra == TRUE) {
      
      div(DT::dataTableOutput("ClustTabIntra"), style = "font-size:10px; height:400px; overflow-X: scroll")
      
    }
    
  })
  
  output$rendSIFPreview <- renderUI({
    
    if (input$PrevSIF == TRUE) {
      
      div(DT::dataTableOutput("SIFPreview"), style = "font-size:10px; height:400px; overflow-X: scroll")
      
    }
    
  })
  
  output$rendSIFPreviewINTER <- renderUI({
    
    if (input$PrevSIFINTER == TRUE) {
      
      div(DT::dataTableOutput("SIFPreviewINTER"), style = "font-size:10px; height:400px; overflow-X: scroll")
      
    }
    
  })
  
  output$rendClustTabInter <- renderUI({
    
    if (input$ViewClustTabInter == TRUE) {
      
      div(DT::dataTableOutput("ClustTabInter"), style = "font-size:10px; height:400px; overflow-X: scroll")
      
    }
    
  })
  
  output$renddnldJaccTabIntra <- renderUI({
    
    req(input$UserPathwayFile)
    downloadButton("dnldJaccTabIntra","Download Jaccard Connectivity Table")
    
  })
  
  output$renddnldJaccTabInter <- renderUI({
    
    req(input$UserPathwayFile)
    downloadButton("dnldJaccTabInter","Download Jaccard Connectivity Table")
    
  })
  
  output$rendJaccDendo <- renderUI({
    
    VisType <- input$ConnView
    if (VisType == "rectangle" | VisType == "phylogenic") {
      
      withSpinner(jqui_resizable(plotlyOutput('JaccDendoIntra', width = "100%", height = "900px")), type = 6)
      
    }
    else if (VisType == "circular") {
      
      withSpinner(jqui_resizable(plotOutput('JaccDendoIntraCirc', width = "100%", height = "900px")), type = 6)
      
    }
    
  })
  
  output$rendPhyloLayout <- renderUI({
    
    if (input$ConnView == "phylogenic") {
      
      selectInput("PhyloLayout","Choose Phylogenic Tree Layout:",
                  choices = c("Auto Layout" = "layout.auto",
                              "Force Directed (DrL) Layout" = "layout_with_drl",
                              "Tree Layout" = "layout_as_tree",
                              "GEM Layout" = "layout.gem",
                              "Mulitdimensional Scaling Layout" = "layout.mds",
                              "Large Graph Layout" = "layout_with_lgl"))
      
    }
    
  })
  
  ####----Reactives----####
  
  ## User Upload gene set list reactive
  user_gs <- reactive({
    
    header_check <- input$HeaderCheckIntra
    gs.u <- input$UserPathwayFile
    ext <- tools::file_ext(gs.u$datapath)
    req(gs.u)
    #validate(need(ext == c("gmt","tsv","txt", "RData"), "Please upload .gmt, .tsv, .txt, or .RData file"))
    validate(need(ext == c("tsv","txt"), "Please upload .tsv or .txt file"))
    
    
    gmt <- as.data.frame(read_delim(gs.u$datapath, delim = '\t', col_names = header_check))
    
    paths <- as.vector(gmt[,1])
    paths <- gsub("[[:punct:]]",".",paths)
    gs_u <- gs[unique(paths)]
    gs_u <- Filter(Negate(is.null), gs_u)
    
    ## If user provides GMT file
    #if (ext == "gmt") {
    #  gmt <- read.gmt(gs.u$datapath)
    #  colnames(gmt) <- c("term","gene")
    #  gs_u <- list()
    #  for (i in unique(gmt[,1])){
    #    gs_u[[i]] <- gmt[gmt[,1] == i,]$gene
    #  }
    #}
    #
    ## If user provides RData list file
    #else if (ext == "RData") {
    #  gs_u <- loadRData(gs.u$datapath)
    #}
    #
    ## If user provides tab-delim two-col file
    #else if (ext == "txt" | ext == "tsv") {
    #  
    #  gmt <- as.data.frame(read_delim(gs.u$datapath, delim = '\t', col_names = header_check))
    #  
    #  if (ncol(gmt) == 1) {
    #    
    #    paths <- as.vector(gmt[,1])
    #    paths <- gsub("[[:punct:]]",".",paths)
    #    gs_u <- gs[unique(paths)]
    #    
    #  }
    #  else if (ncol(gmt) == 2) {
    #    
    #    colnames(gmt) <- c("term","gene")
    #    gs_u <- list()
    #    for (i in unique(gmt[,1])){
    #      gs_u[[i]] <- gmt[gmt[,1] == i,]$gene
    #    }
    #    
    #  }
    #  else if (ncol(gmt) > 2) {
    #    
    #    paths <- as.vector(gmt[,1])
    #    paths <- gsub("[[:punct:]]",".",paths)
    #    gs_u <- gs[unique(paths)]
    #    
    #  }
    #  
    #}
    
    gs_u
    
  })
  
  user_anno <- reactive({
    
    gs.u <- input$UserAnnotationFile
    ext <- tools::file_ext(gs.u$datapath)
    req(gs.u)
    validate(need(ext == c("tsv","txt","csv"), "Please upload .tsv, .txt, or .csv file"))
    if (ext == "txt" | ext == "tsv") {
      
      df <- as.data.frame(read_delim(gs.u$datapath, delim = '\t', col_names = T, comment = "##"))
      
    }
    else if (ext == "csv") {
      
      df <- as.data.frame(read_delim(gs.u$datapath, delim = ',', col_names = T))
      
    }
    
    df
    
  })
  
  ## New Subset data frame based on gene set terms searched
  GeneSetTableSearched <- reactive({
    
    ## Gather Search Terms
    gs_search <- input$GeneSetSerachTerm
    gs_search <- unlist(strsplit(gs_search," "))
    gs_search <- unlist(strsplit(gs_search,"[[:punct:]]"))
    
    ## Subset Gene Set Table for searched terms
    GeneSetTable_searched <- GeneSetTable[apply(GeneSetTable, 1, function(x) any(grepl(paste(gs_search, collapse = "|"), x, ignore.case = T))), ]
    
    GeneSetTable_searched
    
  })
  
  ## Gene Set List Selected
  our_gs <- reactive({
    
    if (input$GeneSetSelection == 1) {
      
      ## Get gene set category name
      gs_category <- GeneSetCatTable[input$GeneSetCategoryTable_rows_selected,2]
      
      ## Get all gene set names from that category
      gs_names <- GeneSetTable[which(GeneSetTable$GeneSet_Sub_Category == gs_category),"GeneSet_Name"]
      
      ## Subset gene set list for specified gene sets
      gs_o <- gs[gs_names]
      gs_o
      
    }
    else if (input$GeneSetSelection == 2) {
      
      if (input$RunSearchedForGeneSets == TRUE) {
        
        ## Reactive gene set table based on search terms
        GeneSetTable_sub <- GeneSetTableSearched()
        
        ## Get gene set names from subset table
        gs_names <- GeneSetTable_sub$GeneSet_Name
        
        ## Subset gene set list for specified gene sets
        gs_o <- gs[gs_names]
        gs_o
        
      }
      
    }
    ## Compare with all gene sets
    else if (input$GeneSetSelection == 3) {
      
      if (input$RunAllGeneSets == TRUE) {
        
        gs_o <- gs
        gs_o
        
      }
      
    }
    
  })
  
  ## Jaccard Matrix Generation
  jacc_react_Intra <- reactive({
    
    if (length(user_gs()) > 0) {
      
      gs1 <- user_gs() #User Gene Sets Input
      gs2 <- user_gs() #User Gene Sets Input
      
      jac_list <- list()
      for (i in 1:length(gs1)) {
        
        jac_vect <- c()
        
        for (j in 1:length(gs2)) {
          
          gs1_genes <- gs1[[i]]
          gs2_genes <- gs2[[j]]
          jac_index <- jaccard(gs1_genes,gs2_genes)
          jac_vect <- c(jac_vect,jac_index)
          
        }
        
        jac_list[[paste(names(gs1)[i],sep = "")]] <- jac_vect
        
      }
      
      jac_df_Intra <- do.call(rbind,jac_list)
      colnames(jac_df_Intra) <- names(gs2)
      
      jac_df_Intra
      
    }
    
  })
  
  ## Jaccard Matrix Generation
  jacc_react_Inter <- reactive({
    
    if (length(our_gs()) > 0 & length(user_gs()) > 0) {
      
      gs1 <- our_gs()  #Refernce Gene Sets Selected
      gs2 <- user_gs() #User Gene Sets Input
      
      jac_list <- list()
      for (i in 1:length(gs1)) {
        
        jac_vect <- c()
        
        for (j in 1:length(gs2)) {
          
          gs1_genes <- gs1[[i]]
          gs2_genes <- gs2[[j]]
          jac_index <- jaccard(gs1_genes,gs2_genes)
          jac_vect <- c(jac_vect,jac_index)
          
        }
        
        jac_list[[paste(names(gs1)[i],sep = "")]] <- jac_vect
        
      }
      
      jac_df_Inter <- do.call(rbind,jac_list)
      colnames(jac_df_Inter) <- names(gs2)
      
      jac_df_Inter
      
    }
    
  })
  
  ####----Data Tables----####
  
  ## Gene Set Category table for tab 1
  output$GeneSetCategoryTable <- DT::renderDataTable({
    
    DT::datatable(GeneSetCatTable,
                  selection = list(mode = 'single', selected = 1),
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")),
                  rownames = F)
    
  })
  
  ## Gene Set Search Table for tab 2
  output$GeneSetSearchedTable <- DT::renderDataTable({
    
    GeneSetTable_searched <- GeneSetTableSearched()
    
    DT::datatable(GeneSetTable_searched,
                  options = list(dom = 'lrtp',
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")),
                  rownames = F)
    
  })
  
  ## Jaccard Matrix Table - INTRA pathway
  output$JaccTableIntra <- DT::renderDataTable({
    
    ## Jaccard Table
    jacc_df <- as.data.frame(jacc_react_Intra())
    colnames(jacc_df) <- gsub("[[:punct:]]", " ",colnames(jacc_df))
    
    DT::datatable(jacc_df,
                  extensions = "FixedColumns",
                  options = list(lengthMenu = c(10,20,50,100,1000,5000,10000),
                                 pageLength = 20,
                                 scrollX = TRUE,
                                 autoWidth = TRUE,
                                 fixedColumns = list(leftColumns = 1)),
                  selection=list(mode = "multiple")) %>%
      formatRound(columns = c(1:ncol(jacc_df)), digits = 4)
    
    
  })
  
  ## Jaccard Matrix Table - Inter Pathway
  output$JaccTableInter <- DT::renderDataTable({
    
    ## Jaccard Table
    jacc_df <- as.data.frame(jacc_react_Inter())
    colnames(jacc_df) <- gsub("[[:punct:]]", " ",colnames(jacc_df))
    
    DT::datatable(jacc_df,
                  extensions = "FixedColumns",
                  options = list(lengthMenu = c(10,20,50,100,1000,5000,10000),
                                 pageLength = 20,
                                 scrollX = TRUE,
                                 autoWidth = TRUE,
                                 fixedColumns = list(leftColumns = 1)),
                  selection=list(mode = "multiple")) %>%
      formatRound(columns = c(1:ncol(jacc_df)), digits = 4)
    
    
  })
  
  ## Cluster Table
  output$ClustTabIntra <- DT::renderDataTable({
    
    
    ## Data matrix
    jacc_df <- jacc_react_Intra()
    jacc_mat <- as.matrix(jacc_df)
    clust_me <- input$ClustMethodIntra
    
    hm_res <- pheatmap::pheatmap(jacc_mat,
                                 clustering_method = clust_me,
                                 silent = T)
    
    jacc_df <- jacc_react_Intra()
    cut_k <- input$NumClusters
    clustTab <- data.frame(Pathways = rownames(jacc_df),
                           Cluster = cutree(hm_res$tree_row,k = cut_k))
    
    DT::datatable(clustTab,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")),
                  rownames = F)
    
  })
  
  output$ClusterTabAnno <- DT::renderDataTable({
    
    if (is.null(input$UserAnnotationFile) == F) {
      
      ## user annotation table
      anno_df <- user_anno()
      
    }
    
    
    ## user gs list
    gs_u <- user_gs()
    
    ## cluster table
    ## Data matrix
    jacc_df <- jacc_react_Intra()
    jacc_mat <- as.matrix(jacc_df)
    clust_me <- input$ClustMethodIntra
    hm_res <- pheatmap::pheatmap(jacc_mat,
                                 clustering_method = clust_me,
                                 silent = T)
    jacc_df <- jacc_react_Intra()
    cut_k <- input$NumClusters
    clustTab <- data.frame(Pathways = rownames(jacc_df),
                           Cluster = cutree(hm_res$tree_row,k = cut_k))
    
    ## table of gene set with genes
    nameVector <- unlist(mapply(function(x,y){ rep(y, length(x)) }, gs_u, names(gs_u)))
    resultDF <- cbind.data.frame(Pathways = nameVector,Genes = unlist(gs_u))
    rownames(resultDF) <- 1:nrow(resultDF)
    
    ## double check punct matches
    clustTab$Pathways <- gsub("[[:punct:]]",".",clustTab$Pathways)
    resultDF$Pathways <- gsub("[[:punct:]]",".",resultDF$Pathways)
    
    ## Merge geneset and cluster table
    df2 <- merge(resultDF,clustTab, by = "Pathways")
    
    if (is.null(input$UserAnnotationFile) == F) {
      
      colnames(anno_df)[1] <- "Genes"
      df3 <- merge(df2,anno_df,by = "Genes", all.x = T)
      
    }
    else if (is.null(input$UserAnnotationFile) == T) {
      
      df3 <- df2
      
    }
    
    
    df4 <- df3 %>%
      relocate(Pathways,Cluster,Genes)
    
    df5 <- df4[order(df4$Pathways),]
    
    DT::datatable(df5,
                  options = list(lengthMenu = c(10,20,50,100,1000,5000,10000),
                                 pageLength = 20,
                                 scrollX = TRUE),
                  rownames = F,
                  selection=list(mode = "multiple"))
    
    
  })
  
  ## Cluster Table
  output$ClustTabInter <- DT::renderDataTable({
    
    hm_res <- JaccHeatInter_react()
    jacc_df <- jacc_react_Inter()
    cut_k <- input$NumClustersInter
    
    clustTab <- data.frame(Pathways = rownames(jacc_df),
                           Cluster = cutree(hm_res$tree_row,k = cut_k))
    
    DT::datatable(clustTab,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")),
                  rownames = F)
    
  })
  
  output$SIFPreview <- DT::renderDataTable({
    
    jacc_df <- jacc_react_Intra()
    dist_cutoff <- input$JaccDistCutoff
    jacc_df <- as.data.frame(jacc_df)
    jacc_df_melt <- melt(jacc_df)
    jacc_cols <- colnames(jacc_df)
    jacc_df_melt$Pathway_B <- rep_len(jacc_cols,length.out = nrow(jacc_df_melt))
    colnames(jacc_df_melt)[c(1,2)] <- c("Pathway_A","Jaccard_Distance")
    jacc_df_melt <- jacc_df_melt[which(round(jacc_df_melt$Jaccard_Distance,4) <= dist_cutoff),]
    jacc_df_melt <- jacc_df_melt[which(jacc_df_melt$Pathway_A != jacc_df_melt$Pathway_B),]
    jacc_df_melt <- jacc_df_melt[!duplicated(t(apply(jacc_df_melt,1,sort))),]
    DT::datatable(jacc_df_melt,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")),
                  rownames = F)
    
  })
  
  output$SIFPreviewINTER <- DT::renderDataTable({
    
    jacc_df <- jacc_react_Inter()
    dist_cutoff <- input$JaccDistCutoffINTER
    jacc_df <- as.data.frame(jacc_df)
    jacc_df_melt <- melt(jacc_df)
    jacc_cols <- colnames(jacc_df)
    jacc_df_melt$Pathway_B <- rep_len(jacc_cols,length.out = nrow(jacc_df_melt))
    colnames(jacc_df_melt)[c(1,2)] <- c("Pathway_A","Jaccard_Distance")
    jacc_df_melt <- jacc_df_melt[which(round(jacc_df_melt$Jaccard_Distance,4) <= dist_cutoff),]
    jacc_df_melt <- jacc_df_melt[which(jacc_df_melt$Pathway_A != jacc_df_melt$Pathway_B),]
    jacc_df_melt <- jacc_df_melt[!duplicated(t(apply(jacc_df_melt,1,sort))),]
    DT::datatable(jacc_df_melt,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")),
                  rownames = F)
    
  })
  
  
  ####----Plots----####
  
  ## Jaccard Heatmap - INTRA
  output$JaccHeatmapIntra <- renderPlot({
    
    ## User inputs
    row_font <- input$HeatRowFontIntra
    col_font <- input$HeatColFontIntra
    row_name <- input$HeatRowNamesIntra
    col_name <- input$HeatColNamesIntra
    col_pall <- input$ColorPaletteIntra
    clust_me <- input$ClustMethodIntra
    row_dend <- input$HeatRowDendHeight
    col_dend <- input$HeatColDendHeight
    
    ## Data matrix
    jacc_df <- jacc_react_Intra()
    jacc_mat <- as.matrix(jacc_df)
    
    ## Color choice
    minimum = min(jacc_mat)
    maximum = max(jacc_mat)
    #if (abs(min(jacc_mat)) > abs(max(jacc_mat))) {
    #  jacc_mat[jacc_mat < -abs(max(jacc_mat))] = -abs(max(jacc_mat))
    #} else {
    #  jacc_mat[jacc_mat > abs(min(jacc_mat))] = abs(min(jacc_mat))
    #}
    bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
    #Heatmap color
    col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
    if (col_pall == "original") {
      HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (col_pall %in% col_sets) {
      HeatMap_Colors <- brewer.pal(n = 5, col_pall)
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (col_pall == "Inferno") {
      hmcols <- inferno(500)
    }
    else if (col_pall == "Viridis") {
      hmcols <- viridis(500)
    }
    else if (col_pall == "Plasma") {
      hmcols <- plasma(500)
    }
    else if (col_pall == "OmniBlueRed") {
      hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
    }
    else if (col_pall == "LightBlueBlackRed") {
      hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
    }
    else if (col_pall == "GreenBlackRed") {
      hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
    }
    pheatmap::pheatmap(jacc_mat,
                             fontsize_row = row_font,
                             fontsize_col = col_font,
                             show_rownames = row_name,
                             show_colnames = col_name,
                             treeheight_col = col_dend,
                             treeheight_row = row_dend,
                             angle_col = 90,
                             border_color = NA,
                             clustering_method = clust_me,
                             color = hmcols)
    
  })
  
  output$JaccHeatmapInter <- renderPlot({
    
    ## User inputs
    row_font <- input$HeatRowFontInter
    col_font <- input$HeatColFontInter
    row_name <- input$HeatRowNamesInter
    col_name <- input$HeatColNamesInter
    col_pall <- input$ColorPaletteInter
    #clust_me <- input$ClustMethodInter
    row_dend <- input$HeatRowDendHeightInter
    col_dend <- input$HeatColDendHeightInter
    
    ## Data matrix
    jacc_df <- jacc_react_Inter()
    jacc_mat <- as.matrix(jacc_df)
    
    ## Color choice
    minimum = min(jacc_mat)
    maximum = max(jacc_mat)
    #if (abs(min(jacc_mat)) > abs(max(jacc_mat))) {
    #  jacc_mat[jacc_mat < -abs(max(jacc_mat))] = -abs(max(jacc_mat))
    #} else {
    #  jacc_mat[jacc_mat > abs(min(jacc_mat))] = abs(min(jacc_mat))
    #}
    bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
    #Heatmap color
    col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
    if (col_pall == "original") {
      HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (col_pall %in% col_sets) {
      HeatMap_Colors <- brewer.pal(n = 5, col_pall)
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (col_pall == "Inferno") {
      hmcols <- inferno(500)
    }
    else if (col_pall == "Viridis") {
      hmcols <- viridis(500)
    }
    else if (col_pall == "Plasma") {
      hmcols <- plasma(500)
    }
    else if (col_pall == "OmniBlueRed") {
      hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
    }
    else if (col_pall == "LightBlueBlackRed") {
      hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
    }
    else if (col_pall == "GreenBlackRed") {
      hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
    }
    pheatmap::pheatmap(jacc_mat,
                             fontsize_row = row_font,
                             fontsize_col = col_font,
                             show_rownames = row_name,
                             show_colnames = col_name,
                             treeheight_col = col_dend,
                             treeheight_row = row_dend,
                             angle_col = 90,
                             border_color = NA,
                             clustering_method = clust_me,
                             color = hmcols)
    
  })
  
  output$JaccDendoIntra <- renderPlotly({
    
    VisType <- input$ConnView
    
    if (VisType == "rectangle" | VisType == "phylogenic") {
    
      jacc_mat <- jacc_react_Intra()
      clust_me <- input$ClustMethodIntra
      nclust <- input$NumClusters
      LableChoice <- input$ShowConnLabels
      FontSize <- input$ConnFontSize
      
      
      hc <- hclust(dist(jacc_mat), method = clust_me)
      hc_dend <- as.dendrogram(hc)
      dend_data <- dendro_data(hc_dend)
      lab_order <- label(dend_data)[3]
      lab_vect <- lab_order$label
      clust_df <- data.frame(Pathways = rownames(jacc_mat),cluster = cutree(hc,k = nclust))
      clust_df2 <- clust_df[match(lab_vect,clust_df$Pathways),]
      colnames(clust_df2)[1] <- "label"
      clust_df2$cluster <- as.factor(clust_df2$cluster)
      
      dend_data[["labels"]] <- merge(dend_data[["labels"]],clust_df2, by="label")
      
      if (VisType == "rectangle") {
        
        FontSize <- FontSize * 6
        g <- ggdendrogram(dend_data,
                          rotate = TRUE,
                          theme_dendro = TRUE,
                          leaf_labels = F,
                          segments = T,
                          labels = F) + 
          labs(title = "") +
          scale_y_continuous(expand = c(0.4, 0))
        if (LableChoice == TRUE) {
          g <- g + geom_text(data=label(dend_data), aes(x, y, label=label, hjust=1, color=cluster),size=FontSize)
        }
        
        gp <- ggplotly(g)
        
      }
      else if (VisType == "phylogenic") {
        
        if (is.null(input$PhyloLayout) == TRUE) {
          PhyloLayout <- "layout.auto"
        }
        else if (is.null(input$PhyloLayout) == FALSE) {
          PhyloLayout <- input$PhyloLayout
        }
        
        
        if (LableChoice == TRUE) {
          k = nclust
          gp <- ggplotly(
            fviz_dend(
              hc_dend,
              cex = FontSize,
              k = nclust,
              show_labels = TRUE,
              type = "phylogenic",
              phylo_layout = PhyloLayout)
          )
          for(i in seq(k, 1)){
            gp$x$data[[i+1]]$text <- gp$x$data[[i+1+k]]$text
          }
          
        }
        else if (LableChoice == FALSE) {
          k = nclust
          gp <- ggplotly(
            fviz_dend(
              hc,
              cex = FontSize,
              k = k,
              show_labels = TRUE,
              type = "phylogenic",
              phylo_layout = PhyloLayout)
          )
          gp$x$data[[1]]$text <- NA 
          for(i in seq(k, 1)){
            gp$x$data[[i+1]]$text <- gp$x$data[[i+1+k]]$text
            gp$x$data[[i + 1 + k]] <- NULL
          }
          
        }
        
      }
      gp
    }
    
  })
  
  output$JaccDendoIntraCirc <- renderPlot({
    
    VisType <- input$ConnView
    
    if (VisType == "circular") {
      
      jacc_mat <- jacc_react_Intra()
      clust_me <- input$ClustMethodIntra
      nclust <- input$NumClusters
      LableChoice <- input$ShowConnLabels
      FontSize <- input$ConnFontSize
      
      hc <- hclust(dist(jacc_mat), method = clust_me)
      hc_dend <- as.dendrogram(hc)
      f <- fviz_dend(hc_dend, k = nclust, cex = FontSize, color_labels_by_k = T, type = "circular", show_labels = LableChoice) +
        theme(legend.position="none") +
        labs(title = "") +
        theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),
              axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank())
      
      
      f
      
    }
    
  })
  
  ####----Download Handlers----####
  
  ## Download Cluster Table - INTRA
  output$dnldClustTabIntra <- downloadHandler(
    filename = function() {
      
      gs.u <- input$UserPathwayFile
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Cluster_Result.txt", sep = '')
      
    },
    content = function(file) {
      
      ## Data matrix
      jacc_df <- jacc_react_Intra()
      jacc_mat <- as.matrix(jacc_df)
      clust_me <- input$ClustMethodIntra
      
      hm_res <- pheatmap::pheatmap(jacc_mat,
                               clustering_method = clust_me,
                               silent = T)
      
      jacc_df <- jacc_react_Intra()
      cut_k <- input$NumClusters
      clustTab <- data.frame(Pathways = rownames(jacc_df),
                             Cluster = cutree(hm_res$tree_row,k = cut_k))
      write_delim(clustTab,file,delim = '\t')
       
    }
  )
  
  output$dnldSIFTabIntra <- downloadHandler(
    filename = function() {
      
      gs.u <- input$UserPathwayFile
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_sif_dend.sif", sep = '')
      
    },
    content = function(file) {
      
      jacc_df <- jacc_react_Intra()
      dist_cutoff <- input$JaccDistCutoff
      jacc_df <- as.data.frame(jacc_df)
      jacc_df_melt <- melt(jacc_df)
      jacc_cols <- colnames(jacc_df)
      jacc_df_melt$Pathway_B <- rep_len(jacc_cols,length.out = nrow(jacc_df_melt))
      colnames(jacc_df_melt)[c(1,2)] <- c("Pathway_A","Jaccard_Distance")
      jacc_df_melt <- jacc_df_melt[which(jacc_df_melt$Jaccard_Distance <= dist_cutoff),]
      jacc_df_melt <- jacc_df_melt[which(jacc_df_melt$Pathway_A != jacc_df_melt$Pathway_B),]
      jacc_df_melt <- jacc_df_melt[!duplicated(t(apply(jacc_df_melt,1,sort))),]
      write_delim(jacc_df_melt,file,delim = '\t')
      
      
    }
  )
  
  output$dnldSIFTabInter <- downloadHandler(
    filename = function() {
      
      gs.u <- input$UserPathwayFile
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_InterPathway_sif_dend.sif", sep = '')
      
    },
    content = function(file) {
      
      jacc_df <- jacc_react_Inter()
      dist_cutoff <- input$JaccDistCutoffINTER
      jacc_df <- as.data.frame(jacc_df)
      jacc_df_melt <- melt(jacc_df)
      jacc_cols <- colnames(jacc_df)
      jacc_df_melt$Pathway_B <- rep_len(jacc_cols,length.out = nrow(jacc_df_melt))
      colnames(jacc_df_melt)[c(1,2)] <- c("Pathway_A","Jaccard_Distance")
      jacc_df_melt <- jacc_df_melt[which(jacc_df_melt$Jaccard_Distance <= dist_cutoff),]
      jacc_df_melt <- jacc_df_melt[which(jacc_df_melt$Pathway_A != jacc_df_melt$Pathway_B),]
      jacc_df_melt <- jacc_df_melt[!duplicated(t(apply(jacc_df_melt,1,sort))),]
      write_delim(jacc_df_melt,file,delim = '\t')
      
      
    }
  )
  
  ## Download Jaccard output - Intra
  output$dnldJaccTabIntra <- downloadHandler(
    filename = function() {
      
      gs.u <- input$UserPathwayFile
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Jaccard_Connectivity.txt", sep = '')
      
    },
    content = function(file) {
      
      jacc_df <- jacc_react_Intra()
      jacc_df <- as.data.frame(jacc_df)
      jacc_df$Pathways <- rownames(jacc_df)
      jacc_df <- jacc_df %>%
        relocate(Pathways)
      write_delim(jacc_df,file,delim = '\t')
      
    }
  )
  
  ## Download Jaccard output - Inter
  output$dnldJaccTabInter <- downloadHandler(
    filename = function() {
      
      gs.u <- input$UserPathwayFile
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Jaccard_Connectivity_Interpathway.txt", sep = '')
      
    },
    content = function(file) {
      
      jacc_df <- jacc_react_Inter()
      jacc_df <- as.data.frame(jacc_df)
      jacc_df$Pathways <- rownames(jacc_df)
      jacc_df <- jacc_df %>%
        relocate(Pathways)
      write_delim(jacc_df,file,delim = '\t')
      
    }
  )
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)



















