# Last Update 07/01/2026

date_to_gene <- function(vec,mouse = FALSE) {
  if (is.null(vec)) exit("Please provide vector argument")
  if (!mouse) {
    vec <- gsub('(\\d+)-(Mar)','MARCH\\1',vec)
    vec <- gsub('(\\d+)-(Sep)','SEPT\\1',vec)
    vec <- gsub('(\\d+)-(Dec)','DEC\\1',vec)
  } else {
    vec <- gsub('(\\d+)-(Mar)','March\\1',vec)
    vec <- gsub('(\\d+)-(Sep)','Sept\\1',vec)
    vec <- gsub('(\\d+)-(Dec)','Dec\\1',vec)
  }
  return(vec)
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

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

highlow2 = function(mat) {
  new_mat = mat;
  new_mat[mat > quantile(as.numeric(mat), na.rm = T)[3]] = "Above_Median";
  new_mat[mat <= quantile(as.numeric(mat), na.rm = T)[3]] = "Below_Median";
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
  string3 <- gsub(" ","_",string3)
  return(string3)
}

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

get_tabData <- function(tab) {
  ## Developed on package "suvival v3.4-0 and survminder v0.4.9
  Variable <- as.character(tab[["terms"]][[3]])
  N <- tab[["n"]]
  out <- capture.output(summary(tab))
  coef_line1 <- strsplit(grep(paste0("^",Variable),out, value = T)[1], "\\s+")[[1]]
  coef_line2 <- strsplit(grep(paste0("^",Variable),out, value = T)[2], "\\s+")[[1]]
  lik_line <- grep("^Likelihood ratio test=",out,value = T)
  lik_line_P <- as.numeric(str_split(str_split(lik_line,", ")[[1]][2],"=")[[1]][2])
  tabData <- c(Variable = Variable, N = as.numeric(N), `Hazard Ratio` = as.numeric(coef_line1[3]), `Standard Error` = as.numeric(coef_line1[4]),
               Low = as.numeric(coef_line2[4]), High = as.numeric(coef_line2[5]), `P.Value` = as.numeric(lik_line_P))
  return(tabData)
}

GetColsOfType <- function(meta,type = c("discrete","continuous"), threshold = 0.75) {
  require(dplyr)
  MetaClass_num <- meta %>%
    type.convert(as.is = TRUE) %>%
    dplyr::select(where(is.numeric)) %>%
    names()
  IntMetaCols <- apply(meta[,MetaClass_num, drop = F],2,function(x) any(round(as.numeric(x)) != as.numeric(x)))
  IntMetaCols <- names(IntMetaCols)[which(IntMetaCols == T)]
  MetaClass_NonNum <- colnames(meta)[which(!colnames(meta) %in% MetaClass_num)]
  ShowOrNot <- apply(meta[,MetaClass_num,drop = F],2,function(x) any(length(levels(as.factor(x)))<(nrow(meta)*threshold)))
  discreteCols <- c(MetaClass_NonNum,names(ShowOrNot)[which(ShowOrNot == T)])
  discreteCols <- discreteCols[which(!discreteCols %in% IntMetaCols)]
  continuousCols <- colnames(meta)[which(!colnames(meta) %in% discreteCols)]
  if (toupper(type) == "DISCRETE") {
    return(discreteCols)
  } else if (toupper(type) == "CONTINUOUS") {
    return(continuousCols)
  } else {
    print("ERROR: Argument 'type' must be 'discrete' or 'continuous'")
  }
}

SubsetSurvData <- function(df,time,id,feat,feat2 = NULL) {
  SampleNameCol <- colnames(df)[1]
  if (is.null(feat2)) {
    df <- df[,c(SampleNameCol,time,id,feat)]
    df[,time] <- as.numeric(df[,time])
    df[,id] <- as.numeric(df[,id])
    if ("TopBottomCutP" %in% c(feat,feat2)) {
      df <- df[which(df$TopBottomCutP != "BetweenCutoff"),]
    }
    colnames(df)[which(colnames(df) == time)] <- "time"
    colnames(df)[which(colnames(df) == id)] <- "ID"
  } else {
    df <- df[,c(SampleNameCol,time,id,feat,feat2)]
    if ("TopBottomCutP" %in% c(feat,feat2)) {
      df <- df[which(df$TopBottomCutP != "BetweenCutoff"),]
    }
    colnames(df)[which(colnames(df) == time)] <- "time"
    colnames(df)[which(colnames(df) == id)] <- "ID"
  }
  return(df)
}

CoxPHobj <- function(df,feat,ref) {
  df[,feat] <- as.factor(df[,feat])
  df[,feat] <- relevel(df[,feat], ref = ref)
  #colnames(df)[which(colnames(df) == feat)] <- gsub("-","_",feat)
  #feat <- gsub("-","_",feat)
  colnames(df)[which(colnames(df) == feat)] <- gsub("[[:punct:]]","_",feat)
  feat <- gsub("[[:punct:]]","_",feat)
  # checks for spaces, first character numerics
  feat <- sprintf(ifelse((grepl(" ", feat) | !is.na(suppressWarnings(as.numeric(substring(feat, 1, 1))))), "`%s`", "%s"), feat)
  tab <- coxph(as.formula(paste0("Surv(time,ID) ~ ",feat)),data = df)
  return(tab)
}

CoxPHtabUni <- function(obj,Feature = NULL) {
  if (!is.null(Feature)) {
    obj <- obj %>%
      gtsummary::tbl_regression(exp = TRUE,
                                label = list(Feature ~ Feature)) %>%
      as_gt()
  } else {
    obj <- obj %>%
      gtsummary::tbl_regression(exp = TRUE) %>%
      as_gt()
  }
  tab_df <- as.data.frame(obj)
  tab_df <- tab_df %>%
    dplyr::select(label,estimate,ci,p.value)
  colnames(tab_df) <- c("Characteristic","Hazard Ratio","95% Confidence Interval","P.Value")
  tab_df <- sapply(tab_df,function(x) { gsub("<br />", "", x) })
  return(tab_df)
}

SurvPlotExpl <- function(CutPlabel,surv_time_col,geneset_name,scoreMethod,metacol_sampletype,SampleTypeSelected,Feature,subFeature,Pval_Tab,HR_Tab) {
  SurvDateType <- sub("\\..*","",surv_time_col)
  pval <- get_lik_pval(Pval_Tab)
  if (as.numeric(pval) < 0.05) {
    pval_char <- "strongly associated"
  } else if (as.numeric(pval) >= 0.05 & as.numeric(pval) < 0.1) {
    pval_char <- "moderately associated"
  } else if (as.numeric(pval) >= 0.1) {
    pval_char <- "not associated"
  }
  if (isTruthy(SampleTypeSelected)) {
    if (SampleTypeSelected != "All Sample Types") {
      if (Feature != "Show all Samples") {
        line1 <- paste0("<li><b>",SurvDateType,"</b> survival analysis ", metacol_sampletype," of <b>",SampleTypeSelected,"</b> Patients.</li>")
        line2 <- paste0("<li>The dataset is filtered by <b>",Feature,"</b> - <b>",subFeature,"</b>.</li>")
      } else {
        line1 <- paste0("<li><b>",SurvDateType,"</b> survival analysis ", metacol_sampletype," of <b>",SampleTypeSelected,"</b> Patients.</li>")
        line2 <- NULL
      }
    } else {
      if (Feature != "Show all Samples") {
        line1 <- paste0("<li><b>",SurvDateType,"</b> survival analysis of all sample types.</li>")
        line2 <- paste0("<li>The dataset is filtered by <b>",Feature,"</b> - <b>",subFeature,"</b>.</li>")
      } else {
        line1 <- paste0("<li><b>",SurvDateType,"</b> survival analysis of all sample types.</li>")
        line2 <- NULL
      }
    }
  } else {
    if (Feature != "Show all Samples") {
      line1 <- paste0("<li><b>",SurvDateType,"</b> survival analysis.</li>")
      line2 <- paste0("<li>The dataset is filtered by <b>",Feature,"</b> - <b>",subFeature,"</b>.</li>")
    } else {
      line1 <- paste0("<li><b>",SurvDateType,"</b> survival analysis.</li>")
      line2 <- NULL
    }
  }
  line3 <- paste0("<li>Kaplan-Meier survival curve categorized by <b>",geneset_name,"</b> <b>",scoreMethod,"</b> ",CutPlabel,".</li>")
  if (nrow(HR_Tab) == 3) {
    HR <- HR_Tab[3,2]
    HR_char <- ifelse(as.numeric(gsub(",","",(HR))) > 1,"high risk","low risk")
    chacteristic <- str_squish(HR_Tab[3,1])
    line4 <- paste("<li>Cox hazard regression analysis finds a Likelihood Ratio P.value of <b>",pval,"</b> and a Hazard Ratio of <b>",HR,"</b>, <b>",chacteristic,"</b> <b>",geneset_name,
                   "</b> is <b>",pval_char,"</b> with <b>",HR_char,"</b> for <b>",SurvDateType,"</b>.</li>",sep = "")
  } else {
    line4 <- paste("<li>Cox hazard regression analysis finds a Likelihood Ratio P.value of <b>",pval,"</b> shows that <b>",geneset_name,
                   "</b> is <b>",pval_char,"</b> with <b>",SurvDateType,"</b>.</li>",sep = "")
  }
  if (is.null(line2)) {
    return(HTML(paste0("<ul>",line1,line3,line4,"</ul>")))
  } else if (!is.null(line2)) {
    return(HTML(paste0("<ul>",line1,line2,line3,line4,"</ul>")))
  }
  
}
SurvPlot <- function(fit,df,title,ylab,pval,conf,legend,median,ylim = c(0,1),xlim,xScale = "Years",xBreaks = 365.25) {
  if (toupper(xScale) == "MONTHS") {
    xScale <- "d_m"
    xLabel <- "Months"
  } else {
    #breakTime <- ifelse(max(df[,"time"]) < 365.25,NULL,365.25)
    xScale <- "d_y"
    xLabel <- "Years"
  }
  surv_summary(fit)
  ggsurv <- survminer::ggsurvplot(fit, data = df, risk.table = TRUE,
                                  title = title,
                                  xscale = c(xScale),
                                  break.time.by=xBreaks,
                                  xlab = xLabel,
                                  ylab = ylab,
                                  ylim = ylim,
                                  submain = "Based on Kaplan-Meier estimates",
                                  caption = "created with survminer",
                                  pval = pval,
                                  conf.int = conf,
                                  ggtheme = theme_bw(),
                                  font.title = c(16, "bold"),
                                  font.submain = c(12, "italic"),
                                  font.caption = c(12, "plain"),
                                  font.x = c(14, "plain"),
                                  font.y = c(14, "plain"),
                                  font.tickslab = c(12, "plain"),
                                  legend = legend,
                                  risk.table.height = 0.20,
                                  surv.median.line = median
  )
  if (median != "none") {
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
  if (!is.null(xlim)) {
    ggsurv$plot$coordinates$limits$x <- c(0,xlim)
    ggsurv$table$coordinates$limits$x <- c(0,xlim)
  }
  ggsurv$table <- ggsurv$table + theme_cleantable()
  ggsurv$plot <- ggsurv$plot + theme(
    legend.key.height = unit(1.5, "cm"), # Adjust height (e.g., to 1.5 cm)
    legend.spacing.y = unit(0.5, "cm")   # Adjust vertical spacing between legend items
  )

  
  #ggsurv$table <- ggsurv$table +
  #  scale_y_discrete(labels = names(fit[["strata"]]))
  
  
  ggsurv
}

#tab <- ggsurvtable(fit, data = df)

SurvPlotTitle <- function(SampleTypeSelected,geneset_name = NULL,scoreMethodLab,Feature,subFeature,CutPLabel,univar = NULL,multivar = NULL) {
  if (isTruthy(geneset_name)) {
    primeFeatLab <- paste0(geneset_name," (",scoreMethodLab,") at ",CutPLabel)
  }
  if (isTruthy(univar)) {
    primeFeatLab <- paste0(univar)
  }
  if (isTruthy(multivar)) {
    primeFeatLab <- paste0(multivar)
  }
  if (isTruthy(SampleTypeSelected)) {
    if (SampleTypeSelected != "All Sample Types") {
      if (Feature != "Show all Samples") {
        PlotTitle <- paste0("Survival Curve of ",primeFeatLab,"\nAcross Patients Classified as ",SampleTypeSelected," - ",Feature," (",subFeature,")")
      } else {
        PlotTitle <- paste0("Survival Curve of ",primeFeatLab,"\nAcross Patients Classified as ",SampleTypeSelected)
      }
    } else {
      if (Feature != "Show all Samples") {
        PlotTitle <- paste0("Survival Curve of ",primeFeatLab,"\nAcross Patients Classified as ",Feature," (",subFeature,")")
      } else {
        PlotTitle <- paste0("Survival Curve of ",primeFeatLab,"\nAcross All Patients")
      }
    }
  } else {
    if (Feature != "Show all Samples") {
      PlotTitle <- paste0("Survival Curve of ",primeFeatLab,"\nAcross Patients Classified as ",Feature," (",subFeature,")")
    } else {
      PlotTitle <- paste0("Survival Curve of ",primeFeatLab,"\nAcross All Patients")
    }
  }
  return(PlotTitle)
}

densPlot <- function(score,quant,xlab,ylab,title,CutPlabel,ShowQuartile = TRUE,user_vline = 0){
  dens_data <- density(score[,1],na.rm = T)
  y_max <- max(dens_data$y)
  y_max_int <- y_max/6
  p <- ggplot(score, aes(x=score[,1])) +
    geom_density(color="darkblue", fill="lightblue", alpha = 0.4) +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(title) +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          plot.title = element_text(size = 20))
  if (!is.null(CutPlabel)) {
    p <- p + geom_vline(data = quant, aes(xintercept = Quantile), linetype = "dashed", color = "darkblue", linewidth = 1)
    if (nrow(quant) == 1) {
      p <- p + geom_text(aes(quant[1,1],y_max-(y_max_int/6),label = paste(as.character(quant[1,1]),CutPlabel),hjust = -0.1,vjust = -0.1),size = 6, check_overlap = T)
    } else if (nrow(quant) == 2) {
      p <- p + geom_text(aes(quant[1,1],y_max-(y_max_int/6),label = paste(as.character(quant[1,1]),CutPlabel[1]),hjust = -0.1,vjust = -0.1),size = 6, check_overlap = T)
      p <- p + geom_text(aes(quant[2,1],y_max-y_max_int,label = paste(as.character(quant[2,1]),CutPlabel[2]),hjust = -0.1,vjust = -0.1),size = 6, check_overlap = T)
    } else if (nrow(quant) == 3) {
      p <- p + geom_text(aes(quant[1,1],y_max-(y_max_int/6),label = paste(as.character(quant[1,1]),CutPlabel[1],sep="\n"),hjust = -0.1,vjust = 0.5),size = 6, check_overlap = T)
      p <- p + geom_text(aes(quant[2,1],y_max-y_max_int,label = paste(as.character(quant[2,1]),CutPlabel[2],sep="\n"),hjust = -0.1,vjust = 0.5),size = 6, check_overlap = T)
      p <- p + geom_text(aes(quant[3,1],y_max-(y_max_int*2),label = paste(as.character(quant[3,1]),CutPlabel[3],sep="\n"),hjust = -0.1,vjust = 0.5),size = 6, check_overlap = T)
    }
  } else {
    if (ShowQuartile == TRUE) {
      p <- p + geom_vline(data = quant, aes(xintercept = Quantile), linetype = "dashed", color = "darkblue", linewidth = 1)
    }
    if (user_vline != 0) {
      p <- p + geom_vline(xintercept = user_vline, linetype = "dashed", color = "darkred", linewidth = 1)
    }
  }
  
  return(p)
}

CoxPHsumm <- function(CoxPHobj,bivarAdd = FALSE,bivarInt = FALSE) {
  out <- capture.output(summary(CoxPHobj))
  xph <- capture.output(cox.zph(CoxPHobj))
  con_line <- grep("^Concordance=",out,value = T)
  lik_line <- grep("^Likelihood ratio test=",out,value = T)
  wal_line <- grep("^Wald test",out,value = T)
  sco_line <- grep("^Score ",out,value = T)
  text <- paste("CoxH Summary:",con_line,lik_line,wal_line,sco_line,"","Proportional Hazards assumption:",xph[1],xph[2],xph[3],sep = "\n")
  if (bivarAdd) {
    text <- paste("CoxH Summary:",con_line,lik_line,wal_line,sco_line,"","Proportional Hazards assumption:",xph[1],xph[2],xph[3],xph[4],sep = "\n")
  }
  if (bivarInt) {
    text <- paste("CoxH Summary:",con_line,lik_line,wal_line,sco_line,"","Proportional Hazards assumption:",xph[1],xph[2],xph[3],xph[4],xph[5], sep = "\n")
  }
  return(cat(text))
}

survFeatRefSelect <- function(meta,Feature,na.rm = TRUE,cont = FALSE,hilo = TRUE) {
  Var_choices <- meta[,Feature]
  if (na.rm == TRUE) {
    Var_choices <- Var_choices[which(is.na(Var_choices) == FALSE)]
    Var_choices <- Var_choices[which(Var_choices != "Inf" & Var_choices != "N/A" & Var_choices != "n/a")]
    Var_choices <- Var_choices[grep("unknown",Var_choices,ignore.case = T, invert = T)]
  }
  if (cont == FALSE) {
    Var_choices <- unique(meta[,Feature])
    Var_choices <- sort(Var_choices, decreasing = T, na.last = T)
  } else if (cont == TRUE) {
    if (hilo == TRUE) {
      Var_choices <- c("Low","High")
    }
  }
  return(Var_choices)
}

#df <- meta_ssgsea_sdf
#colnames(df)[which(colnames(df) == Feature)] <- "Feature"
#labelled::var_label(df) <- list(
#  Feature = Feature #this variable is a numeric -> label works
#)
#forest_model(obj)

forestPlot_Simple <- function(obj,df,Feature,Font) {
  forest <- survminer::ggforest(obj,
                                data = df,
                                main = paste("Hazard Ratio Modeling: ",Feature,sep = ""),
                                fontsize = Font)
  return(forest)
}

linearityPlot <- function(obj,Feature,resid,pred,axisFont,mainFont,tickFont) {
  p <- survminer::ggcoxdiagnostics(obj,
                                   type = resid,
                                   sline = T,
                                   sline.se = T,
                                   ggtheme = theme_minimal(),
                                   ox.scale = pred)
  p <- ggpar(p,
             font.x = axisFont,
             font.y = axisFont,
             font.main = mainFont,
             font.tickslab = tickFont,
             main = paste("Linearity Plot Featuring: ",Feature, sep = ""),
             ylab = paste(str_to_title(resid)," Residuals", sep = "")
  )
  p
}

biVarAnova <- function(obj,obj2) {
  annova_res <- anova(obj,obj2)
  out <- capture.output(annova_res)
  line1 <- out[3]
  line2 <- out[4]
  line3 <- out[5]
  line4 <- out[6]
  line5 <- out[7]
  text <- paste("Model Comparison:",line1,line2,line3,line4,line5,sep = "\n")
  return(cat(text))
}

dnld_ui <- function(id,label) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::downloadButton(ns("dnld"), label = label)
  )
}

dnldPlot_server <- function(id, plot, file, height = 8, width = 8, units = "in", type = "gg") {
  shiny::moduleServer(id, function(input, output, session) {
    output$dnld <- shiny::downloadHandler(
      filename = function() {
        paste(file)
      },
      content = function(file) {
        if (type == "gg") {
          ggplot2::ggsave(filename = file, plot = plot,
                          width = width, height = height,
                          units = units)
        } else if (type == "forest") {
          dims <- get_wh(plot)
          svg(file, width = unname(dims[1]+1), height = unname(dims[2]+1))
          plot(plot)
          dev.off()
        } else if (type == "complex") {
          svg(filename = file, height = height, width = width)
          ComplexHeatmap::draw(plot)
          dev.off()
        }
      }
    )
  })
}

dnldDF_server <- function(id, df, file) {
  shiny::moduleServer(id, function(input, output, session) {
    output$dnld <- shiny::downloadHandler(
      filename = function() {
        paste(file)
      },
      content = function(file) {
        write.table(df,file,sep = '\t', row.names = F)
      }
    )
  })
}

## ============================================================================
## ESTIMATE-based tumor purity scoring + purity correction, no GCT round-trip
## ============================================================================
## Requires only the 'estimate' package (for SI_geneset.gmt and
## common_genes.txt shipped in its extdata/).
##
## install.packages("estimate", repos = "http://r-forge.r-project.org",
##                   dependencies = TRUE)


# =============================================================================
# Internal helpers: load gene set and common genes from estimate extdata
# =============================================================================

#' Locate and validate a file in the estimate package extdata directory
#' @param filename name of the file to locate
#' @return full path to the file
.estimate_extdata_path <- function(filename) {
  path <- system.file("extdata", filename, package = "estimate")
  if (!nzchar(path) || !file.exists(path)) {
    stop(
      "Could not find '", filename, "' in the installed 'estimate' package extdata/. ",
      "Confirm the package is installed:\n",
      "  install.packages('estimate', repos = 'http://r-forge.r-project.org', ",
      "dependencies = TRUE)"
    )
  }
  path
}


#' Read SI_geneset.gmt into a named list of gene symbol vectors
#'
#' @return named list with elements "StromalSignature" and "ImmuneSignature"
.load_gene_sets <- function() {
  gmt_path <- .estimate_extdata_path("SI_geneset.gmt")
  lines <- readLines(gmt_path)
  sets <- lapply(lines, function(l) {
    parts <- strsplit(l, "\t")[[1]]
    list(name = parts[1], genes = parts[-c(1, 2)])
  })
  out        <- lapply(sets, `[[`, "genes")
  names(out) <- vapply(sets, `[[`, character(1), "name")
  out
}


#' Read common_genes.txt and build a primary-symbol + synonym lookup table
#'
#' The returned data.frame has two columns:
#'   GeneSymbol_syn  -- synonym (or the primary symbol itself)
#'   GeneSymbol      -- canonical ESTIMATE common gene symbol
#'
#' Rows where no synonym existed retain the primary symbol in both columns,
#' so a simple left-join on GeneSymbol_syn covers both direct matches and
#' alias matches.
#'
#' @return data.frame with columns GeneSymbol_syn, GeneSymbol
.load_common_genes <- function() {
  path <- .estimate_extdata_path("common_genes.txt")
  cg   <- read.delim(path, stringsAsFactors = FALSE)
  
  # Replace bare "-" with NA so strsplit doesn't produce a spurious synonym
  cg$Synonyms <- gsub("^-$", NA, cg$Synonyms)
  
  syn_list        <- strsplit(cg$Synonyms, "|", fixed = TRUE)
  names(syn_list) <- cg$GeneSymbol
  
  syn_df          <- unique(stack(syn_list))
  colnames(syn_df) <- c("GeneSymbol_syn", "GeneSymbol")
  
  syn_df
}


# =============================================================================
# filter_common_genes: data.frame in, data.frame out
# =============================================================================

#' Filter an expression data.frame to the ESTIMATE common gene panel
#'
#' Replicates the logic of estimate::filterCommonGenes() without any file I/O.
#' Matching is attempted first on primary HGNC symbols, then on synonyms for
#' any genes that were not matched directly -- mirroring the alias-resolution
#' behaviour in tidyestimate::filter_common_genes().
#'
#' @param expr_df A data.frame or matrix, genes x samples. Gene symbols must
#'   be in rownames (default) or a named column specified by \code{gene_col}.
#' @param gene_col Optional. Column name containing gene symbols if they are
#'   not in rownames(expr_df). That column is dropped from the output.
#' @param id Character, either "GeneSymbol" (default) or "EntrezID". Must
#'   match whatever identifier type is used in the rows of expr_df.
#' @param tell_missing Logical. If TRUE (default), prints the number of
#'   common genes not found in expr_df.
#' @param find_alias Logical. If TRUE (default), attempts to resolve missing
#'   genes via the Synonyms column of common_genes.txt. Only unambiguous
#'   1-to-1 alias matches are accepted.
#'
#' @return A data.frame, genes x samples, with rownames set to canonical
#'   ESTIMATE GeneSymbol identifiers and columns matching the input samples.
filter_common_genes <- function(expr_df,
                                gene_col     = NULL,
                                id           = c("GeneSymbol", "EntrezID"),
                                tell_missing = TRUE,
                                find_alias   = TRUE) {
  
  id <- match.arg(id)
  
  # ---- Normalise to matrix with gene-symbol rownames ----------------------
  if (!is.null(gene_col)) {
    if (!gene_col %in% colnames(expr_df))
      stop("gene_col '", gene_col, "' not found in expr_df.")
    gene_symbols <- as.character(expr_df[[gene_col]])
    expr_df      <- expr_df[, setdiff(colnames(expr_df), gene_col), drop = FALSE]
    rownames(expr_df) <- gene_symbols
  }
  
  input_genes <- rownames(expr_df)
  if (is.null(input_genes) || any(is.na(input_genes)) || any(input_genes == ""))
    stop("Gene identifiers not found. Provide them as rownames(expr_df) or via gene_col.")
  
  if (any(duplicated(input_genes))) {
    warning(sum(duplicated(input_genes)), " duplicated gene identifiers; ",
            "keeping the first occurrence of each.")
    expr_df     <- expr_df[!duplicated(input_genes), , drop = FALSE]
    input_genes <- rownames(expr_df)
  }
  
  # ---- Load common genes --------------------------------------------------
  # common_genes.txt column used for matching depends on `id`
  cg_path   <- .estimate_extdata_path("common_genes.txt")
  cg        <- read.delim(cg_path, stringsAsFactors = FALSE)
  cg$Synonyms <- gsub("^-$", NA, cg$Synonyms)
  
  match_col <- id   # "GeneSymbol" or "EntrezID"
  cg_ids    <- as.character(cg[[match_col]])
  
  # ---- Step 1: direct primary-symbol match --------------------------------
  direct_hits <- input_genes %in% cg_ids
  matched_df  <- expr_df[direct_hits, , drop = FALSE]
  
  # Canonical rowname is the GeneSymbol for direct hits
  # (for EntrezID mode, remap rownames to GeneSymbol for downstream use)
  if (id == "EntrezID") {
    id_to_sym           <- setNames(cg$GeneSymbol, cg$EntrezID)
    rownames(matched_df) <- id_to_sym[rownames(matched_df)]
  }
  
  # ---- Step 2: synonym-based resolution for unmatched common genes --------
  if (find_alias && id == "GeneSymbol") {
    missing_common <- setdiff(cg$GeneSymbol, rownames(matched_df))
    
    if (length(missing_common) > 0) {
      # Build synonym -> primary symbol map (only rows for still-missing genes)
      syn_list <- strsplit(cg$Synonyms[cg$GeneSymbol %in% missing_common],
                           "|", fixed = TRUE)
      names(syn_list) <- cg$GeneSymbol[cg$GeneSymbol %in% missing_common]
      syn_df <- unique(stack(syn_list))
      colnames(syn_df) <- c("synonym", "canonical")
      
      # Unambiguous: synonym maps to exactly one canonical AND
      #              that synonym matches exactly one input gene
      syn_in_input   <- syn_df[syn_df$synonym %in% input_genes, ]
      n_canon_per_syn <- table(syn_in_input$synonym)
      n_syn_per_input <- table(syn_in_input$synonym)   # by construction <= 1 per loop
      
      # canonical must appear only once in the filtered set
      n_input_per_canon <- table(syn_in_input$canonical)
      
      unambiguous <- syn_in_input[
        n_canon_per_syn[syn_in_input$synonym] == 1 &
          n_input_per_canon[syn_in_input$canonical] == 1, ]
      
      if (nrow(unambiguous) > 0) {
        alias_df            <- expr_df[unambiguous$synonym, , drop = FALSE]
        rownames(alias_df)  <- unambiguous$canonical
        matched_df          <- rbind(matched_df, alias_df)
        message(length(unique(unambiguous$canonical)),
                " missing gene(s) resolved via synonym matching.")
      }
    }
  }
  
  # ---- Report coverage ----------------------------------------------------
  n_found    <- nrow(matched_df)
  n_common   <- nrow(cg)
  n_missing  <- n_common - n_found
  
  if (tell_missing) {
    message(sprintf(
      "Common genes found: %d / %d  |  mismatched / missing: %d",
      n_found, n_common, n_missing
    ))
  }
  
  matched_df
}


# =============================================================================
# Internal: enrichment scoring engine (ported from estimateScore source)
# =============================================================================

#' Per-sample rank normalization matching estimateScore() exactly
.rank_scale <- function(m) {
  Ng <- nrow(m)
  for (j in seq_len(ncol(m))) {
    m[, j] <- rank(m[, j], ties.method = "average")
  }
  10000 * m / Ng
}


#' Compute weighted KS-style enrichment score for one gene set, all samples.
#' Direct port of the inner loop in estimateScore().
.es_one_geneset <- function(m_ranked, gene.names, gene.set) {
  gene.overlap <- intersect(gene.set, gene.names)
  Ns <- ncol(m_ranked)
  N  <- nrow(m_ranked)
  
  if (length(gene.overlap) == 0) return(rep(NA_real_, Ns))
  
  gene.set2 <- match(gene.overlap, gene.names)
  Nh <- length(gene.set2)
  Nm <- N - Nh
  
  ES.vector <- numeric(Ns)
  
  for (S.index in seq_len(Ns)) {
    gene.list     <- order(m_ranked[, S.index], decreasing = TRUE)
    correl.vector <- m_ranked[gene.list, S.index]
    
    TAG    <- sign(match(gene.list, gene.set2, nomatch = 0))
    no.TAG <- 1 - TAG
    
    correl.vector <- abs(correl.vector)^0.25
    sum.correl    <- sum(correl.vector[TAG == 1])
    
    P0 <- no.TAG / Nm
    F0 <- cumsum(P0)
    Pn <- TAG * correl.vector / sum.correl
    Fn <- cumsum(Pn)
    
    ES.vector[S.index] <- sum(Fn - F0)
  }
  
  ES.vector
}


# =============================================================================
# estimate_score_df: main scoring function with integrated filtering
# =============================================================================

#' Compute ESTIMATE purity scores from a gene expression data.frame,
#' without any GCT file I/O
#'
#' Runs filter_common_genes() as the first step (optional but recommended),
#' then computes stromal, immune, ESTIMATE, and purity scores using the
#' exact enrichment scoring logic from estimate::estimateScore().
#'
#' @param expr_df A numeric data.frame or matrix, genes x samples.
#'   HGNC gene symbols in rownames (default) or a column given by
#'   \code{gene_col}.
#' @param gene_col Optional column name for gene symbols (see
#'   filter_common_genes()).
#' @param platform One of "affymetrix", "agilent", "illumina". The cosine
#'   purity formula was validated on Affymetrix arrays; a warning is issued
#'   for other platforms. Purity is computed regardless and set to NA where
#'   the cosine formula returns < 0.
#' @param filter_genes Logical. If TRUE (default), filter_common_genes() is
#'   run before scoring. Set FALSE only if you have already filtered.
#' @param find_alias Logical, passed to filter_common_genes(). Default TRUE.
#' @param tell_missing Logical, passed to filter_common_genes(). Default TRUE.
#' @param log_input Logical. If TRUE, expr_df is treated as log-transformed
#'   and exponentiated before processing. Rank order is preserved under
#'   monotone transforms so scores are unaffected, but this ensures the
#'   matrix passed to any downstream linear step is on the correct scale.
#' @param log_base,log_pseudocount Used when log_input = TRUE.
#' @param min_overlap Minimum signature gene overlap required per gene set.
#'
#' @return A data.frame: sample, stromal_score, immune_score,
#'   estimate_score, purity. Samples with out-of-bounds ESTIMATE scores
#'   have purity = NA.
estimate_score_df <- function(expr_df,
                              gene_col        = NULL,
                              platform        = c("affymetrix", "agilent", "illumina"),
                              filter_genes    = TRUE,
                              find_alias      = TRUE,
                              tell_missing    = TRUE,
                              log_input       = FALSE,
                              log_base        = 2,
                              log_pseudocount = 1,
                              min_overlap     = 15) {
  
  platform <- match.arg(platform)
  
  # ---- Step 1: filter to common genes -------------------------------------
  if (filter_genes) {
    expr_df <- filter_common_genes(
      expr_df,
      gene_col     = gene_col,
      id           = "GeneSymbol",
      tell_missing = tell_missing,
      find_alias   = find_alias
    )
    gene_col <- NULL   # rownames are now canonical symbols after filtering
  }
  
  # ---- Coerce to numeric matrix with gene rownames -------------------------
  if (!is.null(gene_col)) {
    if (!gene_col %in% colnames(expr_df))
      stop("gene_col '", gene_col, "' not found in expr_df.")
    gene_symbols     <- as.character(expr_df[[gene_col]])
    expr_df          <- expr_df[, setdiff(colnames(expr_df), gene_col), drop = FALSE]
    rownames(expr_df) <- gene_symbols
  }
  
  gene_symbols <- rownames(expr_df)
  if (is.null(gene_symbols) || any(is.na(gene_symbols)) || any(gene_symbols == ""))
    stop("Gene symbols not found. Provide them as rownames(expr_df) or via gene_col.")
  
  m <- as.matrix(expr_df)
  storage.mode(m) <- "numeric"
  rownames(m)      <- gene_symbols
  
  if (any(duplicated(rownames(m)))) {
    warning(sum(duplicated(rownames(m))),
            " duplicated gene symbols; keeping first occurrence.")
    m <- m[!duplicated(rownames(m)), , drop = FALSE]
  }
  
  if (log_input) {
    m <- log_base^m - log_pseudocount
    m[m < 0] <- 0
  }
  
  gene.names   <- rownames(m)
  sample.names <- colnames(m)
  
  # ---- Step 2: load gene sets and rank-scale ------------------------------
  gene_sets <- .load_gene_sets()
  m_ranked  <- .rank_scale(m)
  
  # ---- Step 3: enrichment scores ------------------------------------------
  N.gs     <- length(gene_sets)
  gs.names <- names(gene_sets)
  score.matrix <- matrix(0, nrow = N.gs, ncol = ncol(m))
  
  for (gs.i in seq_len(N.gs)) {
    gene.set     <- gene_sets[[gs.i]]
    gene.overlap <- intersect(gene.set, gene.names)
    message(sprintf("Gene set %d: %-20s  overlap = %d / %d",
                    gs.i, gs.names[gs.i],
                    length(gene.overlap), length(gene.set)))
    
    if (length(gene.overlap) < min_overlap) {
      if (length(gene.overlap) == 0) {
        score.matrix[gs.i, ] <- NA_real_
        next
      }
      warning(sprintf(
        "Gene set '%s': only %d / %d signature genes found (< min_overlap = %d). ",
        gs.names[gs.i], length(gene.overlap), length(gene.set), min_overlap),
        "Scores may be unreliable.")
    }
    
    score.matrix[gs.i, ] <- .es_one_geneset(m_ranked, gene.names, gene.set)
  }
  
  stromal_score  <- score.matrix[1, ]
  immune_score   <- score.matrix[2, ]
  estimate_score <- stromal_score + immune_score
  
  # ---- Step 4: purity -----------------------------------------------------
  if (platform != "affymetrix") {
    message(
      "Note: the cosine purity formula was validated on Affymetrix arrays. ",
      "For RNA-seq (TPM/FPKM) purity values are approximate; relative ",
      "ranking across samples is more reliable than absolute values."
    )
  }
  
  purity            <- cos(0.6049872018 + 0.0001467884 * estimate_score)
  purity[purity < 0] <- NA_real_
  
  data.frame(
    sample         = sample.names,
    stromal_score  = as.numeric(stromal_score),
    immune_score   = as.numeric(immune_score),
    estimate_score = as.numeric(estimate_score),
    purity         = as.numeric(purity),
    row.names      = NULL,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# correct_for_purity
# =============================================================================

#' Correct gene expression matrix for tumor purity
#'
#' Adjusts a genes x samples expression matrix to account for tumor purity.
#'   - mode = "immune" / "stroma": divide by (1 - purity)
#'   - mode = "tumor":             divide by purity
#'
#' Correction is always applied on the linear scale. If log_input = TRUE the
#' function exponentiates, corrects, then re-logs, returning a matrix on the
#' same scale as the input. Samples with NA purity receive NA columns in the
#' output rather than causing an error.
#'
#' @param expr_matrix Numeric matrix or data.frame, genes x samples.
#' @param purity Named numeric vector of purity values (0-1), names matching
#'   colnames(expr_matrix). Unnamed vectors accepted if order matches columns
#'   exactly (a warning is issued). NA values are allowed and propagated.
#' @param mode "immune", "stroma", or "tumor".
#' @param log_input Logical. TRUE if expr_matrix is log-transformed.
#' @param log_base Base used for log transform. Default 2.
#' @param log_pseudocount Pseudocount used in log transform. Default 1.
#' @param clip_purity Logical. Clip purity to purity_floor/purity_ceiling to
#'   avoid division by zero or extreme blow-up. Default TRUE.
#' @param purity_floor,purity_ceiling Bounds applied when clip_purity = TRUE.
#'
#' @return A matrix, same dimensions as expr_matrix, on the same scale as
#'   the input. Columns with NA purity are returned as NA.
correct_for_purity <- function(expr_matrix,
                               purity,
                               mode            = c("immune", "stroma", "tumor"),
                               log_input       = FALSE,
                               log_base        = 2,
                               log_pseudocount = 1,
                               clip_purity     = TRUE,
                               purity_floor    = 0.01,
                               purity_ceiling  = 0.99) {
  
  mode <- match.arg(mode)
  
  if (!is.matrix(expr_matrix) && !is.data.frame(expr_matrix))
    stop("expr_matrix must be a matrix or data.frame.")
  expr_matrix      <- as.matrix(expr_matrix)
  storage.mode(expr_matrix) <- "numeric"
  
  if (!is.numeric(purity)) stop("purity must be a numeric vector.")
  
  # ---- Align purity to sample columns -------------------------------------
  if (!is.null(names(purity))) {
    missing_samples <- setdiff(colnames(expr_matrix), names(purity))
    if (length(missing_samples) > 0)
      stop("Purity values missing for samples: ",
           paste(missing_samples, collapse = ", "))
    purity <- purity[colnames(expr_matrix)]
  } else {
    if (length(purity) != ncol(expr_matrix))
      stop("Unnamed purity vector length must match ncol(expr_matrix).")
    warning("purity vector is unnamed; assuming order matches columns of expr_matrix.")
  }
  
  # ---- Handle NA purity samples gracefully --------------------------------
  na_samples <- is.na(purity)
  if (any(na_samples)) {
    warning(sum(na_samples), " sample(s) have NA purity; those columns will ",
            "be NA in the output: ",
            paste(names(purity)[na_samples], collapse = ", "))
  }
  
  valid <- !na_samples
  if (!any(valid)) stop("All samples have NA purity; cannot correct any samples.")
  
  if (any(purity[valid] < 0 | purity[valid] > 1))
    stop("Non-NA purity values must be between 0 and 1.")
  
  if (clip_purity)
    purity[valid] <- pmin(pmax(purity[valid], purity_floor), purity_ceiling)
  
  denom <- if (mode %in% c("immune", "stroma")) (1 - purity) else purity
  
  if (any(denom[valid] <= 0))
    stop("Denominator is zero for one or more samples after clipping; ",
         "adjust purity_floor / purity_ceiling.")
  
  # ---- Apply correction on linear scale -----------------------------------
  working <- expr_matrix
  if (log_input) {
    working <- log_base^working - log_pseudocount
    working[working < 0] <- 0
  }
  
  corrected              <- working
  corrected[, valid]     <- sweep(working[, valid, drop = FALSE],
                                  MARGIN = 2, STATS = denom[valid], FUN = "/")
  corrected[, na_samples] <- NA_real_
  
  if (log_input)
    corrected[, valid] <- log(corrected[, valid, drop = FALSE] + log_pseudocount,
                              base = log_base)
  
  corrected
}


## ============================================================================
## End-to-end usage
## ============================================================================
## tpm_df: genes (HGNC rownames) x samples, linear TPM/FPKM
##
## # Scoring (filter_common_genes runs automatically as step 1):
## scores <- estimate_score_df(tpm_df, platform = "illumina")
## purity <- setNames(scores$purity, scores$sample)
##
## # Or run filtering separately and pass pre-filtered matrix:
## filtered_df <- filter_common_genes(tpm_df)
## scores      <- estimate_score_df(filtered_df, filter_genes = FALSE)
##
## # Purity correction:
## immune_corrected <- correct_for_purity(tpm_df, purity, mode = "immune")
## tumor_corrected  <- correct_for_purity(tpm_df, purity, mode = "tumor")
##
## # If your matrix is log2(TPM + 1):
## tumor_corrected  <- correct_for_purity(log2_tpm_df, purity, mode = "tumor",
##                                         log_input = TRUE)

