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
  feat <- sprintf(ifelse(grepl(" ", feat), "`%s`", "%s"), feat)
  tab <- coxph(as.formula(paste0("Surv(time,ID) ~ ",feat)),data = df)
  #tab <- coxph(Surv(time,ID) ~ Feature, data = df)
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

SurvPlot <- function(fit,df,title,ylab,pval,conf,legend,median,xlim,xScale = "Years",xBreaks = 365.25) {
  if (toupper(xScale) == "MONTHS") {
    xScale <- "d_m"
    xLabel <- "Months"
  } else {
    #breakTime <- ifelse(max(df[,"time"]) < 365.25,NULL,365.25)
    xScale <- "d_y"
    xLabel <- "Years"
  }
  ggsurv <- survminer::ggsurvplot(fit, data = df, risk.table = TRUE,
                                  title = title,
                                  xscale = c(xScale),
                                  break.time.by=xBreaks,
                                  xlab = xLabel,
                                  ylab = ylab,
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
  ggsurv
}

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














