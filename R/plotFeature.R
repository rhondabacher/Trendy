#' @title Plot features of interest

#' @description plot each feature with (or without) the fitted trend.

#' @inheritParams trendy
#' @param featureNames a list of genes or features to plot
#' @param showFit whether to plot the segmented regression 
#'  fitting (default is TRUE)
#' @param simple if TRUE the plot will not highlight the breakpoints and segments 
#' and will only display a black fitted line. (default is FALSE)
#' @param showLegend if TRUE and simple=FALSE then a legend will be output (default = TRUE)
#' @param trendyOutData segmented regression fitting result from 
#'  running trendy(); if showFit is TRUE and trendyOutData is NULL, then the
#'  segmented regression will be fit for each of the genes and it may 
#'  take longer to run
#' @param cexLegend cex option for sizing of legend text, default is 1.
#' @param legendLocation whether to place the legend to the right 'side' of each plot or at the 'bottom' of a multo-panelled plot (default is 'side').
#' @param xlab x-axis name
#' @param ylab y-axis name
#' @param segColors define colors for the 'breakpoint', and 'up', 'same', 'down' segments (default: segColors = 	c("chartreuse3", "coral1", "black", "cornflowerblue"))
#' @param customTitle default is set the plot title as the name of the feature. Other titles can be specified here. 
#' @param customLabels.x specify x-axis tick labels instead of using the default values from tVectIn.
#' @param spacing.x specify x-axis tick spacing, smaller values give more tick marks.

#' @return plot of gene expression and fitted line
#' @examples d1 <- rbind(c(rep(1,50),seq_len(50)), rev(seq_len(100)))
#'  rownames(d1) <- c("g1","g2")
#'  colnames(d1) <- paste0("time", seq_len(100))
#'  plotFeature(d1, featureNames=c("g1","g2"))
#' @author Ning Leng and Rhonda Bacher
#' @import graphics
#' @import grDevices
#' @inheritParams trendy
#' @export
plotFeature <- 
  function(Data, tVectIn = NULL, featureNames, showFit = TRUE, simple=FALSE,
           showLegend = TRUE, trendyOutData = NULL, cexLegend=1, 
           legendLocation = "side",
           xlab = "Time", 
           ylab = "Normalized Expression",
           segColors = c("chartreuse3", "coral1", "black", "cornflowerblue"),
           customTitle= NULL, customLabels.x=NULL, spacing.x=NULL) 
    
  {
    if (methods::is(Data, "SummarizedExperiment")) {
      if (is.null( SummarizedExperiment::assayNames(Data)) || 
          SummarizedExperiment::assayNames(Data)[1] != "Counts") {
        message("renaming the first element in assays(Data) to 'Counts'")
        SummarizedExperiment::assayNames(Data)[1] <- "Counts"
      }
      Data = Trendy::getCounts(Data)
    }
    if (!(methods::is(Data, "SummarizedExperiment"))) {
      Data <- data.matrix(Data)
    }
    if (anyNA(Data)) {stop("Data contains at least one 
                           value of NA. Unsure how to proceed.")}
    NSample <- ncol(Data)
    if (is.null(tVectIn)) {
      warning(paste0("No values for parameter tVectIn was given. 
                     Trendy will assume data goes from 1:",NSample))
      tVectIn <- seq_len(NSample)
      names(tVectIn) <- colnames(Data)
    }
    if (is.null(names(tVectIn))) {
      names(tVectIn) <- colnames(Data)
    }
    origMF <- par()
    par(mfrow = origMF$mfrow, mar=c(5,4,2,1))
    
    if (simple == FALSE & showLegend == TRUE) {
      if(legendLocation == 'side') {
        par(mar = c(5,4,2,8), mfrow = origMF$mfrow)
      }
      if(legendLocation == 'bottom') {
        par(oma = c(2,.5, .5, .5), mar=c(5,4,2,1), mfrow = origMF$mfrow)
      }
    }
    if(is.null(spacing.x)) {
      spacing.x <- sd(tVectIn)
    }
    ignoreOUT <- lapply(featureNames, function(x) {
      if (is.null(customTitle)) {
        customTitle = x
      }
      if (showFit==FALSE) {
        plot(tVectIn, Data[x,], pch = 20,
             main = customTitle, ylab = ylab, xlab = xlab)
      }
      if (showFit==TRUE){
        if (is.null(trendyOutData)) {tmp.fit <- fitSegBIC(Data = Data[x,],
                                                          tVectIn=tVectIn)
        }
        if (!is.null(trendyOutData)) {tmp.fit <- trendyOutData[[x]]}
        if (simple==TRUE) {
          plot(tVectIn, Data[x,], pch = 20,
               main = customTitle, ylab = ylab, xlab = xlab)
          lines(tVectIn, tmp.fit$Fitted.Values, lwd = 2)
        } else {
          
          plot(tVectIn, Data[x,], pch = 20, col = "#696969", main = customTitle, 
               ylab = ylab, xlab = xlab, xaxt = 'n')
          axis(1, at=round(seq(min(tVectIn),max(tVectIn),spacing.x)), labels = customLabels.x)
          toplot <- which(!duplicated(tVectIn))
          lines(tVectIn[toplot], tmp.fit$Fitted.Values[toplot],lwd = 3, col="#ededed")
          
          abline(v = tmp.fit$Breakpoints, lty = 2, lwd = 3, col=segColors[1])
          ID <- tmp.fit$Trends
          FIT <- tmp.fit$Fitted.Values
          BKS <- c(0, tmp.fit$Breakpoints, max(tVectIn))
          if (length(BKS) > 3 | (length(BKS) == 3 & !is.na(BKS[2]))) {
            for (i in seq_len(length(BKS) - 1)) {
              toCol <- which(tVectIn <= BKS[i + 1] & tVectIn >= BKS[i])
              while (length(unique(tVectIn[toCol]))<2) {
                if (rev(toCol)[1] < length(tVectIn)) {
                  toCol <- c(toCol, rev(toCol)[1] + 1)
                } else {
                  toCol <- c((toCol[1] - 1), toCol)
                }
              }
              
              forCol <- which(tVectIn <= BKS[i + 1] & tVectIn > BKS[i])
              IDseg <- ID[toCol]
              useCol <- switch(names(which.max(table(ID[forCol]))), 
                               `0` = segColors[3], 
                               `-1` = segColors[4], 
                               `1` = segColors[2])
              tvect.seg <- tVectIn[toCol]
              fit.seg <- FIT[toCol]
              toplot <- which(!duplicated(tvect.seg))
              toplot <- toplot[c(1, length(toplot))]
              lines(tvect.seg[toplot], fit.seg[toplot], lwd = 5, col = useCol)
            }
          } else {
            IDseg <- ID[1]
            useCol <- switch(names(which.max(table(IDseg))), 
                             "0" = segColors[3], 
                             "-1" = segColors[4], 
                             "1" = segColors[2])
            toplot <- which(!duplicated(tVectIn))
            lines(tVectIn[toplot], tmp.fit$Fitted.Values[toplot], lwd = 5, col=useCol)
          }
        } 
        if(simple == FALSE & showLegend == TRUE & legendLocation=='side') {
          LEFT = max(tVectIn) + .05*max(tVectIn)
          TOP = max(Data[x,]) + .01*max(Data[x,])
          BOTTOM = seq(min(Data[x,]), max(Data[x,]), length.out = 10)[8]
          
          legend(LEFT, TOP, c("Breakpoint"), xpd = TRUE, horiz = FALSE,
                 inset = c(.02,0), bty = "n", lty = c(2, 1, 1, 1),
                 lwd = c(3,5,5,5),
                 col = segColors[1], cex = cexLegend)
          
          legend(LEFT,BOTTOM, c("Up", "No change", "Down"), xpd = TRUE,
                 horiz = FALSE,
                 bty = "n", title = "Segment trend:",
                 lty = c(1, 1, 1),
                 lwd = c(5,5,5),
                 col = segColors[2:4], cex = cexLegend)
          
        }
        
      }
    })
    if (simple == FALSE & showLegend == TRUE & legendLocation=='bottom') {
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), 
          mar = c(0, 0, 0, 0), new = TRUE)
      plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
      legend('bottom', xpd = NA, 
             c("Breakpoint", "Up", "No change", "Down"), 
             horiz = TRUE, 
             bty = "n", lty = c(2, 1, 1, 1), 
             lwd = c(3,5,5,5), cex=cexLegend, seg.len = 2, 
             col = segColors)
      
      
    }
    
    par(mar=origMF$mar)
    }