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
#' @param condCol color for each condition, names of this vector should be
#'  condition names; if it is NULL (default), no legend will be generated
#' @param condColSample each sample's color. The vector's length should match
#'  number of samples. if it is NULL, expression will be shown in 
#'  black (default).
#' @param condPrefix name to be used in the condition legend (Default is Day, 
#'  but could be hour, minute, etc.).
#' @param xlab x-axis name
#' @param ylab y-axis name
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
				showLegend = TRUE, trendyOutData = NULL, condCol = NULL,
        condColSample = NULL, 
        condPrefix = "Day", xlab = "Time", 
        ylab = "Normalized Expression") 

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
    if (!is.null(condCol)){
        plot(1, type = "n", axes = FALSE, xlab = "", ylab="")
        legend("top", paste(condPrefix, names(condCol)), 
        col = condCol, ncol = 4, lwd = 2)
    }
    if (is.null(condCol)) {condColSample = "black"}
    ignoreOUT <- lapply(featureNames, function(x) {
        
        if (showFit==TRUE){
            if (is.null(trendyOutData)) {tmp.fit <- fitSegBIC(Data = Data[x,],
                tVectIn=tVectIn)
            }
            if (!is.null(trendyOutData)) {tmp.fit <- trendyOutData[[x]]}
						if (simple==TRUE) {
			        plot(tVectIn, Data[x,], pch = 20, col = condColSample, 
			            main = x, ylab = ylab, xlab = xlab)
							lines(tVectIn, tmp.fit$Fitted.Values, lwd = 2)
						} else {
							if(showLegend == TRUE) {par(mfrow=c(1,1), cex=1.5, cex.lab=1, 
									cex.axis=1, cex.main=1.1, 
							        mar=c(5,5,2,2), oma=c(0,.1,.1,6))}
							plot(tVectIn, Data[x,], pch = 20, col = "#696969", main = x, 
										ylab = ylab, xlab = xlab)
							lines(tVectIn, tmp.fit$Fitted.Values, lwd = 3, col="#ededed")
							abline(v = tmp.fit$Breakpoints, lty = 2, lwd = 3, col="chartreuse3")
							ID <- tmp.fit$Trends
							FIT <- tmp.fit$Fitted.Values
							BKS <- c(0, tmp.fit$Breakpoints, max(tVectIn))
							if (length(BKS) > 3 | (length(BKS) == 3 & !is.na(BKS[2]))) {
							   for (i in seq_len(length(BKS) - 1)) {
							       toCol <- which(tVectIn <= BKS[i+1] & tVectIn >= BKS[i])
							       IDseg <- ID[toCol]
							       useCol <- switch(names(which.max(table(IDseg))), 
							       "0" = "black", 
							       "-1" = "cornflowerblue", 
							       "1" = "coral1")
							       lines(tVectIn[toCol], FIT[toCol], lwd = 5, col=useCol)
							   }
							   if(showLegend == TRUE) {
									 par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), 
							     		mar = c(0, 0, 4, 0), new = TRUE)
								   plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
								   legend("topright", c("Breakpoint"), xpd = TRUE, horiz = FALSE,
								       inset = c(.02,0), bty = "n", lty = c(2, 1, 1, 1),
								       lwd = c(3,5,5,5),
								       col = c("chartreuse3" ), cex = 1)
								   par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
								       mar = c(0, 0, 7.5, 0),new = TRUE)
								   plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
								   legend("topright", c("Up", "No change", "Down"), xpd = TRUE,
								   horiz = FALSE,
								   inset = c(.02,0), bty = "n", title = "Segment trend:",
								   lty = c(1, 1, 1),
								   lwd = c(5,5,5),
								   col = c("coral1", "black","cornflowerblue"), cex = 1)
								 }
							
						}
        }
    	}
		})
}