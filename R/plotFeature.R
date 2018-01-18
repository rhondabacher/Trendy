#' @title Plot features of interest

#' @description plot each feature with (or without) the fitted trend.

#' @inheritParams trendy
#' @param featureNames a list of genes or features to plot
#' @param showFit whether to plot the segmented regression 
#'  fitting (default is TRUE)
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
    function(Data, tVectIn = NULL, featureNames, showFit = TRUE, 
        trendyOutData = NULL, condCol = NULL,
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
        plot(tVectIn, Data[x,], pch = 20, col = condColSample, 
            main = x, ylab = ylab, xlab = xlab)
        if (showFit==TRUE){
            if (is.null(trendyOutData)) {tmp.fit <- fitSegBIC(Data = Data[x,],
                tVectIn=tVectIn)
            }
            if (!is.null(trendyOutData)) {tmp.fit <- trendyOutData[[x]]}
            lines(tVectIn, tmp.fit$Fitted.Values, lwd = 2)
        }
    })
}