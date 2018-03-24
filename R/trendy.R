#' @title Trendy

#' @param Data a matrix of normalized expression measurements. Rows are 
#'  genes/features and columns are samples.
#' @param tVectIn a numerical vector indicating the time-points or the 
#'  order of samples. If it is NULL (default), then the time/order 
#'  will be assumed to be equaly spaced from 1:N (N is number of samples).
#' @param saveObject if TRUE then the trendy object produced will be saved 
#'  to use in the Shiny app (default is FALSE).
#' @param fileName the file name (and file path) to save the Trendy object, 
#'  only used if saveObject=TRUE (default name is 
#'  trendyOutputForShiny.RData).
#' @param meanCut genes whose mean is less than MeanCut will not be 
#'  considered, default is 10.
#' @param maxK maximum number of breakpoints to consider. For each gene,
#'  trendy will fit maxK + 1 models containing 0 -> maxK breakpoints
#'  (1 -> (maxK + 1) segments). The model with the lowest BIC 
#'  value will be selected (unless forceRsq = TRUE, see below).
#' @param minNumInSeg minimum number of samples required to be within 
#'  a segment. If a breakpoint model has a segment with fewer than 
#'  minNumInSeg point in any segment, then the model is not 
#'  considered valid.
#' @param pvalCut p-value cutoff. If the p-value of a segment is greater than
#'  PvalCut, then the segment will be called as 'no change'.
#' @param numTry the number of different seeds to try. If all numTry runs 
#'  fail, then the linear regression (no breakpoints, one segment) model 
#'  will be returned.
#' @param keepFit whether to report the fitted object (default is FALSE).
#' @param featureNames optional parameter to specify an explicit subset of 
#'  features/genes to fit the segmented regression model to.
#' @param NCores number of cores to use, default is detectCores() - 1.

#' @description Segmented regression models are fit for each gene. 
#'  The number of model fits is 1 -> maxK. 

#' @return Trend: direction of each sample; -1: down, 0: no change, 1: up
#'  Slope: fitted slopes, Slope.Trend: sign of fitted slopes, 
#'  Slope.Pvalue: p value of each segment, Breakpoint: estimated breakpoints, 
#'  Fitted.Values: fitted values AdjustedR2: adjusted R squared 
#'      value of the model
#'  Fit: fit object
#' @author Ning Leng and Rhonda Bacher

#' @export
#' @importFrom BiocParallel bplapply register MulticoreParam bpparam SnowParam
#' @importFrom parallel detectCores
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment SummarizedExperiment assayNames 
#'  assays colData
#' @importFrom methods is
#' @import stats
#' @import segmented

#' @examples m1 <- rbind(c(rep(1,50),seq_len(50)), rev(seq_len(100)))
#'  rownames(m1) <- c("g1","g2")
#'  colnames(m1) <- paste0("time", seq_len(100))
#'  myTrends <- trendy(m1)

trendy <- 
    function(Data = NULL, tVectIn = NULL, saveObject = FALSE, fileName = NULL,
            meanCut = 10, maxK = 3, minNumInSeg = 5, pvalCut = .1,  
            numTry = 100, keepFit = FALSE, NCores = NULL, featureNames = NULL) 
{
    # Checks
    if (methods::is(Data, "SummarizedExperiment")) {
        if (is.null(SummarizedExperiment::assayNames(Data)) || 
        SummarizedExperiment::assayNames(Data)[1] != "Counts") {
            message("Renaming the first element in assays(Data) to 'Counts'")
            SummarizedExperiment::assayNames(Data)[1] <- "Counts"
            if (is.null(colnames(Trendy::getCounts(Data)))) {
                stop("Must supply sample/cell names!")
            }
        }
    }
    if (!(methods::is(Data, "SummarizedExperiment"))) {
        Data <- data.matrix(Data)
        Data<-SummarizedExperiment::SummarizedExperiment(assays=
            list("Counts"=Data))
    }
    if (anyNA(Trendy::getCounts(Data))) {stop("Data contains at least one 
        value of NA. Unsure how to proceed.")}
    if (is.null(rownames(Trendy::getCounts(Data)))) {
        stop("Must supply feature/gene/row names!")
    }
    if (is.null(colnames(Trendy::getCounts(Data)))) {
        stop("Must supply sample/column names!")
    }
    NSample <- ncol(Trendy::getCounts(Data))
    if (is.null(tVectIn)) {
        warning(paste0("No values for parameter tVectIn were given. 
        Trendy will assume data goes from 1:",NSample))
        tVectIn <- seq_len(NSample)
        names(tVectIn) <- colnames(Trendy::getCounts(Data))
    }
    if (is.null(names(tVectIn))) {
        names(tVectIn) <- colnames(Trendy::getCounts(Data))
    }
    if (is.null(NCores)) {NCores <- max(1, parallel::detectCores() - 1)}
    if (.Platform$OS.type == "windows") {
        param = SnowParam(workers=NCores)
    }

	param = MulticoreParam(workers=NCores)
	
	BiocParallel::register(BPPARAM = param)
    
    if (is.null(featureNames)) {
        featureNames <- rownames(Trendy::getCounts(Data))
    }
    
    SummarizedExperiment::assays(Data)[["Counts"]]<-
        Trendy::getCounts(Data)[featureNames,]
        
    if (length(featureNames) == 1) {

        if (mean(Trendy::getCounts(Data)) >= meanCut) {
            Data.MeanFiltered <- t(data.matrix(Trendy::getCounts(Data)))
            row.names(Data.MeanFiltered) <- featureNames
        } else {
            stop("Gene does not pass the mean cutoff filter!")
        }
    } else {
        toKeep <- which(rowMeans(Trendy::getCounts(Data)) >= meanCut)
        Data.MeanFiltered<-
        Trendy::getCounts(Data)[toKeep,]
        if (sum(rowMeans(Trendy::getCounts(Data)) >= meanCut) == 0)  {
        stop("No genes pass the mean cutoff filter!")
    }
    }
    
    if (NSample < (maxK + 1) * minNumInSeg) {
        maxK <- floor(NSample / minNumInSeg) - 1
        message("Number of samples (", NSample, ") is less than 
        [# segments] * [min number of samples in a segment]. maxK has been
        set to", maxK)
    }
    if (maxK < 1) {
        stop("Invalid value for maxK. Adjust minNumInSeg setting
        in order to run Trendy.")
    }
    
    segAll <- BiocParallel::bplapply(X = seq_len(nrow(Data.MeanFiltered)),
    function(X) { 
        inGene = Data.MeanFiltered[X,]
        fitSegBIC(Data = inGene,
            tVectIn = tVectIn,  maxK = maxK, 
            minNumInSeg = minNumInSeg, 
            pvalCut = pvalCut, numTry = numTry, 
            keepFit = keepFit)
        })
    names(segAll) <- rownames(Data.MeanFiltered)
    
    if (saveObject == TRUE) {
        if (is.null(fileName)){
            fileName <- "trendyForShiny.RData"
        } else {
            fileName <- paste0(fileName, "_trendyForShiny.RData")
        }
        origData <- Trendy::getCounts(Data)
        trendyOut <- segAll
        tVectIn <- tVectIn

        save(trendyOut, origData, tVectIn, file = fileName)}

        S4Vectors::metadata(Data)[["TrendyFits"]] <- segAll
        return(Data)
}