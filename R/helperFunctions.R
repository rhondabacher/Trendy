#' @title getCounts
#' 
#' @usage getCounts(DATA)
#'
#' @param DATA An object of class \code{SummarizedExperiment} that contains 
#' expression data and metadata
#' 
#' @description Convenient helper function to extract the normalized 
#'  expression matrix from the SummarizedExperiment
#' 
#' @return A \code{matrix} which contains the expression data
#'  where genes/features are in rows and samples are in columns
#'
#' @export
#'
#' @importFrom SummarizedExperiment assays
#' @examples 
#'  m1 <- rbind(c(rep(1,50),seq_len(50)), rev(seq_len(100)))
#'  ExampleData <- 
#'  SummarizedExperiment::SummarizedExperiment(assays=list("Counts"=m1))
#'  myData <- getCounts(ExampleData)
    
getCounts <- function(DATA) {
    return(SummarizedExperiment::assays(DATA)[["Counts"]])
}
    
    
#' @title Function to format results for saving.
    
#' @description format data from Trendy which can be saved for later use.
#' @param topTrendyData results from topTrendy() function
#' @param featureNames an optional vector of features 
#'  (if only interested in outputting a subset of features/genes).
#' @return The function will reformat the output from Trendy so
#' that it can be easily save as a .txt or .csv file.
#' If featureNames is supplied then only the information
#' for those features/genes is returned.
#' @author Rhonda Bacher
#' @export
#' @examples data(trendyExampleData)
#'  myTrends <- trendy(Data=trendyExampleData[seq_len(2),])
#'  myTrends <- results(myTrends)
#'  topTrendyRes <- topTrendy(myTrends)
#'  resToSave <- formatResults(topTrendyRes)

formatResults <- function(topTrendyData, featureNames = NULL) {
    if (is.null(featureNames)) {
        featureNames = names(topTrendyData$AdjustedR2)
    }
    ALL_DATA <- data.frame(Feature = featureNames, 
        topTrendyData$Segment.Slopes[featureNames,], 
        topTrendyData$Segment.Trends[featureNames,],
        topTrendyData$Segment.Pvalues[featureNames,],
        topTrendyData$Breakpoints[featureNames,], 
        AdjustedR2 = topTrendyData$AdjustedR2[featureNames],
        topTrendyData$Trends[featureNames,])
    return(ALL_DATA)
}




#' @title internal helper function to format results

#' @description helper function to format result
#' @param IN the object to be formatted
#' @return a formated matrix of results
#' @author Rhonda Bacher
#' @export

formatFunc <- function(IN) {
    
    getL <- vapply(IN, length, integer(1))
    maxLen <- max(getL)
    getNames <- names(IN[which.max(getL)[1]][[1]])
    corrected.list <- lapply(IN, function(x) {c(x, rep(NA,maxLen-length(x)))})
    allSummary <- do.call(rbind, corrected.list)
    colnames(allSummary) <- getNames
    
    return(allSummary)
}

#' @title results
#'   
#' @param DATA An object of class \code{SummarizedExperiment} that contains 
#' normalized expression and other metadata, and the output of the
#' \code{Trendy} function.
#' 
#' @param type A character variable specifying which output is desired, 
#'  with possible values "TrendyFits".
#'  By default results() will
#'  return type="TrendyFits", which is the matrix of normalized 
#'  counts from SCnorm.
#'    
#' @description Convenient helper function to extract the results of running 
#' Trendy. Results data.frames/matrices are stored in the 
#' \code{metadata} slot and can also be accessed without the help of this 
#' convenience function by calling \code{metadata()}.
#'
#' @return A \code{data.frame} containing output as detailed in the
#'  description of the \code{type} input parameter
#'
#' @export
#'
#' @importFrom S4Vectors metadata
#' @examples 
#'  data(trendyExampleData)
#'  Conditions = rep(c(1), each= 90)
#'  trendyOut <- trendy(Data=trendyExampleData[seq_len(2),])
#'  trendyResults <- results(trendyOut)

results <- function(DATA, type=c("TrendyFits")) {
    type <- match.arg(type)
    return(S4Vectors::metadata(DATA)[[type]])
}


#' @title break point fits
#' @param J number of breakpoints in the model


.breakpointFit <- function(J, tVectIn, lmLinear, numTry) {
  lastT <- tVectIn[length(tVectIn)]
  firstT <- tVectIn[1]
  useSeed <- 1
  lmseg.try <- suppressMessages(try(segmented(lmLinear, seg.Z = ~tVectIn,
                                              psi = round(seq(firstT, lastT, length.out = J + 2)[seq_len(J+1)[-1]]), 
                                              control = seg.control(seed = useSeed)), silent = TRUE))
  useSeed2 <- useSeed
  while("try-error" %in% class(lmseg.try) & useSeed2 <= numTry) {
    useSeed2 <- useSeed2 + 1
    lmseg.try <- suppressMessages(try(segmented(lmLinear, 
                                                seg.Z = ~tVectIn, psi = round(seq(firstT, lastT, 
                                                                                  length.out = J + 2)[seq_len(J+1)[-1]]), 
                                                control = seg.control(seed = useSeed2)), silent = TRUE))
    }
  lmseg.try 
}