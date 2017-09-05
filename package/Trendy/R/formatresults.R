#' @title Function to format results for saving.

#' @description format data from Trendy which can be saved for later use.
#' @param Top.Trendy results from topTrendy() function
#' @param Feature.Names an optional vector of features 
#'   (if only interested in outputting a subset of features/genes).
#' @return The function will reformat the output from Trendy so
#' that it can be easily save as a .txt or .csv file.
#' If Feature.Names is supplied then only the information
#' for those features/genes is returned.
#' @author Rhonda Bacher
#' @export

formatResults <- function(Top.Trendy, Feature.Names = NULL){
  
  if(is.null(Feature.Names)) {
    Feature.Names = names(Top.Trendy$AdjustedR2)
  }
  
  ALL_DATA <- data.frame(Feature = Feature.Names, 
                          Top.Trendy$Segment.Slopes[Feature.Names,], 
                          Top.Trendy$Segment.Trends[Feature.Names,],
                          Top.Trendy$Segment.Pvalues[Feature.Names,],
                          Top.Trendy$Breakpoints[Feature.Names,], 
                          AdjustedR2 = Top.Trendy$AdjustedR2[Feature.Names],
                          Top.Trendy$Trends[Feature.Names,])
return(ALL_DATA)
}
