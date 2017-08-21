#' @title Function to format results for saving.

#' @description format data to save.
#' @param Top.Trendy results from toptrendy() function
#' @param Feature.Names an optional vector of genes to use.
#' @return The function will reformat the output from Trendy so
#' that it can be easily save as a .txt or .csv file.
#' If Feature.Names is supplied then only the information
#' for those genes is returned.
#' @author Rhonda Bacher
#' @export

formatresults <- function(Top.Trendy, Feature.Names = NULL){
  
  if(is.null(Feature.Names)) {
    Feature.Names = names(Top.Trendy$radj)
  }
  all.slopes <- formatfunc(Top.Trendy$slp, "slope")
  all.bp <- formatfunc(Top.Trendy$bp, "breakpoint")
  all.radj <- Top.Trendy$radj

  ALL_DATA <- data.frame(feature=Feature.Names, all.slopes[Feature.Names,], 
                         all.bp[Feature.Names,], adjR2= all.radj[Feature.Names])
return(ALL_DATA)
}

#' @title helper function to format results

#' @description helper function to format result
#' @param IN the object to be formatted
#' @param NAME the name of the object
#' @return a formated matrix of results
#' @author Rhonda Bacher

formatfunc <- function(IN, NAME) {
  max.len <- max(sapply(IN, length))
  corrected.list <- lapply(IN, function(x) {c(x, rep(NA, max.len - length(x)))})
  all.summary <- do.call(rbind, corrected.list)
  colnames(all.summary) <- paste0(NAME, seq_len(max.len))

  return(all.summary)
}
