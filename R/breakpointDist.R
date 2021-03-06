#' @title Distribution of breakpoints

#' @description calculates number of breakpoints at each time.
#' @param topTrendyData results from topTrendy() function
#' @param NDigits how many digits to be used when rounding 
#'  (default is 0 (return integers))
#' @return The function takes significant genes called from 
#' the topTrendyData() function. For any time point, this 
#' function calculates how many genes
#' have a breakpoint at this time point. The output is the
#' numbers of genes sorted by time point.

#' @author Ning Leng
#' @export
#' @examples 
#'  m1 <- matrix(c(c(rnorm(50,5,1),sort(rnorm(50, 15, 5))), rnorm(100, 50,10)), 2, 100, TRUE)
#'  rownames(m1) <- c("g1","g2")
#'  colnames(m1) <- paste0("time", seq_len(100))
#'  myTrends <- results(trendy(m1))
#'  topGenes <- topTrendy(myTrends)
#'  bpDist <- breakpointDist(topGenes)

breakpointDist <- 
    function(topTrendyData, NDigits = 0)
{
    tmpbp <- topTrendyData$Breakpoints
    tab <- table(round(unlist(tmpbp), NDigits))
    tabout <- tab[order(as.numeric(names(tab)))]
    gap <- 1/(10^NDigits)
    tabname.num <- as.numeric(names(tabout))
    namesuse <- as.character(seq(min(tabname.num), max(tabname.num), gap))
    tabnames <- names(tabout)
    outv <- rep(0, length(namesuse))
    names(outv) <- namesuse
    outv[tabnames] <- tabout
    return(outv)
}
