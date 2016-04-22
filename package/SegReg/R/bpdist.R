#' @title Distribution of breakpoints
#' @usage bpdist(topres, ndigits=0)
#' @param topres results from topsegreg() function
#' @param ndigits how many digits to be used when rounding. Default is 0 (return integers).
#' @return The function takes significant genes called from 
#' the topsegreg() function. For any time point, this function calculates how many genes
#' have breakpoint at this time point. The output is the numbers of genes sorted
#' by time point
#' @author Ning Leng

bpdist <- function(topres, ndigits=0){
tmpbp <- topres$bp
tab <- table(round(unlist(tmpbp), ndigits))
tabout <- tab[order(as.numeric(names(tab)))]
gap <- 1/(10^ndigits)
tabname.num <- as.numeric(names(tabout))
namesuse <- as.character(seq(min(tabname.num),max(tabname.num),gap))
tabnames <- names(tabout)
outv <- rep(0,length(namesuse))
names(outv) <- namesuse
outv[tabnames] <- tabout
outv
}
