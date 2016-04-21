#' @title Distribution of breakpoints
#' @usage bpdist(topres)
#' @param topres results from topsegreg() function
#' @return The function takes significant genes called from 
#' the topsegreg() function. For any time point, this function calculates how many genes
#' have breakpoint at this time point. The output is the numbers of genes sorted
#' by time point
#' @author Ning Leng

bpdist <- function(topres){
tmpbp <- topres$bp
tab <- table(round(unlist(tmpbp)))
tabout <- tab[order(as.numeric(names(tab)))]
namesuse <- rownames(topres$fitted)
tabnames <- namesuse[as.numeric(names(tabout))]
outv <- rep(0,length(namesuse))
names(outv) <- namesuse
outv[tabnames] <- tabout
outv
}
