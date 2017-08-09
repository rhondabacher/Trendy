#' @title Distribution of breakpoints

#' @description calculates number of breakpoints at each time.
#' @param Top.Trendy results from toptrendy() function
#' @param N.Digits how many digits to be used when rounding 
#'		(default is 0 (return integers))

#' @return The function takes significant genes called from 
#' the toptrendy() function. For any time point, this function calculates how many genes
#' have a breakpoint at this time point. The output is the numbers of genes sorted
#' by time point.

#' @author Ning Leng
#' @export

bpdist <- function(Top.Trendy, N.Digits = 0){
	
	tmpbp <- Top.Trendy$bp
	tab <- table(round(unlist(tmpbp), N.Digits))
	tabout <- tab[order(as.numeric(names(tab)))]
	gap <- 1 / (10^N.Digits)
	tabname.num <- as.numeric(names(tabout))
	namesuse <- as.character(seq(min(tabname.num), max(tabname.num), gap))
	tabnames <- names(tabout)
	outv <- rep(0, length(namesuse))
	names(outv) <- namesuse
	outv[tabnames] <- tabout
	return(outv)
}
