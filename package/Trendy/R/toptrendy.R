#' @title obtain top genes from trendy results

#' @description reformats the list output for genes with a given adjusted R^2 cutoff
#' @param Trendy.Out output from the trendy function
#' @param AdjR2.Cut cutoff for the adjusted R^2. Genes whose adjusted
#'		R^2 is greater than AdjR2.Cut are called as significant. 

#' @return only significant genes will be included in the output. 
#'		The output is reformatted as:
#'		Trend direction of each sample; -1: down, 0: no change, 1: up
#'		Slope: fitted slopes, Slope.Trend: sign of fitted slopes, 
#'		Slope.Pvalue: p value of each segment, 
#'		Breakpoint: estimated breakpoints, Fitted.Values: fitted values 
#'		AdjustedR2: adjusted r value of the model
#'		Fit: fit object
#' @examples d1 <- rbind(c(rep(1,50),1:50), c(100:1))
#'      rownames(d1) <- c("g1","g2")
#'      seg.all <- trendy(d1)
#' 		  top.genes <- topTrendy(seg.all)
#' @export

topTrendy <- function(Trendy.Out, AdjR2.Cut = .5) {
	
	
	seg.all.radj <- sapply(Trendy.Out, function(i) i$AdjustedR2)
	names(seg.all.radj) <- names(Trendy.Out)
	seg.all.radj.sort <- sort(seg.all.radj, decreasing = T)
	seg.all.highr <- names(seg.all.radj.sort)[which(seg.all.radj.sort > AdjR2.Cut)]

	if (length(seg.all.highr) == 0) {
	  stop("No gene pass the adjusted R^2 cutoff filter!")
	}
	
	seg.use <- Trendy.Out[seg.all.highr]	

	reformat <- sapply(1:7, function(i) { 
	            sapply(seg.use, function(j) j[[i]], simplify = FALSE)}, 
	            simplify = FALSE)
	names(reformat) <- names(Trendy.Out[[1]])
	
	reformat[[1]] <- do.call(rbind, reformat[[1]])
	
	reformat[[2]] <- formatFunc(reformat[[2]])
	reformat[[3]] <- formatFunc(reformat[[3]])
	reformat[[4]] <- formatFunc(reformat[[4]])
	reformat[[5]] <- formatFunc(reformat[[5]])
	
	reformat[[6]] <- do.call(rbind, reformat[[6]])
	
	reformat[[7]] <- do.call(c, reformat[[7]])
	
	return(reformat)
}



#' @title helper function to format results

#' @description helper function to format result
#' @param IN the object to be formatted
#' @return a formated matrix of results
#' @author Rhonda Bacher

formatFunc <- function(IN) {
  
  getL <- sapply(IN, length)
  max.len <- max(getL)
  getNames <- names(IN[which.max(getL)[1]][[1]])
  corrected.list <- lapply(IN, function(x) {c(x, rep(NA, max.len - length(x)))})
  all.summary <- do.call(rbind, corrected.list)
  colnames(all.summary) <- getNames

  return(all.summary)
}




