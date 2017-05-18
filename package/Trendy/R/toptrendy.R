#' @title obtain top genes from trendy results

#' @usage toptrendy(Seg.Data, AdjR.Cut = .5)
#' @description reformats the list output for genes with a given adjusted R^2 cutoff
#' @param Seg.Data output from the trendy function
#' @param AdjR.Cut cutoff for the adjusted r^2. Genes whose adjusted
#'		r^2 is greater than adjR.Cut are called as significant. 

#' @return only significant genes will be included in the output. 
#'		The output is reformatted as:
#'		id.sign direction of each sample; -1: down, 0: no change, 1: up
#'		slp: fitted slopes, slp.sign: sign of fitted slopes, 
#'		slp.pval: p value of each segment, 
#'		bp: estimated breakpoints, fitted: fitted values 
#'		radj: adjusted r value of the model
#'		fit: fit object
#' @examples d1 <- rbind(c(rep(1,50),1:50), c(100:1))
#'      rownames(d1) <- c("g1","g2")
#'      seg.all <- trendy(d1)
#' 		top.seg <- toptrendy(seg.all)
#' @export

toptrendy <- function(Seg.Data, AdjR.Cut = .5) {
	
	
	seg.all.radj <- sapply(Seg.Data, function(i) i$radj)
	names(seg.all.radj) <- names(Seg.Data)
	seg.all.radj.sort <- sort(seg.all.radj, decreasing = T)
	seg.all.highr <- names(seg.all.radj.sort)[which(seg.all.radj.sort > AdjR.Cut)]

	seg.use <- Seg.Data[seg.all.highr]	

	reformat <- sapply(1:7, function(i) sapply(seg.use, function(j) j[[i]]))
	names(reformat) <- names(Seg.Data[[1]])
	return(reformat)
}




