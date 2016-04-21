#' @title obtain top genes from segreg results
#' @usage topsegreg{res.segreg, r.cut=.5}
#' @param res.segreg output from the segreg function
#' @param r.cut cutoff for the adjusted r square. genes whose adjusted
#' r square is greater than r.cut are called as significant. 
#' @return only significant genes will be included in the output
#' id.sign direction of each sample; -1: down, 0: no change, 1: up
#' slp: fitted slopes, slp.sign: sign of fitted slopes, slp.pval: p value of each segment, 
#' bp: estimated breakpoints, fitted: fitted values radj: adjusted r value of the model
#' fit: fit object

topsegreg <- function(res.segreg, r.cut=.5){
seg.all <- res.segreg	

seg.all.radj <- sapply(seg.all, function(i)i$radj)
names(seg.all.radj) <- names(seg.all)
seg.all.radj.sort <- sort(seg.all.radj, decreasing=T)
seg.all.highr <- names(seg.all.radj.sort)[which(seg.all.radj.sort> r.cut)]

seg.use <- seg.all[seg.all.highr]	

reformat <- sapply(1:7, function(i)sapply(seg.use, function(j)j[[i]]))
names(reformat) <- names(seg.all[[1]])
reformat	
}




