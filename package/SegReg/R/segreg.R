#' @title segmented regression on a set of genes
#' @usage segreg(data, meancut=10,maxk=3, t.vect=NULL, min.num.in.seg=5, pvalcut=.1,
#'                cutdiff=.1, num.try=100,keepfit=FALSE)
#' @param data normalized expression measure. Rows are genes and columns are samples. The data matrix is expected to be normalized.
#' @param maxk max number of breakpoints to consider. 
#' @param t.vect a numerical vector indicates time points. If it si
#' NULL (default), the time will be assumed to be 1:N in which N is number of samples.
#' @param meancut default is 10. Genes whose mean is less than meancut will not be considered.
#' For each gene, the function will fit maxk+1 models containing 0->maxk breakpoints
#' (1->(maxk+1)) segments. The model with highest adjusted r value will be selected.
#' @param min.num.in.seg min number of samples within a segment
#' @param pvalcut p value cutoff. If the p value of a segment is greater than pvalcut,
#' this segment will be called as 'no change'
#' @param cutdiff default is 0.1. if the difference between r_adj from the k+1 breakpoint model
#' and the r_adj from the k breakpoint model is less than cutdiff, the optimal number
#' of breakpoint will be set as k instead of k+1
#' @param num.try number of different seeds to try. If all num.try runs fails,
#' linear regression results will be returned (which represents one segment case).
#' @param keepfit whether keep the fitted object
#' @return id.sign: direction of each sample; -1: down, 0: no change, 1: up
#' slp: fitted slopes, slp.sign: sign of fitted slopes, slp.pval: p value of each segment, 
#' bp: estimated breakpoints, fitted: fitted values radj: adjusted r value of the model
#' fit: fit object
#' @examples d1 <- rbind(c(rep(1,50),1:50), c(100:1))
#' rownames(d1) <- c("g1","g2")
#' segreg(d1)
#' @author Ning Leng
###################
# fit segmented regression for each gene
# the optimal k is the last one whose Radj is > last one + 0.1
# search from 1 to 5 breakpoints
###################
segreg <- function(data, meancut=10, maxk=3, t.vect=NULL,min.num.in.seg=5, pvalcut=.1,
                    cutdiff=.1, num.try=100,keepfit=FALSE
										){

data.gt10.raw <- data[which(rowMeans(data)>=meancut),]

#data.gt10 <- t(apply(data.gt10.raw,1,scale))
#rownames(data.gt10) <- rownames(data.gt10.raw)
#colnames(data.gt10) <- colnames(data.gt10.raw)

data.gt10 <- data.gt10.raw

nsample <- ncol(data)
if(nsample < (maxk+1)*min.num.in.seg){
	maxk <- floor(nsample/min.num.in.seg)-1
	message("number of samples (", nsample, ") is less than 
	[# segments] * [min number of samples in a segment]. maxk has been
	set to", maxk)
	}

seg.all <- sapply(1:nrow(data.gt10),function(i){
aa=fit.seg(data.gt10, rownames(data.gt10)[i], 
				maxk=maxk, t.vect=t.vect,min.num.in.seg=min.num.in.seg, pvalcut=pvalcut,
                    cutdiff=cutdiff, num.try=num.try,keepfit=keepfit)
},simplify=F)
names(seg.all) <-  rownames(data.gt10)
seg.all
}
