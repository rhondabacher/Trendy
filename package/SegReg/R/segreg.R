#' @title segmented regression on a set of genes
#' @usage segreg(data, meancut=10,maxk=3, t.vect=NULL, min.num.in.seg=5, pvalcut=.1,
#'                cutdiff=.1, num.try=100,keepfit=FALSE, forceRadj = FALSE, NCores=NULL, makeShinyObject=TRUE)
#' @param indata matrix of normalized expression measures. Rows are genes and columns are samples. The data matrix is expected to be normalized.
#' @param maxk max number of breakpoints to consider.
#' @param t.vect a numerical vector indicates time points. If it is
#' NULL (default), the time will be assumed to be 1:N in which N is number of equally spaced samples.
#' @param meancut default is 10. Genes whose mean is less than meancut will not be considered.
#' @param min.num.in.seg min number of samples within a segment.
#' @param pvalcut p-value cutoff. If the p-value of a segment is greater than pvalcut,
#' this segment will be called as 'no change'.
#' @param cutdiff default is 0.1. If the difference between r_sq from the k+1 breakpoint model
#' and the r_sq from the k breakpoint model is less than cutdiff, the optimal number
#' of breakpoint will be set as k instead of k+1.
#' @param num.try number of different seeds to try. If all num.try runs fails,
#' linear regression results will be returned (which represents one segment case).
#' @param keepfit whether to keep the fitted object.
#' @param forceRadj whether to use adjusted r_sq instead of r_sq in model comparisons.
#' @param NCores number of cores to use, default in detectCores() - 1.
#' @param makeShinyObject whether to output the object required for the SegReg shiny app.
#' @param outfile name and path of output of the object required for the shiny app.

#' @description For each gene/feature, the function will fit maxk + 1 segmented regression models containing 0 - maxk breakpoints
#' (1 - (maxk+1)) segments. The model selection step selects the model with k breakpoints as that having both the highest r_sq value in addition to the difference in r_sq's from the k-1 and k breakpoint models being larger than cutdiff, i.e. r_sq(k) - r_sq(k-1) > cutdiff.


#' @return id.sign: direction of each sample based on pvalcut; -1: down, 0: no change, 1: up
#' slp: fitted slopes, slp.sign: sign of fitted slopes, slp.pval: p-value of each segment, 
#' bp: estimated breakpoints, fitted: fitted values radj: adjusted r_sq value of the model
#' fit: fit object
#' @examples d1 <- rbind(c(rep(1,50),1:50), c(100:1))
#' rownames(d1) <- c("g1","g2")
#' segreg(d1)
#' @author  Rhonda Bacher & Ning Leng
###################
# fit segmented regression for each gene
# the optimal k is the last one whose Radj is > last one + 0.1
# search from 1 to 5 breakpoints
###################
#' @export

segreg <- function(indata, meancut = 10, maxk = 3, t.vect = NULL, min.num.in.seg = 5, pvalcut = .1,
                    cutdiff = .1, num.try = 100, keepfit = FALSE, 
					forceRadj = FALSE, NCores = NULL, makeShinyObject = TRUE, outfile = NULL) {


 if(anyNA(indata)) {stop("Data contains at least one value of NA. Unsure how to proceed.")}
 if (is.null(rownames(indata))) {stop("Must supply gene/feature(row) names!")}
 if (is.null(colnames(indata))) {stop("Must supply sample(column) names!")} 
 if (is.null(NCores)) {NCores <- max(1, detectCores() - 1)}
 if (is.null(outfile)) {outfile="rdata_object_forShiny"} 
	 
 data.gt10 <- indata[which(rowMeans(indata) >= meancut),]

	nsample <- ncol(data.gt10)
	if (nsample < (maxk + 1) * min.num.in.seg) {
		maxk <- floor(nsample / min.num.in.seg) - 1
		message("Number of samples (", nsample, ") is less than 
		[# segments] * [min number of samples in a segment]. maxk has been
		set to", maxk)
		}

	NCores <- detectCores() - 1
	
	if (.Platform$OS.type == "windows") {
		NCores = 1
	}	

seg.all <- mclapply(1:nrow(data.gt10), function(x) {
		fit.seg(indata = data.gt10[rownames(data.gt10)[x],], 
			maxk=maxk, t.vect=t.vect, 
			min.num.in.seg=min.num.in.seg, pvalcut=pvalcut, 
			cutdiff=cutdiff, num.try=num.try, 
			keepfit=keepfit, forceRadj = FALSE)}, mc.cores=NCores)

			
names(seg.all) <- rownames(data.gt10)

t.vect = t.vect
orig.data = data.gt10
seg.object = seg.all
save(t.vect, orig.data, seg.object, file=paste0(outfile,".RData"))

currDir <- getwd()
message("SegReg has finished running and the the output object for shiny has been output to ", currDir, "/", outfile,".Rdata")
return(seg.all)
}
