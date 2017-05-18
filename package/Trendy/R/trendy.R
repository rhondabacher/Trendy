#' @title segmented regression on a set of genes

#' @usage trendy(Data = NULL, Mean.Cut=10, Max.K=3, T.Vect=NULL,
#'                Min.Num.In.Seg=5, Pval.Cut=.1, Cut.Diff=.1, Save.Object=FALSE,
#'                File.Name=NULL, Num.Try=100, Keep.Fit=FALSE, Force.Radj = FALSE,
#'                N.Cores=NULL)

#' @param Data matrix of normalized expression measurements. Rows are genes 
#'      and columns are samples.
#' @param Mean.Cut genes whose mean is less than MeanCut will not be 
#'      considered, default is 10.
#' @param Max.K maximum number of breakpoints to consider. For each gene, the
#'      trendy will fit Max.K + 1 models containing 0 -> Max.K breakpoints
#'      (1 -> (Max.K+1) segments). The model with the highest adjusted r^2 
#'      value will be selected.
#' @param T.Vect a numerical vector indicating the time points. If it is NULL
#'      (default), then the time will be assumed to be equaly spaced from 
#'      1:N (N is number of samples).
#' @param Min.Num.In.Seg minimum number of samples required to be within 
#'      a segment.
#' @param Pval.Cut p-value cutoff. If the p-value of a segment is greater than
#'      PvalCut, then the segment will be called as 'no change'.
#' @param Cut.Diff If the difference between the r^2 from the k + 1 breakpoint
#'      model and the r^2 from the k breakpoint model is less than CutFiff, 
#'      then the optimal number of breakpoints is set as k instead of k + 1.
#' @param Save.Object If TRUE then the trendy object produced will be saved 
#'      to use in the shiny app (default is FALSE).
#' @param File.Name the file name (and file path) to save the Trendy object, 
#'      only used if SaveObject=TRUE.
#' @param Num.Try the number of different seeds to try. If all NumTry runs 
#'      fail, then linear regression (no breakpoints, one segment) model will
#'      be returned.
#' @param Keep.Fit whether report the fitted object (default is FALSE).
#' @param Force.Radj whether to use adjusted r^2 in CutDiff evaluation 
#'      instead of r^2.
#' @param N.Cores number of cores to use, default is detectCores() - 1.

#' @description Segmented regression models are fit for each genes. 
#'  The number of model fits is 1 -> Max.K. 

#' @return id.sign: direction of each sample; -1: down, 0: no change, 1: up
#'      slp: fitted slopes, slp.sign: sign of fitted slopes, 
#'      slp.pval: p value of each segment, 
#'      bp: estimated breakpoints, fitted: fitted values 
#'      radj: adjusted r value of the model
#'      fit: fit object
#' @author Ning Leng and Rhonda Bacher

#' @export
#' @importFrom parallel detectCores mclapply
#' @import stats
#' @import segmented
#' @import stats

#' @examples d1 <- rbind(c(rep(1,50),1:50), c(100:1))
#'      rownames(d1) <- c("g1","g2")
#'      trendy(d1)

trendy <- function(Data = NULL, Mean.Cut=10, Max.K=3, T.Vect=NULL,
                Min.Num.In.Seg=5, Pval.Cut=.1, Cut.Diff=.1, Save.Object=FALSE,
                File.Name=NULL, Num.Try=100, Keep.Fit=FALSE, Force.Radj = FALSE,
                N.Cores=NULL) {

  ## Checks
    if (anyNA(Data)) {stop("Data contains at least one value of NA. Unsure how 
            to proceed.")}
    if (is.null(rownames(Data))) {stop("Must supply feature/gene/row names!")}
    if (is.null(colnames(Data))) {stop("Must supply sample/cell names!")}
    if (is.null(N.Cores)) {N.Cores <- max(1, detectCores() - 1)}
    if (.Platform$OS.type == "windows") {
    	N.Cores = 1
    }


    Data.MeanFiltered <- Data[which(rowMeans(Data) >= Mean.Cut),]


    nsample <- ncol(Data)
    if (nsample < (Max.K + 1) * Min.Num.In.Seg) {
    	Max.K <- floor(nsample / Min.Num.In.Seg) - 1
    	message("Number of samples (", nsample, ") is less than 
    	[# segments] * [min number of samples in a segment]. Max.K has been
    	set to", Max.K)}
    if (Max.K < 1) {stop("Invalid value for Max.K. Adjust Min.Num.In.Seg setting
        in order to run Trendy.")}
    
    Features = rownames(Data.MeanFiltered)
    Seg.All <- mclapply(1:nrow(Data.MeanFiltered), function(x) {
    	fit.seg(Data = Data.MeanFiltered[Features[x],], Max.K = Max.K, 
                T.Vect = T.Vect, Min.Num.In.Seg = Min.Num.In.Seg, 
                Pval.Cut = Pval.Cut, Cut.Diff = Cut.Diff, Num.Try = Num.Try, 
                Keep.Fit = Keep.Fit, Force.Radj = FALSE)}, mc.cores = N.Cores)

			
    names(Seg.All) <- Features


    if(Save.Object == TRUE) {
    	if(is.null(T.Vect)) {
    		T.Vect = 1:ncol(Data)
    	}
    	if(is.null(File.Name)){
    		File.Name="MyData"
    	}
    	Orig.Data = Data
    	Seg.Object = Seg.All
    	save(Seg.Object, Orig.Data, T.Vect, file = paste0(File.Name, "_trendyForShiny.RData"))}
    return(Seg.All)
}
