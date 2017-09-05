#' @title segmented regression on a set of genes

#' @param Data matrix of normalized expression measurements. Rows are genes/features 
#'      and columns are samples.
#' @param T.Vect a numerical vector indicating the time-points (or the order of samples). 
#'      If it is NULL (default), then the time/order will be assumed to be equaly spaced from 
#'      1:N (N is number of samples).
#' @param Save.Object if TRUE then the trendy object produced will be saved 
#'      to use in the shiny app (default is FALSE).
#' @param File.Name the file name (and file path) to save the Trendy object, 
#'      only used if SaveObject=TRUE (default name is Data_trendyForShiny.RData).
#' @param Mean.Cut genes whose mean is less than MeanCut will not be 
#'      considered, default is 10.
#' @param Max.K maximum number of breakpoints to consider. For each gene, 
#'      trendy will fit Max.K + 1 models containing 0 -> Max.K breakpoints
#'      (1 -> (Max.K+1) segments). The model with the lowest BIC 
#'      value will be selected (unless Force.Rsq == TRUE, see below).
#' @param Min.Num.In.Seg minimum number of samples required to be within 
#'      a segment.
#' @param Pval.Cut p-value cutoff. If the p-value of a segment is greater than
#'      PvalCut, then the segment will be called as 'no change'.

#' @param Num.Try the number of different seeds to try. If all NumTry runs 
#'      fail, then the linear regression (no breakpoints, one segment) model will
#'      be returned.
#' @param Keep.Fit whether to report the fitted object (default is FALSE).
#' @param Force.Rsq whether to use adjusted R^2 in model evaluation 
#'      instead of BIC, corresponds to very early versions of Trendy. Requires that
#'      each succesive increase in adjusted R squared by larger than Cut.Diff.
#' @param Cut.Diff Only used if Force.Rsq is set to TRUE. 
#'       If the difference between the adjusted R^2 from the k + 1 breakpoint
#'      model and the adjusted R^2 from the k breakpoint model is less than Cut.Diff, 
#'      then the optimal number of breakpoints is set as k instead of k + 1.
#' @param N.Cores number of cores to use, default is detectCores() - 1.

#' @description Segmented regression models are fit for each gene. 
#'  The number of model fits is 1 -> Max.K. 

#' @return Trend: direction of each sample; -1: down, 0: no change, 1: up
#'		Slope: fitted slopes, Slope.Trend: sign of fitted slopes, 
#'		Slope.Pvalue: p value of each segment, Breakpoint: estimated breakpoints, 
#'		Fitted.Values: fitted values AdjustedR2: adjusted r value of the model
#' 		Fit: fit object
#' @author Ning Leng and Rhonda Bacher

#' @export
#' @importFrom parallel detectCores mclapply
#' @import stats
#' @import segmented
#' @import stats

#' @examples d1 <- rbind(c(rep(1,50),1:50), c(100:1))
#'      rownames(d1) <- c("g1","g2")
#'      trendy(d1)

trendy <- function(Data = NULL, 
                    T.Vect = NULL, 
                    Save.Object = FALSE, 
                    File.Name = NULL, 
                    Mean.Cut = 10, 
                    Max.K = 3, 
                    Min.Num.In.Seg = 5, 
                    Pval.Cut = .1,  
                    Num.Try = 100, 
                    Keep.Fit = FALSE, 
                    N.Cores = NULL, 
                    Feature.Names = NULL,
                    Force.Rsq = FALSE, 
                    Cut.Diff = .1) {

  ## Checks
    if (anyNA(Data)) {stop("Data contains at least one value of NA. Unsure how 
            to proceed.")}
    if (is.null(rownames(Data))) {stop("Must supply feature/gene/row names!")}
    if (is.null(colnames(Data))) {stop("Must supply sample/column names!")}
    if (is.null(N.Cores)) {N.Cores <- max(1, parallel::detectCores() - 1)}
    if (.Platform$OS.type == "windows") {
    	N.Cores = 1
    }
    if (is.null(Feature.Names)) {
      Feature.Names <- rownames(Data)
    }
  
    Data <- Data[Feature.Names, ]
    
    if (length(Feature.Names) == 1) {
      
      if (mean(Data) >= Mean.Cut) {
        Data.MeanFiltered <- t(data.matrix(Data))
        row.names(Data.MeanFiltered) <- Feature.Names
      } else {stop("Gene does not pass the mean cutoff filter!")}
    
    } else {
      Data.MeanFiltered <- Data[which(rowMeans(Data) >= Mean.Cut),]
      if (sum(rowMeans(Data) >= Mean.Cut) == 0)  {
        stop("No genes pass the mean cutoff filter!")
      }
   }
        

  nsample <- ncol(Data.MeanFiltered)
    
  if (nsample < (Max.K + 1) * Min.Num.In.Seg) {
    	Max.K <- floor(nsample / Min.Num.In.Seg) - 1
    	message("Number of samples (", nsample, ") is less than 
    	[# segments] * [min number of samples in a segment]. Max.K has been
    	set to", Max.K)}
  if (Max.K < 1) {stop("Invalid value for Max.K. Adjust Min.Num.In.Seg setting
        in order to run Trendy.")}
    
  Features = rownames(Data.MeanFiltered)
  
  if (is.null(T.Vect)) {
    T.Vect <- 1:ncol(Data)
    names(T.Vect) <- colnames(Data)
  }
  if (is.null(names(T.Vect))) {
    names(T.Vect) <- colnames(Data)
  }
    
  if (Force.Rsq == FALSE) {
      Seg.All <- parallel::mclapply(1:length(Features), function(x) {
      	fitSegBIC(Data = Data.MeanFiltered[Features[x],],
                  T.Vect = T.Vect,  Max.K = Max.K, Min.Num.In.Seg = Min.Num.In.Seg, 
                  Pval.Cut = Pval.Cut, Num.Try = Num.Try, 
                  Keep.Fit = Keep.Fit)}, mc.cores = N.Cores)
    }
		if (Force.Rsq == TRUE) {	
      Seg.All <- parallel::mclapply(1:length(Features), function(x) {
      	fitSegRsq(Data = Data.MeanFiltered[Features[x],],
                  T.Vect = T.Vect,  Max.K = Max.K, Min.Num.In.Seg = Min.Num.In.Seg, 
                  Pval.Cut = Pval.Cut, Cut.Diff = Cut.Diff, Num.Try = Num.Try, 
                  Keep.Fit = Keep.Fit)}, mc.cores = N.Cores)
    }                
    
   names(Seg.All) <- Features

 
   
     if(Save.Object == TRUE) {
    	if(is.null(File.Name)){
    		File.Name="Data"
    	}
    	Orig.Data = Data
    	Trendy.Out = Seg.All
    	save(Trendy.Out, Orig.Data, T.Vect, file = paste0(File.Name, "_trendyForShiny.RData"))}
    return(Seg.All)
}
