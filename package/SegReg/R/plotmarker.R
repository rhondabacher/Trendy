#' @title Plot markers of interest
#' @usage plotmarker(indata, t.vect = NULL, outfile = NULL, listfeatures, savePDF = FALSE, plotFit = TRUE, fittedreg=NULL,
#' cond.col=NULL, cond.col.sample=NULL, cond.prefix="day",
#' yname="normalized expression",par.param = c(3,2),
#' pdfheight=15, pdfwidth=10)
#' @param indata matrix of normalized expression measures. Rows are genes and columns are samples. The data matrix is expected to be normalized.
#' @param t.vect a numerical vector indicates time points. If it is
#'  NULL (default), the time will be assumed to be 1:N in which N is number of equally spaced samples.
#' @param outfile output file name; will be ignored if savePDF=F. Default is feature_expressionScatter.
#' @param listfeatures a list of genes of interest
#' @param savePDF whether output figs to a pdf file.
#' @param plotFit whether to plot the segmented regression fitting.
#' @param fittedreg segmented regression fitting result; if fittedres is NULL and seg=T, segmented regression will be fitted for each of the genes and it may take longer runtime
#' @param cond.col color for each condition, names of this vector should be condition names; if it is null, no legends will be generated
#' @param cond.col.sample each sample's color. The vector's length should match number of samples. If it is NULL, expressions will be shown in black.
#' @param cond.prefix prefix to be appended to condition name
#' @param yname y axis name
#' @param par.param mfrow, default c(3,2)
#' @param pdfheight,pdfwidth height and width for the pdf file
#' @return plot of gene expression and fitted lines.
#' @examples d1 <- rbind(c(rep(1,50),1:50), c(100:1))
#' rownames(d1) <- c("g1","g2")
#' plotmarker(d1, listname=c("g1","g2"), pdf = FALSE)
#' @author Ning Leng & Rhonda Bacher
#################
#' @export

plotmarker <- function(indata, t.vect = NULL, outfile = NULL, listfeatures, savePDF = FALSE, plotFit = TRUE, fittedreg = NULL,
	cond.col = NULL, cond.col.sample = NULL, cond.prefix = "Day", 
	yname = "Normalized Expression", par.param = c(3,2),
	pdfheight = 15, pdfwidth = 10){

		data.norm <- indata
		
		if (is.null(t.vect)) { t.vect <- 1:ncol(data.norm) }
		
		out <- vector("list",length(listfeatures))
		names(out) <- listfeatures
		
		if(is.null(outfile)) {outfile = "feature_expressionScatter"}
		if (savePDF) {pdf(paste0(outfile, ".pdf"), height=pdfheight, width=pdfwidth)}
		par(mfrow=par.param)
		
		if (!is.null(cond.col)) {
		plot(1, type="n", axes=FALSE, xlab="", ylab="")
		legend("top", paste(cond.prefix, names(cond.col)),col=cond.col,ncol=4, lwd=2)
		}
		
		if (is.null(cond.col)) {cond.col.sample="black"}
			
		for(i in 1:length(listfeatures)){
			plot(t.vect, data.norm[listfeatures[i],], pch=20, col=cond.col.sample, 
				 main=listfeatures[i], ylab=yname, xlab="")
			if(plotFit){
				if(is.null(fittedreg)) {tmp <- fit.seg(data,listfeatures[i])}
				if(!is.null(fittedreg)) {tmp <- fittedreg[[listfeatures[i]]]}
				
				lines(t.vect, tmp$fitted, lwd=2)
				out[[i]] <- tmp
			}
			box(lwd=2)
		}
		
		if(savePDF) {dev.off()}
		
	return(out)
		
}


