#' @title plot markers of interest
#' @usage plotmarker(data, filename=NULL, listname, pdf=FALSE, seg=TRUE, fittedres=NULL,
#' cond.col=NULL, cond.col.sample=NULL, cond.prefix="day",
#' yname="normalized expression",par.param = c(3,2),
#' pdfheight=15, pdfwidth=10)
#' @param data normalized expression matrix; rows are genes and columns are samples
#' @param filename output file name; will be ignored if pdf=F
#' @param listname a list of genes of interest
#' @param whether output figs to a pdf file
#' @param seg whether plot the segmented regression fitting
#' @param fittedres segmented regression fitting result; if fittedres is NULL and seg=T, segmented regression will be fitted for each of the genes and it may take longer runtime
#' @param cond.col color for each condition, names of this vector should be condition names; if it is null, no legends will be generated
#' @param cond.col.sample each sample's color. The vector's length should match number of samples. if it si NULL, expressions will be shown in black.
#' @param cond.prefix prefix to be appended to condition name
#' @param yname y axis name
#' @param par.param mfrow, default c(3,2)
#' @param pdfheight,pdfwidth height and width for the pdf file
#' @return plot of gene expression and fitted line
#' @examples d1 <- rbind(c(rep(1,50),1:50), c(100:1))
#' rownames(d1) <- c("g1","g2")
#' plotmarker(d1, listname=c("g1","g2"), pdf = FALSE)
#' @author Ning Leng

#################
# check genes
#################
plotmarker <- function(data, x=NULL,filename=NULL, listname, pdf=FALSE, seg=TRUE, fittedres=NULL,
cond.col=NULL, cond.col.sample=NULL, cond.prefix="day", 
yname="normalized expression",par.param = c(3,2),
pdfheight=15, pdfwidth=10){
		data.norm <- data
		if(is.null(x))x <- 1:ncol(data.norm)
		out <- vector("list",length(listname))
		names(out) <- listname
		if(pdf)pdf(filename, height=pdfheight, width=pdfwidth)
		par(mfrow=par.param)
		if(!is.null(cond.col)){
		plot(1, type="n", axes=FALSE, xlab="", ylab="")
		legend("top", paste(cond.prefix, names(cond.col)),col=cond.col,ncol=4, lwd=2)
		}
		if(is.null(cond.col))cond.col.sample="black"
		for(i in 1:length(listname)){
		plot(x, data.norm[listname[i],],pch=20,col=cond.col.sample, 
				 main=listname[i],ylab=yname, xlab="")
		if(seg){
		if(is.null(fittedres))tmp <- fit.seg(data,listname[i])
		if(!is.null(fittedres))tmp <- fittedres[[listname[i]]]
		lines(tmp$fitted,lwd=2)
		out[[i]] <- tmp
		}
		}
		if(pdf)dev.off()
		out
}


