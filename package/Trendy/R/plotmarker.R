#' @title Plot markers of interest
#' @usage plotmarker(Data, T.Vect = NULL, File.Name = NULL, Feature.Names, PDF = FALSE,
#'						Seg.Fit = TRUE, Seg.Data = NULL, Cond.Col = NULL, 
#'						Cond.Col.Sample = NULL, Cond.Prefix = "Day", 
#'						Y.Name = "Normalized Expression", Par.Param = c(3,2),
#'						PDF.Height = 15, PDF.Width = 10)

#' @description plot each feature with (or without) the fitted trend.

#' @inheritParams trendy
#' @param Feature.Names a list of genes or features to plot
#' @param PDF whether to output plots to a pdf file (default is FALSE)
#' @param Seg.Fit whether plot the segmented regression fitting (default is TRUE)
#' @param Seg.Data segmented regression fitting result from running trendy(); i
#'		f Seg.Data is NULL and Seg.Fit = T, then the
#'		segmented regression will be fit for each of the genes and it may take longer to run
#' @param Cond.Col color for each condition, names of this vector should be condition names; 
#'		if it is NULL (default), no legend will be generated
#' @param Cond.Col.Sample each sample's color. The vector's length should match number of samples. 
#'		if it is NULL, expression will be shown in black (default).
#' @param Cond.Prefix prefix to be appended to condition name (Day is default)
#' @param Y.Name y-axis name
#' @param Par.Param mfrow for plot layout(Default is c(3,2))
#' @param PDF.Height,PDF.Width the height and width for the PDF file
#' @return plot of gene expression and fitted line

#' @examples d1 <- rbind(c(rep(1,50),1:50), c(100:1))
#'		rownames(d1) <- c("g1","g2")
#'		plotmarker(d1, Feature.Names=c("g1","g2"), PDF = FALSE)
#' @author Ning Leng and Rhonda Bacher
#' @import graphics
#' @import grDevices
#' @export

plotmarker <- function(Data, T.Vect = NULL, File.Name = NULL, Feature.Names, PDF = FALSE, 
						Seg.Fit = TRUE, Seg.Data = NULL, Cond.Col = NULL, Cond.Col.Sample = NULL,
						Cond.Prefix = "Day", Y.Name = "Normalized Expression", Par.Param = c(3,2), 
						PDF.Height = 15, PDF.Width = 10){
							
	data.norm <- Data
	if(is.null(T.Vect)) {T.Vect <- 1:ncol(data.norm)}
	
	out <- vector("list",length(Feature.Names))
	names(out) <- Feature.Names
	
	if(PDF) {pdf(File.Name, height = PDF.Height, width = PDF.Width)}
	
	par(mfrow = Par.Param)
	if(!is.null(Cond.Col)){
		plot(1, type = "n", axes = FALSE, xlab = "", ylab="")
		legend("top", paste(Cond.Prefix, names(Cond.Col)), col = Cond.Col, ncol = 4, lwd = 2)
	}
	
	if(is.null(Cond.Col)) {Cond.Col.Sample="black"}
		
	for(i in 1:length(Feature.Names)){
		plot(T.Vect, data.norm[Feature.Names[i],], pch = 20, col = Cond.Col.Sample, 
			 main = Feature.Names[i], ylab = Y.Name, xlab = "")
		if(Seg.Fit){
			if(is.null(Seg.Data)) {tmp <- fit.seg(data.norm, Feature.Names[i])}
			if(!is.null(Seg.Data)) {tmp <- Seg.Data[[Feature.Names[i]]]}
			lines(T.Vect, tmp$fitted, lwd = 2)
			out[[i]] <- tmp
		}
	}
	
	if(PDF) dev.off()
	
	return(out)
}


