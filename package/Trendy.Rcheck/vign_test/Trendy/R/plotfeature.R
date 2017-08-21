#' @title Plot features of interest

#' @description plot each feature with (or without) the fitted trend.

#' @inheritParams trendy
#' @param Feature.Names a list of genes or features to plot
#' @param Seg.Fit whether plot the segmented regression fitting (default is TRUE)
#' @param Seg.Data segmented regression fitting result from running trendy(); i
#'		f Seg.Data is NULL and Seg.Fit = T, then the
#'		segmented regression will be fit for each of the genes and it may take longer to run
#' @param Cond.Col color for each condition, names of this vector should be condition names; 
#'		if it is NULL (default), no legend will be generated
#' @param Cond.Col.Sample each sample's color. The vector's length should match number of samples. 
#'		if it is NULL, expression will be shown in black (default).
#' @param Y.Name y-axis name
#' @param X.Name x-axis name
#' @return plot of gene expression and fitted line

#' @examples d1 <- rbind(c(rep(1,50),1:50), c(100:1))
#'		rownames(d1) <- c("g1","g2")
#'		plotmarker(d1, Feature.Names=c("g1","g2"))
#' @author Ning Leng and Rhonda Bacher
#' @import graphics
#' @import grDevices
#' @export

plotfeature <- function(Data, T.Vect = NULL, Feature.Names, Seg.Fit = TRUE, 
                        Seg.Data = NULL, Cond.Col = NULL, Cond.Col.Sample = NULL,
						            Cond.Prefix = "Day", Y.Name = "Normalized Expression", 
                        X.Name = "Time"){
							
	data.norm <- Data
	if(is.null(T.Vect)) {T.Vect <- 1:ncol(data.norm)}
	
	out <- vector("list",length(Feature.Names))
	names(out) <- Feature.Names
	
	
	if(!is.null(Cond.Col)){
		plot(1, type = "n", axes = FALSE, xlab = "", ylab="")
		legend("top", paste(Cond.Prefix, names(Cond.Col)), col = Cond.Col, ncol = 4, lwd = 2)
	}
	
	if(is.null(Cond.Col)) {Cond.Col.Sample="black"}
		
	for(i in 1:length(Feature.Names)){
		plot(T.Vect, data.norm[Feature.Names[i],], pch = 20, col = Cond.Col.Sample, 
			 main = Feature.Names[i], ylab = Y.Name, xlab = X.Name)
		if(Seg.Fit){
			if(is.null(Seg.Data)) {tmp <- fitseg(data.norm, Feature.Names[i])}
			if(!is.null(Seg.Data)) {tmp <- Seg.Data[[Feature.Names[i]]]}
			lines(T.Vect, tmp$fitted, lwd = 2)
			out[[i]] <- tmp
		}
	}

}


