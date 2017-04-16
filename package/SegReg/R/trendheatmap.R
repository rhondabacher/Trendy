#' @title Draw heatmap of gene expression trends
#' @usage trendheatmap(topres, showplot)
#' @param topres results from topsegreg() function
#' @param showplot whether to show the heatmap.
#' @param savePDF save heatmap as a PDF.
#' @param outfile output file name; will be ignored if savePDF=F. Default is trend_heatmap.
#' @param pdfheight,pdfwidth height and width for the pdf file
#' @param sampnames names of samples for the columns of heatmap. Default is null, column names are not shown.
#' @return The function takes significant genes called from 
#' the topsegreg() function. These genes are further grouped into three
#' groups : Up, Down or No Change (Same) in the first segment. Within each group,
#' the genes are sorted by their first break point. The heatmap shows expression trends of these 
#' three groups of genes. In the heatmap, red/blue/black represents up/down/nochange.
#' @author Ning Leng & Rhonda Bacher
#' @export


trendheatmap <- function(topres, showplot=TRUE, savePDF=FALSE, outfile="trend_heatmap", pdfheight= 15, pdfwidth= 10, sampnames=NULL){

	bks.all <- topres$bp
	bks.first <- sapply(bks.all, function(i)i[1])
	names(bks.first) <- colnames(topres[[1]])
	bks.first.sort <- sort(bks.first)

	seg.all.id.highr <- t(topres[[1]])

	bks.sign.sort.list <- list( firstup=sort(bks.first[rownames(seg.all.id.highr)[which(seg.all.id.highr[,1]==1)]]),
								firstdown=sort(bks.first[rownames(seg.all.id.highr)[which(seg.all.id.highr[,1]==-1)]]),
								firstnochange=sort(bks.first[rownames(seg.all.id.highr)[which(seg.all.id.highr[,1]==0)]])
								)

	bks.sign.sort <- c(bks.sign.sort.list[[1]],bks.sign.sort.list[[2]],bks.sign.sort.list[[3]])

	seg.all.id.highr.sign.sort <- seg.all.id.highr[names(bks.sign.sort),]
	if(!is.null(sampnames)) {colnames(seg.all.id.highr.sign.sort) <- sampnames}
	
    if(showplot){
      if (savePDF) {pdf(paste0(outfile, ".pdf"), height=pdfheight, width=pdfwidth)}
      
        heatmap.2(seg.all.id.highr.sign.sort, trace="none",
                  Rowv=FALSE, Colv=FALSE, col=c("slateblue1","black","tomato"),
                  key=FALSE, cexCol=.5,
                  cexRow=.5, dendrogram='none')
               
        legend("top", ncol=3, c("up","same","down"), 
                         fill=rev(c("slateblue1","black","tomato")))
     if (savePDF) { dev.off()}
        
    }

	out <- bks.sign.sort.list
	
	return(out)
}
