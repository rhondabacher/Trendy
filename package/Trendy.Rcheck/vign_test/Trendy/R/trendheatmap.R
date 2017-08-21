#' @title Draw heatmap of gene expression trends

#' @description heatmap of the fitted trends
#' @param Top.Trendy results from toptrendy() function
#' @return The function takes significant genes/features called from 
#'		the toptrendy() function. These genes are further grouped into three
#'		groups : up, down, or no change in the first segment. Within each group,
#'		the genes are sorted by their first break point. The heatmap shows expression 
#'		trends of these three groups of genes. In the heatmap, red/blue/black 
#'		represents up/down/nochange. A list of genes in the heatmap order is returned.
#' @author Ning Leng and Rhonda Bacher
#' @import graphics
#' @import grDevices
#' @importFrom gplots heatmap.2
#' @export

trendheatmap <- function(Top.Trendy){

	bks.all <- Top.Trendy$bp
	bks.first <- sapply(bks.all, function(i) i[1])
	names(bks.first) <- colnames(Top.Trendy[[1]])
	bks.first.sort <- sort(bks.first)

	seg.all.id.highr <- t(Top.Trendy[[1]])

	bks.sign.sort.list <- list(
		firstup = sort(bks.first[rownames(seg.all.id.highr)[which(seg.all.id.highr[,1] == 1)]]),
		firstdown = sort(bks.first[rownames(seg.all.id.highr)[which(seg.all.id.highr[,1] == -1)]]),
		firstnochange = sort(bks.first[rownames(seg.all.id.highr)[which(seg.all.id.highr[,1] == 0)]])
	)

	bks.sign.sort <- c(bks.sign.sort.list[[1]],bks.sign.sort.list[[2]],bks.sign.sort.list[[3]])

	seg.all.id.highr.sign.sort <- seg.all.id.highr[names(bks.sign.sort),]

	gplots::heatmap.2(seg.all.id.highr.sign.sort, trace = "none",
						Rowv = FALSE, Colv = FALSE, col = c("slateblue1", "black", "tomato"),
						key = FALSE, cexRow = .5)
	legend("top", ncol = 3, c("up", "no change", "down"), fill = rev(c("slateblue1", "black", "tomato")))

return(bks.sign.sort.list)
}
