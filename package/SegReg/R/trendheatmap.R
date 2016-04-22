#' @title Draw heatmap of gene expression trends
#' @usage trendheatmap(topres, showplot)
#' @param topres results from topsegreg() function
#' @param showplot whether show the heatmap
#' @return The function takes significant genes called from 
#' the topsegreg() function. These genes are further grouped into three
#' groups : go up/down or no change in the first segment. Within each group,
#' the genes are sorted by their first break point. The heatmap shows expression trends of these 
#' three groups of genes. In the heatmap, red/blue/black represents up/down/nochange.
#' @author Ning Leng


trendheatmap <- function(topres, showplot=TRUE){

bks.all <- topres$bp
bks.first <- sapply(bks.all, function(i)i[1])
names(bks.first) <- colnames(topres[[1]])
bks.first.sort <- sort(bks.first)

seg.all.id.highr <- t(topres[[1]])

bks.sign.sort.list <- list(
									 firstup=sort(bks.first[rownames(seg.all.id.highr)[which(seg.all.id.highr[,1]==1)]]),
									 firstdown=sort(bks.first[rownames(seg.all.id.highr)[which(seg.all.id.highr[,1]==-1)]]),
									 firstnochange=sort(bks.first[rownames(seg.all.id.highr)[which(seg.all.id.highr[,1]==0)]])
									 )

bks.sign.sort <- c(bks.sign.sort.list[[1]],bks.sign.sort.list[[2]],bks.sign.sort.list[[3]])

seg.all.id.highr.sign.sort <- seg.all.id.highr[names(bks.sign.sort),]

if(showplot){heatmap.2(seg.all.id.highr.sign.sort,trace="none",
					Rowv=FALSE,Colv=FALSE,col=c("slateblue1","black","tomato"),
					key=FALSE,
					cexRow=.5)
legend("top",ncol=3, c("up","no change","down"), fill=rev(c("slateblue1","black","tomato")))
									 }
out <- bks.sign.sort.list
}
