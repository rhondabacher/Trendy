#' @title Draw heatmap of gene expression trends

#' @description heatmap of the fitted trends
#' @param topTrendyData results from topTrendyData() function
#' @return The function takes significant genes/features called from 
#'  the topTrendyData() function. These genes are further grouped into three
#'  groups: up, down, or no change in the first segment. Within each group,
#'  the genes are sorted by their first break point. The heatmap shows 
#'  expression trends of these three groups of genes. In the heatmap, 
#'  red/blue/black represents up/down/nochange. A list of genes in the heatmap
#'  order is returned.
#' @author Ning Leng and Rhonda Bacher
#' @import graphics
#' @import grDevices
#' @importFrom gplots heatmap.2
#' @export
#' @examples m1 <- rbind(c(rep(1,50),seq_len(50)), rev(seq_len(100)))
#'  rownames(m1) <- c("g1","g2")
#'  colnames(m1) <- paste0("time", seq_len(100))
#'  myTrends <- results(trendy(m1))
#'  topGenes <- topTrendy(myTrends)
#'  #makeHeat <- trendHeatmap(topGenes)

trendHeatmap <- 
    function (topTrendyData)
    
{
    bks.all <- topTrendyData$Breakpoints
    if (sum(!is.na(bks.all)) == 0 ){
        stop("No breakpoints were found for these genes!")
    }
    
    bks.first <- bks.all[,1]
    names(bks.first) <- rownames(topTrendyData[[1]])
    bks.first.sort <- sort(bks.first)
    
    seg.all.id <- topTrendyData[[1]]
    seg.first <- seg.all.id[,1]
    
    genes <- rownames(seg.all.id)
    bks.sign.sort.list <- list(
        firstup = sort(bks.first[genes[which(seg.first == 1)]]),
        firstdown = sort(bks.first[genes[which(seg.first == -1)]]),
        firstnochange = sort(bks.first[genes[which(seg.first == 0)]]))
        
    bks.sign.sort <- c(bks.sign.sort.list[[1]], 
    bks.sign.sort.list[[2]],bks.sign.sort.list[[3]])
    
    seg.all.id.highr.sign.sort <- seg.all.id[names(bks.sign.sort),]
    
    colnames(seg.all.id.highr.sign.sort) <- gsub(".Fitted.Trend", "", 
    colnames(seg.all.id.highr.sign.sort), fixed=TRUE)
    
    gplots::heatmap.2(seg.all.id.highr.sign.sort, trace = "none",
        Rowv = FALSE, Colv = FALSE, col = c("slateblue1", "black", "tomato"),
        key = FALSE, cexRow = .5, cexCol=.5)
    legend("top", ncol = 3, c("up", "no change", "down"), 
        fill = rev(c("slateblue1", "black", "tomato")))
        
    return(bks.sign.sort.list)
}
