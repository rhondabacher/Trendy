#' @title obtain top genes from trendy results

#' @description reformats the list output for genes with a 
#'  given adjusted R^2 cutoff
#' @param trendyOutData output from the trendy function
#' @param adjR2Cut cutoff for the adjusted R^2. Genes whose adjusted
#'  R^2 is greater than adjR2Cut are called as significant. 

#' @return only significant genes will be included in the output. 
#'  The output is reformatted as:
#'      Trend direction of each sample; -1: down, 0: no change, 1: up
#'		Slope: fitted slopes, Slope.Trend: sign of fitted slopes, 
#'		Slope.Pvalue: p value of each segment, 
#'		Breakpoint: estimated breakpoints, Fitted.Values: fitted values 
#'		AdjustedR2: adjusted r value of the model
#'		Fit: fit object
#' @examples d1 <- rbind(c(rep(1,50),seq_len(50)), rev(seq_len(100)))
#'  rownames(d1) <- c("g1","g2")
#'  colnames(d1) <- paste0("time", seq_len(100))
#'  seg.all <- trendy(d1)
#'  seg.all <- results(seg.all)
#'  top.genes <- topTrendy(seg.all)
#' @export

topTrendy <- 
    function(trendyOutData, adjR2Cut = .5) 

{
    seg.all.radj <- sapply(trendyOutData, function(i) i$AdjustedR2)
    names(seg.all.radj) <- names(trendyOutData)
    seg.all.radj.sort <- sort(seg.all.radj, decreasing = TRUE)
    toKeep <- which(seg.all.radj.sort > adjR2Cut)
    seg.all.highr <- names(seg.all.radj.sort)[toKeep]

    if (length(seg.all.highr) == 0) {
        stop("No gene pass the adjusted R^2 cutoff filter!")
    }

    seg.use <- trendyOutData[seg.all.highr]	

    reformat <- lapply(seq_len(7), function(i) { 
        lapply(seg.use, function(j) j[[i]])
    })
    names(reformat) <- names(trendyOutData[[1]])

    reformat[[1]] <- do.call(rbind, reformat[[1]])
    reformat[[2]] <- formatFunc(reformat[[2]])
    reformat[[3]] <- formatFunc(reformat[[3]])
    reformat[[4]] <- formatFunc(reformat[[4]])
    reformat[[5]] <- formatFunc(reformat[[5]])
    reformat[[6]] <- do.call(rbind, reformat[[6]])
    reformat[[7]] <- do.call(c, reformat[[7]])

    return(reformat)
}







