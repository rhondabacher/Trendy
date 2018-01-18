#' @title Extract pattern from segmented regression

#' @description find dynamic genes that follow a given pattern

#' @param trendyOutData output from trendy() function
#' @param Pattern vector containing pattern to search genes/features 
#'  (e.g, c("up", "down")). If length
#'	is one (e.g c("up")) then it will only consider features with 
#'  constant pattern across the entire time-course.
#' @param adjR2Cut only consider features with adjusted R^2 > adjR2Cut. 
#'	Default = .5.
#' @param Delay search for pattern starting after certain 
#'	time-point (e.g. only genes with a breakpoint > 10).
#' @return Genes: names of genes/features containing pattern 
#'	and the breakpoints corresponding to the pattern.
#' @examples 
#'	myTrends <- trendy(TrendyExampleData[seq_len(5),], tVect=seq_len(40))
#'  myTrends <- results(myTrends)
#'	#extractPattern(myTrends, Pattern = c("up")) #increasing only features
#'  #extractPattern(myTrends, Pattern = c("up", "down")) #features with a peak
#'  #extractPattern(myTrends, Pattern = c("up", "down"), Delay = 20)
#' @author Rhonda Bacher
#' @export


extractPattern <- 
    function(trendyOutData, Pattern = NULL, adjR2Cut = .5, Delay = 0)  
{
    if (is.null(Pattern)) {stop("Must specify a pattern")}
    
    ogpat <- Pattern
    # Restrict to genes that have certain breakpoint pattern
  
    # Get R2_adj for each gene
    segdata.radj <- sort(vapply(trendyOutData, function(i) {
                                    i$AdjustedR2
                                }, numeric(1)), decreasing=TRUE)
    # Get genes with fit >= radjcut
    genes.pass <- names(segdata.radj)[which(segdata.radj >= adjR2Cut)] 
    if (length(genes.pass) == 0) {
        stop("No gene pass the adjusted R^2 cutoff filter!")
    }

    segdata.pass <- trendyOutData[genes.pass]
  
    # Get breakpoints for all genes
    segdata.bks <- lapply(segdata.pass, function(i) i$Breakpoints) 
    segdata.bks[is.na(segdata.bks)] <- 0
    # Get slps for all genes
    segdata.slps <- lapply(segdata.pass, function(i) i$Segment.Trends) 
    Pattern[Pattern == "up"] <- 1
    Pattern[Pattern == "same"] <- 0
    Pattern[Pattern == "down"] <- 2
  
  
    if (length(Pattern) == 1) {
        Pattern = rep(Pattern, length(segdata.pass[[1]]$Trends))
    }
  
    Pattern <- paste(Pattern, collapse="")
  
  
    # Find genes that have the pattern of interest.
    # Hard to initialize since not sure how many will match.
    genes <- list()
    k <- 1
    trackNames <- c()
    for (j in seq_along(genes.pass)) {
        gslp <- segdata.slps[[genes.pass[j]]]
        gslp[gslp == -1] <- 2
    
        if (length(gslp) == 1) { gslp = rep(gslp, nchar(Pattern))}
        slps <- paste(gslp, collapse="")
    
        if (grepl(Pattern, slps) == TRUE) {
            patstart <- gregexpr(Pattern, slps)[[1]]
      
            for (i in seq_along(patstart)) {
      
                if (segdata.bks[[genes.pass[j]]][i] >= Delay) {
                    patend <- patstart[i] + nchar(Pattern)-2
                    brk <- segdata.bks[[genes.pass[j]]][patstart[i]:patend]
                    names(brk) <- paste0(rep("BreakPoint", nchar(Pattern) - 1), 
                                    1:(nchar(Pattern) - 1))
                    genes[[k]] <- brk
                    k = k + 1
                    trackNames <- c(trackNames, genes.pass[j])
                }
            }
        }
    }

    if (length(genes) == 0) {stop("No genes follow this pattern!")}
    genes <- data.frame(Gene = trackNames, 
                        do.call(rbind, genes), stringsAsFactors=FALSE)

    genes <- genes[order(genes$BreakPoint1),]

    return(genes)

}