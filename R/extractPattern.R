#' @title Extract pattern from segmented regression

#' @description find dynamic genes that follow a given pattern

#' @param trendyOutData output from trendy() function
#' @param Pattern vector containing pattern to search genes/features 
#'  (e.g, c("up", "down")), no-change is designated by "same". If length
#'  is one (e.g c("up")) then it will only consider features with 
#'  constant pattern across the entire time-course.
#' @param adjR2Cut only consider features with adjusted R^2 > adjR2Cut. 
#'  Default = .5.
#' @param Delay search for pattern starting after certain 
#'  time-point (e.g. only genes with a breakpoint > 10).
#' @return Genes: names of genes/features containing pattern 
#'  and the breakpoints corresponding to the pattern.
#' @examples 
#'  myTrends <- trendy(trendyExampleData[seq_len(5),], tVect=seq_len(40))
#'  myTrends <- results(myTrends)
#'  #extractPattern(myTrends, Pattern = c("up")) #increasing only features
#'  #extractPattern(myTrends, Pattern = c("same", "down")) 
#'  #extractPattern(myTrends, Pattern = c("up", "down"), Delay = 20)
#' @author Rhonda Bacher
#' @export
extractPattern<-
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
    Pattern[Pattern == "same" | Pattern == "no change"] <- 0
    Pattern[Pattern == "down"] <- 2

    if (length(Pattern) == 1) {
        Pattern = rep(Pattern, length(segdata.pass[[1]]$Trends))
    }

    Pattern <- paste(Pattern, collapse="")



    gslps <- lapply(segdata.slps, function(x) {
        x[x == -1] <- 2
        if (length(x) == 1) {x = rep(x, nchar(Pattern))}
        x <- paste(x, collapse="")
        return(x)
    })
    whichg <- which(grepl(Pattern, gslps))
    if (length(whichg) == 0) {stop("No genes follow this pattern!")}

    gslps <- gslps[whichg]

    patstarts <- lapply(gslps, function(x) {
      patstart <- gregexpr(Pattern, x)[[1]]
      return(patstart)
    })

    outPats <- lapply(seq_along(whichg), function(x) lapply(patstarts[[x]], 
        function(y) {
            patend <- y + nchar(Pattern)-2
            brk <- segdata.bks[[whichg[x]]][seq(y,patend,1)]
            names(brk) <- paste0(rep("BreakPoint", nchar(Pattern) - 1), 
                                    seq(1,(nchar(Pattern) - 1),1))
            RTN <- c(names(gslps[x]), brk)
            return(RTN)
        }))

    outPats <- data.frame(do.call(rbind, do.call(c, outPats)), 
        stringsAsFactors = FALSE)
    outPats[,-1] <- apply(outPats[,-1, drop = FALSE], 2, as.numeric)
    colnames(outPats)[1] <- c("Gene")
    outPats <- outPats[order(outPats$BreakPoint1),]
	
	outPats <- outPats[which(outPats$BreakPoint1 > Delay),]

    return(outPats)
}