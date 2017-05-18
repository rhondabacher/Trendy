#' @title Extract pattern from segmented regression
#' @usage extractpattern(Seg.Data, Pattern = NULL, Radj.Cut = .5, Delay = 0)

#' @description find dynamic genes that follow a given pattern

#' @param Seg.Data output from trendy() function
#' @param Pattern vector containing pattern to search genes/features. If length
#'		is one then it will only consider features with constant pattern 
#'		across the entire time-course.
#' @param Radj.Cut only consider features with adjusted R^2 > Radj.Cut. 
#'		Default = .5.
#' @param Delay search for pattern starting after certain time-point (e.g. only
#'		 genes with a breakpoint > 10).
#' @return genes: names of genes/features containing pattern and the breakpoints
#'		corresponding to the pattern.
#' @examples 
#'	myTrends <- trendy(TrendyExData)
#'	extractpattern(myTrends, Pattern = c("up")) #increasing only features
#'  extractpattern(myTrends, Pattern = c("up", "down")) #features with a peak
#'  extractpattern(myTrends, Pattern = c("up", "down"), Delay = 20)
#' @author Rhonda Bacher
#' @export


extractpattern <- function(Seg.Data, Pattern = NULL, Radj.Cut = .5, Delay = 0)  {
  
  if(is.null(Pattern)) stop("Must specify a pattern")
  ogpat <- Pattern
  ##restrict to genes that have certain breakpoint pattern
  
  segdata.radj <- sort(sapply(Seg.Data, function(i)i$radj), decreasing=TRUE) #get radj for each gene
  
  genes.pass <- names(segdata.radj)[which(segdata.radj >= Radj.Cut)] #get genes with fit >= radjcut
  segdata.pass <- Seg.Data[genes.pass]
  
  
  segdata.bks <- lapply(segdata.pass, function(i) i$bp) #Get breakpoints for all genes
  segdata.bks[is.na(segdata.bks)] <- 0
  segdata.slps <- lapply(segdata.pass, function(i) i$slp.sign) #Get slps for all genes
  
  Pattern[Pattern == "up"] <- 1
  Pattern[Pattern == "same"] <- 0
  Pattern[Pattern == "down"] <- 2
  
  
  if(length(Pattern) == 1) {
    Pattern = rep(Pattern, length(segdata.pass[[1]]$id.sign))
  }
  
  Pattern <- paste(Pattern, collapse="")
  
  
  #find genes that have peaks at ANY breakpoint
  genes <- list()
  k <- 1
  trackNames <- c()
  for (j in 1:length(genes.pass)) {
    
    gslp <- segdata.slps[[genes.pass[j]]]
    gslp[gslp == -1] <- 2
    
    if(length(gslp) == 1) { gslp = rep(gslp, length(segdata.pass[[1]]$id.sign))}
    slps <- paste(gslp, collapse="")
    
    
    if(grepl(Pattern, slps) == TRUE) {
      
      patstart <- gregexpr(Pattern, slps)[[1]]
      
      for(i in 1:length(patstart)) {
      
          if(segdata.bks[[genes.pass[j]]][i] >= Delay) {
            patend <- patstart[i] + nchar(Pattern)-2
            brk <- segdata.bks[[genes.pass[j]]][patstart[i]:patend]
            names(brk) <- paste0(rep("BreakPoint", nchar(Pattern) - 1), 1:(nchar(Pattern) - 1))
            genes[[k]] <- brk
            k = k + 1
            trackNames<- c(trackNames, genes.pass[j])
        }
      
      }
    }
  }

  if(length(genes) == 0) {stop("No genes follow this pattern!")}
	  
  genes <- data.frame(Gene = trackNames, do.call(rbind, genes), stringsAsFactors=FALSE)
  
  genes <- genes[order(genes$BreakPoint1),]

  return(genes)
  
}