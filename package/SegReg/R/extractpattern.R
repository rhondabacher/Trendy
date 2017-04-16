#' @title Extract genes following a given pattern from segmented regression
#' @usage extractpattern(segdata, radjcut = 10, pattern = NULL, delay=0)
#' @param segdata segreg object
#' @param radjcut only consider features with adjusted r_sq > radjcut. Default = .5.
#' @param pattern vector containing pattern to search features. If length is one then only consider features with constant pattern across the entire time-course. 
#' @param delay search for pattern starting after certain time-point.
#' @return Gene: names of features containing pattern.
#' @return Breakpoint: estimated breakpoint connecting the segments of the pattern.
#' @examples extractpattern(segdata, pattern = c("up")) #increasing only features
#' extractpattern(segdata, pattern = c("up", "down")) #features with a peak
#' extractpattern(segdata, pattern = c("up", "down"), delay=20) #features with a peak after 20th time-point
#' extractpattern(segdata, pattern = c("same", "down")) #features which display delaying decreasing behavior.
#' @author Rhonda Bacher
###################
# Enter pattern of interest and return genes having a given pattern.
###################
#' @export

extractpattern <- function (segdata, radjcut = .5, pattern = NULL,  delay = 0) 
{
  
  if (is.null(pattern)) stop("Must specify a pattern")

  ogpat <- pattern

  segdata.radj <- sort(sapply(segdata, function(i)i$radj), decreasing=TRUE) #get radj for each gene
  
  genes.pass <- names(segdata.radj)[which(segdata.radj >= radjcut)] #get genes with fit >= radjcut
  segdata.pass <- segdata[genes.pass]
  
  
  segdata.bks <- lapply(segdata.pass, function(i) i$bp) #Get breakpoints for all genes
  segdata.bks[is.na(segdata.bks)] <- 0
  segdata.slps <- lapply(segdata.pass, function(i) i$slp.sign) #Get slps for all genes
  
  pattern[pattern == "up"] <- 1
  pattern[pattern == "same"] <- 0
  pattern[pattern == "down"] <- 2
  pattern[pattern == -1] <- 2
  
  
  if(length(pattern) == 1) {
    pattern = rep(pattern, length(segdata.pass[[1]]$id.sign))
  }
  
  pattern <- paste(pattern, collapse="")
  
  
  #find genes that have peaks at ANY breakpoint
  genes <- list()
  k <- 1
  trackNames <- c()
  for (j in 1:length(genes.pass)) {
    
    gslp <- segdata.slps[[genes.pass[j]]]
    gslp[gslp == -1] <- 2
    
    if(length(gslp) == 1) { gslp = rep(gslp, length(segdata.pass[[1]]$id.sign))}
    slps <- paste(gslp, collapse="")
    
    
    if(grepl(pattern, slps) == TRUE) {
      
      patstart <- gregexpr(pattern, slps)[[1]]
      
      for(i in 1:length(patstart)) {
      
          if(segdata.bks[[genes.pass[j]]][i] >= delay) {
            patend <- patstart[i] + nchar(pattern)-2
            brk <- segdata.bks[[genes.pass[j]]][patstart[i]:patend]
            names(brk) <- paste0(rep("BreakPoint", nchar(pattern) - 1), 1:(nchar(pattern) - 1))
            genes[[k]] <- brk
            k = k + 1
            trackNames<- c(trackNames, genes.pass[j])
        }
      
      }
    }
  }

  
  genes <- data.frame(Gene = trackNames, do.call(rbind, genes), stringsAsFactors=FALSE)

  genes <- genes[order(genes$BreakPoint1),]

  return(genes)
  
}