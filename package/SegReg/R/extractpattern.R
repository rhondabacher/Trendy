#' @title extract pattern from segmented regression
#' @usage extractpattern(segdata, radjcut = 10, pattern = 0)
#' @param segdata segreg object
#' @param radjcut only consider genes with adjusted R^2 > radjcut. Default = .5.
#' @param pattern vector containing pattern to search genes. If length is one then consider genes with constant pattern.
#' @return genes: name of genes containing pattern
#' @examples extractpattern(segdata, pattern = c("up")) #increasing only genes
#' extractpattern(segdata, pattern = c("up", "down")) #genes with a peak
#' @author Rhonda Bacher
###################
# Enter pattern of interest and return genes having that pattern.
###################
extractpattern <- function (segdata, radjcut = .5, pattern = NULL) 
{
	
	if(is.null(pattern)) stop("Must specify a pattern")
		
	##restrict to genes that have certain breakpoint pattern
	
	segdata.radj <- sapply(segdata, function(i)i$radj) #get radj for each gene
	
	genes.pass <- names(segdata.radj)[which(segdata.radj >= radjcut)] #get genes with fit >= radjcut
	segdata.pass <- segdata[genes.pass]


	segdata.bks <- sapply(segdata.pass, function(i) i$bp) #Get breakpoints for all genes
	segdata.slps <- sapply(segdata.pass, function(i) i$slp.sign) #Get slps for all genes
	
	pattern[pattern == "up"] <- 1
	pattern[pattern == "same"] <- 0
	pattern[pattern == "down"] <- -1
	pattern <- paste(pattern, collapse="")
	
	
	#find genes that have peaks at ANY breakpoint
	genes <- c()
	for (j in 1:length(genes.pass)) {
		
			slps <- paste(segdata.slps[[genes.pass[j]]], collapse="")
			if(grepl(pattern, slps) == TRUE) {
				
				brk <- genes.pass[j]
				genes <- c(genes, brk)
			}

	}
	

	return(genes)
    
}
