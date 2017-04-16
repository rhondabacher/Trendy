#' @title Distribution of breakpoints
#' @usage bpdist(topres, ndigits=0, makeplot=FALSE, savePDF=FALSE, barcolor=NULL, outfile, xlabname=NULL)
#' @param topres results from topsegreg() function
#' @param ndigits how many digits to be used when rounding. Default is 0 (return integers).
#' @param makeplot whether to plot the breakpoint distribution. Default is true.
#' @param savePDF whether to save the plot. Default is FALSE.
#' @param outfile if savePDF is TRUE then this is the name of the file and path of the plot. Default is breakpoint_distribution.pdf saved in the working directory.
#' @param xlabname name of x-axis. Default is Time.
#' @param barcolor color of the barplot.
#' @return The function takes significant genes called from 
#' the topsegreg() function. For any time point, this function calculates how many genes
#' have breakpoint at this time point. The output is the numbers of genes sorted
#' by time point, and optionally a barplot of the breakpoint distribution.
#' @author Rhonda Bacher & Ning Leng
#' @export

bpdist <- function(topres, ndigits = 0, makeplot=TRUE, savePDF=FALSE, outfile=NULL, xlabname=NULL, barcolor=NULL){
	  
	  tmpbp <- topres$bp
	  tab <- table(round(unlist(tmpbp), ndigits))
	  tabout <- tab[order(as.numeric(names(tab)))]
	  gap <- 1/(10^ndigits)
	  tabname.num <- as.numeric(names(tabout))
	  namesuse <- as.character(seq(min(tabname.num),max(tabname.num),gap))
	  tabnames <- names(tabout)
	  outv <- rep(0,length(namesuse))
	  names(outv) <- namesuse
	  outv[tabnames] <- tabout
  
	  if(makeplot)
	  {
	    if (is.null(xlabname)) {xlabname="Time"}
	    if (is.null(barcolor)) {barcolor="royalblue1"} 
	    if (is.null(outfile)) {outfile="breakpoint_distribution"} 
	    if(savePDF == TRUE) {
	      pdf(paste0(outfile,".pdf"),height=5, width=10)
	      par(mar=c(5,5,1,2), cex.lab=2, cex.axis=1.5)
	      barplot(outv, ylab="Number of Breakpoints",
	              col=barcolor, xlab=xlabname,
	              ylim=c(0,max(tab)+max(tab)*.1))
	      box(lwd=2)
	      dev.off()
	    }
    
	    par(mar=c(5,5,1,2), cex.lab=2, cex.axis=1.5)
	    barplot(outv, ylab="Number of Breakpoints",
	            col=barcolor, xlab=xlabname,
	            ylim=c(0,max(tab)+max(tab)*.1))
	    box(lwd=2)
    
    
	  }
  
	  return(outv)
	}
	