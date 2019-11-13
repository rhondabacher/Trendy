#' @title Fit segmented regression models on a feature/gene

#' @description fits segmented regression models
#' @inheritParams trendy
#' @return Trend: direction of each sample; -1: down, 0: no change, 1: up
#'  Slope: fitted slopes, Slope.Trend: sign of fitted slopes, 
#'  Slope.Pvalue: p value of each segment, Breakpoint: estimated breakpoints,
#'  Fitted.Values: fitted values AdjustedR2: adjusted r value of the model
#'  Fit: fit object
#' @author Rhonda Bacher and Ning Leng
#' @export

fitSegBIC <- function(Data, maxK = 2, tVectIn = NULL, 
                      minNumInSeg = 5, pvalCut = .1, 
                      numTry = 5, keepFit = FALSE) 

{
    
    whichFit <- seq_len(maxK)
    
		# If any replicates, jitter a small amount to help with the segmented fitting:
		if (length(unique(tVectIn)) < length(tVectIn)) {tVectIn <- jitter(tVectIn, .1)}
			
    # Start with lm without any breaks
    lmLinear <- lm(Data ~ tVectIn)
    lm.radj <- summary(lmLinear)$adj.r.squared
    lm.rsq <- summary(lmLinear)$r.squared
    lm.slp <- coef(lmLinear)[2]
    names(lm.slp) <- paste0("Segment", seq_len(length(lm.slp)), ".Slope")
    lm.fit <- fitted.values(lmLinear)
    lm.pval <- coef(summary(lmLinear))[2,4]
    lm.sign <- ifelse(lm.slp > 0, 1, -1)
    
    lm.sign[which(lm.pval > pvalCut)] <- 0
    names(lm.sign) <- paste0("Segment", seq_len(length(lm.sign)), ".Trend")
    
    lm.id.sign <- rep(lm.sign, length(tVectIn))
    names(lm.id.sign) <- paste0(names(tVectIn), ".Trend")
    names(lm.fit) <- paste0(names(tVectIn), ".Fitted")
    bic.lm <- BIC(lmLinear)
    
    # Fit all possible breakpoints now:
    fit.bp.all <- lapply(whichFit, breakpointFit, tVectIn, lmLinear, numTry)

    isna <- which(vapply(fit.bp.all, function(i) {
                (class(i)[1] == "character")
              }, logical(1)))

    if (length(isna) == maxK) { #No bp models valid, return linear results
        OUT <- list(Trends = lm.id.sign, Segment.Slopes = lm.slp,
            Segment.Trends = lm.sign,
            Segment.Pvalues = lm.pval, Breakpoints = NA,
            Fitted.Values = lm.fit,
            AdjustedR2 = lm.radj, Fit = lmLinear)
            if (keepFit == FALSE) { OUT <- OUT[seq_len(7)] }
                return(OUT)
            # If it is not solved in 100 trys then return lm
            # results (if numTry=100)
    }


    fit.keep <- fit.bp.all
    if (length(isna) > 0) { # if one of whichFit cant be fitted then remove it.
        fit.keep <- fit.keep[-isna]
        whichFit <- whichFit[-isna]
    }

    # Get info for each fit:
    slp.l <- lapply(fit.keep, function(i) {segmented::slope(i)})
    radj <- vapply(fit.keep, function(i) {summary(i)$adj.r.squared},numeric(1))
    brk.l <- lapply(fit.keep, function(i) {i$psi[,2]})
    id.l <- lapply(fit.keep, function(i) {i$id.group})
    allBIC <- vapply(fit.keep, function(i) {BIC(i)}, numeric(1))


    if (length(whichFit) >= 1) {
        bic.whichmin <- which.min(allBIC)
        if (length(bic.whichmin) > 0) {
            r.choose <- max(bic.whichmin)
        }
        # compare here using bic.
        if (allBIC[r.choose] >= bic.lm) {
            # if linear has smaller BIC then take the linear, return values
            OUT <- list(Trends = lm.id.sign, Segment.Slopes = lm.slp,
                Segment.Trends = lm.sign,
                Segment.Pvalues = lm.pval, Breakpoints = NA,
                Fitted.Values = lm.fit,
                AdjustedR2 = lm.radj, Fit = lmLinear)
                if (keepFit == FALSE) {OUT <- OUT[seq_len(7)]}
                    return(OUT)
            }

            # now make sure that best fit model satisfies having minimum
            # number of segments, if it does not then decrease breakpoints.
            while((min(table(id.l[[r.choose]])) < minNumInSeg | any(brk.l[[r.choose]] < tVectIn[minNumInSeg])) & r.choose > 1) {
                r.choose <- r.choose - 1
            }
						checkMinSeg = 0
						while (any(checkMinSeg < minNumInSeg)) {
							checkMinSeg <- c()
              breaks.temp <- brk.l[[r.choose]]
              if (length(breaks.temp) == 1) {
                checkMinSeg <- c(length(tVectIn[tVectIn < breaks.temp]), length(tVectIn[tVectIn > breaks.temp]))
              }
              if (length(breaks.temp) >= 2) {
                breaks.temp <- c(tVectIn[1], breaks.temp, tVectIn[length(tVectIn)])
                for (i in 1:(length(breaks.temp)-1)) {
							    checkMinSeg <- c(checkMinSeg, length(tVectIn[tVectIn < breaks.temp[i+1] & tVectIn > breaks.temp[i]]))
							  }
              }
							if (any(checkMinSeg < minNumInSeg)) {r.choose <- r.choose - 1}
							if (r.choose == 0) {break}
						}
		        if (r.choose == 0) {
		            # take the linear, return values
		            OUT <- list(Trends = lm.id.sign, Segment.Slopes = lm.slp,
		                Segment.Trends = lm.sign,
		                Segment.Pvalues = lm.pval, Breakpoints = NA,
		                Fitted.Values = lm.fit,
		                AdjustedR2 = lm.radj, Fit = lmLinear)
		                if (keepFit == FALSE) {OUT <- OUT[seq_len(7)]}
		                    return(OUT)
		            }
            if (r.choose == 1 & (min(table(id.l[[r.choose]])) < minNumInSeg | any(brk.l[[r.choose]] < tVectIn[minNumInSeg]))){
                # if 1 bp gives too small segment, then take the linear
                OUT <- list(Trends = lm.id.sign, Segment.Slopes = lm.slp,
                    Segment.Trends = lm.sign,
                    Segment.Pvalues = lm.pval, Breakpoints = NA,
                    Fitted.Values = lm.fit,
                    AdjustedR2 = lm.radj, Fit = lmLinear)
                    if(keepFit == FALSE) {OUT <- OUT[seq_len(7)]}
                        return(OUT)
            }
    }

    # Finally decide if best BP model is better than linear,
    # If not return linear
    if (allBIC[r.choose] >= bic.lm) {
        OUT <- list(Trends = lm.id.sign, Segment.Slopes = lm.slp,
            Segment.Trends = lm.sign,
            Segment.Pvalues = lm.pval, Breakpoints = NA,
            Fitted.Values = lm.fit,
            AdjustedR2 = lm.radj, Fit = lmLinear)
        if (keepFit == FALSE) {OUT <- OUT[seq_len(7)]}
        return(OUT)
        # If lm is better
    }

    # # Actual k (before remove the NA fitting)
    r.choose.ori <- whichFit[r.choose]

    fit.choose <- fit.keep[[r.choose]]
    fv.choose <- fitted.values(fit.choose)
    names(fv.choose) <- paste0(names(tVectIn), ".Fitted")

    bp.choose <- brk.l[[r.choose]]
    if (length(bp.choose) >= 1) {
        names(bp.choose) <- paste0("Breakpoint", seq_along(bp.choose))}
    slp.choose <- slp.l[[r.choose]][[1]][,1]
    names(slp.choose) <- paste0("Segment",seq_len(length(slp.choose)),".Slope")

    slp.t <- slp.l[[r.choose]][[1]][,3]
    slp.pval <- pt(-abs(slp.t), 1)
    names(slp.pval) <- paste0("Segment", seq_len(length(slp.pval)), ".Pvalue")

    slp.sign <- ifelse(slp.t > 0, 1, -1)
    slp.sign[which(slp.pval > pvalCut)] <- 0
    names(slp.sign) <- paste0("Segment", seq_len(length(slp.sign)), ".Trend")

    id.choose <- id.l[[r.choose]]
    id.sign <- slp.sign[id.choose + 1]
    names(id.sign) <- paste0(names(tVectIn), ".Trend")

    OUT = list(Trends = id.sign, Segment.Slopes = slp.choose,
        Segment.Trends = slp.sign,
        Segment.Pvalues = slp.pval, Breakpoints = bp.choose,
        Fitted.Values = fv.choose,
        AdjustedR2 = radj[r.choose], Fit = fit.choose)
    if (keepFit == FALSE) {OUT <- OUT[seq_len(7)]}

    return(OUT)
}
