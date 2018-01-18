#' @title Fit Segmented regression models on a feature/gene

#' @description fits segmented regression models
#' @inheritParams trendy
#' @return Trend: direction of each sample; -1: down, 0: no change, 1: up
#'  Slope: fitted slopes, Slope.Trend: sign of fitted slopes, 
#'  Slope.Pvalue: p value of each segment, Breakpoint: estimated breakpoints,
#'  Fitted.Values: fitted values AdjustedR2: adjusted r value of the model
#'  Fit: fit object
#' @author Ning Leng and Rhonda Bacher


fitSegRsq <- 
    function(Data, tVectIn = NULL, maxK = 5, minNumInSeg = 5, pvalCut = .1, 
        cutDiff = .1, numTry = 100, keepFit = FALSE) 
        
{
    whichFit <- seq_len(maxK)
    
    # Start with lm withOUT bp
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
    
    # Fit all possible breakpoints now:
    fit.bp.all <- lapply(whichFit, .breakpointFit, tVectIn, lmLinear, numTry)
    isna <- which(vapply(fit.bp.all, function(i) {
        "try-error" %in% class(i)}, logical(1)))
        
    # If it is not solved in 100 trys then return lm results (if numTry=100)
    if (length(isna) == maxK) {
        OUT <- list(Trends = lm.id.sign, Segment.Slopes = lm.slp, 
            Segment.Trends = lm.sign, 
            Segment.Pvalues = lm.pval, Breakpoints = NA, 
            Fitted.Values = lm.fit,
            AdjustedR2 = lm.radj, Fit = lmLinear)
        if (keepFit == FALSE) { OUT <- OUT[seq_len(7)] }
        return(OUT) 
        break
    }
    
    
    fit.keep <- fit.bp.all
    # If one of whichFit cant be fitted then remove it.
    if (length(isna) > 0) { 
        fit.keep <- fit.keep[-isna]
        whichFit <- whichFit[-isna]
    }
    
    slp.l <- lapply(fit.keep, segmented::slope)
    radj <- vapply(fit.keep, function(i) {summary(i)$adj.r.squared},numeric(1))
    brk.l <- lapply(fit.keep, function(i) {i$psi[,2]})
    id.l <- lapply(fit.keep, function(i) {i$id.group})
    allBIC <- vapply(fit.keep, function(i) {BIC(i)}, numeric(1))
    
    
    if (length(whichFit) > 1) {
        rsq.diff <- diff(lm.rsq) # get differences in R^2
        rsq.whichmax <- which(rsq.diff > cutDiff) 
        radj.max <- max(radj)
        rsq.max <- max(lm.rsq)
        if (length(rsq.whichmax) > 0) {
            # Start with the largest difference in rsq value
            r.choose <- max(rsq.whichmax) + 1
        } 
        # Compare here using adjusted Rsq. 
        if (length(rsq.whichmax) == 0 & radj[1] > lm.radj) {
            r.choose <- 1
        } # If none larger than diff use smallest bp =1
        if (length(rsq.whichmax) == 0 & radj[1] <= lm.radj) {
            # if none of them satisfy, then take the linear, return values
            OUT <- list(Trends = lm.id.sign, Segment.Slopes = lm.slp, 
                Segment.Trends = lm.sign, 
                Segment.Pvalues = lm.pval, Breakpoints = NA, 
                Fitted.Values = lm.fit,
                AdjustedR2 = lm.radj, Fit = lmLinear)
            if (keepFit == FALSE) {OUT <- OUT[seq_len(7)]}
            return(OUT)
            break
        }
        
        # Now make sure that optimal value satisfies having minimum 
        # Number of segments, if it does not then decrease breakpoints.
        while ((min(table(id.l[[r.choose]])) < minNumInSeg) & r.choose > 1) {
            r.choose <- r.choose - 1
        }
        
        if (r.choose == 1 & (min(table(id.l[[r.choose]])) < minNumInSeg)) {
            # if 1 bp gives small segment, then take the linear
            OUT <- list(Trends = lm.id.sign, Segment.Slopes = lm.slp, 
                Segment.Trends = lm.sign, 
                Segment.Pvalues = lm.pval, Breakpoints = NA, 
                Fitted.Values = lm.fit,
                AdjustedR2 = lm.radj, Fit = lmLinear)
            if (keepFit == FALSE) {OUT <- OUT[seq_len(7)]}
            return(OUT)
            break
        }
    }
    
    if (length(whichFit) == 1) {
        r.choose <- 1
        radj.max <- radj
        rsq.max <- lm.rsq
    }
        
    
    # Finally decide if best BP model is better than linear
    
    if (radj.max < lm.radj || (rsq.max - lm.rsq) < cutDiff || 
        min(table(id.l[[r.choose]])) < minNumInSeg) {
            OUT <- list(Trends = lm.id.sign, Segment.Slopes = lm.slp, 
                Segment.Trends = lm.sign, 
                Segment.Pvalues = lm.pval, Breakpoints = NA, 
                Fitted.Values = lm.fit,
                AdjustedR2 = lm.radj, Fit = lmLinear)
                if (keepFit == FALSE) {OUT <- OUT[seq_len(7)]}
                return(OUT) 
                # If lm is better 
                break
    }
    # Actual k (before remove the NA fitting)
    r.choose.ori <- whichFit[r.choose] 
    
    fit.choose <- fit.keep[[r.choose]]
    fv.choose <- fitted.values(fit.choose)
    names(fv.choose) <- paste0(names(tVectIn), ".Fitted")
    
    bp.choose <- brk.l[[r.choose]]
    if (length(bp.choose) >= 1) {
        names(bp.choose) <- paste0("Breakpoint", seq_along(bp.choose))
    }
    slp.choose <- slp.l[[r.choose]][[1]][,1]
    names(slp.choose)<-paste0("Segment",seq_len(length(slp.choose)), ".Slope")
    
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

