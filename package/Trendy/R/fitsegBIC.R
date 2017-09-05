#' @title Fit Segmented regression models on a feature/gene

#' @description fits segmented regression models
#' @inheritParams trendy
#' @return Trend: direction of each sample; -1: down, 0: no change, 1: up
#'		Slope: fitted slopes, Slope.Trend: sign of fitted slopes, 
#'		Slope.Pvalue: p value of each segment, Breakpoint: estimated breakpoints, 
#'		Fitted.Values: fitted values AdjustedR2: adjusted r value of the model
#' 		Fit: fit object
#' @examples d1 <- rbind(c(rep(1,50),1:50), c(100:1))
#' 		rownames(d1) <- c("g1","g2")
#' 		fitSegBIC(d1, "g1")
#' @author Rhonda Bacher and Ning Leng
#' @import stats
#' @import segmented


fitSegBIC <- function(Data, Max.K = 5, T.Vect = NULL, Min.Num.In.Seg = 5,
                      Pval.Cut = .1, Num.Try = 100, Keep.Fit = FALSE) {

	
  t.l <- T.Vect[length(T.Vect)]
		
	step.r <- 1:Max.K
	dat.tmp <- Data
	seed.use <- 1

	# start with lm withOUT bp
	lm1 <- lm(dat.tmp ~ T.Vect)
	lm.radj <- summary(lm1)$adj.r.squared
	lm.rsq <- summary(lm1)$r.squared
	lm.slp <- coef(lm1)[2]
	names(lm.slp) <- paste0("Segment", seq_len(length(lm.slp)), ".Slope")
	lm.fit <- fitted.values(lm1)
	
	lm.pval <- coef(summary(lm1))[2,4]
	lm.sign <- ifelse(lm.slp > 0, 1, -1)
	
	lm.sign[which(lm.pval > Pval.Cut)] <- 0
	names(lm.sign) <- paste0("Segment", seq_len(length(lm.sign)), ".Trend")
	
	lm.id.sign <- rep(lm.sign, length(T.Vect))
	names(lm.id.sign) <- paste0(names(T.Vect), ".Trend")
	names(T.Vect) <- paste0(names(T.Vect), ".Fitted")
  bic.lm <- BIC(lm1)
  
	# fit all possible breakpoints now:
	fit.l.0 <- sapply(1:length(step.r), function(j) {
		i <- step.r[j]
		lmseg.try <- suppressMessages(try(segmented(lm1, seg.Z = ~T.Vect, 
			            psi = round(seq(1, t.l, length.out = i + 2)[2:(i + 1)]), 
			            control = seg.control(seed = seed.use)), silent = T))
		seed.use2 <- seed.use
		while("try-error" %in% class(lmseg.try) & seed.use2 <= Num.Try) {
	 	   seed.use2 <- seed.use2 + 1
		   lmseg.try <- suppressMessages(try(segmented(lm1, seg.Z = ~T.Vect, 
			               psi = round(seq(1, t.l, length.out = i + 2)[2:(i + 1)]), 
			                control = seg.control(seed = seed.use2)), silent = T))}
	  OUT <- lmseg.try 
  }, simplify = F)

	isna <- which(sapply(fit.l.0, function(i) "try-error" %in% class(i)))

	if (length(isna) == Max.K) {
		OUT <- list(Trends = lm.id.sign, Segment.Slopes = lm.slp, Segment.Trends = lm.sign, 
					          Segment.Pvalues = lm.pval, Breakpoints = NA, Fitted.Values = lm.fit,
                    AdjustedR2 = lm.radj, Fit = lm1)
		if (Keep.Fit == FALSE) { OUT <- OUT[1:7] }
		return(OUT) 
	  # if it is not solved in 100 trys then return lm results (if Num.Try=100)
		break}


	fit.l <- fit.l.0
	if (length(isna) > 0) { # if one of step.r cant be fitted then remove it.
		fit.l <- fit.l.0[-isna]
		step.r <- step.r[-isna]
  }
  
  # Get info for each fit:
	slp.l <- sapply(fit.l, segmented::slope, simplify = F)
	radj <- sapply(fit.l, function(i) {summary(i)$adj.r.squared})
	brk.l <- sapply(fit.l ,function(i) {i$psi[,2]}, simplify = F)
	id.l <- sapply(fit.l, function(i) {i$id.group}, simplify = F)
  allBIC <- sapply(fit.l, function(i) {BIC(i)}, simplify = T)


	if (length(step.r) > 1) {
		bic.whichmin <- which.min(allBIC)
		if (length(bic.whichmin) > 0) {
			r.choose <- max(bic.whichmin)
    } 
		# compare here using bic. 
		if (allBIC[r.choose] >= bic.lm) {
	  	  	# if none of them satisfy, then take the linear, return values
	  		OUT <- list(Trends = lm.id.sign, Segment.Slopes = lm.slp, Segment.Trends = lm.sign, 
					          Segment.Pvalues = lm.pval, Breakpoints = NA, Fitted.Values = lm.fit,
                    AdjustedR2 = lm.radj, Fit = lm1)
	  		if (Keep.Fit == FALSE) {OUT <- OUT[1:7]}
	  		return(OUT)
	  		break
    }

		# now make sure that optimal value satisfies having minimum number of segments, 
		# if it does not then decrease breakpoints.
		while ((min(table(id.l[[r.choose]])) < Min.Num.In.Seg) & r.choose > 1) {
				r.choose <- r.choose - 1}

		if (r.choose == 1 & (min(table(id.l[[r.choose]])) < Min.Num.In.Seg)) {
	  	  	# if 1 bp gives small segment, then take the linear
	  	  	OUT <- list(Trends = lm.id.sign, Segment.Slopes = lm.slp, Segment.Trends = lm.sign, 
					          Segment.Pvalues = lm.pval, Breakpoints = NA, Fitted.Values = lm.fit,
                    AdjustedR2 = lm.radj, Fit = lm1)
	  		if(Keep.Fit == FALSE) {OUT <- OUT[1:7]}
	  		return(OUT)
	  	  	break
    }
	}

	if (length(step.r) == 1) {
		r.choose <- 1}

	# finally decide if best BP model is better than linear
	if(allBIC[r.choose] >= bic.lm){
		OUT <- list(Trends = lm.id.sign, Segment.Slopes = lm.slp, Segment.Trends = lm.sign, 
					          Segment.Pvalues = lm.pval, Breakpoints = NA, Fitted.Values = lm.fit,
                    AdjustedR2 = lm.radj, Fit = lm1)
		if(Keep.Fit == FALSE) {OUT <- OUT[1:7]}
	    return(OUT) 
		# if lm is better 
		break}

	r.choose.ori <- step.r[r.choose] # actual k (before remove the NA fitting)

  fit.choose <- fit.l[[r.choose]]
	fv.choose <- fitted.values(fit.choose)
	names(fv.choose) <- paste0(names(T.Vect), ".Fitted")
	
	bp.choose <- brk.l[[r.choose]]
  if(length(bp.choose) >= 1) {
    names(bp.choose) <- paste0("Breakpoint", 1:length(bp.choose))}
	slp.choose <- slp.l[[r.choose]][[1]][,1]
	names(slp.choose) <- paste0("Segment", seq_len(length(slp.choose)), ".Slope")
	
	slp.t <- slp.l[[r.choose]][[1]][,3]
	slp.pval <- pt(-abs(slp.t), 1)
	names(slp.pval) <- paste0("Segment", seq_len(length(slp.pval)), ".Pvalue")
	
	slp.sign <- ifelse(slp.t > 0, 1, -1)
	slp.sign[which(slp.pval > Pval.Cut)] <- 0
	names(slp.sign) <- paste0("Segment", seq_len(length(slp.sign)), ".Trend")
	
	id.choose <- id.l[[r.choose]]
	id.sign <- slp.sign[id.choose + 1]
	names(id.sign) <- paste0(names(T.Vect), ".Trend")
	
  OUT = list(Trends = id.sign, Segment.Slopes = slp.choose, Segment.Trends = slp.sign, 
				      Segment.Pvalues = slp.pval, Breakpoints = bp.choose, Fitted.Values = fv.choose,
				      AdjustedR2 = radj[r.choose], Fit = fit.choose)
	if(Keep.Fit == FALSE) {OUT <- OUT[1:7]}

	return(OUT)
}
