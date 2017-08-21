#' @title Fit Segmented regression models on a feature/gene

#' @description fits segmented regression models
#' @inheritParams trendy
#' @return id.sign: direction of each sample; -1: down, 0: no change, 1: up
#'		slp: fitted slopes, slp.sign: sign of fitted slopes, 
#'		slp.pval: p value of each segment, bp: estimated breakpoints, 
#'		fitted: fitted values radj: adjusted r value of the model
#' 		fit: fit object
#' @examples d1 <- rbind(c(rep(1,50),1:50), c(100:1))
#' 		rownames(d1) <- c("g1","g2")
#' 		fitseg(d1, "g1")
#' @author Ning Leng and Rhonda Bacher
#' @import stats
#' @import segmented


fitsegRsq <- function(Data,  T.Vect = NULL, Max.K = 5, 
                    Min.Num.In.Seg = 5, Pval.Cut = .1, 
                    Cut.Diff = .1, Num.Try = 100, 
                    Keep.Fit = FALSE) {

	t.use <- 1:length(Data)
	t.l <- length(Data)
	
  if(!is.null(T.Vect)) { 
		t.use <- T.Vect
		t.l <- T.Vect[length(T.Vect)]}
		
	step.r <- c(1:Max.K)
	dat.tmp <- Data
	seed.use <- 1

	# start with lm without bp
	lm1 <- lm(dat.tmp ~ t.use)
	lm.radj <- summary(lm1)$adj.r.squared
	lm.rsq <- summary(lm1)$r.squared
	lm.slp <- coef(lm1)[2]
	lm.fit <- fitted.values(lm1)
	lm.pval <- coef(summary(lm1))[2,4]
	lm.sign <- ifelse(lm.slp > 0, 1, -1)
	lm.sign[which(lm.pval > Pval.Cut)] <- 0
	lm.id.sign <- rep(lm.sign, length(t.use))

	# fit all possible breakpoints now:
	fit.l.0 <- sapply(1:length(step.r), function(j) {
		           i <- step.r[j]
		           lmseg.try <- suppressMessages(try(segmented::segmented(lm1, seg.Z = ~t.use, 
			                        psi = round(seq(1, t.l, length.out = i + 2)[2:(i + 1)]), 
			                         control = segmented::seg.control(seed = seed.use)), silent = T))
		           seed.use2 <- seed.use
		           while("try-error" %in% class(lmseg.try) & seed.use2 <= Num.Try) {
	 	               seed.use2 <- seed.use2 + 1
		               lmseg.try <- suppressMessages(try(segmented::segmented(lm1, seg.Z = ~t.use, 
			                               psi = round(seq(1, t.l, length.out = i + 2)[2:(i + 1)]), 
			                               control = segmented::seg.control(seed = seed.use2)), silent = T))}
		                                 out <- lmseg.try}, 
            simplify = F)

	isna <- which(sapply(fit.l.0, function(i) "try-error" %in% class(i)))

# if it is not solved in 100 trys then return lm results (if Num.Try=100)
	if (length(isna) == Max.K) {
		out <- list(id.sign = lm.id.sign, slp = lm.slp, slp.sign = lm.sign,
				        slp.pval = lm.pval, bp = NA, fitted = lm.fit, radj = lm.radj, fit = lm1)
		if (Keep.Fit == FALSE) { out <- out[1:7] }
		return(out) 
    break}


	fit.l <- fit.l.0
	if (length(isna) > 0) { # if one of step.r cant be fitted then remove it.
		fit.l <- fit.l.0[-isna]
		step.r <- step.r[-isna]}

	slp.l <- sapply(fit.l, slope, simplify = F)
	rsq <- sapply(fit.l, function(i) {summary(i)$r.squared})
	radj <- sapply(fit.l, function(i) {summary(i)$adj.r.squared})
	brk.l <- sapply(fit.l ,function(i) {i$psi[,2]}, simplify = F)
	id.l <- sapply(fit.l, function(i) {i$id.group}, simplify = F)



	if (length(step.r) > 1) {
		rsq.diff <- diff(rsq) # get differences in R^2
		rsq.whichmax <- which(rsq.diff > Cut.Diff) 
		radj.max <- max(radj)
		rsq.max <- max(rsq)
		if (length(rsq.whichmax) > 0) {
			# start with the largest difference in rsq value
			r.choose <- max(rsq.whichmax) + 1} 
		  # compare here using adjusted Rsq. 
		if (length(rsq.whichmax) == 0 & radj[1] > lm.radj) {
			r.choose <- 1} #if none larger than diff use smallest bp =1
		if (length(rsq.whichmax) == 0 & radj[1] <= lm.radj) {
	  	  	# if none of them satisfy, then take the linear, return values
	  		out <- list(id.sign=lm.id.sign, slp=lm.slp, slp.sign=lm.sign, 
					          slp.pval=lm.pval, bp=NA, fitted=lm.fit,radj=lm.radj,fit=lm1)
	  		if (Keep.Fit == FALSE) {out <- out[1:7]}
	  		return(out)
	  		break}

		# now make sure that optimal value satisfies having minimum number of segments, 
		# if it does not then decrease breakpoints.
		while ((min(table(id.l[[r.choose]])) < Min.Num.In.Seg) & r.choose > 1) {
				r.choose <- r.choose - 1}

		if(r.choose == 1 & (min(table(id.l[[r.choose]])) < Min.Num.In.Seg)) {
	  	  	# if 1 bp gives small segment, then take the linear
	  	  	out <- list(id.sign=lm.id.sign, slp=lm.slp, slp.sign=lm.sign, 
						           slp.pval=lm.pval, bp=NA, fitted=lm.fit,radj=lm.radj,fit=lm1)
	  		if(Keep.Fit == FALSE) {out <- out[1:7]}
	  		return(out)
	  	  break}
	}

	if (length(step.r) == 1) {
		r.choose <- 1
		radj.max <- radj
		rsq.max <- rsq}

	# finally decide if best BP model is better than linear

	if (radj.max < lm.radj | (rsq.max - lm.rsq) < Cut.Diff | min(table(id.l[[r.choose]])) < Min.Num.In.Seg){
		out <- list(id.sign=lm.id.sign, slp=lm.slp, slp.sign=lm.sign, 
					        slp.pval=lm.pval, bp=NA, fitted=lm.fit,radj=lm.radj,fit=lm1)
		if(Keep.Fit == FALSE) {out <- out[1:7]}
	    return(out) 
		# if lm is better 
		break}

	r.choose.ori <- step.r[r.choose] # actual k (before remove the NA fitting)

	fit.choose <- fit.l[[r.choose]]
	fv.choose <- fitted.values(fit.choose)
	bp.choose <- brk.l[[r.choose]]
  if(length(bp.choose) >= 1) {
    names(bp.choose) <- paste0("breakpoint", 1:length(bp.choose))}
	slp.choose <- slp.l[[r.choose]][[1]][,1]
	slp.t <- slp.l[[r.choose]][[1]][,3]
	slp.pval <- pt(-abs(slp.t), 1)
  if(length(slp.pval) >= 1) {
    names(slp.pval) <- paste0("pval", 1:length(slp.pval))}
	slp.sign <- ifelse(slp.t > 0, 1, -1)
	slp.sign[which(slp.pval > Pval.Cut)] <- 0
	id.choose <- id.l[[r.choose]]
	id.sign <- slp.sign[id.choose + 1]
	names(id.sign) <- colnames(Data)
	out = list(id.sign=id.sign, slp=slp.choose, slp.sign=slp.sign, 
				slp.pval=slp.pval, bp=bp.choose, fitted=fv.choose,
				radj=radj[r.choose],fit=fit.choose)
	if(Keep.Fit == FALSE) {out <- out[1:7]}

	return(out)
}

