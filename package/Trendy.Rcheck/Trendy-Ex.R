pkgname <- "Trendy"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Trendy')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Orig.Data")
### * Orig.Data

flush(stderr()); flush(stdout())

### Name: Orig.Data
### Title: Example dataset for Trendy Shiny App
### Aliases: Orig.Data
### Keywords: datasets

### ** Examples

data(exampleObject_trendyForShiny)



cleanEx()
nameEx("TrendyExData")
### * TrendyExData

flush(stderr()); flush(stdout())

### Name: TrendyExData
### Title: Example dataset for Trendy
### Aliases: TrendyExData
### Keywords: datasets

### ** Examples

data(TrendyExData)



cleanEx()
nameEx("extractpattern")
### * extractpattern

flush(stderr()); flush(stdout())

### Name: extractpattern
### Title: Extract pattern from segmented regression
### Aliases: extractpattern

### ** Examples

myTrends <- trendy(TrendyExData)
extractpattern(myTrends, Pattern = c("up")) #increasing only features
 extractpattern(myTrends, Pattern = c("up", "down")) #features with a peak
 extractpattern(myTrends, Pattern = c("up", "down"), Delay = 20)



cleanEx()
nameEx("fitSegBIC")
### * fitSegBIC

flush(stderr()); flush(stdout())

### Name: fitSegBIC
### Title: Fit Segmented regression models on a feature/gene
### Aliases: fitSegBIC

### ** Examples

d1 <- rbind(c(rep(1,50),1:50), c(100:1))
		rownames(d1) <- c("g1","g2")
		fit.seg(d1, "g1")



cleanEx()
nameEx("fitsegRsq")
### * fitsegRsq

flush(stderr()); flush(stdout())

### Name: fitsegRsq
### Title: Fit Segmented regression models on a feature/gene
### Aliases: fitsegRsq

### ** Examples

d1 <- rbind(c(rep(1,50),1:50), c(100:1))
		rownames(d1) <- c("g1","g2")
		fitseg(d1, "g1")



cleanEx()
nameEx("plotfeature")
### * plotfeature

flush(stderr()); flush(stdout())

### Name: plotfeature
### Title: Plot features of interest
### Aliases: plotfeature

### ** Examples

d1 <- rbind(c(rep(1,50),1:50), c(100:1))
	rownames(d1) <- c("g1","g2")
	plotmarker(d1, Feature.Names=c("g1","g2"))



cleanEx()
nameEx("toptrendy")
### * toptrendy

flush(stderr()); flush(stdout())

### Name: toptrendy
### Title: obtain top genes from trendy results
### Aliases: toptrendy

### ** Examples

d1 <- rbind(c(rep(1,50),1:50), c(100:1))
     rownames(d1) <- c("g1","g2")
     seg.all <- trendy(d1)
		top.seg <- toptrendy(seg.all)



cleanEx()
nameEx("trendy")
### * trendy

flush(stderr()); flush(stdout())

### Name: trendy
### Title: segmented regression on a set of genes
### Aliases: trendy

### ** Examples

d1 <- rbind(c(rep(1,50),1:50), c(100:1))
     rownames(d1) <- c("g1","g2")
     trendy(d1)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
