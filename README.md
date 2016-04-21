# SegReg
SegReg: Breakpoint analysis of time course expression data

SegReg is a R package that can be used to perform breakpoint analysis on Microarrays or RNA-seq expression data 
with ordered conditions (e.g. time course, spatial course). For each gene or other features, SegReg estimates the optimal number of 
breakpoints as well as the breakpoints by fitting a set of segmented regression models. The top dynamic genes are then identified 
by taking genes that can be well profiled by its gene-specific segmented regression model. SegReg also implements functions to visualize 
the dynamic genes and their trends, to order dynamic genes by their trends, and to compute breakpoint distribution at different 
time points (e.g. detect time points with a large number of expression changes).

Details of SegReg may be found in the vignette:
https://github.com/lengning/SegReg/blob/master/SegReg_vignette.pdf

To download the vignette:
https://github.com/lengning/SegReg/raw/master/SegReg_vignette.pdf
