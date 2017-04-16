# Trendy
Trendy: Breakpoint analysis of time course expression data

Trendy is an R package that can be used to perform breakpoint analysis on microarray or RNA-seq expression data 
with ordered conditions (e.g. time course, spatial course). For each gene or other features, Trendy estimates the optimal number of 
breakpoints, as well as, the breakpoints by fitting a set of segmented regression models. The top dynamic genes are then identified 
by taking genes that can be well profiled by its gene-specific segmented regression model. Trendy also implements functions to visualize 
the dynamic genes and their trends, to order dynamic genes by their trends, and to compute breakpoint distribution at different 
time points (e.g. detect time points with a large number of expression changes).

Details of Trendy may be found in the vignette:
https://github.com/rhondabacher/Trendy/package/Trendy/vignettes/Trendy_vignette.pdf


## Trendy shiny

Trendy shiny assumes you have already run the Trendy function and saved the output as an RData object. 
The app allows you to extract lists of genes/features according to any pattern of interest. 
The patten of interest can also be extracted after a given condition(time-point).

####Make sure packages are installed:

install.packages(c("segmented","gplots")) 

library(devtools)

install_github("rhondabacher/Trendy/package/Trendy")

install.packages("shiny")


####To launch the Shiny app:

runGitHub('rhondabacher/Trendy')

Example of patterns:

"up,down" : Genes/Features which contain a peak through the time-course.

"down,up" : Genes/Features which contain a trough in the time-course.

"up" : Genes/Features which constantly increase throughout the time-course.

