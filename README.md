# Trendy
Trendy: Segmented regression analysis of expression dynamics for high-throughput ordered profiling experiments

Trendy is an R package that can be used to perform breakpoint analysis on microarray or RNA-seq expression data 
with ordered conditions (e.g. time-course, spatial-course). For each gene or other features, Trendy estimates the optimal number of breakpoints, as well as, the breakpoints by fitting a set of segmented regression models. The top dynamic genes are then identified by taking genes that can be well profiled by its gene-specific segmented regression model. Trendy also implements functions to visualize the dynamic genes and their trends, to order dynamic genes by their trends, and to compute breakpoint distribution at different time points (e.g. detect time points with a large number of expression changes).

Details of Trendy may be found in the vignette:
https://github.com/rhondabacher/Trendy/tree/master/package/Trendy/vignettes/Trendy_vignette.pdf


## Trendy R/Shiny

Trendy R/Shiny assumes you have already run the trendy() function in the Trendy package and saved the output as an .RData object (by setting Save.Object = TRUE). The app allows you to extract lists of genes/features according to any pattern of interest. 
The patten of interest can also be extracted after a given condition (time-point) via the delay option.

####  Required packages to install:

install.packages(c("segmented","gplots")) 

install.packages("devtools")

library(devtools)

install_github("rhondabacher/Trendy/package/Trendy")

install.packages("shiny")

install.packages("shinyFiles")

#### To launch the Shiny app:

library(shiny)

library(Trendy)

runGitHub('rhondabacher/Trendy')

Example of patterns:

"up,down" : Genes/Features which contain a peak through the time-course.

"down,up" : Genes/Features which contain a trough in the time-course.


