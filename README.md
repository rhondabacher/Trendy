# Trendy
Trendy: Segmented regression analysis of expression dynamics for high-throughput ordered profiling experiments

Trendy utilizes segmented regression models to simultaneously characterize each gene’s expression pattern and summarize overall dynamic activity in ordered condition experiments. For each gene, Trendy finds the optimal segmented regression model and provides the location and direction of dynamic changes in expression. The top dynamic genes are then identified as genes that can be well profiled by its gene-specific segmented regression model. Trendy also implements functions to visualize the dynamic genes and their trends, to order dynamic genes by their trends, and to compute breakpoint distribution at different time points (e.g. detect time points with a large number of expression changes).

Details of Trendy may be found in the vignette:
https://github.com/rhondabacher/Trendy/tree/master/package/Trendy/vignettes/Trendy_vignette.pdf

and in the bioRxiv preprint:
http://www.biorxiv.org/content/early/2017/09/07/185413

#### The current version of Trendy is in submission to Bioconductor. For previous versions check the release page.


#### To install Trendy:

install.packages("devtools")

library(devtools)

install_github("rhondabacher/Trendy/package/Trendy")


## Trendy R/Shiny Visualization

Trendy R/Shiny assumes you have already run the trendy() function in the Trendy package and saved the output as an .RData object (by setting Save.Object = TRUE). The app allows you to extract lists of genes/features according to any pattern of interest. The patten of interest can also be extracted after a given condition (time-point) via the delay option. The Shiny application also enables the trend and breakpoints of each gene to be explored and visualized.

#### To launch the Shiny app in R:

library(Trendy)

trendyShiny()


Example of patterns:

"up,down" : Genes/Features which contain a peak through the time-course.

"down,up" : Genes/Features which contain a trough in the time-course.


