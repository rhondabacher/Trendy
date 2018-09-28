# Trendy
Trendy: Segmented regression analysis of expression dynamics for high-throughput ordered profiling experiments

Trendy utilizes segmented regression models to simultaneously characterize each gene’s expression pattern and summarize overall dynamic activity in ordered condition experiments. For each gene, Trendy finds the optimal segmented regression model and provides the location and direction of dynamic changes in expression. The top dynamic genes are then identified as genes that can be well profiled by its gene-specific segmented regression model. Trendy also implements functions to visualize the dynamic genes and their trends, to order dynamic genes by their trends, and to compute breakpoint distribution at different time points (e.g. detect time points with a large number of expression changes).

Details of Trendy may be found in the bioRxiv preprint:
http://www.biorxiv.org/content/early/2017/09/07/185413

The vignette for Trendy can be found here:
http://www.bioconductor.org/packages/devel/bioc/vignettes/Trendy/inst/doc/Trendy_vignette.pdf

#### The current version of Trendy is now on Bioconductor: http://www.bioconductor.org/packages/devel/bioc/html/Trendy.html

#### For previous versions check the release page.


#### To install Trendy:

You need to have R version 3.5 installed.

##### Option 1:

library(BiocManager)

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install("Trendy")


##### Option 2:

install.packages("devtools")

library(devtools)

install_github("rhondabacher/Trendy")

##### Option 3 (For R version 3.4):

install.packages("devtools")

library(devtools)

install_github("rhondabacher/Trendy", ref="devel")


## Trendy R/Shiny Visualization


Trendy R/Shiny assumes you have already run the trendy() function in the Trendy package and saved the output as an .RData object (by setting saveObject = TRUE). The app allows you to extract lists of genes/features according to any pattern of interest. The patten of interest can also be extracted after a given condition (time-point) via the Delay option. The Shiny application also enables the trend and breakpoints of each gene to be explored and visualized.


#### To launch the Shiny app in R:

library(Trendy)

trendyShiny()


Example of patterns:

"up,down" : Genes/Features which contain a peak through the time-course.

"down,up" : Genes/Features which contain a trough in the time-course.


