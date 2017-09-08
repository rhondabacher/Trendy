library(Trendy)
library(shiny)
library(DT)
shinyServer(function(input, output, session) {
  
  raVar <- reactiveValues(genes.pass=NULL, outdir=NULL, plotVals=1)
  
  volumes <- c('home'="~")
  shinyDirChoose(input, 'Outdir', roots=volumes, session=session, restrictions=system.file(package='base'))
  output$Dir <- renderPrint({parseDirPath(volumes, input$Outdir)})
   
  In <- eventReactive(input$Submit1, {
    
    load(input$filename$datapath)
     
    if(!is.na(match("Seg.Object", ls()))) {

            Trendy.Out  <- Seg.Object
            Trendy.Out <- lapply(Trendy.Out, function(x) {
                         
                                  names(x[[1]]) <- paste0(names(T.Vect), ".Trend")
                                  names(x[[2]]) <- paste0("Segment", seq_len(length(x[[2]])), ".Slope")
                                  names(x[[3]]) <- paste0("Segment", seq_len(length(x[[2]])), ".Trend")
                                  names(x[[4]]) <- paste0("Segment", seq_len(length(x[[2]])), ".Pvalue")
                                  names(x[[5]]) <-  paste0("Breakpoint", 1:length(x[[5]]))
                                  names(x[[6]]) <- paste0(names(T.Vect), ".Fitted")
                             
                                 return(x) 
                                })
                              
            Trendy.Out <- lapply(Trendy.Out, function(x) {
                                  names(x) <- c("Trends", 
                                  "Segment.Slopes", "Segment.Trends", "Segment.Pvalues",
                                  "Breakpoints", "Fitted.Values", "AdjustedR2")
                                  return(x)})
      }

    all.trendy <- topTrendy(Trendy.Out, -1)
    To.Print <- Trendy::formatResults(all.trendy)
    
    LIST = list(Trendy.Out, Orig.Data, T.Vect, To.Print, all.trendy)
    names(LIST) <- c("Trendy.Out", "Orig.Data", "T.Vect", "To.Print", "TopTrendy")
    return(LIST)
  })
  
  # output$choose_gene <- renderUI({
#
#       IN <- In()
#
#       # Get the data set with the appropriate name
#       allg <- sort(rownames(IN$Orig.Data))
#
#       # Create the checkboxes and select them all by default
#       selectInput("gene", "Gene/Feature Name:",
#                           choices  = allg,
#                           selected = allg[1])
#     })


   observeEvent(input$Submit, {

      withProgress(message = 'Finding gene patterns:', value = 0, {

        IN <- In()

        outdir <- paste0("~", do.call("file.path", input$Outdir[[1]]), "/")

        pattern <- strsplit(input$pattern, ",")[[1]]

        delay <- as.numeric(input$delay)
        rcut <- as.numeric(input$rcut)

        incProgress(0.4, detail = "Extracting genes")

        genes.pass <- extractPattern(IN$Trendy.Out, Pattern = pattern, AdjR2.Cut = rcut, Delay = delay)

        if(input$OutFileName == "") {
          outfilename = paste0("Pattern_", input$pattern)
        }else {outfilename <- input$OutFileName}

        outfileP = paste0(outdir, outfilename, "_Scatter.pdf")
        outfileSS = paste0(outdir, outfilename, "_ShortSummary.csv")
        outfileLS = paste0(outdir, outfilename, "_FullSummary.csv")


        write.table(genes.pass, file = outfileSS, quote=F, row.names=FALSE, sep=",")

        toFormat <- genes.pass$Gene
        fullSummary <- IN$To.Print[toFormat,]

        write.table(fullSummary, file = outfileLS, quote=F, row.names=FALSE, sep=",")


        incProgress(0.4, detail = "Writing genes to output folder")

        if(input$scatterplots == "1") {

          incProgress( 0, detail = "Making scatter plots of patterned genes")

          pdf(outfileP, height=15, width=10)
          par(mar=c(5,5,2,1), mfrow=c(3,2))
          XX <- plotFeature(Data = IN$Orig.Data, T.Vect = IN$T.Vect,
                            Feature.Names = genes.pass$Gene, Show.Fit = TRUE,
                            Trendy.Out = IN$Trendy.Out)
          dev.off()
        }

        incProgress(.2, detail = "Done!")
        raVar$outdir = outdir
        raVar$genes.pass = genes.pass
      })
    })

   output$text1 <- renderText({
    if (is.null(raVar$genes.pass)) return() 
    numG <- length(raVar$genes.pass$Gene)

    MM <- paste(numG, "genes/features with pattern", input$pattern, "have been output to", raVar$outdir )
    return(MM)
     
   })

   observeEvent(input$tab_rows_selected, {
     
       raVar$plotVals <- input$tab_rows_selected
     
     })
       
          
    output$genePlot <- renderPlot({

     IN <- In()

    topg <- as.character(IN$To.Print[raVar$plotVals, 1])

    par(mfrow=c(1,1), cex=1.5, cex.lab=1, cex.axis=1, cex.main=1.1, mar=c(5,5,2,2), oma=c(0,.1,.1,6))

      plot(IN$T.Vect, IN$Orig.Data[topg,], pch=20, col="#696969", main=paste0(topg),
           ylab="Normalized Expr.", xlab="Time", cex.axis=1.2, cex.lab=1.2)
      if(topg %in% names(IN$Trendy.Out)) {
        tmp <- IN$Trendy.Out[[topg]]
        lines(IN$T.Vect, tmp$Fitted.Values, lwd = 3, col="#ededed")
        abline(v = tmp$Breakpoints, lty = 2, lwd = 3, col="chartreuse3")
         ID <- tmp$Trends
         FIT <- tmp$Fitted.Values
          BKS <- c(0, tmp$Breakpoints, max(IN$T.Vect))
          if(length(BKS) > 3 | (length(BKS) == 3 & !is.na(BKS[2]))) {
              for(i in 1:(length(tmp$Breakpoints)+1)) {
                 toCol <- which(IN$T.Vect <= BKS[i+1] & IN$T.Vect >= BKS[i])
                 IDseg <- ID[toCol]
                 useCol <- switch(names(which.max(table(IDseg))), 
                             "0" = "black", 
                             "-1" = "cornflowerblue", 
                             "1" = "coral1")
                lines(IN$T.Vect[toCol], FIT[toCol], lwd = 5, col=useCol)
               }
               par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 4, 0), new = TRUE)
               plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

               legend("topright", c("Breakpoint"), xpd = TRUE, horiz = FALSE, 
                 inset = c(.02,0), bty = "n", 
                 lty = c(2, 1, 1, 1), lwd = c(3,5,5,5), col = c("chartreuse3" ), cex = 1) 
           
      
               par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 7.5, 0), new = TRUE)
               plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

               legend("topright", c("Up", "No change", "Down"), xpd = TRUE, horiz = FALSE, 
                 inset = c(.02,0), bty = "n", title = "Segment trend:",
                 lty = c(1, 1, 1), lwd = c(5,5,5), col = c("coral1", "black","cornflowerblue"  ), cex = 1) 
           }
         }

    }, height=400, width=800)

  # Show the values using an HTML table
  output$tab <- DT::renderDataTable({
      IN <- In()

        
      # whichCol <- !grepl(".Fitted.Trend", colnames(IN$To.Print))
      toprint <- IN$To.Print
      toprint[,-1] <- round(toprint[,-1], 3)

      numSeg <- colnames(toprint)[grepl("^Segment.*Trend", colnames(toprint))]
      numCols <- ncol(toprint)
      mkSmallTable <- c("Feature",  "AdjustedR2", numSeg, colnames(toprint)[grepl("Breakpoint", colnames(toprint))])
      
      
      toprint <- toprint[,mkSmallTable]
     
      COLS <- gsub(".", " ", colnames(toprint), fixed=TRUE)
      DT::datatable(toprint, rownames = FALSE, colnames = COLS, selection = 'single',
        options = list(
          autoWidth = TRUE, scrollX=TRUE,
          columnDefs = list(list(className = 'dt-center', targets = '_all'))
          )) %>% formatStyle(
                  columns = numSeg,
                   valueColumns = numSeg,
                   color = styleEqual(c(-1, 0, 1), c('white', 'white', 'white')),
                   backgroundColor = styleEqual(c(-1, 0, 1), c('cornflowerblue', 'black', 'tomato')),
                   borderRightWidth = '5px', borderStyle = 'solid'
            )
    })


})