library(Trendy)
library(shiny)

shinyServer(function(input, output, session) {
  
  raVar <- reactiveValues(genes.pass=NULL, outdir=NULL)
  
  volumes <- c('home'="~")
  shinyDirChoose(input, 'Outdir', roots=volumes, session=session, restrictions=system.file(package='base'))
  output$Dir <- renderPrint({parseDirPath(volumes, input$Outdir)})
   
  In <- eventReactive(input$Submit1, {
    
    load(input$filename$datapath)
    
    if(exists("Seg.Object")) {
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
  
  output$choose_gene <- renderUI({

      IN <- In()

      # Get the data set with the appropriate name
      allg <- rownames(IN$Orig.Data)

      # Create the checkboxes and select them all by default
      selectInput("gene", "Gene/Feature Name:",
                          choices  = allg,
                          selected = allg[1])
    })


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
      
      
    output$genePlot <- renderPlot({

     IN <- In()

      req(input$gene)

      if(input$gene == "") {
        topg <- rownames(IN$Orig.Data)[1]
        } else {topg <- input$gene}

    
      par(mfrow=c(1,1), cex=1.5, cex.lab=1, cex.axis=1, cex.main=1.1, mar=c(6,6,2,2))

      plot(IN$T.Vect, IN$Orig.Data[topg,], pch=20, col="red", main=paste0(topg),
           ylab="Normalized Expr.", xlab="Time")
      if(topg %in% names(IN$Trendy.Out)) {
        tmp <- IN$Trendy.Out[[topg]]
        lines(IN$T.Vect, tmp$Fitted.Values, lwd = 3)
      }

    }, height=400, width=600)

  # Show the values using an HTML table
  output$tab <- renderDataTable({
      IN <- In()

      req(input$gene)
      if(input$gene == "") {
        topg <- rownames(In$Orig.Data)[1]
        } else {topg <- input$gene}

      show.order <- c(topg, setdiff(rownames(IN$To.Print), topg))
      whichCol <- !grepl(".Fitted.Trend", colnames(IN$To.Print))
      toprint <- IN$To.Print[show.order, whichCol]
      toprint
    })

})