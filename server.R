library(Trendy)
library(shiny)

shinyServer(function(input, output, session) {
  
  volumes <- c('home'="~")
  shinyDirChoose(input, 'Outdir', roots=volumes, session=session, restrictions=system.file(package='base'))
  output$Dir <- renderPrint({parseDirPath(volumes, input$Outdir)})
  
  
  
  
  In <- eventReactive(input$Submit1, {
    
    
    load(input$filename$datapath)
    
    print(ls())
    
    List = c(Seg.Object, Orig.Data, T.Vect)
  })
  
  
  
  output$InCheck <- renderText({
    
    In <- In()
    
    if(!exists("Seg.Object") | !exists("Orig.Data") | !exists("T.Vect")) {
      "Trendy object invalid!"
    } else {"File is uploaded!"}
    
  })
  
  
  OUT1 <- observeEvent(input$Submit, {
    
    
    withProgress(message = 'Finding gene patterns:', value = 0, {
      
      In <- In()
      
      outdir <- paste0("~", do.call("file.path", input$Outdir[[1]]), "/")
      
      
      pattern <- strsplit(input$pattern, ",")[[1]]
      
      delay <- as.numeric(input$delay)
      rcut <- as.numeric(input$rcut)
      
      incProgress(0.4, detail = "Extracting genes")
      
      
      genes.pass <- extractpattern(Seg.Object, Radj.Cut = rcut, Pattern = pattern, Delay = delay)
      # print("X")
      
      
      
      if(input$OutFileName == "") {
        outfilename = paste0("genes_with_", input$pattern)
      }else {outfilename <- input$OutFileName}
      
      
      outfileP = paste0(outdir, outfilename, "_scatter.pdf")
      outfileG = paste0(outdir, outfilename, ".csv")
      
      
      write.table(genes.pass, file=outfileG, quote=F, row.names=FALSE, sep=",")
      
      incProgress(0.4, detail = "Writing genes to output folder")
      
      
      
      if(input$scatterplots == "1") {
        
        incProgress( 0, detail = "Making scatter plots of patterned genes")
        XX <- plotmarker(Data=Orig.Data, T.Vect=T.Vect, File.Name = outfileP, Feature.Names = genes.pass$Gene, PDF = TRUE, Seg.Fit = TRUE,
                         Seg.Data = Seg.Object, Y.Name="Normalized Expression", Par.Param = c(3,2), PDF.Height=15, PDF.Width=10)
        
      }
      
      
      incProgress(.2, detail = "Done!")
    })
    
  })
  
  
  output$genePlot <- renderPlot({
    In <- In()
    req(input$gene)
    if(input$gene == "") {
      topg <- rownames(Orig.Data)[1]
      } else {topg <- input$gene}
    
    par(mfrow=c(1,2), cex=1.5, cex.lab=1, cex.axis=1, cex.main=1.1, mar=c(6,6,2,2))
    
    plot(T.Vect, Orig.Data[topg,], pch=20, col="red", main=paste0(topg),
         ylab="Normalized Expr.", xlab="Time")
    if(topg %in% names(Seg.Object)) {
      tmp <- Seg.Object[[topg]]
      lines(T.Vect, tmp$fitted, lwd = 3)
    }
    
    plot(T.Vect, log2(Orig.Data[topg,]+1), pch=20, col="red", main=paste0(topg),
         ylab="Log2 Normalized Expr.", xlab="Time")
    
    
  }, height=500, width=1000)
  
  
})