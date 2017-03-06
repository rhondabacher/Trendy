library(SegReg)
library(shiny)

shinyServer(function(input, output, session) {
  
	volumes <- c('home'="~")
    shinyDirChoose(input, 'Outdir', roots=volumes, session=session, restrictions=system.file(package='base'))
    output$Dir <- renderPrint({parseDirPath(volumes, input$Outdir)})



	 
    In <- eventReactive(input$Submit1, {
		

      (load(input$filename$datapath))

	  # print(ls())
	  
	  List = c(seg.object, orig.data, t.vect)
	 })
	 

	 
	output$InCheck <- renderText({
	
 		 In <- In()

 	  		if(!exists("seg.object") | !exists("orig.data") | !exists("t.vect")) {
				"Segmented Regression object invalid!"
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
		   
		   
		   genes.pass <- extractpattern(seg.object, radjcut = rcut, pattern = pattern, delay=delay)
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
	  	  XX <- plotmarker(orig.data, t.vect=t.vect, outfile = outfileP, listfeatures = genes.pass$Gene, savePDF = TRUE, 
			  plotFit = TRUE, fittedreg = seg.object, yname="Normalized Expression", par.param = c(3,2), pdfheight=15, pdfwidth=10)
	   
	  }
	  

	   incProgress(.2, detail = "Done!")
  	})

})
   
   
   output$genePlot <- renderPlot({
	   In()
   	  req(input$gene)
	   if(input$gene == "") {topg <- rownames(orig.data)[1]} else {topg <- input$gene}
		   
       par(mfrow=c(1,2), cex=1.5, cex.lab=1, cex.axis=1, cex.main=1.1, mar=c(6,6,2,2))

       plot(t.vect, orig.data[topg,], pch=20, col="red", main=paste0(topg),
            ylab="Normalized Expr.", xlab="Time")
		if(topg %in% names(seg.object)) {
			tmp <- seg.object[[topg]]
			lines(t.vect, tmp$fitted, lwd = 3)
		}
		
       plot(t.vect, log2(orig.data[topg,]+1), pch=20, col="red", main=paste0(topg),
 	             ylab="Log2 Normalized Expr.", xlab="Time")

				 
	}, height=500, width=1000)
	
    
})