#' Trendy shiny app to interactively vizualize results after running trendy().
#'
#'
#' @return Opens a browser window with an interactive \code{shiny} 
#'  app and visualize all precomputed Trendy fits.
#' 
#' @name trendyShiny
#' @import shiny
#' @import utils
#' @import shinyFiles
#' @importFrom magrittr %>%
#' @export

# tVect <- NULL
# origData <- NULL
# trendyOut <- NULL
    
    globalVariables(c("tVect", "origData", "trendyOut"))
    
trendyShiny <- function() {
    
    options(shiny.maxRequestSize = 10000*1024^2) 
    
    
    ui <- fluidPage(
        #  Application title
        headerPanel("Trendy Visualization"),

        fluidRow(
            column(width = 7, 
                tags$div(tags$h4("This shiny app is designed to explore the 
                output from Trendy. First the .RData object 
                output from Trendy must be uploaded.")),
                tags$br(),
                tags$br(),
                fileInput("filename", label = "Input .Rdata from 
                    trendy() run:"),
                actionButton("Submit1","Upload File")
            )
        ),
        
        # Sidebar with sliders that demonstrate various available options
        conditionalPanel(condition = "input.Submit1 > 0",
        fluidRow( 
            column(width = 7, 
                tags$br(),
                tags$br(),
                tags$div(tags$h4("To visualize gene trends one by one, 
                use the 'Visualize genes' tab."),
                tags$h4("To obtain a list of genes according to a specific
                pattern use the 'Obtain gene patterns' tab.")),
                tags$br()
            )
        ),
        fluidRow(
            column(12,
                tabsetPanel(
                    tabPanel("Visualize genes", 
                    tags$br(),
                    column(6, align='left',
                    tags$div(tags$h4("Select a row in the table to update the
                    trend visualization."))
                    ),
                    column(10, align="center",
                    mainPanel(plotOutput('genePlot'), width = "100%"),
                    tags$br(),
                    DT::dataTableOutput("tab"),
                    tags$br(),
                    tags$br()
                    )
                    ),
                tabPanel("Obtain gene patterns",
                tags$br(),
                column(5,
                    tags$div(tags$b("Please select a folder for output :")),
                    shinyDirButton('Outdir', label ='Select Output Folder', 
                    title = 'Please select a folder'),
                    tags$br(),
                    tags$br(),
                    textInput("pattern", "Enter pattern (separate by 
                        comma, no spaces):", "up,down"),
                    textInput("rcut", "Only consider genes with an 
                        adjusted R squared greater than: ", ".5", 
                        placeholder="Must be between 0 and 1."),
                    textInput("delay", "Only consider genes with the pattern
                        occurring after time-point: ", "0"),
                    radioButtons("scatterplots",
                        label = "Output a plot of patterned genes?",
                        choices = list("Yes" = 1,"No" = 2),
                        selected = 1),
                    textInput("OutFileName", 
                        label = "Name of output files (default will
                            be the specified pattern)", value = ""),
                    tags$br(),
                    tags$br(),
                    actionButton("Submit","Submit for processing"),
                    tags$br(),
                    tags$br()),
                    column(6, textOutput("text1"),
                    tags$head(tags$style("#text1{color: red;font-size: 
                        26px;font-style: bold;}")))
                )
                )
            )
        )),
    fluidRow(mainPanel(""))
    )
    
    
server <- function(input, output, session) {
    
    raVar <- reactiveValues(genes.pass=NULL, outdir=NULL, plotVals=1)
    
    volumes <- c('home'="~")
    shinyDirChoose(input, 'Outdir', roots=volumes, session=session, 
    restrictions=system.file(package='base'))
    output$Dir <- renderPrint({parseDirPath(volumes, input$Outdir)})
    
    In <- eventReactive(input$Submit1, {
        load(input$filename$datapath)

        allTrendy <- Trendy::topTrendy(trendyOut, -1)
        ToPrint <- Trendy::formatResults(allTrendy)

        LIST = list(trendyOut, origData, tVect, ToPrint, allTrendy)
        names(LIST) <- c("trendyOut", "origData", "tVect", 
                "ToPrint", "TopTrendy")
        return(LIST)
    })
    
    observeEvent(input$Submit, {
        
        withProgress(message = 'Finding gene patterns:', value = 0, {
            
            IN <- In()
            outdir <- paste0("~", do.call("file.path", input$Outdir[[1]]), "/")
            
            pattern <- strsplit(input$pattern, ",")[[1]]
            delay <- as.numeric(input$delay)
            rcut <- as.numeric(input$rcut)
            
            incProgress(0.4, detail = "Extracting genes")
            genes.pass <- Trendy::extractPattern(IN$trendyOut, 
                Pattern = pattern, adjR2Cut = rcut, Delay = delay)
                
            if (input$OutFileName == "") {
                outfilename = paste0("Pattern_", input$pattern)
            } else {outfilename <- input$OutFileName}
            
            outfileP = paste0(outdir, outfilename, "_Scatter.pdf")
            outfileSS = paste0(outdir, outfilename, "_ShortSummary.csv")
            outfileLS = paste0(outdir, outfilename, "_FullSummary.csv")
            
            write.table(genes.pass, file = outfileSS, quote = FALSE,
                row.names = FALSE, sep = ",")
                
            toFormat <- genes.pass$Gene
            fullSummary <- IN$ToPrint[toFormat,]
            
            write.table(fullSummary, file = outfileLS, quote=FALSE, 
                row.names=FALSE, sep=",")
            incProgress(0.4, detail = "Writing genes to output folder")
            if (input$scatterplots == "1") {
                
                incProgress( 0, detail = "Making scatter plots of 
                    patterned genes")
                    
                pdf(outfileP, height=15, width=10)
                par(mar=c(5,5,2,1), mfrow=c(3,2))
                XX <- Trendy::plotFeature(Data = IN$origData, 
                    tVectIn = IN$tVect,
                    featureNames = genes.pass$Gene, showFit = TRUE,
                    trendyOutData = IN$trendyOut)
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
        
        MM <- paste(numG, "genes/features with pattern", input$pattern, 
            "have been output to", raVar$outdir )
        return(MM)
    })
    
    observeEvent(input$tab_rows_selected, {
        raVar$plotVals <- input$tab_rows_selected
    })
    
    output$genePlot <- renderPlot({
        IN <- In()
        topg <- as.character(IN$ToPrint[raVar$plotVals, 1])
        par(mfrow=c(1,1), cex=1.5, cex.lab=1, cex.axis=1, cex.main=1.1, 
        mar=c(5,5,2,2), oma=c(0,.1,.1,6))
        plot(IN$tVect, IN$origData[topg,], pch=20, col="#696969", 
            main=paste0(topg), ylab="Normalized Expr.", xlab="Time", 
            cex.axis=1.2, cex.lab=1.2)
        if (topg %in% names(IN$trendyOut)) {
            tmp <- IN$trendyOut[[topg]]
            lines(IN$tVect, tmp$Fitted.Values, lwd = 3, col="#ededed")
            abline(v = tmp$Breakpoints, lty = 2, lwd = 3, col="chartreuse3")
            ID <- tmp$Trends
            FIT <- tmp$Fitted.Values
            BKS <- c(0, tmp$Breakpoints, max(IN$tVect))
            if (length(BKS) > 3 | (length(BKS) == 3 & !is.na(BKS[2]))) {
                for (i in seq_len(length(BKS)+1)) {
                    toCol <- which(IN$tVect <= BKS[i+1] & IN$tVect >= BKS[i])
                    IDseg <- ID[toCol]
                    useCol <- switch(names(which.max(table(IDseg))), 
                    "0" = "black", 
                    "-1" = "cornflowerblue", 
                    "1" = "coral1")
                    lines(IN$tVect[toCol], FIT[toCol], lwd = 5, col=useCol)
                }
                par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), 
                    mar = c(0, 0, 4, 0), new = TRUE)
                plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
                legend("topright", c("Breakpoint"), xpd = TRUE, horiz = FALSE,
                    inset = c(.02,0), bty = "n", lty = c(2, 1, 1, 1), 
                    lwd = c(3,5,5,5), 
                    col = c("chartreuse3" ), cex = 1)
                par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), 
                    mar = c(0, 0, 7.5, 0),new = TRUE)
                plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
                legend("topright", c("Up", "No change", "Down"), xpd = TRUE, 
                horiz = FALSE, 
                inset = c(.02,0), bty = "n", title = "Segment trend:", 
                lty = c(1, 1, 1), 
                lwd = c(5,5,5), 
                col = c("coral1", "black","cornflowerblue"), cex = 1) 
            }
        }
        
    }, height=400, width=800)
    
    # Show the values using an HTML table
    output$tab <- DT::renderDataTable({
        IN <- In()
        toprint <- IN$ToPrint
        toprint[,-1] <- round(toprint[,-1], 3)
        
        numSeg <-colnames(toprint)[grepl("^Segment.*Trend",colnames(toprint))]
        numCols <- ncol(toprint)
        mkSmallTable <- c("Feature", "AdjustedR2", numSeg, 
        colnames(toprint)[grepl("Breakpoint", colnames(toprint))])
        
        toprint <- toprint[,mkSmallTable]
        
        COLS <- gsub(".", " ", colnames(toprint), fixed=TRUE)
        DT::datatable(toprint, rownames = FALSE, colnames = COLS, 
            selection = 'single',
            options = list(
                autoWidth = TRUE, scrollX=TRUE,
                columnDefs = list(list(className='dt-center',targets = '_all'))
                )) %>% DT::formatStyle(
                    columns = numSeg,
                    valueColumns = numSeg,
                    color = DT::styleEqual(c(-1, 0, 1), c('white', 'white', 
                    'white')),
                    backgroundColor = DT::styleEqual(c(-1, 0, 1), 
                    c('cornflowerblue', 'black', 'tomato')),
                    borderRightWidth = '5px', borderStyle = 'solid')
    })
    }
    runApp(shinyApp(ui = ui, server = server))
}