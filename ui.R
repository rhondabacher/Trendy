library(shiny)
library(shinyFiles)
library(DT)
#library(gdata)
options(shiny.maxRequestSize=10000*1024^2) 
# Define UI for slider demo application
shinyUI(fluidPage(
  #  Application title
  headerPanel("Trendy Visualization"),
  
 fluidRow(
		 column(width = 7, 
				tags$div(tags$h4("This shiny app is designed to explore the output from Trendy. 
                          First the .RData object output from Trendy must be uploaded."),
                 tags$br()
                 ),
     tags$br(),
	
     fileInput("filename", label = "Input .Rdata from trendy() run:"),
	   actionButton("Submit1","Upload File")
	   
  )),
	
  # Sidebar with sliders that demonstrate various available options
 conditionalPanel(condition = "input.Submit1 > 0",
	 fluidRow( 
    column(width = 7, 
      tags$br(),
      tags$br(),
  		tags$div(
                tags$h4("To obtain a list of genes according to a specific  
                         pattern use the 'Obtain gene patterns' tab."),
                tags$h4("To visualize gene trends one by one, use the 'Visualize genes' tab.")
               ),
        tags$br()
      )),
  fluidRow(    
     column(12,
     tabsetPanel(
	       tabPanel("Obtain gene patterns", #width=22,height=20,
               # file
		             tags$br(),
			  
   				column(5,	
            tags$div(tags$b("Please select a folder for output :")),
            shinyDirButton('Outdir', label ='Select Output Folder', title = 'Please select a folder'),
            tags$br(),
            tags$br(),
					
            
            textInput("pattern", "Enter pattern (separate by comma, no spaces):", "up,down"),
   				  textInput("rcut", "Only consider genes with an adjusted R squared greater than: ", ".5", 
                                      placeholder="Must be between 0 and 1."),
					  textInput("delay", "Only consider genes with the pattern occurring after time-point: ", "0"),
                      
            
            radioButtons("scatterplots",
                         label = "Output a plot of patterned genes?",
                         choices = list("Yes" = 1,
                                        "No" = 2),
                         selected = 1),

            textInput("OutFileName", 
                      label = "Name of output files (default will be the specified pattern)", 
                      value = ""),
               br(),
               br(),
               actionButton("Submit","Submit for processing"),
               tags$br(),
               tags$br()
              ),
              column(6, textOutput("text1"),
                  tags$head(tags$style("#text1{color: red;
                                 font-size: 26px;
                                 font-style: bold;
                                 }"
                         )
              )
                )
        ),
		
	tabPanel("Visualize genes", 
	   tags$br(),
     # br(),
   
      # uiOutput("choose_gene"),
      column(6, align='left',
  		tags$div(
                tags$h4("Select a row in the table to update the trend visualization."))
      ),
	    column(10, align="center",
	           mainPanel(plotOutput('genePlot'), width = "100%"),
             tags$br(),
             DT::dataTableOutput("tab"),
             tags$br(),
             tags$br()
  ))
)))),
fluidRow(mainPanel(""))
))
