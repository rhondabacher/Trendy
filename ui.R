library(shiny)
library(shinyFiles)
#library(gdata)
options(shiny.maxRequestSize=10000*1024^2) 
# Define UI for slider demo application
shinyUI(fluidPage(
  #  Application title
  headerPanel("SegReg"),
  
  	 fluidRow(
		 column(width = 7, 
			 tags$br(),
				tags$div(tags$b("Upload .RData object first, then obtain list of genes according to some pattern or visualize genes one by one.")),
                tags$br(),
                tags$br(),
	
    fileInput("filename", label = "Input .Rdata from segreg() run:"),
	actionButton("Submit1","Upload File"),
	    tags$br(),
	    tags$br(),
		
	mainPanel(textOutput("InCheck")),
	br(),
	br()
	)),
	
  # Sidebar with sliders that demonstrate various available options
 conditionalPanel(condition = "input.Submit1 > 0",
	 fluidRow( column(8,
  tabsetPanel(
	  tabPanel("Obtain gene patterns", #width=22,height=20,
               # file
			tags$br(),
			  
   					tags$div(tags$b("Please select a folder for output :")),
                    shinyDirButton('Outdir', label ='Select Output Folder', title = 'Please select a folder'),
                    tags$br(),
                    tags$br(),
					
                      # Num permutation
                      textInput("pattern", "Enter pattern (separate by comma, no spaces):", "up,down"),
   				      textInput("rcut", "Only consider genes with adjusted R squared greater than: ", ".5", placeholder="Must be between 0 and 1."),
					  textInput("delay", "Only consider genes with pattern after time-point: ", "0"),
                      
                      # Normalization
                      radioButtons("scatterplots",
                                   label = "Output a plot of patterned genes?",
                                   choices = list("Yes" = 1,
                                                  "No" = 2),
                                   selected = 1),

                      textInput("OutFileName", 
                                label = "Output file name (will default to pattern)", 
                                value = ""),
               br(),
               br(),
               actionButton("Submit","Submit for processing"),
               tags$br(),
               tags$br()
		
  		),
		
		tabPanel("Visualize genes", 
		   tags$br(),
		   br(),
		   
		    textInput("gene", "Gene/Feature Name:", ""),
			
		    column(8,
		           # Show a plot of the generated distribution
		           mainPanel(plotOutput('genePlot'), width = "100%")))
	)))),
	 fluidRow(mainPanel(""))
))
