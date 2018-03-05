library(shiny)

# Define server logic required to generate and plot a random distribution
server <- function(input, output) {

  dataInput <- reactive({
    sessionEnvir <- sys.frame()
    if (!is.null(input$f1)) load(input$f1$datapath, sessionEnvir)
  })
  output$datastr <- renderPrint({
    if (is.null(dataInput())){return("Nothin to show")} else{dataInput(); str(Wheat_BFR) }
  })

  output$contents <- renderTable({

    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.

    req(input$file1)

    df <- read.csv(input$file1$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote)

    if(input$disp == "head") {
      return(head(df))
    }
    else {
      return(df)
    }

  })


  observeEvent(input$plot,{
    if ( is.null(input$file)) return(NULL)
    inFile <- input$file
    file <- inFile$datapath
    # load the file into new environment and get it from there
    e = new.env()
    name <- load(file, envir = e)
    data <- e[[name]]

    # Plot the data
    output$hist <- renderPlot({
      hist(data)
    })
  })

}
