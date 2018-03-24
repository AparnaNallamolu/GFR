
# Define UI for data upload app ----
fluidPage(
  # App title ----
  titlePanel("Bayesian Functional Regression"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----

    sidebarPanel(
      # # Input: Select a file ----
      fileInput('File1', 'Choose a RDS File with the Dataset ', accept = c('.RDS', '.Rds', '.rds')),

      # fileInput('File2', '(Optional) Choose a Genomic Matrix RDS File', accept = c('.RDS', '.Rds', '.rds')),
      #
      # fileInput('File3', '(Optional) Choose a Bands RDS File', accept = c('.RDS', '.Rds', '.rds')),
      #
      # fileInput('File4', '(Optional) Choose a Wavelengths RDS File', accept = c('.RDS', '.Rds', '.rds')),

      # Horizontal line ----
      tags$hr(),

      # Select the LineID
      selectInput('DatasetID', 'Choose the datasetID:', c(''), multiple = FALSE, selectize = FALSE),

      # Select the LineID
      selectInput('response', 'Response Type:', c('gaussian', 'ordinal'), multiple = FALSE, selectize = FALSE),

      # Select numbers
      numericInput("iter", "Iterations:", min = 1, max = 50000, value = 1000, step = 1),
      numericInput("burn", "Burning:", min = 1, max = 49999, value = 700, step = 1),
      numericInput("thin", "Thin:", min = 1, max = 1000, value = 1, step = 1),

      # Select seed
      numericInput("seed", "Select seed:", min = 1, max = 1000, value = 123, step = 1),

      tags$hr(),

      # Select Crossvalidation
      #
      selectInput('CrossV', 'CrossValidation Type', c('KFold', 'RandomPartition'), multiple = FALSE, selectize = TRUE),
      numericInput("nRep", "Folds:", min = 1, max = 20, value = 3, step = 1),
      # This outputs the dynamic UI component
      wellPanel(
      uiOutput("ui")),

      actionButton("proButton", "Process the data")
      # submitButton("Process data")
    ),

    # Main panel for displaying outputs ----
    mainPanel(
      h3("Design used: "),
      uiOutput('Design'),
      h3("Analysis Results"),
      plotOutput("plotres"),
      br(),
      dataTableOutput('res'),

      # Output: Data file ----
      h3("Dataset Information"),
      dataTableOutput('contents')
    )
  )
)

