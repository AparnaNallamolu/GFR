library(shiny)
# Define server logic required to generate and plot a random distribution
function(input, output, session) {

  dataSet <- NULL
  pm <- NULL
  ETA <- NULL

  design <- list('Bayes-Single' = '$$\\eta = NULL$$',
                 'Bayes-SingleBands' = '$$\\eta = Bands$$',
                 'Bayes-Env' = '$$\\eta = Env + Line + Env \\times Line$$',
                 'Bayes-EnvBands' = '$$\\eta = Env + Line + Env \\times Line + Bands$$',
                 'Bayes-Trait' = '$$\\eta = Trait + Line + Trait \\times Line$$',
                 'Bayes-TraitBands' = '$$\\eta = Trait + Line + Trait \\times Line + Bands$$',
                 'Bayes-Multi' = '$$\\eta = Env + Trait + Line + Line \\times Env + Line \\times Trait + Env\\times Line\\times Trait $$',
                 'Bayes-MultiBands' = '$$\\eta = Env + Trait + Line + Line \\times Env + Line \\times Trait + Env\\times Line\\times Trait + Bands$$')

  ############# SHINY APP STARTS HERE
  #############
  output$ui <- renderUI({
    # Depending on input$input_type, we'll generate a different
    # UI component and send it to the client.
    switch(input$CrossV,
           "KFold" = {updateNumericInput(session, 'nRep' , label = "Folds:")
             NULL},
           "RandomPartition" = {
             updateNumericInput(session, 'nRep' ,label = "Partitions:")
             numericInput("nPer", "Testing %:", min = 0, max = 1, value = .25, step = .05)
             }
    )
  })

  output$CrossV_text <- renderText({
    input$CrossV
  })
  output$AnalysisType_text <- renderText({
    input$AnalysisType
  })

  output$contents <- renderDataTable({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.

    req(input$File1)

    dataSet <<- readRDS(input$File1$datapath)
    dataSetLabels <- colnames(dataSet)
    s_options <- as.list(dataSetLabels)

    updateSelectInput(session, "DatasetID",
                    choices = s_options)
    return(dataSet)
  })

  observeEvent(input$proButton,{
    output$res = renderDataTable({
      if (input$CrossV == 'KFold') {
        Crossvalidation_list <- list(Type = input$CrossV, nFolds = input$nRep)
      } else if (input$CrossV == 'RandomPartition') {
        Crossvalidation_list <- list(Type = input$CrossV, NPartitions = input$nRep, PTesting = input$nPer)
      }

      ETA <<- GFR::ETAGenerate(dataSet, datasetID = input$DatasetID, priorType = 'BRR')
      pm <<- GFR::BFR(ETA = ETA, nIter = input$iter, burnIn = input$burn, set_seed = input$seed, CrossValidation = Crossvalidation_list, verbose = T)
      cleanDat(forceClean = TRUE)
      summary(pm)
    })
    output$Design <- renderUI({
      print(ETA$Design)
      withMathJax(design[[which(names(design) == ETA$Design)]])
    })
    output$plotres <- renderPlot({boxplot(pm)
      })
  })
}
