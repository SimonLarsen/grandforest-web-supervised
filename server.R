library(shiny)
library(shinyjs)
library(shinysky)
library(visNetwork)
library(data.table)
library(magrittr)
library(grandforest)
library(igraph)
library(ComplexHeatmap)
library(circlize)
library(org.Hs.eg.db)
library(ggplot2)

source("grandforest-web-common/get_network.R")
source("grandforest-web-common/enrichment.R")
source("grandforest-web-common/feature_graph.R")
source("evaluation.R")

shinyServer(function(input, output, session) {
  currentEdges <- reactiveVal()
  currentData <- reactiveVal()
  currentModel <- reactiveVal()
  currentEnrichmentTable <- reactiveVal()
  currentEvaluation <- reactiveVal()
  currentPredictions <- reactiveVal()
  
  output$hasModel <- reactive({
    req(currentModel())
    return(TRUE)
  })
  outputOptions(output, "hasModel", suspendWhenHidden=FALSE)
  
  output$hasEnrichmentTable <- reactive({
    D <- req(currentEnrichmentTable())
    return(nrow(D) > 0)
  })
  outputOptions(output, "hasEnrichmentTable", suspendWhenHidden=FALSE)
  
  output$hasPredictions <- reactive({
    req(currentPredictions())
    return(TRUE)
  })
  outputOptions(output, "hasPredictions", suspendWhenHidden=FALSE)
  
  observeEvent(input$submitButton, {
    if(!input$useExampleData && !isTruthy(input$file)) {
      showNotification("Please select an expression data file and wait for it upload before submitting", type="error")
      return()
    }
    if(input$ntrees < MIN_NUM_TREES || input$ntrees > MAX_NUM_TREES) {
      showNotification(paste0("Number of trees must be >= ", MIN_NUM_TREES, " and <= ", MAX_NUM_TREES, "."), type="error")
      return()
    }
    if(!input$useExampleData && input$depvar == "") {
      showNotification("Missing dependent variable name.", type="error")
      return()
    }
    
    withProgress(value=0, message="Training grand forest", {
      # Read data file
      setProgress(value=0, detail="Reading data file")
      if(input$useExampleData) {
        D <- readRDS(EXAMPLE_DATA_PATH)
        depvar <- EXAMPLE_DATA_DEPVAR
        modelType <- EXAMPLE_DATA_MODELTYPE
      } else {
        D <- fread(input$file$datapath, header=TRUE, sep=",")
        depvar <- input$depvar
        modelType <- input$modelType
      }
      
      if((!depvar %in% colnames(D))) {
        showNotification("Dependent variable name does not match any column name in data file.", type="error")
        return()
      }
      
      # Scale and center data
      setProgress(value=0.1, detail="Normalizing data")
      depvar_col <- which(colnames(D) == depvar)
      D <- tryCatch({
        data.table(D[[depvar]], scale(D[,-depvar_col,with=FALSE], center=TRUE, scale=TRUE))
      }, error = function(e) {
        showNotification("Normalization failed. Not all columns are numeric.", type="error")
        req(FALSE)
      })
      colnames(D)[1] <- depvar
      
      # convert dependent variable to correct type
      if(modelType == "classification" || modelType == "probability") {
        D[[depvar]] <- as.factor(D[[depvar]])
      } else {
        D[[depvar]] <- as.numeric(D[[depvar]])
      }
      
      # Read graph file
      setProgress(value=0.2, detail="Preparing network")
      path <- get_network_file(input$graph)
      edges <- fread(path, header=FALSE, sep="\t", colClasses=rep("character", 2))
      colnames(edges) <- c("from","to")
      currentEdges(edges)
      
      # Train grand forest model
      setProgress(value=0.3, detail="Training model")
      fit <- grandforest(
        data=D, graph_data=edges,
        dependent.variable.name=depvar,
        probability=(modelType=="probability"),
        num.trees=input$ntrees,
        importance="impurity"
      )
      
      setProgress(value=0.9, detail="Finishing up")
      currentData(D)
      currentModel(list(fit=fit, depvar=depvar, type=modelType))
      currentEvaluation(NULL)
      currentPredictions(NULL)
    })
  })
    
  observeEvent(input$predictButton, {
    if(!isTruthy(input$predictFile)) {
      showNotification("Please upload a data set for prediction.", type="error")
      return()
    }
    
    withProgress(value=0, message="Predicting", {
      setProgress(value=0.1, detail="Preparing data")
      req(input$predictFile)
      req(currentModel())
      
      D <- fread(input$predictFile$datapath, header=TRUE, sep=",")
      fit <- currentModel()$fit
      
      setProgress(value=0.4, detail="Computing predictions")
      preds <- predict(fit, data=D)
      
      setProgress(value=0.9, detail="Preparing output")
      currentPredictions(preds$predictions)
    })
  })
  
  featureTable <- reactive({
    req(currentModel())
    imp <- importance(currentModel()$fit)
    
    genes <- names(imp)
    symbols <- mapIds(org.Hs.eg.db, genes, "SYMBOL", "ENTREZID")
    
    D <- data.table(gene=genes, name=symbols, importance=imp)
    D[order(D$importance, decreasing=TRUE),]
  })

  output$summary <- renderUI({
    model <- req(currentModel())
    fit <- model$fit
    depvar <- model$depvar
    tag("dl", list(class="dl-horizontal",
      tag("dt", "Model type"), tag("dd", fit$forest$treetype),
      tag("dt", "Dep. var. name"), tag("dd", depvar),
      tag("dt", "No. samples"), tag("dd", fit$num.samples),
      tag("dt", "No. features"), tag("dd", fit$num.independent.variables),
      tag("dt", "No. trees"), tag("dd", fit$num.trees),
      tag("dt", "OOB prediction error"), tag("dd", fit$prediction.error)
    ))
  })

  output$featureGraph <- renderVisNetwork({
    model <- req(currentModel())
    imp <- importance(model$fit)
    features <- head(featureTable(), n=input$nfeatures)
    
    D <- currentData()
    edges <- currentEdges()
    depvar <- model$depvar
    groups <- D[[depvar]]
    labels <- if(input$featureGraphGeneSymbols) features$name else features$gene

    feature_graph(D, edges, features$gene, labels, groups)
  })
  
  output$featureTable <- renderDataTable({
    head(featureTable(), n=input$nfeatures)
  }, options=list(
    scrollX = TRUE,
    searching = FALSE,
    pageLength = 10
  ))

  featureHeatmapPlot <- reactive({
    depvar <- currentModel()$depvar
    features <- head(featureTable(), n=input$nfeatures)
    D <- currentData()
    groups <- as.factor(D[[depvar]])
    D <- D[,features$gene,with=FALSE]
    colnames(D) <- paste0(features$gene, " (", features$name, ")")
    
    col.ramp <- colorRamp2(c(-2, 0, 2), c("magenta", "black", "green"))
    Heatmap(D, name="expression", split=groups, col=col.ramp)
  })
  
  output$featureHeatmap <- renderPlot({
    featureHeatmapPlot()
  })
  
  observeEvent(input$enrichmentButton, {
    req(featureTable())
    
    withProgress(message="Performing gene set enrichment", {
      setProgress(value=0.1, detail="Preparing data")
      features <- head(featureTable(), input$nfeatures)
      genes <- as.character(features$gene)
      depvar <- currentModel()$depvar
      D <- currentData()
      universe <- setdiff(colnames(D), depvar)
      
      setProgress(value=0.1, detail="Computing enrichment")
      out <- gene_set_enrichment(genes, universe, input$enrichmentType, input$enrichmentPvalueCutoff, input$enrichmentQvalueCutoff)
      
      setProgress(value=0.1, detail="Finishing up")
      currentEnrichmentTable(out)
    })
  })
  
  output$enrichmentTable <- renderDataTable({
    D <- req(as.data.frame(currentEnrichmentTable()))
    D$ID <- gene_set_enrichment_get_links(D$ID, isolate(input$enrichmentType))
    return(D)
  }, options = list(
    pageLength = 10,
    scrollX = TRUE
  ), escape = FALSE)
  
  output$enrichmentPlot <- renderPlot({
    D <- req(currentEnrichmentTable())
    DOSE::dotplot(D, showCategory=20)
  })
  
  output$predictionsTable <- renderDataTable({
    req(currentPredictions())
    preds <- currentPredictions()
    data.frame(sample=1:length(preds), prediction=preds)
  }, options = list(
    searching = FALSE,
    pageLength = 10
  ))
  
  observeEvent(input$evaluationButton, {
    if(input$cvTrees < MIN_NUM_TREES || input$cvTrees > MAX_NUM_TREES) {
      showNotification(paste0("Number of trees must be >= ", MIN_NUM_TREES, " and <= ", MAX_NUM_TREES, "."), type="error")
      return()
    }
    D <- currentData()
    edges <- currentEdges()
    depvar <- currentModel()$depvar
    modelType <- currentModel()$type
    
    out <- evaluation_cv(D, depvar, modelType, edges, input$cvTrees, input$cvFolds, input$cvRepetitions)
    currentEvaluation(out)
  })
  
  output$evaluationPerformance <- renderPlot({
    req(currentEvaluation())
    D <- currentEvaluation()$performance
    dummy <- data.frame(x = rep(1,4), y = c(0,1,-1,1), repetition=1, measure = c(rep("F1-score", 2), rep("Matthew's correlation coefficient", 2)))
    
    ggplot(aes(x=fold, y=value, color=as.factor(repetition)), data=D) +
      geom_line() +
      geom_point() +
      geom_blank(aes(x=x, y=y), data=dummy) +
      scale_x_continuous(breaks=1:max(D$fold)) +
      facet_wrap(~measure, ncol=1) +
      theme_minimal() +
      theme(
        text = element_text(size=18),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_text(size=14),
        axis.title.y = element_blank(),
        legend.position = "none"
      )
  })
  
  output$evaluationStability <- renderPlot({
    req(currentEvaluation())
    D <- currentEvaluation()$stability
    
    ggplot(aes(x=size, y=value), data=D) +
      geom_boxplot() +
      theme_minimal() +
      facet_wrap(~measure, ncol=1) +
      scale_y_continuous(limits=c(0,1)) +
      theme(
        text = element_text(size=18),
        axis.title.x = element_text(size=14),
        axis.title.y = element_blank()
      ) +
      labs(x="Top % features")
  })
  
  output$dlFeatureGraph <- downloadHandler(
    filename = function() {
      paste0("network_top", input$nfeatures, ".sif")
    },
    content = function(file) {
      features <- head(featureTable(), n=input$nfeatures)
      edges <- currentEdges()
      edges <- subset(edges, from %in% features$gene & to %in% features$gene)
      D <- data.frame(from=edges$from, type=".", to=edges$to)
      if(input$featureGraphGeneSymbols) {
        D$from <- mapIds(org.Hs.eg.db, as.character(D$from), "SYMBOL", "ENTREZID")
        D$to <- mapIds(org.Hs.eg.db, as.character(D$to), "SYMBOL", "ENTREZID")
      }
      write.table(D, file, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
    }
  )
  
  output$dlFeatureHeatmap <- downloadHandler(
    filename = "heatmap.pdf",
    content = function(file) {
      pdf(file=file, width=10, height=10)
      print(featureHeatmapPlot())
      dev.off()
    }
  )

  output$dlFeatureTable <- downloadHandler(
    filename = function() {
      paste0("genes_top", input$nfeatures, ".csv")
    },
    content = function(file) {
      write.csv(head(featureTable(), input$nfeatures), file, row.names=FALSE)
    }
  )
  
  output$dlFeatureTableFull <- downloadHandler(
    filename = "genes.csv",
    content = function(file) {
      write.csv(featureTable(), file, row.names=FALSE)
    }
  )
  
  output$dlPredictionsTable <- downloadHandler(
    filename = "predictions.csv",
    content = function(file) {
      write.csv(currentPredictions(), file, row.names=FALSE)
    }
  )
  
  output$dlEnrichmentTable <- downloadHandler(
    filename = "enrichment.csv",
    content = function(file) {
      write.csv(as.data.frame(currentEnrichmentTable()), file, row.names=FALSE)
    }
  )
  
  hide(id="loading-content", anim=TRUE, animType="fade")
})
