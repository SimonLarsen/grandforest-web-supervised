library(shiny)
library(shinyjs)
library(data.table)
library(magrittr)
library(grandforest)

source("grandforest-web-common/read_data.R")
source("grandforest-web-common/get_network.R")
source("grandforest-web-common/enrichment.R")
source("grandforest-web-common/targets.R")
source("grandforest-web-common/feature_graph.R")
source("grandforest-web-common/mapping.R")
source("grandforest-web-common/make_links.R")
source("evaluation.R")

shinyServer(function(input, output, session) {
  currentEdges <- reactiveVal()
  currentData <- reactiveVal()
  currentSpecies <- reactiveVal()
  currentModel <- reactiveVal()
  currentEnrichmentTable <- reactiveVal()
  currentTargetsTable <- reactiveVal()
  currentEvaluation <- reactiveVal()
  currentPredictions <- reactiveVal()

  output$hasModel <- reactive({
    req(currentModel())
    return(TRUE)
  })
  outputOptions(output, "hasModel", suspendWhenHidden=FALSE)

  output$hasEnrichmentTable <- reactive({
    req(currentEnrichmentTable())
    return(TRUE)
  })
  outputOptions(output, "hasEnrichmentTable", suspendWhenHidden=FALSE)

  output$hasTargetsTable <- reactive({
    req(currentTargetsTable())
    return(TRUE)
  })
  outputOptions(output, "hasTargetsTable", suspendWhenHidden=FALSE)

  output$hasPredictions <- reactive({
    req(currentPredictions())
    return(TRUE)
  })
  outputOptions(output, "hasPredictions", suspendWhenHidden=FALSE)

  output$currentSpecies <- reactive({
    currentSpecies()
  })
  outputOptions(output, "currentSpecies", suspendWhenHidden=FALSE)

  output$graphSelect <- renderUI({
    selectInput("graph", "Network", get_network_options(input$species))
  })

  output$enrichmentTypeSelect <- renderUI({
    selectInput("enrichmentType", "Enrichment type", gene_set_enrichment_types(currentSpecies()))
  })

  observeEvent(input$submitButton, {
    if(!input$useExampleData && !isTruthy(input$file)) {
      alert("Please select an expression data file and wait for it upload before submitting")
      return()
    }
    if(input$useExampleData && input$species != "human") {
      alert("Please set species to \"Homo sapiens\" when using example data.")
      return()
    }
    if(input$ntrees < MIN_NUM_TREES || input$ntrees > MAX_NUM_TREES) {
      alert(paste0("Number of trees must be >= ", MIN_NUM_TREES, " and <= ", MAX_NUM_TREES, "."))
      return()
    }
    if(!input$useExampleData && input$depvar == "") {
      alert("Missing dependent variable name.")
      return()
    }
    if(input$modelType == "survival" && input$statusvar == "") {
      alert("Missing status variable for survival model.")
      return()
    }

    withProgress(value=0, message="Training grand forest", {
      # Read data file
      setProgress(value=0.1, detail="Reading data file")
      if(input$useExampleData) {
        D <- readRDS(EXAMPLE_DATA_PATH)
        depvar <- EXAMPLE_DATA_DEPVAR
        modelType <- EXAMPLE_DATA_MODELTYPE
        statusvar <- NULL
      } else {
        D <- tryCatch(
          read_expression_file(input$file$datapath),
          error = function(e) { alert(e$message); return(NULL) }
        )
        if(is.null(D)) return()
        depvar <- input$depvar
        modelType <- input$modelType
        statusvar <- if(modelType == "survival") input$statusvar else NULL
      }

      if(!(depvar %in% colnames(D))) {
        alert("Dependent variable name does not match any column in data file.")
        return()
      }
      if(modelType == "survival" && !(statusvar %in% colnames(D))) {
        alert("Status variable name does not match any column in data file.")
        return()
      }

      # Read graph file
      setProgress(value=0.2, detail="Preparing network")
      edges <- readRDS(get_network_file(input$graph))

      all_nodes <- unique(c(edges$from, edges$to))
      found_genes <- intersect(colnames(D), all_nodes)
      missing_genes <- setdiff(setdiff(colnames(D), c(depvar,statusvar)), all_nodes)

      if(length(found_genes) == 0) {
        alert("No expression was found for any genes in the network. Make sure you chose the right species, and that gene IDs are NCBI Entrez ID.")
        return()
      }

      # Extract valid features
      D <- D[,c(depvar,statusvar,found_genes),with=FALSE]

      # convert dependent variable to correct type
      if(modelType == "classification" || modelType == "probability") {
        D[[depvar]] <- as.factor(D[[depvar]])
        if(length(levels(D[[depvar]])) > MAX_NUM_CLASSES) {
          alert(sprintf("Dependent variable contains too many different classes. Only up to %d allowed.", MAX_NUM_CLASSES))
          return()
        }
      } else if(modelType == "survival") {
        D[[depvar]] <- as.numeric(D[[depvar]])
        D[[statusvar]] <- as.numeric(D[[statusvar]])
        if(any(D[[statusvar]] < 0 | D[[statusvar]] > 1, na.rm=TRUE)) {
          alert("Survival status variable could not be coerced to 0 and 1.")
          return()
        }
      } else if(modelType == "regression") {
        D[[depvar]] <- as.numeric(D[[depvar]])
      } else {
        alert("Encountered unknown model type."); return()
      }

      # Train grand forest model
      setProgress(value=0.3, detail="Training model")
      fit <- grandforest(
        data=D, graph_data=edges,
        dependent.variable.name=depvar,
        status.variable.name=statusvar,
        probability=(modelType=="probability"),
        num.trees=input$ntrees,
        importance="impurity"
      )

      setProgress(value=0.9, detail="Finishing up")
      currentEdges(edges)
      currentData(D)
      currentSpecies(input$species)
      currentModel(list(
        fit=fit,
        type=modelType,
        depvar=depvar,
        statusvar=statusvar,
        found_genes=found_genes,
        missing_genes=missing_genes
      ))
      currentEnrichmentTable(NULL)
      currentTargetsTable(NULL)
      currentEvaluation(NULL)
      currentPredictions(NULL)
    })
  })

  observeEvent(input$predictButton, {
    if(!isTruthy(input$predictFile)) {
      alert("Please upload a data set for prediction.")
      return()
    }
    if(currentModel()$type == "survival") {
      alert("Prediction not supported for survival data.")
      return()
    }

    withProgress(value=0, message="Predicting", {
      setProgress(value=0.1, detail="Preparing data")
      req(input$predictFile)
      req(currentModel())

      D <- tryCatch(
        read_expression_file(input$predictFile$datapath),
        error = function(e) { alert(e$message); return(NULL) }
      )
      if(is.null(D)) return()

      setProgress(value=0.4, detail="Computing predictions")
      fit <- currentModel()$fit
      preds <- predict(fit, data=D)

      setProgress(value=0.9, detail="Preparing output")
      currentPredictions(preds$predictions)
    })
  })

  featureTable <- reactive({
    req(currentModel())
    imp <- importance(currentModel()$fit)

    genes <- names(imp)
    symbols <- map_ids_fallback(genes, "SYMBOL", "ENTREZID", currentSpecies())
    D <- data.table(gene=genes, name=symbols, importance=imp)
    D[order(D$importance, decreasing=TRUE),]
  })

  output$summary <- renderUI({
    model <- req(currentModel())
    fit <- model$fit
    depvar <- model$depvar
    found_pct <- length(model$found_genes) / (length(model$found_genes)+length(model$missing_genes)) * 100
    tag("dl", list(class="dl-horizontal",
      tag("dt", "Model type"), tag("dd", fit$forest$treetype),
      tag("dt", "Dep. var. name"), tag("dd", depvar),
      tag("dt", "Genes in network"), tag("dd", list(
        sprintf("%.2f %%", found_pct),
        downloadLink("dlMissingGenes", sprintf("(%d missing)", length(model$missing_genes)))
      )),
      tag("dt", "No. samples"), tag("dd", fit$num.samples),
      tag("dt", "No. features"), tag("dd", fit$num.independent.variables),
      tag("dt", "No. trees"), tag("dd", fit$num.trees),
      tag("dt", "OOB prediction error"), tag("dd", sprintf("%.3f", fit$prediction.error))
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
    D <- head(featureTable(), n=input$nfeatures)
    D$gene <- make_links(D$gene, "ncbi_gene")
    return(D)
  }, options=list(pageLength=10, searching=FALSE, scrollX = TRUE), escape=FALSE)

  featureHeatmapPlot <- reactive({
    depvar <- currentModel()$depvar
    features <- head(featureTable(), n=input$nfeatures)
    D <- currentData()
    anno <- NULL
    gaps_row <- NULL
    
    if(currentModel()$type == "classification" || currentModel()$type == "probability") {
      anno <- data.frame(as.factor(D[[depvar]]))
      colnames(anno) <- depvar
      row_order <- order(anno[[depvar]])
    
      gaps_row <- head(cumsum(table(anno[[depvar]])), -1)
    }
    else if(currentModel()$type == "regression") {
      anno <- data.frame(as.numeric(D[[depvar]]))
      colnames(anno) <- depvar
      row_order <- order(anno[[depvar]])
    }
    else if(currentModel()$type == "survival") {
      anno <- data.frame(
        time = as.numeric(D[[depvar]]),
        status = as.factor(D[[currentModel()$statusvar]])
      )
      row_order <- order(anno$time)
    }
    rownames(anno) <- paste0("p", seq_len(nrow(anno)))
    
    D <- scale(D[,features$gene,with=FALSE])
    rownames(D) <- paste0("p", seq_len(nrow(D)))
    if(input$featureHeatmapGeneSymbol) {
      colnames(D) <- features$name
    }
    
    D[D >  2] <-  2
    D[D < -2] <- -2
    
    pheatmap::pheatmap(
      D[row_order,],
      color = colorRampPalette(c("magenta","black","green"))(100),
      annotation_row=anno,
      cluster_rows=FALSE,
      show_rownames=FALSE,
      gaps_row = gaps_row
    )
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

      setProgress(value=0.2, detail="Computing enrichment")
      out <- gene_set_enrichment(genes, currentSpecies(), universe, input$enrichmentType, input$enrichmentPvalueCutoff, input$enrichmentQvalueCutoff)

      setProgress(value=0.9, detail="Finishing up")
      currentEnrichmentTable(out)
    })
  })

  output$enrichmentTable <- renderDataTable({
    D <- as.data.frame(req(currentEnrichmentTable()))
    validate(need(nrow(D) > 0, "No significantly enriched gene sets found."))
    D <- get_gene_set_enrichment_links(D, isolate(input$enrichmentType))
    return(D)
  }, options = list(pageLength=10, scrollX=TRUE), escape=FALSE)

  output$enrichmentPlot <- renderPlot({
    D <- req(currentEnrichmentTable())
    validate(need(nrow(D) > 0, "No significantly enriched gene sets found."))
    DOSE::dotplot(D, showCategory=20)
  })

  observeEvent(input$targetsButton, {
    req(featureTable())
    features <- head(featureTable(), input$nfeatures)
    genes <- as.character(features$gene)
    out <- get_gene_targets(genes, input$targetsType)
    currentTargetsTable(out)
  })

  output$targetsTable <- renderDataTable({
    D <- req(currentTargetsTable())
    validate(need(nrow(D) > 0, "No targets found."))
    get_gene_targets_table(D, isolate(input$targetsType))
  }, options = list(pageLength=10, scrollX=TRUE), escape=FALSE)

  output$targetsNetwork <- renderVisNetwork({
    D <- req(currentTargetsTable())
    validate(need(nrow(D) > 0, "No targets found."))
    get_gene_targets_network(D, isolate(input$targetsType), input$targetsNetworkSymbols)
  })

  output$predictionsTable <- renderDataTable({
    req(currentPredictions())
    preds <- currentPredictions()
    data.frame(sample=1:length(preds), prediction=preds)
  }, options = list(searching = FALSE, pageLength = 10))

  observeEvent(input$evaluationButton, {
    if(input$cvTrees < MIN_NUM_TREES || input$cvTrees > MAX_NUM_TREES) {
      alert(paste0("Number of trees must be >= ", MIN_NUM_TREES, " and <= ", MAX_NUM_TREES, "."))
      return()
    }
    if(currentModel()$type == "survival") {
      alert("Model evaluation not supported for survival data.")
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
    dummy <- subset(rbind(
      data.frame(x=1, repetition=1, y=c(0,1),  measure="F1 score"),
      data.frame(x=1, repetition=1, y=c(0,1),  measure="F1 score micro avg."),
      data.frame(x=1, repetition=1, y=c(0,1),  measure="F1 score macro avg."),
      data.frame(x=1, repetition=1, y=c(-1,1), measure="Matthews correlation coefficient"),
      data.frame(x=1, repetition=1, y=c(0,1),  measure="R squared"),
      data.frame(x=1, repetition=1, y=0,       measure="Root-mean-squared error")
    ), measure %in% levels(D$measure))

    library(ggplot2)
    ggplot(aes(x=fold, y=value, color=as.factor(repetition)), data=D) +
      geom_line() +
      geom_point() +
      geom_blank(aes(x=x, y=y), data=dummy) +
      scale_x_continuous(breaks=1:max(D$fold)) +
      facet_wrap(~measure, ncol=1, scales="free") +
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
    library(ggplot2)
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
        D$from <- map_ids_fallback(as.character(D$from), "SYMBOL", "ENTREZID", currentSpecies())
        D$to <- map_ids_fallback(as.character(D$to), "SYMBOL", "ENTREZID", currentSpecies())
      }
      write.table(D, file, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
    }
  )

  output$dlFeatureHeatmap <- downloadHandler(
    filename = "heatmap.pdf",
    content = function(file) {
      p <- featureHeatmapPlot()
      pdf(file=file, width=9, height=9)
      print(p)
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
    filename = function() { paste0("enrichment.", input$enrichmentType, ".csv") },
    content = function(file) { write.csv(as.data.frame(currentEnrichmentTable()), file, row.names=FALSE) }
  )

  output$dlTargetsTable <- downloadHandler(
    filename = function() { paste0("targets.", input$targetsType, ".csv") },
    content = function(file) { write.csv(currentTargetsTable(), file, row.names=FALSE) }
  )

  output$dlMissingGenes <- downloadHandler(
    filename = "missing_genes.txt",
    content = function(file) {
      write(currentModel()$missing_genes, file=file)
    }
  )

  hide(id="loading-content", anim=TRUE, animType="fade")
})
