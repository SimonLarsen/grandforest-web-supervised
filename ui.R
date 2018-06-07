library(shiny)
library(shinyjs)
library(visNetwork)
library(shinycssloaders)

source("grandforest-web-common/enrichment.R")
source("grandforest-web-common/targets.R")

tooltip_label <- function(text, tooltip) {
  HTML(sprintf('<span data-toggle="tooltip" data-placement="top" title="%s">%s <i class="fa fa-question-circle"></i></span>', tooltip, text))
}

shinyUI(tagList(
  tags$head(
    tags$link(rel="stylesheet", type="text/css", href="style.css"),
    tags$link(rel="stylesheet", type="text/css", href="loader.css"),
    tags$script(type="text/javascript", '/* load tooltips */ $(function () { $(\'[data-toggle="tooltip"]\').tooltip() })')
  ),
  useShinyjs(),
  div(id="loading-content", h2("Loading..."), div(class="loader", "Loading")),
  navbarPage("Grand Forest • Supervised", footer=column(width=12, hr(), p(paste0("Grand Forest • Supervised workflow • Version ", APP_VERSION))), inverse=TRUE,
    tabPanel(HTML("Analysis</a></li><li><a href=\"https://grandforest.compbio.sdu.dk/guide/supervised\" target=_blank>User guide</a></li><li><a href=\"https://grandforest.compbio.sdu.dk/#cite\" target=_blank>Cite"),
      sidebarLayout(
        sidebarPanel(width=3,
          tags$h3("Train model", class="sidebar-top-heading"),
          checkboxInput("useExampleData", "Use example data"),
          conditionalPanel("input.useExampleData == false",
            fileInput("file", tooltip_label("Expression data", "See `User guide` for a description of supported file formats.")),
            selectInput("modelType", "Model type", list(
              "Classification" = "classification",
              "Survival" = "survival",
              "Regression" = "regression",
              "Probability" = "probability"
            )),
            textInput("depvar", tooltip_label("Dependent variable name", "Name of column containing dependent variable. Must match exactly, including capitalization.")),
            conditionalPanel("input.modelType == 'survival'",
              textInput("statusvar", tooltip_label("Status variable name", "Name of column containing survival status variable. 0 = alive, 1 = death."))
            )
          ),
          numericInput("ntrees", tooltip_label("Number of decision trees", "How many decision trees to train. Larger number produces better results but takes longer to compute."), DEFAULT_NUM_TREES, min = MIN_NUM_TREES, max = MAX_NUM_TREES),
          selectInput("species", "Species", list(
            "Homo sapiens" = "human",
            "Mus musculus" = "mouse"
          )),
          uiOutput("graphSelect"),
          actionButton("submitButton", "Train grand forest", class = "btn-primary"),
          conditionalPanel("output.hasModel == true",
            h3("Model summary"),
            uiOutput("summary"),
            h3("Parameters"),
            sliderInput("nfeatures", tooltip_label("Module size", "Number of features to extract, ranked by importance."), min=MIN_NUM_FEATURES, max=MAX_NUM_FEATURES, value=DEFAULT_NUM_FEATURES, step=1, width = "400px")
          )
        ),
        mainPanel(width=9,
          conditionalPanel("output.hasModel == true",
            tabsetPanel(id="mainTabs", type="tabs",
              tabPanel("Analysis",
                tags$div(class="page-header", h2("Analysis")),
                fluidRow(
                  column(width=6,
                    h3("Feature subnetwork"),
                    wellPanel(
                      withSpinner(visNetworkOutput("featureGraph", height=450)),
                      fluidRow(
                        column(width=4, downloadButton("dlFeatureGraph", "Download network", class="btn-sm")),
                        column(width=4, checkboxInput("featureGraphGeneSymbols", "Show gene symbols", value=TRUE))
                      )
                    )
                  ),
                  column(width=6,
                    h3("Heatmap"),
                    wellPanel(
                      withSpinner(plotOutput("featureHeatmap", height=500)),
                      fluidRow(
                        column(width=4, downloadButton("dlFeatureHeatmap", "Download heatmap", class="btn-sm")),
                        column(width=4, checkboxInput("featureHeatmapGeneSymbol", "Show gene symbols", value=TRUE))
                      )
                    )
                  )
                ),
                h3("Genes"),
                div(class="body-tabs", tabsetPanel(type="tabs",
                  tabPanel("Feature table",
                    dataTableOutput("featureTable"),
                    downloadButton("dlFeatureTable", "Download table", class="btn-sm"),
                    downloadButton("dlFeatureTableFull", "Download full table", class="btn-sm")
                  ),
                  tabPanel("Gene set enrichment",
                    fluidRow(
                      column(width=4, uiOutput("enrichmentTypeSelect")),
                      column(width=4, numericInput("enrichmentPvalueCutoff", "p-value cutoff", value=0.05, min=0, max=1, step=0.01)),
                      column(width=4, numericInput("enrichmentQvalueCutoff", "q-value cutoff", value=0.2, min=0, max=1, step=0.01))
                    ),
                    actionButton("enrichmentButton", "Run enrichment analysis", class="btn-primary"),
                    conditionalPanel("output.hasEnrichmentTable == true",
                      hr(),
                      tabsetPanel(
                        tabPanel("Table",
                          dataTableOutput("enrichmentTable"),
                          downloadButton("dlEnrichmentTable", "Download table", class="btn-sm")
                        ),
                        tabPanel("Dot plot", withSpinner(plotOutput("enrichmentPlot")))
                      )
                    )
                  ),
                  tabPanel("Drug/miRNA targets",
                    conditionalPanel("output.currentSpecies == 'human'", 
                      selectInput("targetsType", "Database", gene_target_sources()),
                      actionButton("targetsButton", "Get gene targets", class="btn-primary"),
                      conditionalPanel("output.hasTargetsTable == true",
                        hr(),
                        tabsetPanel(type="tabs",
                          tabPanel("Table",
                            withSpinner(dataTableOutput("targetsTable")),
                            downloadButton("dlTargetsTable", "Download table", class="btn-sm")
                          ),
                          tabPanel("Network",
                            withSpinner(visNetworkOutput("targetsNetwork")),
                            checkboxInput("targetsNetworkSymbols", "Show gene symbols", value=TRUE)
                          )
                        )
                      )
                    ),
                    conditionalPanel("output.currentSpecies != 'human'",
                      div(class="alert alert-info", "Drug/miRNA target search only supported for Homo sapiens.")
                    )
                  )
                ))
              ),
              tabPanel("Evaluation",
                tags$div(class="page-header",
                  h2("Evaluation")
                ),
                h3("Parameters"),
                wellPanel(
                  fluidRow(
                    column(width=4, numericInput("cvFolds", "Number of cross-validation folds", min=2, max=10, value=5)),
                    column(width=4, numericInput("cvRepetitions", "Number of repetitions", min=1, max=5, value=1)),
                    column(width=4, numericInput("cvTrees", "Number of decision trees", min=MIN_NUM_TREES, max=MAX_NUM_TREES, value=DEFAULT_NUM_TREES))
                  ),
                  actionButton("evaluationButton", "Run evaluation", class="btn-primary")
                ),
                h3("Results"),
                fluidRow(
                  column(width=6,
                    h4("Predictive performance"),
                    wellPanel(withSpinner(plotOutput("evaluationPerformance")))
                  ),
                  column(width=6,
                    h4("Feature stability"),
                    wellPanel(withSpinner(plotOutput("evaluationStability")))
                  )
                )
              )
            )
          )
        )
      )
    )
  )
))
