install.packages(c('devtools','shiny', 'shinythemes', 'visNetwork', 'data.table', 'igraph', 'magrittr', 'circlize', 'ggplot2', 'e1071', 'caret', 'mccr'))

source('https://bioconductor.org/biocLite.R')
biocLite('qvalue', ask=FALSE)
biocLite('ComplexHeatmap', ask=FALSE)
biocLite('org.Hs.eg.db', ask=FALSE)
biocLite('DOSE', ask=FALSE)
biocLite('clusterProfiler', ask=FALSE)

devtools::install_github('AnalytixWare/ShinySky')
devtools::install_github('SimonLarsen/grandforest', auth_token='edb93ddf877c24e3e71a66fdb96092c2e96bed82')
