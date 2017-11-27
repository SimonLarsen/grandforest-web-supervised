library(caret)
library(mccr)

evaluation_stability <- function(fits, sizes=c(1,5,10,20,40)) {
  features <- lapply(fits, function(fit) {
    imp <- importance(fit)
    names(imp[order(imp, decreasing=TRUE)])
  })
  
  total.features <- length(features[[1]])
  do.call(rbind, lapply(sizes, function(size) {
    nfeatures <- round(total.features * (size/100))
    jaccard <- sapply(combn(features, 2, simplify=FALSE), function(pair) {
      a <- head(pair[[1]], nfeatures)
      b <- head(pair[[2]], nfeatures)
      length(intersect(a,b)) / length(union(a,b))
    })
    
    data.frame(measure="Jaccard index", value=jaccard, size=paste0(size,"%"))
  }))
}

evaluation_classification <- function(preds, labels) {
  list(
    f.scores = confusionMatrix(preds, labels)$byClass["F1"],
    mat.scores = mccr(as.numeric(preds)-1, as.numeric(labels)-1)
  )
}

evaluation_regression <- function(preds, labels) {
  
}

evaluation_probability <- function(preds, labels) {
  
}

evaluation_cv <- function(D, depvar, modelType, edges, ntrees, cvFolds, cvRepetitions) {
  withProgress(value=0, max=cvFolds*cvRepetitions, message="Performing evaluation", {
    f.scores <- list()
    mat.scores <- list()
    stability <- list()
    
    for(rep in 1:cvRepetitions) {
      fits <- c()
      rep.f.scores <- c()
      rep.mat.scores <- c()
      folds <- createFolds(D[[depvar]], k=cvFolds, returnTrain=TRUE, list=TRUE)
      
      for(i in 1:cvFolds) {
        setProgress(value=(rep-1)*cvFolds+(i-1), detail=paste("Repetition ", rep, " Fold", i))
        fits[[i]] <- grandforest(data=D[folds[[i]],], graph_data=edges, dependent.variable.name=depvar, num.trees=ntrees, importance="impurity")
        preds <- predict(fits[[i]], data=D[-folds[[i]],])$predictions
        labels <- D[-folds[[i]],][[depvar]]
        
        eval <- evaluation_classification(preds, labels)
        
        rep.f.scores[[i]] <- eval$f.scores
        rep.mat.scores[[i]] <- eval$mat.scores
      }
      
      f.scores[[rep]] <- rep.f.scores
      mat.scores[[rep]] <- rep.mat.scores
      
      rep.stability <- evaluation_stability(fits)
      rep.stability$repetition <- rep
      
      stability[[rep]] <- rep.stability
    }
  })
  
  f.scores <- do.call(rbind, lapply(1:length(f.scores), function(i) {
    data.frame(fold=1:length(f.scores[[i]]), value=f.scores[[i]], repetition=i, measure="F1-score")
  }))
  mat.scores <- do.call(rbind, lapply(1:length(mat.scores), function(i) {
    data.frame(fold=1:length(mat.scores[[i]]), value=mat.scores[[i]], repetition=i, measure="Matthew's correlation coefficient")
  }))
  
  return(list(
    performance = rbind(f.scores, mat.scores),
    stability = do.call(rbind, stability)
  ))
}