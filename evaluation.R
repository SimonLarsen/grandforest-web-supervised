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
  if(length(levels(labels)) == 2) {
    data.frame(
      `F1 score` = caret::confusionMatrix(preds, labels)$byClass["F1"],
      `Matthews correlation coefficient` = mccr::mccr(as.numeric(preds)-1, as.numeric(labels)-1),
      check.names = FALSE
    )
  } else {
    cm <- caret::confusionMatrix(preds, labels)$byClass
    mprec <- mean(cm[,"Precision"])
    mrec <- mean(cm[,"Recall"])
    
    data.frame(
      `F1 score macro avg.` = mean(cm[,"F1"]),
      `F1 score micro avg.` = 2 * (mprec * mrec) / (mprec + mrec),
      check.names = FALSE
    )
  }
}

evaluation_regression <- function(preds, labels) {
  ym <- mean(labels)
  SSres <- sum((labels - preds)^2)
  SStot <- sum((labels - ym)^2)
  
  R2 <- 1 - SSres/SStot
  RMSE <- sqrt(mean((preds - labels)^2))
  data.frame(
    `R squared` = R2,
    `Root-mean-squared error` = RMSE,
    check.names = FALSE
  )
}

evaluation_probability <- function(preds, labels) {
  preds_f <- colnames(preds)[apply(preds, 1, which.max)]
  evaluation_classification(preds_f, labels)
}

evaluation_cv <- function(D, depvar, modelType, edges, ntrees, cvFolds, cvRepetitions) {
  if(modelType == "survival") {
    alert("Model evaluation not supported for survival data.")
    return()
  }
  withProgress(value=0, max=cvFolds*cvRepetitions, message="Performing evaluation", {
    scores <- data.frame()
    stability <- data.frame()
    
    for(rep in 1:cvRepetitions) {
      fits <- c()
      rep.scores <- data.frame()
      folds <- caret::createFolds(D[[depvar]], k=cvFolds, returnTrain=TRUE, list=TRUE)
      
      for(i in 1:cvFolds) {
        setProgress(value=(rep-1)*cvFolds+(i-1), detail=paste("Repetition ", rep, " Fold", i))
        fits[[i]] <- grandforest(data=D[folds[[i]],], graph_data=edges, dependent.variable.name=depvar, num.trees=ntrees, importance="impurity")
        preds <- predict(fits[[i]], data=D[-folds[[i]],])$predictions
        labels <- D[-folds[[i]],][[depvar]]
        
        if(modelType == "classification") {
          eval <- evaluation_classification(preds, labels)
        } else if(modelType == "regression") {
          eval <- evaluation_regression(preds, labels)
        } else if(modelType == "probability") {
          eval <- evaluation_probability(preds, labels)
        }
        
        rep.scores <- rbind(rep.scores, eval)
      }
      
      rep.scores$fold <- 1:cvFolds
      rep.scores$repetition <- rep
      scores <- rbind(scores, rep.scores)
      
      rep.stability <- evaluation_stability(fits)
      rep.stability$repetition <- rep
      
      stability <- rbind(stability, rep.stability)
    }
  })
  
  list(
    performance = reshape2::melt(scores, id.vars=c("repetition","fold"), variable.name="measure", value.name="value"),
    stability = stability
  )
}
