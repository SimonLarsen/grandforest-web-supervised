# This file is part of caret.
# Obtained from: https://github.com/topepo/caret/blob/master/pkg/caret/R/confusionMatrix.R.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

confusionMatrix <-
  function(data, ...){
    UseMethod("confusionMatrix")
  }

confusionMatrix.default <- function(data, reference,
                                    positive = NULL,
                                    dnn = c("Prediction", "Reference"),
                                    prevalence = NULL,
                                    mode = "sens_spec",
                                    ...) {
  if(!(mode %in% c("sens_spec", "prec_recall", "everything")))
    stop("`mode` should be either 'sens_spec', 'prec_recall', or 'everything'")
  if(!is.factor(data) | !is.factor(reference)) {
    stop("`data` and `reference` should be factors with the same levels.", call. = FALSE)
  }
  if(!is.character(positive) & !is.null(positive)) stop("positive argument must be character")
  
  if(length(levels(data)) > length(levels(reference)))
    stop("the data cannot have more levels than the reference")
  
  if(!any(levels(data) %in% levels(reference))){
    stop("The data must contain some levels that overlap the reference.")
  }
  
  if(!all(levels(data) %in% levels(reference))){
    badLevel <- levels(data)[!levels(data) %in% levels(reference)]
    if(sum(table(data)[badLevel]) > 0){
      stop("The data contain levels not found in the data.")
    } else{
      warning("The data contains levels not found in the data, but they are empty and will be dropped.")
      data <- factor(as.character(data))
    }
  }
  
  if(any(levels(reference) != levels(data))) {
    warning("Levels are not in the same order for reference and data. Refactoring data to match.")
    data <- as.character(data)
    data <- factor(data, levels = levels(reference))
  }
  classLevels <- levels(data)
  numLevels <- length(classLevels)
  if(numLevels < 2)
    stop("there must be at least 2 factors levels in the data")
  
  if(numLevels == 2 & is.null(positive))  positive <- levels(reference)[1]
  
  classTable <- table(data, reference, dnn = dnn, ...)
  
  getFromNamespace("confusionMatrix.table", "caret")(classTable, positive, prevalence = prevalence, mode = mode)
}

confusionMatrix.table <- function(data, positive = NULL,
                                  prevalence = NULL, mode = "sens_spec", ...){
  requireNamespaceQuietStop("e1071")
  if(!(mode %in% c("sens_spec", "prec_recall", "everything")))
    stop("`mode` should be either 'sens_spec', 'prec_recall', or 'everything'")
  if(length(dim(data)) != 2) stop("the table must have two dimensions")
  if(!all.equal(nrow(data), ncol(data))) stop("the table must nrow = ncol")
  if(!all.equal(rownames(data), colnames(data))) stop("the table must the same classes in the same order")
  if(!is.character(positive) & !is.null(positive)) stop("positive argument must be character")
  
  classLevels <- rownames(data)
  numLevels <- length(classLevels)
  if(numLevels < 2)
    stop("there must be at least 2 factors levels in the data")
  
  if(numLevels == 2 & is.null(positive))  positive <- rownames(data)[1]
  
  
  if(numLevels == 2 & !is.null(prevalence) && length(prevalence) != 1)
    stop("with two levels, one prevalence probability must be specified")
  
  if(numLevels > 2 & !is.null(prevalence) && length(prevalence) != numLevels)
    stop("the number of prevalence probability must be the same as the number of levels")
  
  if(numLevels > 2 & !is.null(prevalence) && is.null(names(prevalence)))
    stop("with >2 classes, the prevalence vector must have names")
  
  propCI <- function(x) {
    res <- try(binom.test(sum(diag(x)), sum(x))$conf.int, silent = TRUE)
    if(inherits(res, "try-error"))
      res <- rep(NA, 2)
    res
  }
  
  propTest <- function(x){
    res <- try(
      binom.test(sum(diag(x)),
                 sum(x),
                 p = max(apply(x, 2, sum)/sum(x)),
                 alternative = "greater"),
      silent = TRUE)
    res <- if(inherits(res, "try-error"))
      c("null.value.probability of success" = NA, p.value = NA)
    else
      res <- unlist(res[c("null.value", "p.value")])
    res
  }
  
  overall <- c(unlist(e1071::classAgreement(data))[c("diag", "kappa")],
               propCI(data),
               propTest(data),
               mcnemar.test(data)$p.value)
  
  names(overall) <- c("Accuracy", "Kappa", "AccuracyLower", "AccuracyUpper", "AccuracyNull", "AccuracyPValue", "McnemarPValue")
  
  if(numLevels == 2) {
    if(is.null(prevalence)) prevalence <- sum(data[, positive])/sum(data)
    negative <- classLevels[!(classLevels %in% positive)]
    tableStats <- c(sensitivity.table(data, positive),
                    specificity.table(data, negative),
                    posPredValue.table(data, positive, prevalence = prevalence),
                    negPredValue.table(data, negative, prevalence = prevalence),
                    precision.table(data, relevant = positive),
                    recall.table(data, relevant = positive),
                    F_meas.table(data, relevant = positive),
                    prevalence,
                    sum(data[positive, positive])/sum(data),
                    sum(data[positive, ])/sum(data))
    names(tableStats) <- c("Sensitivity", "Specificity",
                           "Pos Pred Value", "Neg Pred Value",
                           "Precision", "Recall", "F1",
                           "Prevalence", "Detection Rate",
                           "Detection Prevalence")
    tableStats["Balanced Accuracy"] <- (tableStats["Sensitivity"]+tableStats["Specificity"])/2
    
  } else {
    
    tableStats <- matrix(NA, nrow = length(classLevels), ncol = 11)
    
    for(i in seq(along = classLevels)) {
      pos <- classLevels[i]
      neg <- classLevels[!(classLevels %in% classLevels[i])]
      prev <- if(is.null(prevalence)) sum(data[, pos])/sum(data) else prevalence[pos]
      tableStats[i,] <- c(sensitivity.table(data, pos),
                          specificity.table(data, neg),
                          posPredValue.table(data, pos, prevalence = prev),
                          negPredValue.table(data, neg, prevalence = prev),
                          precision.table(data, relevant = pos),
                          recall.table(data, relevant = pos),
                          F_meas.table(data, relevant = pos),
                          prev,
                          sum(data[pos, pos])/sum(data),
                          sum(data[pos, ])/sum(data), NA)
      tableStats[i,11] <- (tableStats[i,1] + tableStats[i,2])/2
    }
    rownames(tableStats) <- paste("Class:", classLevels)
    colnames(tableStats) <- c("Sensitivity", "Specificity",
                              "Pos Pred Value", "Neg Pred Value",
                              "Precision", "Recall", "F1",
                              "Prevalence", "Detection Rate",
                              "Detection Prevalence", "Balanced Accuracy")
  }
  
  structure(
    list(positive = positive,
         table = data,
         overall = overall,
         byClass = tableStats,
         mode = mode,
         dots = list(...)),
    class = "confusionMatrix")
}