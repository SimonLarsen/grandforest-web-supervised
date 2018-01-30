## Overview

## Expression data format

Expression data must be uploaded in CSV format (comma-separated values) with comma (`","`) as column separator and new line as row separator.

* Each row corresponds to a sample/patient.
* One column must contain the dependent variable.
* The remaining columns (independent variables) must be numeric.
* The column name of each independent variable should be an **Entrez gene ID**. All columns not matching a gene in the network will be disregarded.

Example file:

```
label,2099,351,7534,8452
case,0.40,0.50,0.88,0.81
control,0.42,0.95,0.31,0.15
control,0.91,0.99,0.88,0.57
case,0.72,0.19,0.32,0.16
case,0.19,0.85,0.03,0.37
control,0.971,0.78,0.96,0.18
```

If the dependent variable values contain commas, they must be enclosed in quotes:

```
group,2099,351,7534,8452
"case,day1",0.40,0.50,0.88,0.81
"case,day1",0.42,0.95,0.31,0.15
"case,day7",0.91,0.99,0.88,0.57
"case,day7",0.72,0.19,0.32,0.16
"control",0.19,0.85,0.03,0.37
"control",0.971,0.78,0.96,0.18
```

A valid file can be exported from R using the `write.csv` function:

```r
write.csv(x, row.names=FALSE)
```

## Training a model

## Model evaluation

## Prediction
