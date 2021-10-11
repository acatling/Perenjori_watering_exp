---
title: "README.md"
author: "Alexandra Catling"
date: "06/10/2021"
output: html_document
---

## README
Repository for WA Perenjori Watering Experiment
PhD Research
Experiment conducted in 2020

data_preparation sheet contains all dataframes needed to run plots and models in other sheets
WA_questions_analysis is the most up to date analysis of the paper's questions
dataall is the dataframe containing most of the information: row per individual with Species, Site, Plot, neighbour abundance, survival, abiotic environment info, trait data.



```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r cars}
summary(cars)
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
