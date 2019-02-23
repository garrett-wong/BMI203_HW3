setwd("~/Documents/winterY1/algorithms/homeworks/homework3/HW3_due_02_23")
myData <- read.table("optimizeGapPenalties.txt")
colnames(myData) <- c("openGapCost", "extendGapCost", "falsePositiveRate")
library(ggplot2)
ggplot(data=myData, aes(x=openGapCost, y=extendGapCost, fill=falsePositiveRate)) + geom_tile()

myData[myData$falsePositiveRate == min(myData$falsePositiveRate),]
