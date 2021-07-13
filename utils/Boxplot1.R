#---- BoxPlot ------
library(ggplot2)
load("Boxdata.Rdata")
ggplot(Boxdata, aes(x = methods, y = AveLoss, fill=Initial)) + geom_boxplot() + 
    ylab("Average error") + theme(legend.key.size = unit(1, 'cm')) + ylim(0.1, 1.65) + 
    facet_wrap(~DataType) + theme(legend.position="top") 
