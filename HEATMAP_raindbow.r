library(ggplot2)
library(reshape)
library(gplots)
library(RColorBrewer)
library(grid)
library(colorRamps)


setwd("/Users/davidgae/Desktop/L12_SIM/")
data <- read.csv("./heatmap.csv", comment.char="#")

n_matrix <- data.matrix(data)
pdf("SEQ_SIGheatmap.pdf")
y <- melt(n_matrix)

y$value1 <- cut(y$value,breaks = c(0,0.3,0.6,0.9,1.0),include.lowest=TRUE) #breaks the
myPalette <- colorRampPalette(matlab.like2(256), space="Lab")


p <- ggplot(na.omit(y),aes(x=X1, y=X2),reverse=TRUE) #omit NA values (na.omit).

p + geom_tile(aes(fill=value1))+scale_y_discrete(limits=unique(y$X2))+scale_x_discrete(expand= c(0,0),limits=c(0,50,100,150,200,250,300,350,400,450,500,550,600))+scale_fill_manual(labels = c("0","0.1-0.3","0.31-0.60","0.61-0.9","0.91-1"),values=c("#FCFCFC","#00FF00","#FFFF00","#FF0000"),guide = guide_legend(reverse=TRUE, title=NULL, direction = "vertical",label.hjust = 0.5, label.vjust = 0.5, label.theme = element_text(angle = 90),size=0.1))+ xlab("") + ylab("") +theme(legend.position="right")+ theme(legend.text=element_text(size=1))+theme(axis.text.x=element_text(angle=90),axis.text.y=element_text(angle=180))

dev.off()
