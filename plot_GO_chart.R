###############################################################################
###############################################################################
#script to design the pie charts
###############################################################################
###############################################################################

#this script is largelly inspired by this post by Arsalvacion: 
#http://r-nold.blogspot.fi/2012/10/nscb-sexy-stats-version-2.html


library(RColorBrewer)


op<-par(mfrow=c(1,2),mar=c(1,4.1,2,0.1))

level2<-read.table(file="chart_data_20121107_1851_level2.txt", header=T, sep="\t")
pie(level2$Score, label=paste(level2$GO.term," (",level2$Score,")"),
    init.angle=270,radius=0.7, col=brewer.pal(6,"Set1"), cex=0.6)
par(new=TRUE)
pie(c(1), labels=NA, radius=0.4)
par(new=TRUE)
pie(c(1), labels=NA, border='white', radius=0.39)
text(0,0,labels="Gene Ontology\nBiological Process\nLevel 2", cex=1, font=2)

level3<-read.table(file="modified_level3.txt", header=T, sep="\t")
pie(level3$Score, label=paste(level3$GO.term," (",level3$Score,")"),
    init.angle=5,radius=0.7, col=brewer.pal(7,"Set1"), cex=0.6)
par(new=TRUE)
pie(c(1), labels=NA, radius=0.4)
par(new=TRUE)
pie(c(1), labels=NA, border='white', radius=0.39)
text(0,0,labels="Gene Ontology\nBiological Process\nLevel 3", cex=1, font=2)

par(op)

