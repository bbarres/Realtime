##############################################################################/
##############################################################################/
#Plot the output of BayeScan software
##############################################################################/
##############################################################################/

#set the working directory



##############################################################################/
#Define the function####
##############################################################################/

#this function was retrieve from the BayeScan software package and a few 
#modification were made

plot_bayescan<-function(res,FDR=0.05,size=1,pos=0.35,highlight=NULL,
                        name_highlighted=F,add_text=T) {
  if (is.character(res))
    res=read.table(res)
  
  colfstat=5
  colq=colfstat-2
  
  highlight_rows=which(is.element(as.numeric(row.names(res)),highlight))
  non_highlight_rows=setdiff(1:nrow(res),highlight_rows)
  
  outliers=as.integer(row.names(res[res[,colq]<=FDR,]))
  
  ok_outliers=TRUE
  if (sum(res[,colq]<=FDR)==0)
    ok_outliers=FALSE;
  res[res[,colq]<=0.0001,colq]=0.0001
  
  # plot
  plot(log10(res[,colq]),res[,colfstat],xlim=rev(range(log10(res[,colq]))),
       xlab="log10(q value)",ylab="Fst",type="n",las=1,cex.axis=0.8)
  points(log10(res[non_highlight_rows,colq]),res[non_highlight_rows,
                                                 colfstat],pch=19,cex=size)
  
  if (name_highlighted) {
    if (length(highlight_rows)>0) {
      text(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],
           row.names(res[highlight_rows,]),col="red",cex=size*1.2,font=2)
    }
    
  }
  
  else {
    points(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],
           col="red",pch=19,cex=size)
    # add names of loci over p and vertical line
    
    if (ok_outliers & add_text) {
      text(log10(res[res[,colq]<=FDR,][,colq])+
             pos*(round(runif(nrow(res[res[,colq]<=FDR,]),1,2))*2-3),
           res[res[,colq]<=FDR,][,colfstat],row.names(res[res[,colq]<=FDR,]),
           cex=size)
    }
  }
  lines(c(log10(FDR),log10(FDR)),c(-1,1),lwd=2)
  return(list("outliers"=outliers,"nb_outliers"=length(outliers)))
}


##############################################################################/
#Computation and plotting####
##############################################################################/

plot_bayescan("data/lim_vs_natrenf_fst.txt",pos=0.35,FDR=0.77,
              name_highlighted=FALSE,add_text=FALSE,
              highlight=97)
# if you save the output in a variable, you can recall the different results:
# results<-plot_bayescan("output_fst.txt",0,FDR=0.05)
# results$outliers
# results$nb_outliers


##############################################################################/
#END
##############################################################################/