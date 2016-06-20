setwd("~/CoralProject")
DELTA = 0.00033333333333

#This function saves png files of barplots for binned data segments of a given "greasieness" for each sample
#It provides x axis labels for any data samples whose indices (from among 1 through ncol(df)) are specified in bottomPages
saveBinnedBarplots <- function(df, mat, binWidth, bottomPages, solvent){
  setwd("~/CoralProject")
  for(i in 1:nrow(mat)){
    png(filename=paste("Binned_Areas_",as.character(binWidth), "_PAL", as.character(df[i,1]),"_", solvent, "_", as.character(df[i,3]),".png", sep=""), height=100, width=400, bg="white")
    par(mar=c(.2,0,1,.01), family="serif")
    barplot(mat[i,],col="skyblue",ylab = "Area",xlab="Bin", names.arg=c(1:ncol(mat)), bty="n")
    mtext(text=as.character(c(1:ncol(mat))), at=seq(.55,.55+ncol(mat)*1.205, 1.205), line=-11, cex=.5)
    legend(x="topright", legend=paste("PAL", as.character(df[i,1]), sep=""), bty="n")
    dev.off()
  }
  for(j in bottomPages){
    png(filename=paste("Binned_Areas_",as.character(binWidth), "_PAL", as.character(df[i,1]),"_", solvent, "_", as.character(df[i,3]),".png", sep=""), height=100, width=400, bg="white")
    par(mar=c(3,0,1,.01), family="serif")
    barplot(mat[i,],col="skyblue",ylab = "Area",xlab="Bin", names.arg=c(1:ncol(mat)), axes=FALSE, bty="n")
    legend(x="topright", legend=paste("PAL", as.character(df[i,1]), sep=""), bty="n")
    dev.off()
  }
}

saveBinnedBarplots(dfAllT_DCM, binned_matrix_DCM, 1100, c(7, 15, 23, 24, 31, 34, 41, 45), "DCM")
saveBinnedBarplots(dfAllT_HEX, binned_matrix_HEX, 1100, c(7, 15, 18 , 25, 26, 33), "HEX")

#This function saves png files of gas chromatogram plots for all data samples of a given "greasiness" (solvent).
#It also provides x axis labels for any data samples whose indices (from among 1 through ncol(df)) are specified in bottomPageIndices
saveGasChromatogramPlot <- function(df, bottomPageIndices, solvent){
  setwd("~/CoralProject")
  for(i in 1:nrow(df)){
    png(filename=paste("Chromatogram_PAL_",as.character(df[i,1]),"_", solvent,"_",as.character(df[i,3]) ,".png", sep=""), height=80, width=400, bg="white")
    par(mar=c(.5,.4,.4,.5), family="serif")
    plot(x=seq(1,(1 + DELTA*(ncol(df)-4)),DELTA), y=df[i,4:ncol(df)],col="black", type="l", xlab=NA, ylab=NA, bty="n", axes=FALSE, xaxs="i", yaxs="i")
    legend(x="topright", legend=paste("PAL", as.character(df[i,1]), sep=""), bty="n")
    dev.off()
  }
  for(i in bottomPageIndices){
    png(filename=paste("Chromatogram_PAL_",as.character(df[i,1]),"_", solvent,"_",as.character(df[i,3]) ,".png", sep=""), height=80, width=400, bg="white")
    par(mar=c(2.6,.4,.4,.5), family="serif")
    plot(x=seq(1,(1 + DELTA*(ncol(df)-4)),DELTA), y=df[i,4:ncol(df)],col="black", type="l", ylab=NA, bty="n", xaxs="i", yaxs="i")
    legend(x="topright", legend=paste("PAL", as.character(df[i,1]), sep=""), bty="n")
    dev.off()
  }
}


saveGasChromatogramPlot(dfAllT_DCM, c(7, 15, 23, 24, 31, 34, 41, 45), "DCM")
saveGasChromatogramPlot(dfAllT_HEX, c(7, 15, 18 , 25, 26, 33), "HEX")





#plots of original time data with bin lines shown
#xadj is the amount by which we push the label of the bin over from the line 
#yheight is the heigh at which the bin labels are displayed
#binWidth
drawPlotBinLines <- function(df, row, binWidth, xadj=(0.000175*binWidth), yheight=floor(max(as.numeric(df[row, TRUNC_BEGIN:TRUNC_END]))), charSize=1, binLabels=FALSE){
  par(mar=c(.1,.1,.4,.1), family="serif")
  times = seq(1,26.6667, DELTA)
  a=seq(6, 6+(TRUNC_END-TRUNC_BEGIN)*DELTA,binWidth*DELTA)
  a=a[-length(a)]
  b=rep(yheight,length(a))
  plot(df[row,TRUNC_BEGIN:TRUNC_END],pch=20,col="white",xlim=c(5,21.6667),frame=FALSE)
  lines(times[TRUNC_BEGIN:TRUNC_END],df[row, TRUNC_BEGIN:TRUNC_END],lwd=2)
  abline(v=seq(6, 6+(TRUNC_END-TRUNC_BEGIN)*DELTA,binWidth*DELTA), col="steelblue1")
  if(binLabels){
    text(x=a+xadj, y=b,labels=as.character(c(1:length(b))),col="grey42", cex=charSize)
  }
  legend(x="topright", legend=paste("PAL", as.character(df[i,1]), sep=""), bty="n")
}


binWidth = 1100
#DCM truncated and binned chromatograms
for(i in 1:nrow(dfAllT_DCM)){
  png(filename=paste("Bin_Trunc_Chrom_",as.character(binWidth), "_PAL_",as.character(dfAllT_DCM[i,1]),"_DCM_",as.character(dfAllT_DCM[i,3]) ,".png", sep=""), height=80, width=400, bg="white")
  drawPlotBinLines(dfAllT_DCM, row=i, binWidth=1100, binLabels = TRUE, charSize = .5)
  dev.off()
}
  
#HEX truncated and binned chromatograms
for(i in 1:nrow(dfAllT_HEX)){
  png(filename=paste("Bin_Trunc_Chrom_", as.character(binWidth), "_PAL_",as.character(dfAllT_HEX[i,1]),"_HEX_",as.character(dfAllT_HEX[i,3]) ,".png", sep=""), height=80, width=400, bg="white")
  drawPlotBinLines(dfAllT_HEX, row=i, binWidth=1100, binLabels = TRUE, charSize = .5)
  dev.off()
}



BINWIDTH = 2000
#DCM truncated and binned chromatograms
for(i in 1:nrow(dfAllT_DCM)){
  png(filename=paste("Bin_Trunc_Chrom_",as.character(BINWIDTH), "_PAL_",as.character(dfAllT_DCM[i,1]),"_DCM_",as.character(dfAllT_DCM[i,3]) ,".png", sep=""), height=80, width=400, bg="white")
  drawPlotBinLines(dfAllT_DCM, row=i, binWidth=BINWIDTH, binLabels = TRUE, charSize = .5)
  dev.off()
}

#HEX truncated and binned chromatograms
for(i in 1:nrow(dfAllT_HEX)){
  png(filename=paste("Bin_Trunc_Chrom_", as.character(binWidth), "_PAL_",as.character(dfAllT_HEX[i,1]),"_HEX_",as.character(dfAllT_HEX[i,3]) ,".png", sep=""), height=80, width=400, bg="white")
  drawPlotBinLines(dfAllT_HEX, row=i, binWidth=BINWIDTH, binLabels = TRUE, charSize = .5)
  dev.off()
}

saveBinnedBarplots(dfAllT_DCM, binned_matrix_DCM, BINWIDTH, c(7, 15, 23, 24, 31, 34, 41, 45), "DCM")
saveBinnedBarplots(dfAllT_HEX, binned_matrix_HEX, BINWIDTH, c(7, 15, 18 , 25, 26, 33), "HEX")


BINWIDTH = 3000
#DCM truncated and binned chromatograms
for(i in 1:nrow(dfAllT_DCM)){
  png(filename=paste("Bin_Trunc_Chrom_",as.character(BINWIDTH), "_PAL_",as.character(dfAllT_DCM[i,1]),"_DCM_",as.character(dfAllT_DCM[i,3]) ,".png", sep=""), height=80, width=400, bg="white")
  drawPlotBinLines(dfAllT_DCM, row=i, binWidth=BINWIDTH, binLabels = TRUE, charSize = .5)
  dev.off()
}

#HEX truncated and binned chromatograms
for(i in 1:nrow(dfAllT_HEX)){
  png(filename=paste("Bin_Trunc_Chrom_", as.character(binWidth), "_PAL_",as.character(dfAllT_HEX[i,1]),"_HEX_",as.character(dfAllT_HEX[i,3]) ,".png", sep=""), height=80, width=400, bg="white")
  drawPlotBinLines(dfAllT_HEX, row=i, binWidth=BINWIDTH, binLabels = TRUE, charSize = .5)
  dev.off()
}

saveBinnedBarplots(dfAllT_DCM, binned_matrix_DCM, BINWIDTH, c(7, 15, 23, 24, 31, 34, 41, 45), "DCM")
saveBinnedBarplots(dfAllT_HEX, binned_matrix_HEX, BINWIDTH, c(7, 15, 18 , 25, 26, 33), "HEX")

