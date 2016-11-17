library(data.table)
library(R.matlab)

source('../../Func/plotFFT.R')

plotEEG <- function(max.files=150, train=TRUE, patient=1) {
  
  path = getwd()
  
  patient = as.numeric(unlist(strsplit(path,'_'))[2])
  file_list = list.files(path=path, pattern = '[1-9]')
  
  x <- seq(0, 0.5*pi, 0.01)
  
  if ( train ) {
    file_list = paste(patient,'_', c(1:150), "_",c(0,1), ".mat", sep = "")
  }
  file_list=file_list[1:max.files]
  
  j = 0
  for ( afile in file_list ) {
    j = j + 1
    
    cat('File ',j,': ', afile, '\n')
    # Read.
    mat <- readMat(afile)
    dt_tmp = data.table(mat[[1]][[1]])
    #parts = unlist(strsplit(afile,'_'))
    
    #N = as.numeric(parts[[3]])
    #filename=paste0(afile,'.png')
    
    png(file=paste0('Plots/',afile,'.png'), width = 40, height = 20, units = 'cm',re=600)
    par(mfrow=c(8,4))
    par(mar=c(1,1,1,1))
    for ( i in seq(1,16) ) {
      plot(unlist(dt_tmp[,i,with=FALSE]),type='l',ylab = '',xlab = '')
      res <- plotFFT(x, unlist(dt_tmp[,i,with=FALSE]), 400)
    }
    dev.off()
    
    #png(file=paste0(filename,'fft.png'), width = 6, height = 6, units = 'cm')
    #par(mfrow=c(2,2))
    #for ( i in seq(1,16) ) {
    #  res <- plotFFT(x, unlist(dt_tmp[,i,with=FALSE]), 400)
    #}
    #dev.off()
    
    rm(list=c('mat','dt_tmp'))
    gc()
  }
}
