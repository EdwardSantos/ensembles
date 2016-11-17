require(data.table)
require(R.matlab)

merge_files <- function(sample=NULL) {
  if ( is.null(sample) ) {
    cat('Input needed: sample must be train or test.')
  }
  
  path = getwd()
  
  patient = as.numeric(unlist(strsplit(path,'_'))[2])
  file_list = list.files(path=paste0(path,'/',sample,'/matlab'), full.names=T, pattern = '[1-9]')
  
  i = 0
  dt = data.table()
  dt_info = data.table()
  for ( afile in file_list ) {
    i = i + 1

    cat('File ',i,': ', afile, '\n')
    # Read.
    mat <- readMat(afile)
    dt_tmp = data.table(mat[[1]][[1]])

    parts = unlist(strsplit(afile, '_'))
    
    target_train = NA
    if ( sample=='train' ) {
      N = as.numeric(unlist(strsplit(parts[[3]], '.mat'))[1])
      
      target_train = as.numeric(substr(parts[[4]],0,1))
      dt_tmp[,target:=target_train]
      dt_tmp[,N:=N]
    } else {
      N = as.numeric(unlist(strsplit(parts[4], '[.]'))[1])
      dt_tmp[,N:=N]
    }
    
    
    dt = rbind(dt_tmp, dt)
    
    if ( i%%20==0 ) {
      cat('Saving data.table ', i/20,'\n')
      save(dt, file = sprintf('%s/dt_%d.RData',sample,i/20))
      rm(list='dt')
      gc()
      dt = data.table()
    }
    
    sampling_rate   = mat[[1]][[2]][[1]]
    nSamplesSegment = mat[[1]][[3]][[1]]
    
    tmp = data.table(file=N, target=target_train, sRate=sampling_rate, nSamplesSegment=nSamplesSegment)
    
    dt_info = rbind(dt_info, tmp)
    
    
    rm(list=c('mat','dt_tmp'))
    gc()
  }
  save(dt, dt_info, file=paste0(sample,'/dt_last.RData'))
}



