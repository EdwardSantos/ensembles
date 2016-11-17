# Filter duplicates.


# Loop over all ids. For each new id, sample 20 numbers from V2. Compare with  

#patient = as.numeric(unlist(strsplit(path,'_'))[2])
set.seed(42)
rows_selected = sort(round(runif(2000, 1,120000)))
path = getwd()

for ( sample in c('train','test') ) {
  dt_path       = paste0(path,'/',sample,'/data_tables')
  file_list = list.files(path=dt_path,full.names = T)
  
  i=0
  for ( afile in file_list ) {
    i=i+1
    cat('Processing file ', afile, '\n')
    load(afile)
    
    # V13 seems to have highest variance.
    
    dt = dt[,list(N,target,V13)]
    dt[,V13:=c(0.0,diff(V13)), by=N]
    dt = dt[,.SD[rows_selected], by=N]
    setkeyv(dt, 'N')
    
    dt[,var(.SD)]
    
    cols = c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V13','V14','V15','V16')
    
    dt_sds   = dt[,lapply(.SD, sd), by=N, .SDcols=cols]
    
    which.max(colMeans(dt_sds)[seq(2,17)])
    setnames(dt_sds, c('id',paste0('sd',cols)))
    
    
} 


dt_features = data.table()

i=0

  
  
  
}





