
require(data.table)
require(seewave)
require(reshape2)
require(DescTools)
#require(Hmisc)

source('../../Func/plotFFT.R')
source('../../Func/summaries.R')
source('../../Func/utils.R')

#patient=2
#sample='test'

build_features <-function(patient, sample) {
  #patient = as.numeric(unlist(strsplit(path,'_'))[2])
  path = getwd()
  features_path = paste0(path,'/',sample,'/features/')
  dt_path       = paste0(path,'/',sample,'/data_tables')
  
  file_list = list.files(path=dt_path,full.names = T)
  
  dt_features = data.table()
  
  i=0
  for ( afile in file_list ) {
    i=i+1
    cat('Processing file ', afile, '\n')
    
    load(afile)
    
    # remove rows with all zeros.
    rsum = rowSums(dt[,seq(1,16),with=FALSE], na.rm = T)
    all_zero_rows = rsum==0
    dt = dt[!all_zero_rows]
    
    if ( nrow(dt)==0 ) {
      next
    }
    
    cols = c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V13','V14','V15','V16')
    
    if ( 'target'%in%colnames(dt) ) {
      dt[,id:=paste0(N,'_',target)]
      dt[,N:=NULL]
      dt[,target:=NULL]
    } else {
      dt[,id:=N]
      dt[,N:=NULL]
    }
    
    # Winsorize.
    dt = dt[,lapply(.SD, Winsorize, probs=c(0.01,0.99)), by=id, .SDcols=cols]
    
    
    dt_means = dt[,lapply(.SD, mean), by=id, .SDcols=cols]
    setnames(dt_means, c('id',paste0('mean',cols)))
    
    dt_sds   = dt[,lapply(.SD, sd), by=id, .SDcols=cols]
    setnames(dt_sds, c('id',paste0('sd',cols)))
    
    dt_rob = dt[,as.list(unlist(lapply(.SD, rob_kurt_skew))), by=id, .SDcols=cols]
    robnames = c('id',paste0(rep(paste0(paste0('V',seq(1,16))),each=2),c('_skew','_kurt')))
    setnames(dt_rob,robnames)
    
    
    dt_cacf = dt[,lapply(.SD, cacf_calculator), by=id, .SDcols=cols]
    cacfnames = c('id',paste0(rep(paste0(paste0('V',seq(1,16)))),c('_cacf')))
    setnames(dt_cacf,cacfnames)
    
    
    # Till now we have 16*4 parameters.
    # for each electrode, build FFT spectrum. Average Power in each wave.
    #cat('Computing FFT...\n')
    ids = unique(dt$id)
    
    fft_summaries = data.table()
    for ( idi in ids ) {
      #cat(idi,'\n')
      fft_summary = dt[id==idi,as.list(unlist(lapply(.SD,fft_summary))), .SDcols=cols]
      fft_summary[,id:=idi]
      fft_summaries = rbind(fft_summary, fft_summaries)
      #rm(list=)
    }
    #fft_summaries = dt[,as.list(unlist(lapply(.SD,fft_summary))), by=id, .SDcols=cols]
    fftnames = c(paste0(rep(paste0(paste0('V',seq(1,16))),each=6),paste0('_band_',seq(1,6))), 'id')
    setnames(fft_summaries,fftnames)
    
     cat('Computing coherence...\n')
     coh_summaries = data.table()
     for ( idi in ids ) {
       cat(idi,'\n')
       coh_summary=dt[id==idi, as.list(unlist(coh_summary(.SD))), .SDcols=cols]
       coh_summary[,id:=idi]
       coh_summaries = rbind(coh_summary, coh_summaries)
     }
     #cohnames = c(paste0('coh_mean_',seq(1,6)),paste0('coh_sd_',seq(1,6)),'id')
     setnames(coh_summaries,cohnames)
    
    #cat('Computing Cross-Correlations...\n')
    # Cross-correlations.
    all_corr=data.table()
    for ( idi in ids ) {
      one_id = data.table(id=idi)
      for ( coli in seq(1,14) ) {
        for ( colj in seq(coli+2,16,2) ) {
          dt_mini = dt[id==idi]
          rows = seq(1,nrow(dt_mini),100)
          
          #ken  = cor(dt_mini[rows,coli,with=FALSE], dt_mini[rows,colj,with=FALSE], use="complete.obs", method="kendall")
          per  = cor(dt_mini[rows,coli,with=FALSE], dt_mini[rows,colj,with=FALSE], use="complete.obs", method="pearson")
          #sper = cor(dt_mini[rows,coli,with=FALSE], dt_mini[rows,colj,with=FALSE], use="complete.obs", method="spearman")
          
          #kendal_name   = paste0('ken_',coli,'_',colj)
          #person_name   = paste0('per_',coli,'_',colj)
          #spearman_name = paste0('sper_',coli,'_',colj)
          #one_id[,(paste0('ken_',coli,'_',colj)):=ken]
          one_id[,(paste0('per_',coli,'_',colj)):=per]
          #one_id[,(paste0('sper_',coli,'_',colj)):=sper]
        }
      }
      all_corr = rbind(one_id,all_corr)
    }
    
    #ccr_pears = rcorr(as.matrix(dt[,seq(1,16),with=FALSE]), type='pearson')
    #ccr_pears = ccr_pears$r[lower.tri(ccr_pears$r)]
    #ccr_spearman = rcorr(as.matrix(dt), type='spearman')
    #ccr_spearman = ccr_spearman$r[lower.tri(ccr_spearman$r)]
    
    
    #cat('Computing vol statistics...\n')
    # build mean(vol), var(vol), mean(volvol), var(volvol)
    dt[,t:=seq(1,nrow(.SD)),by=id]
    melted = melt(dt, id.vars = c('id','t'))
    setkeyv(melted,c('id','t'))
    #x= melted[id=='1',value]
    #melted[,lagged:=shift(value), by=id]
    #melted
    melted[,vol   :=as.numeric(volEWMA(value, N=2000)), by=c('id','variable')]
    melted[,volvol:=as.numeric(volEWMA(vol,   N=2000)), by=c('id','variable')]
    vol_sum = melted[,list(mean(vol),var(vol),mean(volvol),var(volvol)), by=c('id','variable')]
    setnames(vol_sum,c('id','variable','mean_vol','var_vol','mean_volvol','var_volvol'))
    
    t1=dcast(vol_sum[,list(id,variable,mean_vol)], id~variable,    value.var = 'mean_vol')
    setnames(t1,c('id',paste0(seq(1,16),'_mean_vol')))
    
    t2=dcast(vol_sum[,list(id,variable,var_vol)], id~variable,     value.var = 'var_vol')
    setnames(t2,c('id',paste0(seq(1,16),'_var_vol')))
    
    t3=dcast(vol_sum[,list(id,variable,mean_volvol)], id~variable, value.var = 'mean_volvol')
    setnames(t3,c('id',paste0(seq(1,16),'_mean_volvol')))
    
    t4=dcast(vol_sum[,list(id,variable,var_volvol)], id~variable,  value.var = 'var_volvol')
    setnames(t4,c('id',paste0(seq(1,16),'_var_volvol')))
    
    merged_vol_summary = merge(t1,t2, by='id')
    merged_vol_summary = merge(merged_vol_summary,t3, by='id')
    merged_vol_summary = merge(merged_vol_summary,t4, by='id')
    
    #vols = apply(dt,2,volEWMA)
    #vols = dt[,apply(.SD, volEWMA,), by=id, .SDcols=cols]
    #cat('Now merging everything...\n')
    #targets = unique(dt[,target,by=id])
    dt_tmp = merge(dt_means, dt_sds,   by='id')
    #dt_tmp = merge(dt_tmp, dt_mads,     by='id')
    #dt_tmp = merge(dt_tmp, dt_iqrs,     by='id')
    dt_tmp = merge(dt_tmp, dt_rob,     by='id')
    dt_tmp = merge(dt_tmp, dt_cacf,    by='id')
    dt_tmp = merge(dt_tmp, fft_summaries,      by='id')
    dt_tmp = merge(dt_tmp, coh_summaries,      by='id')
    dt_tmp = merge(dt_tmp, merged_vol_summary, by='id')
    dt_tmp = merge(dt_tmp, all_corr,           by='id')
    
    #dt_features = merge(dt_features, targets,  by='id')
    
    dt_features = rbind(dt_features,dt_tmp)
    
    
    if ( i%%2==0 ) {
      save(dt_features, file=paste0(features_path,i,'.RData'))
    }
    
    rm(list=c('dt','dt_cacf','dt_means','dt_sds','dt_tmp','merged_vol_summary','melted',
              'all_corr','dt_mini','coh_summaries',
              'dt_rob','fft_summaries'))
    gc()
    
    #fft_summaries[,i:=1]
    #fft_summaries[,i:=cumsum(i),by=id]
    #fft_summaries[,fft_var:=paste0(variable,'_fftwave_',i)]
    #fft_summaries = melt(fft_summaries,id.vars=c('id','i'))
    #fft_summaries = dcast(fft_summaries, id~fft_var)
    
    #coh_summaries = dt[,as.list(apply(.SD,coh_summary)), by=id, .SDcols=cols]
    
    #fft_summaries[,i:=1]
    #fft_summaries[,i:=cumsum(i),by=id]
    #fft_summaries = melt(fft_summaries,id.vars=c('id','i'))
    #fft_summaries[,fft_var:=paste0(variable,'_fftwave_',i)]
    #fft_summaries = dcast(fft_summaries, id~fft_var)
    
    #y = dt[id=='10',V1]
    #fft_summary(y)
    
    
    #plot(FFTdata[1:nrow(FFTdata)/2,])
    
    #cacf_names = paste0(rep(paste0('cacf',cols),each=4),rep(c('_1','_2','_3','_4')))
    #dt_cacf = dt[,lapply(.SD, cacf_calculator), by=N, .SDcols=cols]
    #setnames(dt_cacf, cacf_names)
    #dt[,id:=paste0(N,target)]
    #dt_mini = copy(dt[N%in%c(1,10)])
    #library(signal)
    #dt2 = dt_mini[,lapply(.SD, specgram, Fs=240000/600, n=2048), by='id', .SDcols=cols]
    
    #test = specgram(as.matrix(dt[id=='10', .SD, .SDcols=cols]), Fs=240000/600, n=2048)
    #test2=spectrum(as.matrix(dt[id=='10', .SD, .SDcols=cols]), )
    #test2$coh
    
    
    # first col is coherence in kHz
    #c=coh(dt[id==10,V1],dt[id==10,V2],f=400)
    # pick 50 bands? average? sum?
    
    
    #ts = dt[N==1,V16]
    #out = specgram(ts, Fs=240000/600, n=2048)
    #plot(out$S)
    
    # Cross-correlations.
    #ccr_pears = rcorr(as.matrix(dt), type='pearson')
    #ccr_pears = ccr_pears$r[lower.tri(ccr_pears$r)]
    #ccr_spearman = rcorr(as.matrix(dt), type='spearman')
    #ccr_spearman = ccr_spearman$r[lower.tri(ccr_spearman$r)]
  }
  #dt_features[,N:=as.numeric(unlist(strsplit(id,'_'))[1]),      by=id]  
  #dt_features[,target:=as.numeric(unlist(strsplit(id,'_'))[2]), by=id]
  #dt_features[,id:=NULL]
  
  save(dt_features,file = paste0(features_path,'features.RData'))
}


