rm(list=ls())
gc()
require(data.table)
require(xgboost)
require(MASS)
#require(rnn)
require(randomForest)
#require(mice)
require(caret)
require(plyr)

source('~/local/kaggle/eeg/model_utils.R')


build_global_model <- function() {
  cat('Building global models...\n')
  data = load_all_data()
  dt_train = data[['train']]
  dt_test  = data[['test']]
  
  cat('Removing Kendal and Spearman features...\n')
  cols_2_from_test = grep("ken", colnames(dt_test))
  cols_2_from_train = grep("ken", colnames(dt_train))
  dt_train[,c(cols_2_from_train):=NULL]
  dt_test[,c(cols_2_from_test):=NULL]
  
  cols_2_from_test = grep("sper", colnames(dt_test))
  cols_2_from_train = grep("sper", colnames(dt_train))
  dt_train[,c(cols_2_from_train):=NULL]
  dt_test[,c(cols_2_from_test):=NULL]
  
  simple_models = c('xgboost','rf','nnet')
  #simple_models = c('nnet')
  for ( model in simple_models ) {
    standalone_models(dt_train, dt_test, model, suffix='cor_reduced')
  }
  
  # dont like to use channel specific variables in a global model... 
  # summarise things in 10 deciles...
  #meanvolvol_cols_test  = grep("mean_volvol", colnames(dt_test))
  #dt_train[,i:=seq(1,nrow(dt_train))]
  #meanvolvol_cols_train = grep("mean_volvol",colnames(dt_train))
  #ord = order(dt_train[4,c(meanvolvol_cols_train),with=FALSE])
}

build_balanced_models <- function(N_models=3) {
  cat('Building patient specific models...')
  data = load_data()
  dt_train = data[['train']]
  dt_test  = data[['test']]
  
  cat('Removing Kendal and Spearman features...\n')
  cols_2_from_test = grep("ken", colnames(dt_test))
  cols_2_from_train = grep("ken", colnames(dt_train))
  dt_train[,c(cols_2_from_train):=NULL]
  dt_test[,c(cols_2_from_test):=NULL]
  
  cols_2_from_test = grep("sper", colnames(dt_test))
  cols_2_from_train = grep("sper", colnames(dt_train))
  dt_train[,c(cols_2_from_train):=NULL]
  dt_test[,c(cols_2_from_test):=NULL]
  
  # remove vol and volvol features...
  #cols_2_from_test = grep("vol", colnames(dt_test))
  #cols_2_from_train = grep("vol", colnames(dt_train))
  #dt_train[,c(cols_2_from_train):=NULL]
  #dt_test[,c(cols_2_from_test):=NULL]
  
  #for ( model in simple_models ) {
  #simple_models = c('xgboost','rf','nnet')
  #  standalone_models(dt_train, dt_test, model, suffix='all')
  #}
  
  #features_2_rm = find_bad_features(colnames(dt_test))
  #cat('Removing features ', features_2_rm,'\n')
  
  #dt_train[,c(features_2_rm):=NULL]
  #dt_test[,c(features_2_rm):=NULL]
  #
  #simple_models = c('xgboost','rf','nnet')
  #for ( model in simple_models ) {
  #  standalone_models(dt_train, dt_test, model, suffix='trimmed')
  #}
  
  #return(0)
  
  # From the above, decided:
  # - use probability/raw,
  # - use all features/trim
  ######################################
  ##              WELCOME             ##
  ##      to the BEAUTIFUL world      ##
  ##  of bagged & stacked estimation  ##
  ######################################
  cat('Bagging and Stacking starts...')
  
  ###################################
  # Handle the imbalance.
  train_1 = copy(dt_train[target=='1'])
  train_0 = copy(dt_train[target=='0'])
  N = round(nrow(train_0)/nrow(train_1))
  # break train_0 into N parts:
  all_rows = seq(1,nrow(train_0))
  dts = list()
  L = round(nrow(train_0)/N)
  for ( i in 1:(N-1) ) {
    sel = sample(all_rows, L)
    dts[[i]] <- rbind(train_1,train_0[sel])
    dts[[i]] = dts[[i]][sample(nrow(dts[[i]]))]
    all_rows = all_rows[!(all_rows%in%sel)]
  }
  dts[[N]] = rbind(train_1,train_0[all_rows])
  dts[[N]] = dts[[N]][sample(nrow(dts[[N]]))]
  
  ###################################
  #rnn_model = model_rnn(dt_train, dt_test)
  
  ###################################
  # Build N models of each desired type. 
  models = c('xgboost','rf','nnet','svm','gbm')[1:N_models]
  all_raw_models = list()
  for ( i in 1:(N-2) ) {
    cat('Building Raw Models...\n')
    raw_model = model_caret(dts[[i]], models, i)
    all_raw_models[[i]] = raw_model
  }
  
  ###################################
  N_samples = length(all_raw_models)
  N_models  = length(all_raw_models[[1]])
  model_names = names(all_raw_models[[1]])
  
  all_models = list()
  for ( i in 1:N_models ) {
    model_name = model_names[i]
    all_models[model_name] <- NA
    for ( j in 1:N_samples ) {
      all_models[[model_name]] = c(all_models[[model_name]], list(all_raw_models[[j]][[model_name]]))
    }
    all_models[[model_name]][1] <-NULL # remove NAs
  }
  
  saved = copy(all_models)
  ###################################
  ## Merge all the raw models.
  all_models  = bag_models(all_models, dts[[N-1]])
  
  all_models  = stack_models(all_models, dts[[N]])
  
  # Predict test using bagged models (i.e., 1 prediction from xgboost, another from RF,...
  ###################################
  ## Impute non-finite values.
  for (j in 1:ncol(dt_test)) set(dt_test, which(is.infinite(dt_test[[j]])), j, NA)
  for (j in names(dt_test)) {
    set(dt_test,which(is.na(dt_test[[j]])),j,mean(dt_test[[j]],na.rm=T))
  }
  
  path = getwd()
  patient = as.numeric(unlist(strsplit(path,'_'))[2])
  filename = paste0('new_',patient,'_',dt_test$id,'.mat')
  ##################################################
  N_models   = length(all_models)-1
  N_samples  = length(all_models[[1]])-1
  model_names = names(all_models)
  all_predictions = data.table(dummy=NA)
  for ( m in 1:N_models ) { # For each model...
    model_name = model_names[m]
    bag_predictions = data.table(dummy=NA) # ... will have a dt with predictions...
    for ( i in 1:N_samples ) {
      prediction = predict(all_models[[model_name]][[i]], newdata=dt_test, type='prob')$S
      bag_predictions = cbind(bag_predictions, prediction)
    }
    bag_predictions[,dummy:=NULL]
    setnames(bag_predictions,paste0('B',1:N_samples))
    prediction = predict(all_models[[model_name]][['bagged']], newdata=bag_predictions, type='response')
    prediction = data.table(prediction)
    setnames(prediction,model_name)
    all_predictions = cbind(all_predictions, prediction)
    
    sub_bagged <- data.table(File=filename, Class=prediction)
    write.csv(sub_bagged, file=sprintf('models/bagged_prediction_%s.csv',model_name), row.names=F, quote=F)
  }
  all_predictions[,dummy:=NULL]

  stacked_prediction <- predict(all_models[['stacked']], newdata=all_predictions,type='response')
  ##################################################
  sub_final <- data.table(File=filename, Class=stacked_prediction)
  write.csv(sub_final, file='models/stacked_prediction.csv', row.names=F, quote=F)
}


stack_models <- function(all_models, train) {
  train_bagged = copy(train)
  for (j in 1:ncol(train_bagged)) set(train_bagged, which(is.infinite(train_bagged[[j]])), j, NA)
  for (j in names(train_bagged)) {
    set(train_bagged,which(is.na(train_bagged[[j]])),j,mean(train_bagged[[j]],na.rm=T))
  }
  
  train_bagged[target=='0',targetF:='N']
  train_bagged[target=='1',targetF:='S']
  train_bagged[,targetF:=as.factor(targetF)]
  targets = train_bagged[,targetF]
  
  train_bagged[,target:=NULL]
  train_bagged[,targetF:=NULL]
  
  ## Here we go...
  #stacked_model = list()
  N_models   = length(all_models)
  N_samples  = length(all_models[[1]])-1
  model_names = names(all_models)
  all_predictions = data.table(dummy=NA)
  for ( m in 1:N_models ) { # For each model...
    model_name = model_names[m]
    bag_predictions = data.table(dummy=NA) # ... will have a dt with predictions...
    for ( i in 1:N_samples ) {
      prediction = predict(all_models[[model_name]][[i]], newdata=train_bagged, type='prob')$S
      bag_predictions = cbind(bag_predictions, prediction)
    }
    bag_predictions[,dummy:=NULL]
    setnames(bag_predictions,paste0('B',1:N_samples))
    prediction = predict(all_models[[model_name]][['bagged']], newdata=bag_predictions, type='response')
    prediction = data.table(prediction)
    setnames(prediction,model_name)
    all_predictions = cbind(all_predictions, prediction)
  }
  all_predictions[,dummy:=NULL]
  all_predictions[,target:=targets]

  stacked_model <- glm(target~., 
                       family=binomial(link='logit'), 
                       data=all_predictions,
                       control = list(maxit = 100))
  
  tmp = summary(stacked_model)
  capture.output(tmp,file = sprintf('models/glm_stacking_summary.txt'))
  cat('Finished stacking!\n') 
  print(tmp)
  
  all_models[['stacked']] = stacked_model
  
  return(all_models)
}

bag_models <- function(all_models, train) {
  train_bagged = copy(train)
  for (j in 1:ncol(train_bagged)) set(train_bagged, which(is.infinite(train_bagged[[j]])), j, NA)
  for (j in names(train_bagged)) {
    set(train_bagged,which(is.na(train_bagged[[j]])),j,mean(train_bagged[[j]],na.rm=T))
  }
  
  train_bagged[target=='0',targetF:='N']
  train_bagged[target=='1',targetF:='S']
  train_bagged[,targetF:=as.factor(targetF)]
  targets = train_bagged[,targetF]
  
  train_bagged[,target:=NULL]
  train_bagged[,targetF:=NULL]
  
  ## Here we go...
  bagged_models = list()
  N_models  = length(all_models)
  N_samples = length(all_models[[1]])
  model_names = names(all_models)
  for ( m in 1:N_models ) { # For each model...
    model_name = model_names[m]
    bag_predictions = data.table(dummy=NA) # ... will have a dt with predictions...
      for ( i in 1:N_samples ) {
        prediction = predict(all_models[[model_name]][[i]], newdata=train_bagged, type='prob')$S
        bag_predictions = cbind(bag_predictions, prediction)
      }
    bag_predictions[,dummy:=NULL]
    setnames(bag_predictions,paste0('B',1:N_samples))
    bag_predictions[,target:=targets]
      
    model <- glm(target~., 
                 family=binomial(link='logit'), 
                 data=bag_predictions,
                 control = list(maxit = 1000))
      
    tmp = summary(model)
    capture.output(tmp,file = sprintf('models/glm_bagging_summary_%s.txt', model_name))
    cat('Finished bagging ', model_name, '\n') 
    cat('Summary:\n')
    print(tmp)
      
    all_models[[model_name]][['bagged']] = model
  }
  return(all_models)
}
