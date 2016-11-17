library('NeuralNetTools')

##################################################################################
#model_caret <- function(dt_train, dt_log_train, dt_test, models=c('xgboost'), partition=0) {
model_caret <- function(dt_train, models=c('xgboost'), partition=0) {
  train  = copy(dt_train)
  #train2 = copy(dt_log_train)
  #test   = copy(dt_test)
  
  ###################################
  ## Impute non-finite values.
  for (j in 1:ncol(train)) set(train, which(is.infinite(train[[j]])), j, NA)
  for (j in names(train)) {
    set(train,which(is.na(train[[j]])),j,mean(train[[j]],na.rm=T))
  }
  
  #train2[,target:=NULL]
  train[target=='0',targetF:='N']
  train[target=='1',targetF:='S']
  train[,targetF:=as.factor(targetF)]
  train[,target:=NULL]
  
  raw_models = list()
  
  trcontrol <- trainControl(method="cv",
                            number = 3,
                            verboseIter = FALSE,
                            returnData=FALSE,
                            returnResamp = "all",
                            allowParallel = FALSE,
                            summaryFunction=twoClassSummary,
                            classProbs=T)
  ###################################
  if ( 'gbm' %in% models ) {
    cat('   Training gbm...\n')
    #gbm_grid_1 <- expand.grid(n.trees=c(500,1000),interaction.depth=c(6,15), shrinkage = 0.01,n.minobsinnode=c(3,6,10))
    gbm_grid_1 <- expand.grid(n.trees=c(500),interaction.depth=c(6), shrinkage = 0.01,n.minobsinnode=c(6))
    
    gbm_model <- train(targetF~., 
                       data=train, 
                       method="gbm", 
                       trControl = trcontrol,
                       tuneGrid  = gbm_grid_1,
                       metric="ROC")
    
    raw_models[['gbm']] = gbm_model
  }
  
  ###################################
  if ( 'nnet' %in% models ) {
    cat('   Training nnet...\n')
    nnet_grid_1 <- expand.grid(size=2,decay=0.5)
    #nnet_grid_1 <- expand.grid(size=c(2,5,7),decay=c(0.5,0.1))
    
    nnet_model <- train(targetF~., 
                        data=train, 
                        method="nnet", 
                        trControl = trcontrol,
                        tuneGrid  = nnet_grid_1,
                        preProcess = c("center", "scale"),
                        metric="ROC")
    
    raw_models[['nnet']] = nnet_model
  }
  
  ###################################
  if ( 'svmRadial' %in% models ) {
    cat('   Training SVM...\n')
    svm_grid_1 <- expand.grid(sigma=2^(-20), C=2^(2))
    #svm_grid_1 <-expand.grid(sigma= 2^c(-25, -20, -15,-10, -5, 0), C= 2^c(0:5))
    
    svm_model <- train(targetF~., 
                       data=train, 
                       method="svmRadial", 
                       trControl = trcontrol,
                       tuneGrid  = svm_grid_1,
                       preProc = c("center", "scale"),
                       metric="ROC")
    
    raw_models[['svm']] = svm_model
  }
  
  ###################################
  if ( 'rf' %in% models ) {
    cat('   Training RF...\n')
    #rf_grid_1 <- expand.grid(mtry=c(5,10,20))
    rf_grid_1 <- expand.grid(mtry=c(5))
    
    #names(getModelInfo())
    # what is strata?=factor(...)
    
    rf_model <- train(targetF ~., 
                      data=train, 
                      method="rf", 
                      ntree=1000,
                      trControl = trcontrol,
                      tuneGrid  = rf_grid_1,
                      metric="ROC")
    
    raw_models[['rf']] = rf_model
  }
  
  ###################################
  if ( 'xgboost' %in% models ) {
    cat('   Training xgboost...\n')
    xgb_grid_1 <- expand.grid(nrounds= 100,
                              max_depth=c(6),
                              gamma=c(0),
                              min_child_weight=c(0.5),
                              colsample_bytree=c(0.5),
                              eta=c(0.01))
    #xgb_grid_1 <- expand.grid(nrounds= 10000,
    #                          max_depth=c(6,10,15),
    #                          gamma=c(0,1),
    #                          min_child_weight=c(0.1,0.5),
    #                          colsample_bytree=c(0.1,0.5),
    #                          eta=c(0.01))
    
    xgb_model <- train(targetF~.,
                       data=train,
                       method="xgbTree",
                       trControl = trcontrol,
                       tuneGrid  = xgb_grid_1,
                       metric="ROC")
    
    raw_models[['xgboost']] = xgb_model
  }
  return(raw_models)
}
##################################################################################

##################################################################################
load_data <- function() {
  path = getwd()
  patient = as.numeric(unlist(strsplit(path,'_'))[2])
  
  set.seed(42)
  
  safe = fread('../train_and_test_data_labels_safe.csv')
  safe[,patient:=as.numeric(substr(x = image, 0,1))]
  safe[,N:=lapply(image, f<-function(x) as.numeric(unlist(strsplit(x,'_'))[2]))]
  safe = safe[!is.na(N)]
  safe = safe[patient==patient,list(safe,N)]
  safe[,N:=as.numeric(N)]
  
  load('train/features/features.RData')
  dt_train=copy(dt_features)
  
  dt_train[,N:=as.numeric(unlist(strsplit(id,'_'))[1]),      by=id] 
  dt_train[,target:=as.numeric(unlist(strsplit(id,'_'))[2]), by=id]
  Ns_2_rm = safe[safe==0,N]
  dt_train = dt_train[!(N%in%Ns_2_rm)]
  dt_train[,id:=NULL]
  dt_train[,N:=NULL]
  
  # Load test data.
  load('test/features/features.RData')
  dt_test = copy(dt_features)
  return(list(train=dt_train, test=dt_test))
}
##################################################################################

##################################################################################
save_prediction <- function(method, ids, prediction_test, prediction_train) {
  path = getwd()
  patient = as.numeric(unlist(strsplit(path,'_'))[2])
  
  filename = paste0('new_',patient,'_',ids,'.mat')
  sub <- data.table(File=filename, Class=prediction_test)
  write.csv(sub, file = sprintf("models/%s.csv",method), row.names=F, quote=F)
  
  write.csv(prediction_train, file=sprintf("models/%s_train.csv",method), row.names=F, quote=F)
}
##################################################################################

##################################################################################
model_randomForest <- function(dt_train, dt_test) {
    train = copy(dt_train)
    test  = copy(dt_test)
    
    for (j in 1:ncol(train)) set(train, which(is.infinite(train[[j]])), j, NA)
    for (j in 1:ncol(test )) set(test,  which(is.infinite(test[[j]])),  j, NA)
    
    for (j in names(train)) {
      set(train,which(is.na(train[[j]])),j,mean(train[[j]],na.rm=T))
    }
    for (j in names(test)) {
      set(test,which(is.na(test[[j]])),j,mean(test[[j]],na.rm=T))
    }
    
    imbalance = round(nrow(train[target==0])/nrow(train[target==1]))
    new_rows = do.call("rbind", replicate((imbalance-1), train[target==1], simplify = FALSE))
    train_imbalancefixed = rbind(train,new_rows)
    
    X1_fixed = as.matrix(train_imbalancefixed[,-ncol(train), with=FALSE])
    X1 = as.matrix(train[,-ncol(train), with=FALSE])
    X2 = as.matrix(test)
    Y1 = as.factor(train_imbalancefixed$target)
    
    #rf_model = randomForest(x=X1, y=Y1, ntree=100)  
    
    # Algorithm Tune (tuneRF)
    rf_model <- tuneRF(X1_fixed, Y1, stepFactor=1.5, improve=1e-5, ntree=500, doBest=TRUE)
    #print(bestmtry)
    
    train_prediction = predict(rf_model, X1)
    test_prediction  = predict(rf_model, X2)
    
    save_prediction(method='randomforest', test$id, test_prediction, train_prediction)
    return(rf_model)
  }
  ##################################################################################
  
  ##################################################################################
  standalone_models <- function(dt_train, dt_test, model_name, suffix='0') {
    path = getwd()
    patient = as.numeric(unlist(strsplit(path,'_'))[2])
    
    filename = paste0('new_',patient,'_',dt_test$id,'.mat')
    
    train = copy(dt_train)
    test  = copy(dt_test)
    #new_rows = do.call("rbind", replicate((imbalance-1), train[target==1], simplify = FALSE))
    #train=rbind(train,new_rows)
    #   Suggestions:
    #   eta = c(0.01,0.2)
    #   min_child_weight=1
    #   max_depth=c(3,10)
    #   gamma=c(0,1)
    #   #max_delta_step=
    #   subsample=c(0.5,1)
    #   colsample_bytree=c(0.5,1)
    #   lambda=1
    #   alpha=0
    #   scale_pos_weight=1
    for (j in 1:ncol(train)) set(train, which(is.infinite(train[[j]])), j, NA)
    for (j in 1:ncol(test )) set(test,  which(is.infinite(test[[j]])),  j, NA)
    
    for (j in names(train)) {
      set(train,which(is.na(train[[j]])),j,mean(train[[j]],na.rm=T))
    }
    for (j in names(test)) {
      set(test,which(is.na(test[[j]])),j,mean(test[[j]],na.rm=T))
    }
    
    train[target=='0',targetF:='N']
    train[target=='1',targetF:='S']
    train[,targetF:=as.factor(targetF)]
    train[,target:=NULL]
    
    # Create balanced training data.
    imbalance = round(nrow(train[targetF=='N'])/nrow(train[targetF=='S']))
    new_rows = do.call("rbind", replicate((imbalance-1), train[targetF=='S'], simplify = FALSE))
    train_balanced = rbind(train,new_rows)
    train_balanced = train_balanced[sample(nrow(train),replace = F)]
    
    trcontrol <- trainControl(method="cv",
                              number = 5,
                              verboseIter  = TRUE,
                              returnData   = FALSE,
                              returnResamp = "all",
                              allowParallel = FALSE,
                              summaryFunction=twoClassSummary,
                              classProbs=T)
    
    if ( model_name=='xgboost' ) {
      cat('   Training xgboost\n')
      
      #xgb_grid_1 <- expand.grid(nrounds= 100,
      #                          max_depth=c(6),
      #                          gamma=c(0),
      #                          min_child_weight=c(0.5),
      #                          colsample_bytree=c(0.5),
      #                          eta=c(0.01))
      
      xgb_grid_1 <- expand.grid(nrounds= 10000,
                                max_depth=c(6,10,15),
                                gamma=c(0,1),
                                min_child_weight=c(0.1,0.5),
                                colsample_bytree=c(0.1,0.5),
                                eta=c(0.01))
      
      model_imbalanced <- train(targetF~.,
                                data=train,
                                method="xgbTree",
                                trControl = trcontrol,
                                tuneGrid  = xgb_grid_1,
                                metric="ROC")
      
      model_balanced <- train(targetF~.,
                              data=train_balanced,
                              method="xgbTree",
                              trControl = trcontrol,
                              tuneGrid  = xgb_grid_1,
                              metric="ROC")
      
      xgb_balanced   = model_balanced$finalModel
      xgb_imbalanced = model_imbalanced$finalModel
      
      importance_balanced   <- xgb.importance(model = xgb_balanced)
      importance_imbalanced <- xgb.importance(model = xgb_imbalanced)
      
      write.csv(importance_balanced,   file='models/xgb_importance_balanced.csv')
      write.csv(importance_imbalanced, file='models/xgb_importance_imbalanced.csv')
    }
    
    if ( model=='rf' ) {
      cat('   Training RF...\n')
      rf_grid_1 <- expand.grid(mtry=c(5,10,20))
      #rf_grid_1 <- expand.grid(mtry=c(5))
      
      model_balanced <- train(targetF ~., 
                              data=train_balanced, 
                              method="rf", 
                              ntree=1000,
                              trControl = trcontrol,
                              tuneGrid  = rf_grid_1,
                              metric="ROC")
      
      model_imbalanced <- train(targetF ~., 
                        data=train, 
                        method="rf", 
                        ntree=1000,
                        trControl = trcontrol,
                        tuneGrid  = rf_grid_1,
                        metric="ROC")
      
      rf_imbalanced = model_imbalanced$finalModel
      rf_balanced   = model_balanced$finalModel
      
      importance_imbalanced = importance(rf_imbalanced)
      importance_balanced   = importance(rf_balanced)
      
      write.csv(importance_balanced,   file='models/rf_importance_balanced.csv')
      write.csv(importance_imbalanced, file='models/rf_importance_imbalanced.csv')
      
    }
    
    if ( model_name=='nnet') {
      cat('   Training nnet...\n')
      #nnet_grid_1 <- expand.grid(size=2,decay=0.5)
      nnet_grid_1 <- expand.grid(size=c(2,5,7),decay=c(0.5,0.1))
      
      model_imbalanced <- train(targetF~., 
                          data=train, 
                          method="nnet", 
                          trControl = trcontrol,
                          tuneGrid  = nnet_grid_1,
                          preProcess = c("center", "scale"),
                          metric="ROC")
      
      model_balanced <- train(targetF~., 
                                data=train_balanced, 
                                method="nnet", 
                                trControl = trcontrol,
                                tuneGrid  = nnet_grid_1,
                                preProcess = c("center", "scale"),
                                metric="ROC")
      
      nnet_imbalanced = model_imbalanced$finalModel
      nnet_balanced   = model_balanced$finalModel
      
      importance_imbalanced = garson(nnet_imbalanced, bar_plot=F)
      importance_balanced   = garson(nnet_balanced, bar_plot=F)
      
      write.csv(importance_balanced,   file='models/nnet_importance_balanced.csv')
      write.csv(importance_imbalanced, file='models/nnet_importance_imbalanced.csv')
    }
    #png(filename = 'models/xgb_importance_matrix.png',height=100,width=20,units='cm',res=900)
    #xgb.plot.importance(importance_matrix = importance_matrix)
    #dev.off()
    
    ###########################################
    ###### Build stand-alone submissions.######
    ###########################################
    # Probabilities.
    prediction_balanced   = predict(model_balanced,   newdata=test, type='prob')
    prediction_imbalanced = predict(model_imbalanced, newdata=test, type='prob')
    sub_balanced   <- data.table(File=filename, Class=prediction_balanced$S)
    sub_imbalanced <- data.table(File=filename, Class=prediction_imbalanced$S)
    write.csv(sub_balanced,   file=sprintf('models/standalone_%s_imbalanced_probs_%s.csv',model_name,suffix), row.names=F, quote=F)
    write.csv(sub_imbalanced, file=sprintf('models/standalone_%s_balanced_probs_%s.csv',model_name,suffix),   row.names=F, quote=F)
    
    # Raw.
    prediction_balanced   = predict(model_balanced,   newdata=test, type='raw')
    prediction_imbalanced = predict(model_imbalanced, newdata=test, type='raw')    
    sub_balanced   <- data.table(File=filename, Class2=prediction_balanced)
    sub_imbalanced <- data.table(File=filename, Class2=prediction_imbalanced)
    sub_balanced[Class2=='N',Class:=0]
    sub_imbalanced[Class2=='N',Class:=0]
    sub_balanced[Class2=='S',Class:=1]
    sub_imbalanced[Class2=='S',Class:=1]
    sub_balanced[,Class2:=NULL]
    sub_imbalanced[,Class2:=NULL]
    write.csv(sub_balanced,   file=sprintf('models/standalone_%s_balanced_raw_%s.csv',model_name,suffix),   row.names=F, quote=F)
    write.csv(sub_imbalanced, file=sprintf('models/standalone_%s_imbalanced_raw_%s.csv',model_name,suffix), row.names=F, quote=F)
    
    #return(xgb_model)
  }
  ##################################################################################
  
  
  find_bad_features <- function(feature_names) {
    feature_names = feature_names[!(feature_names%in%c('id'))]
    
    f1 = fread('models/nnet_importance_balanced.csv')
    f2 = fread('models/nnet_importance_imbalanced.csv')
    f3 = fread('models/rf_importance_balanced.csv')
    f4 = fread('models/rf_importance_imbalanced.csv')
    f5 = fread('models/xgb_importance_balanced.csv')
    f5[,f:=as.integer(Feature)]
    f6 = fread('models/xgb_importance_imbalanced.csv')
    f6[,f:=as.integer(Feature)]
    
    map = data.table(f=seq(1,length(feature_names)), V1=feature_names)
    f5 = merge(f5, map, by='f')
    f6 = merge(f6, map, by='f')
    f5=f5[,list(V1.y,Gain)]
    f6=f6[,list(V1.y,Gain)]
    
    setorder(f1,-rel_imp)
    setorder(f2,-rel_imp)
    setorder(f3,-MeanDecreaseGini)
    setorder(f4,-MeanDecreaseGini)
    setorder(f5,-Gain)
    setorder(f6,-Gain)
    
    frac = 0.3
    N2=392
    N1=round(N2*frac)
    s1 = f1[N1:N2,V1]
    s2 = f2[N1:N2,V1]
    s3 = f3[N1:N2,V1]
    s4 = f4[N1:N2,V1]
    s5 = f5[round(365*frac):365,V1.y]
    s6 = f6[round(365*frac):365,V1.y]
    
    gsub("[^[:alnum:][:blank:]+?&/_]", "", s1)
    gsub("[^[:alnum:][:blank:]+?&/_]", "", s2)
    gsub("[^[:alnum:][:blank:]+?&/_]", "", s3)
    gsub("[^[:alnum:][:blank:]+?&/_]", "", s4)
    gsub("[^[:alnum:][:blank:]+?&/_]", "", s5)
    gsub("[^[:alnum:][:blank:]+?&/_]", "", s6)
    
    bad_features = Reduce(intersect, list(s1,s2,s3,s4,s5,s6))

    return(bad_features)
}
  
  
  