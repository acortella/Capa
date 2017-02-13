## Initisation du chemin du dossier contenant les données
path <- "C:/Users/A/Desktop/CP_data/"
setwd(path)

## Installation des librairies si besoin
if(!isTRUE(require(data.table))) install.packages("data.table")
if(!isTRUE(require(caTools))) install.packages("caTools")
if(!isTRUE(require(plotly))) install.packages("plotly")
if(!isTRUE(require(e1071))) install.packages("e1071")
if(!isTRUE(require(caret))) install.packages("caret")
if(!isTRUE(require(doSNOW))) install.packages("doSNOW")

## Chargement des librairies
library(data.table)
library(caTools)
library(plotly)
library(e1071)
library(caret)
library(doSNOW)

load("data_processing.RData")

####################################################################################################################
###################################### Classification ##############################################################
####################################################################################################################

## Divise les volumétries utilisées du training set et du test set par la dernière volumétrie limite
div <- function(my_data){
    my_data <- cbind(my_data[,"ind", with=F], 
                          as.matrix(my_data[,2:95,with=F])/as.vector(as.matrix(my_data[,92,with=F])))
    my_data <- my_data[,':='(
        "y10"=apply(my_data[,"y10",with=F],1,function(x){ifelse(unlist(x)>=1, "YES", "NO")}),
        "y30"=apply(my_data[,"y30",with=F],1,function(x){ifelse(unlist(x)>=1, "YES", "NO")}),
        "y60"=apply(my_data[,"y60",with=F],1,function(x){ifelse(unlist(x)>=1, "YES", "NO")}))][,-c("limit"),with=F]
    return(my_data)
}
training_set <- div(training_set)
test_set <- div(test_set)


## Vrai positif dans le test set
## Cas 10 - 72
P_10_test <- nrow(test_set[y10=="YES"])
## Cas 30 - 225
P_30_test <- nrow(test_set[y30=="YES"])
## Cas 60 - 580
P_60_test <- nrow(test_set[y60=="YES"])


### Filter training set
filtre_coefficient <- function(my_vector, my_case=10, my_limit=1){
    my_vector <- unlist(my_vector)
    if(length(my_vector[!is.na(my_vector)>2]) && !is.na(my_vector[length(my_vector)])){
        l_last_value <- my_vector[length(my_vector)]
        my_vector <- diff(my_vector, 1)[!is.na(my_vector)]
        ifelse(any(lapply(my_vector, function(x){x*my_case+l_last_value}))>=my_limit, return(1), return(0))
    }else{
        return(0)
    }
}

## Filter max limit
filtre_max_limit <- function(my_vector, my_limit=1){
    ifelse(any(unlist(my_vector)>=my_limit*0.90), 1, 0)
}

## Filtre du training_set
training_set <- training_set[,":="("fc"=apply(training_set[,2:91,with=F],1,function(x){filtre_coefficient(x, 10, 1)}),
                                       "fm"=apply(training_set[,2:91,with=F],1,filtre_max_limit))][fc==1][fm==1][,-c("fc", "fm"),with=F]


### Feature engineering
get_features <- function(my_training=training_class, my_test=test_class, my_case){
  if(my_case==10){
      my_training <- my_training[,':='("response"=y10)][,-c("y10", "y30", "y60"),with=F]
      my_test <- my_test[,':='("response"=y10)][,-c("y10", "y30", "y60"),with=F]
  }else if(my_case==30){
      my_training <- my_training[,':='("response"=y30)][,-c("y10", "y30", "y60"),with=F]
      my_test <- my_test[,':='("response"=y30)][,-c("y10", "y30", "y60"),with=F]
  }else if(my_case==60){
      my_training <- my_training[,':='("response"=y60)][,-c("y10", "y30", "y60"),with=F]
      my_test <- my_test[,':='("response"=y60)][,-c("y10", "y30", "y60"),with=F]
  }
  
  ### Génération de features mathématique simple pour chaque dimension d%%10 de 10 à 90
  my_training <- my_training[,':='("max_90"=(apply(my_training[,2:91,with=F],1,max)),
                                   "min_90"=(apply(my_training[,2:91,with=F],1,min)),
                                   "median_90"=(apply(my_training[,2:91,with=F],1,median)),
                                   "mean_90"=(apply(my_training[,2:91,with=F],1,mean)),
                                   "sd_90"=(apply(my_training[,2:91,with=F],1,sd)),
                                   
                                   "max_80"=(apply(my_training[,12:91,with=F],1,max)),
                                   "min_80"=(apply(my_training[,12:91,with=F],1,min)),
                                   "median_80"=(apply(my_training[,12:91,with=F],1,median)),
                                   "mean_80"=(apply(my_training[,12:91,with=F],1,mean)),
                                   "sd_80"=(apply(my_training[,12:91,with=F],1,sd)),
                                   
                                   "max_70"=(apply(my_training[,22:91,with=F],1,max)),
                                   "min_70"=(apply(my_training[,22:91,with=F],1,min)),
                                   "median_70"=(apply(my_training[,22:91,with=F],1,median)),
                                   "mean_70"=(apply(my_training[,22:91,with=F],1,mean)),
                                   "sd_70"=(apply(my_training[,22:91,with=F],1,sd)),
                                   
                                   "max_60"=(apply(my_training[,32:91,with=F],1,max)),
                                   "min_60"=(apply(my_training[,32:91,with=F],1,min)),
                                   "median_60"=(apply(my_training[,32:91,with=F],1,median)),
                                   "mean_60"=(apply(my_training[,32:91,with=F],1,mean)),
                                   "sd_60"=(apply(my_training[,32:91,with=F],1,sd)),
                                   
                                   "max_50"=(apply(my_training[,42:91,with=F],1,max)),
                                   "min_50"=(apply(my_training[,42:91,with=F],1,min)),
                                   "median_50"=(apply(my_training[,42:91,with=F],1,median)),
                                   "mean_50"=(apply(my_training[,42:91,with=F],1,mean)),
                                   "sd_50"=(apply(my_training[,42:91,with=F],1,sd)),
                                   
                                   "max_40"=(apply(my_training[,52:91,with=F],1,max)),
                                   "min_40"=(apply(my_training[,52:91,with=F],1,min)),
                                   "median_40"=(apply(my_training[,52:91,with=F],1,median)),
                                   "mean_40"=(apply(my_training[,52:91,with=F],1,mean)),
                                   "sd_40"=(apply(my_training[,52:91,with=F],1,sd)),
                                   
                                   "max_30"=(apply(my_training[,62:91,with=F],1,max)),
                                   "min_30"=(apply(my_training[,62:91,with=F],1,min)),
                                   "median_30"=(apply(my_training[,62:91,with=F],1,median)),
                                   "mean_30"=(apply(my_training[,62:91,with=F],1,mean)),
                                   "sd_30"=(apply(my_training[,62:91,with=F],1,sd)),
                                   
                                   "max_20"=(apply(my_training[,72:91,with=F],1,max)),
                                   "min_20"=(apply(my_training[,72:91,with=F],1,min)),
                                   "median_20"=(apply(my_training[,72:91,with=F],1,median)),
                                   "mean_20"=(apply(my_training[,72:91,with=F],1,mean)),
                                   "sd_20"=(apply(my_training[,72:91,with=F],1,sd)),
                                   
                                   "max_10"=(apply(my_training[,82:91,with=F],1,max)),
                                   "min_10"=(apply(my_training[,82:91,with=F],1,min)),
                                   "median_10"=(apply(my_training[,82:91,with=F],1,median)),
                                   "mean_10"=(apply(my_training[,82:91,with=F],1,mean)),
                                   "sd_10"=(apply(my_training[,82:91,with=F],1,sd)))]
  
  
    #### Supprime les colonnes qui sont trop corrélées
    correlationMatrix <- cor(my_training[,-c("ind", "response"),with=F])
    highlyCorrelated <- sort(findCorrelation(correlationMatrix, cutoff=0.9))
    temp_my_training <- my_training[,-c("ind", "response"),with=F]
    while(length(highlyCorrelated)>0){
        temp_my_training <- temp_my_training[,-c(highlyCorrelated[1]),with=F]
        correlationMatrix <- cor(temp_my_training)
        highlyCorrelated <- sort(findCorrelation(correlationMatrix, cutoff=0.9))
        print(highlyCorrelated)
        if(length(highlyCorrelated)>=2){
            highlyCorrelated <- highlyCorrelated[2:length(highlyCorrelated)]
        }else{
            break
        }
    }
    my_training <- cbind(my_training[,c("ind"),with=F], temp_my_training, my_training[,c("response"),with=F])
  
    #### Feature selection avec les forêts aléatoires
    ## Créer une liste de "seeds" pour la reproduction des résultats
    set.seed(1)
    #length is = (n_repeats*nresampling)+1
    seeds <- vector(mode = "list", length = 11)
    #(3 is the number of tuning parameter, mtry for rf, here equal to ncol(iris)-2)
    for(i in 1:10) seeds[[i]]<- sample.int(n=1000, ncol(my_training[,-c("ind", "response"),with=F]))
    #for the last model
    seeds[[11]]<-sample.int(1000, 1)
    
    ## Transform en facteur pour le modèle
    my_training <- my_training[,response:=as.factor(response)]
    
    ## Prend un échantillon du set d'entrainement (tous les yes + échantillon de 385 non)
    nb_random <- nrow(my_training[response!="YES"])-nrow(my_training[response=="YES"])
    nb_random <- sample(nb_random, 385, replace=F)
    temp_table <- cbind(my_training[response!="YES"], data.table("random"=1:nrow(my_training[response!="YES"])))
    my_training <- rbindlist(list(temp_table[random %in% nb_random][,-c("random"),with=F], my_training[response=="YES"]))
    my_training <- as.data.table(as.data.frame(my_training)[sample(nrow(my_training)),])
    temp_train <- as.data.frame(my_training[,-c("ind", "response"),with=F])
    temp_test <- unlist((my_training[,"response",with=F]))
    
    preProcValues <- preProcess(temp_train, method = c("center", "scale"))
    trainTransformed <- predict(preProcValues, temp_train)
    
    ## Feature selection algorithm
    control <- rfeControl(functions=rfFuncs, method="cv", number=10, seeds=seeds)
    cl <- makeCluster(3)
    registerDoSNOW(cl)
  
    
    results <- rfe(temp_train,  temp_test, sizes=c(1:68), rfeControl=control)
    registerDoSNOW(cl)
    stopCluster(cl)
    
    col_keep <- names(results$fit$forest$ncat)[1:5]
    
    return(col_keep)
}



best_score <- function(data, mu=10^(-10), deltat=0.01, alpha=2, beta=1, seuil=7){
    
}





feature_10 <- get_features(my_training=training_class, my_test=test_class, 10)
feature_30 <- get_features(my_training=training_class, my_test=test_class, 30)
feature_60 <- get_features(my_training=training_class, my_test=test_class, 60)


### Machine learning
ml_cl <- function(my_training=training_class, my_test=test_class, my_case, my_feature, my_algo){
  my_feature <- feature_10
  my_training <- training_class
  my_test <- test_class
  my_case <- 10
  if(my_case==10){
    my_training <- my_training[,':='("response"=y10)][,-c("y10", "y30", "y60"),with=F]
    my_test <- my_test[,':='("response"=y10)][,-c("y10", "y30", "y60"),with=F]
  }else if(my_case==30){
    my_training <- my_training[,':='("response"=y30)][,-c("y10", "y30", "y60"),with=F]
    my_test <- my_test[,':='("response"=y30)][,-c("y10", "y30", "y60"),with=F]
  }else if(my_case==60){
    my_training <- my_training[,':='("response"=y60)][,-c("y10", "y30", "y60"),with=F]
    my_test <- my_test[,':='("response"=y60)][,-c("y10", "y30", "y60"),with=F]
  }
  
  my_training <- my_training[,':='("max_90"=(apply(my_training[,2:91,with=F],1,max)),
                                   "min_90"=(apply(my_training[,2:91,with=F],1,min)),
                                   "median_90"=(apply(my_training[,2:91,with=F],1,median)),
                                   "mean_90"=(apply(my_training[,2:91,with=F],1,mean)),
                                   "sd_90"=(apply(my_training[,2:91,with=F],1,sd)),
                                   
                                   "max_80"=(apply(my_training[,12:91,with=F],1,max)),
                                   "min_80"=(apply(my_training[,12:91,with=F],1,min)),
                                   "median_80"=(apply(my_training[,12:91,with=F],1,median)),
                                   "mean_80"=(apply(my_training[,12:91,with=F],1,mean)),
                                   "sd_80"=(apply(my_training[,12:91,with=F],1,sd)),
                                   
                                   "max_70"=(apply(my_training[,22:91,with=F],1,max)),
                                   "min_70"=(apply(my_training[,22:91,with=F],1,min)),
                                   "median_70"=(apply(my_training[,22:91,with=F],1,median)),
                                   "mean_70"=(apply(my_training[,22:91,with=F],1,mean)),
                                   "sd_70"=(apply(my_training[,22:91,with=F],1,sd)),
                                   
                                   "max_60"=(apply(my_training[,32:91,with=F],1,max)),
                                   "min_60"=(apply(my_training[,32:91,with=F],1,min)),
                                   "median_60"=(apply(my_training[,32:91,with=F],1,median)),
                                   "mean_60"=(apply(my_training[,32:91,with=F],1,mean)),
                                   "sd_60"=(apply(my_training[,32:91,with=F],1,sd)),
                                   
                                   "max_50"=(apply(my_training[,42:91,with=F],1,max)),
                                   "min_50"=(apply(my_training[,42:91,with=F],1,min)),
                                   "median_50"=(apply(my_training[,42:91,with=F],1,median)),
                                   "mean_50"=(apply(my_training[,42:91,with=F],1,mean)),
                                   "sd_50"=(apply(my_training[,42:91,with=F],1,sd)),
                                   
                                   "max_40"=(apply(my_training[,52:91,with=F],1,max)),
                                   "min_40"=(apply(my_training[,52:91,with=F],1,min)),
                                   "median_40"=(apply(my_training[,52:91,with=F],1,median)),
                                   "mean_40"=(apply(my_training[,52:91,with=F],1,mean)),
                                   "sd_40"=(apply(my_training[,52:91,with=F],1,sd)),
                                   
                                   "max_30"=(apply(my_training[,62:91,with=F],1,max)),
                                   "min_30"=(apply(my_training[,62:91,with=F],1,min)),
                                   "median_30"=(apply(my_training[,62:91,with=F],1,median)),
                                   "mean_30"=(apply(my_training[,62:91,with=F],1,mean)),
                                   "sd_30"=(apply(my_training[,62:91,with=F],1,sd)),
                                   
                                   "max_20"=(apply(my_training[,72:91,with=F],1,max)),
                                   "min_20"=(apply(my_training[,72:91,with=F],1,min)),
                                   "median_20"=(apply(my_training[,72:91,with=F],1,median)),
                                   "mean_20"=(apply(my_training[,72:91,with=F],1,mean)),
                                   "sd_20"=(apply(my_training[,72:91,with=F],1,sd)),
                                   
                                   "max_10"=(apply(my_training[,82:91,with=F],1,max)),
                                   "min_10"=(apply(my_training[,82:91,with=F],1,min)),
                                   "median_10"=(apply(my_training[,82:91,with=F],1,median)),
                                   "mean_10"=(apply(my_training[,82:91,with=F],1,mean)),
                                   "sd_10"=(apply(my_training[,82:91,with=F],1,sd)))]
  
  
  my_test <- my_test[,':='("max_90"=(apply(my_test[,2:91,with=F],1,max)),
                           "min_90"=(apply(my_test[,2:91,with=F],1,min)),
                           "median_90"=(apply(my_test[,2:91,with=F],1,median)),
                           "mean_90"=(apply(my_test[,2:91,with=F],1,mean)),
                           "sd_90"=(apply(my_test[,2:91,with=F],1,sd)),
                           
                           "max_80"=(apply(my_test[,12:91,with=F],1,max)),
                           "min_80"=(apply(my_test[,12:91,with=F],1,min)),
                           "median_80"=(apply(my_test[,12:91,with=F],1,median)),
                           "mean_80"=(apply(my_test[,12:91,with=F],1,mean)),
                           "sd_80"=(apply(my_test[,12:91,with=F],1,sd)),
                           
                           "max_70"=(apply(my_test[,22:91,with=F],1,max)),
                           "min_70"=(apply(my_test[,22:91,with=F],1,min)),
                           "median_70"=(apply(my_test[,22:91,with=F],1,median)),
                           "mean_70"=(apply(my_test[,22:91,with=F],1,mean)),
                           "sd_70"=(apply(my_test[,22:91,with=F],1,sd)),
                           
                           "max_60"=(apply(my_test[,32:91,with=F],1,max)),
                           "min_60"=(apply(my_test[,32:91,with=F],1,min)),
                           "median_60"=(apply(my_test[,32:91,with=F],1,median)),
                           "mean_60"=(apply(my_test[,32:91,with=F],1,mean)),
                           "sd_60"=(apply(my_test[,32:91,with=F],1,sd)),
                           
                           "max_50"=(apply(my_test[,42:91,with=F],1,max)),
                           "min_50"=(apply(my_test[,42:91,with=F],1,min)),
                           "median_50"=(apply(my_test[,42:91,with=F],1,median)),
                           "mean_50"=(apply(my_test[,42:91,with=F],1,mean)),
                           "sd_50"=(apply(my_test[,42:91,with=F],1,sd)),
                           
                           "max_40"=(apply(my_test[,52:91,with=F],1,max)),
                           "min_40"=(apply(my_test[,52:91,with=F],1,min)),
                           "median_40"=(apply(my_test[,52:91,with=F],1,median)),
                           "mean_40"=(apply(my_test[,52:91,with=F],1,mean)),
                           "sd_40"=(apply(my_test[,52:91,with=F],1,sd)),
                           
                           "max_30"=(apply(my_test[,62:91,with=F],1,max)),
                           "min_30"=(apply(my_test[,62:91,with=F],1,min)),
                           "median_30"=(apply(my_test[,62:91,with=F],1,median)),
                           "mean_30"=(apply(my_test[,62:91,with=F],1,mean)),
                           "sd_30"=(apply(my_test[,62:91,with=F],1,sd)),
                           
                           "max_20"=(apply(my_test[,72:91,with=F],1,max)),
                           "min_20"=(apply(my_test[,72:91,with=F],1,min)),
                           "median_20"=(apply(my_test[,72:91,with=F],1,median)),
                           "mean_20"=(apply(my_test[,72:91,with=F],1,mean)),
                           "sd_20"=(apply(my_test[,72:91,with=F],1,sd)),
                           
                           "max_10"=(apply(my_test[,82:91,with=F],1,max)),
                           "min_10"=(apply(my_test[,82:91,with=F],1,min)),
                           "median_10"=(apply(my_test[,82:91,with=F],1,median)),
                           "mean_10"=(apply(my_test[,82:91,with=F],1,mean)),
                           "sd_10"=(apply(my_test[,82:91,with=F],1,sd)))]
  
  my_training <- my_training[,c("ind", my_feature, "response"),with=F][,response:=as.factor(response)]
  my_test <- my_test[,c("ind", my_feature, "response"),with=F][,response:=as.factor(response)]
  
  ## Prend un échantillon du set d'entrainement (tous les yes + échantillon de 385 non)
  nb_random <- nrow(my_training[response!="YES"])-nrow(my_training[response=="YES"])
  nb_random <- sample(nb_random, 385, replace=F)
  temp_table <- cbind(my_training[response!="YES"], data.table("random"=1:nrow(my_training[response!="YES"])))
  my_training <- rbindlist(list(temp_table[random %in% nb_random][,-c("random"),with=F], my_training[response=="YES"]))
  my_training <- as.data.table(as.data.frame(my_training)[sample(nrow(my_training)),])
  temp_train <- as.data.frame(my_training[,-c("ind", "response"),with=F])
  temp_test <- unlist((my_training[,"response",with=F]))
  
  set.seed(1)
  #length is = (n_repeats*nresampling)+1
  seeds <- vector(mode = "list", length = 11)
  #(3 is the number of tuning parameter, mtry for rf, here equal to ncol(iris)-2)
  for(i in 1:10) seeds[[i]]<- sample.int(n=1000, ncol(my_training[,-c("ind", "response"),with=F]))
  #for the last model
  seeds[[11]]<-sample.int(1000, 1)
  
  control <- rfeControl(functions=rfFuncs, method="cv", number=10, seeds=seeds)
  my_model <- caret::train()
  
  system.time({results <- rfe(temp_train,  temp_test, sizes=c(1:68), rfeControl=control)})
  col_keep_1 <- names(results$fit$forest$ncat)
  
  
  
  
  
  
  #####################
  # training_class_2 <- training_class
  # ### 
  # training_class_2 <- as.data.table(as.data.frame(training_class_2)[sample(nrow(training_class_2)),])
  # temp_train <- as.data.frame(training_class_2[,-c("ind", "response"),with=F])
  # temp_test <- as.factor(unlist((training_class_2[,"response",with=F])))
  # system.time({results <- rfe(temp_train,  temp_test, sizes=c(1:ncol(temp_train)), rfeControl=control)})
  # col_keep_2 <- names(results$fit$forest$ncat)
  
  
  preProcValues <- preProcess(training_class, method = c("center", "scale"))
  
  trainTransformed <- predict(preProcValues, training_class)
  testTransformed <- predict(preProcValues, test_call)
  
  
  control <- rfeControl(functions=rfFuncs, method="cv", number=10, seeds=seeds, index=createFolds(training_class$y10))
  
  library(doSNOW)
  
  cl <- makeCluster(detectCores()-1)
  registerDoSNOW(cl)
  model1 <- train(Species~., iris, method='rf', trControl=myControl)
  stopCluster(cl)
  
  predict(model1, testTransformed, type='prob')
  
  
  
  training_class_2 <- training_class[, col_keep, with=F]
  
  correlationMatrix <- cor(training_class_2[,1:8])
  
  
  # summarize the correlation matrix
  print(correlationMatrix)
  # find attributes that are highly corrected (ideally >0.75)
  highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.5)
  # print indexes of highly correlated attributes
  print(highlyCorrelated)
  
  na.omit(temp_train)
  nrow(temp_train)
  nrow(temp_test)
  
  ?rfe
  unique(names(temp_train))
  names(temp_test)
  
  library(mlbench)
  
  
  
  data(PimaIndiansDiabetes)
  # define the control using a random forest selection function
  control <- rfeControl(functions=rfFuncs, method="cv", number=10)
  # run the RFE algorithm
  results <- rfe(PimaIndiansDiabetes[,1:8], PimaIndiansDiabetes[,9], sizes=c(1:8), rfeControl=control)
  
  ?quantile   
}






feature_coefficient <- function(my_vector, my_case, my_limit=1){
  my_vector <- unlist(my_vector)
  if(length(my_vector[!is.na(my_vector)>2]) && !is.na(my_vector[length(my_vector)])){
    l_last_value <- my_vector[lenght(my_vector)]
    my_vector <- diff(my_vector, 1)[!is.na(my_vector)]
    my_vector <- lapply(my_vector, function(x){x*my_case+l_last_value})
    return(length(my_vector[which(my_vector>1)]))
  }else{
    return(0)
  }
}


feature_max_limit <- function(my_vector){
  my_vector <- unlist(my_vector)
  my_vector <- length(my_vector[which(my_vector>=1)])
}





