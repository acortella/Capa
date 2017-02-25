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
if(!isTRUE(require(glmnet))) install.packages("glmnet")
if(!isTRUE(require(party))) install.packages("party")
if(!isTRUE(require(mboost))) install.packages("mboost")
if(!isTRUE(require(plyr))) install.packages("plyr")
if(!isTRUE(require(obliqueRF))) install.packages("obliqueRF")
if(!isTRUE(require(ranger))) install.packages("ranger")
if(!isTRUE(require(rotationForest))) install.packages("rotationForest")
if(!isTRUE(require(RSNNS))) install.packages("RSNNS")
if(!isTRUE(require(caTools))) install.packages("caTools")
if(!isTRUE(require(kernlab))) install.packages("kernlab")
if(!isTRUE(require(arm))) install.packages("arm")


## Chargement des librairies
library(data.table)
library(caTools)
library(plotly)
library(e1071)
library(caret)
library(doSNOW)
library(glmnet)
library(party)
library(mboost)
library(plyr)
library(obliqueRF)
library(ranger)
library(rotationForest)
library(RSNNS)
library(caTools)
library(kernlab)
library(arm)

load("data_processing.RData")

####################################################################################################################
###################################### Regression ##############################################################
####################################################################################################################


## Vrai positif dans le test set
## Cas 10 - 72
P_10_test <- nrow(test_set[y10=="YES"])
## Cas 30 - 225
P_30_test <- nrow(test_set[y30=="YES"])
## Cas 60 - 580
P_60_test <- nrow(test_set[y60=="YES"])


ind <- training_set[,ind]
test_set <- test_set[,':='(
    "y10"=apply(test_set[,"y10",with=F],1,function(x){ifelse(unlist(x)>=0, "YES", "NO")}),
    "y30"=apply(test_set[,"y30",with=F],1,function(x){ifelse(unlist(x)>=0, "YES", "NO")}),
    "y60"=apply(test_set[,"y60",with=F],1,function(x){ifelse(unlist(x)>=0, "YES", "NO")}))]






### Filter training set
filtre_coefficient <- function(my_vector, my_case=60){
    my_vector <- unlist(my_vector)
    my_limit <-  my_vector[length(my_vector)]
    my_vector <- my_vector[1:(length(my_vector)-1)]
    l_last_value <- my_vector[length(my_vector)]
    my_vector <- diff(my_vector, 1)
    my_vector <- my_vector[!is.na(my_vector)]
    ifelse(any(unlist(lapply(my_vector, function(x){x*my_case+l_last_value}))>=my_limit), return(1), return(0))
}

## Filter max limit
filtre_max_limit <- function(my_vector){
    my_vector <- unlist(my_vector)
    my_limit <-  my_vector[length(my_vector)]
    my_vector <- my_vector[1:(length(my_vector)-1)]
    ifelse(any(unlist(my_vector)>=my_limit*0.90), 1, 0)
}

## filtre training_set
filtre_training <- function(my_training=training_set){
    vect <- sort(unique(my_training[,ind]))
    filtered <- NULL
    for(g in 1:length(vect)){
        temp_train <- my_training[ind==vect[g]]
        temp_train <- temp_train[,":="("fc"=apply(temp_train[,2:92,with=F],1,function(x){filtre_coefficient(x, 60)}),
                                           "fm"=apply(temp_train[,2:92,with=F],1,filtre_max_limit))]
        if(all(temp_train[,fc]==0) && all(temp_train[,fm]==0)){
            filtered <- c(filtered, vect[g])
        }
        print(c(g, length(filtered)))
    }
    my_training <- my_training[!ind %in% filtered]
    return(my_training)
}
training_set <- filtre_training(training_set)


## Soustrait les max du training set et du test set par la dernière volumétrie limite
sous <- function(my_data){
    my_data <- cbind(my_data[,1:91, with=F],  as.matrix(my_data[,93:95,with=F])-as.vector(as.matrix(my_data[,92,with=F])))
    return(my_data)
}
training_set <- sous(training_set)
test_set <- sous(test_set)

## Optimise le score d'après la métrique de score
best_score <- function(pred, obs, th=NULL, mu=10^(-10), delta=0.01, alpha=2, beta=1, stand=7){
    result <- data.table("obs"=obs, "pred"=pred)
    if(!is.null(th)){
        result <- result[,':='("pred"=apply(result[,"pred",with=F],1,function(x){ifelse(x>=th,1,0)}))]
    }
    result <- result [,':='("detail"=apply(result[,c("obs", "pred"),with=F],1,get_confusion))]
    TP <- nrow(result[detail=="TP"])
    FP <- nrow(result[detail=="FP"])
    FN <- nrow(result[detail=="FN"])
    
    qual <- ifelse(is.na(TP/(TP+FP)), 0, TP/(TP+FP))
    quant <- ifelse(is.na(TP/(TP+FN)), 0, TP/(TP+FN))
    if(qual>=(stand/(stand+1))){
        score <- ((alpha*qual)+(beta*quant))/(alpha+beta)
    }else{
        score <- (0.5*(qual^2))+(mu*quant)/((qual^2)+delta)
    }
    return(score)
}

## Donne la matrice de confusion
get_confusion <- function(x){
    x <- unlist(x)
    if(x[1]==1){
        if(x[2]==1){return("TP")
        }else if(x[2]==0){return("FN")}
    }else if(x[1]==0){
        if(x[2]==1){return("FP")
        }else if(x[2]==0){return("TN")}
    }
}

## Métrique définie (contient deux autres fonctions pour les intégrer dans les clusters de DoSnow)
quaqua_Summary_reg <- function (data, lev = NULL, model = NULL) {
    # saveRDS(data, "data.rds")
    # data <- readRDS("data.rds")
    
    obs <- unlist(data$obs)
    
    # if("YES" %in% obs){
    #     obs[obs=="YES"] <- 1
    # }
    # 
    # if("NO" %in% obs){
    #     obs[obs=="NO"] <- 0
    # }
    obs[obs>=0] <- 2
    obs[obs<0] <- 1
    # obs <- unlist(lapply(obs, function(x){ifelse(x>=0, 2, 1)}))
    
    if("pred"  %in% names(data)){
        pred <- unlist(data$pred)
        pred[is.na(pred)] <- -1
        # pred[x>=0] <- 1
        # pred[x<0] <- 0
        pred <- unlist(lapply(pred, function(x){ifelse(is.na(x), -1, x)}))
        pred <- unlist(lapply(pred, function(x){ifelse(x>=0, 1, 0)}))
        score <- pred + obs
        TP <- length(score[which(score==3)])
        FP <- length(score[which(score==1)])
        FN <- length(score[which(score==2)])
        mu <- 10^(-10)
        delta <- 0.01
        alpha <- 2
        beta <- 1 
        stand <- 7
        qual <- ifelse(is.na(TP/(TP+FP)), 0, TP/(TP+FP))
        quant <- ifelse(is.na(TP/(TP+FN)), 0, TP/(TP+FN))
        if(qual>=(stand/(stand+1))){
            score <- ((alpha*qual)+(beta*quant))/(alpha+beta)
        }else{
            score <- (0.5*(qual^2))+(mu*quant)/((qual^2)+delta)
        }
    }else{
        score <- NA
    }
    out <- score
    names(out) <- "qq_metric"
    out
}

#####################################################
#### #####################  FE #########
#####################################################
## Ajoute de nouvelles features
## Ajoute de nouvelles features
add_feature <- function(my_training, my_case){
    my_training <- my_training[,':='(
        
        
        "max_10"=(apply(my_training[,81:90,with=F],1,max)),
        "min_10"=(apply(my_training[,81:90,with=F],1,min)),
        "median_10"=(apply(my_training[,81:90,with=F],1,median)),
        "mean_10"=(apply(my_training[,81:90,with=F],1,mean)),
        "sd_10"=(apply(my_training[,81:90,with=F],1,sd)),
        
        "max_20"=(apply(my_training[,71:90,with=F],1,max)),
        "min_20"=(apply(my_training[,71:90,with=F],1,min)),
        "median_20"=(apply(my_training[,71:90,with=F],1,median)),
        "mean_20"=(apply(my_training[,71:90,with=F],1,mean)),
        "sd_20"=(apply(my_training[,71:90,with=F],1,sd)),
        
        "max_30"=(apply(my_training[,61:90,with=F],1,max)),
        "min_30"=(apply(my_training[,61:90,with=F],1,min)),
        "median_30"=(apply(my_training[,61:90,with=F],1,median)),
        "mean_30"=(apply(my_training[,61:90,with=F],1,mean)),
        "sd_30"=(apply(my_training[,61:90,with=F],1,sd)),
        
        "max_40"=(apply(my_training[,51:90,with=F],1,max)),
        "min_40"=(apply(my_training[,51:90,with=F],1,min)),
        "median_40"=(apply(my_training[,51:90,with=F],1,median)),
        "mean_40"=(apply(my_training[,51:90,with=F],1,mean)),
        "sd_40"=(apply(my_training[,51:90,with=F],1,sd)),
        
        "max_50"=(apply(my_training[,41:90,with=F],1,max)),
        "min_50"=(apply(my_training[,41:90,with=F],1,min)),
        "median_50"=(apply(my_training[,41:90,with=F],1,median)),
        "mean_50"=(apply(my_training[,41:90,with=F],1,mean)),
        "sd_50"=(apply(my_training[,41:90,with=F],1,sd)),
        
        "max_60"=(apply(my_training[,31:90,with=F],1,max)),
        "min_60"=(apply(my_training[,31:90,with=F],1,min)),
        "median_60"=(apply(my_training[,31:90,with=F],1,median)),
        "mean_60"=(apply(my_training[,31:90,with=F],1,mean)),
        "sd_60"=(apply(my_training[,31:90,with=F],1,sd)),
        
        "max_70"=(apply(my_training[,21:90,with=F],1,max)),
        "min_70"=(apply(my_training[,21:90,with=F],1,min)),
        "median_70"=(apply(my_training[,21:90,with=F],1,median)),
        "mean_70"=(apply(my_training[,21:90,with=F],1,mean)),
        "sd_70"=(apply(my_training[,21:90,with=F],1,sd)),
        
        "max_80"=(apply(my_training[,11:90,with=F],1,max)),
        "min_80"=(apply(my_training[,11:90,with=F],1,min)),
        "median_80"=(apply(my_training[,11:90,with=F],1,median)),
        "mean_80"=(apply(my_training[,11:90,with=F],1,mean)),
        "sd_80"=(apply(my_training[,11:90,with=F],1,sd)),
        
        "max_90"=(apply(my_training[,1:90,with=F],1,max)),
        "min_90"=(apply(my_training[,1:90,with=F],1,min)),
        "median_90"=(apply(my_training[,1:90,with=F],1,median)),
        "mean_90"=(apply(my_training[,1:90,with=F],1,mean)),
        "sd_90"=(apply(my_training[,1:90,with=F],1,sd))
    )]
    
    return(my_training)
}
selected_feature <- function(my_dataset, my_case, my_feature){
    if(max_90 %in% my_feature){my_training <- my_training[,':='("max_90"=(apply(my_training[,1:90,with=F],1,max)))]}
    if(min_90 %in% my_feature){my_training <- my_training[,':='("min_90"=(apply(my_training[,1:90,with=F],1,min)))]}
    if(median_90 %in% my_feature){my_training <- my_training[,':='("median_90"=(apply(my_training[,1:90,with=F],1,median)))]}
    if(mean_90 %in% my_feature){my_training <- my_training[,':='("mean_90"=(apply(my_training[,1:90,with=F],1,mean)))]}
    if(sd_90 %in% my_feature){my_training <- my_training[,':='("sd_90"=(apply(my_training[,1:90,with=F],1,sd)))]}
    if(fc_90 %in% my_feature){my_training <- my_training[,':='("fc_90"=(apply(my_training[,1:90,with=F],1,function(x){feature_coefficient(x,my_case)})))]}
    if(fm_90 %in% my_feature){my_training <- my_training[,':='("fm_90"=(apply(my_training[,1:90,with=F],1,feature_max_limit)))]}
    
    if(max_80 %in% my_feature){my_training <- my_training[,':='("max_80"=(apply(my_training[,11:90,with=F],1,max)))]}
    if(min_80 %in% my_feature){my_training <- my_training[,':='("min_80"=(apply(my_training[,11:90,with=F],1,min)))]}
    if(median_80 %in% my_feature){my_training <- my_training[,':='("median_80"=(apply(my_training[,11:90,with=F],1,median)))]}
    if(mean_80 %in% my_feature){my_training <- my_training[,':='("mean_80"=(apply(my_training[,11:90,with=F],1,mean)))]}
    if(sd_80 %in% my_feature){my_training <- my_training[,':='("sd_80"=(apply(my_training[,11:90,with=F],1,sd)))]}
    if(fc_80 %in% my_feature){my_training <- my_training[,':='("fc_80"=(apply(my_training[,11:90,with=F],1,function(x){feature_coefficient(x,my_case)})))]}
    if(fm_80 %in% my_feature){my_training <- my_training[,':='("fm_80"=(apply(my_training[,11:90,with=F],1,feature_max_limit)))]}
     
    if(max_70 %in% my_feature){my_training <- my_training[,':='("max_70"=(apply(my_training[,21:90,with=F],1,max)))]}
    if(min_70 %in% my_feature){my_training <- my_training[,':='("min_70"=(apply(my_training[,21:90,with=F],1,min)))]}
    if(median_70 %in% my_feature){my_training <- my_training[,':='("median_70"=(apply(my_training[,21:90,with=F],1,median)))]}
    if(mean_70 %in% my_feature){my_training <- my_training[,':='("mean_70"=(apply(my_training[,21:90,with=F],1,mean)))]}
    if(sd_70 %in% my_feature){my_training <- my_training[,':='("sd_70"=(apply(my_training[,21:90,with=F],1,sd)))]}
    if(fc_70 %in% my_feature){my_training <- my_training[,':='("fc_70"=(apply(my_training[,21:90,with=F],1,function(x){feature_coefficient(x,my_case)})))]}
    if(fm_70 %in% my_feature){my_training <- my_training[,':='("fm_70"=(apply(my_training[,21:90,with=F],1,feature_max_limit)))]}
    
    if(max_60 %in% my_feature){my_training <- my_training[,':='("max_60"=(apply(my_training[,31:90,with=F],1,max)))]}
    if(min_60 %in% my_feature){my_training <- my_training[,':='("min_60"=(apply(my_training[,31:90,with=F],1,min)))]}
    if(median_60 %in% my_feature){my_training <- my_training[,':='("median_60"=(apply(my_training[,31:90,with=F],1,median)))]}
    if(mean_60 %in% my_feature){my_training <- my_training[,':='("mean_60"=(apply(my_training[,31:90,with=F],1,mean)))]}
    if(sd_60 %in% my_feature){my_training <- my_training[,':='("sd_60"=(apply(my_training[,31:90,with=F],1,sd)))]}
    if(fc_60 %in% my_feature){my_training <- my_training[,':='("fc_60"=(apply(my_training[,31:90,with=F],1,function(x){feature_coefficient(x,my_case)})))]}
    if(fm_60 %in% my_feature){my_training <- my_training[,':='("fm_60"=(apply(my_training[,31:90,with=F],1,feature_max_limit)))]}
    
    if(max_50 %in% my_feature){my_training <- my_training[,':='("max_50"=(apply(my_training[,41:90,with=F],1,max)))]}
    if(min_50 %in% my_feature){my_training <- my_training[,':='("min_50"=(apply(my_training[,41:90,with=F],1,min)))]}
    if(median_50 %in% my_feature){my_training <- my_training[,':='("median_50"=(apply(my_training[,41:90,with=F],1,median)))]}
    if(mean_50 %in% my_feature){my_training <- my_training[,':='("mean_50"=(apply(my_training[,41:90,with=F],1,mean)))]}
    if(sd_50 %in% my_feature){my_training <- my_training[,':='("sd_50"=(apply(my_training[,41:90,with=F],1,sd)))]}
    if(fc_50 %in% my_feature){my_training <- my_training[,':='("fc_50"=(apply(my_training[,41:90,with=F],1,function(x){feature_coefficient(x,my_case)})))]}
    if(fm_50 %in% my_feature){my_training <- my_training[,':='("fm_50"=(apply(my_training[,41:90,with=F],1,feature_max_limit)))]}
    
    if(max_40 %in% my_feature){my_training <- my_training[,':='("max_40"=(apply(my_training[,41:90,with=F],1,max)))]}
    if(min_40 %in% my_feature){my_training <- my_training[,':='("min_40"=(apply(my_training[,41:90,with=F],1,min)))]}
    if(median_40 %in% my_feature){my_training <- my_training[,':='("median_40"=(apply(my_training[,41:90,with=F],1,median)))]}
    if(mean_40 %in% my_feature){my_training <- my_training[,':='("mean_40"=(apply(my_training[,41:90,with=F],1,mean)))]}
    if(sd_40 %in% my_feature){my_training <- my_training[,':='("sd_40"=(apply(my_training[,41:90,with=F],1,sd)))]}
    if(fc_40 %in% my_feature){my_training <- my_training[,':='("fc_40"=(apply(my_training[,41:90,with=F],1,function(x){feature_coefficient(x,my_case)})))]}
    if(fm_40 %in% my_feature){my_training <- my_training[,':='("fm_40"=(apply(my_training[,41:90,with=F],1,feature_max_limit)))]}
    
    if(max_30 %in% my_feature){my_training <- my_training[,':='("max_30"=(apply(my_training[,41:90,with=F],1,max)))]}
    if(min_30 %in% my_feature){my_training <- my_training[,':='("min_30"=(apply(my_training[,41:90,with=F],1,min)))]}
    if(median_30 %in% my_feature){my_training <- my_training[,':='("median_30"=(apply(my_training[,41:90,with=F],1,median)))]}
    if(mean_30 %in% my_feature){my_training <- my_training[,':='("mean_30"=(apply(my_training[,41:90,with=F],1,mean)))]}
    if(sd_30 %in% my_feature){my_training <- my_training[,':='("sd_30"=(apply(my_training[,41:90,with=F],1,sd)))]}
    if(fc_30 %in% my_feature){my_training <- my_training[,':='("fc_30"=(apply(my_training[,41:90,with=F],1,function(x){feature_coefficient(x,my_case)})))]}
    if(fm_30 %in% my_feature){my_training <- my_training[,':='("fm_30"=(apply(my_training[,41:90,with=F],1,feature_max_limit)))]}
    
    if(max_20 %in% my_feature){my_training <- my_training[,':='("max_20"=(apply(my_training[,41:90,with=F],1,max)))]}
    if(min_20 %in% my_feature){my_training <- my_training[,':='("min_20"=(apply(my_training[,41:90,with=F],1,min)))]}
    if(median_20 %in% my_feature){my_training <- my_training[,':='("median_20"=(apply(my_training[,41:90,with=F],1,median)))]}
    if(mean_20 %in% my_feature){my_training <- my_training[,':='("mean_20"=(apply(my_training[,41:90,with=F],1,mean)))]}
    if(sd_20 %in% my_feature){my_training <- my_training[,':='("sd_20"=(apply(my_training[,41:90,with=F],1,sd)))]}
    if(fc_20 %in% my_feature){my_training <- my_training[,':='("fc_20"=(apply(my_training[,41:90,with=F],1,function(x){feature_coefficient(x,my_case)})))]}
    if(fm_20 %in% my_feature){my_training <- my_training[,':='("fm_20"=(apply(my_training[,41:90,with=F],1,feature_max_limit)))]}
    
    if(max_10 %in% my_feature){my_training <- my_training[,':='("max_10"=(apply(my_training[,41:90,with=F],1,max)))]}
    if(min_10 %in% my_feature){my_training <- my_training[,':='("min_10"=(apply(my_training[,41:90,with=F],1,min)))]}
    if(median_10 %in% my_feature){my_training <- my_training[,':='("median_10"=(apply(my_training[,41:90,with=F],1,median)))]}
    if(mean_10 %in% my_feature){my_training <- my_training[,':='("mean_10"=(apply(my_training[,41:90,with=F],1,mean)))]}
    if(sd_10 %in% my_feature){my_training <- my_training[,':='("sd_10"=(apply(my_training[,41:90,with=F],1,sd)))]}
    if(fc_10 %in% my_feature){my_training <- my_training[,':='("fc_10"=(apply(my_training[,41:90,with=F],1,function(x){feature_coefficient(x,my_case)})))]}
    if(fm_10 %in% my_feature){my_training <- my_training[,':='("fm_10"=(apply(my_training[,41:90,with=F],1,feature_max_limit)))]}
    
    return(my_training)
}
    

save.image(file="regression.RData")
load(file="regression.RData")

## Apprentissage statistique
ml_reg <- function(my_training, my_test, my_case, my_algo, my_eGrid, my_metric, my_type="raw"){
    my_training <- training_set
    my_test <- test_set
    my_case <- 30
    my_algo <- "rpart2"
    my_eGrid <- expand.grid(.maxdepth = c(10, 20, 30))
    my_metric <- "qq_metric"
    my_type <- "raw"

    grid_1 <- expand.grid(.alpha = 1, .lambda = 0)
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

    vect <- sort(unique(my_training[,ind]))
    all_result <- NULL

    begin <- Sys.time()
    
    for(i in 1:length(vect)){
        # i <- 1
        time_fe <- Sys.time()
        
        temp_test_set <- my_test[ind==vect[i]][,-c("ind"),with=F]
        temp_train <- my_training[ind==vect[i]][,-c("ind"),with=F]
        
        if(nrow(temp_train)>=95){
            temp_fold <- cbind(data.frame("fold"=1:nrow(temp_train)), temp_train)
            folds_in <- createFolds(temp_fold[,"fold"], k = 5)
            folds_out <- folds_in
            row <- nrow(temp_train)
            folds_in$Fold1 <- 1:trunc(row*90/95)
            folds_in$Fold2 <- 1:trunc(row*91/95)
            folds_in$Fold3 <- 1:trunc(row*92/95)
            folds_in$Fold4 <- 1:trunc(row*93/95)
            folds_in$Fold5 <- 1:trunc(row*94/95)
            
            folds_out$Fold1 <- (trunc(row*90/95)+1):trunc(row*91/95)
            folds_out$Fold2 <- (trunc(row*91/95)+1):trunc(row*91/95)
            folds_out$Fold3 <- (trunc(row*92/95)+1):trunc(row*91/95)
            folds_out$Fold4 <- (trunc(row*93/95)+1):trunc(row*91/95)
            folds_out$Fold5 <- (trunc(row*94/95)+1):trunc(row*91/95)
            
        }else{
            print(c("1", i))
        }
      
        ## Génère les graines (reproduction des résultats) #length is = (n_repeats*nresampling)+1
        set.seed(1)
        seeds <- vector(mode = "list", length = 200)
        for(j in 1:199) seeds[[j]]<- sample.int(n=1000, 199)
        seeds[[200]]<-sample.int(1000, 200)

        temp_resp <- temp_train[, "response", with=F]
        fe_sel <- add_feature(temp_train[, -c("response"), with=F], my_case)
 
        tryCatch({
            #### Supprime les colonnes qui sont trop corrélées
            correlationMatrix <- cor(fe_sel)
            highlyCorrelated <- sort(findCorrelation(correlationMatrix, cutoff=0.95))
            temp_my_training <- as.data.table(fe_sel)
            while(length(highlyCorrelated)>0){
                temp_my_training <- temp_my_training[,-c(highlyCorrelated[1]),with=F]
                correlationMatrix <- cor(temp_my_training)
                highlyCorrelated <- sort(findCorrelation(correlationMatrix, cutoff=0.95))
                # print(highlyCorrelated)
                if(length(highlyCorrelated)>=2){
                    highlyCorrelated <- highlyCorrelated[2:length(highlyCorrelated)]
                }else{
                    break
                }
            }
            training <- cbind(fe_sel[,c(highlyCorrelated),with=F], temp_resp)
            if(ncol(training)<=2){
                training <- cbind(fe_sel, temp_resp)
            }
            
        }, error = function(e){
            training <- cbind(fe_sel, temp_resp)
        })
        

        control <- trainControl(seeds=seeds)
        cl <- makeCluster(3)
        registerDoSNOW(cl)
        model <- caret::train(response ~ .,
                              data=training,
                              method = "glmnet",
                              metric= "RMSE",
                              preProcess = c("center","scale"),
                              tuneGrid= grid_1,
                              trControl = control)
        stopCluster(cl)
        importance <- varImp(model, scale=FALSE)
        vec_t_value <- importance$importance$Overall
        p_value <-  abs(qt(c(.090), df=(nrow(temp_train)-1)))
        significant_variable <- data.table("value"=importance$importance$Overall, "names"=names(training[,-c("response"),with=F]))
        significant_variable <- significant_variable[value>=p_value][,names]
        # print(significant_variable)
        if(length(significant_variable)<=1){
            significant_variable <- names(fe_sel)
        }
        
        ## Paramètre la cross validation et la nature des résultats
        if(my_metric=="qq_metric"){
            if(my_type=="prob"){
                ## Paramètre la cross validation et la nature des résultats
                control <- trainControl(method="cv", number=5, seeds=seeds, index=folds_in, indexOut = folds_out, summaryFunction = quaqua_Summary_reg)
            }else if(my_type=="raw"){
                control <- trainControl(method="cv", number=5, seeds=seeds, index=folds_in, indexOut = folds_out, summaryFunction = quaqua_Summary_reg)
            }
        }else{
            if(my_type=="prob"){
                ## Paramètre la cross validation et la nature des résultats
                control <- trainControl(method="cv", number=5, seeds=seeds, index=folds_in, indexOut = folds_out)
            }else if(my_type=="raw"){
                control <- trainControl(method="cv", number=5, seeds=seeds, index=folds_in, indexOut = folds_out)
            }
        }
        # cl <- makeCluster(3)
        # registerDoSNOW(cl)
        # tryCatch({
        #     model <- caret::train(response ~ .,
        #                           data=  training,
        #                           method = "glm",
        #                           metric= my_metric,
        #                           preProcess = c("center","scale"),
        #                           trControl = control)
        # 
        #     
        #     significant_variable <- summary(model)$coeff[-1,4] < 0.05
        #     significant_variable <- unlist(lapply(significant_variable, function(x){ifelse(isTRUE(x), 1, 0)}))
        #     significant_variable <-  names(significant_variable[which(significant_variable==1)])
        #     if(length(significant_variable)==0){
        #         significant_variable <- names(fe_sel)
        #     }
        # }, error = function(e) {
        #     print(c(i, "1", e))
        #     significant_variable <- names(fe_sel)
        # })
        # stopCluster(cl)
        
        # time_fe <- as.numeric(-(time_fe - Sys.time()))
        print(c(i, time_fe - Sys.time()))
        train_rep <- temp_train[, "response", with=F]
        test_rep <- temp_test_set[, "response", with=F]
        
        temp_train <- add_feature(temp_train[, -c("response"), with=F], my_case)
        temp_test_set <- add_feature(temp_test_set[, -c("response"), with=F], my_case)
        
        temp_train  <- cbind(temp_train[,c(significant_variable),with=F], train_rep)
        temp_test_set  <- cbind(temp_test_set[,c(significant_variable),with=F], test_rep)
      
        cl <- makeCluster(3)
        registerDoSNOW(cl)
        
        tryCatch({
        model <- caret::train(response ~ .,
                data=temp_train,
                method = my_algo,
                tuneGrid = my_eGrid,
                metric= my_metric,
                preProcess = c("center","scale"),
                trControl = control)
        temp_predict <- data.table("ind"=vect[i], "predict"=predict(model, temp_test_set), 
                                   temp_test_set[,"response",with=F])
        
        }, error = function(e) {
            print(c(i, "2", e))
            temp_predict <- data.table("ind"=vect[i], "predict"=(-1),
                                       temp_test_set[,"response",with=F])
        })
        
        stopCluster(cl)
        
        
        all_result <- rbindlist(list(all_result, temp_predict))
    }
    
    if(my_metric=="qq_metric"){
        all <-  copy(all_result)[,':='("predict"=apply(all_result[,"predict",with=F],1,function(x){ifelse(x>=0, 1, 0)}),
                                           "response"=apply(all_result[,"response",with=F],1,function(x){ifelse(x>=0, 1, 0)})
        )]
        temp_pred <- unlist(all[,predict])
        temp_obs <- unlist(all[,response])
        score <- best_score(temp_pred, temp_obs)
    }else if(my_metric=="RMSE"){
        temp_pred <- unlist(all_result[,predict])
        temp_obs <- unlist(all_result[,response])
        score <- sqrt(sum((temp_pred - temp_obs)^2 , na.rm = TRUE)/length(temp_obs))
    }
    all <-  copy(all_result)[,':='("predict"=apply(all_result[,"predict",with=F],1,function(x){ifelse(x>=0, 1, 0)}),
                                   "response"=apply(all_result[,"response",with=F],1,function(x){ifelse(x>=0, 1, 0)})
    )]
    temp_pred <- unlist(all[,predict])
    temp_obs <- unlist(all[,response])
    obs_rest <- unlist(my_test[! ind %in% vect][,response])
    obs_rest <- unlist(lapply(obs_rest,function(x){ifelse(x>=0,1,0)}))
    temp_obs <- c(temp_obs, obs_rest)
    temp_pred <- c(temp_pred, rep(0, (length(temp_obs)-length(temp_pred))))
    score_final_entire <- best_score(temp_pred, temp_obs)
    time <- as.numeric(-(begin - Sys.time()))
    return(list("score"=score, "temps"=time, "case"=my_case, "algo"=my_algo, my_eGrid, "metric"=my_metric, "type"=my_type, "final"=score_final_entire))
}


## Modèle à tester
############## TREE #################
### CART
algo <- "rpart2"
eGrid <- expand.grid(.maxdepth = c(10, 20, 30))
simple_tree_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw")
simple_tree_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw")
simple_tree_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw")
simple_tree_10_rm <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw")
simple_tree_30_rm <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw")
simple_tree_60_rm <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw")


### Bagged CART
algo <- "treebag"
eGrid <- NULL
Bagged_tree_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw")
Bagged_tree_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw")
Bagged_tree_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw")
Bagged_tree_10_rm <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw")
Bagged_tree_30_rm <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw")
Bagged_tree_60_rm <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw")


### Boosted Tree
algo <- "blackboost"
eGrid <- expand.grid(.mstop = c(50, 100, 150, 200), .maxdepth = c(10, 20, 30))
Boosted_tree_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw")
Boosted_tree_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw")
Boosted_tree_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw")
Boosted_tree_10_rm <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw")
Boosted_tree_30_rm <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw")
Boosted_tree_60_rm <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw")


### Random Forest
algo <- "ranger"
eGrid <- expand.grid(.mtry = c(1,2,3,4,5,6,7))
rf_tree_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw")
rf_tree_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw")
rf_tree_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw")
rf_tree_10_rm <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw")
rf_tree_30_rm <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw")
rf_tree_60_rm <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw")



############## Linear Model #################
## Generalized Linear Model
algo <- "glmnet"
eGrid <- expand.grid(.alpha = 1, .lambda = 0)
# glm_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw")
# glm_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw")
# glm_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw")
# glm_10_rm <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw")
# glm_30_rm <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw")
# glm_60_rm <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw")


## L1 Logistic regression (lasso)
# algo <- "glmnet"
# eGrid <- expand.grid(.alpha = 1, .lambda = 10^-(-5:5))
# Lone_glm_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric")
# Lone_glm_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric")
# Lone_glm_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric")
# Lone_glm_10_rm <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE")
# Lone_glm_30_rm <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE")
# Lone_glm_60_rm <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw")


## L2 Logistic regression (ridge)
# algo <- "glmnet"
# eGrid <- expand.grid(.alpha = 0, .lambda = 10^-(-5:5))
# Ltwo_glm_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw")
# Ltwo_glm_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw")
# Ltwo_glm_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw")
# Ltwo_glm_10_rm <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw")
# Ltwo_glm_30_rm <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw")
# Ltwo_glm_60_rm <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw")


############### Plus proches voisins ###########
# algo <- "knn"
# eGrid <- expand.grid(.k = c(3, 6, 9, 12))
# knn_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw")
# knn_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw")
# knn_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw")
# knn_10_rm <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw")
# knn_30_rm <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw")
# knn_60_rm <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw")

############### Réseaux de neurones ###########
# algo <- "mlp"
# eGrid <- expand.grid(.size = 1:5)
# rn_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw")
# rn_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw")
# rn_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw")
# rn_10_rm <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw")
# rn_30_rm <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw")
# rn_60_rm <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw")


############### SVM ###########
### Linear Kernel
algo <- "svmLinear"
eGrid <- expand.grid(.C = c(10^(-5), 10^(-3), 10, 10^(3), 10^(6), 10^(9)))
svmLinear_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw")
svmLinear_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw")
svmLinear_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw")
svmLinear_10_rm <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw")
svmLinear_30_rm <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw")
svmLinear_60_rm <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw")


### Radial Kernel
algo <- "svmRadial"
eGrid <- expand.grid(.C = c(10^(-5), 10^(-3), 10, 10^(3), 10^(6), 10^(9)), sigma=c(10^(-5), 10^(-3), 10^(-1), 10))
svmRadial_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw")
svmRadial_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw")
svmRadial_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw")
svmRadial_10_rm <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw")
svmRadial_30_rm <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw")
svmRadial_60_rm <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw")


# # ### Radial Kernel
# algo <- "rvmLinear"
# eGrid <- NULL
# rvmLinear_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw")
# rvmLinear_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw")
# rvmLinear_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw")
# rvmLinear_10_rm <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw")
# rvmLinear_30_rm <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw")
# rvmLinear_60_rm <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw")
# # 
# # 
# # 
# # ### Radial Kernel
# algo <- "rvmRadial"
# eGrid <- expand.grid(.sigma=c(10^(-5), 10^(-3), 10^(-1), 10))
# rvmRadial_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw")
# rvmRadial_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw")
# rvmRadial_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw")
# rvmRadial_10_rm <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw")
# rvmRadial_30_rm <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw")
# rvmRadial_60_rm <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw")


### Bayesian Generalized Linear Model
algo <- "bayesglm"
eGrid <- NULL
# arm_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw")
# arm_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw")
# arm_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw")
# arm_10_rm <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw")
# arm_30_rm <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw")
# arm_60_rm <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw")


save.image(file="regression_3.RData")
