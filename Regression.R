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

### Filter training set
filtre_coefficient <- function(my_vector, my_case=60){
    my_vector <- temp_train[1,2:92,with=F]
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
    for(i in 1:length(vect)){
        temp_train <- my_training[ind==vect[i]]
        temp_train <- temp_train[,":="("fc"=apply(temp_train[,2:92,with=F],1,function(x){filtre_coefficient(x, 60)}),
                                           "fm"=apply(temp_train[,2:92,with=F],1,filtre_max_limit))]
        if(all(temp_train[,fc]==0) && all(temp_train[,fm]==0)){
            filtered <- c(filtered, vect[i])
        }
        print(c(i, length(filtered)))
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
    obs <- unlist(lapply(obs, function(x){ifelse(x>=0, 2, 1)}))
    
    if("pred"  %in% names(data)){
        pred <- unlist(data$pred)
        pred <- unlist(lapply(pred, function(x){ifelse(is.na(x), 1, x)}))
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
    # 
    # if (any(is.nan(out))) 
    #     out[is.nan(out)] <- NA
    # out
    
    # }, error = function(e) {
    #     saveRDS(data, "data.rds")
    #     print(e)
    # 
    # })
}

save.image(file="regression.RData")
load(file="regression.RData")


## depth_selection
depth_selection <- function(my_training, my_test, cases=c(10,30,60)){
    # my_training <- training_set
    # my_test <- test_set
    # cases <- c(10,30,60)
    grid_1 <- expand.grid(.alpha = 1, .lambda = 0)
    for(k in 1:length(cases)){
        # k <- 1
        my_case <- cases[k]
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
        
        vect_table <- NULL
        for(i in 1:length(vect)){
            # i<- 1
            print(i)
            temp_test_set <- my_test[ind==vect[i]][,-c("ind"),with=F]
            temp_train <- my_training[ind==vect[i]][,-c("ind"),with=F]
            temp_fold <- cbind(data.frame("fold"=1:nrow(temp_train)), temp_train)
            folds <- createFolds(temp_fold[,"fold"], k = 6)
            
            row <- nrow(temp_train)
            folds$Fold1 <- 1:trunc((row/6))
            folds$Fold2 <- 1:(trunc((2*row/6)))
            folds$Fold3 <- 1:(trunc((3*row/6)))
            folds$Fold4 <- 1:(trunc((4*row/6)))
            folds$Fold5 <- 1:(trunc((5*row/6)))
            folds$Fold6 <- 1:row
            
            ## Génère les graines (reproduction des résultats) #length is = (n_repeats*nresampling)+1
            set.seed(1)
            seeds <- vector(mode = "list", length = 200)
            for(i in 1:199) seeds[[i]]<- sample.int(n=1000, 199)
            seeds[[200]]<-sample.int(1000, 200)
            
            ## Détermine rapidement la meilleure longueur
            depth_score <- NULL
            
            control <- trainControl(method="cv", number=5, seeds=seeds, index=folds[1:5], indexOut = folds[2:6])
            for(j in 1:9){
                depth_training <- temp_train[,(1+(j-1)*10):91,with=F]
                
                cl <- makeCluster(3)
                registerDoSNOW(cl)
                model <- caret::train(response ~ .,
                                      data=depth_training,
                                      method = "glmnet",
                                      metric="RMSE",
                                      tuneGrid =grid_1,
                                      trControl = control)
                
                depth_score <- c(depth_score, model$results$RMSE)
            }
            names(depth_score) <- c(1:9)
            my_max <- max(depth_score)
            depth_score <- depth_score[which(depth_score==my_max)]
            depth_score <- as.numeric(names(depth_score[length(depth_score)]))
            depth_score <- 1+(depth_score-1)*10
            vect_table <- rbindlist(list(vect_table, data.table("case"=my_case, "ind"==vect[i], "depth"=depth_score)))
        }
    }
    return(vect_table)
}

table_depth <- depth_selection(training_set, test_set, c(10,30,60))



## Apprentissage statistique
ml_reg <- function(my_training, my_test, my_case, my_algo, my_eGrid, my_metric, my_type="raw", my_table_depth){
    # my_training <- training_set
    # my_test <- test_set
    # my_case <- 10
    # my_algo <- "rpart2"
    # my_eGrid <- expand.grid(.maxdepth = c(10, 20, 30))
    # my_metric <- "qq_metric"
    # my_type <- "raw"
    # grid_1 <- expand.grid(.alpha = 0, .lambda = 0)
    # my_table_depth <- table_depth
    
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
        # i<- 1
        print(i)
        temp_test_set <- my_test[ind==vect[i]][,-c("ind"),with=F]
        temp_train <- my_training[ind==vect[i]][,-c("ind"),with=F]
        temp_fold <- cbind(data.frame("fold"=1:nrow(temp_train)), temp_train)
        folds <- createFolds(temp_fold[,"fold"], k = 6)

        row <- nrow(temp_train)
        folds$Fold1 <- 1:trunc((row/6))
        folds$Fold2 <- 1:(trunc((2*row/6)))
        folds$Fold3 <- 1:(trunc((3*row/6)))
        folds$Fold4 <- 1:(trunc((4*row/6)))
        folds$Fold5 <- 1:(trunc((5*row/6)))
        folds$Fold6 <- 1:row
      
        ## Génère les graines (reproduction des résultats) #length is = (n_repeats*nresampling)+1
        set.seed(1)
        seeds <- vector(mode = "list", length = 200)
        for(i in 1:199) seeds[[i]]<- sample.int(n=1000, 199)
        seeds[[200]]<-sample.int(1000, 200)


        depth_score <- my_table_depth[ind==vect[i]][case==my_case][,depth]
        temp_train <- temp_train[,depth_score:91,with=F]
        temp_test_set <- temp_test_set[,depth_score:91,with=F]
        
        
        ## Paramètre la cross validation et la nature des résultats
        if(my_metric=="qq_metric"){
            if(my_type=="prob"){
                ## Paramètre la cross validation et la nature des résultats
                control <- trainControl(method="cv", number=5, seeds=seeds, index=folds[1:5], indexOut = folds[2:6], summaryFunction = quaqua_Summary_reg)
            }else if(my_type=="raw"){
                control <- trainControl(method="cv", number=5, seeds=seeds, index=folds[1:5], indexOut = folds[2:6], summaryFunction = quaqua_Summary_reg)
            }
        }else{
            if(my_type=="prob"){
                ## Paramètre la cross validation et la nature des résultats
                control <- trainControl(method="cv", number=5, seeds=seeds, index=folds[1:5], indexOut = folds[2:6])
            }else if(my_type=="raw"){
                control <- trainControl(method="cv", number=5, seeds=seeds, index=folds[1:5], indexOut = folds[2:6])
            }
        }
        
        cl <- makeCluster(3)
        registerDoSNOW(cl)
        
        model <- caret::train(response ~ .,
                data=temp_train,
                method = my_algo,
                tuneGrid = my_eGrid,
                metric= my_metric,
                trControl = control)

        stopCluster(cl)
        
        temp_predict <- data.table("ind"=vect[i], "predict"=predict(model, temp_test_set), 
                                   temp_test_set[,"response",with=F])
        all_result <- rbindlist(list(all_result, temp_predict))
    }
    
    if(my_metric=="qq_metric"){
        all_result <-  all_result[,':='("predict"=apply(all_result[,"predict",with=F],1,function(x){ifelse(x>=0, 1, 0)}),
                                           "response"=apply(all_result[,"response",with=F],1,function(x){ifelse(x>=0, 1, 0)})
        )]
        temp_pred <- unlist( all_result[,pred])
        temp_obs <- unlist(all_result[,obs])
        score <- best_score(temp_pred, temp_obs)
    }else if(my_metric=="RMSE"){
        temp_pred <- unlist(all_result[,predict])
        temp_obs <- unlist(all_result[,response])
        score <- sqrt(sum((temp_pred - temp_obs)^2 , na.rm = TRUE)/length(temp_obs))
    }
    temp_test <- unlist(test_set[! ind %in% vect][,response])
    temp_test <- unlist(lapply(test_set,function(x){ifelse(x>=0,1,0)}))
    temp_obs <- c(temp_obs, temp_test)
    temp_pred <- c(temp_obs, rep(0,length(temp_test)))
    score_final_entire <- best_score(temp_pred, temp_obs)
    time <- as.numeric(-(begin - Sys.time()))
    return(list("score"=score, "temps"=time, "case"=my_case, "algo"=my_algo, my_eGrid, "metric"=my_metric, "type"=my_type, "final"=score_final_entire))
}


## Modèle à tester
############## TREE #################
### CART
algo <- "rpart2"
eGrid <- expand.grid(.maxdepth = c(10, 20, 30))
simple_tree_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw", table_depth)
simple_tree_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw", table_depth)
simple_tree_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw", table_depth)
simple_tree_10_acc <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw", table_depth)
simple_tree_30_acc <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw", table_depth)
simple_tree_60_acc <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw", table_depth)


### Bagged CART
algo <- "treebag"
eGrid <- NULL
Bagged_tree_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", , table_depth)
Bagged_tree_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", , table_depth)
Bagged_tree_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", , table_depth)
Bagged_tree_10_acc <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", , table_depth)
Bagged_tree_30_acc <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", , table_depth)
Bagged_tree_60_acc <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", , table_depth)


### Boosted Tree
algo <- "blackboost"
eGrid <- expand.grid(.mstop = c(50, 100, 150, 200), .maxdepth = c(10, 20, 30))
Boosted_tree_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", , table_depth)
Boosted_tree_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", , table_depth)
Boosted_tree_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", , table_depth)
Boosted_tree_10_acc <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", , table_depth)
Boosted_tree_30_acc <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", , table_depth)
Boosted_tree_60_acc <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", , table_depth)


### Random Forest
algo <- "ranger"
eGrid <- expand.grid(.mtry = floor((1:floor(length(feature_10)/10))*10))
rf_tree_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", , table_depth)
rf_tree_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", , table_depth)
rf_tree_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", , table_depth)
rf_tree_10_acc <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", , table_depth)
rf_tree_30_acc <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", , table_depth)
rf_tree_60_acc <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", , table_depth)



############## Linear Model #################
## Generalized Linear Model
algo <- "glm"
eGrid <- NULL
glm_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", , table_depth)
glm_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", , table_depth)
glm_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", , table_depth)
glm_10_acc <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", , table_depth)
glm_30_acc <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", , table_depth)
glm_60_acc <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", , table_depth)


## L1 Logistic regression (lasso)
algo <- "glmnet"
eGrid <- expand.grid(.alpha = 1, .lambda = 10^-(-5:5))
Lone_glm_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", , table_depth)
Lone_glm_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", , table_depth)
Lone_glm_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", , table_depth)
Lone_glm_10_acc <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", , table_depth)
Lone_glm_30_acc <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", , table_depth)
Lone_glm_60_acc <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw", table_depth)


## L2 Logistic regression (ridge)
algo <- "glmnet"
eGrid <- expand.grid(.alpha = 0, .lambda = 10^-(-5:5))
Ltwo_glm_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw", table_depth)
Ltwo_glm_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw", table_depth)
Ltwo_glm_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw", table_depth)
Ltwo_glm_10_acc <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw", table_depth)
Ltwo_glm_30_acc <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw", table_depth)
Ltwo_glm_60_acc <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw", table_depth)


############### Plus proches voisins ###########
algo <- "knn"
eGrid <- expand.grid(.k = c(3, 6, 9, 12))
knn_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw", table_depth)
knn_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw", table_depth)
knn_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw", table_depth)
knn_10_acc <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw", table_depth)
knn_30_acc <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw", table_depth)
knn_60_acc <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw", table_depth)

############### Réseaux de neurones ###########
algo <- "mlp"
eGrid <- expand.grid(.size = 1:5)
rn_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw", table_depth)
rn_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw", table_depth)
rn_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw", table_depth)
rn_10_acc <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw", table_depth)
rn_30_acc <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw", table_depth)
rn_60_acc <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw", table_depth)


############### SVM ###########
### Linear Kernel
algo <- "svmLinear"
eGrid <- expand.grid(.C = c(10^(-5), 10^(-3), 10, 10^(3), 10^(6), 10^(9)))
svmLinear_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw", table_depth)
svmLinear_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw", table_depth)
svmLinear_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw", table_depth)
svmLinear_10_acc <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw", table_depth)
svmLinear_30_acc <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw", table_depth)
svmLinear_60_acc <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw", table_depth)


### Radial Kernel
algo <- "svmRadial"
eGrid <- expand.grid(.C = c(10^(-5), 10^(-3), 10, 10^(3), 10^(6), 10^(9)), sigma=c(10^(-5), 10^(-3), 10^(-1), 10))
svmRadial_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw", table_depth)
svmRadial_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw", table_depth)
svmRadial_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw", table_depth)
svmRadial_10_acc <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw", table_depth)
svmRadial_30_acc <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw", table_depth)
svmRadial_60_acc <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw", table_depth)


### Relevance Vector Machine wiLinear Kernel
algo <- "rvmLinear"
eGrid <- expand.grid(.C = c(10^(-5), 10^(-3), 10, 10^(3), 10^(6), 10^(9)))
svmLinear_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw", table_depth)
svmLinear_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw", table_depth)
svmLinear_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw", table_depth)
svmLinear_10_acc <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw", table_depth)
svmLinear_30_acc <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw", table_depth)
svmLinear_60_acc <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw", table_depth)


### Radial Kernel
algo <- "rvmRadial"
eGrid <- expand.grid(.C = c(10^(-5), 10^(-3), 10, 10^(3), 10^(6), 10^(9)), sigma=c(10^(-5), 10^(-3), 10^(-1), 10))
svmRadial_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw", table_depth)
svmRadial_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw", table_depth)
svmRadial_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw", table_depth)
svmRadial_10_acc <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw", table_depth)
svmRadial_30_acc <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw", table_depth)
svmRadial_60_acc <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw", table_depth)


### Bayesian Generalized Linear Model
algo <- "arm"
eGrid <- NULL
arm_10_qq <- ml_reg(training_set, test_set, 10, algo, eGrid, "qq_metric", "raw", table_depth)
arm_30_qq <- ml_reg(training_set, test_set, 30, algo, eGrid, "qq_metric", "raw", table_depth)
arm_60_qq <- ml_reg(training_set, test_set, 60, algo, eGrid, "qq_metric", "raw", table_depth)
arm_10_kap <- ml_reg(training_set, test_set, 10, algo, eGrid, "RMSE", "raw", table_depth)
arm_30_kap <- ml_reg(training_set, test_set, 30, algo, eGrid, "RMSE", "raw", table_depth)
arm_60_kap <- ml_reg(training_set, test_set, 60, algo, eGrid, "RMSE", "raw", table_depth)



