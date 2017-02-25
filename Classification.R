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
if(!isTRUE(require(ipred))) install.packages("ipred")


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
library(ipred)


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

#####################################################
#### ##################### FIltres ##################
#####################################################
### Filter  coefficient
filtre_coefficient <- function(my_vector, my_case=10, my_limit=1){
    my_vector <- unlist(my_vector)
    if(length(my_vector[!is.na(my_vector)>2]) && !is.na(my_vector[length(my_vector)])){
        l_last_value <- my_vector[length(my_vector)]
        my_vector <- diff(my_vector, 1)
        my_vector <- my_vector[!is.na(my_vector)]
        ifelse(any(unlist(lapply(my_vector, function(x){x*my_case+l_last_value}))>=my_limit), return(1), return(0))
    }else{
        return(0)
    }
}

## Filter max limit
filtre_max_limit <- function(my_vector, my_limit=1){
    ifelse(any(unlist(my_vector)>=my_limit*0.90), 1, 0)
}

## Filtre du training_set
training_set <- training_set[,":="("fc"=apply(training_set[,2:91,with=F],1,function(x){filtre_coefficient(x, 60, 1)}),
                                       "fm"=apply(training_set[,2:91,with=F],1,filtre_max_limit))][fc==1][fm==1][,-c("fc", "fm"),with=F]


#####################################################
#### #####################  Score functions #########
#####################################################
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
quaqua_Summary_cl <- function (data, lev = NULL, model = NULL) {
    library(data.table)
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
    
    
    # tryCatch({
      
    # saveRDS(data, "data.rds")
    # data <- readRDS("data.rds")
    
    obs <- as.character(unlist(data$obs))
    
    if("YES" %in% obs){
        obs[obs=="YES"] <- 1
    }
    
    if("NO" %in% obs){
        obs[obs=="NO"] <- 0
    }
    # obs <- unlist(lapply(obs, function(x){ifelse(x=="YES", 1,0)}))
    # obs[obs=="YES"] <- 1
    # obs[obs=="NO"] <- 0
    
    if("YES" %in% names(data)){
        pred <- unlist(data$YES)
        pred[is.na(pred)] <- 0
        # pred <- unlist(lapply(pred, function(x){ifelse(is.na(x), 0, x)}))
        if(length(pred)==0){
            score <- NA
        }else{
            value_max <- 1
            value_min <- 0
            init_seuil <- 10:0/10
            score <- 0
            th <- sort(unique(c(0, pred, 1)),decreasing = TRUE)
            temp_th <- th
            for(i in 1:9){
                # i <- 1
                puis <- 10^(i)
                temp_seuil <- unique(trunc(temp_th*(puis))/puis)
                temp_seuil <- temp_seuil[which(temp_seuil<=value_max)][which(temp_seuil>=value_min)]
                temp_seuil <- temp_seuil[!is.na(temp_seuil)]
                if(length(temp_seuil)==2){
                    values <- temp_seuil[1:2]
                    s_1 <- best_score(pred, obs, values[1])
                    s_2 <- best_score(pred, obs, values[2])
                    score <- max(s_1, s_2)
                    break
                }else{
                    for(j in 1:(length(temp_seuil)-2)){
                        # j <- 1
                        # print(j)
                        values <- temp_seuil[j:(j+2)]
                        s_1 <- best_score(pred, obs, values[1])
                        s_2 <- best_score(pred, obs, values[2])
                        s_3 <- best_score(pred, obs, values[3])
                        # print(c(s_1, s_2, s_3))
                        if(s_1==1 || s_2==1 || s_3==1){
                            score <- 1
                            break
                        }
                        if(s_3<s_2){
                            value_max <- values[1]
                            value_min <- values[3]
                            break
                        }
                    }
                }
                if(i==6){score <- max(s_1, s_2, s_3)}
                if(score==1){break}
            }
            out <- score  
        }
    }else if("pred"  %in% names(data)){
        pred <- as.character(unlist(data$pred))
        pred[is.na(pred)] <- 0
        # pred <- unlist(lapply(pred, function(x){ifelse(is.na(x), 0, x)}))
        if(length(pred)==0){
            score <- NA
        }else{
            pred[pred=="YES"] <- 1
            pred[pred=="NO"] <- 0
            # pred <- unlist(lapply(pred, function(x){ifelse(x=="YES", 1,0)}))
            score <- best_score(pred, obs, 0.5)
        }
    }else{
        score <- NA
    }
    out <- score
    names(out) <- "qq_metric"
    out
}

## Métrique définie (contient deux autres fonctions pour les intégrer dans les clusters de DoSnow)
ROC_Summary <- function (data, lev = NULL, model = NULL) {
    library(pROC)
        # saveRDS(data, "data.rds")
        # data <- readRDS("data.rds")

        obs <- as.character(unlist(data$obs))
        obs[obs=="YES"] <- 1
        obs[obs=="NO"] <- 0
 
        if("YES" %in% names(data)){
            pred <- unlist(data$YES)
            pred[is.na(pred)] <- 0
            # pred <- unlist(lapply(pred, function(x){ifelse(is.na(x), 0, x)}))
            score <- roc(response=obs,predictor= pred)
        }else if("pred"  %in% names(data)){
            pred <- as.character(unlist(data$pred))
            pred[is.na(pred)] <- 0
            # pred[pred=="NO"] <- 0
            # pred[pred=="YES"] <- 0
            pred <- unlist(lapply(pred, function(x){ifelse(is.na(x) || x=="NO", 0, 1)}))
            score <- roc(response=obs,predictor= pred)
        }
        out <- as.numeric(score$auc)
        names(out) <- "my_ROC"

        # if (any(is.nan(out)))
        #     out[is.nan(out)] <- NA
        out
}


#####################################################
#### #####################  FE #########
#####################################################

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

## Trouve les meilleures features pour chaque cas
get_features <- function(my_training=training_set, my_test=test_set, case=c(10, 30, 60)){
    # my_training <- training_set
    # my_test <- test_set
    table_feature <- NULL
    init_training <- my_training
    init_test <- my_test
    grid_1 <- expand.grid(.alpha = 1, .lambda = 0)
    for(c in 1:length(case)){
        # c <- 2
        my_case <- case[c]
        if(my_case==10){
            my_training <- init_training[,':='("response"=y10)][,-c("y10", "y30", "y60", "ind"),with=F]
            my_test <- init_test [,':='("response"=y10)][,-c("y10", "y30", "y60"),with=F]
        }else if(my_case==30){
            my_training <- init_training[,':='("response"=y30)][,-c("y10", "y30", "y60", "ind"),with=F]
            my_test <-  init_test[,':='("response"=y30)][,-c("y10", "y30", "y60"),with=F]
        }else if(my_case==60){
            my_training <- init_training[,':='("response"=y60)][,-c("y10", "y30", "y60", "ind"),with=F]
            my_test <- init_test[,':='("response"=y60)][,-c("y10", "y30", "y60"),with=F]
        }
          
        ### Génération de features mathématique simple pour chaque dimension d%%10 de 10 à 90
        temp_resp <- my_training[, "response", with=F]
        fe_sel <- add_feature(my_training[, -c("response"), with=F], my_case)
        training <- cbind(fe_sel, temp_resp)
      
        ## Génère les graines (reproduction des résultats) #length is = (n_repeats*nresampling)+1
        set.seed(1)
        seeds <- vector(mode = "list", length = 200)
        for(j in 1:199) seeds[[j]]<- sample.int(n=1000, 199)
        seeds[[200]]<-sample.int(1000, 200)
            
        my_metric <- "my_ROC"
        ## Paramètre la cross validation et la nature des résultats
        control <- trainControl(seeds=seeds)
       
        training <- training[,response:=as.factor(response)]
        temp_train <- as.data.frame(training[,-c("response"),with=F])
        temp_test <- unlist((training[,"response",with=F]))
            
        #### Supprime les colonnes qui sont trop corrélées
        correlationMatrix <- cor(temp_train)
        highlyCorrelated <- sort(findCorrelation(correlationMatrix, cutoff=0.95))
        temp_my_training <- as.data.table(temp_train)
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
        temp_train <- as.data.frame(temp_my_training)
        
        ## Test de student
        cl <- makeCluster(3)
        registerDoSNOW(cl)
        model <- caret::train(x=temp_train,
                                  y=temp_test,
                                  method = "glmnet",
                                  preProcess = c("center","scale"),
                                  tuneGrid= grid_1,
                                  trControl = control)
        stopCluster(cl)
        importance <- varImp(model, scale=FALSE)
        vec_t_value <- importance$importance$Overall
        p_value <- abs(qt(c(.095), df=(nrow(temp_train)-1)))
        significant_variable <- data.table("value"=importance$importance$Overall, "names"=names(temp_train))
        significant_variable <- significant_variable[value>=p_value][,names]
        if(length(significant_variable)<=1){
            significant_variable <- names(fe_sel)
        }
        table_feature <- rbindlist(list(table_feature, data.table("case"=my_case, "var"=significant_variable)))
    }
    return(table_feature)
}
table_feature <- get_features(my_training=training_set, my_test=test_set, case=c(10, 30, 60))

## Applique les features
tr_feature <- function(my_training=training_set, my_test=test_set, table_feature, my_case, my_metric){
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
  
    feature <- unlist(table_feature[case==my_case][,var])
    
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
    
  
    my_training <- cbind(data.table("type"="training"), my_training[,c("ind", feature, "response"),with=F][,response:=as.factor(response)])
    my_test <- cbind(data.table("type"="test"), my_test[,c("ind", feature, "response"),with=F][,response:=as.factor(response)])
    return(rbindlist(list(my_training, my_test)))
}
data_10 <- tr_feature(training_set, test_set, table_feature, 10)
data_30 <- tr_feature(training_set, test_set, table_feature, 30)
data_60 <- tr_feature(training_set, test_set, table_feature, 60)


save.image(file="classification.RData")
load(file="classification.RData")

big_result_table <- NULL
### Machine learning
ml_cl <- function(my_data=data_10, my_case, my_algo, my_eGrid, my_metric, my_type="prob"){
    # my_algo <- "rpart2"
    # my_eGrid <- expand.grid(.maxdepth = c(10, 20, 30))
    # my_data <- data_10
    # my_case <- 10
    # my_metric <- "my_ROC"
    # my_type <- "raw"

    # #
    my_training <- my_data[type=="training"][,-c("type"),with=F]
    my_test <- my_data[type=="test"][,-c("type"),with=F]
    temp_train <- as.data.frame(my_training[,-c("ind", "response"),with=F])
    temp_test <- unlist((my_training[,"response",with=F]))
    
    ## Créé les k-folds
    temp_fold <- cbind(data.frame("fold"=1:nrow(temp_train)), temp_train)
    nb_row <- 1:nrow(my_training)
    folds_in <- createFolds(temp_fold[,"fold"], k = 5)
    folds_out <- folds_in
    folds_in$Fold1 <- nb_row[! nb_row %in% folds_in$Fold1]
    folds_in$Fold2 <- nb_row[! nb_row %in% folds_in$Fold2]
    folds_in$Fold3 <- nb_row[! nb_row %in% folds_in$Fold3]
    folds_in$Fold4 <- nb_row[! nb_row %in% folds_in$Fold4]
    folds_in$Fold5 <- nb_row[! nb_row %in% folds_in$Fold5]
    
    ## Génère les séquences pour choisies pour chaque itération
    ## (reproduction des résultats) #length is = (n_repeats*nresampling)+1
    ## Ces graines sont beaucoup plus nombreuses que nécessaire pour éviter de les changer à chaque modification
    ## d'algorithme et de choix d'yhperparamètre
    set.seed(1)
    seeds <- vector(mode = "list", length = 200)
    for(i in 1:199) seeds[[i]]<- sample.int(n=1000, 199)
    seeds[[200]]<-sample.int(1000, 200)
    
    if(my_metric=="qq_metric"){
        if(my_type=="prob"){
            ## Paramètre la cross validation et la nature des résultats
            control <- trainControl(method="cv", number=5, seeds=seeds, index=folds_in, indexOut=folds_out, classProbs=TRUE, summaryFunction = quaqua_Summary_cl)
        }else if(my_type=="raw"){
            control <- trainControl(method="cv", number=5, seeds=seeds, index=folds_in, indexOut=folds_out, classProbs=FALSE, summaryFunction = quaqua_Summary_cl)
        }
    }else if(my_metric=="my_ROC"){
        if(my_type=="prob"){
            ## Paramètre la cross validation et la nature des résultats
            control <- trainControl(method="cv", number=5, seeds=seeds, index=folds_in, indexOut=folds_out, summaryFunction=ROC_Summary, classProbs=TRUE)
        }else if(my_type=="raw"){
            control <- trainControl(method="cv", number=5, seeds=seeds, index=folds_in, indexOut=folds_out, summaryFunction=ROC_Summary, classProbs=FALSE)
        }
    }else{
        if(my_type=="prob"){
            ## Paramètre la cross validation et la nature des résultats
            control <- trainControl(method="cv", number=5, seeds=seeds, index=folds_in, indexOut=folds_out, classProbs=TRUE)
        }else if(my_type=="raw"){
            control <- trainControl(method="cv", number=5, seeds=seeds, index=folds_in, indexOut=folds_out, classProbs=FALSE)
        }
    }

    ## Initialise le cluster
    cl <- makeCluster(3)
    registerDoSNOW(cl)
    
    ## Enregistre le temps
    begin <- proc.time()
    
    ## Entrainement du modèle
    model <-  caret::train(
            x = temp_train, 
            y =  temp_test,
            method = my_algo,
            tuneGrid = my_eGrid,
            trControl = control,
            preProcess = c("center","scale"),
            metric=my_metric)

    ## arrêt du cluster
    stopCluster(cl)
    
    ## Un fois qu'on a le meilleure modèle, on l'applique sur le test set
    
    if(my_metric=="qq_metric"){
        if(my_type=="prob"){
            temp_test <- my_test[,-c("ind", "response"),with=F]
            sc_test <- cbind(as.data.table(predict(model, newdata = temp_test, type=my_type)$YES), my_test[,"response",with=F])
            setnames(sc_test, c("YES", "obs"))
            score <- quaqua_Summary_cl(data=sc_test)
        }else if(my_type=="raw"){
            temp_test <- my_test[,-c("ind", "response"),with=F]
            sc_test <- cbind(as.data.table(predict(model, newdata = temp_test, type=my_type)), my_test[,"response",with=F])
            setnames(sc_test, c("pred", "obs"))
            score <- quaqua_Summary_cl(data=sc_test)
        }
        score <- round(score, digits=3)
        time <- round(begin - proc.time(), digits=0)
        big_result_table <<- rbindlist(list(big_result_table, data.table(my_case, my_algo, my_metric, score, time[3])))
        return(list("score"=score, "temps"=time[3], "cas"=my_case, "algo"=my_algo, "metric"=my_metric))
    }else if(my_metric=="Accuracy"){
        temp_test <- my_test[,-c("ind", "response"),with=F]
        sc_test <- cbind(as.data.table(my_test[,"response",with=F]), predict(model, newdata = temp_test, type="raw"))
        setnames(sc_test, c("my_obs", "my_pred"))
        sc_test <- sc_test[,':='("my_obs"=apply(sc_test[,"my_obs",with=F],1,function(x){ifelse(x=="YES", 1, 0)}))]
        sc_test <- sc_test[,':='("my_pred"=apply(sc_test[,"my_pred",with=F],1,function(x){ifelse(x=="NO"|| is.na(x), 0, 1)}))]
        sc_test <- sc_test[,':='("detail"=apply(sc_test[,1:2,with=F],1,function(x){get_confusion(x)}))]
        TP <- nrow(sc_test[detail=="TP"])
        FP <- nrow(sc_test[detail=="FP"])
        FN <- nrow(sc_test[detail=="FN"])
        TN <- nrow(sc_test[detail=="TN"])
        score <- (TP+TN)/nrow(sc_test)
        score <- round(score, digits=3)
        time <- round(begin - proc.time(), digits=0)
        big_result_table <<- rbindlist(list(big_result_table, data.table(my_case, my_algo, my_metric, score, time[3])))
        return(list("score"=score, "temps"=time[3], "cas"=my_case, "algo"=my_algo, "metric"=my_metric))
    }else if(my_metric=="Kappa"){
        temp_test <- my_test[,-c("ind", "response"),with=F]
        sc_test <- cbind(as.data.table(my_test[,"response",with=F]), as.data.table(predict(model, newdata = temp_test, type="raw")))
        setnames(sc_test, c("my_obs", "my_pred"))
        sc_test <- sc_test[,':='("my_obs"=apply(sc_test[,"my_obs",with=F],1,function(x){ifelse(x=="YES", 1, 0)}))]
        sc_test <- sc_test[,':='("my_pred"=apply(sc_test[,"my_pred",with=F],1,function(x){ifelse(x=="NO"|| is.na(x), 0, 1)}))]
        sc_test <- sc_test[,':='("detail"=apply(sc_test[,1:2,with=F],1,function(x){get_confusion(x)}))]
        TP <- nrow(sc_test[detail=="TP"])
        FP <- nrow(sc_test[detail=="FP"])
        FN <- nrow(sc_test[detail=="FN"])
        TN <- nrow(sc_test[detail=="TN"])
        po <- (TP+TN)/(TP+FN+FP+TN)
        ma <- ((TP+FN)*(TP+FP))/(TP+FN+FP+TN)
        mb <- ((TN+FN)*(TN+FP))/(TP+FN+FP+TN)
        pe <- (ma+mb)/(TP+FN+FP+TN)
        score <- round((po - pe)/(1 - pe), digits=3)
        time <- round(begin - proc.time(), digits=0)
        big_result_table <<- rbindlist(list(big_result_table, data.table(my_case, my_algo, my_metric, score, time[3])))
        return(list("score"=score, "temps"=time[3], "cas"=my_case, "algo"=my_algo, "metric"=my_metric))
    }else if(my_metric=="my_ROC"){
        if(my_type=="prob"){
            temp_test <- my_test[,-c("ind", "response"),with=F]
            sc_test <- cbind(as.data.table(my_test[,"response",with=F]), as.data.table(predict(model, newdata = temp_test, type="prob")$YES))
            setnames(sc_test, c("my_obs", "my_pred"))
            sc_test <- sc_test[,':='("my_obs"=apply(sc_test[,"my_obs",with=F],1,function(x){ifelse(x=="YES", 1, 0)}))]
            # sc_test <- sc_test[,':='("my_pred"=apply(sc_test[,"my_pred",with=F],1,function(x){ifelse(x=="NO"|| is.na(x), 0, 1)}))]
            score <- roc(response= unlist(sc_test[,my_obs]), predictor=unlist(sc_test[,my_pred]))
            score <- round(as.numeric(score$auc), digits=3)
        }else if(my_type=="raw"){
            temp_test <- my_test[,-c("ind", "response"),with=F]
            sc_test <- cbind(as.data.table(my_test[,"response",with=F]), as.data.table(predict(model, newdata = temp_test, type="raw")))
            setnames(sc_test, c("my_obs", "my_pred"))
            sc_test <- sc_test[,':='("my_obs"=apply(sc_test[,"my_obs",with=F],1,function(x){ifelse(x=="YES", 1, 0)}))]
            sc_test <- sc_test[,':='("my_pred"=apply(sc_test[,"my_pred",with=F],1,function(x){ifelse(x=="NO"|| is.na(x), 0, 1)}))]
            score <- roc(response= unlist(sc_test[,my_obs]), predictor=unlist(sc_test[,my_pred]))
            score <- round(as.numeric(score$auc), digits=3)
        }
        time <- round(begin - proc.time(), digits=0)
        big_result_table <<- rbindlist(list(big_result_table, data.table(my_case, my_algo, my_metric, score, time[3])))
        return(list("score"=score, "temps"=time[3], "cas"=my_case, "algo"=my_algo, "metric"=my_metric))
    }
}


# begin <- Sys.time()
# time <- as.numeric(-(begin - Sys.time()))
## Modèle à tester
############## TREE #################
### CART
algo <- "rpart2"
eGrid <- expand.grid(.maxdepth = c(10, 20, 30))
simple_tree_10_qq <- ml_cl(data_10, 10, algo, eGrid, "qq_metric", "prob")
simple_tree_30_qq <- ml_cl(data_30, 30, algo, eGrid, "qq_metric", "prob")
simple_tree_60_qq <- ml_cl(data_60, 60, algo, eGrid, "qq_metric", "prob")
simple_tree_10_acc <- ml_cl(data_10, 10, algo, eGrid, "Accuracy", "prob")
simple_tree_30_acc <- ml_cl(data_30, 30, algo, eGrid, "Accuracy", "prob")
simple_tree_60_acc <- ml_cl(data_60, 60, algo, eGrid, "Accuracy", "prob")
simple_tree_10_kap <- ml_cl(data_10, 10, algo, eGrid, "Kappa", "prob")
simple_tree_30_kap <- ml_cl(data_30, 30, algo, eGrid, "Kappa", "prob")
simple_tree_60_kap <- ml_cl(data_60, 60, algo, eGrid, "Kappa", "prob")
simple_tree_10_my_ROC <- ml_cl(data_10, 10, algo, eGrid, "my_ROC", "prob")
simple_tree_30_my_ROC <- ml_cl(data_30, 30, algo, eGrid, "my_ROC", "prob")
simple_tree_60_my_ROC <- ml_cl(data_60, 60, algo, eGrid, "my_ROC", "prob")

### Bagged CART
algo <- "treebag"
eGrid <- NULL
Bagged_tree_10_qq <- ml_cl(data_10, 10, algo, eGrid, "qq_metric", "prob")
Bagged_tree_30_qq <- ml_cl(data_30, 30, algo, eGrid, "qq_metric", "prob")
Bagged_tree_60_qq <- ml_cl(data_60, 60, algo, eGrid, "qq_metric", "prob")
Bagged_tree_10_acc <- ml_cl(data_10, 10, algo, eGrid, "Accuracy", "prob")
Bagged_tree_30_acc <- ml_cl(data_30, 30, algo, eGrid, "Accuracy", "prob")
Bagged_tree_60_acc <- ml_cl(data_60, 60, algo, eGrid, "Accuracy", "prob")
Bagged_tree_10_kap <- ml_cl(data_10, 10, algo, eGrid, "Kappa", "prob")
Bagged_tree_30_kap <- ml_cl(data_30, 30, algo, eGrid, "Kappa", "prob")
Bagged_tree_60_kap <- ml_cl(data_60, 60, algo, eGrid, "Kappa", "prob")
Bagged_tree_10_my_ROC <- ml_cl(data_10, 10, algo, eGrid, "my_ROC", "prob")
Bagged_tree_30_my_ROC <- ml_cl(data_30, 30, algo, eGrid, "my_ROC", "prob")
Bagged_tree_60_my_ROC <- ml_cl(data_60, 60, algo, eGrid, "my_ROC", "prob")

### Boosted Tree
algo <- "blackboost"
eGrid <- expand.grid(.mstop = c(50, 100, 150, 200), .maxdepth = c(10, 20, 30))
Boosted_tree_10_qq <- ml_cl(data_10, 10, algo, eGrid, "qq_metric", "prob")
Boosted_tree_30_qq <- ml_cl(data_30, 30, algo, eGrid, "qq_metric", "prob")
Boosted_tree_60_qq <- ml_cl(data_60, 60, algo, eGrid, "qq_metric", "prob")
Boosted_tree_10_acc <- ml_cl(data_10, 10, algo, eGrid, "Accuracy", "prob")
Boosted_tree_30_acc <- ml_cl(data_30, 30, algo, eGrid, "Accuracy", "prob")
Boosted_tree_60_acc <- ml_cl(data_60, 60, algo, eGrid, "Accuracy", "prob")
Boosted_tree_10_kap <- ml_cl(data_10, 10, algo, eGrid, "Kappa", "prob")
Boosted_tree_30_kap <- ml_cl(data_30, 30, algo, eGrid, "Kappa", "prob")
Boosted_tree_60_kap <- ml_cl(data_60, 60, algo, eGrid, "Kappa", "prob")
Boosted_tree_10_my_ROC <- ml_cl(data_10, 10, algo, eGrid, "my_ROC", "prob")
Boosted_tree_30_my_ROC <- ml_cl(data_30, 30, algo, eGrid, "my_ROC", "prob")
Boosted_tree_60_my_ROC <- ml_cl(data_60, 60, algo, eGrid, "my_ROC", "prob")

### Random Forest
algo <- "ranger"
eGrid <- expand.grid(.mtry = c(1,2,3,4,5,6,7))
rf_tree_10_qq <- ml_cl(data_10, 10, algo, eGrid, "qq_metric", "prob")
rf_tree_30_qq <- ml_cl(data_30, 30, algo, eGrid, "qq_metric", "prob")
rf_tree_60_qq <- ml_cl(data_60, 60, algo, eGrid, "qq_metric", "prob")
rf_tree_10_acc <- ml_cl(data_10, 10, algo, eGrid, "Accuracy", "prob")
rf_tree_30_acc <- ml_cl(data_30, 30, algo, eGrid, "Accuracy", "prob")
rf_tree_60_acc <- ml_cl(data_60, 60, algo, eGrid, "Accuracy", "prob")
rf_tree_10_kap <- ml_cl(data_10, 10, algo, eGrid, "Kappa", "prob")
rf_tree_30_kap <- ml_cl(data_30, 30, algo, eGrid, "Kappa", "prob")
rf_tree_60_kap <- ml_cl(data_60, 60, algo, eGrid, "Kappa", "prob")
rf_tree_10_my_ROC <- ml_cl(data_10, 10, algo, eGrid, "my_ROC", "prob")
rf_tree_30_my_ROC <- ml_cl(data_30, 30, algo, eGrid, "my_ROC", "prob")
rf_tree_60_my_ROC <- ml_cl(data_60, 60, algo, eGrid, "my_ROC", "prob")

############## Linear Model #################
## Generalized Linear Model
algo <- "glmnet"
eGrid <- expand.grid(.alpha = 1, .lambda = 0)
glm_10_qq <- ml_cl(data_10, 10, algo, eGrid, "qq_metric", "prob")
glm_30_qq <- ml_cl(data_30, 30, algo, eGrid, "qq_metric", "prob")
glm_60_qq <- ml_cl(data_60, 60, algo, eGrid, "qq_metric", "prob")
glm_10_acc <- ml_cl(data_10, 10, algo, eGrid, "Accuracy", "prob")
glm_30_acc <- ml_cl(data_30, 30, algo, eGrid, "Accuracy", "prob")
glm_60_acc <- ml_cl(data_60, 60, algo, eGrid, "Accuracy", "prob")
glm_10_kap <- ml_cl(data_10, 10, algo, eGrid, "Kappa", "prob")
glm_30_kap <- ml_cl(data_30, 30, algo, eGrid, "Kappa", "prob")
glm_60_kap <- ml_cl(data_60, 60, algo, eGrid, "Kappa", "prob")
glm_10_my_ROC <- ml_cl(data_10, 10, algo, eGrid, "my_ROC", "prob")
glm_30_my_ROC <- ml_cl(data_30, 30, algo, eGrid, "my_ROC", "prob")
glm_60_my_ROC <- ml_cl(data_60, 60, algo, eGrid, "my_ROC", "prob")

## L1 Logistic regression (lasso)
algo <- "glmnet"
eGrid <- expand.grid(.alpha = 1, .lambda = 10^-(-5:5))
Lone_glm_10_qq <- ml_cl(data_10, 10, algo, eGrid, "qq_metric", "prob")
Lone_glm_30_qq <- ml_cl(data_30, 30, algo, eGrid, "qq_metric", "prob")
Lone_glm_60_qq <- ml_cl(data_60, 60, algo, eGrid, "qq_metric", "prob")
Lone_glm_10_acc <- ml_cl(data_10, 10, algo, eGrid, "Accuracy", "prob")
Lone_glm_30_acc <- ml_cl(data_30, 30, algo, eGrid, "Accuracy", "prob")
Lone_glm_60_acc <- ml_cl(data_60, 60, algo, eGrid, "Accuracy", "prob")
Lone_glm_10_kap <- ml_cl(data_10, 10, algo, eGrid, "Kappa", "prob")
Lone_glm_30_kap <- ml_cl(data_30, 30, algo, eGrid, "Kappa", "prob")
Lone_glm_60_kap <- ml_cl(data_60, 60, algo, eGrid, "Kappa", "prob")
Lone_glm_10_my_ROC <- ml_cl(data_10, 10, algo, eGrid, "my_ROC", "prob")
Lone_glm_30_my_ROC <- ml_cl(data_30, 30, algo, eGrid, "my_ROC", "prob")
Lone_glm_60_my_ROC <- ml_cl(data_60, 60, algo, eGrid, "my_ROC", "prob")

## L2 Logistic regression (ridge)
algo <- "glmnet"
eGrid <- expand.grid(.alpha = 0, .lambda = 10^-(-5:5))
Ltwo_glm_10_qq <- ml_cl(data_10, 10, algo, eGrid, "qq_metric", "prob")
Ltwo_glm_30_qq <- ml_cl(data_30, 30, algo, eGrid, "qq_metric", "prob")
Ltwo_glm_60_qq <- ml_cl(data_60, 60, algo, eGrid, "qq_metric", "prob")
Ltwo_glm_10_acc <- ml_cl(data_10, 10, algo, eGrid, "Accuracy", "prob")
Ltwo_glm_30_acc <- ml_cl(data_30, 30, algo, eGrid, "Accuracy", "prob")
Ltwo_glm_60_acc <- ml_cl(data_60, 60, algo, eGrid, "Accuracy", "prob")
Ltwo_glm_10_kap <- ml_cl(data_10, 10, algo, eGrid, "Kappa", "prob")
Ltwo_glm_30_kap <- ml_cl(data_30, 30, algo, eGrid, "Kappa", "prob")
Ltwo_glm_60_kap <- ml_cl(data_60, 60, algo, eGrid, "Kappa", "prob")
Ltwo_glm_10_my_ROC <- ml_cl(data_10, 10, algo, eGrid, "my_ROC", "prob")
Ltwo_glm_30_my_ROC <- ml_cl(data_30, 30, algo, eGrid, "my_ROC", "prob")
Ltwo_glm_60_my_ROC <- ml_cl(data_60, 60, algo, eGrid, "my_ROC", "prob")

## Boosted Logistic Regression
algo <- "LogitBoost"
eGrid <- expand.grid(.nIter = c(50, 100, 150, 200))
b_lg_10_qq <- ml_cl(data_10, 10, algo, eGrid, "qq_metric", "prob")
b_lg_30_qq <- ml_cl(data_30, 30, algo, eGrid, "qq_metric", "prob")
b_lg_60_qq <- ml_cl(data_60, 60, algo, eGrid, "qq_metric", "prob")
b_lg_10_acc <- ml_cl(data_10, 10, algo, eGrid, "Accuracy", "prob")
b_lg_30_acc <- ml_cl(data_30, 30, algo, eGrid, "Accuracy", "prob")
b_lg_60_acc <- ml_cl(data_60, 60, algo, eGrid, "Accuracy", "prob")
b_lg_10_kap <- ml_cl(data_10, 10, algo, eGrid, "Kappa", "prob")
b_lg_30_kap <- ml_cl(data_30, 30, algo, eGrid, "Kappa", "prob")
b_lg_60_kap <- ml_cl(data_60, 60, algo, eGrid, "Kappa", "prob")
b_lg_10_my_ROC <- ml_cl(data_10, 10, algo, eGrid, "my_ROC", "prob")
b_lg_30_my_ROC <- ml_cl(data_30, 30, algo, eGrid, "my_ROC", "prob")
b_lg_60_my_ROC <- ml_cl(data_60, 60, algo, eGrid, "my_ROC", "prob")


############### Plus pmy_ROChes voisins ###########
algo <- "knn"
eGrid <- expand.grid(.k = c(1, 2, 3, 4,5, 6))
knn_10_qq <- ml_cl(data_10, 10, algo, eGrid, "qq_metric", "raw")
knn_30_qq <- ml_cl(data_30, 30, algo, eGrid, "qq_metric", "raw")
knn_60_qq <- ml_cl(data_60, 60, algo, eGrid, "qq_metric", "raw")
knn_10_acc <- ml_cl(data_10, 10, algo, eGrid, "Accuracy", "raw")
knn_30_acc <- ml_cl(data_30, 30, algo, eGrid, "Accuracy", "raw")
knn_60_acc <- ml_cl(data_60, 60, algo, eGrid, "Accuracy", "raw")
knn_10_kap <- ml_cl(data_10, 10, algo, eGrid, "Kappa", "raw")
knn_30_kap <- ml_cl(data_30, 30, algo, eGrid, "Kappa", "raw")
knn_60_kap <- ml_cl(data_60, 60, algo, eGrid, "Kappa", "raw")
knn_10_my_ROC <- ml_cl(data_10, 10, algo, eGrid, "my_ROC", "raw")
knn_30_my_ROC <- ml_cl(data_30, 30, algo, eGrid, "my_ROC", "raw")
knn_60_my_ROC <- ml_cl(data_60, 60, algo, eGrid, "my_ROC", "raw")

############### Réseaux de neurones ###########
algo <- "mlp"
eGrid <- expand.grid(.size = 1:5)
rn_10_qq <- ml_cl(data_10, 10, algo, eGrid, "qq_metric", "prob")
rn_30_qq <- ml_cl(data_30, 30, algo, eGrid, "qq_metric", "prob")
rn_60_qq <- ml_cl(data_60, 60, algo, eGrid, "qq_metric", "prob")
rn_10_acc <- ml_cl(data_10, 10, algo, eGrid, "Accuracy", "prob")
rn_30_acc <- ml_cl(data_30, 30, algo, eGrid, "Accuracy", "prob")
rn_60_acc <- ml_cl(data_60, 60, algo, eGrid, "Accuracy", "prob")
rn_10_kap <- ml_cl(data_10, 10, algo, eGrid, "Kappa", "prob")
rn_30_kap <- ml_cl(data_30, 30, algo, eGrid, "Kappa", "prob")
rn_60_kap <- ml_cl(data_60, 60, algo, eGrid, "Kappa", "prob")
rn_10_my_ROC <- ml_cl(data_10, 10, algo, eGrid, "my_ROC", "prob")
rn_30_my_ROC <- ml_cl(data_30, 30, algo, eGrid, "my_ROC", "prob")
rn_60_my_ROC <- ml_cl(data_60, 60, algo, eGrid, "my_ROC", "prob")

############### SVM ###########
### Linear Kernel
algo <- "svmLinear"
eGrid <- expand.grid(.C = c(10^(-5), 10^(-3), 10, 10^(3), 10^(6), 10^(9)))
svmLinear_10_qq <- ml_cl(data_10, 10, algo, eGrid, "qq_metric", "raw")
svmLinear_30_qq <- ml_cl(data_30, 30, algo, eGrid, "qq_metric", "raw")
svmLinear_60_qq <- ml_cl(data_60, 60, algo, eGrid, "qq_metric", "raw")
svmLinear_10_acc <- ml_cl(data_10, 10, algo, eGrid, "Accuracy", "raw")
svmLinear_30_acc <- ml_cl(data_30, 30, algo, eGrid, "Accuracy", "raw")
svmLinear_60_acc <- ml_cl(data_60, 60, algo, eGrid, "Accuracy", "raw")
svmLinear_10_kap <- ml_cl(data_10, 10, algo, eGrid, "Kappa", "raw")
svmLinear_30_kap <- ml_cl(data_30, 30, algo, eGrid, "Kappa", "raw")
svmLinear_60_kap <- ml_cl(data_60, 60, algo, eGrid, "Kappa", "raw")
svmLinear_10_my_ROC <- ml_cl(data_10, 10, algo, eGrid, "my_ROC", "raw")
svmLinear_30_my_ROC <- ml_cl(data_30, 30, algo, eGrid, "my_ROC", "raw")
svmLinear_60_my_ROC <- ml_cl(data_60, 60, algo, eGrid, "my_ROC", "raw")


### Radial Kernel
algo <- "svmRadial"
eGrid <- expand.grid(.C = c(10^(-5), 10^(-3), 10, 10^(3), 10^(6), 10^(9)), .sigma=c(10^(-5), 10^(-3), 10^(-1), 10))
svmRadial_10_qq <- ml_cl(data_10, 10, algo, eGrid, "qq_metric", "raw")
svmRadial_30_qq <- ml_cl(data_30, 30, algo, eGrid, "qq_metric", "raw")
svmRadial_60_qq <- ml_cl(data_60, 60, algo, eGrid, "qq_metric", "raw")
svmRadial_10_acc <- ml_cl(data_10, 10, algo, eGrid, "Accuracy", "raw")
svmRadial_30_acc <- ml_cl(data_30, 30, algo, eGrid, "Accuracy", "raw")
svmRadial_60_acc <- ml_cl(data_60, 60, algo, eGrid, "Accuracy", "raw")
svmRadial_10_kap <- ml_cl(data_10, 10, algo, eGrid, "Kappa", "raw")
svmRadial_30_kap <- ml_cl(data_30, 30, algo, eGrid, "Kappa", "raw")
svmRadial_60_kap <- ml_cl(data_60, 60, algo, eGrid, "Kappa", "raw")
svmRadial_10_my_ROC <- ml_cl(data_10, 10, algo, eGrid, "my_ROC", "raw")
svmRadial_30_my_ROC <- ml_cl(data_30, 30, algo, eGrid, "my_ROC", "raw")
svmRadial_60_my_ROC <- ml_cl(data_60, 60, algo, eGrid, "my_ROC", "raw")

### Bayesian Generalized Linear Model
algo <- "bayesglm"
eGrid <- NULL
arm_10_qq <- ml_cl(data_10, 10, algo, eGrid, "qq_metric", "prob")
arm_30_qq <- ml_cl(data_30, 30, algo, eGrid, "qq_metric", "prob")
arm_60_qq <- ml_cl(data_60, 60, algo, eGrid, "qq_metric", "prob")
arm_10_acc <- ml_cl(data_10, 10, algo, eGrid, "Accuracy", "prob")
arm_30_acc <- ml_cl(data_30, 30, algo, eGrid, "Accuracy", "prob")
arm_60_acc <- ml_cl(data_60, 60, algo, eGrid, "Accuracy", "prob")
arm_10_kap <- ml_cl(data_10, 10, algo, eGrid, "Kappa", "prob")
arm_30_kap <- ml_cl(data_30, 30, algo, eGrid, "Kappa", "prob")
arm_60_kap <- ml_cl(data_60, 60, algo, eGrid, "Kappa", "prob")
arm_10_my_ROC <- ml_cl(data_10, 10, algo, eGrid, "my_ROC", "prob")
arm_30_my_ROC <- ml_cl(data_30, 30, algo, eGrid, "my_ROC", "prob")
arm_60_my_ROC <- ml_cl(data_60, 60, algo, eGrid, "my_ROC", "prob")

save.image(file="classification_3.RData")
