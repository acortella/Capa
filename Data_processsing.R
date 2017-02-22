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

## Chargement des données
load("data_cleaning.RData")


#################################################################
###################################### Data processing #########
################################################################
#### Cette function permet, A partir d'un data.table d'historique de données, de calculer les maximums successifs pour les 3 différents
#### cas : 10 jours, 30 jours et 60 jours
get_max <- function(my_dataset){
    # ##
    # my_dataset <- g_used_data
    # ##
    ref <- my_dataset[,1:2,with=F]
    my_dataset <- copy(my_dataset[,3:length(my_dataset),with=F])
    se_dates <- names(my_dataset)
    
    max_10 <- matrix(data=0, nrow=nrow(my_dataset), ncol=(ncol(my_dataset)-60))
    max_30 <- matrix(data=0, nrow=nrow(my_dataset), ncol=(ncol(my_dataset)-60))
    max_60 <- matrix(data=0, nrow=nrow(my_dataset), ncol=(ncol(my_dataset)-60))
    
    for(i in 1:nrow(my_dataset)){
        # print(i)
        temp_vect <- unlist(my_dataset[i])
        temp_10 <- c(runmax(temp_vect, 10, endrule="NA", align="left")[2:length(temp_vect)], NA)
        temp_30 <- c(runmax(temp_vect, 30, endrule="NA", align="left")[2:length(temp_vect)], NA)
        temp_60 <- c(runmax(temp_vect, 60, endrule="NA", align="left")[2:length(temp_vect)], NA)
        
        max_10[i,] <- temp_10[1:(length(temp_10)-60)]
        max_30[i,] <- temp_30[1:(length(temp_10)-60)]
        max_60[i,] <- temp_60[1:(length(temp_10)-60)]
    }
    max_10 <- as.data.table(max_10)
    max_30 <- as.data.table(max_30)
    max_60 <- as.data.table(max_60)
    names(max_10) <- names(max_30) <- names(max_60) <- se_dates[1:(length(se_dates)-60)]
    max_10 <- cbind(ref, max_10)
    max_30 <- cbind(ref, max_30)
    max_60 <- cbind(ref, max_60)
    
    return(list(max_10, max_30, max_60))
}

all_max <- get_max(used_data)
max_10 <- as.data.table(all_max[1])
max_30 <- as.data.table(all_max[2])
max_60 <- as.data.table(all_max[3])
remove(all_max)

### Generate the training and the test set
real_used <- used_data[,(ncol(used_data)-59):ncol(used_data), with=F]

training_used <- cbind(data.table("ind"=1:nrow(used_data)), used_data[,3:(ncol(used_data)-60-100-89), with=F])
training_limit <- cbind(data.table("ind"=1:nrow(used_data)), limit_data[,3:(ncol(limit_data)-60-100-89), with=F])
training_max_10 <- cbind(data.table("ind"=1:nrow(used_data)), max_10[,3:(ncol(max_10)-100-89), with=F])
training_max_30 <- cbind(data.table("ind"=1:nrow(used_data)), max_30[,3:(ncol(max_30)-100-89), with=F])
training_max_60 <- cbind(data.table("ind"=1:nrow(used_data)), max_60[,3:(ncol(max_60)-100-89), with=F])

test_used <- cbind(data.table("ind"=1:nrow(used_data)), used_data[,(ncol(used_data)-59-100-89):(ncol(used_data)-60), with=F])
test_limit <-  cbind(data.table("ind"=1:nrow(used_data)), limit_data[,(ncol(limit_data)-59-100-89):(ncol(limit_data)-60), with=F])
test_max_10 <-  cbind(data.table("ind"=1:nrow(used_data)), max_10[,(ncol(max_10)-100-88):ncol(max_10), with=F])
test_max_30 <-  cbind(data.table("ind"=1:nrow(used_data)), max_30[,(ncol(max_30)-100-88):ncol(max_30), with=F])
test_max_60 <-  cbind(data.table("ind"=1:nrow(used_data)), max_60[,(ncol(max_60)-100-88):ncol(max_60), with=F])

# ## Supprime les données non-utiles
# remove(max_10, max_30, max_60, real_used, limit_data, used_data, diff_table)

## Fonctions générant le training set et le test set qui sera pris en entrée pour la classification et la régression, supprime les NA
get_set <- function(my_dataset_used, my_dataset_limit, my_dataset_max_10, my_dataset_max_30, my_dataset_max_60){
    # my_dataset_used <- training_used
    # my_dataset_limit <- training_limit
    # my_dataset_max_10 <- training_max_10
    # my_dataset_max_30 <- training_max_30
    # my_dataset_max_60 <- training_max_60
    
    training_cl <- NULL
    for(i in 1:(length(my_dataset_used)-90)){
        print(i)
        temp_ref <- my_dataset_used[,"ind",with=F]
        temp_used <- my_dataset_used[,(i+1):(90+i), with=F]
        temp_limit <- my_dataset_limit[,(90+i), with=F]
        temp_max_10 <- my_dataset_max_10[,(90+i), with=F]
        temp_max_30 <- my_dataset_max_30[,(90+i), with=F]
        temp_max_60 <- my_dataset_max_60[,(90+i), with=F]
        
        temp_training_cl <- cbind(temp_ref, temp_used, temp_limit, temp_max_10, temp_max_30, temp_max_60)
        temp_training_cl <- na.omit(temp_training_cl)
        training_cl <- rbindlist(list(training_cl, temp_training_cl))
    }
    setnames(training_cl, c("ind", paste0("x",1:90), "limit", "y10", "y30", "y60"))
    return(training_cl)
}

training_set <- get_set(training_used, training_limit, training_max_10, training_max_30, training_max_60) 
test_set <- get_set(test_used, test_limit, test_max_10, test_max_30, test_max_60) 

keep <- unique(training_set[,ind])
test_set <- test_set[ind %in% keep]
# limit_data <- limit_data[ind %in% keep]
# max_60 <- max_60[ind %in% keep]

## Vérifie que les id correspondent
length(sort(unique(training_set[,ind])))==length(sort(unique(test_set[,ind]), decreasing=F))


#####################################################################################################################
###################################### Echantillonnage ##############################################################
#####################################################################################################################

## Choix du nombe de bases de données retenues
p <- 0.5
s <- 95/100
e <- 5/100
t <- 1.96
i_ech <- (t*t*p*(1-p))/(e*e)

## Population ajustée
nb_bd <- nrow(test_set)/100
n_ech <- i_ech/(1 + (i_ech-1)/nb_bd)
n_ech <- ceiling(n_ech)

## Génération de nombre aléatoires
set.seed(1)
rand_row <- base::sample(1:nrow(limit_data), n_ech, replace = FALSE)

## Détermination des bases de données restantes et suppressions des autres
temp <- cbind(data.table("rel_ind"=1:nrow(limit_data)), limit_data)
temp <- temp [rel_ind %in% rand_row]
temp  <- unlist(temp [,ind])

test_set <-  test_set[ind %in% temp]
training_set <- training_set[ind %in% temp]


## Supprime les données non-utilisées
remove(test_limit, test_used, test_max_10, test_max_30, test_max_60, training_max_10, training_max_30,
       training_max_60, training_limit, training_used, e, n_ech, p, rand_row, s, t, temp, limit_data,
       max_10, max_30, max_60, used_data, keep, real_used)

## Sauvegarde des données
save.image(file="data_processing.RData")
