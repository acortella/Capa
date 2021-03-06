## Initisation du chemin du dossier contenant les donn�es
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

## Chargement des donn�es
limit_data <- fread(paste0(path, "limit_data.csv"))
used_data <- fread(paste0(path, "used_data.csv"))


######################################################################
###################################### Nettoyage des donn�es #########
######################################################################

#### Etape 1
#### Supprime les bases de donn�es qui ont au moins une fois une 
#### volum�trie utilis�e strictement sup�rieure � la volum�tire
#### limite pour une date donn�e.
vol_sup_lim <- function(my_data_limit, my_data_used){
    matrix_my_data_limit <- as.matrix(my_data_limit[,3:ncol(my_data_limit),with=F])
    matrix_my_data_used <- as.matrix(my_data_used[,3:ncol(my_data_used),with=F])
    diff_table <- as.data.table(matrix_my_data_used/matrix_my_data_limit)
    
    diff_table <- as.data.table(matrix_my_data_limit-matrix_my_data_used)
    diff_table <- cbind(my_data_limit[,2,with=F],
        diff_table[, ':='("sup"=apply(diff_table[,1:ncol(diff_table),
                     with=F],1,function(x){ifelse(any(x<0),0,1)}))])
    ind_delete <- diff_table[sup==0][,ind]
    
    my_data_limit <- my_data_limit[! ind %in% ind_delete]
    my_data_used <- my_data_used[! ind %in% ind_delete]
    return(rbindlist(list(my_data_limit, my_data_used)))
}
temp <- vol_sup_lim(limit_data, used_data)
limit_data <- temp[1:(nrow(temp)/2),]
used_data <- temp[((nrow(temp)/2)+1):nrow(temp),]
## Print le nombre de base de donn�es restantes
print(c("sup_lim", nrow(limit_data)))


#### Etape 2
#### Analyse du nombre de bases de donn�es avec les volum�tries 
##enti�res en fonction du temps et suppression des dates
#### dont la stabilit� de collecte n'est pas assur�e.
vol_collecte <- function(my_data){
    days <- NULL
    nb_entire <- NULL
    for(i in 3:1460){
      # print(c(i, nrow(na.omit(limit_data[,i:1460,with=F]))))
      days <- c(days, i)
      nb_entire <- c(nb_entire, unlist(nrow(na.omit(my_data[,i:1460,with=F]))))
    }
    p_days_data <- data.frame("days"=days, "nb_bases"=nb_entire)
    p_collecte <- plot_ly(x=~days, y=~nb_entire, data=p_days_data, mode="lines", type="scatter")
    study_table <- data.table("day"=days, "nb"=nb_entire)
    first_day <- unlist(study_table[nb!=0][1])[1]
    return(list(p_collecte, first_day))
}
temp <- vol_collecte(limit_data)
plot_db_entire <- temp[1]
limit_data <- limit_data[,-c(3:(unlist(temp[2])-1)),with=F]
used_data <- used_data[,-c(3:(unlist(temp[2])-1)),with=F]


#### Etape 3
#### On fixe un test set de 90 jours sur 100 jours et on garde 60 jours
#### pour les pr�visions, on garde aussi au moins 5 jours sur
#### le training set pour faire de la cross-validation avec la regression
#### On garde aussi 90 jours pour pouvoir calculer les �quations � 90 inconnues
#### avec la probl�matique de r�gression, on garde aussi 90 en plus car l'on ne
#### va prendre en entr�es que des valeurs diff�rentes de NA
#### Suppression des bases de donn�es qui n'ont pas assez de volum�trie
#### connue pour faire partie de l'�tude
vol_minimum <- function(my_data, d_prev=60, d_feat_test=90,
                        d_nb_test=99, d_training_cv=5, d_inconnues=90, d_na=90){
    d_min <- d_prev + d_feat_test + d_nb_test + d_training_cv + d_inconnues + d_na
    rel_ind_l <- unlist(na.omit(my_data[,-c(1, 3:(ncol(my_data)-d_min+1)),
                                        with=F])[,ind])
    my_data <- my_data[ind %in% rel_ind_l]
    return(my_data)
}
limit_data <- vol_minimum(limit_data, d_prev=60, d_feat_test=90, d_nb_test=100)
used_data <- vol_minimum(used_data, d_prev=60, d_feat_test=90, d_nb_test=100)
print(c("vol_min", nrow(limit_data)))


#### Etape 4
#### On choisit de ne pas travailler avec des valeurs imput�es, 
#### on supprime donc les bases de donn�es qui ont des
#### valeurs manquantes. On consid�re en effet que ces bases de donn�es
#### ont rencontr�es des probl�mes dans la collecte
#### et que les autres donn�es sont potentiellement source d'erreurs
vol_entire <- function(my_data){
  my_data <- my_data[, ':='("entire"=apply(my_data[,3:ncol(my_data),with=F],1,function(x){
    x <- unlist(x)
    i <- 1
    while(is.na(x[i])){
      i <- i + 1
    }
    if(i==(length(x)+1)){
      return(0)
    }else{
      x <- x[i:length(x)] 
      ifelse(length(x)==length(x[!is.na(x)]),return(1),return(0))
    }
  }))]
  my_data <- my_data[entire==1][,-c("entire"),with=F]
}
limit_data <- vol_entire(limit_data)
used_data <- vol_entire(used_data)
print(c("vol_ent", nrow(limit_data)))
#### V�rifie si les num�ros d'indices correspondent entre le set de
#### donn��s de volum�trie utilis�e et de volum�trie limite
all((unlist(limit_data[,ind])==unlist(used_data[,ind]))==TRUE)


#### Etape 5
#### V�rifie la r�elle valeur de la limite (valeur pour laquelle on
#### consid�re que on a d�pass� la limite)
real_l <- function(my_data_limit, my_data_used){
    matrix_limit_data <- as.matrix(my_data_limit[,3:ncol(my_data_limit),with=F])
    matrix_used_data <- as.matrix(my_data_used[,3:ncol(my_data_used),with=F])
    diff_table <- as.data.table(matrix_used_data/matrix_limit_data)
    temp_table <- NULL
    for(i in 1:1061){
        loop_table <- na.omit(diff_table[,i,with=F])
        names(loop_table) <- "V1"
        temp_table<- rbindlist(list(temp_table, loop_table))
    }
    temp_vect <- sort(unname(unlist(temp_table[,1,with=F])), decreasing = FALSE)
    return(data.table("value"=temp_vect))
}
temp <- real_l(limit_data, used_data)

## Plot 0 � 1
temp <- temp[value>=0][value<=1]
temp <- cbind(data.table("nb"=1:nrow(temp)), temp)
temp_vect <- temp[, ':='("modulo"=apply(temp[,c("nb"),with=F], 1, 
       function(x){ifelse(x%%100==0,1,0)}))][modulo==1]
p_0_to_1 <- plot_ly(x=temp_vect[,nb], y=temp_vect[,value], 
     xaxis="nb_entry", yaxis="threshold", mode="lines", type="scatter")

## Plot 0.95 � 1
temp <- temp[value>=0.95]
p_0_95_to_1 <- plot_ly(x=temp[,nb], y=temp[,value], 
    xaxis="nb_entry", yaxis="threshold", mode="lines", type="scatter")
s1 <- temp[1,][,nb]:5920964
s2 <- 5920964:temp[nrow(temp),][,nb]
m1 <- lm(unlist(temp[nb>=s1[1]][nb<=s1[length(s1)]][,value])~s1)
m2 <- lm(unlist(temp[nb>=s2[1]][nb<=s2[length(s2)]][,value])~s2)

#  Aggrandit les intervalles pour la visibilit�e sur la courbe
s1 <- c(s1, (s1[length(s1)]+1):(s1[length(s1)]+2000))
s2 <- c((s2[1]-4000):(s2[1]-1),s2)
temp <- temp[,':='("pred_1"=apply(temp[,"nb",with=F],1,function(x){
        ifelse(x %in% s1, m1$coefficients[2]*x+m1$coefficients[1],NA)}),
                   "pred_2"=apply(temp[,"nb",with=F],1,function(x){
        ifelse(x %in% s2, m2$coefficients[2]*x+m2$coefficients[1],NA)})
                   )]
p_0_95_to_1_bis <- plot_ly(x=temp[,nb], y=temp[,value], name ="courbe_ref",
            xaxis="nb_entry", yaxis="threshold", mode="lines", type="scatter") %>%
            add_trace(y=temp[,pred_1], name = 'tangent 1',mode = 'lines') %>%
            add_trace(y =temp[,pred_2], name = 'tangent 2',mode = 'lines')

#### Le seuil est dertermin� comme l'ordonn� de la courbe de r�f�rence 
#### dont l'abscisse est l'abscisse d'intersection des deux tangentes
abs <- (m2$coefficients[1]-m1$coefficients[1])/(m1$coefficients[2]-m2$coefficients[2])

## Le nouveau seuil est donc
threshold <- round(temp[,':='("eucl_d"=apply(temp[,"nb",with=F],1,function(x){
                  sqrt((x-abs)^2)}))][which.min(eucl_d)][,value],digits=4)

## Transformation de limit en threshold*limit
limit_data <- cbind(limit_data[,1:2,with=F], 
        as.data.table(as.matrix(limit_data[,3:ncol(limit_data),with=F])*threshold))

## Supprime les donn�es inutiles pour la suite
remove(abs, s1, s2, m1, m2, temp_vect, temp)

## Sauvegarde des donn�es
save.image(file="data_cleaning.RData")


