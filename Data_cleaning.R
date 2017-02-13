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

## Chargement des données
limit_data <- fread(paste0(path, "limit_data.csv"))
used_data <- fread(paste0(path, "used_data.csv"))


###########################################################################################################################
###################################### Nettoyage des données ##############################################################
###########################################################################################################################

#### Etape 1
#### Supprime les bases de données qui ont au moins une fois une volumétrie utilisée strictement supérieure à la volumétire
#### limite pour une date donnée.
vol_sup_lim <- function(my_data_limit, my_data_used){
    matrix_my_data_limit <- as.matrix(my_data_limit[,3:ncol(my_data_limit),with=F])
    matrix_my_data_used <- as.matrix(my_data_used[,3:ncol(my_data_used),with=F])
    diff_table <- as.data.table(matrix_my_data_used/matrix_my_data_limit)
    
    diff_table <- as.data.table(matrix_my_data_limit-matrix_my_data_used)
    diff_table <- cbind(my_data_limit[,2,with=F],
                        diff_table[, ':='("sup"=apply(diff_table[,1:ncol(diff_table),with=F],1,function(x){ifelse(any(x<0),0,1)}))])
    ind_delete <- diff_table[sup==0][,ind]
    
    my_data_limit <- my_data_limit[! ind %in% ind_delete]
    my_data_used <- my_data_used[! ind %in% ind_delete]
    return(rbindlist(list(my_data_limit, my_data_used)))
}
temp <- vol_sup_lim(limit_data, used_data)
limit_data <- temp[1:(nrow(temp)/2),]
used_data <- temp[((nrow(temp)/2)+1):nrow(temp),]
## Print le nombre de base de données restantes
print(c("sup_lim", nrow(limit_data)))


#### Etape 2
#### Analyse du nombre de bases de données avec les volumétries entières en fonction du temps et suppression des dates
#### dont la stabilité de collecte n'est pas assurée.
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
#### On choisit de ne pas travailler avec des valeurs imputées, on supprime donc les bases de données qui ont des
#### valeurs manquantes. On considère en effet que ces bases de données ont rencontrées des problèmes dans la collecte
#### et que les autres données sont potentiellement source d'erreurs
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
#### Vérifie si les numéros d'indices correspondent entre le set de donnéés de volumétrie utilisée et de volumétrie limite
all((unlist(limit_data[,ind])==unlist(used_data[,ind]))==TRUE)


#### Etape 4
#### On fixe un test set de 90 jours sur 100 jours et on garde 60 jours pour les prévisions
#### Suppression des bases de données qui n'ont pas assez de volumétrie connue pour faire partie de l'étude
vol_minimum <- function(my_data, d_prev=60, d_feat_test=90, d_nb_test=100){
  d_min <- d_prev + d_feat_test + d_nb_test
  rel_ind_l <- unlist(na.omit(my_data[,-c(1, 3:(ncol(my_data)-d_min+1)),with=F])[,ind])
  my_data <- my_data[ind %in% rel_ind_l]
  return(my_data)
}
limit_data <- vol_minimum(limit_data, d_prev=60, d_feat_test=90, d_nb_test=100)
used_data <- vol_minimum(used_data, d_prev=60, d_feat_test=90, d_nb_test=100)
print(c("vol_min", nrow(limit_data)))


#### Etape 5
#### Vérifie la réelle valeur de la limite (valeur pour laquelle on considère que on a dépassé la limite)
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

## Plot 0 à 1
temp <- temp[value>=0][value<=1]
temp <- cbind(data.table("nb"=1:nrow(temp)), temp)
temp_vect <- temp[, ':='("modulo"=apply(temp[,c("nb"),with=F], 1, function(x){ifelse(x%%100==0,1,0)}))][modulo==1]
p_0_to_1 <- plot_ly(x=temp_vect[,nb], y=temp_vect[,value], xaxis="nb_entry", yaxis="threshold", mode="lines", type="scatter")

## Plot 0.95 à 1
temp <- temp[value>=0.95]
p_0_95_to_1 <- plot_ly(x=temp[,nb], y=temp[,value], xaxis="nb_entry", yaxis="threshold", mode="lines", type="scatter")
s1 <- temp[1,][,nb]:6349220
s2 <- 6349220:temp[nrow(temp),][,nb]
m1 <- lm(unlist(temp[nb>=s1[1]][nb<=s1[length(s1)]][,value])~s1)
m2 <- lm(unlist(temp[nb>=s2[1]][nb<=s2[length(s2)]][,value])~s2)

#  Aggrandit les intervalles pour la visibilitée sur la courbe
s1 <- c(s1, (s1[length(s1)]+1):(s1[length(s1)]+2000))
s2 <- c((s2[1]-4000):(s2[1]-1),s2)
temp <- temp[,':='("pred_1"=apply(temp[,"nb",with=F],1,function(x){
                          ifelse(x %in% s1, m1$coefficients[2]*x+m1$coefficients[1],NA)}),
                   "pred_2"=apply(temp[,"nb",with=F],1,function(x){
                          ifelse(x %in% s2, m2$coefficients[2]*x+m2$coefficients[1],NA)})
                   )]
p_0_95_to_1 <- plot_ly(x=temp[,nb], y=temp[,value], name ="courbe_ref", xaxis="nb_entry", yaxis="threshold", mode="lines", type="scatter") %>%
                  add_trace(y=temp[,pred_1], name = 'tangent 1',mode = 'lines') %>%
                  add_trace(y =temp[,pred_2], name = 'tangent 2',mode = 'lines')

## Le seuil est derterminé comme l'ordonné de la courbe de référence dont l'abscisse est l'abscisse d'intersection des deux tangentes
abs <- (m2$coefficients[1]-m1$coefficients[1])/(m1$coefficients[2]-m2$coefficients[2])

## Le nouveau seuil est donc
threshold <- temp[,':='("eucl_d"=apply(temp[,"nb",with=F],1,function(x){
                  sqrt((x-abs)^2)}))][which.min(eucl_d)][,value]

## Transformation de limit en threshold*limit
limit_data <- cbind(limit_data[,1:2,with=F], as.data.table(as.matrix(limit_data[,3:ncol(limit_data),with=F])*threshold))

## Supprime les données inutiles pour la suite
remove(abs, s1, s2, m1, m2, temp_vect, temp)

## Sauvegarde des données
save.image(file="data_cleaning.RData")















# Copy échantillonnage non nécessaire


# #### Etude du nombre de date de rupture pour le cas 10 de façon à supprimer les bases de données in-interessantes
# matrix_limit_data <- as.matrix(limit_data[,3:ncol(limit_data),with=F])
# ncol(matrix_limit_data)
# ncol(max_60[,3:ncol(max_60),with=F])
# matrix_max <- cbind(as.matrix(max_60[,3:ncol(max_60),with=F]), matrix(0, nrow=nrow(max_60), ncol=60))
# diff_table <- as.data.table(matrix_limit_data-matrix_max)
# remove(matrix_limit_data, matrix_max, matrix_used_data, p_days_data)
# diff_table <- cbind(limit_data[,2,with=F],
#                     # diff_table[, ':='("sup"=apply(diff_table[,1:ncol(diff_table),with=F],1,function(x){ifelse(any(x<=0),0,1)}))])[sup==1]
#                     diff_table[, ':='("sup"=apply(diff_table[,1:ncol(diff_table),with=F],1,function(x){
#                       x <- unlist(x)
#                       length(x[which(x<=0)])}))])
# 
# table(diff_table[,sup])
# diff_table <- diff_table[sup>0]

# n_0 <- nrow(diff_table[sup<=0])
# n_1 <- nrow(diff_table[sup>0])
# 
# # p <- n_1 / nrow(diff_table)
# s <- 95/100
# e <- 5/100
# t <- 1.96
# 
# # n_ech <- (t*t*p*(1-p))/(e*e)
# 
# p <- 0.5
n_ech <- (t*t*p*(1-p))/(e*e)
n_ech <- ceiling(n_ech)

set.seed(1)
rand_row <- base::sample(1:nrow(diff_table), n_ech, replace = FALSE)
# rand_row <- base::sample(1:nrow(diff_table), n_ech, replace = FALSE) changement
diff_table <- cbind(data.table("rel_ind"=1:nrow(diff_table)), diff_table)
diff_table <- diff_table[rel_ind %in% rand_row]

kept_row <-  unlist(diff_table[,ind])

limit_data <-  limit_data[ind %in% kept_row]
used_data <- used_data[ind %in% kept_row]
max_10 <- max_10[ind %in% kept_row]
max_30 <- max_30[ind %in% kept_row]
max_60 <- max_60[ind %in% kept_row]