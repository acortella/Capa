setwd("C:/Users/A/Desktop")

library(data.table)
library(caTools)
library(plotly)

## Chargement des données
limit_data <- fread("C:/Users/A/Desktop/CP_data/limit_data.csv")
used_data <- fread("C:/Users/A/Desktop/CP_data/used_data.csv")

## Data processing
#### Delete data which have at least one used > one limit
matrix_limit_data <- as.matrix(limit_data[,3:ncol(limit_data),with=F])
matrix_used_data <- as.matrix(used_data[,3:ncol(used_data),with=F])
diff_table <- as.data.table(matrix_used_data/matrix_limit_data)

diff_table <- as.data.table(matrix_limit_data-matrix_used_data)
diff_table <- cbind(limit_data[,2,with=F],
                    diff_table[, ':='("sup"=apply(diff_table[,1:ncol(diff_table),with=F],1,function(x){ifelse(any(x<0),0,1)}))])

ind_delete <- diff_table[sup==0][,ind]
length(ind_delete)

limit_data <- limit_data[! ind %in% ind_delete]
used_data <- used_data[! ind %in% ind_delete]
remove(diff_table, ind_delete, matrix_limit_data, matrix_used_data)
print(c("sup_lim", nrow(limit_data)))


# days <- NULL
# nb_entire <- NULL
# for(i in 3:1460){
#     # print(c(i, nrow(na.omit(limit_data[,i:1460,with=F]))))
#     days <- c(days, i)
#     nb_entire <- c(nb_entire, unlist(nrow(na.omit(limit_data[,i:1460,with=F]))))
# }
# p_days_data <- data.frame("days"=days, "nb_bases"=nb_entire)
# p_days <- plot_ly(x=~days, y=~nb_entire, data=p_days_data)
# study_table <- data.table("day"=days, "nb"=nb_entire)
# first_day <- unlist(study_table[nb!=0][1])[1]
# print(c("stability", nrow(limit_data)))
# 
# limit_data <- limit_data[,-c(3:399),with=F]
# used_data <- used_data[,-c(3:399),with=F]
# remove(study_table)


## On fixe un test set de 90 jours sur 100 jours et on garde 60 jours pour les prévisions
rel_ind_l <- unlist(na.omit(limit_data[,-c(1, 3:814),with=F])[,ind])
rel_ind_u <- unlist(na.omit(used_data[,-c(1, 3:814),with=F])[,ind])
limit_data <- limit_data[ind %in% rel_ind_l]
used_data <- used_data[ind %in% rel_ind_u]
print(c("pred_test", nrow(limit_data)))

## Supprime les bases de données qui ont des valeurs manquantes
limit_data <- limit_data[, ':='("entire"=apply(limit_data[,3:ncol(limit_data),with=F],1,function(x){
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
limit_data <- limit_data[entire==1][,-c("entire"),with=F]

used_data <- used_data[, ':='("entire"=apply(used_data[,3:ncol(used_data),with=F],1,function(x){
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
used_data <- used_data[entire==1][,-c("entire"),with=F]

## Vérifie si les numéros d'indices correspondent
all((unlist(limit_data[,ind])==unlist(used_data[,ind]))==TRUE)



matrix_limit_data <- as.matrix(limit_data[,3:ncol(limit_data),with=F])
matrix_used_data <- as.matrix(used_data[,3:ncol(used_data),with=F])
diff_table <- as.data.table(matrix_used_data/matrix_limit_data)

temp_table <- NULL
for(i in 1:1061){
  loop_table <- na.omit(diff_table[,i,with=F])
  names(loop_table) <- "V1"
  temp_table<- rbindlist(list(temp_table, loop_table))
}
temp_vect <- sort(unname(unlist(temp_table[,1,with=F])), decreasing = FALSE)

## Plot 0 à 1
temp_vect <- temp_vect[which(temp_vect>=0)]
temp_vect <- temp_vect[which(temp_vect<1)]
names(temp_vect) <- 1:length(temp_vect)
p1 <- plot_ly(x=names(temp_vect), y=temp_vect)

## Plot 0.93 à 1
temp_vect <- temp_vect[which(temp_vect>=0.94)]
temp_vect <- temp_vect[which(temp_vect<1)]
begin <- length(temp_vect[which(temp_vect<0.94)])+1
names(temp_vect) <- begin:length(temp_vect)+begin
p2 <- plot(names(temp_vect), temp_vect)
d <- c(0,0, diff(temp_vect, differences = 2))
names(d) <- begin:length(temp_vect)+begin
p3 <- plot_ly(x = names(d), y = d, type = 'scatter', mode = 'lines')

p4 <- plot_ly(x = 70000:length(d), y = d[70000:length(d)], type = 'scatter', mode = 'lines')

## filter inf to 10 micro
d_filt <- d[70000:length(d)][which(abs(d)>10E-6)]
p4 <- plot_ly(x = 1:length(d_filt), y = d_filt, type = 'scatter', mode = 'lines')












# limit_data <- fread("C:/Users/A/Desktop/CP_data/limit_data.csv")[,3:1460,with=F]
# used_data <- fread("C:/Users/A/Desktop/CP_data/used_data.csv")[,3:1460,with=F]


## Cette function permet, A partir d'un data.table d'historiue de donnÃ©es, de calucler les maximums successifs
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
        print(i)
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

matrix_limit_data <- as.matrix(limit_data[,3:ncol(limit_data),with=F])
matrix_max <- cbind(as.matrix(max_10[,3:ncol(max_60),with=F]), matrix(0, nrow=nrow(max_60), ncol=60))
diff_table <- as.data.table(matrix_limit_data-matrix_max)
diff_table <- cbind(limit_data[,2,with=F],
                    # diff_table[, ':='("sup"=apply(diff_table[,1:ncol(diff_table),with=F],1,function(x){ifelse(any(x<=0),0,1)}))])[sup==1]
                    diff_table[, ':='("sup"=apply(diff_table[,1:ncol(diff_table),with=F],1,function(x){
                          x <- unlist(x)
                          length(x[which(x<=0)])}))])
hist(diff_table[,sup])


## stratification
nb_sup <- nrow()
nb_inf <- 












### Generate the training and the test set
real_used <- used_data[,(ncol(used_data)-59):ncol(used_data), with=F]

test_used <- used_data[,(ncol(used_data)-60-98):(ncol(used_data)-60), with=F]
test_limit <- limit_data[,(ncol(limit_data)-60-98):(ncol(limit_data)-60), with=F]
test_max_10 <- max_10[,(ncol(max_10)-98):ncol(max_10), with=F]
test_max_30 <- max_30[,(ncol(max_30)-98):ncol(max_30), with=F]
test_max_60 <- max_60[,(ncol(max_60)-98):ncol(max_60), with=F]

training_used <- used_data[,1:(ncol(used_data)-60-99), with=F]
training_limit <- limit_data[,1:(ncol(limit_data)-60-99), with=F]
training_max_10 <- max_10[,1:(ncol(max_10)-99), with=F]
training_max_30 <- max_30[,1:(ncol(max_30)-99), with=F]
training_max_60 <- max_60[,1:(ncol(max_60)-99), with=F]
  
remove(max_10, max_30, max_60, real_used, limit_data, used_data)


### Classification
get_cl_training <- function(my_dataset_used, my_dataset_limit, my_dataset_max_10, my_dataset_max_30, my_dataset_max_60){
    my_dataset_used <- training_used
    my_dataset_limit <- training_limit
    my_dataset_max_10 <- training_max_10
    my_dataset_max_30 <- training_max_30
    my_dataset_max_60 <- training_max_60
    
    for(i in 1:(length(my_dataset_used)-89)){
        i <- 1
        temp_used <- my_dataset_used[,i:(89+i), with=F]
        temp_limit <- my_dataset_limit[,i:(89+i), with=F]
        temp_max_10 <- my_dataset_max_10[,(89+i), with=F]
        temp_max_30 <- my_dataset_max_30[,(89+i), with=F]
        temp_max_60 <- my_dataset_max_60[,(89+i), with=F]
      
        temp_used <- as.matrix(temp_used)
        temp_limit <- as.matrix(temp_limit)
        
        temp_dataset <- as.data.table(temp_used - temp_limit)
        temp_dataset[,]

        
        
    }
}











### Regression 























  
  
  
  
  





filtre_coefficient <- function(my_vector, my_case, my_limit=1){
    my_vector <- unlist(my_vector)
    if(length(my_vector[!is.na(my_vector)>2]) && !is.na(my_vector[length(my_vector)])){
        l_last_value <- my_vector[lenght(my_vector)]
        my_vector <- diff(my_vector, 1)[!is.na(my_vector)]
        ifelse(any(lapply(my_vector, function(x){x*my_case+l_last_value}))>=my_limit, return(1), return(0))
    }else{
        return(0)
    }
}

filtre_max_limit <- function(my_vector, my_limit=1){
    ifelse(any(my_vector>=my_limit), return(1), return(0))
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





