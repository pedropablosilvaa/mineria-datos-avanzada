# Taller de Minería de datos avanzada
# Profesor: Max Chacón, ¿
# Alumno: Pedro Pablo Silva Antilef
# Magister en Ingniería Informática - Universidad de Santiago de Chile

# get path of current file
library("rstudioapi") 
# plot graphs
library("ggplot2")
library("tidyverse")
library("hrbrthemes")
# model-based clustering
library("mclust")
# correlation plot
library(corrplot)
library(randomForest)
library(MASS)

#--------------------------------------------------------------------
#-------------------------- Carga de datos   ------------------------
#--------------------------------------------------------------------


wd_path = dirname(getSourceEditorContext()$path)
setwd(wd_path)
setwd('..')
parent_path = getwd()
data_path = file.path(getwd(),
                      "data",
                      "seeds_dataset.txt")

#data_path
data = read.table(data_path,
                   header = FALSE,
                   sep = "")
#data
column_names_list = c("area",
                      "perimeter",
                      "compactness",
                      "l_kernel",
                      "w_kernel",
                      "asymmetry",
                      "l_groove",
                      "type")

colnames(data) = column_names_list 

anyNA(data)


#--------------------------------------------------------------------
#-------------------------- Random Forest   -------------------------
#--------------------------------------------------------------------



#------------------ Dataset Completo - Experimento 1  ---------------


data = read.table(data_path,
                  header = FALSE,
                  sep = "")
colnames(data) = column_names_list 

set.seed(6) #Semilla usada para la generación del proceso aleatorio.
data$type = as.factor(data$type)
data.rf = randomForest(type ~ .,
                       data=data,
                       ntree = 10000,
                       mtry=5,
                       importance=TRUE,
                       proximity=TRUE) # se genera un árbol en base a la Formula “type ~.” que determina como salida el atributo type y como entrada todas las variables del dataSet

print(data.rf) #muestra un resumen de los resultados de RF.

plot(data.rf)


round(importance(x = data.rf), 2)#Entrega la importancia de cada atributo sobre las instancias de la clase.
varImpPlot(data.rf)



data.mds = cmdscale(1 - data.rf$proximity, eig=TRUE) #escalamiento clásico multidimensional, usando valores propios
op = par(pty="s") #A character specifying the type of plot region to be used; "s" generates a square plotting region and "m" generates the maximal plotting region.
pairs(cbind(data[,1:7], data.mds$points), cex=0.6, gap=0,
      col=c("red", "green", "blue")[as.numeric(data$type)],
      main="Seeds Data: Predictors and MDS of Proximity Based on RandomForest") # se agregan los puntos del escalamiento multidimensional con los valores originales de las variables
par(op)

print(data.mds$GOF) #a numeric vector of length 2, equal to say (g.1,g.2), where g.i = (sum{j=1..k} λ[j]) / (sum{j=1..n} T.i(λ[j])), where λ[j] are the eigenvalues (sorted in decreasing order), T.1(v) = abs(v), and T.2(v) = max(v, 0)
MDSplot(data.rf, data$type)

plot(data.rf)



parcoord(data[,1:7],var.label = TRUE,col=c("red", "green", "blue")[as.numeric(data$type)])
legend("bottomright",legend = c("Kama", "Rosa", "Canadian"),cex=0.5, fill=2:4)


#------------------ Dataset sin "compactness" - Experimento 2  ---------------


data = read.table(data_path,
                  header = FALSE,
                  sep = "")
colnames(data) = column_names_list 

data_exp2 = subset(data, select = -c(compactness))



set.seed(1) 
data_exp2$type = as.factor(data_exp2$type)
data_exp2.rf = randomForest(type ~ .,
                       data=data_exp2,
                       ntree = 3000,
                       mtry=3,
                       importance=TRUE,
                       proximity=TRUE) # se genera un árbol en base a la Formula “type ~.” que determina como salida el atributo type y como entrada todas las variables del dataSet


print(data_exp2.rf) #muestra un resumen de los resultados de RF.

plot(data_exp2.rf)


round(importance(x = data_exp2.rf), 2)#Entrega la importancia de cada atributo sobre las instancias de la clase.
varImpPlot(data_exp2.rf)



data_exp2.mds = cmdscale(1 - data_exp2.rf$proximity, eig=TRUE) #escalamiento clásico multidimensional, usando valores propios
op = par(pty="s") #A character specifying the type of plot region to be used; "s" generates a square plotting region and "m" generates the maximal plotting region.
pairs(cbind(data_exp2[,1:6], data_exp2.mds$points), cex=0.6, gap=0,
      col=c("red", "green", "blue")[as.numeric(data_exp2$type)],
      main="Seeds Data: Predictors and MDS of Proximity Based on RandomForest") # se agregan los puntos del escalamiento multidimensional con los valores originales de las variables
par(op)

print(data_exp2.mds$GOF) #a numeric vector of length 2, equal to say (g.1,g.2), where g.i = (sum{j=1..k} λ[j]) / (sum{j=1..n} T.i(λ[j])), where λ[j] are the eigenvalues (sorted in decreasing order), T.1(v) = abs(v), and T.2(v) = max(v, 0)
MDSplot(data_exp2.rf, data_exp2$type)



parcoord(data_exp2[,1:6],var.label = TRUE,col=c("red", "green", "blue")[as.numeric(data_exp2$type)])
legend("bottomright",legend = c("Kama", "Rosa", "Canadian"),cex=0.5, fill=2:4)




#------------------ Dataset sólo con l_groove, area, perimeter - Experimento 3  ---------------



data = read.table(data_path,
                  header = FALSE,
                  sep = "")
colnames(data) = column_names_list 

data_exp3 = subset(data, select = c(l_groove, area, perimeter, type))

set.seed(1) 
data_exp3$type = as.factor(data_exp3$type)
data_exp3.rf = randomForest(type ~ .,
                            data=data_exp3,
                            ntree = 3000,
                            mtry=1,
                            importance=TRUE,
                            proximity=TRUE) # se genera un árbol en base a la Formula “type ~.” que determina como salida el atributo type y como entrada todas las variables del dataSet


print(data_exp3.rf) #muestra un resumen de los resultados de RF.

plot(data_exp3.rf)


#------------------ Dataset sólo con l_groove, area, asymmetry - Experimento 4  ---------------



data = read.table(data_path,
                  header = FALSE,
                  sep = "")
colnames(data) = column_names_list 

data_exp4 = subset(data, select = c(l_groove, area, asymmetry, type))

set.seed(1) 
data_exp4$type = as.factor(data_exp4$type)
data_exp4.rf = randomForest(type ~ .,
                            data=data_exp4,
                            ntree = 5000,
                            mtry=2,
                            importance=TRUE,
                            proximity=TRUE) # se genera un árbol en base a la Formula “type ~.” que determina como salida el atributo type y como entrada todas las variables del dataSet


print(data_exp4.rf) #muestra un resumen de los resultados de RF.

plot(data_exp4.rf)
legend("bottomright",legend = c("Kama", "Rosa", "Canadian"),cex=0.5, fill=2:4)

round(importance(x = data_exp4.rf), 2)#Entrega la importancia de cada atributo sobre las instancias de la clase.
varImpPlot(data_exp4.rf)



data_exp4.mds = cmdscale(1 - data_exp4.rf$proximity, eig=TRUE) #escalamiento clásico multidimensional, usando valores propios
op = par(pty="s") #A character specifying the type of plot region to be used; "s" generates a square plotting region and "m" generates the maximal plotting region.
pairs(cbind(data_exp4[,1:3], data_exp4.mds$points), cex=0.6, gap=0,
      col=c("red", "green", "blue")[as.numeric(data_exp4$type)],
      main="Seeds Data: Predictors and MDS of Proximity Based on RandomForest") # se agregan los puntos del escalamiento multidimensional con los valores originales de las variables
par(op)

print(data_exp4.mds$GOF) #a numeric vector of length 2, equal to say (g.1,g.2), where g.i = (sum{j=1..k} λ[j]) / (sum{j=1..n} T.i(λ[j])), where λ[j] are the eigenvalues (sorted in decreasing order), T.1(v) = abs(v), and T.2(v) = max(v, 0)
MDSplot(data_exp4.rf, data_exp4$type)
legend("bottomright",legend = c("Kama", "Canadian", "Rosa"),cex=0.5, fill=2:4)



parcoord(data_exp4[,1:3],var.label = TRUE,col=c("red", "green", "blue")[as.numeric(data_exp4$type)])
legend("bottomright",legend = c("Kama", "Rosa", "Canadian"),cex=0.5, fill=2:4)




#------------------ Dataset sólo con l_groove, area - Experimento 5  ---------------



data = read.table(data_path,
                  header = FALSE,
                  sep = "")
colnames(data) = column_names_list 

data_exp5 = subset(data, select = c(l_groove, area, type))

set.seed(1) 
data_exp5$type = as.factor(data_exp5$type)
data_exp5.rf = randomForest(type ~ .,
                            data=data_exp5,
                            ntree = 1000,
                            mtry=1,
                            importance=TRUE,
                            proximity=TRUE) # se genera un árbol en base a la Formula “type ~.” que determina como salida el atributo type y como entrada todas las variables del dataSet


print(data_exp5.rf) #muestra un resumen de los resultados de RF.

plot(data_exp5.rf)




#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#------------------------------------- Experimentos finales   ---------------------------------
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------


#---------------------------------- Experimento 1  -------------------------------------------



data = read.table(data_path,
                  header = FALSE,
                  sep = "")
colnames(data) = column_names_list 

data_exp4 = subset(data, select = c(l_groove, area, asymmetry, type))

set.seed(15) 
data_exp4$type = as.factor(data_exp4$type)
data_exp4.rf = randomForest(type ~ .,
                            data=data_exp4,
                            ntree = 400,
                            mtry=2,
                            importance=TRUE,
                            proximity=TRUE) # se genera un árbol en base a la Formula “type ~.” que determina como salida el atributo type y como entrada todas las variables del dataSet


print(data_exp4.rf) #muestra un resumen de los resultados de RF.

plot(data_exp4.rf)
legend("bottomright",legend = c("Kama", "Rosa", "Canadian"),cex=0.5, fill=2:4)

round(importance(x = data_exp4.rf), 2)#Entrega la importancia de cada atributo sobre las instancias de la clase.
varImpPlot(data_exp4.rf)



data_exp4.mds = cmdscale(1 - data_exp4.rf$proximity, eig=TRUE) #escalamiento clásico multidimensional, usando valores propios
op = par(pty="s") #A character specifying the type of plot region to be used; "s" generates a square plotting region and "m" generates the maximal plotting region.
pairs(cbind(data_exp4[,1:3], data_exp4.mds$points), cex=0.6, gap=0,
      col=c("red", "green", "blue")[as.numeric(data_exp4$type)],
      main="Seeds Data: Predictors and MDS of Proximity Based on RandomForest") # se agregan los puntos del escalamiento multidimensional con los valores originales de las variables
par(op)

print(data_exp4.mds$GOF) #a numeric vector of length 2, equal to say (g.1,g.2), where g.i = (sum{j=1..k} λ[j]) / (sum{j=1..n} T.i(λ[j])), where λ[j] are the eigenvalues (sorted in decreasing order), T.1(v) = abs(v), and T.2(v) = max(v, 0)
MDSplot(data_exp4.rf, data_exp4$type)
legend("bottomright",legend = c("Kama", "Canadian", "Rosa"),cex=0.5, fill=2:4)



parcoord(data_exp4[,1:3],var.label = TRUE,col=c("red", "green", "blue")[as.numeric(data_exp4$type)])
legend("bottomright",legend = c("Kama", "Rosa", "Canadian"),cex=0.5, fill=2:4)



#---------------------------------- Experimento 2  -------------------------------------------


data = read.table(data_path,
                  header = FALSE,
                  sep = "")
colnames(data) = column_names_list 

set.seed(28) #Semilla usada para la generación del proceso aleatorio.
data$type = as.factor(data$type)
data.rf = randomForest(type ~ .,
                       data=data,
                       ntree = 600,
                       mtry=5,
                       importance=TRUE,
                       proximity=TRUE) # se genera un árbol en base a la Formula “type ~.” que determina como salida el atributo type y como entrada todas las variables del dataSet

print(data.rf) #muestra un resumen de los resultados de RF.

plot(data.rf)


round(importance(x = data.rf), 2)#Entrega la importancia de cada atributo sobre las instancias de la clase.
varImpPlot(data.rf)



data.mds = cmdscale(1 - data.rf$proximity, eig=TRUE) #escalamiento clásico multidimensional, usando valores propios
op = par(pty="s") #A character specifying the type of plot region to be used; "s" generates a square plotting region and "m" generates the maximal plotting region.
pairs(cbind(data[,1:7], data.mds$points), cex=0.6, gap=0,
      col=c("red", "green", "blue")[as.numeric(data$type)],
      main="Seeds Data: Predictors and MDS of Proximity Based on RandomForest") # se agregan los puntos del escalamiento multidimensional con los valores originales de las variables
par(op)

print(data.mds$GOF) #a numeric vector of length 2, equal to say (g.1,g.2), where g.i = (sum{j=1..k} λ[j]) / (sum{j=1..n} T.i(λ[j])), where λ[j] are the eigenvalues (sorted in decreasing order), T.1(v) = abs(v), and T.2(v) = max(v, 0)
MDSplot(data.rf, data$type)

plot(data.rf)



parcoord(data[,1:7],var.label = TRUE,col=c("red", "green", "blue")[as.numeric(data$type)])
legend("bottomright",legend = c("Kama", "Rosa", "Canadian"),cex=0.5, fill=2:4)

