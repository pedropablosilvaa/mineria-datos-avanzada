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

set.seed(1) #Semilla usada para la generación del proceso aleatorio.
data$type = as.factor(data$type)
data.rf = randomForest(type ~ .,
                       data=data,
                       ntree = 1500,
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
legend("bottomright",legend = c("setosa", "versicolor", "virginica"),cex=0.5, fill=2:4)


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
                       ntree = 1000,
                       mtry=6,
                       importance=TRUE,
                       proximity=TRUE) # se genera un árbol en base a la Formula “type ~.” que determina como salida el atributo type y como entrada todas las variables del dataSet


print(data_exp2.rf) #muestra un resumen de los resultados de RF.

plot(data_exp2.rf)


#------------------ Dataset sin "compactness" - Experimento 3  ---------------
