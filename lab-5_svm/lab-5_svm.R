# Taller de Minería de datos avanzada
# Profesor: Max Chacón, Felipe Bello Robles
# Alumno: Pedro Pablo Silva Antilef
# Magister en Ingniería Informática - Universidad de Santiago de Chile

# get path of current file
library("rstudioapi") 
library(e1071)
library(RWeka)
library(caTools)

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
#-------------  Selección de características ------------------------
#--------------------------------------------------------------------


data$type <- as.factor(data$type)
ranking <- InfoGainAttributeEval(formula=type ~  ., data = data)
print(ranking)



#--------------------------------------------------------------------
#-------------------- Support Vector Machine   ----------------------
#--------------------------------------------------------------------


#------------------------- Dataset completo   -----------------------

x <- subset(data, select = -type)
y <- type

# Kernel lineal considerando dataset completo
set.seed(40)
obj <- tune(svm, type~., data = data, kernel = "linear",ranges = list(cost = 2^(-3:6)), tunecontrol = tune.control(sampling = "cross", cross = 2 ))
summary(obj)
# Se plotea los errores con los costos
plot(obj)
summary(obj$best.model)

# se genera la prediccion con el mejor modelo
pred <- predict(obj$best.model, x)
table(pred, type)

# kernel radial
obj <- tune(svm, type~., data = data, kernel = "radial", ranges = list(gamma = 2^(-4:8), cost = 2^(-3:6), tunecontrol = tune.control(sampling = "cross", cross = 2 )))
summary(obj)
plot(obj$gamma)

pred <- predict(obj$best.model, x)
table(pred, type)



#------------------------- Dataset editado   -----------------------

#seleccion de caracteristicas
data_edited = data[,c('area','l_kernel', 'l_groove','type')]

# kernel lineal
obj <- tune(svm, type~., data = data_edited, kernel = "linear",ranges = list(cost = 2^(-3:6)), tunecontrol = tune.control(sampling = "cross", cross = 2 ))
summary(obj)
plot(obj)

pred <- predict(obj$best.model, x)
table(pred, type)


# kernel radial
obj <- tune(svm, type~., data = data_edited, kernel = "radial", ranges = list(gamma = 2^(-4:8), cost = 2^(-3:6), tunecontrol = tune.control(sampling = "cross", cross = 2 )))
summary(obj)
pred <- predict(obj$best.model, x)
table(pred, type)


#set.seed(40)
#spliting dataset
#split = sample.split(data$type, SplitRatio = 0.65)
#training_set = subset(data, split == TRUE)
#test_set = subset(data, split == FALSE)

#svm on trianing set
#classifier = svm(formula = type ~ .,
#                 data = training_set,
#                 type = 'C-classification',
#                 kernel = 'linear')

#y_pred = predict(classifier, newdata = test_set[-8])
#cm = table(test_set[, 8], y_pred)
#cm




#s





attach(data)

x <- subset(data, select = -type)
y <- type
model <- svm(x, y)
# test with train data
pred <- predict(model, x)
table(pred, y)

# compute decision values and probabilities:
pred <- predict(model, x, decision.values = TRUE)
attr(pred, "decision.values")[1:4,]

# visualize (classes by color, SV by crosses):
plot(cmdscale(dist(data[,-8])), col = as.integer(data[,8]), pch = c("o","+")[1:150 %in% model$index + 1])

