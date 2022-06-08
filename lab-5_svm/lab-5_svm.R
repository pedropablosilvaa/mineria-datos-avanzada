# Taller de Minería de datos avanzada
# Profesor: Max Chacón, Felipe Bello Robles
# Alumno: Pedro Pablo Silva Antilef
# Magister en Ingniería Informática - Universidad de Santiago de Chile

library(e1071)
library(RWeka)


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


attach(data)

formula=type ~ .
model <- svm(formula, data = data)
print(model)
summary(model)


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
plot(cmdscale(dist(data[,-7])), col = as.integer(data[,7]), pch = c("o","+")[1:150 %in% model$index + 1])


obj <- tune(svm, type~., data = data, kernel = "linear",ranges = list(cost = 2^(-1:4)), tunecontrol = tune.control(sampling = "cross", cross = 2 ))
summary(obj)

plot(obj)
summary(obj$best.model)


pred <- predict(obj$best.model, x)
table(pred, type)


obj <- tune(svm, type~., data = data, kernel = "radial", ranges = list(gamma = 2^(-2:4), cost = 2^(-1:4), tunecontrol = tune.control(sampling = "cross", cross = 2 )))
summary(obj)


pred <- predict(obj$best.model, x)
table(pred, type)

obj <- tune(svm, type~., data = data, kernel = "radial", ranges = list(gamma = 2^(-7:12), cost = 2^(-7:14) , tunecontrol = tune.control(sampling = "cross", cross = 2)))
pred <- predict(obj$best.model, x)
table(pred, type)
