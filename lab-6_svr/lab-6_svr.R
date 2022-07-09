# Taller de Minería de datos avanzada
# Profesor: Max Chacón, Felipe Bello Robles
# Alumno: Pedro Pablo Silva Antilef
# Magister en Ingniería Informática - Universidad de Santiago de Chile

# get path of current file
library("rstudioapi") 
library(e1071)
library(doParallel)
library(stats)
library(babynames)
library(dplyr)
#library(RWeka)
#library(caTools)


#----------------------------------------------------------------------
#----------------------------- Funciones----- -------------------------
#----------------------------------------------------------------------


#funcion que genera los retardos multiples
retardos_multi <- function(
    signalData,
    lags
    )
{
  
  signal.uni <- signalData
  max.lag <- max(unlist(lags)) + 1
  indices <- 1:nrow(signal.uni)
  lag.mat <- embed(indices, max.lag)
  
  col.names <- list("PAMn","VFSCn")
  columns <- NULL
  lagged.columns.names <- c()
  for(colname in col.names){
    
    lag.order <- lags[[colname]]
    columns[[colname]] <- signal.uni[lag.mat[, 1], colname]
    if(!is.null(lag.order) && lag.order > 0){
      for(i in 1:lag.order){
        new.colname <- paste(colname, paste0("lag", i), sep = ".")
        lagged.columns.names <- c(lagged.columns.names, new.colname)
        columns[[new.colname]] <- signal.uni[lag.mat[, i+1], colname]
      }
    
    }
    
  }
  folded.signal <- data.frame(columns)
  
  sorting <- order(lag.mat[, 1])
  folded.signal <- folded.signal[sorting, ]
  list(folded.signal = folded.signal, lagged.columns.names = lagged.columns.names)
}


#Funcion para dividir dataset en partes iguales conservando el orden

split_df_in_half <- function(df) {
  test <- split(df, rep(1:2,each=length(df)/2))
  return(test)
}


#Funcion para obtener entrenar y testear modelos considerando una grilla de hiperparametros
generate_model <- function(train_df, test_df, parameters){
  start_time <- Sys.time()
  print(start_time)
  salida <- (c( foreach(i = 1:nrow(parameters),  combine = rbind, .inorder = FALSE) %dopar% {
    c <- parameters[i, ]$cost
    n <- parameters[i, ]$nu
    g <- parameters[i, ]$gamma
    l <- parameters[i, ]$lagsList
    lag<-list(PAMn = 1,VFSCn = 0)
    #aca deberia ir training data
    signal.train <- retardos_multi(train_df, lag)
    #se guardan los retardos
    retDatos=signal.train$folded.signal
    x=subset(retDatos, select = -VFSCn)
    y=retDatos$VFSCn
    #se genera el modelo con los parametros definidos anteriormente
    modelo <- e1071::svm(x, y, type = "nu-regression", kernel = "radial", cost = c, nu = n, gamma=g)
    
    #aca generar los retardos de igual manera como se hizo anteriormente
    signal.test <- retardos_multi(test_df, lag)
    retDatos.test = signal.test$folded.signal
    
    x_test = subset(retDatos.test, select = -VFSCn)
    y_test = retDatos.test$VFSCn
    
    #se genera la prediccion
    pred <- predict(modelo, x_test)
    corr_pred <- cor(pred, y_test, method = "pearson")
    
    
    
    c(l, c, n, g, corr_pred)
    
  }))
  
  
  output <- matrix(unlist(salida), ncol = 5, byrow = TRUE)
  mejoresModelos<-output[order(output[,5], decreasing = TRUE),]
  print(mejoresModelos)
  
  df_resultados <- as.data.frame(mejoresModelos)
  colnames(df_resultados) <- c("lag", "cost", "nu", "gamma", "corr_pred")
  df_resultados <- df_resultados[order(df_resultados$corr_pred, decreasing = TRUE), ]
  
  #write.csv(df_resultados, name_df, row.names = FALSE)
  end_time <- Sys.time()
  print(end_time - start_time)
  return(df_resultados)
  
  
}


# Funcion que simula el escalon de presion y plotea el resultado

autoregulacion <- function(best_results_df){
  start_time <- Sys.time()
  for (i in 1:200){
    
    PAMn<-(data_1$PAM-min(data_1$PAM))/(max(data_1$PAM)-min(data_1$PAM))
    VFSCn<-(data_1$VFSC-min(data_1$VFSC))/(max(data_1$VFSC)-min(data_1$VFSC))
    data <- data.frame(PAMn,VFSCn)
    lag<-list(PAMn = best_results_df[i,1],VFSCn = 0)
    signal.train <- retardos_multi(data, lag)
    retDatos=signal.train$folded.signal
    
    x=subset(retDatos, select = -VFSCn)
    y=retDatos$VFSCn
    mejorModelo <- svm(x, y, kernel = "radial",type = "nu-regression", cost = best_results_df[i,2], nu = best_results_df[i,3], gamma=best_results_df[i,4])
    
    PAMn=inverseStep
    VFSCn=inverseStep 
    data <- data.frame(PAMn,VFSCn)
    lag<-list(PAMn = best_results_df[i,1],VFSCn = 0)
    signal.train <- retardos_multi(data, lag)
    retDatos=signal.train$folded.signal
    x=subset(retDatos, select = -VFSCn)
    y=retDatos$VFSCn
    
    stepTime=seq(Ts,(length(retDatos$PAMn))*Ts,Ts)
    stepResponse <- predict(mejorModelo, x ) 
    plot(stepTime,retDatos$PAMn,type="l", col="red")
    lines(stepTime,stepResponse, col = "blue")
    legend("topright", c("Escalon de presión", "respuesta al escalón"), title = "autorregulación", pch = 1, col=c("red","blue"),lty=c(1,1))
    print(paste("corr=",best_results_df[i,5]))
    readline(prompt="Press [enter] to continue")
  }

  end_time <- Sys.time()
  end_time - start_time
  
  
}

# constantes

Ts=0.2
Tiempo=seq(Ts,(length(VFSC))*Ts,Ts) 


#se define el numero de nucleos
registerDoParallel(cores = 8)
#se definen los parametros de la grilla
cost <- 2^seq(-4, 12, 1)
nu <- seq(0.1, 0.9, 0.1)
gamma<-2^seq(-6, 12, 1)
lagsList<-seq(1,5,1)


#--------------------------------------------------------------------
#-------------------------- Carga de datos   ------------------------
#--------------------------------------------------------------------


wd_path = dirname(getSourceEditorContext()$path)
setwd(wd_path)
setwd('..')
parent_path = getwd()
data_path_1 = file.path(getwd(),
                      "data",
                      "G1_001.csv")

data_path_2 = file.path(getwd(),
                        "data",
                        "G1_002.csv")

#data paciente 1
data_1 = read.csv(data_path_1,
                  header = TRUE
                )
#data paciente 2
data_2 = read.csv(data_path_2,
                  header = TRUE
)

data_2 = head(data_2, - 6)              # Apply head function
data_2                  

anyNA(data_1)
anyNA(data_2)

wd_path = dirname(getSourceEditorContext()$path)
setwd(wd_path)


#--------------------------------------------------------------------
#------------------------  Paciente A -------------------------------
#--------------------------------------------------------------------




#--------------------------------------------
#-----------  Paciente 1  -------------------
#--------------------------------------------

set.seed(100)
attach(data_1)


#formula del modelo
formula=VFSC ~ PAM


#se cargan los parametros de la grilla
parms <- expand.grid(lagsList=lagsList, cost = cost, nu = nu, gamma=gamma)


#Normalizacion de los datos
PAMn<-(data_1$PAM-min(data_1$PAM))/(max(data_1$PAM)-min(data_1$PAM))
VFSCn<-(data_1$VFSC-min(data_1$VFSC))/(max(data_1$VFSC)-min(data_1$VFSC))
#dataframe a normalizado
data_1_n <- data.frame(PAMn,VFSCn)
#se genera los pasos de tiempo
Ts=0.2
Tiempo=seq(Ts,(length(VFSC))*Ts,Ts) 

# Se divide el dataset en dos
#ind <- sample(2, nrow(data_1_n), replace = TRUE,prob = c(0.5, 0.5))

length= as.numeric(length(data_1_n[, 1])/2)
data_1_A =data_1_n[1:length, ]
data_1_B =tail(data_1_n, length)
#data_1_n




#gamma -> 
#costo -> ok
#nu -> 
model_train_data_1A = generate_model(data_1_A, data_1_B, parms)
#8:37PM start
write.csv(model_train_data_1A, './output_pedro/train_1A.csv', row.names = FALSE)




model_train_data_1B = generate_model(data_1_B, data_1_A, parms)
write.csv(model_train_data_1B, './output/train_1B_refine.csv', row.names = FALSE)




inverseStep=matrix(1,300/Ts,1)
inverseStep[(150/Ts):(300/Ts),1]=0


#train_data_1A <- train_data_1A[order(train_data_1A$corr_pred, decreasing = TRUE), ]


autoregulacion(model_train_data_1A)

autoregulacion(model_train_data_1B)


#paciente1_trainA


#--------------------------------------------
#-----------  Paciente 2  -------------------
#--------------------------------------------

set.seed(100)
attach(data_2)


#formula del modelo
formula=VFSC ~ PAM

#se define el numero de nucleos
registerDoParallel(cores = 8)


#Normalizacion de los datos
PAMn<-(data_2$PAM-min(data_2$PAM))/(max(data_2$PAM)-min(data_2$PAM))
VFSCn<-(data_2$VFSC-min(data_2$VFSC))/(max(data_2$VFSC)-min(data_2$VFSC))



#dataframe a normalizado
data_2_n <- data.frame(PAMn,VFSCn)
#se genera los pasos de tiempo
Ts=0.2
Tiempo=seq(Ts,(length(VFSC))*Ts,Ts) 

# Se divide el dataset en dos
#ind <- sample(2, nrow(data_1_n), replace = TRUE,prob = c(0.5, 0.5))

length= as.numeric(length(data_2_n[, 1])/2)
data_2_A =data_2_n[1:length, ]
data_2_B =tail(data_2_n, length)
#data_1_n

#se definen los parametros de la grilla
#cost <- 2^seq(-5, -3, 1)
#nu <- seq(0.1, 0.3, 0.1)
#gamma<-2^seq(-8, -4, 1)
#lagsList<-seq(1,5,1)
#se cargan los parametros de la grilla
parms <- expand.grid(lagsList=lagsList, cost = cost, nu = nu, gamma=gamma)


model_train_data_2A = generate_model(data_2_A, data_2_B, parms)
Sys.time()
write.csv(model_train_data_2A, './output/train_2A_refine.csv', row.names = FALSE)



model_train_data_2B = generate_model(data_2_B, data_2_A, parms)
Sys.time()
write.csv(model_train_data_2B, './output/train_2B_refine.csv', row.names = FALSE)






#-----------------------------------------------------------------------
#---------------------- extraccion de informacion  ----------------------
#-----------------------------------------------------------------------
# Paciente 1


best_model_1 <- read.csv("./output/train_1B_refine.csv")
#best_model_1 <- best_model_1[order(best_model_1[,5], decreasing = TRUE),]
#test_best_model_1 <- best_model_1[best_model_1$corr_pred > 0.35, ]

PAMn <- (data_1$PAM-min(data_1$PAM))/(max(data_1$PAM)-min(data_1$PAM))
VFSCn <- (data_1$VFSC-min(data_1$VFSC))/(max(data_1$VFSC)-min(data_1$VFSC))
data_1_n <- data.frame(PAMn,VFSCn)


lag<-list(PAMn = best_model_1[5,1],VFSCn = 0)
signal.train <- retardos_multi(data_1_n, lag)
retDatos=signal.train$folded.signal

x=subset(retDatos, select = -VFSCn)
y=retDatos$VFSCn

best_svm_1 <- svm(x, y, kernel = "radial",type = "nu-regression", cost = best_model_1[5,2],
                        nu = best_model_1[5,3], gamma=best_model_1[5,4])

print(best_svm_1)

inverseStep=matrix(1,(300/Ts),1)
inverseStep[(150/Ts):(300/Ts),1]=0

PAMn=inverseStep
VFSCn=inverseStep 

data <- data.frame(PAMn, VFSCn)
lag <- list(PAMn = best_model_1[5,1],VFSCn = 0)

signal.train <- retardos_multi(data, lag)

retDatos = signal.train$folded.signal
x = subset(retDatos, select = -VFSCn)


stepTime=seq(Ts,(length(retDatos$PAMn))*Ts,Ts)

stepResponse <- predict(best_svm_1, x)

plot(stepTime,retDatos$PAMn,type="l", col="red")
lines(stepTime,stepResponse, col = "blue")
legend("topright", 
       c("Escalon de presión", 
        "respuesta al escalon"), 
       title = "autorregulacion", 
       pch = 1, 
       col=c("red","blue"),
       lty=c(1,1),
       inset = 0.01,
       cex = 0.7)

#----------------

formula <- data_1$VFSC ~ data_1$PAM

best_model_1_unnormalized <- svm(formula, data_1, kernel = "radial",type = "nu-regression", cost = best_model_1[5,2],
                                      nu = best_model_1[5,3], gamma=best_model_1[5,4])

VFSC_tunedModel <- predict(best_model_1_unnormalized, data_1$PAM)

plot(Tiempo, data_1$VFSC, type="l", col = "blue")
lines(Tiempo, VFSC_tunedModel, col = "red")
legend("topright", c("VFSC","VFSC_estimated"), col = c("blue", "red"),
       title = "VFSC vs VFSC estimado",  lty=1:2, cex=0.5)



#-----------------------------------



#plot(Tiempo,data_1$PAM, type="l", col = "black", ylim=range( c(data_1$PAM, data_1$VFSC) ))
#lines(Tiempo, data_1$VFSC,  col = "red")
#legend("right",
#      col=c("black","red"),
#       legend = c("VFSC","PAM"),
#       title = "Datos Paciente 1",
#       pch = 1,
#       
#       lty=c(2,1),
#       inset = 0.01
#       )




#-----------------------------------------------------------------------
#---------------------- extraccion de informacion  ----------------------
#-----------------------------------------------------------------------
# Paciente 2



best_model_2 <- read.csv("./output/train_2A_refine.csv")

PAMn <- (data_2$PAM-min(data_2$PAM))/(max(data_2$PAM)-min(data_2$PAM))
VFSCn <- (data_2$VFSC-min(data_2$VFSC))/(max(data_2$VFSC)-min(data_2$VFSC))

data_2_n <- data.frame(PAMn,VFSCn)

lag<-list(PAMn = best_model_2[5,1],VFSCn = 0)
signal.train <- retardos_multi(data_2_n, lag)
retDatos=signal.train$folded.signal

x=subset(retDatos, select = -VFSCn)
y=retDatos$VFSCn

best_svm_2 <- svm(x, y, kernel = "radial",type = "nu-regression", cost = best_model_2[5,2],
                  nu = best_model_2[5,3], gamma=best_model_2[5,4])

print(best_svm_2)

inverseStep=matrix(1,(300/Ts),1)
inverseStep[(150/Ts):(300/Ts),1]=0

PAMn=inverseStep
VFSCn=inverseStep 

data <- data.frame(PAMn, VFSCn)
lag <- list(PAMn = best_model_2[5,1],VFSCn = 0)

signal.train <- retardos_multi(data, lag)

retDatos = signal.train$folded.signal
x = subset(retDatos, select = -VFSCn)


stepTime=seq(Ts,(length(retDatos$PAMn))*Ts,Ts)

stepResponse <- predict(best_svm_2, x)

plot(stepTime,retDatos$PAMn,type="l", col="red")
lines(stepTime,stepResponse, col = "blue")
legend("topright", 
       c("Escalon de presión", 
         "Respuesta al escalón"), 
       title = "Autorregulación", 
       pch = 1, 
       col=c("red","blue"),
       lty=c(1,1),
       inset = 0.01,
       cex = 0.7)

#----------------

formula <- data_2$VFSC ~ data_2$PAM

best_model_2_unnormalized <- svm(formula, data_2, kernel = "radial",type = "nu-regression", cost = best_model_2[5,2],
                                 nu = best_model_2[5,3], gamma=best_model_2[5,4])

VFSC_tunedModel <- predict(best_model_2_unnormalized, data_2$PAM)

plot(Tiempo, data_2$VFSC, type="l", col = "blue")
lines(Tiempo, VFSC_tunedModel, col = "red")
legend("topright", c("VFSC","VFSC_estimated"), col = c("blue", "red"),
       title = "VFSC vs VFSC estimado",  lty=1:2, cex=0.5)




#plot(Tiempo,data_2$PAM, type="l", col = "black", ylim=range( c(data_2$PAM, data_2$VFSC) ))
#lines(Tiempo, data_2$VFSC,  col = "red")
#legend("right",
#       col=c("black","red"),
#       legend = c("VFSC","PAM"),
#       title = "Datos Paciente 2",
#       pch = 1,
#       
#       lty=c(2,1),
#       inset = 0.01
#)


