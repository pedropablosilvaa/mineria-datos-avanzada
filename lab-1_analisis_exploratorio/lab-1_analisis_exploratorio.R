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

#--------------------------------------------------------------------
#--------------------- Análisis exploratorio ------------------------
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


# basic histogram

hist_plot = function(df, variable) {
  p <- ggplot(df, aes_string(x=variable)) +
    theme()+
    geom_histogram(fill="#69b3a2",
                   color="#e9ecef",
                   alpha=0.9) +
    ggtitle(capture.output(cat("Histograma de variable:",variable))) + 
    theme_ipsum() +
    theme(
      plot.title = element_text(size=15), plot.background = element_rect(fill = 'white', colour = 'black')
    )
  path_file = capture.output(cat(wd_path,"/images", sep = ""))
  setwd(path_file)
  file_name = capture.output(cat(variable,"_plot.jpeg", sep = ""))
  ggsave(file_name, device = "jpeg")
  setwd(wd_path)
  return(p)
}

numcerical_var <- column_names_list[1:(length(column_names_list)-1)]

for (value in numcerical_var)
{
  hist_plot(data, value)
}


corr_data = cor(data)
head(round(corr_data,2))

# visualizing correlogram
# as circle
corrplot(corr_data, method="pie")
#ggsave('data/corrplot.jpeg', device = "jpeg")





#--------------------------------------------------------------------
#--------------------------- Clustering  ----------------------------
#--------------------------------------------------------------------


#mod1 = Mclust(data[,1:7]) #DEFAULT 
#summary(mod1)

#mod2 = Mclust(data[,1:7], G = 3)  #Numero de grupos = 3.
#summary(mod2, parameters = TRUE)

#mod6 = Mclust(data[,1:7], prior = priorControl(functionName="defaultPrior", shrinkage=0.1), modelNames ="EII")  
#"EII" = spherical, equal volume # Using prior #The function priorControl is used to specify a conjugate prior for EM within MCLUST. 
#summary(mod6,parameter = TRUE)
#plot(mod6, what = "classification")

class<-data$type

#------------------------- Experimento 1 --------------------------

#------------------------------ BIC 1 -----------------------------
BIC<-mclustBIC(data[,1:7], prior = priorControl(functionName="defaultPrior", shrinkage=0.1))
plot(BIC)  #se grafican los BIC por configuración de parámetros
summary(BIC)  # se presentan los mejores valores BIC

mod11=Mclust(data[,1:7],x=BIC) # en base al mejor valor BIC se realiza el mclust
summary(mod11)#se muestra resultado y tabla de clustering


plot(mod11, what = "classification")  #se grafica la configuración de agrupamientos.
legend("bottomright", legend = 1:3,
       col = mclust.options("classPlotColors"),
       pch = mclust.options("classPlotSymbols"),title = "Class labels:")

table(class, mod11$classification) #distribución de clases por cada grupo.

#------------------------------ ICL 1 -----------------------------

ICL<-mclustICL(data[,1:7])
plot(ICL)  #se grafican los ICL por configuración de parámetros
summary(ICL)  # se presentan los mejores valores ICL

mod12=Mclust(data[,1:7], G=4, prior = priorControl(functionName="defaultPrior", shrinkage=0.1), modelNames ="EEV") # en base al mejor valor ICL se realiza el mclust
summary(mod12)#se muestra resultado y tabla de clustering


plot(mod12, what = "classification")  #se grafica la configuración de agrupamientos.
legend("bottomright", legend = 1:4,
       col = mclust.options("classPlotColors"),
       pch = mclust.options("classPlotSymbols"),title = "Class labels:")

table(class, mod12$classification) #distribución de clases por cada grupo.









#------------------------- Experimento 2 --------------------------


#------------------------------ BIC 2 -----------------------------

BIC<-mclustBIC(data[,2:7], prior = priorControl(functionName="defaultPrior", shrinkage=0.1))
plot(BIC)  #se grafican los BIC por configuración de parámetros
summary(BIC)  # se presentan los mejores valores BIC

mod11=Mclust(data[,2:7], G=3, modelNames ="EEE") # en base al mejor valor BIC se realiza el mclust
summary(mod11)#se muestra resultado y tabla de clustering


plot(mod11, what = "classification")  #se grafica la configuración de agrupamientos.
legend("bottomright", legend = 1:3,
       col = mclust.options("classPlotColors"),
       pch = mclust.options("classPlotSymbols"),title = "Class labels:")

table(class, mod11$classification) #distribución de clases por cada grupo.

#------------------------------ ICL 2 -----------------------------

ICL<-mclustICL(data[,2:7])
plot(ICL)  #se grafican los ICL por configuración de parámetros
summary(ICL)  # se presentan los mejores valores ICL
data_2 = data.frame(data[,2:6])
mod12=Mclust(data_2, G=3, modelNames ="VEE") # en base al mejor valor ICL se realiza el mclust
summary(mod12)#se muestra resultado y tabla de clusteringw


plot(mod12, what = "classification")  #se grafica la configuración de agrupamientos.
legend("bottomright", legend = 1:3,
       col = mclust.options("classPlotColors"),
       pch = mclust.options("classPlotSymbols"),title = "Class labels:")

table(class, mod12$classification) #distribución de clases por cada grupo.







#------------------------- Experimento 3 --------------------------


#------------------------------ BIC 3 -----------------------------

BIC<-mclustBIC(data[,3:7], prior = priorControl(functionName="defaultPrior", shrinkage=0.1))
plot(BIC)  #se grafican los BIC por configuración de parámetros
summary(BIC)  # se presentan los mejores valores BIC

mod11=Mclust(data[,3:7],x=BIC) # en base al mejor valor BIC se realiza el mclust
summary(mod11)#se muestra resultado y tabla de clustering


plot(mod11, what = "classification")  #se grafica la configuración de agrupamientos.
legend("bottomright", legend = 1:3,
       col = mclust.options("classPlotColors"),
       pch = mclust.options("classPlotSymbols"),title = "Class labels:")

table(class, mod11$classification) #distribución de clases por cada grupo.

#------------------------------ ICL 3 -----------------------------

ICL<-mclustICL(data[,3:7])
plot(ICL)  #se grafican los ICL por configuración de parámetros
summary(ICL)  # se presentan los mejores valores ICL
data_3 = data.frame(data[,3:7])
mod12=Mclust(data_3, G=3, modelNames ="EVE") # en base al mejor valor ICL se realiza el mclust
summary(mod12)#se muestra resultado y tabla de clusteringw


plot(mod12, what = "classification")  #se grafica la configuración de agrupamientos.
legend("bottomright", legend = 1:3,
       col = mclust.options("classPlotColors"),
       pch = mclust.options("classPlotSymbols"),title = "Class labels:")

table(class, mod12$classification) #distribución de clases por cada grupo.


