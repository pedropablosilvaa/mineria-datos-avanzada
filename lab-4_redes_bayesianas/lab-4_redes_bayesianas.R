# Taller de Minería de datos avanzada
# Profesor: Max Chacón, Felipe 
# Alumno: Pedro Pablo Silva Antilef
# Magister en Ingniería Informática - Universidad de Santiago de Chile

# get path of current file
library("rstudioapi") 
# plot graphs
library("ggplot2")
library("tidyverse")
library("hrbrthemes")
#Max entropy
library("maxent")
#clean corpus
library("SnowballC")
# word cloud
library("wordcloud")
library("stringi")
library("caret")
library("bnlearn")
library('dplyr')
library("data.table")
library("stringr")
library("bnlearn")
#--------------------------------------------------------------------
#-------------------------- Carga de datos   ------------------------
#--------------------------------------------------------------------








#---------------------------------------------------------------------------------------------
#-------------------------- Experimento 1 - Con todos las relaciones   ----------------------
#-------------------------------------------------------------------------------------------


data(hailfinder)

#----------- Hill-Climbing---------

#relaciones = se aplica algoritmo al df
res_1_hc <- hc(hailfinder) #Algoritmo Hill-Climbing
#plot(res_1_hc)
graphviz.plot(res_1_hc, layout = "fdp")
sc_1_hc<-score(res_1_hc,hailfinder) # BIC por default
print(sc_1_hc)

#----------- Max-Min Hill-Climbing (mmhc) ---------

res_1_mmhc <- mmhc(hailfinder) #Algoritmo Hill-Climbing
#plot(res_1_mmhc)
graphviz.plot(res_1_mmhc, layout = "fdp")
sc_1_mmhc<-score(res_1_mmhc,hailfinder) # BIC por default
print(sc_1_mmhc)


#----------- Max-Min Parents and Children (mmpc)---------

res_1_mmpc <- mmpc(hailfinder) #Algoritmo Hill-Climbing
#plot(res_1_mmhc)
graphviz.plot(res_1_mmpc, layout = "fdp")
sc_1_mmpc<-score(res_1_mmpc,hailfinder) # BIC por default
#print(sc_1_mmpc)



#---------------------------------------------------------------------------------------------
#-------------------------- Experimento 2 - Con blacklist------------   ----------------------
#---------------------------------------------------------------------------------------------

#origin_list = list("IRCloudCover","VISCloudCov")
#destination_list = list("CombClouds","CombClouds")
#bl<-data.frame(unlist(origin_list),
#               unlist(destination_list)
#               ) #Lista negra de relaciones (Origen, Destino)

hailfinder_changed <- subset(hailfinder, select = -c(IRCloudCover, VISCloudCov,N07muVerMo,
                                                     SubjVertMo,QGVertMotion, SatContMoist,
                                                     RaoContMoist, VISCloudCov, IRCloudCover,
                                                     LowLLapse, MeanRH, MidLLapse ))
# Se aplica el algoritmo

res_2_hc <- hc(hailfinder_changed)
graphviz.plot(res_2_hc, layout = "fdp")
print(res_2_hc)
sc_2_hc<-score(res_2_hc,hailfinder_changed) # BIC por default
print(sc_2_hc)


res_2_mmhc <- mmhc(hailfinder_changed)
graphviz.plot(res_2_mmhc, layout = "fdp")
print(res_2_mmhc)
sc_2_mmhc<-score(res_2_mmhc,hailfinder_changed) # BIC por default
print(sc_2_mmhc)


res_2_mmpc <- mmpc(hailfinder_changed)
graphviz.plot(res_2_mmpc, layout = "fdp")
print(res_2_mmpc)
sc_2_mmpc<-score(res_2_mmpc,hailfinder_changed) # BIC por default
print(sc_2_mmpc)



fittedbn_2 <- bn.fit(res_2_mmhc, data = hailfinder_changed) # Se obtiene la tabla de probabilidades condicionales mediante EM. (Máxima Expectación, propagación de la evidencia)
print(fittedbn_2$R5Fcst) #se obtiene la información respecto del nodo Proteins


cpquery(fittedbn_2, event = (R5Fcst=="SVR"), evidence = (CombMoisture=="VeryWet")) 

cpquery(fittedbn_2, event = (R5Fcst=="SVR"), evidence = (CombVerMo=="StrongUp")) 

cpquery(fittedbn_2, event = (R5Fcst=="SVR"), evidence = (CombClouds=="Cloudy")) 

cpquery(fittedbn_2, event = (R5Fcst=="SVR"), evidence = (LLIW=="Strong")) 

cpquery(fittedbn_2, event = (R5Fcst=="SVR"), evidence = (CldShadeConv=="Marked")) 

cpquery(fittedbn_2, event = (R5Fcst=="SVR"), evidence = (WndHodograph=="DCVZFavor")) 
