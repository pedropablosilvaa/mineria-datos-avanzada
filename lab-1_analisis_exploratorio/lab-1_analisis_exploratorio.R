# Taller de Minería de datos avanzada
# Profesor: Max Chacón, Felipe-Andrés Bello
# Alumno: Pedro Pablo Silva Antilef
# Magister en Ingniería Informática - Universidad de Santiago de Chile

getwd()
wd_path = getwd()
setwd('..')
parent_path = getwd()
data_path = file.path(getwd(),
                      "data",
                      "seeds_dataset.txt")
setwd(wd_path)
#data_path
data = read.table(data_path,
                   header = FALSE,
                   sep = "")
#data
colnames(data) <- c("area",
                    "perimeter",
                    "compactness"
                    ,"l_kernel",
                    "l_width",
                    "asymmetry",
                    "l_groove",
                    "type")

