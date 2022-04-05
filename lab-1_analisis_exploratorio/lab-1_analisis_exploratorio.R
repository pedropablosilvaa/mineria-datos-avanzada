# Taller de Minería de datos avanzada
# Profesor: Max Chacón, Felipe-Andrés Bello
# Alumno: Pedro Pablo Silva Antilef
# Magister en Ingniería Informática - Universidad de Santiago de Chile

# to get path of current file
library("rstudioapi") 
# to plot graphs
library("ggplot2")
library("tidyverse")
library("hrbrthemes")



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
                      "l_width",
                      "asymmetry",
                      "l_groove",
                      "type")

colnames(data) = column_names_list 

anyNA(data)


# basic histogram

hist_plot = function(df, variable) {
  p <- ggplot(df, aes_string(x=variable)) +
    geom_histogram(aes(y =..density..),
                       fill="#69b3a2",
                       color="#e9ecef",
                       alpha=0.9) +
    geom_density(col = "red")
    ggtitle(capture.output(cat("Histograma de variable:",variable))) + 
    theme_ipsum() +
    theme(
      plot.title = element_text(size=15)
    )
  path_file = capture.output(cat(wd_path,"/images", sep = ""))
  setwd(path_file)
  file_name = capture.output(cat(variable,"_plot.png", sep = ""))
  ggsave(file_name, device = "png")
  setwd(wd_path)
  return(p)
}

numcerical_var <- column_names_list[1:(length(column_names_list)-1)]

for (value in numcerical_var)
{
  hist_plot(data, value)
}


