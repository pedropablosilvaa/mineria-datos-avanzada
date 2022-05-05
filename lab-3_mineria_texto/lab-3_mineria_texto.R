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


#--------------------------------------------------------------------
#-------------------------- Carga de datos   ------------------------
#--------------------------------------------------------------------


wd_path = dirname(getSourceEditorContext()$path)
setwd(wd_path)
setwd('..')
parent_path = getwd()
data_path = file.path(getwd(),
                      "data",
                      "rt_reviews.csv")

#data_path
data = read.csv(data_path,  encoding = "UTF-8", stringsAsFactors=FALSE)
data = subset(data, stri_enc_isascii(data$Review))

ggplot(data, aes(x = Freshness)) +
  geom_bar()

data = sample_n(data, 40000)
#--------------------------------------------------------------------
#-------------------------- Preprocesamiento   ----------------------
#--------------------------------------------------------------------






corpus = Corpus(VectorSource(data$Review[1:40000]))
#print(corpus)
#summary(corpus)
#inspect(corpus[1])
#corpus[[1]]$content
for (i in 1:10) print (corpus[[i]]$content)
corpus = tm_map(corpus, 
                 content_transformer(removePunctuation))

corpus = tm_map(corpus, 
                 content_transformer(removeWords),
                 stopwords("english"))

corpus = tm_map(corpus,
                 content_transformer(tolower))

corpus = tm_map(corpus, 
                content_transformer(removeWords), 
                stopwords("english"))

corpus = tm_map(corpus, stemDocument)

corpus = tm_map(corpus, stripWhitespace) 

corpus = tm_map(corpus, content_transformer(removeNumbers))


matrix = DocumentTermMatrix(corpus)
sparse = as.compressed.matrix(matrix)


#--------------------------------------------------------------------
#---------------------- Sentimental analysis   ----------------------
#--------------------------------------------------------------------



#------------------ Dataset Completo - Experimento 1  ---------------


f = tune.maxent(sparse[1:30000,],
                data$Freshness[1:30000],
                nfold=3,
                showall=TRUE,
                verbose=TRUE)
print(f)

model = maxent(sparse[1:30000,],
               data$Freshness[1:30000], 
               l1_regularizer=0,
               l2_regularizer=1.0, 
               use_sgd=FALSE, 
               set_heldout=0, 
               verbose=TRUE)

results = predict(model,sparse[30001:40000,]) 

freshness_list = data$Freshness[30001:40000]

predicted = results[,1]

df_final = do.call(rbind, Map(data.frame, A=predicted, B=freshness_list))

example <- confusionMatrix(data=as.factor(predicted), reference = as.factor(freshness_list))


#------------------ Dataset Completo - Experimento 2  ---------------


f_2 = tune.maxent(sparse[1:30000,],
                data$Freshness[1:30000],
                nfold=5,
                showall=TRUE,
                verbose=TRUE)
print(f_2)

model_2 = maxent(sparse[1:30000,],
               data$Freshness[1:30000], 
               l1_regularizer=0,
               l2_regularizer=1.0, 
               use_sgd=FALSE, 
               set_heldout=0, 
               verbose=TRUE)

results = predict(model,sparse[30001:40000,]) 

freshness_list = data$Freshness[30001:40000]

predicted = results[,1]

df_final = do.call(rbind, Map(data.frame, A=predicted, B=freshness_list))

example <- confusionMatrix(data=as.factor(predicted), reference = as.factor(freshness_list))







#----------------------- words cloud -------------------




dtm <- TermDocumentMatrix(corpus) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
df <- data.frame(word = names(words),freq=words)


set.seed(1234) # for reproducibility 
wordcloud(words = df$word, 
          freq = df$freq, 
          min.freq = 1,
          max.words=200, 
          random.order=FALSE, 
          rot.per=0.35,            
          colors=brewer.pal(8, "Dark2")
          )




