## ----setup, include=FALSE-----------------------------------------------------
options(knitr.duplicate.label = "allow")

# Establecemos la configuración por defecto del documento
knitr::opts_chunk$set(echo = TRUE, 
               comment = NA, 
               prompt = TRUE, 
               tidy = FALSE, 
               fig.width = 7, 
               fig.height = 7, 
               message = FALSE, 
               warning = FALSE, 
               cache = FALSE)
Sys.setlocale("LC_TIME", "C")


## ----paquetes, include=FALSE--------------------------------------------------
# Instalamos los paquetes necesarios para resolver los ejercicios
if(!(require(BiocManager, quietly = TRUE))) install.packages("BiocManager")
if(!(require(ggtext, ))) install.packages("ggtext")
if(!(require(FactoMineR, ))) install.packages("FactoMineR")
if(!(require(factoextra, ))) install.packages("factoextra")
if(!(require(Rcpp, ))) install.packages("Rcpp")
if(!(require(printr))) {
  install.packages(
    'printr',
    type = 'source',
    repos = c('http://yihui.name/xran', 'http://cran.rstudio.com')
  )
}
BiocManager::install("SummarizedExperiment")
BiocManager::install("POMA")
BiocManager::install("Biobase")


## -----------------------------------------------------------------------------
# Valores clínicos y metabolómicos
data_values <- read.csv("C:\\Users\\juanm\\Desktop\\MASTER_BIOINFORMATICA\\Analisis_Datos_Omicos\\UOC_metabodata\\Datasets\\2018-MetabotypingPaper\\DataValues_S013.csv", row.names = 1) 

# Metadatos de columnas (de las variables)
data_info <- read.csv("C:\\Users\\juanm\\Desktop\\MASTER_BIOINFORMATICA\\Analisis_Datos_Omicos\\UOC_metabodata\\Datasets\\2018-MetabotypingPaper\\DataInfo_S013.csv", row.names = 1)    

# Metadatos adicionales de los metabolitos
aa_information <- read.csv("C:\\Users\\juanm\\Desktop\\MASTER_BIOINFORMATICA\\Analisis_Datos_Omicos\\UOC_metabodata\\Datasets\\2018-MetabotypingPaper\\AAInformation_S006.csv")


## -----------------------------------------------------------------------------
library(SummarizedExperiment)


## -----------------------------------------------------------------------------
# Dimensiones de los datos del assay
cat("Dimensiones de data_values (assays):", dim(data_values), "\n")

# Dimensiones de colData
cat("Dimensiones de data_info (colData):", dim(data_info), "\n")

# Dimensiones de los metadatos de los metabolitos
cat("Dimensiones de aa_information (metadatos$metabolitos):", dim(aa_information), "\n")


## -----------------------------------------------------------------------------
head(data_values[, 1:9]) # Seleccionamos solo las 9 primeras columnas


## -----------------------------------------------------------------------------
library(dplyr)
# Transponemos la matriz
data_values <- as.data.frame(t(data_values))

# Dataset de metadatos de la componente "colData" (primeras 9 columnas)
metadatos_colData <- data_values[1:9, ]

# Dataset de la componente assays (resto de las columnas) 
assay_Data <- as.data.frame(data_values[-c(1:9),])
# Transformamos sus variables de nuevo a tipo "double"
assay_Data <- assay_Data %>% mutate_if(is.character, as.numeric)

# Comprobamos que los datasets se han dividido correctamente
head(metadatos_colData[, 1:9]) # Seleccionamos los 9 primeros sujetos
head(assay_Data[, 1:9])


## -----------------------------------------------------------------------------
# Dataset de metadatos adicionales para las primeras 9 filas
metadatos_fila9_info <- data_info[1:9, ]

# Dataset de colData (resto de las columnas)
rowData_info <- data_info[-c(1:9), ]

# Comprobamos los datasets
head(metadatos_fila9_info)
head(rowData_info)


## -----------------------------------------------------------------------------
                                                    # Matriz de datos 
se <- SummarizedExperiment(assays = list(counts = as.matrix(assay_Data)),
                                                    # Metadatos de columnas (sujetos)
                           colData = t(metadatos_colData),
                                                    # Metadatos de filas (variables)     
                           rowData = as.data.frame(rowData_info))           
 
# Añadimos los metadatos adicionales                       
metadata(se)$metabolitos <- aa_information
metadata(se)$metadata_adicional <- metadatos_fila9_info

# Visualizamos el objeto `SummarizedExperiment`
print(se)


## -----------------------------------------------------------------------------
summary(se)


## -----------------------------------------------------------------------------
dim(se)


## -----------------------------------------------------------------------------
dimnames(se[1:10, 1:10])


## -----------------------------------------------------------------------------
colData(se)


## -----------------------------------------------------------------------------
rowData(se)


## -----------------------------------------------------------------------------
head(metadata(se)$metabolitos)


## -----------------------------------------------------------------------------
metadata(se)$metadata_adicional


## -----------------------------------------------------------------------------
# Dimensiones de la Matriz
dim(se)

# Datos de las primeras 10 variables de los tres primeros pacientes
assay(se[1:10, 1:3])

# Datos de las primeras 5 variables del segundo paciente
assay(se[1:5, 2])



## -----------------------------------------------------------------------------
boxplot(assay(se, "counts")[, 1:10], 
        main = "Boxplot de los primeros 10 sujetos", 
        ylab = "Valores Metabólicos", las = 2)


## -----------------------------------------------------------------------------
library("POMA")
library("Rcpp")
library("ggtext")
library("magrittr")
library("dplyr")


## -----------------------------------------------------------------------------
se_imputado <- se %>% PomaImpute(method = "median", 
                                 zeros_as_na = FALSE, 
                                 remove_na = TRUE, 
                                 cutoff = 60)


## -----------------------------------------------------------------------------
correlation_matrix <- cor(assay(se_imputado))
heatmap(correlation_matrix, 
        main = "Mapa de calor de la correlación entre muestras")


## -----------------------------------------------------------------------------
PomaBoxplots(se_imputado, 
             x = "samples")


## -----------------------------------------------------------------------------
se_normal <- se_imputado %>% 
  PomaNorm(method = "auto_scaling")


## -----------------------------------------------------------------------------
PomaBoxplots(se_normal, 
             x = "samples")


## -----------------------------------------------------------------------------
library(ggplot2)


## -----------------------------------------------------------------------------
# Reemplazamos valores negativos en el assay del objeto `se_imputado` por 0
assay(se_imputado)[assay(se_imputado) < 0] <- 0

# Comprobamos que ya no hay observaciones menores que 0
(negative_rows <- rownames(assay(se_imputado))[rowSums(assay(se_imputado) < 0) > 0])


## -----------------------------------------------------------------------------
# Normalizamos los datos
se_normal <- se_imputado %>% 
  PomaNorm(method = "auto_scaling")

# Extraemos la matriz del assay (datos metabolómicos)
se_matrix <- t(assay(se_normal)) # Trasponemos para poder graficar por "SURGERY" 

# Realizamos el PCA
pca_result <- prcomp(se_matrix, center = TRUE, scale. = TRUE)


## -----------------------------------------------------------------------------
library(ggplot2)

# Convertimos los resultados del PCA a un data frame para poder graficarlos
pca_data <- as.data.frame(pca_result$x) # Contiene las componentes principales
# Añadimos la condición para colorear según tipo de cirugía
pca_data$Condicion <- colData(se)$SURGERY  

# Calculamos la proporción de varianza explicada por cada CP
varianza_explicada <- pca_result$sdev^2 / sum(pca_result$sdev^2)
# Convertimos a porcentaje
porcent_var_expl <- round(varianza_explicada * 100, 1)

# Graficamos las primeras dos componentes principales
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condicion)) +
  geom_point() +
  labs(title = "PCA de datos metabolómicos",
       x = paste0("PC1 (", porcent_var_expl[1], "% varianza)"),
       y = paste0("PC2 (", porcent_var_expl[2], "% varianza)")) +
  theme_minimal()



## -----------------------------------------------------------------------------
porcent_var_expl


## -----------------------------------------------------------------------------
library("FactoMineR")
library("factoextra")
se_pca <- PCA(se_matrix, graph = FALSE)
fviz_eig(se_pca, addlabels = TRUE, ylim = c(0, 15))


## -----------------------------------------------------------------------------
save(se, file = "C:/Users/juanm/Desktop/MASTER_BIOINFORMATICA/Analisis_Datos_Omicos/ADO_PEC1/datos_se.Rda")


## -----------------------------------------------------------------------------
library(knitr)
purl("C:/Users/juanm/Desktop/MASTER_BIOINFORMATICA/Analisis_Datos_Omicos/ADO_PEC1/ADO_PEC1.Rmd", 
     output = "codigo_exploracion_datos.R", documentation = 1)

