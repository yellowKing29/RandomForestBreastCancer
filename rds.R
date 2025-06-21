library(Matrix)
library(tidyr)
library(ggplot2)
library(readr)
library(dplyr)
library(caret)
library(ranger)
library(glmnet)


# CARGA DE DATOS
matriz <- readRDS("BreastCancer/1863-counts_cells_cohort1.Rds")
metadata <- read_csv("BreastCancer/1872-BIOKEY_metaData_cohort1_web.csv")


#VISUALIZACION GENERAL DE DATOS
#str(matriz)
#summary(matriz)

cell_names <- colnames(matriz)


# IDENTIFICAR COLUMNAS CON "_Pre_" (TRUE) y "_On_" (FALSE)
keep_cells <- grepl("_Pre_", cell_names, fixed = TRUE)


# FILTRAR LA MATRIZ
matriz_filtrada <- matriz[, keep_cells]
#str(matriz_filtrada)

Columnas <- colnames(matriz_filtrada)
#print(Columnas)


#NORMALIZACION A 10,000
norm1 <- 1e4
Xnorm <- matriz_filtrada %*% Diagonal(x = norm1 / colSums(matriz_filtrada))
#str(Xnorm)


#NORMALIZACION LOG
Xlog <- log1p(Xnorm)
#head(Xlog)
#str(Xlog)

#minVal <- min(Xlog)
#maxVal <- max(Xlog)
#print(minVal)
#print(maxVal)

#min_val_nz <- min(Xlog@x)  # Solo valores no cero
#max_val_nz <- max(Xlog@x)
#print(min_val_nz)
#print(max_val_nz)


#CALCULAR EXPRESION POR GEN
gene_sums <- rowSums(Xlog)


#SEPARACION DE GENENES ≥ 1
genes_to_keep <- which(gene_sums >= 1)
X_filtered <- Xlog[genes_to_keep, ]
#str(X_filtered)


# NUMERO MINIMO DE CELULAS DE UN GEN
n_cells_min <- 10


# CALCULAR EN CUANTAS COLUMNAS CADA GEN TIENE EXPRESION NO NULA
genes_expressed_cells <- rowSums(X_filtered != 0)


# FILTRAR GENES CON AL MENOS `n_cells_min` CELULAS CON EXPRESION
genes_to_keep <- which(genes_expressed_cells >= n_cells_min)
X_filtered2 <- X_filtered[genes_to_keep, ]
#str(X_filtered2)


# FILTRAR METADATA CUIDANDO QUE COINCIDA CON SPARSE MATRIX
metadata_filtered <- metadata %>%
  filter(Cell %in% Columnas) %>%
  arrange(match(Cell, Columnas))  # Asegura el mismo orden


# EXTRAER VECTOR ALINEADO
expansion <- metadata_filtered$expansion


# REEMPLAZAR n/a por NE
expansion[expansion == "n/a"] <- NA
expansion[is.na(expansion)] <- "NE"
#print(expansion)


#str(X_filtered2)
#str(expansion)

set.seed(123)

# Crear índice estratificado (usa etiquetas "expansion")
train_idx <- createDataPartition(expansion, p = 0.7, list = FALSE)

# Normalizar tipo de dato
train_idx <- as.vector(train_idx)  # o también: train_idx[, 1]


# Dividir matriz de expresión (dgCMatrix)
X_train <- X_filtered2[, train_idx]
X_test  <- X_filtered2[, -train_idx]

# Dividir vector de etiquetas
y_train <- expansion[train_idx]
y_test  <- expansion[-train_idx]

#str(X_train)
#str(X_test)
#str(y_train)
#str(y_test)

#INICIA RANDOM FOREST

# Transponer para que cada fila sea una célula (observación)
X_train_t <- Matrix::t(X_train)
X_test_t  <- Matrix::t(X_test)

# Asegura que y esté en formato factor binario o multiclase
y_train_factor <- as.factor(y_train)
y_test_factor <- as.factor(y_test)

# Entrenamiento (multinomial si hay más de 2 clases)
set.seed(123)
modelo_glmnet <- cv.glmnet(
  x = X_train_t,
  y = y_train_factor,
  family = ifelse(length(unique(y_train_factor)) > 2, "multinomial", "binomial"),
  alpha = 1,             # LASSO (usa alpha = 0 para Ridge)
  type.measure = "class" # usa clasificación como métrica
)

# Predicción
predicciones <- predict(modelo_glmnet, newx = X_test_t, s = "lambda.min", type = "class")

# Evaluación
confusionMatrix(as.factor(predicciones), y_test_factor)












#cell_totals <- colSums(X_filtered2)
# Filtrar células con al menos cierta expresión total
#min_counts <- 100
#cells_to_keep <- which(cell_totals >= min_counts)
# Aplicar el filtro
#X_final <- X_filtered2[, cells_to_keep]
#str(X_final)

get_nonzero_per_column <- function(sparse_mat) {
  result <- vector("list", ncol(sparse_mat))
  for (j in seq_len(ncol(sparse_mat))) {
    idx_start <- sparse_mat@p[j] + 1
    idx_end   <- sparse_mat@p[j + 1]
    if (idx_end >= idx_start) {
      result[[j]] <- sparse_mat@x[idx_start:idx_end]
    } else {
      result[[j]] <- numeric(0)  # columna con solo ceros
    }
  }
  result
}

#BOXPLOT

# Extraer datos no cero por columna
#col_data <- get_nonzero_per_column(Xlog)
#quartz()
# Ahora hacer boxplot solo con columnas seleccionadas (por ejemplo las 1:50)
#boxplot(col_data[1:50],
#        outline = FALSE,
#        las = 2,
#        main = "Distribución por columna (valores ≠ 0)",
#        ylab = "Expresión")
#Sys.sleep(120)




#Suma filas despreciando los na 
#rowSums(Xlog, na.rm = TRUE)


