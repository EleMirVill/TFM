## DATOS DE RNA-Seq

# Instalación de Bioconductor y librerías necesarias
install.packages("BiocManager") # Gestiona paquetes de Bioconductor
BiocManager::install("maftools", force=TRUE)
BiocManager::install("TCGAbiolinks", force = TRUE)
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("annaffy")
BiocManager::install("ipd.Db", force = TRUE)
install.packages("htmlTable")
install.packages("pheatmap")

# Cargar librerías
library(BiocManager)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(limma)
library(edgeR)
library(annaffy)
library(htmlTable)
library(ipdDb)
library(pheatmap)


## RECOPILACION DE DATOS DE EXPRESION (COUNTS) DESDE LA BASE DE DATOS GDC

# Obtener lista de los proyectos disponibles en GDC
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-BRCA') # Resumen del proyecto TCGA-BRCA (CÁNCER DE MAMA)

# Lanzar consulta para recuperar datos de expresión génica
query_TCGA <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Transcriptome Profiling",
                       experimental.strategy = "RNA-Seq",
                       workflow.type = "STAR - Counts",
                       access = "open",
                       barcode = c("TCGA-BH-A18M-11A-33R-A12D-07", "TCGA-BH-A18R-11A-42R-A12D-07", 
                                   "TCGA-BH-A1FE-11B-14R-A13Q-07", "TCGA-BH-A0BZ-11A-61R-A12P-07",
                                   "TCGA-BH-A0DH-11A-31R-A089-07", "TCGA-BH-A1F0-11B-23R-A137-07",
                                   "TCGA-BH-A18V-11A-52R-A12D-07", "TCGA-BH-A18L-11A-42R-A12D-07",
                                   "TCGA-BH-A0DL-11A-13R-A115-07", "TCGA-BH-A18S-11A-43R-A12D-07",
                                   "TCGA-PL-A8LV-01A-21R-A41B-07", "TCGA-BH-A0BC-01A-22R-A084-07",
                                   "TCGA-AR-A1AX-01A-11R-A12P-07", "TCGA-AC-A2FO-01A-11R-A180-07",
                                   "TCGA-AQ-A0Y5-01A-11R-A14M-07", "TCGA-AC-A3EH-01A-22R-A22K-07",
                                   "TCGA-A8-A099-01A-11R-A00Z-07", "TCGA-D8-A140-01A-11R-A115-07",
                                   "TCGA-A8-A097-01A-11R-A034-07", "TCGA-AQ-A04L-01B-21R-A10J-07"))
getResults(query_TCGA) # Muestra la consulta por pantalla

# Descargar datos
GDCdownload(query_TCGA)

# Preprocesamiento de los datos descargados (counts)
SKCM.counts <- GDCprepare(query = query_TCGA, summarizedExperiment = TRUE)
rm(query_TCGA)

# Subir información de los datos descargados (tipos de muestras)
colData <- read.delim("C:/Users/villa/Downloads/id_infoStudy.csv", stringsAsFactors = TRUE, sep = ',')

# Matriz de expresión - counts
counts_data <- assay(SKCM.counts)

# Obtener el nombre del gen id
gene_names <- SKCM.counts@rowRanges@elementMetadata@listData$gene_name
rownames(counts_data) <- gene_names

# Contar el número de genes originales
num_genes_originales <- nrow(counts_data)
cat("Número de genes originales: ", num_genes_originales, "\n")


# Convertir counts a objeto DGEList para limma
dge <- DGEList(counts = counts_data)
dge$samples$group <- colData$sample_type

# Prefiltrar: eliminar filas con pocos conteos
# Definir el umbral de conteo (por ejemplo, un mínimo de conteos totales en todas las muestras)
min_counts <- 10

# Convertir conteos a CPM (counts per million)
cpm_data <- cpm(dge)

# Mantener genes con al menos min_counts en al menos una muestra
keep <- rowSums(cpm_data >= min_counts) >= 1

# Filtrar el objeto DGEList utilizando el índice keep
dge <- dge[keep, , keep.lib.sizes = FALSE]  # Esta opción recalcula el tamaño de las bibliotecas

# Después de la filtración
num_genes_filtrados <- nrow(dge$counts)
cat("Número de genes después de la filtración: ", num_genes_filtrados, "\n")

# Paso 1: Calcular los valores de CPM antes de la normalización
cpm_before_norm <- cpm(dge, log = TRUE)

summary(cpm_before_norm)

# Paso 2: Normalizar los datos usando TMM
dge <- calcNormFactors(dge)

# Paso 3: Calcular los valores de CPM después de la normalización
cpm_after_norm <- cpm(dge, log = TRUE)

summary(cpm_after_norm)

# Paso 4: Crear un dataframe combinado para antes y después de la normalización
cpm_data_combined <- data.frame(
  Sample = rep(colnames(cpm_before_norm), each = nrow(cpm_before_norm) * 2),
  CPM = c(as.vector(cpm_before_norm), as.vector(cpm_after_norm)),
  Condition = rep(c("Antes de Normalización", "Después de Normalización"), times = ncol(cpm_before_norm) * nrow(cpm_before_norm))
)

# Acortar los nombres de las muestras eliminando el prefijo "TCGA-" y parte del sufijo
cpm_data_combined$Sample <- gsub("^TCGA-|-.{4}$", "", cpm_data_combined$Sample)

# Paso 5: Graficar las cajas para cada muestra, mostrando alternadamente Antes y Después
ggplot(cpm_data_combined, aes(x = Sample, y = CPM, fill = Condition)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.75)) +  # Dodge las cajas para separarlas un poco
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),  # Ajuste para que las etiquetas se vean mejor
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  ggtitle("Distribución de CPM Antes y Después de la Normalización por Muestra") +
  ylab("CPM (log scale)") +
  xlab("Muestras") +
  theme_minimal() +
  scale_fill_manual(values = c("Antes de Normalización" = "#00be67", "Después de Normalización" = "#c77cff")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  # Muestra las etiquetas solo una vez por muestra
  coord_flip()  # Girar el gráfico

# Estandarizar los datos en log2 para reducir la variabilidad
countdata_log2 <- log2(cpm_after_norm+ 1)

dim(countdata_log2)

# Ajustar el grupo de referencia para el análisis
# Re-nivelar los factores para establecer el tejido "Normal" como referencia
# Esto permite que el análisis compare la expresión génica en tejidos tumorales contra el tejido normal.
dge$samples$group <- relevel(dge$samples$group, ref = "Solid Tissue Normal")

# Crear el diseño experimental para 10 muestras normales y 10 tumorales
# Se utiliza el factor de condición para distinguir entre muestras "Normal" y "Tumor"
experimental.design <- model.matrix(~ -1 + factor(c(rep("Normal", 10), rep("Tumor", 10))))
colnames(experimental.design) <- c("Normal", "Tumor")

# Ajustar el modelo lineal para cada gen en el conjunto de datos
# Este modelo usa la matriz de diseño para estimar la expresión en "Normal" y "Tumor"
expression.level <- countdata_log2
linear.fit<-lmFit(expression.level, experimental.design)

# Crear una matriz de contraste para comparar expresiones en "Normal" vs "Tumor"
# El contraste "Normal - Tumor" compara directamente los niveles de expresión entre ambos grupos
contrast.matrix <- makeContrasts( Normal - Tumor, levels = c("Normal", "Tumor") )
contrast.matrix

# Aplicar ajuste de Bayes empírico para obtener significancia estadística de los contrastes
contrast.linear.fit<-contrasts.fit(linear.fit, contrast.matrix)
contrast.results<-eBayes(contrast.linear.fit)
nrow(expression.level)

# Extraer los resultados de los genes más diferencialmente expresados
# Se selecciona el coeficiente del contraste "Normal - Tumor" y se ordena por log Fold Change (logFC)
results <- topTable(contrast.results, number = nrow(expression.level), coef = 1, sort.by = "logFC")
head(results)
summary(results)


# Obtener datos de fold-change y p-valor ajustado
fold.change.Breast <- results$logFC
adj.Pval.Breast <- results$adj.P.Val
genes.ids.Breast <- rownames(results)

# Filtrar genes con valores faltantes en fold-change o p-valor ajustado
valid_genes <- !is.na(fold.change.Breast) & !is.na(adj.Pval.Breast)

# Filtrar para obtener datos solo de los genes válidos (con datos completos)
fold.change.Breast_valid <- fold.change.Breast[valid_genes]
adj.Pval.Breast_valid <- adj.Pval.Breast[valid_genes]
genes.ids.Breast_valid <- genes.ids.Breast[valid_genes]

# Calcular el número total de genes válidos con datos completos
num_genes_valid <- length(genes.ids.Breast_valid)

# Ajuste de log2FoldChange y adj. p-valor para obtener menos de 1800 genes
fold_change_thresholds <- c(1.5, 2, 2.5)
p_value_thresholds <- c(0.01, 0.005, 0.001)

# Iterar sobre diferentes umbrales para encontrar una combinación que dé menos de 1800 genes
best_thresholds <- NULL
for (fc in fold_change_thresholds) {
  for (pval in p_value_thresholds) {
    selected_genes <- (abs(fold.change.Breast_valid) > fc) & (adj.Pval.Breast_valid < pval)
    num_selected_genes <- sum(selected_genes, na.rm = TRUE)
    
    if (num_selected_genes <= 1800) {
      best_thresholds <- list(fold_change = fc, p_value = pval, num_genes = num_selected_genes)
      break
    }
  }
  if (!is.null(best_thresholds)) break
}

if (!is.null(best_thresholds)) {
  cat("Mejores umbrales encontrados:\n")
  cat("Fold-change:", best_thresholds$fold_change, "\n")
  cat("P-valor:", best_thresholds$p_value, "\n")
  cat("Número de genes seleccionados:", best_thresholds$num_genes, "\n")
} else {
  cat("No se encontraron umbrales que den menos de 1800 genes.\n")
}

# Crear el objeto breast.all.DEG con los genes seleccionados
selected_genes <- (abs(fold.change.Breast_valid) > best_thresholds$fold_change) & (adj.Pval.Breast_valid < best_thresholds$p_value)
breast.all.DEG <- genes.ids.Breast_valid[selected_genes]

num_genes_selected <- length(breast.all.DEG)

# Obtener expresión promedio
ave_expr <- rowMeans(countdata_log2)

# Crear la tabla de resultados usando topTable
results <- topTable(contrast.results, number = Inf, coef = 1, sort.by = "logFC")
results$AveExpr <- ave_expr[rownames(results)]  # Agregar la expresión promedio a los resultados

# Crear la tabla final
de_table <- data.frame(
  Gene = rownames(results),           # Identificadores de genes
  logFC = results$logFC,              # Fold-change
  AveExpr = results$AveExpr,          # Expresión promedio
  t = results$t,                      # Estadístico t
  p_val = results$P.Value,            # p-valor
  p_val_adj = results$adj.P.Val,      # p-valor ajustado
  B = results$B                       # Log-odds de estar diferencialmente expresado
)

# Número de filas en la tabla de resultados
num_rows <- nrow(de_table)


# Filtrar la tabla de resultados para mantener solo los genes seleccionados
de_table_filtered <- de_table[de_table$Gene %in% breast.all.DEG, ]

# Contar el número de filas en la tabla filtrada
num_rows_filtered <- nrow(de_table_filtered)
cat("Número de filas en de_table_filtered:", num_rows_filtered, "\n")

# Guardar la tabla filtrada en un archivo CSV
write.csv(de_table_filtered, file = "C:/Users/villa/Downloads/filtered_DEGs_table.csv", row.names = FALSE)

# Filtrar la tabla para mantener solo los genes válidos
de_table <- de_table[valid_genes,]
cat("Número de filas en de_table:", num_rows, "\n")

# Guardar la tabla en un archivo CSV
write.csv(de_table, file = "C:/Users/villa/Downloads/limma_results.csv", row.names = FALSE)

# Visualización de los datos: volcano plot
log.padj.breast <- -log10(adj.Pval.Breast_valid)
plot(fold.change.Breast_valid, log.padj.breast, ylab = "-log10(p value)", 
     xlab = "log2 fold change", pch = 19, cex = 0.5, col = "grey", xlim = c(-6, 6))

# Genes activados (puntos rojos)
activated_genes <- (fold.change.Breast_valid > best_thresholds$fold_change) & 
  (adj.Pval.Breast_valid < best_thresholds$p_value)

# Verificar cuántos genes activados hay
cat("Número de genes activados:", sum(activated_genes), "\n")

# Dibujar los puntos rojos para genes activados
#points(fold.change.Breast_valid[activated_genes], log.padj.breast[activated_genes], col = "red", pch = 19, cex = 0.5)

# Dibujar los puntos rojos para genes activados
points(fold.change.Breast_valid[activated_genes], log.padj.breast[activated_genes], 
       pch = 19, cex = 1, col = "red")

# Genes reprimidos (puntos azules)
repressed_genes <- (fold.change.Breast_valid < -best_thresholds$fold_change) & 
  (adj.Pval.Breast_valid < best_thresholds$p_value)

# Verificar cuántos genes reprimidos hay
cat("Número de genes reprimidos:", sum(repressed_genes), "\n")

# Dibujar los puntos azules para genes reprimidos
points(fold.change.Breast_valid[repressed_genes], log.padj.breast[repressed_genes], 
       pch = 19, cex = 1, col = "blue")

# Líneas de referencia para p-valor y fold-change
abline(h = -log10(best_thresholds$p_value), col = "blue", lty = 2)  # Línea horizontal en -log10(p_value)
abline(v = c(-best_thresholds$fold_change, best_thresholds$fold_change), col = "blue", lty = 2)  # Líneas verticales en fold-change

# Seleccionar los genes más separados (más extremos en ambos ejes)
num_genes_to_label <- 5  # Número de genes que quieres etiquetar

# Para genes activados
if (sum(activated_genes) > 0) {
  # Ordenar genes activados por la combinación de fold-change y -log10(p-value)
  extreme_activated_indices <- order(
    abs(fold.change.Breast_valid[activated_genes]) * log.padj.breast[activated_genes], 
    decreasing = TRUE
  )[1:min(sum(activated_genes), num_genes_to_label)]
  
  # Verificar si hay genes para etiquetar
  if (length(extreme_activated_indices) > 0) {
    # Etiquetar los genes activados más extremos
    text(
      fold.change.Breast_valid[activated_genes][extreme_activated_indices] + 0.2, 
      log.padj.breast[activated_genes][extreme_activated_indices], 
      labels = genes.ids.Breast_valid[activated_genes][extreme_activated_indices], 
      pos = 4, cex = 0.7, col = "red", srt = 45  # Ángulo de rotación para mejor visibilidad
    )
  }
}

# Para genes reprimidos
if (sum(repressed_genes) > 0) {
  # Ordenar genes reprimidos por la combinación de fold-change y -log10(p-value)
  extreme_repressed_indices <- order(
    abs(fold.change.Breast_valid[repressed_genes]) * log.padj.breast[repressed_genes], 
    decreasing = TRUE
  )[1:min(sum(repressed_genes), num_genes_to_label)]
  
  # Verificar si hay genes para etiquetar
  if (length(extreme_repressed_indices) > 0) {
    # Etiquetar los genes reprimidos más extremos
    text(
      fold.change.Breast_valid[repressed_genes][extreme_repressed_indices] - 0.2, 
      log.padj.breast[repressed_genes][extreme_repressed_indices], 
      labels = genes.ids.Breast_valid[repressed_genes][extreme_repressed_indices], 
      pos = 2, cex = 0.7, col = "blue", srt = 45  # Ángulo de rotación para mejor visibilidad
    )
  }
}

#Volcano2 mejor este que el anterior
# Crear un dataframe para el Volcano Plot
library(dplyr)
library(ggplot2)

volcano_data <- data.frame(
  logFC = results$logFC,
  negLogPval = -log10(results$adj.P.Val),
  Significance = ifelse(
    results$adj.P.Val < 0.01 & abs(results$logFC) > 1.5,
    ifelse(results$logFC > 1.5, "Activado", "Reprimido"),
    "No significativo"
  ),
  Gene = rownames(results)  # Añadir los nombres de los genes como una columna
)

# Generar el Volcano Plot con ggplot2
base_plot <- ggplot(volcano_data, aes(x = logFC, y = negLogPval, color = Significance)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(
    values = c("Activado" = "red", "Reprimido" = "blue", "No significativo" = "grey"),
    name = "Significancia"  # Cambia el título de la leyenda
  ) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot de Genes Diferencialmente Expresados",
    x = "Log2 Fold Change",
    y = "-log10(P-valor ajustado)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Seleccionar los genes activados y reprimidos más significativos para etiquetar
top_activated_genes <- volcano_data %>%
  filter(Significance == "Activado") %>%
  arrange(desc(negLogPval)) %>%
  head(10)

top_repressed_genes <- volcano_data %>%
  filter(Significance == "Reprimido") %>%
  arrange(desc(negLogPval)) %>%
  head(10)

# Añadir las etiquetas al gráfico
base_plot +
  geom_text(
    data = top_activated_genes,
    aes(label = Gene),
    color = "red",
    angle = 40,  # Rotación de las etiquetas para mejorar legibilidad
    hjust = 1,
    vjust = 1,
    size = 2.5
  ) +
  geom_text(
    data = top_repressed_genes,
    aes(label = Gene),
    color = "blue",
    angle = 30,  # Rotación de las etiquetas para mejorar legibilidad
    hjust = 1,
    vjust = 1,
    size = 2.5
  )

# Contar el número de genes en cada categoría
num_genes_by_category <- table(volcano_data$Significance)
print(num_genes_by_category)

#heatmap

# Cargar datos de expresión para muestras normales desde el archivo TXT
#normal_data <- read.delim("C:/Users/villa/Downloads/NormalBreast_DEGs_table2.txt", sep = "\t", row.names = 1)

# Cargar datos de expresión para muestras tumorales desde el archivo TXT
#tumor_data <- read.delim("C:/Users/villa/Downloads/TumorBreast_DEGs_table2.txt", sep = "\t", row.names = 1)

# Eliminar genes con valores NA en ambos datasets
#normal_data_clean <- na.omit(normal_data)
#tumor_data_clean <- na.omit(tumor_data)

# Obtener los nombres de los genes que están presentes en ambos conjuntos de datos
#common_genes <- intersect(rownames(normal_data_clean), rownames(tumor_data_clean))

# Filtrar ambos conjuntos de datos para incluir solo los genes comunes
#normal_data_clean <- normal_data_clean[common_genes, ]
#tumor_data_clean <- tumor_data_clean[common_genes, ]

# Combinar los datos limpios y filtrados de normal y tumor
#combined_data_clean <- cbind(normal_data_clean, tumor_data_clean)

# Cargar la librería pheatmap
library(pheatmap)

# Generar el heatmap para los datos combinados
#pheatmap(
#  combined_data_clean,
#  scale = "row",
#  clustering_distance_rows = "euclidean",
#  clustering_distance_cols = "euclidean",
#  clustering_method = "complete",
#  show_rownames = FALSE,
#  show_colnames = TRUE,
#  color = colorRampPalette(c("blue", "white", "red"))(50),
#  main = "Mapa de Calor de Expresión Génica - Genes Comunes en Muestras Normales y Tumorales"
#)



# Definir umbrales de filtrado
fc_threshold <- 1.5
pval_threshold <- 0.01

# Genes activados
activated_genes <- results %>%
  filter(logFC > fc_threshold & adj.P.Val < pval_threshold) %>%
  arrange(desc(logFC)) %>%
  head(10)  # Los 10 genes más activados

# Genes reprimidos
repressed_genes <- results %>%
  filter(logFC < -fc_threshold & adj.P.Val < pval_threshold) %>%
  arrange(logFC) %>%
  head(10)  # Los 10 genes más reprimidos


# Crear tabla para genes activados
activated_table <- data.frame(
  Gene = rownames(activated_genes),
  `Log2 Fold Change` = activated_genes$logFC,
  `-log10(P-valor ajustado)` = -log10(activated_genes$adj.P.Val)
)

# Crear tabla para genes reprimidos
repressed_table <- data.frame(
  Gene = rownames(repressed_genes),
  `Log2 Fold Change` = repressed_genes$logFC,
  `-log10(P-valor ajustado)` = -log10(repressed_genes$adj.P.Val)
)

# Visualizar ambas tablas
print("Tabla de los 10 genes más activados:")
print(activated_table)

print("Tabla de los 10 genes más reprimidos:")
print(repressed_table)

#segunda tabla
# Definir umbrales de filtrado
fc_threshold <- 1.5
pval_threshold <- 0.01

# Genes activados
activated_genes <- results %>%
  filter(logFC > fc_threshold & adj.P.Val < pval_threshold) %>%
  arrange(desc(logFC)) %>%
  head(10)  # Los 10 genes más activados

# Genes reprimidos
repressed_genes <- results %>%
  filter(logFC < -fc_threshold & adj.P.Val < pval_threshold) %>%
  arrange(logFC) %>%
  head(10)  # Los 10 genes más reprimidos

# Asegurar que los nombres de genes se obtienen correctamente
activated_genes$Gene <- rownames(activated_genes)
repressed_genes$Gene <- rownames(repressed_genes)

# Crear tabla para genes activados
activated_table <- data.frame(
  Gene = activated_genes$Gene,
  `Log2 Fold Change` = activated_genes$logFC,
  `-log10(P-valor ajustado)` = -log10(activated_genes$adj.P.Val)
)

# Crear tabla para genes reprimidos
repressed_table <- data.frame(
  Gene = repressed_genes$Gene,
  `Log2 Fold Change` = repressed_genes$logFC,
  `-log10(P-valor ajustado)` = -log10(repressed_genes$adj.P.Val)
)

# Visualizar ambas tablas
cat("Tabla de los 10 genes más activados:\n")
print(activated_table)

cat("\nTabla de los 10 genes más reprimidos:\n")
print(repressed_table)


# Mostrando el número total de genes seleccionados
cat("Número total de genes seleccionados: ", length(breast.all.DEG), "\n")






# Cambio nombres columnas
colnames(countdata_log2) <- colData$name_id
normal.breast.DEG.table <- countdata_log2[, c("Normal_breast_1", "Normal_breast_2", 
                                              "Normal_breast_3", "Normal_breast_4",
                                              "Normal_breast_5", "Normal_breast_6",
                                              "Normal_breast_7", "Normal_breast_8", 
                                              "Normal_breast_9", "Normal_breast_10")]


normal.breast.DEG.table <- normal.breast.DEG.table[rownames(normal.breast.DEG.table) %in% breast.all.DEG,]
normal.breast.DEG.table <- cbind(attr_name = rownames(normal.breast.DEG.table), normal.breast.DEG.table)

tumor.breast.DEG.table <- countdata_log2[, c("Tumor_breast_1", "Tumor_breast_2", 
                                             "Tumor_breast_3", "Tumor_breast_4", 
                                             "Tumor_breast_5", "Tumor_breast_6",
                                             "Tumor_breast_7", "Tumor_breast_8",
                                             "Tumor_breast_9", "Tumor_breast_10")]
tumor.breast.DEG.table <- tumor.breast.DEG.table[rownames(tumor.breast.DEG.table) %in% breast.all.DEG,]
tumor.breast.DEG.table <- cbind(attr_name = rownames(tumor.breast.DEG.table), tumor.breast.DEG.table)


# Guardar la tabla filtrada para muestras normales
write.table(normal.breast.DEG.table, file="C:/Users/villa/Downloads/NormalBreast_DEGs_table2.csv", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
html_normal_breast_DEGs <- htmlTable(normal.breast.DEG.table, 
                                     caption = "Normal Breast DEGs",
                                     rnames = FALSE)
write(html_normal_breast_DEGs, file = "C:/Users/villa/Downloads/NormalBreast_DEGs_table2.html")

write.table(normal.breast.DEG.table, file = "C:/Users/villa/Downloads/NormalBreast_DEGs_table2.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

# Guardar la tabla filtrada para muestras tumorales
write.table(tumor.breast.DEG.table, file="C:/Users/villa/Downloads/TumorBreast_DEGs_table2.csv", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(tumor.breast.DEG.table, file="C:/Users/villa/Downloads/TumorBreast_DEGs_table2.txt", sep = "\t", 
           quote = FALSE, row.names = FALSE, col.names = TRUE)



####analisis de enriquecimiento3#####

# Cargar los genes desde el archivo "TumorBreast_DEGs_table2.csv"
tumor_breast_deg_table <- read.csv("C:/Users/villa/Downloads/TumorBreast_DEGs_table2.csv", sep = "\t", header = TRUE)

# Asegúrate de extraer la columna con los nombres de los genes
genes <- as.character(tumor_breast_deg_table$attr_name)  # Ajusta 'attr_name' al nombre exacto de la columna que contiene los IDs de genes

# Verificar los genes cargados
head(genes)  # Muestra los primeros genes para confirmar la carga correcta

# Cargar las librerías necesarias para el análisis de enriquecimiento funcional
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
library(clusterProfiler)
library(org.Hs.eg.db)

# Realizar análisis de enriquecimiento de GO y KEGG
# Convertir los IDs de genes a un formato compatible
genes <- as.character(breast.all.DEG)  # Asegúrate de que los genes estén en formato de texto

# Análisis de enriquecimiento GO
enrich_go <- enrichGO(gene = genes,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",  # Ajusta esto según el tipo de ID en breast.all.DEG
                      ont = "BP",          # Ontología de procesos biológicos
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

# Análisis de enriquecimiento KEGG
genes_entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Extrae solo los Entrez IDs para el enriquecimiento
entrez_ids <- genes_entrez$ENTREZID
# Realizar el enriquecimiento en KEGG usando los Entrez IDs
enrich_kegg <- enrichKEGG(gene = entrez_ids,     # Aquí usamos entrez_ids en lugar de genes
                          organism = "hsa",      # Homo sapiens
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)

# Revisa los resultados
print(enrich_kegg)


# Guardar los resultados en archivos CSV
write.csv(as.data.frame(enrich_go), file = "C:/Users/villa/Downloads/GO_enrichment_results.csv", row.names = FALSE)
write.csv(as.data.frame(enrich_kegg), file = "C:/Users/villa/Downloads/KEGG_enrichment_results.csv", row.names = FALSE)

# Generar visualizaciones
library(ggplot2)
dotplot(enrich_go, showCategory = 10) + ggtitle("GO Enrichment (Biological Processes)")
dotplot(enrich_kegg, showCategory = 10) + ggtitle("KEGG Pathway Enrichment")


#####NO SE PARA QUE SIRVE#####

# Instalar Bioconductor después de actualizar R
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Instalar el paquete ipdDb
BiocManager::install("ipdDb", force=TRUE)


# Cargar el paquete
library(ipdDb)

tumor.breast.DEG.table_df <- as.data.frame(tumor.breast.DEG.table)
# Aquí se supone que necesitas usar un objeto de anotación
# Supongamos que `aaf.handler()` es una función para manipular datos
annotated_data <- aaf.handler(tumor.breast.DEG.table_df, "ipd.Db")
annotated_data <- aaf.handler(tumor.breast.DEG.table_df, name = "ipdDb")

head(tumor.breast.DEG.table_df)

# Extraer los identificadores de las sondas
probeids <- tumor.breast.DEG.table_df$attr_name  # Esto asumiendo que 'attr_name' tiene los IDs que necesitas

# Definir el nombre del chip (asegúrate de que sea correcto)
chip_name <- "ipdDb"  # Ajusta esto según el nombre correcto del chip que estás utilizando

# Llamar a aaf.handler() con los argumentos correctos
annotated_data <- aaf.handler(probeids, chip = chip_name, name = "ipdDb")

# Verificar el resultado
if (is.null(annotated_data)) {
  print("annotated_data es NULL. Verifica la entrada o la función.")
} else {
  print(head(annotated_data))  # Muestra las primeras filas de los datos anotados
}


# Crear la tabla HTML con los datos anotados
html_tumor_breast_DEGs <- htmlTable(annotated_data, 
                                    caption = "Tumor Breast DEGs",
                                    rnames = FALSE)


html_tumor_breast_DEGs<-aafTableAnn(tumor.breast.DEG.table, 
                           caption = "Tumor Breast DEGs",
                           rnames = FALSE)
write(html_tumor_breast_DEGs, file = "C:/Users/villa/Downloads/TumorBreast_DEGs_table2.html")

write.table(tumor.breast.DEG.table, file = "C:/Users/villa/Downloads/TumorBreast_DEGs_table2.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)


# Verificar los nombres de genes en countdata_log2
genes_in_data <- rownames(countdata_log2)

# Verificar los nombres de genes en breast.all.DEG
missing_genes <- setdiff(breast.all.DEG, genes_in_data)
if (length(missing_genes) > 0) {
  cat("Genes en breast.all.DEG que no están en countdata_log2:\n")
  print(missing_genes)
} else {
  cat("Todos los genes en breast.all.DEG están presentes en countdata_log2.\n")
}

# Mostrar el número de filas y las primeras filas de la tabla
cat("Número de filas en tumor.breast.DEG.table:", nrow(tumor.breast.DEG.table), "\n")
if (nrow(tumor.breast.DEG.table) > 0) {
  head(tumor.breast.DEG.table)
} else {
  cat("La tabla tumor.breast.DEG.table está vacía.\n")
}

# Verificar el número de genes en breast.all.DEG
cat("Número total de genes en breast.all.DEG:", length(breast.all.DEG), "\n")

# Comprobar si los nombres de los genes se están alineando correctamente
head(breast.all.DEG)
head(rownames(countdata_log2))


if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("biomaRt")  # Asegúrate de tener biomaRt instalado
library(biomaRt)
library(htmlTable)  # Necesario para crear tablas HTML

# Conectar a Ensembl
mart <- useMart("ensembl")
dataset <- useDataset("hsapiens_gene_ensembl", mart)

# Asumiendo que tumor.breast.DEG.table tiene una columna 'attr_name' con los gene IDs
tumor.breast.DEG.table_df <- as.data.frame(tumor.breast.DEG.table)
probeids <- tumor.breast.DEG.table_df$attr_name

# Obtener anotaciones
annotated_data <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),
                        filters = "ensembl_gene_id",  # Cambia si 'probeids' son otros identificadores
                        values = probeids, 
                        mart = dataset)

# Verificar el resultado de las anotaciones
if (nrow(annotated_data) == 0) {
  print("annotated_data está vacío. Verifica los identificadores o el dataset.")
} else {
  print(head(annotated_data))  # Muestra las primeras filas de los datos anotados
}

# Crear la tabla HTML con los datos anotados
html_tumor_breast_DEGs <- htmlTable(annotated_data, 
                                    caption = "Tumor Breast DEGs",
                                    rnames = FALSE)

# Guardar la tabla HTML
write(html_tumor_breast_DEGs, file = "C:/Users/villa/Downloads/TumorBreast_DEGs_table4.html")

# Guardar los datos anotados en un archivo TXT
write.table(annotated_data, file = "C:/Users/villa/Downloads/TumorBreast_DEGs_table4.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

# Guardar los datos anotados en un archivo CSV
write.csv(annotated_data, file = "C:/Users/villa/Downloads/TumorBreast_DEGs_table4.csv", 
          row.names = FALSE)

# Verificar los nombres de genes en countdata_log2
genes_in_data <- rownames(countdata_log2)

# Verificar los nombres de genes en breast.all.DEG
missing_genes <- setdiff(breast.all.DEG, genes_in_data)
if (length(missing_genes) > 0) {
  cat("Genes en breast.all.DEG que no están en countdata_log2:\n")
  print(missing_genes)
} else {
  cat("Todos los genes en breast.all.DEG están presentes en countdata_log2.\n")
}

# Mostrar el número de filas y las primeras filas de la tabla original
cat("Número de filas en tumor.breast.DEG.table:", nrow(tumor.breast.DEG.table), "\n")
if (nrow(tumor.breast.DEG.table) > 0) {
  head(tumor.breast.DEG.table)
} else {
  cat("La tabla tumor.breast.DEG.table está vacía.\n")
}

# Verificar el número de genes en breast.all.DEG
cat("Número total de genes en breast.all.DEG:", length(breast.all.DEG), "\n")

# Comprobar si los nombres de los genes se están alineando correctamente
head(breast.all.DEG)
head(rownames(countdata_log2))

# Verificar qué genes de breast.all.DEG no están en la tabla de anotación
missing_genes <- setdiff(breast.all.DEG, annotated_data$external_gene_name)
cat("Genes en breast.all.DEG que no están en la tabla de anotación:\n")
print(missing_genes)

# Verificar cuántos genes de la tabla de anotación no están en breast.all.DEG
not_in_deg <- setdiff(annotated_data$external_gene_name, breast.all.DEG)
cat("Genes en la tabla de anotación que no están en breast.all.DEG:\n")
print(not_in_deg)

# Primero, obtener los IDs de Ensembl usando external_gene_name como filtro
probeids <- breast.all.DEG  # Usar los nombres de los genes en breast.all.DEG

# Obtener los IDs de Ensembl
ensembl_data <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = breast.all.DEG,
                      mart = dataset)

# Asegúrate de revisar qué IDs has obtenido
print(ensembl_data)

# Luego, puedes utilizar esos IDs para la anotación que deseas
if (nrow(ensembl_data) > 0) {
  # Ahora puedes usar estos IDs de Ensembl para obtener la descripción
  annotated_data <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),
                          filters = "ensembl_gene_id",
                          values = ensembl_data$ensembl_gene_id,
                          mart = dataset)
  
  # Verificar resultados
  if (nrow(annotated_data) == 0) {
    print("annotated_data está vacío.")
  } else {
    print(head(annotated_data))  # Mostrar las primeras filas de los datos anotados
  }
} else {
  print("No se encontraron IDs de Ensembl.")
}


# Genes en breast.all.DEG que no están en la tabla de anotación
missing_genes <- setdiff(breast.all.DEG, annotated_data$external_gene_name)
cat("Genes en breast.all.DEG que no están en la tabla de anotación:\n")
print(missing_genes)

# Genes en la tabla de anotación que no están en breast.all.DEG
not_in_deg <- setdiff(annotated_data$external_gene_name, breast.all.DEG)
cat("Genes en la tabla de anotación que no están en breast.all.DEG:\n")
print(not_in_deg)



