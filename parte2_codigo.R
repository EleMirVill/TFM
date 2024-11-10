# Cargar las librerías necesarias solo si no están instaladas
if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("infotheo")
if (!requireNamespace("impute", quietly = TRUE)) {
  BiocManager::install("impute")
}
install.packages("pROC")


# Cargar las librerías
library(clusterProfiler)
library(org.Hs.eg.db)
library(igraph)
library(ggplot2)
library(infotheo)  # Para NMI
library(randomForest)
library(e1071)  # Para SVM
library(impute)
library(pROC)

# Paso 1: Cargar los datos de expresión
tumor_data <- read.csv("C:/Users/villa/Downloads/TumorBreast_DEGs_table2.csv", sep = "\t", header = TRUE, row.names = 1)


# Verificar si los datos se han cargado correctamente
if (is.null(tumor_data)) {
  stop("Error al cargar el archivo de datos. Verifica la ruta y el formato del archivo.")
}

# Imputación de valores faltantes usando KNN

# Contar valores faltantes por columna
missing_per_column <- colSums(is.na(tumor_data))
print(missing_per_column)

# Contar valores faltantes por fila
missing_per_row <- rowSums(is.na(tumor_data))
print(missing_per_row)

# Crear un data frame para los valores faltantes por columna
missing_data_col <- data.frame(Column = names(missing_per_column), MissingValues = missing_per_column)

# Filtrar para mostrar solo columnas con valores faltantes
missing_data_col <- missing_data_col[missing_data_col$MissingValues > 0, ]

# Gráfico de barras
ggplot(missing_data_col, aes(x = reorder(Column, -MissingValues), y = MissingValues)) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  labs(title = "Número de Valores Faltantes por Columna", x = "Columnas", y = "Valores Faltantes")

# Crear un data frame para los valores faltantes por fila
missing_data_row <- data.frame(Row = seq_along(missing_per_row), MissingValues = missing_per_row)

# Filtrar para mostrar solo filas con valores faltantes
missing_data_row <- missing_data_row[missing_data_row$MissingValues > 0, ]

# Gráfico de barras
ggplot(missing_data_row, aes(x = Row, y = MissingValues)) +
  geom_bar(stat = "identity", fill = "red") +
  labs(title = "Número de Valores Faltantes por Fila", x = "Filas", y = "Valores Faltantes")

tumor_data_imputed <- impute.knn(as.matrix(tumor_data))$data

# Paso 5: Verificar la imputación. Comprobar si los valores faltantes han sido imputados
print("Datos después de la imputación:")
print(sum(is.na(tumor_data_imputed)))

# Paso 2: Transponer los datos para que las filas sean muestras y las columnas sean genes
expression_data <- t(tumor_data_imputed)



# ===================== # APRENDIZAJE NO SUPERVISADO =====================

# === Optimización del umbral basado en la modularidad ===
# Define los umbrales a probar y calcula modularidad y métricas adicionales

# Función para calcular la modularidad de una red en función del umbral
calculate_modularity <- function(expression_data, threshold) {
  # Crear la red aplicando el umbral
  adj_matrix <- cor(expression_data)
  adj_matrix[abs(adj_matrix) < threshold] <- 0
  graph_net <- graph.adjacency(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # Detectar comunidades usando el algoritmo de Walktrap
  communities <- cluster_walktrap(graph_net)
  
  # Calcular la modularidad
  modularity(communities)
}

# Generar umbrales a evaluar
thresholds <- seq(0.5, 0.95, by = 0.05)

# Inicializar una lista para almacenar los resultados
results <- data.frame(threshold = thresholds, modularity = NA, avg_degree = NA, num_edges = NA)

# Iterar sobre los umbrales para calcular la modularidad y métricas adicionales
for (i in 1:length(thresholds)) {
  threshold <- thresholds[i]
  
  # Calcular la matriz de correlación y aplicar el umbral usando valores absolutos
  adj_matrix <- abs(cor(expression_data))  # Tomamos el valor absoluto de la correlación
  adj_matrix[adj_matrix < threshold] <- 0  # Aplicamos el umbral para crear la red
  
  # Crear la red sin pesos negativos
  graph_net <- graph.adjacency(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # Calcular modularidad usando el algoritmo Walktrap
  communities <- cluster_walktrap(graph_net)
  modularity_value <- modularity(communities)
  
  # Calcular métricas adicionales: grado promedio y número de aristas
  avg_degree <- mean(degree(graph_net))
  num_edges <- ecount(graph_net)
  
  # Guardar resultados
  results$modularity[i] <- modularity_value
  results$avg_degree[i] <- avg_degree
  results$num_edges[i] <- num_edges
}

# Seleccionar el umbral óptimo basado en la modularidad máxima
optimal_threshold <- results$threshold[which.max(results$modularity)]
cat("El umbral óptimo basado en modularidad es:", optimal_threshold, "\n")

# Visualizar los resultados para explorar el comportamiento de la red con diferentes umbrales
plot(results$threshold, results$modularity, type = "b", col = "blue", xlab = "Umbral", ylab = "Modularidad", main = "Modularidad vs Umbral")
points(optimal_threshold, max(results$modularity), col = "red", pch = 19)  # Marcar el óptimo
plot(results$threshold, results$avg_degree, type = "b", col = "darkgreen", xlab = "Umbral", ylab = "Grado Promedio", main = "Grado Promedio vs Umbral")
plot(results$threshold, results$num_edges, type = "b", col = "purple", xlab = "Umbral", ylab = "Número de Aristas", main = "Número de Aristas vs Umbral")


# === Fin de Optimización de Umbral ===

# Crear la red final usando el umbral óptimo
# Definir los dos umbrales que quieres comparar
threshold_1 <- 0.8
threshold_2 <- 0.9

# Crear la red para el umbral de 0.8
adj_matrix_0_8 <- cor(expression_data)
adj_matrix_0_8[abs(adj_matrix_0_8) < threshold_1] <- 0
network_0_8 <- graph.adjacency(adj_matrix_0_8, mode = "undirected", weighted = TRUE, diag = FALSE)



# Crear la red para el umbral de 0.9
adj_matrix_0_9 <- cor(expression_data)
adj_matrix_0_9[abs(adj_matrix_0_9) < threshold_2] <- 0
network_0_9 <- graph.adjacency(adj_matrix_0_9, mode = "undirected", weighted = TRUE, diag = FALSE)

# Establecer el espacio de gráficos en dos columnas para ver ambas redes en la misma salida
par(mfrow = c(1, 2))

# Visualización de la red con umbral 0.8
plot(network_0_8, vertex.label = NA, vertex.size = log(degree(network_0_8) + 1) * 5,
     main = paste("Red con umbral 0.8"),cex.main = 1.5)

# Visualización de la red con umbral 0.9
plot(network_0_9, vertex.label = NA, vertex.size = log(degree(network_0_9) + 1) * 5,
     main = paste("Red con umbral 0.9"),cex.main = 1.5)

# Restablecer el espacio de gráficos a un solo panel
par(mfrow = c(1, 1))

# Función para calcular hubs
calculate_hubs <- function(expression_data, threshold) {
  pearson_corr <- cor(expression_data, method = "pearson")
  spearman_corr <- cor(expression_data, method = "spearman")
  kendall_corr <- cor(expression_data, method = "kendall")
  
  # Crear redes binarias para el umbral dado
  pearson_net <- (abs(pearson_corr) > threshold) * 1
  spearman_net <- (abs(spearman_corr) > threshold) * 1
  kendall_net <- (abs(kendall_corr) > threshold) * 1
  
  # Ensamblar la red final con umbrales flexibles
  ensemble_net <- (pearson_net + spearman_net + kendall_net) >= 2
  graph_net <- graph_from_adjacency_matrix(ensemble_net, mode = "undirected", diag = FALSE)
  
  # Calcular centralidad por grado (hubs)
  degree_values <- degree(graph_net)
  hubs <- names(degree_values[degree_values > quantile(degree_values, 0.9)])  # Genes en el top 10% por conexiones
  
  return(list(graph_net = graph_net, hubs = hubs, degree_values = degree_values))
}

# Generar las características para la red con el umbral 0.8
final_network_8 <- calculate_hubs(expression_data, 0.8) 

# Generar las características para la red con el umbral 0.9
final_network_9 <- calculate_hubs(expression_data, 0.9)


# Función para calcular NMI entre dos genes
calculate_nmi <- function(gene1_data, gene2_data) {
  joint_data <- data.frame(gene1_data, gene2_data)
  nmi_value <- mutinformation(discretize(joint_data[,1]), discretize(joint_data[,2]))
  return(nmi_value)
}


# Nombres de genes en el gráfico 0.9 y 0.8
graph_genes_9 <- V(final_network_9$graph_net)$name
print(graph_genes_9)

graph_genes_8 <- V(final_network_8$graph_net)$name
print(graph_genes_8)

# Nombres de columnas en expression_data
expression_genes <- colnames(expression_data)
print(expression_genes)

# Realizar dos comparaciones, con el 0.8 y 0.9

common_genes_8 <- intersect(graph_genes_8, expression_genes)
print(common_genes_8)
num_common_genes_8 <- length(common_genes_8)

common_genes_9 <- intersect(graph_genes_9, expression_genes)
print(common_genes_9)
num_common_genes_9 <- length(common_genes_9)

# Comparar número de aristas
num_edges_0_8 <- ecount(network_0_8)
num_edges_0_9 <- ecount(network_0_9)

# Comparar grados de hubs
hubs_0_8 <- final_network_8$hubs
hubs_0_9 <- final_network_9$hubs

# Imprimir resultados
cat("Número de aristas (Umbral 0.8):", num_edges_0_8, "\n")
cat("Número de aristas (Umbral 0.9):", num_edges_0_9, "\n")

# Comparar hubs
cat("Hubs (Umbral 0.8):", paste(hubs_0_8, collapse = ", "), "\n")
cat("Hubs (Umbral 0.9):", paste(hubs_0_9, collapse = ", "), "\n")


# ===================== # APRENDIZAJE SUPERVISADO =====================

# Modificar la función extract_features2 para incluir métricas de centralidad
extract_features2 <- function(graph_net, expression_data) {
  edges <- as.data.frame(as_edgelist(graph_net))
  colnames(edges) <- c("Gene1", "Gene2")
  
  # Filtrar solo las interacciones cuyos genes están presentes en expression_data
  edges <- edges[edges$Gene1 %in% colnames(expression_data) & edges$Gene2 %in% colnames(expression_data), ]
  
  # Calcular correlaciones y NMI para cada interacción
  edges$pearson_corr <- apply(edges, 1, function(x) cor(expression_data[, x[1]], expression_data[, x[2]], method = "pearson"))
  edges$spearman_corr <- apply(edges, 1, function(x) cor(expression_data[, x[1]], expression_data[, x[2]], method = "spearman"))
  edges$kendall_corr <- apply(edges, 1, function(x) cor(expression_data[, x[1]], expression_data[, x[2]], method = "kendall"))
  edges$nmi <- apply(edges, 1, function(x) calculate_nmi(expression_data[, x[1]], expression_data[, x[2]]))
  
  # Calcular grado y centralidad de los genes involucrados
  degree_values <- degree(graph_net)
  edges$degree_Gene1 <- degree_values[edges$Gene1]
  edges$degree_Gene2 <- degree_values[edges$Gene2]
  centrality_values <- betweenness(graph_net)
  edges$centrality_Gene1 <- centrality_values[edges$Gene1]
  edges$centrality_Gene2 <- centrality_values[edges$Gene2]
  
  return(edges)
}

# Generar las características para la red con el umbral óptimo
features_optim_9 <- extract_features2(final_network_9$graph_net, expression_data)
features_optim_8 <- extract_features2(final_network_8$graph_net, expression_data)

# Asegúrate de que features_optim contiene la columna pearson_corr
if ("pearson_corr" %in% colnames(features_optim_8)) {
  features_optim_8$label <- factor(ifelse(abs(features_optim_8$pearson_corr) > 0.8, 1, 0))
} else {
  stop("pearson_corr no encontrado en features_optim.")
}

# Asegúrate de que features_optim contiene la columna pearson_corr
if ("pearson_corr" %in% colnames(features_optim_9)) {
  features_optim_9$label <- factor(ifelse(abs(features_optim_9$pearson_corr) > 0.8, 1, 0))
} else {
  stop("pearson_corr no encontrado en features_optim.")
}

# Mostrar las primeras filas de las características generadas
print(head(features_optim_9))
print(head(features_optim_8))

# Resumen estadístico
summary(features_optim_9)
summary(features_optim_8)

# Histograma de las correlaciones para 0.9
hist(features_optim_9$pearson_corr, main = "Distribución de Correlaciones (Umbral 0.9)", xlab = "Correlaciones", col = "lightblue")

# Histograma de las correlaciones para 0.8
hist(features_optim_8$pearson_corr, main = "Distribución de Correlaciones (Umbral 0.8)", xlab = "Correlaciones", col = "lightgreen")

# Correlación entre grado y centralidad para el umbral 0.9
cor_9 <- cor(features_optim_9$degree_Gene1, features_optim_9$centrality_Gene1)
cat("Correlación entre grado y centralidad (Umbral 0.9):", cor_9, "\n")

# Correlación entre grado y centralidad para el umbral 0.8
cor_8 <- cor(features_optim_8$degree_Gene1, features_optim_8$centrality_Gene1)
cat("Correlación entre grado y centralidad (Umbral 0.8):", cor_8, "\n")

# Dibujar la red para el umbral 0.8

par(mfrow = c(1, 2))
plot(network_0_8, 
     vertex.size = log(degree(network_0_8) + 1) * 5,  # Tamaño de los vértices basado en el grado
     vertex.color = ifelse(V(network_0_8)$name %in% hubs_0_8, "red", "orange"),  # Hubs en rojo, otros en naranja
     main = "Hubs en la Red (Umbral 0.8)",
     vertex.label = NA  # No mostrar etiquetas
)

# Dibujar la red para el umbral 0.9
plot(network_0_9, 
     vertex.size = log(degree(network_0_9) + 1) * 5,  # Tamaño de los vértices basado en el grado
     vertex.color = ifelse(V(network_0_9)$name %in% hubs_0_9, "red", "orange"),  # Hubs en rojo, otros en naranja
     main = "Hubs en la Red (Umbral 0.9)",
     vertex.label = NA  # No mostrar etiquetas
)

#####ANALISIS DE ENRIQUECIMIENTO SOBRE LOS HUBS:
# Cargar las librerías necesarias
library(clusterProfiler)
library(org.Hs.eg.db)

# Convertir los nombres de los hubs a IDs de ENTREZ (si es necesario) para KEGG
# Esto asume que `hubs_0_8` y `hubs_0_9` contienen nombres en formato de símbolos de genes (SYMBOL)
hubs_0_8_entrez <- bitr(hubs_0_8, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
hubs_0_9_entrez <- bitr(hubs_0_9, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

# Realizar el análisis de enriquecimiento GO para los hubs de la red con umbral 0.8
enrich_go_0_8 <- enrichGO(gene = hubs_0_8, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                          ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)

# Realizar el análisis de enriquecimiento KEGG para los hubs de la red con umbral 0.8
enrich_kegg_0_8 <- enrichKEGG(gene = hubs_0_8_entrez, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

# Realizar el análisis de enriquecimiento GO para los hubs de la red con umbral 0.9
enrich_go_0_9 <- enrichGO(gene = hubs_0_9, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                          ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)

# Realizar el análisis de enriquecimiento KEGG para los hubs de la red con umbral 0.9
enrich_kegg_0_9 <- enrichKEGG(gene = hubs_0_9_entrez, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

# Visualizar los resultados de enriquecimiento
print("Enriquecimiento GO para hubs con umbral 0.8:")
print(enrich_go_0_8)
print("Enriquecimiento KEGG para hubs con umbral 0.8:")
print(enrich_kegg_0_8)

print("Enriquecimiento GO para hubs con umbral 0.9:")
print(enrich_go_0_9)
print("Enriquecimiento KEGG para hubs con umbral 0.9:")
print(enrich_kegg_0_9)

###Combino los genes y los mapeo###
# Extraer genes de los términos enriquecidos en GO para umbral 0.8
go_genes_0_8 <- unique(unlist(strsplit(enrich_go_0_8@result$geneID, "/")))
go_genes_0_8 <- intersect(go_genes_0_8, hubs_0_8)

# Extraer genes de los términos enriquecidos en KEGG para umbral 0.8
kegg_genes_0_8 <- unique(unlist(strsplit(enrich_kegg_0_8@result$geneID, "/")))
kegg_genes_0_8 <- intersect(kegg_genes_0_8, hubs_0_8)

# Combinar los genes de GO y KEGG, eliminando duplicados
combined_genes_0_8 <- union(go_genes_0_8, kegg_genes_0_8)

# Mostrar los genes combinados para el umbral 0.8
print("Genes combinados enriquecidos para umbral 0.8:")
print(combined_genes_0_8)

# Separar entre IDs numéricos (ENTREZ) y nombres de genes (símbolos)
entrez_ids_8 <- combined_genes_0_8[grepl("^[0-9]+$", combined_genes_0_8)]
symbol_ids_8 <- combined_genes_0_8[!grepl("^[0-9]+$", combined_genes_0_8)]

# Convertir ENTREZ IDs a símbolos de genes
entrez_to_symbol_8 <- bitr(entrez_ids_8, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

# Combinar todos los símbolos en un solo conjunto
combined_symbols_8 <- unique(c(symbol_ids_8, entrez_to_symbol_8$SYMBOL))

# Mostrar el conjunto homogéneo de símbolos de genes
print("Conjunto de genes combinado en formato de símbolos de genes:")
print(combined_symbols_8)


# Extraer genes de los términos enriquecidos en GO para umbral 0.9
go_genes_0_9 <- unique(unlist(strsplit(enrich_go_0_9@result$geneID, "/")))
go_genes_0_9 <- intersect(go_genes_0_9, hubs_0_9)

# Extraer genes de los términos enriquecidos en KEGG para umbral 0.9
kegg_genes_0_9 <- unique(unlist(strsplit(enrich_kegg_0_9@result$geneID, "/")))
kegg_genes_0_9 <- intersect(kegg_genes_0_9, hubs_0_9)

# Combinar los genes de GO y KEGG, eliminando duplicados
combined_genes_0_9 <- union(go_genes_0_9, kegg_genes_0_9)

# Mostrar los genes combinados para el umbral 0.9
print("Genes combinados enriquecidos para umbral 0.9:")
print(combined_genes_0_9)

# Separar entre IDs numéricos (ENTREZ) y nombres de genes (símbolos)
entrez_ids_9 <- combined_genes_0_9[grepl("^[0-9]+$", combined_genes_0_9)]
symbol_ids_9 <- combined_genes_0_9[!grepl("^[0-9]+$", combined_genes_0_9)]

# Convertir ENTREZ IDs a símbolos de genes
entrez_to_symbol_9 <- bitr(entrez_ids_9, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

# Combinar todos los símbolos en un solo conjunto
combined_symbols_9 <- unique(c(symbol_ids_9, entrez_to_symbol_8$SYMBOL))

# Mostrar el conjunto homogéneo de símbolos de genes
print("Conjunto de genes combinado en formato de símbolos de genes:")
print(combined_symbols_9)


# Crear y visualizar la subred de genes enriquecidos para el umbral 0.8
subgraph_enriched_0_8 <- induced_subgraph(network_0_8, 
                                          V(network_0_8)$name %in% combined_symbols_8)

# Visualizar la subred enriquecida en el umbral 0.8
plot(subgraph_enriched_0_8,
     vertex.color = ifelse(V(subgraph_enriched_0_8)$name %in% combined_symbols_8, "blue", "orange"),
     vertex.size = log(degree(subgraph_enriched_0_8) + 1) * 5,
     vertex.label = NA,
     main = "Subred de Genes Enriquecidos (Umbral 0.8)")

# Visualización de la red completa con umbral 0.8, resaltando hubs y hubs enriquecidos
plot(network_0_8, 
     vertex.color = ifelse(V(network_0_8)$name %in% combined_symbols_8 & V(network_0_8)$name %in% hubs_0_8, "blue", 
                           ifelse(V(network_0_8)$name %in% hubs_0_8, "red", "orange")),
     vertex.size = log(degree(network_0_8) + 1) * 5,  # Tamaño basado en el grado
     vertex.label = NA,  # Sin etiquetas para una visualización más limpia
     main = "Red Completa con Hubs (Umbral 0.8)")

# Crear y visualizar la subred de genes enriquecidos para el umbral 0.9
subgraph_enriched_0_9 <- induced_subgraph(network_0_9, 
                                          V(network_0_9)$name %in% combined_symbols_9)

# Visualizar la subred enriquecida en el umbral 0.9
plot(subgraph_enriched_0_9,
     vertex.color = ifelse(V(subgraph_enriched_0_9)$name %in% combined_symbols_9, "blue", "orange"),
     vertex.size = log(degree(subgraph_enriched_0_9) + 1) * 5,
     vertex.label = NA,
     main = "Subred de Genes Enriquecidos (Umbral 0.9)")

plot(network_0_9, 
     vertex.color = ifelse(V(network_0_9)$name %in% combined_symbols_9 & V(network_0_9)$name %in% hubs_0_9, "blue", 
                           ifelse(V(network_0_9)$name %in% hubs_0_9, "red", "orange")),
     vertex.size = log(degree(network_0_9) + 1) * 5,  # Tamaño basado en el grado
     vertex.label = NA,  # Sin etiquetas para una visualización más limpia
     main = "Red Completa con Hubs (Umbral 0.9)")



# Exportar la red completa de umbral 0.8
write_graph(network_0_8, file = "network_0_8.graphml", format = "graphml")

# Exportar la red completa de umbral 0.9
write_graph(network_0_9, file = "network_0_9.graphml", format = "graphml")

# Crear archivo de atributos para la red con umbral 0.8
node_attributes_0_8 <- data.frame(
  id = V(network_0_8)$name,
  color = ifelse(V(network_0_8)$name %in% combined_symbols_8 & V(network_0_8)$name %in% hubs_0_8, "blue",
                 ifelse(V(network_0_8)$name %in% hubs_0_8, "red", "orange"))
)

# Guardar como archivo CSV para Cytoscape
write.csv(node_attributes_0_8, file = "node_attributes_0_8.csv", row.names = FALSE)

# Repetir para la red con umbral 0.9
node_attributes_0_9 <- data.frame(
  id = V(network_0_9)$name,
  color = ifelse(V(network_0_9)$name %in% combined_symbols_9 & V(network_0_9)$name %in% hubs_0_9, "blue",
                 ifelse(V(network_0_9)$name %in% hubs_0_9, "red", "orange"))
)

write.csv(node_attributes_0_9, file = "node_attributes_0_9.csv", row.names = FALSE)

# Dotplot para GO (umbral 0.8)
dotplot(enrich_go_0_8, showCategory = 20, title = "Dotplot de GO para Hubs (Umbral 0.8)")

# Dotplot para KEGG (umbral 0.8)
dotplot(enrich_kegg_0_8, showCategory = 20, title = "Dotplot de KEGG para Hubs (Umbral 0.8)")

# Dotplot para GO (umbral 0.9)
dotplot(enrich_go_0_9, showCategory = 20, title = "Dotplot de GO para Hubs (Umbral 0.9)")

# Dotplot para KEGG (umbral 0.9)
dotplot(enrich_kegg_0_9, showCategory = 20, title = "Dotplot de KEGG para Hubs (Umbral 0.9)")




##########AQUÍ TERMINA######


# División en conjunto de entrenamiento (70%) y prueba (30%)
set.seed(123)
train_indices_8 <- sample(1:nrow(features_optim_8), 0.7 * nrow(features_optim_8))
train_data_8 <- features_optim_8[train_indices_8, ]
test_data_8 <- features_optim_8[-train_indices_8, ]

train_indices_9 <- sample(1:nrow(features_optim_9), 0.7 * nrow(features_optim_9))
train_data_9 <- features_optim_9[train_indices_9, ]
test_data_9 <- features_optim_9[-train_indices_9, ]

# Asegúrate de que no haya NA en train_data
if (any(is.na(train_data_8))) {
  train_data_8 <- na.omit(train_data_8)
}

if (any(is.na(train_data_9))) {
  train_data_9 <- na.omit(train_data_9)
}

cat("Número de muestras en el conjunto de entrenamiento:", nrow(train_data_8), "\n")
cat("Número de muestras en el conjunto de prueba:", nrow(test_data_8), "\n")

cat("Número de muestras en el conjunto de entrenamiento:", nrow(train_data_9), "\n")
cat("Número de muestras en el conjunto de prueba:", nrow(test_data_9), "\n")


# Random Forest

rf_model_8 <- randomForest(label ~ pearson_corr + spearman_corr + kendall_corr + nmi + degree_Gene1 + degree_Gene2 + centrality_Gene1 + centrality_Gene2, 
                           data = train_data_8, importance = TRUE)

rf_model_9 <- randomForest(label ~ pearson_corr + spearman_corr + kendall_corr + nmi + degree_Gene1 + degree_Gene2 + centrality_Gene1 + centrality_Gene2, 
                           data = train_data_9, importance = TRUE)

# Verificar la distribución de clases en el conjunto de entrenamiento 8 
cat("Distribución de clases en train_data_8:\n") 
print(table(train_data_8$label))
# Verificar la distribución de clases en el conjunto de prueba 8
cat("Distribución de clases en test_data_8:\n") 
print(table(test_data_8$label))

# Verificar la distribución de clases en el conjunto de entrenamiento 9 
cat("Distribución de clases en train_data_9:\n") 
print(table(train_data_9$label))
# Verificar la distribución de clases en el conjunto de prueba 9 
cat("Distribución de clases en test_data_9:\n") 
print(table(test_data_9$label))

# Imprimir resultados del modelo
print(rf_model_8)
rf_predictions_8 <- predict(rf_model_8, test_data_8)

print(rf_model_9)
rf_predictions_9 <- predict(rf_model_9, test_data_9)


# Matriz de confusión para Random Forest
confusion_matrix_rf_8 <- table(test_data_8$label, rf_predictions_8)
print("Confusion Matrix for Random Forest para umbral 0.8:")
print(confusion_matrix_rf_8)

confusion_matrix_rf_9 <- table(test_data_9$label, rf_predictions_9)
print("Confusion Matrix for Random Forest para umbral 0.9:")
print(confusion_matrix_rf_9)


library(caret)
confusionMatrix(rf_predictions_8, test_data_8$label)
confusionMatrix(rf_predictions_9, test_data_9$label)

# Cargar la librería ggplot2
library(ggplot2)

# Crear la matriz de confusión como un data frame
cm <- data.frame(
  Prediction = c("0", "0", "1", "1"),
  Actual = c("0", "1", "0", "1"),
  Freq = c(2, 0, 2, 17)
)

# Reorganizar el data frame para que ggplot lo reconozca correctamente
cm <- matrix(cm$Freq, nrow = 2)
rownames(cm) <- c("0", "1")
colnames(cm) <- c("0", "1")

# Convertir la matriz a un formato adecuado para ggplot
cm_df <- as.data.frame(as.table(cm))

# Crear el gráfico
ggplot(cm_df, aes(x = Var1, y = Var2, fill = Freq)) + 
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "blue") +
  geom_text(aes(label = Freq), vjust = 1) + 
  labs(x = "Predicted", y = "Actual") + 
  theme_minimal() +
  ggtitle("Confusion Matrix Heatmap")

# Visualización de la importancia de las variables en Random Forest

varImpPlot(rf_model_8)
varImpPlot(rf_model_9)

#curva ROC

# Asegúrate de que las predicciones sean numéricas
rf_predictions_numeric_8 <- as.numeric(rf_predictions_8) - 1  # Convertir 1 a 1 y 0 a 0 (restando 1)
rf_predictions_numeric_9 <- as.numeric(rf_predictions_9) - 1  # Convertir 1 a 1 y 0 a 0 (restando 1)

# Calcular la curva ROC
roc_obj_8 <- roc(test_data_8$label, rf_predictions_numeric_8)
roc_obj_9 <- roc(test_data_9$label, rf_predictions_numeric_9)

# Plotear la curva ROC
plot(roc_obj_8, main = "Curva ROC para Random Forest con umbral 0.8", col = "blue")
plot(roc_obj_9, main = "Curva ROC para Random Forest con umbral 0.9", col = "blue")

auc_value_8 <- auc(roc_obj_8)
cat("Área bajo la curva (AUC):", auc_value_8, "\n")
auc_value_9 <- auc(roc_obj_9)
cat("Área bajo la curva (AUC):", auc_value_9, "\n")

library(pROC)

# Asegúrate de que las predicciones sean numéricas
rf_predictions_numeric_8 <- as.numeric(rf_predictions_8) - 1
rf_predictions_numeric_9 <- as.numeric(rf_predictions_9) - 1

# Calcular la curva ROC
roc_obj_rf_8 <- roc(test_data_8$label, rf_predictions_numeric_8)
roc_obj_rf_9 <- roc(test_data_9$label, rf_predictions_numeric_9)

# Graficar la curva ROC
plot(roc_obj_rf_8, main = "Curva ROC para Random Forest con umbral 0.8", col = "blue")
abline(a=0, b=1, col="red", lty=2)  # Línea de referencia

plot(roc_obj_rf_9, main = "Curva ROC para Random Forest con umbral 0.9", col = "blue")
abline(a=0, b=1, col="red", lty=2)  # Línea de referencia



#### SVM####
svm_model_8 <- svm(label ~ pearson_corr + spearman_corr + kendall_corr + nmi + degree_Gene1 + degree_Gene2 + centrality_Gene1 + centrality_Gene2, 
                   data = train_data_8)
svm_model_9 <- svm(label ~ pearson_corr + spearman_corr + kendall_corr + nmi + degree_Gene1 + degree_Gene2 + centrality_Gene1 + centrality_Gene2, 
                   data = train_data_9)


svm_predictions_8 <- predict(svm_model_8, test_data_8)
svm_predictions_9 <- predict(svm_model_9, test_data_9)

# Matriz de confusión para SVM

confusion_matrix_svm_8 <- table(test_data_8$label, svm_predictions_8)
print("Confusion Matrix for SVM para umbral 0.8:")
print(confusion_matrix_svm_8)

confusion_matrix_svm_9 <- table(test_data_9$label, svm_predictions_9)
print("Confusion Matrix for SVM para umbral 0.9:")
print(confusion_matrix_svm_9)

# Cargar la librería ggplot2
library(ggplot2)

# Convertir la matriz de confusión a un data frame
cm_svm_8 <- as.data.frame(confusion_matrix_svm_8)

# Asegúrate de que las columnas se llaman correctamente
colnames(cm_svm_8) <- c("Actual", "Predicted", "Freq")

# Crear el gráfico
ggplot(cm_svm_8, aes(x = Predicted, y = Actual, fill = Freq)) + 
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "blue") +
  geom_text(aes(label = Freq), vjust = 1) + 
  labs(x = "Predicted", y = "Actual") + 
  theme_minimal() +
  ggtitle("Confusion Matrix Heatmap for SVM para el umbral 0.8")

# Convertir la matriz de confusión a un data frame
cm_svm_9 <- as.data.frame(confusion_matrix_svm_9)

# Asegúrate de que las columnas se llaman correctamente
colnames(cm_svm_9) <- c("Actual", "Predicted", "Freq")

# Crear el gráfico
ggplot(cm_svm_9, aes(x = Predicted, y = Actual, fill = Freq)) + 
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "blue") +
  geom_text(aes(label = Freq), vjust = 1) + 
  labs(x = "Predicted", y = "Actual") + 
  theme_minimal() +
  ggtitle("Confusion Matrix Heatmap for SVM para el umbral 0.8")



# Asegúrate de que las predicciones sean numéricas
svm_predictions_numeric_8 <- as.numeric(svm_predictions_8) - 1
svm_predictions_numeric_9 <- as.numeric(svm_predictions_9) - 1

# Calcular la curva ROC
roc_svm_8 <- roc(test_data_8$label, svm_predictions_numeric_8)
roc_svm_9 <- roc(test_data_9$label, svm_predictions_numeric_9)


# Cargar la librería pROC si aún no está cargada
library(pROC)

# Graficar la curva ROC
plot(roc_svm_8, main = "Curva ROC para SVM umbral 0.8", col = "green")
abline(a = 0, b = 1, col = "red", lty = 2)  # Línea de referencia

plot(roc_svm_9, main = "Curva ROC para SVM umbral 0.9", col = "green")
abline(a = 0, b = 1, col = "red", lty = 2)  # Línea de referencia

# Calcular el AUC
auc_svm_8 <- auc(roc_svm_8)
cat("Área bajo la curva (AUC) para SVM para el umbral 0.8:", auc_svm_8, "\n")
# Calcular el AUC
auc_svm_9 <- auc(roc_svm_9)
cat("Área bajo la curva (AUC) para SVM para el umbral 0.9:", auc_svm_9, "\n")


#Compración de curvas ROC
plot(roc_obj_rf_8, main = "Comparación de Curvas ROC umbral 0.8", col = "blue")
plot(roc_svm_8, add = TRUE, col = "green")
legend("bottomright", legend = c("RF (Azul)", "SVM (Verde)"), col = c("blue", "green"), lty = 1)

plot(roc_obj_rf_9, main = "Comparación de Curvas ROC umbral 0.9", col = "blue")
plot(roc_svm_9, add = TRUE, col = "green")
legend("bottomright", legend = c("Random Forest", "SVM"), col = c("blue", "green"), lty = 1)



####GENES HUBS####

# Filtrar los datos para incluir solo interacciones de genes hubs (umbral 0.8)
hubs_features_8 <- features_optim_8[
  features_optim_8$Gene1 %in% hubs_0_8 & features_optim_8$Gene2 %in% hubs_0_8, ]

# Filtrar los datos para incluir solo interacciones de genes hubs (umbral 0.9)
hubs_features_9 <- features_optim_9[
  features_optim_9$Gene1 %in% hubs_0_9 & features_optim_9$Gene2 %in% hubs_0_9, ]

# División en conjunto de entrenamiento (70%) y prueba (30%)
set.seed(123)
train_indices_8 <- sample(1:nrow(hubs_features_8), 0.7 * nrow(hubs_features_8))
train_data_8 <- hubs_features_8[train_indices_8, ]
test_data_8 <- hubs_features_8[-train_indices_8, ]

train_indices_9 <- sample(1:nrow(hubs_features_9), 0.7 * nrow(hubs_features_9))
train_data_9 <- hubs_features_9[train_indices_9, ]
test_data_9 <- hubs_features_9[-train_indices_9, ]

# Eliminar filas con valores NA
train_data_8 <- na.omit(train_data_8)
test_data_8 <- na.omit(test_data_8)
train_data_9 <- na.omit(train_data_9)
test_data_9 <- na.omit(test_data_9)

cat("Número de muestras en el conjunto de entrenamiento (umbral 0.8):", nrow(train_data_8), "\n")
cat("Número de muestras en el conjunto de prueba (umbral 0.8):", nrow(test_data_8), "\n")

cat("Número de muestras en el conjunto de entrenamiento (umbral 0.9):", nrow(train_data_9), "\n")
cat("Número de muestras en el conjunto de prueba (umbral 0.9):", nrow(test_data_9), "\n")

# Random Forest para genes hubs
library(randomForest)

rf_model_8 <- randomForest(label ~ pearson_corr + spearman_corr + kendall_corr + nmi +
                             degree_Gene1 + degree_Gene2 + centrality_Gene1 + centrality_Gene2,
                           data = train_data_8, importance = TRUE)

rf_model_9 <- randomForest(label ~ pearson_corr + spearman_corr + kendall_corr + nmi +
                             degree_Gene1 + degree_Gene2 + centrality_Gene1 + centrality_Gene2,
                           data = train_data_9, importance = TRUE)

# Imprimir resultados del modelo
print(rf_model_8)
print(rf_model_9)

# Predicciones y matrices de confusión
rf_predictions_8 <- predict(rf_model_8, test_data_8)
rf_predictions_9 <- predict(rf_model_9, test_data_9)

confusion_matrix_rf_8 <- table(test_data_8$label, rf_predictions_8)
confusion_matrix_rf_9 <- table(test_data_9$label, rf_predictions_9)

print("Matriz de Confusión para Random Forest (umbral 0.8):")
print(confusion_matrix_rf_8)

print("Matriz de Confusión para Random Forest (umbral 0.9):")
print(confusion_matrix_rf_9)

# Evaluar el rendimiento del modelo con métricas de caret
library(caret)
metrics_rf_8 <- confusionMatrix(rf_predictions_8, test_data_8$label)
metrics_rf_9 <- confusionMatrix(rf_predictions_9, test_data_9$label)

print(metrics_rf_8)
print(metrics_rf_9)

# Visualización de la importancia de las variables en Random Forest
varImpPlot(rf_model_8)
varImpPlot(rf_model_9)

# Curva ROC y AUC
library(pROC)

rf_predictions_numeric_8 <- as.numeric(rf_predictions_8) - 1
rf_predictions_numeric_9 <- as.numeric(rf_predictions_9) - 1

roc_obj_rf_8 <- roc(test_data_8$label, rf_predictions_numeric_8)
roc_obj_rf_9 <- roc(test_data_9$label, rf_predictions_numeric_9)

plot(roc_obj_rf_8, main = "Curva ROC para Random Forest (Genes Hubs, Umbral 0.8)", col = "blue")
abline(a = 0, b = 1, col = "red", lty = 2)

plot(roc_obj_rf_9, main = "Curva ROC para Random Forest (Genes Hubs, Umbral 0.9)", col = "blue")
abline(a = 0, b = 1, col = "red", lty = 2)

auc_value_8 <- auc(roc_obj_rf_8)
auc_value_9 <- auc(roc_obj_rf_9)

cat("Área bajo la curva (AUC) para umbral 0.8:", auc_value_8, "\n")
cat("Área bajo la curva (AUC) para umbral 0.9:", auc_value_9, "\n")

####SVM####

# SVM para genes hubs
library(e1071)

# Entrenar modelos SVM para genes hubs
svm_model_8 <- svm(label ~ pearson_corr + spearman_corr + kendall_corr + nmi + 
                     degree_Gene1 + degree_Gene2 + centrality_Gene1 + centrality_Gene2, 
                   data = train_data_8)
svm_model_9 <- svm(label ~ pearson_corr + spearman_corr + kendall_corr + nmi + 
                     degree_Gene1 + degree_Gene2 + centrality_Gene1 + centrality_Gene2, 
                   data = train_data_9)

# Realizar predicciones para los conjuntos de prueba
svm_predictions_8 <- predict(svm_model_8, test_data_8)
svm_predictions_9 <- predict(svm_model_9, test_data_9)

# Matrices de confusión para SVM
confusion_matrix_svm_8 <- table(test_data_8$label, svm_predictions_8)
print("Matriz de Confusión para SVM (genes hubs, umbral 0.8):")
print(confusion_matrix_svm_8)

confusion_matrix_svm_9 <- table(test_data_9$label, svm_predictions_9)
print("Matriz de Confusión para SVM (genes hubs, umbral 0.9):")
print(confusion_matrix_svm_9)

# Cargar la librería ggplot2 para visualizar los gráficos
library(ggplot2)

# Visualización de la matriz de confusión para umbral 0.8
cm_svm_8 <- as.data.frame(confusion_matrix_svm_8)
colnames(cm_svm_8) <- c("Actual", "Predicted", "Freq")

ggplot(cm_svm_8, aes(x = Predicted, y = Actual, fill = Freq)) + 
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "blue") +
  geom_text(aes(label = Freq), vjust = 1) + 
  labs(x = "Predicted", y = "Actual") + 
  theme_minimal() +
  ggtitle("Heatmap de Matriz de Confusión para SVM (genes hubs, umbral 0.8)")

# Visualización de la matriz de confusión para umbral 0.9
cm_svm_9 <- as.data.frame(confusion_matrix_svm_9)
colnames(cm_svm_9) <- c("Actual", "Predicted", "Freq")

ggplot(cm_svm_9, aes(x = Predicted, y = Actual, fill = Freq)) + 
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "blue") +
  geom_text(aes(label = Freq), vjust = 1) + 
  labs(x = "Predicted", y = "Actual") + 
  theme_minimal() +
  ggtitle("Heatmap de Matriz de Confusión para SVM (genes hubs, umbral 0.9)")

# Calcular la curva ROC para los modelos SVM
library(pROC)

# Convertir predicciones a numéricas para calcular la curva ROC
svm_predictions_numeric_8 <- as.numeric(svm_predictions_8) - 1
svm_predictions_numeric_9 <- as.numeric(svm_predictions_9) - 1

roc_svm_8 <- roc(test_data_8$label, svm_predictions_numeric_8)
roc_svm_9 <- roc(test_data_9$label, svm_predictions_numeric_9)

# Graficar la curva ROC para SVM (umbral 0.8)
plot(roc_svm_8, main = "Curva ROC para SVM (genes hubs, umbral 0.8)", col = "green")
abline(a = 0, b = 1, col = "red", lty = 2)

# Graficar la curva ROC para SVM (umbral 0.9)
plot(roc_svm_9, main = "Curva ROC para SVM (genes hubs, umbral 0.9)", col = "green")
abline(a = 0, b = 1, col = "red", lty = 2)

# Calcular el área bajo la curva (AUC) para ambos umbrales
auc_svm_8 <- auc(roc_svm_8)
auc_svm_9 <- auc(roc_svm_9)

cat("Área bajo la curva (AUC) para SVM (genes hubs, umbral 0.8):", auc_svm_8, "\n")
cat("Área bajo la curva (AUC) para SVM (genes hubs, umbral 0.9):", auc_svm_9, "\n")

# Comparación de curvas ROC entre Random Forest y SVM
plot(roc_obj_rf_8, main = "Comparación de Curvas ROC (genes hubs, umbral 0.8)", col = "blue")
plot(roc_svm_8, add = TRUE, col = "green")
legend("bottomright", legend = c("Random Forest (Azul)", "SVM (Verde)"), col = c("blue", "green"), lty = 1)

plot(roc_obj_rf_9, main = "Comparación de Curvas ROC (genes hubs, umbral 0.9)", col = "blue")
plot(roc_svm_9, add = TRUE, col = "green")
legend("bottomright", legend = c("Random Forest", "SVM"), col = c("blue", "green"), lty = 1)

