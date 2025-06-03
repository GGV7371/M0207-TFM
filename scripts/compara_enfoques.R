#######################################################################
# Archivo: compara_enfoques.R
# Descripción: Código de la opción 8 del menú
# Autora: Gemma Gariglio Viejo
# Fecha: 02-03-2025
# Dependencias: graphics.R
#######################################################################

source("R/graphics.R")

#==============================================================
# Función: throw_compara
# Descripción: Lanza el análisis DEG basado en el método lmdme
#==============================================================

#' Lanza todos los pasos para realizar la compración de los resultados 
#' entre lmdme y DESEq2
#' 
#' Consiste en realizar:
#' - Diagrama de venn de DEG 
#' 
#' @return No devuelve ningún valor. Una vez terminado el análisis
#' regresa a al proceso que lo llamó


throw_compara<-function(){
  
  cat("Generando diagrama de Venn de DEGs...\n")
  
  # Genera una lista de DEG con ambos enfoques
  venn_list <- list(lmdme=unique(table_cluster_lmdme$gene),
                    DESEq2=unique(table_cluster_deseq$gene))
  
  # LLamada a función que genera diagrama de Venn
  get_vennplot(venn_list)
  
  cat("Comparando Agrupamiento...\n")
  
  # Resumen de la distribución de genes por cluster para lmdme
  lmdme_counts <- table_cluster_lmdme %>%
    dplyr::select(gene, cluster) %>%
    distinct() %>%
    dplyr::count(cluster, name = "n_genes") %>%
    mutate(method = "LMDME",
           cluster = as.character(cluster),
           porc= round(100*n_genes/sum(n_genes),2))
    #rename(cluster_lmdme = cluster) %>%
    #mutate(por_lmdme = round(100*n_lmdme/sum(n_lmdme),2))
  
  # Resumen de la distribución de genes por cluster para deseq2
  deseq_counts <- table_cluster_deseq %>%
    dplyr::select(gene, cluster) %>%
    distinct() %>%
    dplyr::count(cluster, name = "n_genes") %>%
    mutate(method = "DESeq2",
           cluster = as.character(cluster),
           porc= round(100*n_genes/sum(n_genes),2))
    #rename(cluster_deseq = cluster) %>%
    #mutate(por_deseq = round(100*n_deseq/sum(n_deseq),2))
  
  # Union de tablas
  combined_counts <- bind_rows(lmdme_counts, deseq_counts)
  
  # Barplot de las distribuciones
  ggplot(combined_counts, aes(x = cluster, y = n_genes, fill = method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    geom_text(aes(label = paste0(porc, "%")), 
              position = position_dodge(width = 0.8),
              vjust = -0.5, size = 4) +
    labs(title = "Distribución de genes por clúster en ambos métodos",
         x = "Clúster",
         y = "Número de genes",
         fill = "Método") +
    theme_minimal(base_size = 14)
  
  # Guarda fichero con barplot 
  ggsave("results/images/barplot_comparacion_clustering.png", width = 8, height = 6)
  
  # Tabla con genes-cluster para lmdme
  genes_lmdme <- table_cluster_lmdme %>%
    distinct(gene, cluster) %>%
    rename(cluster_lmdme = cluster)
  
  # Tabla con genes-cluster para deseq2
  genes_deseq <- table_cluster_deseq %>%
    distinct(gene, cluster) %>%
    rename(cluster_deseq = cluster)
  
  # Unir tablas
  coincidencias <- inner_join(genes_lmdme, genes_deseq, by = "gene")
  
  heatmap_data <- coincidencias %>%
    dplyr::count(cluster_lmdme, cluster_deseq, name = "n_genes")
  
  # Crear heatmap
  ggplot(heatmap_data, aes(x = cluster_deseq, y = cluster_lmdme, fill = n_genes)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_genes), 
              color = "black", 
              size = 4) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(title = "Coincidencia de genes por clúster (LMDME vs DESeq2)",
         x = "Clúster DESeq2", 
         y = "Clúster LMDME", 
         fill = "Nº de genes") +
    theme_minimal(base_size = 14)
 
   # Guardar
  ggsave("results/images/heatmap_comparacion_clustering.png", width = 8, height = 6)
  
  cat("Comparando enriquecimiento funcional términos GO...\n")
  
  get_compara_enrich("GO")
  
  cat("Comparando enriquecimiento funcional términos KEGG...\n")
  
  get_compara_enrich("KEGG")
  
  
   return()
  
}
  