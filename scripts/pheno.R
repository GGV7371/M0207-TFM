#######################################################################
# Archivo: pheno.R
# Descripción: Código de la opción 9 del menú 
# Autr: Gemma Gariglio Viejo
# Fecha: 02/06/2025
# Dependencias: utils.R, procedrures.R, lmdme.R, deseq2.R
#######################################################################

source("scripts/lmdme.R")
source("scripts/deseq2.R")
source("R/utils.R")
source("R/procedures.R")

#==============================================================
# Función: throw_pheno
# Descripción: Carga los ficheros con los fenotipos
# y llama a las funciones para generar la correlación entre los 
# terminos GO/KEGG y los fenotipos
#==============================================================

#' Carga el fichero con los metadatos de fenotipos
#' y llama a los funciones que generan los heatmap de correlación
#' 
#' @return No devuelve ningún valor. Una vez terminado el análisis
#' regresa a al proceso que lo llamó

throw_pheno<-function(){
  # Lee el fichero con la matriz de recuentos
  if (!exists("pheno", envir = .GlobalEnv)) {
    cat("Introduce la ubicación del fichero con la información de fenotipos\n")
    file_counts<-readline(prompt = "Metadatos fenotípicos (.csv, .xlsx):" )
    pheno<-load_data(file_counts)
    assign("pheno",pheno, envir = .GlobalEnv )
  }
  
  cat("Analisis fenotipico para LMDME\n")
  
  cat("Generando heatmap de fenotipos frente a cluster...\n")
  heatmap_pheno_cluster(cluster_avg_lmdme, "lmdme")
  
  cat("Generando heatmap de fenotipos frente a Terminos GO/KEGG...\n")
  heatmap_pheno_enrich(enrich_conditions_GO_lmdme, "GO", "lmdme")
  heatmap_pheno_enrich(enrich_conditions_KEGG_lmdme, "KEGG", "lmdme")
  heatmap_pheno_enrich_cluster(enrich_conditions_GO_lmdme, "GO", "lmdme")
  
  cat("Analisis fenotipico para DESeq2\n")
  
  cat("Generando heatmap de fenotipos frente a cluster...\n")
  heatmap_pheno_cluster(cluster_avg_deseq, "DESeq2")
  
  cat("Generando heatmap de fenotipos frente a Terminos GO/KEGG...\n")
  heatmap_pheno_enrich(enrich_conditions_GO_deseq, "GO", "DESeq2")
  heatmap_pheno_enrich(enrich_conditions_KEGG_deseq, "KEGG", "DESeq2")
  heatmap_pheno_enrich_cluster(enrich_conditions_KEGG_deseq, "GO", "DESeq2")
  
}