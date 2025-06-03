#######################################################################
# Archivo: enrichers_deseq2.R
# Descripción: Funciones de análisis funcional (GO/KEGG) sobre clusters
# Dependencias: utils.R
#######################################################################

source("R/utils.R")

#==============================================================
# Función: throw_enricher_deseq2
# Descripción: Lanza el análisis funcional basado en resultados 
# de clustering DESEq2
#==============================================================

#' Lanza el análisis funcional basado en resultados de clustering deseq2
#'
#' Verifica que se ha realizado la identificación de DEG y clustering
#' por el método deseq2. Para ello comprueba que existen en el entorno global:
#' - 'table_cluster_deseq': tabla de asignación de genes a clusters
#' - 'annondb': base de anotaciones funcionales (GO, KEGG)
#'
#' Si existen lanza las funciones:
#' - get_compare_cluster_go() para comparar términos GO entre clusters
#' - get_compare_cluster_kegg() para comparar rutas metabólicas entre clusters
#' - get_enricher_tukey_go() para encontrar términos GO en las comparaciones 
#'   significativas
#' - get_enricher_tukey_kegg() para encontrar rutas metabólicas en las 
#'   comparaciones significativas
#' @return No devuelve ningún valor. Imprime mensajes de progreso y ejecuta 
#' funciones de análisis.

throw_enricher_deseq2<-function(){
  # Comprueba que se ha realizado el análisis de DEG y que se han cargado
  # los archivos de anotaciones
  if (!("table_cluster_deseq" %in% ls(envir = .GlobalEnv))) {
    cat("DEG no encontrados.Por favor, realize primero el análisis DESEq2...\n")
  }else if (!("annondb" %in% ls(envir = .GlobalEnv))){ 
    cat("No se han cargado los ficheros para Analisis Funcional. Seleccione la opción 4...\n")
  }else {
    
    #cat("Comparando terminos GO entre clusters...\n")
    
    # Compara términos GO entre clusters
    #get_compare_cluster_go(table_cluster_deseq,"DESeq2")
    
    #cat("Comparando términos KEGG entre clusters...\n")
    
    # Compara términos KEGG entre clusters
    #get_compare_cluster_kegg(table_cluster_deseq,"DESeq2")
    
    #cat("Buscando términos GO comparando con el control...\n")
    
    # Busca términos GO para las comparaciones significativas
    #get_enricher_tukey_go(table_cluster_deseq, 
    #                      tukey_deseq$tukey_table, 
    #                      "CTR:embryo", 
    #                      0.05, 
    #                      "DESeq2")
    #get_enricher_tukey_go(table_cluster_deseq, 
    #                      tukey_deseq$tukey_table, 
    #                      "CTR:unprovisioned_egg", 
    #                      0.05, 
    #                      "DESeq2")
    #get_enricher_tukey_go(table_cluster_deseq, 
    #                      tukey_deseq$tukey_table, 
    #                      "CTR:provisioned_egg", 
    #                      0.05, 
    #                      "DESeq2")
    
    #cat("Buscando terminos KEGG comparando con el control...\n")
    
    # Busca rutas metabólicas para las comparaciones significativas
   
    #get_enricher_tukey_kegg(table_cluster_deseq, 
    #                        tukey_deseq$tukey_table, 
    #                        "CTR:embryo", 
    #                        0.05, 
    #                        "DESeq2")
    #get_enricher_tukey_kegg(table_cluster_deseq, 
    #                        tukey_deseq$tukey_table, 
    #                        "CTR:unprovisioned_egg", 
    #                        0.05, 
    #                        "DESeq2")
    #get_enricher_tukey_kegg(table_cluster_deseq, 
    #                        tukey_deseq$tukey_table, 
    #                        "CTR:provisioned_egg", 
    #                        0.05, 
    #                        "DESeq2")
    
    #cat("Buscando terminos KEGG para la condicion dominante en cada ...\n")
    
    #get_kegg_by_dominant_condition(table_cluster_deseq, 
    #                               tukey_deseq$dom_table, 
    #                               "DESeq2",
    #                               "KEGG_dominant")
    
    cat("Buscando terminos GO por condicion significativa respecto al control...\n")
    enrich_conditions_GO_deseq<-get_enrichment_by_significant_conditions(genes_list_deseq,
                                                                         unique(annondb$go_db$Gene_ID),
                                                                         annondb$go_db %>% 
                                                                           dplyr::select(GO_ID, Gene_ID) %>% 
                                                                           drop_na(),
                                                                         annondb$go_db %>% 
                                                                           dplyr::select(GO_ID, TERM) %>% 
                                                                           drop_na(),
                                                                         "DESEq2",
                                                                         0.05,
                                                                         "GO",
                                                                         "enrich_conditions")
    
    assign("enrich_conditions_GO_deseq",enrich_conditions_GO_deseq, envir = .GlobalEnv )
    
    cat("Buscando terminos KEGG por condicion significativa respecto al control...\n")
    enrich_conditions_KEGG_deseq<-get_enrichment_by_significant_conditions(genes_list_deseq,
                                                                           unique(annondb$kegg_db$Gene_ID),
                                                                           annondb$kegg_db %>% 
                                                                             dplyr::select(KEGG_ID, Gene_ID) %>% 
                                                                             drop_na(),
                                                                           annondb$kegg_db %>% 
                                                                             dplyr::select(KEGG_ID, NAME) %>% 
                                                                             drop_na(),
                                                                           "DESEq2",
                                                                           0.05,
                                                                           "KEGG",
                                                                           "enrich_conditions")
    
    assign("enrich_conditions_KEGG_deseq",enrich_conditions_KEGG_deseq, envir = .GlobalEnv )
    
    #enrich_gsea_go_deseq<-get_enrich_gsea(genes_model_deseq,
    #                                   "deseq",
    #                                   annondb$go_db %>% 
    #                                     dplyr::select(GO_ID, Gene_ID) %>% 
    #                                     drop_na(),
    #                                   annondb$go_db %>% 
    #                                     dplyr::select(GO_ID, TERM) %>% 
    #                                     drop_na(),
    #                                   minSize = 10, 
    #                                   maxSize = 500, 
    #                                   pCutoff = 0.05)
    #assign("enrich_gsea_go_deseq",enrich_gsea_go_deseq, envir = .GlobalEnv )
  }
}
