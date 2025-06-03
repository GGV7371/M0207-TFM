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
#' Si existen lanza la funcion:
#' - get_enrichment_by_significant_conditions() para encontrar términos GO/KEGG 
#'   en las comparaciones significativas en función de los parametros
#'   que se le pasan
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
