#######################################################################
# Archivo: enrichers_lmdme.R
# Descripción: Funciones de análisis funcional (GO/KEGG) sobre clusters
# Dependencias: utils.R
#######################################################################

source("R/utils.R")

#==============================================================
# Función: throw_enricher_lmdme
# Descripción: Lanza el análisis funcional basado en resultados 
# de clustering lmdme
#==============================================================

#' Lanza el análisis funcional basado en resultados de clustering lmdme
#'
#' Verifica que se ha realizado la identificación de DEG y clustering
#' por el método lmdme. Para ello comprueba que existen en el entorno global:
#' - 'table_cluster_lmdme': tabla de asignación de genes a clusters
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
#' 

throw_enricher_lmdme<-function(){
  # Comprueba que se ha realizado el análisis de DEG y que se han cargado
  # los archivos de anotaciones
  if (!("table_cluster_lmdme" %in% ls(envir = .GlobalEnv))) {
    cat("DEG no encontrados.Por favor, realize primero el análisis lmdme...\n")
  }else if (!("annondb" %in% ls(envir = .GlobalEnv))){ 
    cat("No se han cargado los ficheros para Analisis Funcional. Seleccione la opción 4...\n")
  }else { 
    
    cat("Buscando terminos GO por condicion significativa respecto al control...\n")
    enrich_conditions_GO_lmdme<-get_enrichment_by_significant_conditions(genes_list_lmdme,
                                                                         unique(annondb$go_db$Gene_ID),
                                                                         annondb$go_db %>% 
                                                                           dplyr::select(GO_ID, Gene_ID) %>% 
                                                                           drop_na(),
                                                                         annondb$go_db %>% 
                                                                           dplyr::select(GO_ID, TERM) %>% 
                                                                           drop_na(),
                                                                         "lmdme",
                                                                         0.05,
                                                                         "GO",
                                                                         "enrich_conditions")
    
    assign("enrich_conditions_GO_lmdme",enrich_conditions_GO_lmdme, envir = .GlobalEnv )
    
    cat("Buscando terminos KEGG por condicion significativa respecto al control...\n")
    enrich_conditions_KEGG_lmdme<-get_enrichment_by_significant_conditions(genes_list_lmdme,
                                                                           unique(annondb$kegg_db$Gene_ID),
                                                                           annondb$kegg_db %>% 
                                                                             dplyr::select(KEGG_ID, Gene_ID) %>% 
                                                                             drop_na(),
                                                                           annondb$kegg_db %>% 
                                                                             dplyr::select(KEGG_ID, NAME) %>% 
                                                                             drop_na(),
                                                                           "lmdme",
                                                                           0.05,
                                                                           "KEGG",
                                                                           "enrich_conditions")
    
    assign("enrich_conditions_KEGG_lmdme",enrich_conditions_KEGG_lmdme, envir = .GlobalEnv )

    #enrich_gsea_go_lmdme<-get_enrich_gsea(genes_model_lmdme,
    #                                      "lmdme",
    #                                      annondb$go_db %>% 
    #                                        dplyr::select(GO_ID, Gene_ID) %>% 
    #                                        drop_na(),
    #                                      annondb$go_db %>% 
    #                                        dplyr::select(GO_ID, TERM) %>% 
    #                                        drop_na(),
    #                                      minSize = 10, 
    #                                      maxSize = 500, 
    #                                      pCutoff = 0.05)
    #assign("enrich_gsea_go_lmdme",enrich_gsea_go_lmdme, envir = .GlobalEnv )
  }
}
