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
    
    #cat("Comparando términos GO entre clústeres...\n")
    
    # Comparar términos GO entre clústeres
    #get_compare_cluster_enrichment(table_cluster_lmdme,
    #                               term2gene = annondb$go_db %>% 
    #                                 dplyr::select(GO_ID, Gene_ID) %>% 
    #                                 drop_na(),
    #                               term2name = annondb$go_db %>% 
    #                                 dplyr::select(GO_ID, TERM) %>% 
    #                                 drop_na(),
    #                               "GO",
    #                               "lmdme",
    #                               0.05,
    #                              "compare_cluster_enrichment")
    
    #cat("Comparando términos GO entre condiciones...\n")
    #get_compare_condition_enrichment(table_cluster_lmdme,
    #                                 term2gene = annondb$go_db %>% 
    #                                   dplyr::select(GO_ID, Gene_ID) %>%
    #                                   drop_na(),
    #                                 term2name = annondb$go_db %>% 
    #                                   dplyr::select(GO_ID, TERM) %>%
    #                                   drop_na(),
    #                                 "GO",
    #                                 "lmdme",
    #                                 0.05,
    #                                 "compare_condition_enrichment")
    
    #cat("Comparando terminos KEGG entre clusters...\n")
    
    # Comparar términos KEGG entre clústeres
    #get_compare_cluster_kegg(table_cluster_lmdme,"lmdme")
    #get_compare_cluster_enrichment(table_cluster_lmdme,
    #                               term2gene = annondb$kegg_db %>% 
    #                                 dplyr::select(KEGG_ID, Gene_ID) %>% 
    #                                 drop_na(),
    #                               term2name = annondb$kegg_db %>% 
    #                                 dplyr::select(KEGG_ID, NAME) %>% 
    #                                 drop_na(),
    #                               "KEGG",
    #                               "lmdme",
    #                               0.05,
    #                               "compare_cluster_enrichment")
    
    #cat("Comparando terminos KEGG entre condiciones...\n")
    
    #get_compare_condition_enrichment(table_cluster_lmdme,
    #                                 term2gene = annondb$kegg_db %>% 
    #                                 dplyr::select(KEGG_ID, Gene_ID) %>% 
    #                                   drop_na(),
    #                                 term2name = annondb$kegg_db %>% 
    #                                   dplyr::select(KEGG_ID, NAME) %>% 
    #                                   drop_na(),
    #                                 "KEGG",
    #                                 "lmdme",
    #                                 0.05,
    #                                 "compare_condition_enrichment")
    
    #cat("Buscando términos GO para las comparaciones significativas...\n")
    
    # Busca términos GO para las comparaciones significativas
    #get_enricher_tukey_go(table_cluster_lmdme, 
    #                      tukey_lmdme$tukey_table, 
    #                      "CTR:embryo", 
    #                      0.05, 
    #                      "lmdme")
    
    #get_enricher_tukey_go(table_cluster_lmdme, 
    #                      tukey_lmdme$tukey_table, 
    #                      "CTR:unprovisioned_egg", 
    #                      0.05, 
    #                      "lmdme")
    
    #get_enricher_tukey_go(table_cluster_lmdme, 
    #                      tukey_lmdme$tukey_table, 
    #                      "CTR:provisioned_egg", 
    #                      0.05, 
    #                      "lmdme")
    
    #cat("Buscando rutas metabólicas para las comparaciones significativas...\n")
    
    # Busca rutas metabólicas para las comparaciones significativas
    #get_enricher_tukey_kegg(table_cluster_lmdme, 
    #                        tukey_lmdme$tukey_table, 
    #                        "CTR:embryo", 
    #                        0.05, 
    #                       "lmdme")
    
    #get_enricher_tukey_kegg(table_cluster_lmdme, 
    #                        tukey_lmdme$tukey_table, 
    #                        "CTR:unprovisioned_egg", 
    #                        0.05, 
    #                        "lmdme")
    
    #get_enricher_tukey_kegg(table_cluster_lmdme, 
    #                        tukey_lmdme$tukey_table, 
    #                        "CTR:provisioned_egg", 
    #                        0.05, 
    #                        "lmdme")
    
    #cat("Buscando rutas metabólicas para las condiciones dominantes...\n")
    
    #get_enrich_by_dominant_condition(table_cluster_lmdme, 
    #                                 tukey_lmdme, 
    #                                 unique(annondb$go_db$Gene_ID),
    #                                 annondb$go_db %>% 
    #                                   dplyr::select(GO_ID, Gene_ID) %>% 
    #                                   drop_na(),
    #                                 annondb$go_db %>% 
    #                                   dplyr::select(GO_ID, TERM) %>% 
    #                                   drop_na(),
    #                                 "lmdme",
    #                                 0.05,
    #                                 "GO_dominant")
    
    #get_enrich_by_dominant_condition(table_cluster_lmdme, 
    #                                 tukey_lmdme,
    #                                 unique(annondb$kegg_db$Gene_ID),
    #                                 annondb$kegg_db %>% 
    #                                   dplyr::select(KEGG_ID, Gene_ID) %>% 
    #                                   drop_na(),
    #                                 annondb$kegg_db %>% 
    #                                   dplyr::select(KEGG_ID, NAME) %>% 
    #                                   drop_na(),
    #                                 "lmdme",
    #                                 0.05,
    #                                 "KEGG_dominant")
    
    #for (s in unique(metadata$stage)){
    #  print(s)
    #  get_enrichment_by_significant_conditions(table_cluster_lmdme, 
    #                                           tukey_lmdme$tukey_table,
    #                                           paste0("CTR:",s),
    #                                           unique(annondb$go_db$Gene_ID),
    #                                           annondb$go_db %>% 
    #                                             dplyr::select(GO_ID, Gene_ID) %>% 
    #                                             drop_na(),
    #                                           annondb$go_db %>% 
    #                                             dplyr::select(GO_ID, TERM) %>% 
    #                                             drop_na(),
    #                                           "lmdme",
    #                                           0.05,
    #                                           "GO",
    #                                           "enrich_conditions")
    #}
    
    
    #cat("Obteniendo DEG por condicion significativa respecto al control...\n")
    #genes_list<-get_significant_genes_by_condition(table_cluster_lmdme, 
    #                                               tukey_lmdme$tukey_table)
    
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

    #for (s in unique(metadata$stage)){
    #  get_enrichment_by_significant_conditions(table_cluster_lmdme, 
    #                                           tukey_lmdme$tukey_table,
    #                                           paste0("CTR:",s),
    #                                           unique(annondb$kegg_db$Gene_ID),
    #                                           annondb$kegg_db %>% 
    #                                             dplyr::select(KEGG_ID, Gene_ID) %>% 
    #                                             drop_na(),
    #                                           annondb$kegg_db %>% 
    #                                             dplyr::select(KEGG_ID, NAME) %>% 
    #                                             drop_na(),
    #                                           "lmdme",
    #                                           0.05,
    #                                           "KEGG",
    #                                           "enrich_conditions")
    #}
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
