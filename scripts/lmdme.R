#######################################################################
# Archivo: lmdme.R
# Descripción: Código de la opción 2 del menú
# Autr: Gemma Gariglio Viejo
# Fecha: 
# Dependencias: utils.R, lmdmefun.R, procedrures.R, graphics.R
#######################################################################

source("R/lmdmefun.R")
source("R/procedures.R")
source("R/utils.R")
source("R/graphics.R")

#==============================================================
# Función: throw_lmdme
# Descripción: Lanza el análisis DEG basado en el método lmdme
#==============================================================

#' Lanza todos los pasos para realizar el análisis de expresión
#' diferencial basado en el modelo lmdme
#' 
#' Consiste en realizar:
#' - Transformación, Filtrado y Normalización de los datos de recuento
#' - Ajuste del modelo lmdme
#' - Identificación de DEG
#' - Visualización de resultados: Genera gráficos y los guarda
#' - Clustering
#' - Comparación ANOVA-Tukey
#' 
#' @return No devuelve ningún valor. Una vez terminado el análisis
#' regresa a al proceso que lo llamó

 
throw_lmdme<-function(){
  if (!("counts" %in% ls(envir = .GlobalEnv)) || !("metadata" %in% ls(envir = .GlobalEnv))) {
    cat("❌ Datos no encontrados. Cargue los ficheros de trabajo...\n")
  }else{  
    
    #==================================================
    # LOESS + lmdme
    #==================================================
    #--------------------------------------------------
    # Transformación y normalización
    #--------------------------------------------------
    cat("Transformación y normalización...\n")
    
    logcounts <- get_logCounts()
    logcountsloess <- normalizeCyclicLoess(logcounts)
    
    logCountsloessLong<-get_logCountsLong(logcountsloess)
    
    get_box_recuentos(logCountsloessLong, "loess")
    
    # Agrupa los transcritos si existen más de uno y filtra
    normdata<-join_transcripts(logcountsloess)
    
    # Normalización de muestras frente al control
    #matrix_lmdme <- get_centrado_ctr(normdata, 
    #                                 metadata,
    #                                 "iden",
    #                                 "treatment",
    #                                 "CTR",
    #                                 "stage")
    
    matrix_lmdme<-normdata
    
    #--------------------------------------------------
    # Ajuste del modelo
    #--------------------------------------------------
    cat("Ajustando el modelo...\n")
    design<-metadata %>% 
      dplyr::select(c(treatment, stage, iden)) 
    
    conds<-c("treatment", "stage")
    model <- get_formula_modelo(conds, "*")
    
    fit <- lmdme(model= model, 
                 data = matrix_lmdme, 
                 design = design)
    
    #--------------------------------------------------
    # Identificación de DEG
    #--------------------------------------------------
    cat("Identificando DEG...\n")
    
    # Define el termino a analizar
    term=paste0(conds[1],":",conds[2])
    
    # Extrae F.p.valores, p.values y coef del modelo
    fvals<-F.p.values(fit, term = term)
    
    id<-which(fvals< 0.05)
    
    fit_plsr<-fit
    
    # Realiza plsr sobre el modelo
    decomposition(fit_plsr,
                         decomposition = "plsr", 
                         type = "coefficient", 
                         #type = "residual",
                         term = term, 
                         subset = id, 
                         scale = "row")
    
    
    # Extrae p valores después de realizar plsr
    pvals<-p.values(fit_plsr,term = term)
    
    minPVal<-apply(pvals,1,min)
   
    
    minPVal_df <- data.frame(gene = names(minPVal),
                             pvalue = as.vector(minPVal))
    minPVal_df <- minPVal_df[!is.na(minPVal_df$pvalue), ]
   
    # obtiene un vector con los DEGs segun p.value <0.05
    sig_genes_pval_lmdme <- names(minPVal[minPVal < 0.05])
    
    # Calcula el p ajustado
    padj <- apply(pvals, 1, function(x) p.adjust(x, method = "BH"))
    padj <- t(padj)  
    minPadj <- apply(padj, 1, min)
    
    minPadj_df <- data.frame(gene = names(minPadj),
                             padj = as.vector(minPadj))
    minPadj_df <- minPadj_df[!is.na(minPadj_df$padj), ]
    
    genes_model_lmdme<-inner_join(minPVal_df, minPadj_df, by = "gene")
    assign("genes_model_lmdme",genes_model_lmdme, envir = .GlobalEnv )
    
    
    # obtiene un vector con los DEGs segun p.adj <0.05
    sig_genes_padj_lmdme <- names(minPadj[minPadj < 0.05])
    
    cat("DEG by pval LRT:", length(sig_genes_pval_lmdme), "\n")
    cat("DEG by padj por LRT:", length(sig_genes_padj_lmdme), "\n")
    
    cat("Guardando listado de DEGS...\n")
    # Guarda un fichero con el listado de DEGs con p.value<0.05
    write.csv(as.data.frame(sig_genes_pval_lmdme), 
              "results/datafiles/sigGenesPval_lmdme.csv")
    # Guarda un fichero con el listado de DEGs con p.adj<0.05
    write.csv(as.data.frame(sig_genes_padj_lmdme), 
              "results/datafiles/sigGenesPadj_lmdme.csv")
    write.csv(as.data.frame(genes_model_lmdme), 
              "results/datafiles/genes_model_lmdme.csv")
    
    #--------------------------------------------------
    # Visualización de resultados
    #--------------------------------------------------
    cat("Generando visualizaciones...\n")
    
    # Filtrado y escalado de la matriz de expresión para los DEGs
    scaled_matrix_lmdme<-get_matrix_scaled(matrix_lmdme, 
                                           sig_genes_padj_lmdme)
    
    # crea un variable con las condiciones
    anno_col<-anota_col("stage", "treatment", NULL)
    # ordena las muestras según las condiciones
    anno_col_ord<-anota_col("stage", "treatment", "orden")
    
    # ordena la matriz escalada según las condiciones
    scaled_matrix_lmdme_ord<-get_order_samples(scaled_matrix_lmdme, 
                                               anno_col)
    cat("Guardando Heatmap de DEGS...\n") 
    
    heat_sig_lmdme<-get_heatmap_sig(scaled_matrix_lmdme_ord,
                                    anno_col_ord,
                                    NULL,
                                    "lmdme")
    
    # Clustering
    
    cat("Realizando clustering...\n")
    
    # Calcula el número óptimo de clusteres
    k<-get_kmedoid_optimal_clusters(scaled_matrix_lmdme, 10, "silhouette", "lmdme")
    
    # Asigna a cada gen su cluster
    genes_clusters_lmdme<-get_clustering(scaled_matrix_lmdme, k)
    
    anno_row<-anota_row(genes_clusters_lmdme)
    
    # LLama a la función que genera el heatmap
    heat_sig_lmdme<-get_heatmap_sig(scaled_matrix_lmdme_ord,
                                    anno_col_ord,
                                    anno_row,
                                    "lmdme")

    # Genera un gráfico PCA
    cat("Guardando gráfico PCA...\n") 
    
    get_pca(scaled_matrix_lmdme, genes_clusters_lmdme, "lmdme")
    
    # tabla expresión de DEGs agrupados por cluster
    table_cluster_lmdme <- get_table_cluster(scaled_matrix_lmdme, 
                                             genes_clusters_lmdme)
    
    cat("Guardando la tabla de DEGs agrupados por cluster...\n")
    
    write.csv(table_cluster_lmdme, 
              "results/datafiles/table_gene_cluster_lmdme.csv")
    
    # Genera un heatmap de los cluster
    
    cat("Guardando heatmap de los clusters...\n") 
    
    # Genera un variable global
    assign("table_cluster_lmdme",table_cluster_lmdme, envir = .GlobalEnv )
    
    cat("Guardando tabla de expresion media by cluster...\n")
    
    cluster_avg_lmdme <- table_cluster_lmdme %>%
      group_by(cluster, trt_st, iden) %>%
      summarise(mean_expr = mean(expression), .groups = "drop")
    
    cluster_avg_lmdme <- cluster_avg_lmdme %>%
      left_join(table_cluster_lmdme %>%
                  distinct(iden, trt_st, treatment, stage),
                by = c("iden", "trt_st"))
    
    assign("cluster_avg_lmdme",cluster_avg_lmdme, envir = .GlobalEnv )
    
    cat("Realizando comparaciones basadas en ANOVA-Tukey...\n")
    
    #tukey_lmdme<-get_anova_tukey(cluster_avg_lmdme, "lmdme")
    tukey_lmdme<-get_anova_tukey_by_state(cluster_avg_lmdme, "lmdme")
    assign("tukey_lmdme",tukey_lmdme, envir = .GlobalEnv )
    
    cat("Obteniendo DEG por condicion significativa respecto al control...\n")
    genes_list_lmdme<-get_significant_genes_by_condition(table_cluster_lmdme, 
                                                         tukey_lmdme$tukey_table)
    assign("genes_list_lmdme",genes_list_lmdme, envir = .GlobalEnv )
    
    cat("Fin del análisis\n")
  }
  return()
}


