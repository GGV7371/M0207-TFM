#######################################################################
# Archivo: deseq2.R
# Descripción: Código de la opción 3 del menú
# Dependencias: utils.R, procedrures.R, graphics.R
#######################################################################

source("R/utils.R")
source("R/procedures.R")
source("R/graphics.R")

#==============================================================
# Función: throw_deseq2
# Descripción: Lanza el análisis DEG basado en el método DESEq2
#==============================================================

#' Lanza todos los pasos para realizar el análisis de expresión
#' diferencial basado en el método DESEq2+LRT
#' 
#' Consiste en realizar:
#' - Crear el objeto DESEqDataSet
#' - Filtrar genes con baja expresión
#' - Identificación de DEG
#' - Visualización de resultados: Genera gráficos y los guarda
#' - Clustering
#' 
#' @return No devuelve ningún valor. Una vez terminado el análisis
#' regresa a al proceso que lo llamó

throw_deseq2<-function(){
  if (!("counts" %in% ls(envir = .GlobalEnv)) || !("metadata" %in% ls(envir = .GlobalEnv))) {
    cat("❌ Datos no encontrados. Cargue los ficheros de trabajo...\n")
    return()
  }else{
  
    #==================================================
    # DESEQ2
    #==================================================
    cat("Creando el objeto DESEqDataSet...\n")
    
    conds<-c("treatment", "stage")
    design_formula <- get_formula_modelo(conds, "*")
    
    # Create DeseqDataSet Object
    ddsM <- DESeqDataSetFromMatrix(countData = counts,
                                   colData = metadata,
                                   design = design_formula)
    
    #--------------------------------------------------
    # Filtration
    #--------------------------------------------------
    #cat("Filtering out genes with low expression...\n")
    keep <- rowSums(counts(ddsM) >= 10) >= 3
    dds <- ddsM[keep, ]
    #dds<-ddsM
    print(nrow(dds))
    #--------------------------------------------------
    # Normalization VST for visualizations
    #--------------------------------------------------
    vsd <- vst(ddsM, blind = TRUE)
    vsdcounts <- assay(vsd)
    longvsddata <- get_logCountsLong(vsdcounts)
    get_box_recuentos(longvsddata, "vst")

    #--------------------------------------------------
    # DEG Analysis
    #--------------------------------------------------
    cat("Identificando DEG...\n")
    
    dds <- DESeq(dds)
    
    print(table(mcols(dds)$replace, useNA = "ifany"))
    
    desig_LRT<-get_formula_modelo(conds, "+")
    dds <- DESeq(dds, test = "LRT", reduced = desig_LRT)
    
    res_lrt <- results(dds)
    
    #assign("genes_model_deseq",res_lrt, envir = .GlobalEnv )
    table(is.na(res_lrt$padj))
    table(res_lrt$pvalue < 0.05)
    
    res_lrt <- res_lrt[!is.na(res_lrt$padj), ]
    res_lrt_df<-as.data.frame(res_lrt) %>%
      rownames_to_column(var="gene") %>%
      dplyr::select(gene, stat,pvalue,padj)
    
    assign("genes_model_deseq",res_lrt_df, envir = .GlobalEnv )
    
    sig_genes_pval_deseq <- rownames(res_lrt)[res_lrt$pvalue < 0.05]
    sig_genes_padj_deseq <- rownames(res_lrt)[res_lrt$padj < 0.05]
    
    cat("DEG by p-value LRT:", length(sig_genes_pval_deseq), "\n")
    cat("DEG by padj LRT:", length(sig_genes_padj_deseq), "\n")
    
    write.csv(as.data.frame(sig_genes_pval_deseq), 
              "results/datafiles/sigGenesPval_DESeq.csv")
    
    write.csv(as.data.frame(sig_genes_padj_deseq), 
              "results/datafiles/sigGenesPadj_DESeq.csv")
    
    #--------------------------------------------------
    # Results Visualization
    #--------------------------------------------------
    cat("Generando visualizaciones...\n")
    
    vsd<-vst(dds, blind = FALSE)
    matrix_deseq<-assay(vsd)
    
    scaled_matrix_deseq<-get_matrix_scaled(matrix_deseq, 
                                           sig_genes_padj_deseq)
    
    anno_col<-anota_col("stage", "treatment", NULL)
    anno_col_ord<-anota_col("stage", "treatment", "orden")
    
    scaled_matrix_deseq_ord<-get_order_samples(scaled_matrix_deseq, 
                                               anno_col)
    
    cat("Guardando Heatmap de DEGS...\n") 
    
    heat_sig_deseq<-get_heatmap_sig(scaled_matrix_deseq_ord,
                                    anno_col_ord,
                                    NULL,
                                    "DESEq")
    
    cat("Realizando clustering...\n")
    
    k<-get_kmedoid_optimal_clusters(scaled_matrix_deseq, 10, "silhouette", "DESeq2")
    
    genes_clusters_deseq<-get_clustering(scaled_matrix_deseq, k)
    
    anno_row<-anota_row(genes_clusters_deseq)
    
    cat("Guardando heatmap de los clusters...\n") 
    
    heat_sig_deseq<-get_heatmap_sig(scaled_matrix_deseq_ord,
                                    anno_col_ord,
                                    anno_row,
                                    "DESEq2")
    
    cat("Guardando gráfico PCA...\n") 
    
    get_pca(scaled_matrix_deseq, genes_clusters_deseq, "DESeq2")
    
    table_cluster_deseq <- get_table_cluster(scaled_matrix_deseq, 
                                             genes_clusters_deseq)
    
    cat("Guardando la tabla de DEGs agrupados por cluster...\n")
    
    write.csv(table_cluster_deseq, 
              "results/datafiles/table_gene_cluster_deseq.csv")
    assign("table_cluster_deseq",table_cluster_deseq, envir = .GlobalEnv )
    
    cat("Guardando tabla de expresion media by cluster...\n")
    
    cluster_avg_deseq <- table_cluster_deseq %>%
      group_by(cluster, trt_st, iden) %>%
      summarise(mean_expr = mean(expression), .groups = "drop")
    
    cluster_avg_deseq <- cluster_avg_deseq %>%
      left_join(table_cluster_deseq %>%
                  distinct(iden, trt_st, treatment, stage),
                by = c("iden", "trt_st"))
    
    assign("cluster_avg_deseq",cluster_avg_deseq, envir = .GlobalEnv )
    
    
    
    cat("Realizando comparaciones basadas en ANOVA-Tukey...\n")
    
    #tukey_deseq<-get_anova_tukey(cluster_avg_deseq, "DESEq2")
    tukey_deseq<-get_anova_tukey_by_state(cluster_avg_deseq, "DESEq2")
    assign("tukey_deseq",tukey_deseq, envir = .GlobalEnv )
    
    cat("Obteniendo DEG por condicion significativa respecto al control...\n")
    genes_list_deseq<-get_significant_genes_by_condition(table_cluster_deseq, 
                                                   tukey_deseq$tukey_table)
    assign("genes_list_deseq",genes_list_deseq, envir = .GlobalEnv )
    
    cat("Fin del análisis\n")
  }
    return()
}

