############################################################
# Archivo: graphics.R
# Descripción: Funciones para grear las visualizaciones
# Autor: Gemma Gariglio Viejo
# Fecha: [Fecha]
############################################################

source("R/utils.R")
source("R/procedures.R")

#==============================================================
# Función: get_dim_biblio
# Descripción:  Genera gráficos de barras del total de recuentos
# por muestra y del numero de genes por muestra
#==============================================================

#' Genera dos barplot: 
#' - Total de recuentos por muestra (tamaño de la biblioteca) 
#' - Numero de genes detectados por muestra 
#' y los guarda en un fichero en el directorio de resultados
#' @return No devuelve nada

get_dim_biblio<-function(){
  
  #Suma los recuentos en cada muestra
  r<-colSums(counts)
  
  # Crea un data frame de los recuentos y condiciones experimentales
  r_d <- data.frame(
    sample = names(r),
    totalReads = r,
    treatment = metadata$treatment,
    stage = metadata$stage
  )
  
  # Genera un barplot de los recuentos por muestra
  bplotr<-ggplot(r_d, aes(x = sample, y = r)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "Total de Recuentos por Muestra", 
         x = "Muestra", 
         y = "Total de recuentos") +
    theme_minimal() +
    theme(axis.text.x = element_text(size=6, angle = 90, hjust = 1))
  
  # suprime los mensasjes que R muestra por consola
  suppressMessages(ggsave(paste0("results/images/barplot_libreria.png"), plot=bplotr))
  
  dg <- colSums(counts > 0)
  
  dg_d <- data.frame(
    sample = names(dg),
    detected_genes = dg
  )
  
  bplordg<-ggplot(dg_d, aes(x = sample, y = dg)) +
    geom_bar(stat = "identity", fill = "tomato") +
    labs(title = "Número de Genes Detectados por Muestra", x = "Muestra", y = "Genes detectados (>0)") +
    theme_minimal() +
    theme(axis.text.x = element_text(size=6, angle = 90, hjust = 1))
  
  suppressMessages(ggsave(paste0("results/images/barplot_genes_muestra.png"), plot=bplordg))
}

#==============================================================
# Función: get_boxplot_recuentos
# Descripción:  Genera gráficas con el tamaño de la biblioteca
# y con el numero de recuentos por gen
#==============================================================

#' Genera dos gráficos: 
#' - barplot con la distribución de recuentos transformados 
#' - density plot de los recuentos recuentos transformados 
#' y los guarda en un fichero en el directorio de resultados
#' @return No devuelve nada. Guarda un fichero con el gráfico

get_box_recuentos<-function(lc_d, t, cond){

  #l_c<-get_logCounts("log")
  
  #lc_d<-get_logCountsLong(lc_d)
  library(scales)
  
  num_grupos <- length(unique(lc_d$trt_st))
  colores <- hue_pal()(num_grupos)

  bxplot<-ggplot(lc_d, aes(x = iden, 
                           y = logcount,
                           fill = trt_st)) +
    geom_boxplot(width = 0.8,
                 position = position_dodge(0.9)) +
    #geom_boxplot(fill = "lightblue") +
    scale_fill_manual(values = colores) +
    scale_x_discrete(limits = unique(lc_d$iden)) +
    labs(title = "A. Distribución de Recuentos", 
         x = "Muestra", 
         y = t) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, 
                                     hjust = 0.5, 
                                     size = rel(0.5)),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 10),
          legend.position = "bottom")
   
  
  suppressMessages(ggsave(filename=paste0("results/images/boxplot_recuentos_",
                                          t,".png"), 
                          plot=bxplot))
  
  dnplot<-ggplot(lc_d, aes(x = logcount, 
                           group = iden, 
                           fill=trt_st)) +
    geom_density(alpha = 0.2, linewidth = 0.25) +
    scale_fill_manual(values = colores) +
    facet_wrap(~trt_st, dir = 'v') +
    labs(title = "B. Densidad de Recuentos", 
         x = t, 
         y = "Densidad") +
    theme_minimal() +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8))
  
  suppressMessages(ggsave(filename=paste0("results/images/density_recuentos_",
                                          t,".png"), 
                          plot=dnplot))
}

#==============================================================
# Función: get_pca_brutos
# Descripción:  Genera un gráfico PCA de los recuentos brutos
#==============================================================

#' Genera un grafico PCA de los recuentos transformados 
#' y los guarda en un fichero en el directorio de resultados
#' @return No devuelve nada

get_pca_brutos<-function(){
  library(FactoMineR)
  library(factoextra)
  
  l_c<-get_logCounts()
  lc_d<-get_logCountsLong(l_c)
  
  pca_res <- prcomp(t(l_c), scale. = TRUE)
  
  # Crear un data frame para ggplot
  pca_df <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    sample = colnames(pca_res$x),
    treatment = metadata$treatment,
    stage = metadata$stage
  )
  
  pcaplot<-ggplot(pca_df, aes(x = PC1, y = PC2, color = treatment, shape = stage)) +
    geom_point(size = 3) +
    labs(title = "PCA de Muestras RNA-seq", x = "PC1", y = "PC2") +
    theme_minimal()
  invisible(ggsave(filename=paste0("results/images/pca_recuentos.png"), plot=pcaplot))
}


#==============================================================
# Función: get_heatmap_sig
# Descripción:  Genera un heatmap de los genes significativos
#==============================================================

#' Genera un heatmap de los genes significativos 
#' y los guarda en un fichero en el directorio de resultados
#' si el parámetro a_r no es nulo el grafico se anota con
#' los clusters asignados a cada gen
#' @param ms_o matriz escalada y ordenada
#' @param a_c_o
#' @param a_r
#' @param m
#' @return No devuelve nada

get_heatmap_sig<-function(ms_o, a_c_o, a_r, m){
  # Si están anotadas las filas crea heatmap por cluster
  if (!is.null(a_r)){
    heat<-pheatmap(ms_o,
                   annotation_col = a_c_o,
                   annotation_row = a_r,
                   clustering_method = "ward.D2",
                   cluster_cols = FALSE,       
                   cluster_rows = TRUE,        
                   show_rownames = TRUE,    
                   fontsize_row = 3,
                   fontsize_col = 5,
                   color = colorRampPalette(c("blue", "white", "red"))(299),
                   main = paste(m, "Cluster de genes significativos"),
                   filename = paste0("results/images/heatmap_cluster_",m,".png"))
    
    clusters <- sort(unique(a_r$cluster))
    for (cl in clusters) {
      genes_cl <- rownames(a_r)[a_r$cluster == cl]
      ms_cl <- ms_o[rownames(ms_o) %in% genes_cl, , drop = FALSE]
      
      if (nrow(ms_cl) > 1) {
        pheatmap(ms_cl,
                 annotation_col = a_c_o,
                 clustering_method = "ward.D2",
                 cluster_cols = FALSE,
                 cluster_rows = TRUE,
                 show_rownames = TRUE,
                 fontsize_row = 3,
                 fontsize_col = 5,
                 color = colorRampPalette(c("blue", "white", "red"))(299),
                 main = paste(m, "- Cluster", cl, "(", nrow(ms_cl), ")"),
                 filename = paste0("results/images/heatmap_cluster_", m, "_cl", cl, ".png"))
      }
    }
  }else {
    heat<-pheatmap(ms_o,
                   annotation_col = a_c_o,
                   scale = "none",
                   cluster_cols = FALSE, 
                   show_rownames = FALSE,
                   fontsize_col = 5,
                   color = colorRampPalette(c("blue", "white", "red"))(299),
                   main = paste(m, "Heatmap -  genes significativos"),
                   filename = paste0("results/images/heatmap_sig_",m,".png"))
  }

}

#==============================================================
# Función: get_pca
# Descripción:  Genera un gráfico PCA de genes agrupados 
# por cluster
#==============================================================

#' Genera un gráfico PCA de genes agrupados por cluster
#' y los guarda en un fichero en el directorio de resultados
#' @param ms matriz escalada
#' @param gc data frame de genes y clusters
#' @param m tipo de método
#' @return No devuelve nada

get_pca<-function(ms, gc, m){
  # Transponer: genes como observaciones, muestras como variables
  pca_res <- prcomp(ms, center = TRUE, scale. = FALSE)
  
  var_exp <- pca_res$sdev^2 / sum(pca_res$sdev^2)
  percent_var <- round(100 * var_exp[1:2], 1)
  
  # Extraer componentes principales
  pca_df <- as.data.frame(pca_res$x[, 1:2])  # PC1 y PC2

  pca_df$cluster <- factor(gc$cluster[match(rownames(pca_df), rownames(gc))])
  
  # Visualizar
  library(ggplot2)
  
  cl <- sort(unique(pca_df$cluster))
  
  col_cl <- get_colors(cl)
  
  ggpca<-ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster,
                            fill = cluster)) +
    stat_ellipse(type = "norm", alpha = 0.2, geom = "polygon") +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(values = col_cl) +
    labs(title = "PCA de genes agrupados por clúster de expresión",
         subtitle = m,
         x = paste0("PC1 (", percent_var[1], "%)"),
         y = paste0("PC2 (", percent_var[2], "%)"), 
         color = "Clúster") +
    theme_minimal()
  
  ggsave(paste0("results/images/pca_cluster_",m,".png"), plot=ggpca)
  
}

#==============================================================
# Función: get_boxplot
# Descripción:  Genera un boxplot de la expresión media 
# por cluster
#==============================================================

#' Genera un boxplot de la expresión media  por cluster
#' y lo guarda en un fichero en el directorio de resultados
#' @param c_a matriz escalada
#' @param m tipo de método
#' @return No devuelve nada

get_boxplot<- function(c_a, m){
  ggboxplot<-ggplot(c_a, aes(x = trt_st, y = mean_expr, fill =trt_st)) +
    geom_boxplot() +
    facet_wrap(~ cluster, scales = "free_y") +
    theme_minimal() +
    labs(title = "Expresión media por clúster de genes",
         x = "Grupo experimental", y = "Expresión media (escalada)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  ggsave(paste0("results/images/boxplot_cluster_",m,".png"), plot=ggboxplot)
}

#==============================================================
# Función: get_boxplot_anova
# Descripción:  Genera un boxplot de la expresión media 
# por cluster
#==============================================================

#' Genera un boxplot de la expresión media  por cluster
#' despues de hacer un ANOVA - Tukey y lo guarda en un fichero
#' en el directorio de resultados
#' @param c_a matriz escalada
#' @param a_t  
#' @param m tipo de método
#' @return No devuelve nada

get_boxplot_anova<-function(c_a, a_t, m){
  library("ggpubr")

  a_t<-a_t %>%
    mutate(group1 = factor(group1, levels = levels(c_a$trt_st)),
           group2 = factor(group2, levels = levels(c_a$trt_st)))
  plot_anova_tukey <- ggplot(c_a, 
                             aes(x = trt_st, 
                                 y = mean_expr, 
                                 fill = trt_st)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
    facet_wrap(~ cluster, scales = "free_y") +
    theme_minimal() +
    labs(title = "Expresión media por clúster (ANOVA + Tukey HSD)",
         x = "Tratamiento_Etapa", 
         y = "Expresión media (Mean Expression)") +
    theme(axis.text.x = element_text(angle = 45, 
                                     hjust = 1, 
                                     size = 8),
          plot.title = element_text(hjust = 0.5)) +
    stat_pvalue_manual(data = a_t,
                       x ="group1",
                       xend = "group2",
                       y.position = "y.position",
                       label = "p.adj",
                       bracket.size = 0.4,
                       tip.length = 0.01,
                       size = 3,
                       group_by = "cluster")
  
  ggsave(paste0("results/images/boxplot_anova_tukey_",m,".png"), plot=plot_anova_tukey)
  
}

#==============================================================
# Función: get_boxplot_tukey_cld
# Descripción:  Genera un boxplot de la expresión media 
# por cluster con letras de significancia
#==============================================================

#' Genera un boxplot de la expresión media  por cluster
#' después de hacer comparaciones basadas en ANOVA - Tukey 
#' y lo guarda en un fichero en el directorio de resultados
#' asigna letras de significancia
#' 
#' @param d_c matriz escalada
#' @param cl identificador de cluster  
#' @param cld tipo de método
#' @param m string con el tipo de enfoque
#' @return No devuelve nada

get_boxplot_tukey_cld<-function(d_c, cl, cld, m){
  
  cld_data <- cld %>%
    as.data.frame() %>%
    dplyr::select(trt_st, .group)
  
  d_c <- d_c %>%
    dplyr::left_join(cld_data, by = "trt_st")
  
  # Graficar con letras de significancia
  p <- ggplot(d_c, aes(x = trt_st, y = mean_expr, fill = trt_st)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    geom_text_repel(data = cld_data, aes(x = trt_st, 
                                         y = max(d_c$mean_expr, 
                                                 na.rm = TRUE) + 0.2, 
                                         label = .group),
                    inherit.aes = FALSE,
                    size = 4,
                    fontface = "bold",
                    show.legend = FALSE) +
    labs(title = paste("Clúster", cl, "- Comparaciones por condición"),
         subtitle = m,
         x = "Condición experimental",
         y = "Expresión") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  ggsave(paste0("results/images/boxplot_cluster", cl, "_", m, ".png"),
         plot = p, width = 12, height = 6)
  

    #p <- ggplot(d_c, aes(x = trt_st, 
    #                   y = mean_expr,
    #                   fill = trt_st)) +
    #geom_boxplot(outlier.shape = NA) +
    #geom_jitter(width = 0.2, alpha = 0.5) +
    #geom_text(data = cld, aes(x = trt_st, 
    #                          y = max(d_c$mean_expr) * 1.1, 
    #                          label = .group),
    #          inherit.aes = FALSE, size = 3) +
    #theme_minimal() +
    #labs(title = paste("Expresión media en clúster ", cl, " (TUKEY + CLD)"),
    #     subtitle = m,
    #     x = "Tratamiento_Etapa", 
    #     y = "Expresión media (Mean Expression)") +
    #theme(axis.text.x = element_text(angle = 90, 
    #                                 hjust = 1, 
    #                                 size = 8),
    #      plot.title = element_text(hjust = 0.5))
  
  # Guardar el gráfico
  #ggsave(paste0("results/images/boxplot_tukey_cluster_", cl, "_", m, ".png"),
  #       plot = p,
  #       width = 12, height = 8, dpi = 300)
}

get_boxplot_tukey_cld_all_states <- function(d_c_all, cl, cld_all, m) {
  library(ggplot2)
  library(dplyr)
  
  # Añadir letras de significancia
  cld_df <- cld_all %>%
    rename(trt_st = trt_st, letra = .group) %>%
    distinct(trt_st, letra)
  
  d_c_all <- d_c_all %>%
    left_join(cld_df, by = "trt_st")
  
  # Boxplot
  p <- ggplot(d_c_all, aes(x = trt_st, y = mean_expr, fill = treatment)) +
    geom_boxplot(alpha = 0.7) +
    geom_text(aes(label = letra), stat = "summary", fun = mean, vjust = -0.8, size = 4) +
    facet_wrap(~ stage, scales = "free_x") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = paste("Expresión por condición en clúster", cl),
         x = "Condición", y = "Expresión")
  
  ggsave(paste0("results/images/boxplot_allstates_cluster", cl, "_", m, ".png"),
         plot = p, width = 12, height = 6)
}


get_dotplot <- function(e_g, cl,m){
  
  gdotplot<-dotplot(e_g,
                    showCategory = 20,
                    title = paste0("Cluster ",cl, "-GO enrichment"),
                    #subtitle = m,
                    font.size = 10,
                    label_format = 50) +
    theme_minimal()

  ggsave(paste0("results/images/dotplot_cluster_", cl,"-",m,".png"), plot=gdotplot)
  
}

get_barplot <- function(e_g, cl, m){
  gbarplot<-barplot(e_g,
                    showCategory = 15,
                    title = paste0("Cluster ",cl, "-GO enrichment", "(", m, ")"),
                    font.size =10)+
    theme_minimal()
  ggsave(paste0("results/images/barplot_cluster_", cl,"-",m,".png"), plot=gbarplot)
}

get_vennplot <- function(v_l){
  venn_plot <- venn.diagram(
    x = v_l,
    filename = "results/images/venn_lmdme_deseq2.png",
    imagetype = "png",
    output = TRUE,
    main = "DEGs comparados entre lmdme y DESeq2",
    col = "black",
    fill = c("lightblue", "lightgreen"),
    alpha = 0.5,
    cex = 1.5,
    fontface = "bold",
    cat.cex = 1.2,
    cat.pos = c(-20, 20),
    cat.dist = 0.05
  )
}

get_heatmap_enricher<-function(e_g, cl, m, tipo){
  
  matriz_heatmap <- e_g %>%
    filter(p.adjust<0.05) %>%
    mutate(log_padj = -log10(p.adjust)) %>%
    dplyr::select(comparison, Description, log_padj) %>%
    pivot_wider(names_from = comparison, values_from = log_padj, values_fill = 0) %>%
    column_to_rownames("Description") %>%
    as.matrix()
  
  if (nrow(matriz_heatmap) < 2 || ncol(matriz_heatmap) < 2) {
    cat("⚠️  La matriz del heatmap no tiene suficientes filas o columnas para generar el gráfico (mínimo 2).\n")
    return()
  }
  
  # Dibujar el heatmap
  pheatmap(matriz_heatmap,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "none",
           fontsize_row = 8,
           fontsize_col = 10,
           main = paste0("Heatmap de términos ",tipo, " enriquecidos",
                         " en clúster ", cl,
                         " (",m,")"),
           filename = paste0("results/images/heatmap_enrich_cluster",cl, "_", tipo, "_", m,".png"))
}

plot_dotplot_enrichment <- function(e_g, 
                                    cl,
                                    tipo = "GO", 
                                    top_terms = 15, 
                                    m,
                                    output_name = "dotplot_compare") {
  
  
  # Selecciona los términos más frecuentes (por número de comparaciones)
  top_ids <- e_g %>%
    filter(p.adjust < 0.05) %>%
    dplyr::count(ID, sort = TRUE) %>%
    slice_head(n = top_terms) %>%
    pull(ID)
  

  if (length(top_ids) == 0) {
    cat("⚠️  La tabla de resumen está vacía. No se generará el dotplot.\n")
    return()
  }
  
  dt_plot <- e_g %>%
    filter(ID %in% top_ids)
    
  
  p <- ggplot(dt_plot, aes(x = comparison, y = Description, size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = paste("Dotplot términos ", tipo, " en cluster", cl, " (", m, ")"),
      x = "Comparación",
      y = "Término enriquecido"
    )
  
  ggsave(paste0("results/images/", output_name, "_", "cluster", cl, "_",tipo, "-", m,".png"), plot = p, width = 12, height = 7)
  cat("✅ Dotplot conjunto de enriquecimiento guardado correctamente.\n")
}

