############################################################
# Archivo: procedures.R
# Descripción: Funciones del proyecto
# Autor: Gemma Gariglio Viejo
# Fecha: 02/06/2025
############################################################
source("R/utils.R")

library(readxl)
library(here)

#==============================================================
# Función: prepare_counts
# Descripción:  Prepara los datos de recuentos para el análisis
#==============================================================

#' Prepara los datos de recuentos para el análisis
#' eliminando las columnas que no son necesarias
#' y estableciendo el nombre de las filas con
#' como el identificador de genes
#' @param x matriz con los recuentos brutos
#' @return data.frame con recuentos

prepare_counts<-function(x){
  x<- x %>% 
    as.data.frame()
  #readxl::read_excel("data/gene_count.xlsx")
  x<-x %>% 
    dplyr::select(-c(gene_name, 
                     gene_start, 
                     gene_end, 
                     gene_strand, 
                     gene_length,
                     gene_biotype, 
                     gene_description,
                     tf_family,
                     gene_chr)) %>%
    column_to_rownames("gene_id") #%>%
  
  return(x)
}

#==============================================================
# Función: prepare_metadata
# Descripción: Prepara los metadatos para el análisis
#==============================================================

#' Prepara los metadatos para el análisis
#' convirtiendo las variable categóricas en factores 
#' cambiando el contenido de alguna variable
#' creando nuevas variables para análisis posteriores
#' @param md data.frame de metadatos original
#' @return data.frame con metadatos procesados y nuevas variables

prepare_metadata<-function(md){
  
  md<-md %>%
    as.data.frame()
  
  # Factorización de variable categóricas
  
  md$stage<-factor(md$stage)
  md$treatment<-factor(md$treatment)
  
  # Corrección de valores NA en variable dose
  md$dose[is.na(md$dose)]<-"0"
  md$dose<-factor(md$dose)
  
  # Creación de variables
  # variable rep con el identificador de replica
  md$rep<- ave(1:nrow(md), 
               md$stage, 
               md$treatment, 
               md$dose, 
                     FUN = seq_along)
  
  # variable con una versión corta del estado
  md <- md %>%
    mutate(st = case_when(
      stage == "embryo"  ~ "ee",
      stage == "provisioned_egg" ~ "pe",
      stage == "unprovisioned_egg"  ~ "ue"
    ))
  
  # variable trt_st que combina el contaminante con el estado
  md <- md %>%
    mutate(trt_st = paste(treatment, stage, sep = ":"))
  md$trt_st<-factor(md$trt_st)
  
  # variable trt_dose que combina el contaminante con la dosis
  md <- md %>%
    mutate(trt_dose = paste(treatment, dose, sep = "_"))
  md$trt_dose<-factor(md$trt_dose)
  md$trt_dose<-droplevels(md$trt_dose)
  
  # variable iden contiene un identificador de la muestra mas 
  # representativo
  md <- md %>%
    mutate(iden = paste(treatment, st, dose, rep, sep = "_"))

  return(md)
}

#==============================================================
# Función: get_table_muestras
# Descripción: Crea una tabla de las muestras según contaminante
# y estado de desarrollo
#==============================================================

#' Crea una tabla de las muestras según contaminante
#' y estado de desarrollo y la guarda en un fichero
#' con formato word 
#' 
#' @return no devuelve nada

get_table_muestras<-function(){
  tabla_contingencia <- as.data.frame(base::table(metadata$treatment, metadata$stage))
  colnames(tabla_contingencia) <- c("Tratamiento", "Etapa", "Nº Muestras")
  tabla_flex <- flextable(tabla_contingencia) %>%
    autofit()
  
  # Exportar a Word
  doc <- read_docx() %>%
    body_add_par("Tabla de Muestras: Contaminante vs Etapa", style = "heading 1") %>%
    body_add_flextable(tabla_flex)
  
  print(doc, target = "results/intermediateData/tabla_muestras.docx")
}

#==============================================================
# Función: get_table_calidad
# Descripción: Crea una tabla con el tamaño de las bibliotecas
# de RNA-seq
#==============================================================

#' Crea una tabla con el tamaño de las bibliotecas
#' de RNA-seq y lo guarda en un fichero con formato word 
#' 
#' @return no devuelve nada

get_table_calidad<-function(){
  library(gtsummary)
  
  # Transponer matriz para que las muestras sean filas
  counts_t <- as.data.frame(t(counts))
  counts_t <- counts_t %>%
    rownames_to_column("sample")
  
  # Calcular tamaño de biblioteca y genes detectados
  summary_df <- counts_t %>%
    rowwise() %>%
    mutate(
      total_counts = sum(c_across(-sample)),
      detected_genes = sum(c_across(-sample) > 0)
    ) %>%
    ungroup() %>%
    dplyr::select(sample, total_counts, detected_genes)
  
  library(officer)
  library(flextable)
  
  # Crear tabla con flextable
  tabla_flex <- flextable(summary_df) %>%
    set_header_labels(
      sample = "Muestra",
      total_counts = "Tamaño de biblioteca",
      detected_genes = "Genes detectados"
    ) %>%
    autofit()
  
  tabla_estadisticas <- summary_df %>%
    dplyr::select(-sample) %>%  # solo variables numéricas
    tbl_summary(
      type = all_continuous() ~ "continuous2",  # permite múltiples estadísticas
      statistic = all_continuous() ~ c(
        "Media (±DE)" = "{mean} (±{sd})",
        "Mediana" = "{median}",
        "Mínimo" = "{min}",
        "Máximo" = "{max}"
      ),
      digits = all_continuous() ~ 0,
      label = list(
        total_counts ~ "Tamaño de biblioteca",
        detected_genes ~ "Genes detectados"
      )
    ) %>%
    bold_labels()
  
  tabla_estadisticas_flex <- as_flex_table(tabla_estadisticas)
  
  doc <- read_docx() %>%
    body_add_par("Resumen por muestra", style = "heading 1") %>%
    body_add_flextable(tabla_flex) %>%
    body_add_par("", style = "Normal") %>%
    body_add_par("Estadísticas descriptivas", style = "heading 1") %>%
    body_add_flextable(tabla_estadisticas_flex)
  
  # Guardar el archivo
  print(doc, target = "results/intermediateData/tabla_recuentos.docx")
}

#==============================================================
# Función: filter_0
# Descripción: Elimina los genes cuya suma de recuentos esta por 
# debajo en al menos dos muestras
#==============================================================

#' Elimina los genes cuya media de recuentos esta por debajo
#' de diez 
#' @return matriz filtrada

filter_0<-function(){
  
  keep <- apply(counts,1,mean)
  x <- counts[keep>=10, ]

  return(x)
}

#==============================================================
# Función: getLogCounts
# Descripción: Transforma por pseudoconteo la matriz de recuentos
#==============================================================

#' Transforma logaritmicamente por pseudoconteo la matriz
#' de recuento
#' @return matriz transformada logaritmicamente

get_logCounts<-function(){

    lgc <- log2(counts + 1)

  return(lgc)
}

#==============================================================
# Función: getLogCountsLong
# Descripción: Devuelve el contenido de una matriz en formato largo
#==============================================================

#' Devuelve los recuentos en formato largo
#' @param lc matriz de genes x muestra
#' @return data.frame en formato largo

get_logCountsLong<-function(lc){
  
  #convierte una matriz en data.frame
  lc_d<-as.data.frame(lc)
  
  # crea un columna con el nombre de las filas
  #lc_d$gene<-rownames(lc)
  
  # transforma el data.frame a formato largo y asigna nombre
  # a las columnas
  lc_d <- lc_d %>%
    pivot_longer(cols = everything(), #-gene,
                 names_to = "iden",
                 values_to = "logcount")
  
  # añade las columnas treatment y stage de metadata
  lc_d <- lc_d %>%
    left_join(metadata %>% dplyr::select(iden, 
                                         treatment, 
                                         stage, 
                                         trt_st), 
              by="iden") %>%
    arrange(iden)
  
  
  write.csv(lc_d, paste0("results/datafiles/long_","" , ".csv"), row.names = FALSE)
  return(lc_d)
}

#==============================================================
# Función: get_matrix_scaled
# Descripción: Escala la matriz de recuentos de los genes
# significativos
#==============================================================

#' Escala la matriz de recuentos
#' @param m matriz con genes form filas y muestras por columna
#' @param g_s vector con el listado de genes significativos

get_matrix_scaled<-function(x, g_s){
  
  # filtrado de los valores de recuentos para los genes significativos
  x<-x[rownames(x) %in% g_s,]
  # escalado
  ms<-t(scale(t(x)))
  
  return(ms)
}

#==============================================================
# Función: anota_col
# Descripción: crea un data frame con las anotaciones
# de las muestras de las condiciones que se le pasan
#==============================================================

#' crea un data frame con las anotaciones de muestras
#' @param cond1 string con el nombre de la condición
#' @param cond2 string con el nombre de la condición 
#' @return a_c data frame con anotaciones de muestras

anota_col<-function(cond1, cond2, orden) {
  # crea un data frame con las condiciones a estudiar
  if (is.null(orden)){
    a_c<-metadata %>%
      dplyr::select(cond1,cond2)
  } else {
    a_c<-metadata %>%
      arrange(cond1,cond2) %>%
      dplyr::select(cond1,cond2)
  }
  

  return(a_c)
}

#==============================================================
# Función: get_order_samples
# Descripción: ordena la matriz de expresión en función de las
# anotaciones de muestras y las condiciones que se le pasan y 
# el data frame de anotaciones de muestras
#==============================================================

#' crea un data frame ordenado
#' @param ms matriz de expresión escalada
#' @param a_c data frame de anotaciones de muestras
#' @return data frame 
 
get_order_samples<-function(ms, a_c){

  o_s<-rownames(a_c)
  ms_o<-ms[,o_s]

  return(ms_o)
}

#==============================================================
# Función: get_kmedoid_optimal_clusters
# Descripción: Obtiene el nº de clusters óptimos aplicando 
# método silhouette
#==============================================================

#' Obtiene el nº de clusters óptimos aplicando el método
#' silhouette
#' 
#' @param ms matriz escalada de expresión
#' @param max_clusters número máximo de cluster a comprobar
#' @param method método estadístico a aplicar
#' @return cluster óptimo

get_kmedoid_optimal_clusters <- function(ms, 
                                         max_clusters = 10, 
                                         method = "silhouette",
                                         m) {
  
  sil_width <- numeric(max_clusters-1)
  
  for (k in 2:max_clusters) {
    pam_fit <- pam(ms, k)
    sil_width[k-1] <- pam_fit$silinfo$avg.width
  }
  
  # Mejor número de clusters: el que maximiza el ancho promedio de silhouette
  best_k <- which.max(sil_width) + 1  # +1 porque el índice 1 corresponde a k=2
  
  library(fpc)
  pamk.best <- pamk(ms, krange=1:5, critout=TRUE)
  cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")
  
  nk=as.integer(pamk.best$nc)
  print(nk)
  
  prefix<-"mejor_nclus_lmdme"
  
  elb<-fviz_nbclust(ms, kmeans, method = "wss") +
    geom_vline(xintercept = 3, linetype = 2) +
    labs(title = "Método del Codo")
  
  ggsave(filename = paste0("results/images/",prefix, "_elbow", m,".png"), 
         plot = elb, width = 6, height = 4)
  
  sil<-fviz_nbclust(ms, kmeans, method = "silhouette") +
    labs(title = "Método del Silhouette")
  
  ggsave(filename = paste0("results/images/",prefix, "_sil", m, ".png"), 
         plot = sil, width = 6, height = 4)
  
  gap_stat <- clusGap(ms, FUN = kmeans, nstart = 25, K.max = 10, B = 50)
  gap<-fviz_gap_stat(gap_stat) +
    labs(title = "Método del Gap Statistic")
  
  ggsave(filename = paste0("results/images/", prefix, "_gap.png"), 
         plot = gap, width = 6, height = 4)
  
  return(nk)
}

#==============================================================
# Función: get_clustering
# Descripción: asigna a cada gen a uno de los clusteres calculados
# por el método silhouette
#==============================================================

#' Asigna a cada gen a uno de los clusteres calculados por el
#' método silhouette
#' @param ms matriz escalada de expresión
#' @param k numero optimo de cluster a asignar
#' @return data.frame con los identificadores de genes y el cluster asignado

get_clustering<-function(ms, k=3){
  
  # metodo pam para asignar clusteres
  clust<-pam(ms, k)
  g_c<-data.frame(cluster = factor(clust$clustering))
  
  return(g_c)
}

#==============================================================
# Función: anota_row
# Descripción: crea un data frame con el identificador de gene 
# y su cluster
#==============================================================

#' crea un data frame
#' @param g_c data frame 
#' @return data frame con identificador de gen y su cluster

anota_row<-function(g_c) {
  # crea un data frame con el identificador de gene y su cluster

  a_r<-g_c["cluster"]

  return(a_r)
}

#==============================================================
# Función: get_table_cluster
# Descripción: Crea un data.frame en formato largo con los datos
# de expresión, el cluster asignado y los metadatos
#==============================================================

#' Crea un data.frame en formato largo con los datos
#' de expresión, el cluster asignado y los metadatos
#' 
#' @param ms matriz escalada de expresión
#' @param g_c data.frame con los identificadores de genes y cluster
#' @return data.frame en formato largo

get_table_cluster<-function(ms, g_c){
  ms_df<-as.data.frame(ms) %>%
    rownames_to_column("gene")
  
  e_long<- ms_df %>%
    pivot_longer(-gene,
                 names_to = "iden",
                 values_to = "expression") %>%
    inner_join(metadata, by ="iden") %>%
    inner_join(g_c %>% rownames_to_column("gene"), by = "gene")
  
  t_c <- e_long %>%
    dplyr::select(gene,iden,expression, treatment, stage, trt_st, cluster)
  
  return(t_c)
}

#==============================================================
# Función: get_anova_tukey_by_state
# Descripción: 
#==============================================================

#  Realiza pruebas ANOVA y Sidak para las compraciones con respecto 
#' al control y les asigna letras de significancia
#' @param c_a matriz escalada de expresión
#' @param m cadena con el nombre del método utilizado
#' @return lista de data frames

get_anova_tukey_by_state <- function(c_a, m) {
  library("emmeans")
  library("multcomp")
  
  # define listas para utilizar más adelante
  tukey_list <- list()
  dominant_conditions <- list()
  d_c_all <- list()
  cld_all <- list()
  
  # Extraer las etapas únicas (parte después de ":")
  estados <- unique(c_a$stage)
  
   # Para cada etapa
   for (estado in estados) {
    # filtra por etapa
    st_data <- c_a %>% filter(stage == estado)
    # obtiene los identificadores unicos de los cluster
    clusters <- unique(st_data$cluster)
    
    # para cada cluster
    for (cl in clusters) {
      #filtra por cluster
      d_c <- st_data %>% filter(cluster == cl)
      
      # Si no hay más de dos condiciones se sale del bucle
      if (length(unique(d_c$trt_st)) < 2) next
      
      # calcula anova y medias marginales
      aov_m <- aov(mean_expr ~ trt_st, data = d_c)
      emm <- emmeans(aov_m, ~ trt_st)
      
      # Determinar el control correspondiente a ese estado
      ref_ctrl <- paste0("CTR:", estado)
      if (!(ref_ctrl %in% levels(emm@grid$trt_st))) next
      
      #realiza los contrastes
      contr <- contrast(emm, method = "trt.vs.ctrl", 
                        ref = which(levels(emm@grid$trt_st) == ref_ctrl))

      #Crea el dataframe con los resultados
      tukey_r <- as.data.frame(contr) %>%
        mutate(group1 = ref_ctrl,
               group2 = str_trim(str_extract(contrast, "^[^-]+")),
               cluster = cl,
               estado = estado) %>%
        dplyr::select(cluster, estado, group1, group2, estimate, p.value)
      
      # crea una lista para cada cluster
      tukey_list[[paste(cl, estado, sep = "_")]] <- tukey_r
      
      # Haya las letras de significancia para boxplot
      cld_r <- cld(emm, Letters = letters, adjust = "tukey")
      
      # crea un dataframe con la comparación dominante
      dominant <- cld_r %>%
        as.data.frame() %>%
        arrange(desc(emmean)) %>%
        slice(1) %>%
        dplyr::select(trt_st, emmean)
      
      dominant_conditions[[paste(cl, estado, sep = "_")]] <- data.frame(
        cluster = cl,
        estado = estado,
        dominant_condition = dominant$trt_st,
        max_emmean = dominant$emmean
      )
      
      d_c_all[[as.character(cl)]] <- bind_rows(d_c_all[[as.character(cl)]], 
                                               d_c)
      cld_all[[as.character(cl)]] <- bind_rows(cld_all[[as.character(cl)]], 
                                               as.data.frame(cld_r))
    }
  }
  
  # une todos los dataframe en una unica lista
  tukey_table <- bind_rows(tukey_list)
  dom_table <- bind_rows(dominant_conditions)
 
  #colnames(tukey_table) <- c("Cluster", "estado", "Control", "Contaminante", "estimate", "p.value")
  
  show_table_word(tukey_table, m, 
                  "Prueba Sidak de los contrastes frente al control",
                  "contrastes")
  
  write.csv(tukey_table, paste0("results/datafiles/tukey_contrasts_by_state_", m, ".csv"), row.names = FALSE)
  write.csv(dom_table, paste0("results/datafiles/dominant_conditions_by_state_", m, ".csv"), row.names = FALSE)
  
  # genera un boxplot con todos los resultados
  for (cl in names(d_c_all)){
    get_boxplot_tukey_cld_all_states(d_c_all[[cl]], cl, cld_all[[cl]], m)
  }
  
  # genera un dataframe con las letras asignadas
  letras_table <- bind_rows(cld_all, .id = "cluster") %>%
    dplyr::select(cluster, trt_st, emmean, SE, df, .group) %>%
    rename(letra = .group)
  
  # Guardar tabla de letras asignadas en un fichero
  write.csv(letras_table,
            paste0("results/datafiles/letras_significancia_clusters_", m, ".csv"),
            row.names = FALSE)
  
  cat("✅ Comparaciones Tukey por estado completadas\n")
  
  return(list(tukey_table = tukey_table, dom_table = dom_table))
}


#==============================================================
# Función: get_significant_genes_by_condition
# Descripción: 
#==============================================================

#' Identifica genes diferencialmente expresados según condiciones 
#' experimentales comparadas con el grupo control, para cada clúster 
#' de genes y nivel de una condición.
#' Solo se consideran comparaciones donde el grupo control y la condición 
#' pertenecen al mismo estado. 
#' 
#' @param t_c data frame con expresión
#' @param t_t data frame con resultados de comparaciones múltiples
#' @param p_cutoff umbral de significancia para el valor p (default = 0.05)
#' @param ctr_p cadena que identifica el grupo control(default = "CTR:")
#' @param m cadena con el nombre del método utilizado ("lmdme", "DESEq2")
#' @return lista de data frames donde cada uno contiene los genes relevantes
#'         para la comparación  

get_significant_genes_by_condition <- function(t_c, 
                                               t_t, 
                                               p_cutoff = 0.05, 
                                               ctr_p = "CTR:",
                                               m = "lmdme") {
 
  genes_by_condition <- list()
  
  # obtiene un vector de cluster unicos
  clusters <- unique(t_t$cluster)
  
  for (cl in clusters) {
    # Filtra comparaciones significativas respecto del control
    # dentro del clúster
    tukey_cl <- t_t %>%
      filter(cluster == cl, p.value < p_cutoff) %>%
      filter(grepl(ctr_p, group1) | grepl(ctr_p, group2))
    
    # condición para salir si no hay comparaciones
    if (nrow(tukey_cl) == 0) next
    
    for (i in seq_len(nrow(tukey_cl))) {
      row <- tukey_cl[i, ]
      
      # Detecta cuál es la condición experimental y cual el control
      cond <- if (grepl(ctr_p, row$group1)) row$group2 else row$group1
      ctrl <- if (grepl(ctr_p, row$group1)) row$group1 else row$group2
      
      # Extrae los estados
      cond_stage <- str_split(cond, ":", simplify = TRUE)[, 2]
      ctrl_stage <- str_split(ctrl, ":", simplify = TRUE)[, 2]
      
      # condición para controlar si estamos en el mismo estado
      if (cond_stage != ctrl_stage) next
      
      # Calcula la media de expresión por gen en cada condición
      mean_expr <- t_c %>%
        filter(cluster == cl, trt_st %in% c(ctrl, cond)) %>%
        group_by(gene, trt_st) %>%
        summarise(mean_exp = mean(expression, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = trt_st, values_from = mean_exp)
      
      # Condición para comprobar que existen las condiciones
      if (!(cond %in% colnames(mean_expr)) || !(ctrl %in% colnames(mean_expr))) next
      
      # Selecciona los genes con mayor expresión en la condición experimental
      genes_sig <- mean_expr %>%
        filter(.data[[cond]] > .data[[ctrl]]) %>%
        mutate(cluster = cl,
               condition = cond,
               control = ctrl,
               estimate = row$estimate)
        #pull(gene)
      
      # condición para guardar si hay al menos 5 genes significativos
      if (length(genes_sig) >= 5) {
        key <- paste0("cl", cl, "_", gsub(":", "_", cond))
        genes_by_condition[[key]] <- genes_sig
        write.csv(genes_sig, paste0("results/datafiles/", "genes_relev_", m, "_", key, ".csv"), row.names = FALSE)
      }
      
    }
  }
  
  return(genes_by_condition)
}

#==============================================================
# Función: get_enrichment_by_significant_conditions
# Descripción: Ejecuta enriquecimiento funcional para un conjunto de genes
# por comparación y genera visualizaciones
#==============================================================

#' Ejecuta enriquecimiento funcional para un conjunto de genes
#' por comparación y genera visualizaciones.
#' 
#' Este bloque realiza enriquecimiento funcional usando la función `enricher()` 
#' de clusterProfiler para una lista de genes por comparación 
#' (e.g., condición vs control). Guarda los resultados, genera gráficos 
#' (dotplot y cnetplot), y construye un resumen general.
#' 
#' @genes_list vector con genes significativos de la comparación actual
#' @uni  universo de genes (todos los posibles)
#' @t2g tabla TERM2GENE (asociaciones entre términos y genes)
#' @t2n tabla TERM2NAME (nombres descriptivos de los términos)
#' @m cadena con el nombre del método utilizado ("lmdme", "DESEq2")
#' @p_cutoff umbral de significancia para el valor p (default = 0.05)
#' @tipo cadena con el tipo de termino ("GO", "KEGG")
#' @output_prefix cadena para darle nombre al fichero 
#' @return devuelve un data frame global con todos los términos enriquecidos.

get_enrichment_by_significant_conditions <- function(genes_list,
                                                     uni,
                                                     t2g,
                                                     t2n,
                                                     m = "lmdme",
                                                     p_cutoff = 0.05,
                                                     tipo = "GO",
                                                     output_prefix = "enrich_significant") {
  
  enrich_list<-list()
  
  for (name in names(genes_list)) {
    t_g <- genes_list[[name]]
    genes<-unique(t_g$gene)
    
    if (length(genes) < 5) {
      cat("⚠️  Pocos genes para", name, "- omitiendo...\n")
      next
    }
    
    enr <- enricher(gene = genes,
                    universe = uni,
                    TERM2GENE = t2g, 
                    TERM2NAME = t2n, 
                    pvalueCutoff = p_cutoff)
    
    if (!is.null(enr) && nrow(enr@result %>% filter(p.adjust < p_cutoff)) > 0) {
      res <- enr@result %>%
        mutate(comparison = name,
               cluster = unique(t_g$cluster),
               condition = gsub("^cl\\d+_", "", name),
               control = unique(t_g$control))
      
      cat("✅  Encontrado enriquecimeinto significativo para", name,"\n")
      
      enrich_list[[name]] <- res
      
      
      dplot <- dotplot(enr, showCategory = 15, font.size = 10) +
        ggtitle(name)
      ggsave(paste0("results/images/", 
                    m, "_", output_prefix, "_", 
                    tipo, "_dotplot_",name,".png"),
             plot = dplot, width = 10, height = 6)
     
      cnet<-cnetplot(enr, categorySize="pvalue",
                     circular = TRUE) +
        ggtitle(paste0("Cnetplot - ", name, "(", m ,")"))
      
      ggsave(paste0("results/images/", 
                    m, "_", output_prefix, "_", 
                    tipo, "_cnet_",name, ".png"),
             plot = cnet, width = 10, height = 6)
      
    } else {
        cat("⚠️  Sin enriquecimiento significativo en ", name, "\n")
    }
  }
  
  cat("Unificando resultados...\n")
  
  enrich_global <- bind_rows(enrich_list)
  
  if (nrow(enrich_global) != 0) {
    
    #show_table_word(enrich_global, m, 
    #                paste0("Terminos ", tipo, " enriquecidos - ", m),
    #                paste0("contrastes_", tipo, "_", m))
   
    write.csv(enrich_global, paste0("results/datafiles/", m, "_", 
                                    output_prefix, "_", 
                                    tipo, "_resumen.csv"), 
              row.names = FALSE)
    
    
    
    #get_heatmap_enricher(enrich_global, m, tipo)
    #plot_dotplot_enrichment(enrich_global,
    #                        tipo, 
    #                        10, 
    #                        m,
    #                        "dotplot_compare")
    
    cat("Generando heatmap y dotplot para cada cluster...\n")
    for (cl in unique(enrich_global$cluster)){
      df_cl <- filter(enrich_global, cluster == cl)
      
      if (nrow(df_cl) == 0) {
        cat("⚠️  La tabla de resumen está vacía. No se generará el dotplot.\n")
        next
      }
      
      get_heatmap_enricher(df_cl, cl, m, tipo)
      
      plot_dotplot_enrichment(df_cl, 
                              cl, 
                              tipo = tipo, 
                              top_terms = 15, 
                              m = m,
                              output_name = paste0("dotplot_compare_cluster", cl))
      
      
    }
    
    summ<-enrich_global %>%
      group_by(ID) %>%
      summarise(n_condiciones = n_distinct(condition)) %>%
      arrange(desc(n_condiciones)) %>%
      filter(n_condiciones > 1)  # términos comunes
    
    return(enrich_global)
  } else {
      return()
  }
  
  
}

#==============================================================
# Función: get_compara_enrich
# Descripción: compara los resultados de los dos enfoques
#==============================================================

get_compara_enrich<-function(tipo){
  
  if (tipo == "GO"){
    enrichment_lmdme <- enrich_conditions_GO_lmdme %>%
      mutate(method = "lmdme")
    
    enrichment_deseq <- enrich_conditions_GO_deseq %>%
      mutate(method = "DESeq2")
  } else{
    enrichment_lmdme <- enrich_conditions_KEGG_lmdme %>%
      mutate(method = "lmdme")
    
    enrichment_deseq <- enrich_conditions_KEGG_deseq %>%
      mutate(method = "DESeq2")
  }

  enrichment_all <- bind_rows(enrichment_lmdme, enrichment_deseq)
  
  enrichment_all_filtered <- enrichment_all %>%
    filter(p.adjust < 0.05)
  
  term_comparison <- enrichment_all_filtered %>%
    group_by(ID, Description, method) %>%
    summarise(n = n(), .groups = "drop")
  
  term_comparison_wide <- term_comparison %>%
    tidyr::pivot_wider(names_from = method, values_from = n, values_fill = 0)
  
  term_comparison_wide <- term_comparison_wide %>%
    mutate(shared = ifelse(lmdme > 0 & DESeq2 > 0, "Sí", "No"))
  
  ggplot(term_comparison, aes(x = method, 
                              y = Description, 
                              size = n, 
                              fill = method)) +
    geom_point(alpha = 0.7, shape = 21, color = "black") +
    scale_size_continuous(range = c(3, 10)) +
    labs(title = "Comparación de términos enriquecidos por método",
         x = "Método", 
         y = "Término GO/KEGG", 
         size = "Frecuencia") +
    theme_minimal(base_size = 14) +
    theme(axis.text.y = element_text(size = 6))
  
  ggsave(paste0("results/images/enrich_comparation_",tipo, ".png"), width = 8, height = 6)
  
  
  terms_lmdme <- unique(enrichment_lmdme$ID[enrichment_lmdme$p.adjust < 0.05])
  terms_deseq <- unique(enrichment_deseq$ID[enrichment_deseq$p.adjust < 0.05])
  
  venn.plot <- venn.diagram(
    x = list(LMDME = terms_lmdme, DESeq2 = terms_deseq),
    category.names = c("lmdme", "DESeq2"),
    filename = NULL,
    output = FALSE,
    fill = c("skyblue", "salmon"),
    alpha = 0.5,
    cex = 1.5,
    cat.cex = 1.2,
    cat.fontfamily = "Courier"
  )
  
  grid::grid.draw(venn.plot)
  
  png(paste0("results/images/venn_enrichment_comparison_",tipo,".png"), 
      width = 800, 
      height = 700)
  grid.draw(venn.plot)
  
  # Añadir título
  grid.text("Comparación de términos enriquecidos", 
            x = 0.5, 
            y = 0.95, 
            gp = gpar(fontsize = 18, 
                      fontface = "bold"))
  
  # Añadir subtítulo
  grid.text(paste0("Términos ", tipo, " significativos: LMDME vs DESeq2"), 
            x = 0.5, 
            y = 0.91, 
            gp = gpar(fontsize = 14, 
                      fontface = "italic"))
  
  dev.off()
  
  
  return()
  
}

get_enrich_gsea<-function(g_m,
                          m,
                          t2g,
                          t2n,
                          minSize = 10, 
                          maxSize = 500, 
                          pCutoff = 0.05) {
  
  #terms_name <- colnames(g_m)
  #terms_name <- terms_name[grepl("\\.x$", terms_name )]
  #terms_name <- gsub("\\.x$", "", terms_name )  # quitar sufijo
  
  #gsea_lists<-list()
  
  #for (name in terms_name) {
    #print(name)
    #rank_vector <- g_m %>%
    #  dplyr::select(gene, coef = all_of(paste0(name, ".x"))) %>%
    #  filter(!is.na(coef)) %>%
    #  arrange(desc(coef)) %>%
    #  deframe()
    
    #if (length(rank_vector) < 100) {
    #  warning("Termino omitido por tamaño insuficiente: ", name)
    #  next
    #}
    rank_vector<- -log10(g_m$padj)
    names(rank_vector)<-g_m$gene
    rank_vector<-rank_vector[!is.na(rank_vector)]
    rank_vector <- sort(rank_vector, decreasing = TRUE)
    gsea_result<-GSEA(geneList = rank_vector,
                      TERM2GENE = t2g,
                      TERM2NAME = t2n,
                      minGSSize = minSize,
                      maxGSSize = maxSize,
                      pvalueCutoff = pCutoff,
                      verbose = FALSE)
    
    #gsea_lists[[name]] <- gsea_result
  #}
  
 dot<-dotplot(gsea_result, showCategory = 20)
 
 ggsave(paste0("results/images/dotplot_gsea_go_",m,".png"), plot=dot)
  
  return(gsea_result)
  
}

heatmap_pheno_cluster<-function(c_a, m){
  # selecciona las columnas sample e iden de metadata
  iden<-metadata %>%
    dplyr::select(sample, iden)
  
  # añade a la tabla de fenotipos la variable iden y convierte
  # a factor las variables stage,treatment, trt_st
  pheno<-pheno %>%
    left_join(iden, by = "sample") %>%
    mutate(across(c(stage,treatment, trt_st), as.factor))
  
  # Añade al dataframe con los expresiones medias por cluster
  # los fenotipos correspondientes a cada muestra
  expr_pheno<-c_a %>%
    left_join(pheno, by = c("iden" = "iden"))
  
  # Elimina del dataframe de expresión las variable 
  # que no son necesarias, convierte a formato largo los valores 
  # de los fenotipos, agrupa por cluster y fenotipo
  # y calcula la correlación entre la expresión media y el fenotipo
  cor_results <- expr_pheno %>%
    dplyr::select(-ends_with("_SE")) %>%
    pivot_longer(cols = ends_with("_Mean"), names_to = "phenotype", values_to = "value") %>%
    group_by(cluster, phenotype) %>%
    dplyr::summarise(cor_test = list(cor.test(mean_expr, value, method = "pearson")),
                     .groups = "drop") %>%
    mutate(correlation = map_dbl(cor_test, ~.x$estimate),
           p_value = map_dbl(cor_test, ~.x$p.value))
  
  # Selecciona la variables y convierte a formato ancho
  cor_matrix <- cor_results %>%
    dplyr::select(cluster, phenotype, correlation) %>%
    pivot_wider(names_from = cluster, values_from = correlation) %>%
    column_to_rownames("phenotype")
  
  # genera el heatmap por cluster
  pheatmap(as.matrix(cor_matrix),
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           filename = paste0("results/images/heatmap_pheno_", m, ".png"),
           main = paste0("Correlación expresión (clusters) vs. fenotipos - " , m))
}

heatmap_pheno_enrich<-function(e_c, tipo, m){
  
  zscore_matrix <- e_c %>%
    group_by(ID, condition) %>%
    summarise(zScore = mean(zScore, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = condition, values_from = zScore, values_fill = list(zScore = 0)) %>%
    column_to_rownames("ID")
  
  colnames(zscore_matrix)<-sub("_", ":", colnames(zscore_matrix))
  
  pheno_summary <- pheno %>%
    group_by(trt_st) %>%
    summarise(across(ends_with("_Mean"), mean, na.rm = TRUE)) %>%
    ungroup()
  
  pheno_matrix <- pheno_summary %>%
    column_to_rownames("trt_st") %>%
    t()  # transponer: filas = fenotipos, columnas = condición
  
  common_conds <- intersect(colnames(zscore_matrix), colnames(pheno_matrix))
  
  zscore_matrix <- zscore_matrix[, common_conds]
  
  pheno_matrix <- pheno_matrix[, common_conds]
  
  cor_matrix <- cor(t(zscore_matrix), t(pheno_matrix), method = "pearson")
  
  # Convertir a formato largo
  cor_df <- as.data.frame(as.table(cor_matrix)) %>%
    rename(GO_term = Var1, Phenotype = Var2, Correlation = Freq)
  
  # Filtrar los top 5 términos más correlacionados (+/-) por fenotipo
  top_go_terms <- cor_df %>%
    group_by(Phenotype) %>%
    arrange(desc(abs(Correlation))) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    distinct(GO_term) %>%
    pull()
  
  # Filtrar matriz original
  cor_matrix_top <- cor_matrix[rownames(cor_matrix) %in% top_go_terms, ]
  
  go_desc <- e_c %>%
    distinct(ID, Description)
  desc_lookup <- go_desc$Description
  names(desc_lookup) <- go_desc$ID
  
  # reemplaza los nombre de las filas por las descripciones si existen
  # y si no deja el ID
  rownames(cor_matrix_top) <- ifelse(!is.na(desc_lookup[rownames(cor_matrix_top)]),
                                     desc_lookup[rownames(cor_matrix_top)],
                                     rownames(cor_matrix_top))
  
  # recorta 50 caraceteres las descripciones largas
  rownames(cor_matrix_top) <- stringr::str_trunc(rownames(cor_matrix_top), width = 50)
  
  pheatmap(cor_matrix_top,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           main = paste0("Correlación Términos ", tipo, " vs. Fenotipos -", m),
           fontsize_row = 8,
           fontsize_col = 9,
           display_numbers = TRUE,
           number_color = "black",
           filename = paste0("results/images/heatmap_correlacion_",tipo,
                             "_fenotipos_", m, ".png"))
}

#==============================================================
# Función: heatmap_pheno_enrich_cluster
# Descripción: genera un heatmap de la correlaciónentre los 
# terminos GO/KEGG y las variables fenotipicas
#==============================================================

#  Genera un heatmap de la correlaciónentre los 
#  terminos GO/KEGG y las variables fenotipicas
#' @param e_c dataframe con resultado de terminos enriquecidos
#' @param tipo string con el tipo de termino (GO o KEGG)
#' @param m cadena con el nombre del método utilizado
#' @return lista de data frames
#' 
heatmap_pheno_enrich_cluster<- function(e_c, tipo, m) {
  pheno_summary <- pheno %>%
    group_by(trt_st) %>%
    summarise(across(ends_with("_Mean"), mean, na.rm = TRUE)) %>%
    ungroup()
  
  pheno_matrix <- pheno_summary %>%
    column_to_rownames("trt_st") %>%
    t()
  
  go_desc <- e_c %>%
    distinct(ID, Description)
  desc_lookup <- go_desc$Description
  names(desc_lookup) <- go_desc$ID
  
  # Dividir por cluster
  cluster_list <- split(e_c, e_c$cluster)
  
  for (cl in names(cluster_list)) {
    cluster_data <- cluster_list[[cl]]
    
    # Dividir por condición dentro del cluster
    cluster_data$condition <- sub("_", ":", cluster_data$condition)
    condition_list <- split(cluster_data, cluster_data$condition)
    
    # Asegurar unicidad y generar matriz zScore
    zscore_matrix <- cluster_data %>%
      group_by(ID, condition) %>%
      summarise(zScore = mean(zScore, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = condition, values_from = zScore, values_fill = list(zScore = 0)) %>%
      column_to_rownames("ID")
      
      #colnames(zscore_matrix) <- sub("_", ":", colnames(zscore_matrix))
      
      common_conds <- intersect(colnames(zscore_matrix), colnames(pheno_matrix))

      zscore_matrix <- zscore_matrix[, common_conds, drop = FALSE]
      pheno_sub <- pheno_matrix[, common_conds, drop = FALSE]
      
      if (ncol(zscore_matrix) < 2 || ncol(pheno_sub) < 2) next
      
      cat("Pase por la condicion if\n")
      
      cor_matrix <- cor(t(zscore_matrix), t(pheno_sub), method = "pearson")
      
      # Seleccionar top términos correlacionados
      cor_df <- as.data.frame(as.table(cor_matrix)) %>%
        rename(GO_term = Var1, Phenotype = Var2, Correlation = Freq)
      
      top_go_terms <- cor_df %>%
        group_by(Phenotype) %>%
        arrange(desc(abs(Correlation))) %>%
        slice_head(n = 5) %>%
        ungroup() %>%
        distinct(GO_term) %>%
        pull()
      
      cor_matrix_top <- cor_matrix[rownames(cor_matrix) %in% top_go_terms, , drop = FALSE]
      
      rownames(cor_matrix_top) <- ifelse(!is.na(desc_lookup[rownames(cor_matrix_top)]),
                                         desc_lookup[rownames(cor_matrix_top)],
                                         rownames(cor_matrix_top))
      rownames(cor_matrix_top) <- stringr::str_trunc(rownames(cor_matrix_top), width = 50)
      
      asterisks_matrix <- ifelse(abs(cor_matrix_top) >= 0.85, "*", "")
      
      # Evita errores con matrices vacías
      if (nrow(cor_matrix_top) > 1 && ncol(cor_matrix_top) > 1) {
        pheatmap(cor_matrix_top,
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 main = paste0("Terminos ", tipo, " vs Fenotipos - Cluster ", cl,
                               " - ", m),
                 fontsize_row = 8,
                 fontsize_col = 9,
                 display_numbers = asterisks_matrix,
                 number_color = "black",
                 filename = paste0("results/images/heatmap_correlacion_", tipo,
                                   "_fenotipos_", m, "_cluster_", cl, ".png"))
      }
    }
}



