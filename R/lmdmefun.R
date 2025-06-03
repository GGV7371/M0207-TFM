############################################################
# Archivo: lmdmefun.R
# Descripción: Funciones propias del enfoque LMDME
# Autora: Gemma Gariglio Viejo
# Fecha: 02-03-2025
############################################################

#==============================================================
# Función: join_transcripts
# Descripción:  Une los transcritos en un unico gen
# en caso de que existan 
#==============================================================

#' Une los transcritos en un único gen
#' en caso de que existan 
#'
#' @param logcl matriz de recuentos transformada
#' @return dataframe 
#' 
join_transcripts <- function(logcl){
  ids<-as.character(rownames(logcl))
  ld <- 2^logcl #use the selected normalization
  AggData<-(aggregate(ld~(ids), FUN=sum))
  df<-as.data.frame(AggData[,-1])
  rownames(df)<-as.character(AggData[,1])

  print(dim(df))
  
  return(df)
}


# FUNCIONES NO UTILIZADAS EN LA ULTIMA VERSION
get_centrado_ctr <- function(df, m, id_col = "iden",
                             cond1 = "treatment", byc = "CTR",
                             cond2 = "stage") {
  
  # Transpone y une con metadata
  dfm <- df %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = id_col) %>%
    left_join(m, by = id_col)  
  
    #left_join(m, by = id_col)

  m <- as.matrix(df)
  
  # Calcula medias por grupo dentro del grupo de control
  niveles<-unique(dfm[[cond2]])

  for (nivel in niveles) {
    idx_control <- dfm[[cond1]] == byc & dfm[[cond2]] == nivel
    medc <- apply(m[, idx_control], 1, mean)

    idx_total <- dfm[[cond2]] == nivel
    m[, idx_total] <- m[, idx_total] / medc
  }
  
  # Log transform
  ml <- log2(m + 0.001)
  
  return(ml)
}


get_loadings<-function(){
  
  # Extraer loadings (pesos de los genes en la componente principal asociada al efecto)
  loadings <- modelo_decomp@loadings  # o modelo_decomp$loadings si es S3
  
  # Extraer p-values de cada gen
  pvalues_genes <- modelo_decomp@p.value  # o modelo_decomp$p.value
  
  # Definir genes significativos (p-value < 0.05)
  genes_significativos <- names(pvalues_genes)[pvalues_genes < 0.05]
  
  # Crear un data frame con información
  tabla_resultados <- data.frame(
    gen = genes_significativos,
    pvalue = pvalues_genes[genes_significativos],
    loading = loadings[genes_significativos],
    direccion = ifelse(loadings[genes_significativos] > 0, "Sobrerregulado", "Infrarregulado")
  )
  
  
}