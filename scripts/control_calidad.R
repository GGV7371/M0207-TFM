#######################################################################
# Archivo: control_calidad.R
# Descripción: Código de la opción 2 del menú 
# Autr: Gemma Gariglio Viejo
# Fecha: 
# Dependencias: utils.R, procedrures.R, lmdme.R, deseq2.R
#######################################################################

source("R/graphics.R")
source("R/procedures.R")

throw_control_calidad<-function(){
  
  counts<-counts_origin
  assign("counts",counts, envir = .GlobalEnv)
  
  cat("Generando análisis exploratorio de datos brutos...\n")
  
  get_table_calidad()
  
  cat("Generando barplot de recuentos por muestra...\n")
  
  get_dim_biblio()
  
  cat("Generando boxplot y density plot de recuentos log2 por muestra...\n")
  
  logcounts<-get_logCounts()
  
  logCountsLong<-get_logCountsLong(logcounts)
  
  get_box_recuentos(logCountsLong,"log")
  
  cat("Generando PCA ...\n")
  #get_pca_brutos()
  
  cat("Filtrando los genes con recuentos nulos en todas las muestras...\n")
  print(dim(counts))
  counts<-filter_0()
  assign("counts",counts, envir = .GlobalEnv)
  print(dim(counts))
  
  #logcounts<-get_logCounts()
  
  #logcountsloess <- normalizeCyclicLoess(logcounts)
  #assign("logcountsloess",logcountsloess, envir = .GlobalEnv)
  
  #logCountsloessLong<-get_logCountsLong(logcountsloess)
  
  #get_box_recuentos(logCountsloessLong, "loess")
  
}