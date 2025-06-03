#######################################################################
# Archivo: control_calidad.R
# Descripción: Código de la opción 2 del menú 
# Autora: Gemma Gariglio Viejo
# Fecha: 02-06-2025
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
  
  cat("Filtrando los genes con recuentos nulos en todas las muestras...\n")
  print(dim(counts))
  counts<-filter_0()
  assign("counts",counts, envir = .GlobalEnv)
  print(dim(counts))
  
  
}