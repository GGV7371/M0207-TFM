#######################################################################
# Archivo: load.R
# Descripción: Código de la opción 1 del menú 
# Autr: Gemma Gariglio Viejo
# Fecha: 
# Dependencias: utils.R, procedrures.R, lmdme.R, deseq2.R
#######################################################################

source("scripts/lmdme.R")
source("scripts/deseq2.R")
source("R/utils.R")
source("R/procedures.R")

#==============================================================
# Función: throw_load
# Descripción: Carga los ficheros de trabajo y los prepara 
#==============================================================

#' Instala y carga los paquetes de R y de Bioconductor
#' Carga los ficheros de recuentos y metadatos y los prepara
#' para el posterior análisis
#' 
#' Consiste en realizar:
#' - Instala y carga los paquetes de R y de Bioconductor
#' - Solicitar al usuario la ubicación de los ficheros de:
#' - - recuentos
#' - - metadatos
#' - Preparar los datos para el análisis
#' 
#' @return No devuelve ningún valor. Una vez terminado el análisis
#' regresa a al proceso que lo llamó

throw_load<-function(){
  # Lee el fichero con la matriz de recuentos
  cat("Introduce la ubicación de los ficheros de trabajo\n")
  file_counts<-readline(prompt = "Matriz de recuentos (.csv, .xlsx):" )
  cd<-load_data(file_counts)
  
  # prepara los datos de recuento para el análisis
  counts<-prepare_counts(cd)
  
  # Lee el fichero con las condiones experimentales por muestra
  file_metadata<-readline(prompt = "Metadatos (.csv, .xlsx):" )
  md<-load_data(file_metadata)
  
  # prepara los metadatos para el análisis
  metadata<-prepare_metadata(md)
  
  # Cambia los nombres de las filas en metadata y columnas en
  # counts
  rownames(metadata)<-metadata$iden
  colnames(counts)<-rownames(metadata)
  
  #comprueba que los nombres de las columnas en la data.frame de recuentos
  # es igual al nombre de las filas en metadata
  all.equal(colnames(counts),rownames(metadata))
  
  # convierte el data.frame de recuentos en matriz para posterior
  # procesos
  counts<-as.matrix(counts)
  
  # Crea variables globales para los recuentos y los metadatos
  assign("counts",counts, envir = .GlobalEnv)
  assign("counts_origin",counts, envir = .GlobalEnv)
  assign("metadata",metadata, envir = .GlobalEnv)
  
  cat("Preparando estadisticas de los datos...\n")
  get_table_muestras()
  
  cat("✅  Ficheros cargados con éxito.\n")

}



