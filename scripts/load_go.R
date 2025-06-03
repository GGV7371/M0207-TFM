#######################################################################
# Archivo: load_go.R
# Descripción: Código de la opción 4 del menú 
# Autr: Gemma Gariglio Viejo
# Fecha: 
# Dependencias: utils.R, procedrures.R, lmdme.R, deseq2.R
#######################################################################

source("scripts/lmdme.R")
source("scripts/DESEq2.R")
source("R/utils.R")
source("R/procedures.R")

throw_loadAnnon<-function(){

  cat("Introduce la ubicación de los ficheros de trabajo\n")
  file_genes<-readline(prompt = "Genes-Proteinas (.csv, .xlsx, .tsv):" )
  t_gen<-load_data(file_genes)
  
  file_annon<-readline(prompt = "Anotaciones Funcionales (.csv, .xlsx, .tsv):" )
  t_annon<-load_data(file_annon)
  
  cat("Generando anotaciones GO/KEGG...\n")
  annondb<-get_annon_db(t_gen, t_annon)
  
  assign("annondb",annondb, envir = .GlobalEnv)

}



