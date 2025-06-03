#######################################################################
# Archivo: main.R
# Descripción: Menú principal
# Autor: Gemma Gariglio Viejo
# Fecha: 02/06/2025
# Dependencias: utils.R, lmdmefun.R, procedrures.R, graphics.R
#######################################################################

source("scripts/load.R")
source("scripts/control_calidad.R")
source("scripts/lmdme.R")
source("scripts/deseq2.R")
source("scripts/load_go.R")
source("scripts/enricher_lmdme.R")
source("scripts/enricher_deseq2.R")
source("scripts/compara_enfoques.R")
source("scripts/pheno.R")

# Define el repositorio de paquetes
repository<-c('http://yihui.name/xran',
              'http://cran.rstudio.com')

# Define los paquetes de R para instalar
required_packages<-c("BiocManager",
                     "tcltk",
                     "readxl",
                     "dplyr",
                     "tidyverse",
                     "lmdme",
                     "limma",
                     "pheatmap",
                     "gplots",
                     "RColorBrewer",
                     "cluster",
                     "VennDiagram",
                     "ggrepel",
                     "flextable",
                     "gtsummary",
                     "officer",
                     "factoextra"
)

# Defines los paquetes de BioConducor para instalar
required_BiocManager<- c("DESeq2",
                         "clusterProfiler",
                         "GO.db"
                         )

# Llamada a la función "install_load_packages" function
install_load_packages(required_packages, repository)

# Llamada a la función "install_load_packBioC"
install_load_packBioC(required_BiocManager)

# Crea los directorios de trabajo
create_work_directory()

#===================================================
# Menú
#===================================================
# Opciones a mostrar en pantalla
repeat {
  cat("=== Selecciona una opción ===\n")
  cat("0: Salir\n")
  cat("1: Cargar Ficheros de Análisis DEG\n")
  cat("2: Exploración visual datos brutos\n")
  cat("3: Análisis LOESS+lmdme\n")
  cat("4: Análisis DESeq2\n")
  cat("5: Cargar Ficheros para Enriquecimiento\n")
  cat("6: Enriquecimiento funcional para lmdme\n")
  cat("7: Enriquecimiento funcional para DESEq2\n")
  cat("8: Comparación de enfoques\n")
  cat("9: Análisis phenotípico\n")
  
  # Lee la opción seleccionada
  opcion <- as.integer(readline(prompt = "Opción: "))
  
  # Acciones a lanzar
  if (opcion == 0) {
    cat("👋 Saliendo del programa...\n")
    break
  } else if (opcion == 1) {
    throw_load()
  } else if (opcion == 2) {
    throw_control_calidad()
  } else if (opcion == 3) {
    throw_lmdme()
  } else if (opcion == 4) {
    throw_deseq2()
  } else if (opcion == 5) {
    throw_loadAnnon()
  } else if (opcion == 6) {
    throw_enricher_lmdme()
  } else if (opcion == 7) {
    throw_enricher_deseq2()
  } else if (opcion == 8) {
    throw_compara()
  } else if (opcion == 9) {
    throw_pheno()
  } else {
    cat("❌ Opción no válida. Intentalo de nuevo.\n")
  }
}
