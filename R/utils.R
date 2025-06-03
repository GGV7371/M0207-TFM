############################################################
# Archivo: utils.R
# Descripción: Funciones auxiliares del proyecto
# Autor: Gemma Gariglio Viejo
# Fecha: 02/06/2025
############################################################

#===========================================================
# Función: install_load_packages
# Descripción: Instala y/o carga paquetes de R
#===========================================================

#' Instala y/o carga paquetes de R
#' comprobando si los paquetes ya están cargados y si no 
#' lo están, los instala y luego los carga
#' 
#' @param packages Vector de nombres de paquetes
#' @param repos URL del repositorio CRAN para la instalación
#' @return No devuelve nada

install_load_packages<-function(packages, 
                                repos = "https://cloud.r-project.org"){
  
  for (pkg in packages){
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg, type = 'source', repos = repos)
      library(pkg, character.only = TRUE)
    } else {
      library(pkg, character.only = TRUE)
    }
  }
}

#===========================================================
# Función: install_load_packBioC
# Descripción: Instala y/o carga paquetes de Bioconductor
#===========================================================

#' Instala y/o carga paquetes de Bioconductor
#' comprobando si ya están cargados y si no lo están
#' los instala
#' 
#' @param biocPackages Vector de nombres de paquetes
#' @return No devuelve nada

install_load_packBioC <- function (bioCPackages){
  for (pkg in bioCPackages){
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      BiocManager::install(pkg)
      library(pkg, character.only = TRUE)
    } else {
      library(pkg, character.only = TRUE)
    }  
  }
}

#===========================================================
# Función: create_work_directory
# Descripción: Crea los directorios de trabajo
#===========================================================

#' Crea el directorio de resultados y sus subdirectorios
#'  
#' sin parámetros
#' @return No devuelve nada

create_work_directory<-function(){
  trabajoDir<-getwd()
  
  directories <- c("results/")
  
  #file.remove(
  #  # Create a character vector of relative paths
  #  # to all files in the variable directories
  #  list.files(path = directories,
  #             all.files = TRUE,
  #             full.names = TRUE,
  #             recursive = TRUE)
  #)
  
  # Si no existe el directorio results lo crea
  if (!(dir.exists("./results"))) {
    dir.create("results")
  }
  
  # se posiciona en el directorio results
  setwd("./results")
  
  # Define los nombre de los subdirectorio
  directories <- c("datafiles", "intermediateData", "images")
  
  # Crea los subdirectorios dentro de results
  lapply(directories, function(x){
    if (!(dir.exists(x))){
      dir.create(x)
    }
  })
  
  # vuelve al directorio raiz
  setwd("..")
  
}

#===========================================================
# Función: load_data
# Descripción: Carga un archivo de datos en formato CSV, 
# Excel (.xlsx) o tsv
#===========================================================

#' Carga un archivo de datos en formato CSV, Excel (.xlsx) o tsv
#' 
#' Se puede escribir ruta/nombre del fichero o bien si se pulsa enter, 
#' abre el explorador de archivos en el diretorio de trabajo para seleccionar uno.
#' Valida que el archivo exista y que el formato sea el esperado.
#' 
#' @param archivo Ruta al archivo
#' @return data.frame o tibble con los datos que contiene el archivo

load_data <- function(archivo){  
  if (archivo == ""){
    archivo<-rstudioapi::selectFile()
  }
  
  if (!file.exists(archivo)) {
    stop("❌ El archivo no existe. Revisa el nombre o la ruta.")
  }
  
  # Lectura de datos en función de la extensión del archivo
  if (grepl("\\.csv$", archivo, ignore.case = TRUE)) {
    datos <- read.csv(archivo)
  } else if (grepl("\\.xlsx?$", archivo, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("📦 El paquete 'readxl' es necesario para leer archivos Excel.")
    }
    datos <- readxl::read_excel(archivo)
  } else if (grepl("\\.tsv?$", archivo, ignore.case = TRUE)) {
      datos <- read.table(archivo, sep = "\t", header = TRUE)
  }else{
    stop("Formato no soportado. Usa .csv, .xlsx o .tsv")
  }
  return(datos)
}

#===========================================================
# Función: get_formula_modelo
# Descripción: Modela la formula para emplear en modelos 
# estadísticos
#===========================================================

#' Modela la formula para emplear en modelos estadísticos
#' 
#' @param vars vector con nombres de variables (character vector)
#' @param operador espera un "+", ":" o "*" (string)
#' @return formula para usar en modelos estadisticos
#' 
#' @examples 
#' get_formula_deseq(C("treatment","stage"), "*")
#' # devuelve ~ treatment * stage

get_formula_modelo <- function(vars, operador) {
  
  if (length(vars) == 0) {
    stop("Debes proporcionar al menos una variable.")
  }
  
  if (!operador %in% c("+", ":", "*")) {
    stop("Operador inválido. Usa '+', ':', o '*'.")
  }
  
  # Pegar las variables con "*" (interacción completa)
  formula_text <- paste("~", paste(vars, collapse = operador))
  
  # Convertir el texto en fórmula
  return(as.formula(formula_text))
}

#===========================================================
# Función: get_colors
# Descripción: Asigna colores de la paleta especificada 
# al vector que se le pase como parámetro
#===========================================================

#' Asigna colores de la paleta especificada al vector 
#' que se le pase como parámetro
#' 
#' @param cl vector al que asignar colores
#' @param paleta nombre de la paleta de colores (string)
#' @return vector de colores

get_colors<-function(cl){
  col_cl <- setNames(brewer.pal(n = length(cl), 
                                name = "Set2"),
                     cl)
  return(col_cl)
}

#===========================================================
# Función: get_annon_db
# Descripción: Crea dos tabla de anotaciones una para 
# los términos GO y otra para las rutas metabolicas KEGG
#===========================================================

#' Crea dos tabla de anotaciones una para los términos
#' GO y otra las rutas metabólicas KEGG
#' 
#' @param t_gen data.frame que asocia identificadores de genes
#' con identificadores de proteínas
#' @param t_annon data.frame que asocia identificadores de 
#' proteinas con términos GO, KEGG
#' @return lista con las dos tablas de anotaciones

get_annon_db<-function(t_gen, t_annon){
  
  # une las tablas 
  genes_annon <- t_annon %>%
    inner_join(t_gen, by = c("query" ="Protein.accession"))

  cat("Creando tabla de anotaciones GO\n")
  go_db <- genes_annon %>%
    dplyr::select(Gene.ID, GOs) %>%
    rename(Gene_ID = Gene.ID,
           GO_ID = GOs) %>%
    mutate(Gene_ID = as.character(Gene_ID),
           GO_ID = str_split(GO_ID, ",")) %>%
    unnest(GO_ID) %>%
    distinct()
    #mutate(across(everything(), ~na_if(., "-")))
  
  # Busca las descripciones de los términos GO
  go_desc <- AnnotationDbi::select(GO.db,
                                   keys = unique(go_db$GO_ID),
                                   columns = "TERM",
                                   keytype = "GOID")
  
  # Une la tabla de identificadores GO con la de descripciones GO
  go_db <- left_join(go_db, go_desc, by = c("GO_ID" = "GOID"))
  
  
  cat("Creando tabla de anotaciones KEGG\n")
  kegg_db <- genes_annon %>%
    dplyr::select(Gene.ID, KEGG_Pathway) %>%
    rename(Gene_ID = Gene.ID,
           KEGG_ID = KEGG_Pathway) %>%
    mutate(Gene_ID = as.character(Gene_ID),
           KEGG_ID = str_split(KEGG_ID, ",")) %>%
    unnest(KEGG_ID) %>%
    distinct()
  
  library(KEGGREST)
  
  # Obtener lista completa de pathways KEGG
  pathways <- keggList("pathway")
  
  # Convertir a data.frame
  term2name <- data.frame(
    KEGG_ID = sub("path:", "", names(pathways)),
    NAME = as.character(pathways),
    stringsAsFactors = FALSE
  )
  
 kegg_db<-kegg_db %>%
   left_join(term2name, by = "KEGG_ID")
 
 
  return(list(go_db=go_db, kegg_db=kegg_db))
}

#===========================================================
# Función: show_table_word
# Descripción: Da formato word a una tabla y lo guarda 
# en un fichero
#===========================================================

#' Da formato word a una tabla y lo guarda en un fichero tipo word
#' 
#' @param ft data.frame o tibble
#' @param m string con el tipo de enfoque
#' @param title string con el titulo de la tabla
#' @param fichero string con el nombre del fichero
#' @return no devuelve nada
 
show_table_word<-function(ft, m, title, fichero) {
  
  # Formatea la tabla
  ft <- flextable(ft) %>%
    font(fontname = "Calibri") %>%
    align(align = "center", part = "all") %>%
    fontsize(size =10, part = "all") %>%
    bold(part = "header") %>%
    border_remove() %>%
    border_outer(border = fp_border(color = "cyan", width = 2), part = "header") %>%
    border_outer(border = fp_border(color = "cyan", width = 0.5), part = "body") %>%
    border_inner(border = fp_border(color = "cyan", width = 0.5), part = "all") %>%
    autofit()
  
  # Exportar a Word
  doc <- read_docx() %>%
    body_add_par(title, 
                 style = "heading 1") %>%
    body_add_flextable(ft)
  
  print(doc, target = paste0("results/intermediateData/", fichero, "_", m,".docx"))
}