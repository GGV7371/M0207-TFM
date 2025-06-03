# M0207-TFM
Repositorio M0207-TFM (Master Bioinformática, Bioestadística UOC)
Autora: Gemma Gariglio Viejo
Proyecto en R con el pipeline para el analisis de DEGs, análisis
funcional y correlación fenotípica correspondiente al trabajo fin
de máster.
Una vez cargado el proyecto en RStudio, para ejecutar el pipeline 
y acceder al menu de opciones hay que escribir en la consola
source("main.R"). La primera vez comenzará a cargar las librerias 
necesarias para la ejecución y luego mostrará el menú

=== Selecciona una opción ===
0: Salir
1: Cargar Ficheros de Análisis DEG
2: Exploración visual datos brutos
3: Análisis LOESS+lmdme
4: Análisis DESeq2
5: Cargar Ficheros para Enriquecimiento
6: Enriquecimiento funcional para lmdme
7: Enriquecimiento funcional para DESEq2
8: Comparación de enfoques
9: Análisis phenotípico
Opción:

Si se elige la opción 1: Cargar Ficheros de Análisis DEG o la opción 
5: Cargar Ficheros para Enriquecimiento pedirá que se introduzca la
ubicación de los ficheros de trabajo. Por ejemplo:

Introduce la ubicación de los ficheros de trabajo
Matriz de recuentos (.csv, .xlsx):

Se puede escribir la ruta o pulsar intro para que muestre el explorador 
de archivos para seleccionar el fichero.

En la carpeta /data contiene los ficheros con la matriz de recuentos, metadatos
y anotaciones que son la base de este trabajo. La carpeta /results contiene los ficheros 
que generan los distintos procesos que se ejecutan.
