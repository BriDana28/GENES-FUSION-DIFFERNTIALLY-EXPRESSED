# CREAR MATRIZ CON LOS DATOS DE LOS GENES FUSIONADOS Y EDITAR LOS DATOS PARA PODER REALIZAR UN ANALISIS DE
# ENREQUECIMIENTO PARA PODER ANALIZAR SI LAS VARIABLES SON INDEPENDIENTES O DEPENDIENTES 

# Cargar los paquetes necesarios
library(VennDiagram)
library(readxl)
library(openxlsx)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(data.table)
library(dlpyr)

# Cargar base de datos de fusion del TumorFusionPortal
data <- read_excel("C:/Users/brida/Documents/SERVICIO SOCIAL/nar-02671-data-e-2017-File007.xlsx") 

#Filtro de tumores BRCA en el documento de TumorFusionPortal
tumores_BRCA <- subset(data, `Tissue` == "BRCA")

# Calcular la columna sum_reads
tumores_BRCA$sum_reads <- tumores_BRCA$`Discordant Read Pairs` + tumores_BRCA$`Junction Spanning Reads`

# Filtro de tumores BRCA en dos categorías 
tumores_BRCA_10 <- subset(tumores_BRCA, sum_reads > 10)
tumores_BRCA_50 <- subset(tumores_BRCA, sum_reads > 50)

# Crear histograma 
hist(tumores_BRCA$sum_reads, xlim = c(1, 500), breaks = 500)

tumores_BRCA_seleccionados <- tumores_BRCA[, c("Sample","Gene_A", "Gene_B")]

# Crear matriz de genes fusinados en el contexto BRCA 
matriz_labjack <- as.matrix(tumores_BRCA_seleccionados)

##########################################################################################################################
# CAMBIO DE NOMBRES DE ID DE GENES A FORMATO ENSEMBL EN LA MATRIZ DE GENES FUSIONADOS 
# Aplicar funcion "unique" a las columnas "Gene_A" y "Gene_B"
nombres_genes_originales <- unique(c(tumores_BRCA_seleccionados$Gene_A, tumores_BRCA_seleccionados$Gene_B))

# Crear una lista de genes fusinados para cambiarlos a ENSEMBL
genes_fusionados <- list()

# Crear un objeto para la base de datos org.Hs.eg.db
orgdb <- org.Hs.eg.db

# Mapear los nombres en la base de datos org.Hs.eg.db de genes originales a identificadores Ensembl 
for (columna in c("Gene_A", "Gene_B")) {
  genes_originales <- unique(tumores_BRCA[[columna]])
  mapeo <- select(orgdb, keys = genes_originales, keytype = "SYMBOL", columns = "ENSEMBL")
  genes_fusionados[[columna]] <- mapeo$ENSEMBL[match(genes_originales, mapeo$SYMBOL)]
}

# Reemplazar los nombres de genes originales en la matriz con los identificadores Ensembl
matriz_labjack <- as.matrix(tumores_BRCA[, c("Sample", "Gene_A", "Gene_B")])
matriz_labjack[, "Gene_A"] <- genes_fusionados[["Gene_A"]][match(matriz_labjack[, "Gene_A"], unique(tumores_BRCA$Gene_A))]
matriz_labjack[, "Gene_B"] <- genes_fusionados[["Gene_B"]][match(matriz_labjack[, "Gene_B"], unique(tumores_BRCA$Gene_B))]

# Función para extraer el identificador en formato "ENSG (ensable)" a partir de identificadores actuales
extraer_identificador_ENS <- function(identificador_actual) {
  identificador_ENS <- sub(".*(ENSG\\d+\\.\\d+).*", "\\1", identificador_actual)
  return(identificador_ENS)
}

# Aplicar la función para obtener los identificadores en formato Ensembl y agregarlo a las columnas 
matriz_labjack[, "Gene_A"] <- sapply(matriz_labjack[, "Gene_A"], extraer_identificador_ENS)
matriz_labjack[, "Gene_B"] <- sapply(matriz_labjack[, "Gene_B"], extraer_identificador_ENS)

# Imprimir las primeras filas de matriz_labjack con los genes en formato Ensembl
head(matriz_labjack)
#########################################################################################################
#CREAR VECTORES PARA PODER REALIZAR UNA COMPARARCION DE GENES FUSIOANDOS Y GENES DIFERENCIALMENTE EXPRESADOS 

# Obtener los genes de fusión de ambas columnas "Gene_A" y "Gene_B" en matriz_labjack
genes_fusion <- union(matriz_labjack[, "Gene_A"], matriz_labjack[, "Gene_B"])

# Crear un vector con los nombres de los genes diferencialmente expresados con cambios sustanciales
genes_dif_exp_cambios_sustanciales <- de.genes

# Crear un vector con todos los genes utilizados en el análisis de expresión diferencial
todos_los_genes <- rownames(dge)

##########################################################################################
#ELIMINAR EL NUMERO DE VESION DE CADA GEN PARA SOLO TENER SU IDENTIFICADOR ENSEMBL 

# Eliminar la versión y conservar solo el identificador Ensembl principal
genes_dif_exp_cambios_sustanciales_sinversion <- sub("\\.\\d+", "", genes_dif_exp_cambios_sustanciales)

# Ver los identificadores resultantes
print(genes_dif_exp_cambios_sustanciales_sinversion)

# Eliminar la versión y conservar solo el identificador Ensembl principal
todos_los_genes_sin_version <- sub("\\.\\d+", "", todos_los_genes)


#######################################################################################
#CREAR GRAFICO DE VEEN PARA ANALIZAR LA INTERSECCION DE GENES DIFEERNCIALMENTE EXPRESADOS Y GENES DE FUSION EN TODO EL UNIVERSO DE GENES PREVIAMENTE OBTENIDO

#Función Unique
# Aplicar unique a tus vectores de nombres de genes para eliminar cualquier repetecion o inconsistencia en los datos 
genes_fusion <- unique(genes_fusion)
genes_dif_exp_cambios_sustanciales_sinversion <- unique(genes_dif_exp_cambios_sustanciales_sinversion)
todos_los_genes_sin_version <- unique(todos_los_genes_sin_version)

# Encontrar la intersección entre los genes de fusión y los genes diferencialmente expresados con cambios sustanciales
genes_interseccion <- intersect(genes_fusion, genes_dif_exp_cambios_sustanciales_sinversion)

# Encontrar los genes de fusión que no están en la intersección
genes_fusion_no_interseccion <- setdiff(genes_fusion, genes_interseccion)

# Encontrar los genes diferencialmente expresados con cambios sustanciales que no están en la intersección
genes_dif_exp_cambios_sustanciales_no_interseccion <- setdiff(genes_dif_exp_cambios_sustanciales_sinversion, genes_interseccion)

# Encontrar los genes que no están en ninguna de las dos categorías
genes_fuera <- setdiff(todos_los_genes, union(genes_fusion, genes_dif_exp_cambios_sustanciales_sinversion))

# Eliminar valores NA de la lista genes_fusion
genes_fusion <- genes_fusion[!is.na(genes_fusion)]


# Crear el diagrama de Venn
venn.plot <- venn.diagram(
  x = list(
    Todos = todos_los_genes_sin_version,
    Fusión = genes_fusion,
    Dif_Exp = genes_dif_exp_cambios_sustanciales_sinversion
  ),
  category.names = c("Todos", "Fusión", "Dif_Exp"),
  filename = NULL,
  output = TRUE
)

# Mostrar el diagrama de Venn
grid.draw(venn.plot)

# Crear una interseccion en los genes de fusion con el universo de genes ya que se identificaron 117 genes que no se sabe si estan expresados 
genes_fusion_comunes <- intersect(genes_fusion, todos_los_genes_sin_version)

# Crear el diagrama de Venn
venn.plot <- venn.diagram(
  x = list(
    Todos = todos_los_genes_sin_version,
    Fusión = genes_fusion_comunes,
    Dif_Exp = genes_dif_exp_cambios_sustanciales_sinversion
  ),
  category.names = c("Todos", "Fusión", "Dif_Exp"),
  filename = NULL,
  output = TRUE
)

# Mostrar el diagrama de Venn
grid.draw(venn.plot)
##########################################################################################################################################
# REALIZAR MATRIZ DE CONTIGENCIA PARA PODER APLICAR UNA PRUEBA FISHER Y CHI-Cuadrada 

# Obtener los nombres de los genes en la intersección
genes_interseccion <- intersect(genes_fusion_comunes, genes_dif_exp_cambios_sustanciales_sinversion)

# Calcular los tamaños de las intersecciones y conjuntos
int <- length(genes_interseccion)
A <- length(setdiff(genes_fusion_comunes, genes_dif_exp_cambios_sustanciales_sinversion))
B <- length(setdiff(genes_dif_exp_cambios_sustanciales_sinversion, genes_fusion_comunes))
C <- length(setdiff(todos_los_genes_sin_version, union(genes_fusion_comunes, genes_dif_exp_cambios_sustanciales_sinversion)))

# Crear y transponer la matriz 2x2
matriz_2x2_prueba <- matrix(c(int, A, B, C), nrow = 2, byrow = TRUE, dimnames = list(c("SI.DE", "NO.DE"), c("SI.FUS", "NO.FUS")))

# Imprimir la matriz
print(matriz_2x2_prueba)

#Realizar prueba Fisher y ji-cuadrada 
resultado_fisher <- fisher.test(matriz_2x2_prueba)
resultado_jicuadrada <- chisq.test(matriz_2x2_prueba)
# Realizar prueba de proporciones
prop_prueba <- prop.test(matriz_2x2_prueba)

######################################################################################################################################################
#OBTENER TABLA DE LOS NUMEROS MARGINALES 

# Calcula los totales marginales
totals <- margin.table(matriz_2x2_prueba)

# Imprime los totales
print(totals)

# Calcular los totales de filas y columnas
total_filas <- rowSums(matriz_2x2_prueba)
total_columnas <- colSums(matriz_2x2_prueba)

# Calcular el Universo
universo <- totals

# Crear una nueva matriz con Total y Universo agregados
nueva_matriz <- rbind(matriz_2x2_prueba, total_columnas)
nueva_matriz <- rbind(nueva_matriz, c(universo, NA))

# Agregar los nombres de fila y columna
row.names(nueva_matriz) <- c("SI.DE", "NO.DE", "Total", "Universo")
colnames(nueva_matriz) <- c("SI.FUS", "NO.FUS")

# Imprimir la nueva matriz
print(nueva_matriz)


# Calcular las proporciones 
proporcion_SI_FUS_SI_DE <- (nueva_matriz[1, 1] / universo) * (nueva_matriz[3, 1] / universo) * universo
proporcion_NO_FUS_SI_DE <- (nueva_matriz[1, 2] / universo) * (nueva_matriz[3, 2] / universo) * universo
proporcion_SI_FUS_NO_DE <- (nueva_matriz[2, 1] / universo) * (nueva_matriz[3, 1] / universo) * universo
proporcion_NO_FUS_NO_DE <- (nueva_matriz[2, 2] / universo) * (nueva_matriz[3, 2] / universo) * universo

# Crear una nueva matriz con los cálculos
nueva_matriz_calculos <- matrix(c(proporcion_SI_FUS_SI_DE, proporcion_NO_FUS_SI_DE, proporcion_SI_FUS_NO_DE, proporcion_NO_FUS_NO_DE), nrow = 2, byrow = TRUE, dimnames = list(c("SI.DE", "NO.DE"), c("SI.FUS", "NO.FUS")))

# Imprimir la nueva matriz
print(nueva_matriz_calculos)











