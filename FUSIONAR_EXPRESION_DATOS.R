#CREAR MATRIZ DE EXPRESION QUE CONTIENE LOS NIVELES DE EXPRESION DE CADA GEN PARA CADA MUESTRA 

# Leer metadatos descargados en TCGA (GDC)
metadata <- read.table("TCGA-BRCA_metadata.txt", sep="\t", header=TRUE)

# Folder de datos de expresion de proyecto BRCA descargados en TCGA (GDC)
folder_gdc <- "gdc_download_20230329_213220.612876"

# Procesar la primera fila 
x <- 1

# Leer los datos de expresion 
id_sample <- metadata$id[x]
name <- paste(folder_gdc, paste(id_sample, metadata$file_name[x],  sep="/" ), sep = "/")
base.data <- read.table(name, sep="\t", header=TRUE)

# Eliminar lecturas ambiguas y no alineadas 
base.data <- base.data[!grepl("N_", base.data$gene_id ),]

# Mantener solo los recuentos variados 
base.data <- base.data[, c("gene_id", "unstranded")]
names(base.data) <- c("gene_id", id_sample)

# Fusionar los archivos restantes 
for (x in 2:nrow(metadata)){
  id_sample <- metadata$id[x]
  name <- paste(folder_gdc, paste(id_sample, metadata$file_name[x],  sep="/" ), sep = "/")
  
  exp.data <- read.table(name, sep="\t", header=TRUE)
  
  exp.data <- exp.data[!grepl("N_", exp.data$gene_id ),]
  exp.data <- exp.data[, c("gene_id", "unstranded")]
  names(exp.data) <- c("gene_id", id_sample)
  
  merge.data <- merge(base.data, exp.data, by = "gene_id", all=TRUE)
  base.data <- merge.data
}

# Guardar la matriz en un archivo tipo tsv 
write.table(base.data, file="TCGA-BRCA.all.expression.tsv", row.names = FALSE, col.names = TRUE, sep=",")


