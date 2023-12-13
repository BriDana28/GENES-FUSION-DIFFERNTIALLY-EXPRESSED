#CREAR MATRIZ DE DISENO DE CLASIFICACION DE MUESTRAS "tumor" y "normal" PARA ASIGNARLA A UN OBJETO DGEList 
# CON EL FIN DE REALIZAR UN ANALISIS DE EXPRESIÓN 

# Cargar paquetes necesarios
library(edgeR)
library(ggplot2)
library(limma)  
library(statmod)
library(pheatmap)

# Cargar los metadatos y la matriz de diseno con los niveles de expresion relacionados al proyecto TCGA-BRCA
metadata <- read.delim("TCGA-BRCA_metadata.txt", sep = "\t")

base.data <- read.delim("TCGA-BRCA.all.expression.tsv",
                        sep = ",")

# Unificar muestras entre metadatos y datos de expresion
samples.expr <- colnames(base.data[, 2:ncol(base.data)])
samples.names <- gsub(".", "-", samples.expr, fixed = TRUE)
samples.names <- gsub("^X", "", samples.names)

gene.id <- base.data$gene_id
base.data <- base.data[, -1]
colnames(base.data) <- samples.names

metadata <- metadata[metadata$id %in% samples.names,]

# Eliminar muestras metastaticas
sample.metastatic <- metadata$id[metadata$cases.0.samples.0.sample_type == "Metastatic"]
base.data <- base.data[, !colnames(base.data) %in% sample.metastatic]
metadata <- metadata[metadata$cases.0.samples.0.sample_type != "Metastatic",]

# Crear matriz de diseno para clasificar las muestras en grupos ("tumor" y "normal")
group <- ifelse(metadata$cases.0.samples.0.sample_type == "Primary Tumor", "tumor", "normal")
group <- factor(group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Crear objeto DGEList a partir de la matriz de diseno
base.data.mat <- as.matrix(base.data)
dge <- DGEList(counts = base.data.mat, group = group)
rownames(dge) <- gene.id

# Filtrar los genes que no tienen suficiente expresión en al menos una muestra
keep <- filterByExpr(dge, design)
dge <- dge[keep,, keep.lib.sizes = FALSE]

# Normalizar los datos
dge <- calcNormFactors(dge, method = "TMM")

# Estimar la dispersión
dge <- estimateDisp(dge, design, robust = TRUE)
fit <- glmQLFit(dge, design, robust = TRUE)

# Seleccionar un contraste
contrast <- makeContrasts(tumor - normal, levels = design)
res.de <- glmQLFTest(fit, contrast = contrast)

# Obtener los genes diferencialmente expresados
tags <- topTags(res.de, n= 5000, p.value = 0.0001)
de.genes <- rownames(tags$table)[tags$table$FDR < 0.05 & abs(tags$table$logFC) > 2]


##########################################################################################################
#FILTRAR EL OBJETO DGEList PARA MANTENER SOLO LOS GENES DIFERENCIALMENTE EXPRESADOS 

# Filtrar el objeto DGEList 
dge.de_genes <- dge[de.genes, ]

# Normalizar su expresion
dge.norm <- cpm(dge.de_genes, normalized.lib.sizes = TRUE, log = TRUE)
dge.norm.vec <- as.vector(dge.norm)

n.samples <- ncol(dge.norm)
n.genes <- nrow(dge.norm)

genes.vec <- rep(de.genes, times = n.samples)

samples.vec <- rep(dge.de_genes$samples$group, each = n.genes)

##########################################################################################################
#CREAR MARCO DE DATOS PARA PODER EJECUTAR UN BOXPLOT PARA COMPARAR LAS MUESTRAS TUMORALES VS NORMALES POR CADA GEN

# Crear el marco de datos para el boxplot
df_gene_expr <- data.frame(Sample_Type = samples.vec, logCPM =dge.norm.vec, Gene_ID = genes.vec )

# Crear una carpeta para los boxplots
if (!file.exists("boxplots_per_gene3")) {
  dir.create("boxplots_per_gene3")
}

# Recorrer a traves de los genes y crear boxplots 
for (gene in de.genes) {
  gene_data <- df_gene_expr[df_gene_expr$Gene_ID == gene, ]
  
  boxplot <- ggplot(gene_data, aes(x = factor(Sample_Type, levels = c("normal", "tumor")), y = logCPM, fill = Sample_Type)) +
    geom_boxplot() +
    labs(title = paste("Boxplot para", gene), x = "Sample Type", y = "Log CPM") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Guardar el boxplot en la carpeta correspondiente al gen
  file_path <- file.path("boxplots_per_gene3", paste(gene, ".png", sep = ""))
  ggsave(file_path, plot = boxplot, width = 6, height = 4)
}

#########################################################################################################
# COMPARAR LOS GRUPOS "tumor" y "normal" EN CUANTO A SU EXPRESION DIFERENCIAL 

# Filtrar el marco de datos para obtener solo los datos de tumor y normal
tumor_data <- df_gene_expr[df_gene_expr$Sample_Type == "tumor", ]
normal_data <- df_gene_expr[df_gene_expr$Sample_Type == "normal", ]

# Crear una matriz de expresión para los datos de tumor y normal
expression_matrix <- matrix(
  c(tumor_data$logCPM, normal_data$logCPM),
  nrow = length(unique(df_gene_expr$Gene_ID)),
  byrow = FALSE
)

# Obtener nombres de genes para etiquetas de filas
gene_labels <- unique(df_gene_expr$Gene_ID)

# Crear un vector de colores para el heatmap
colors <- colorRampPalette(c("blue", "white", "red"))(50)

# Crear el heatmap usando la función pheatmap()
pheatmap(expression_matrix,
         color = colors,
         scale = "row",
         clustering_distance_rows = "euclidean",  # Puedes ajustar el método de distancia
         clustering_method = "ward.D",  # Puedes ajustar el método de clustering
         labels_row = gene_labels,  # Etiquetas de filas
         labels_col = c(rep("Tumor", ncol(tumor_data)), rep("Normal", ncol(normal_data))),  # Etiquetas de columnas
         main = "Heatmap Comparando Tumor vs Normal")


