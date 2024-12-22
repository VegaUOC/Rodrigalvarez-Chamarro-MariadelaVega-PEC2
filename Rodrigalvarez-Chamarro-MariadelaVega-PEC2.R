## Variables
data_folder <- "data"
results_folder <- "results"

##### Carga de librerias #####

# Carga de librerías necesarias para la realización de la PEC

# Librería Bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Lectura de CELFiles
if (!require(oligo)){
  BiocManager::install(oligo)
}
library(oligo)

# Paquete de calidad
if (!require(arrayQualityMetrics)){
  BiocManager::install(arrayQualityMetrics)
}
library(arrayQualityMetrics)

# Paquete de anotación
if (!require(mouse4302.db)){
  BiocManager::install("mouse4302.db")
}
library("mouse4302.db")


# Filtrado y anotación de genes
if (!require(genefilter)){
  BiocManager::install(genefilter)
}
library(genefilter)

if (!require("limma")){
  BiocManager::install("limma")
}
library("limma")

if (!require(annotate)){
  BiocManager::install(annotate)
}
library(annotate)

if (!require(biomaRt)){
  BiocManager::install(biomaRt)
}
library(biomaRt)

BiocManager::install(c("clusterProfiler", "org.Mm.eg.db"))
library("clusterProfiler")
library(org.Mm.eg.db)

# Directorio de trabajo
getwd()

##### Carga y Preparación de los datos #####

# Preparación de los datos. Filtrado de los datos

# Función para realizar el filtrado de los datos
filter_microarray <- function(allTargets, seed = 123) {
  # Configurar la semilla aleatoria
  print(seed)
  set.seed(seed)
  
  # Filtrar las filas donde 'time' no sea 'hour 2'
  filtered <- subset(allTargets, time != "hour 2")
  
  # Dividir el dataset por grupos únicos de 'infection' + 'agent'
  filtered$group <- interaction(filtered$infection, filtered$agent)
  
  # Seleccionar 4 muestras al azar de cada grupo
  selected <- do.call(rbind, lapply(split(filtered, filtered$group), function(group_data) {
    if (nrow(group_data) > 4) {
      group_data[sample(1:nrow(group_data), 4), ]
    } else {
      group_data
    }
  }))
  
  # Obtener los índices originales como nombres de las filas seleccionadas
  original_indices <- match(selected$sample, allTargets$sample)
  
  # Modificar los rownames usando 'sample' y los índices originales
  rownames(selected) <- paste0(selected$sample, ".", original_indices)
  
  # Eliminar la columna 'group' y devolver el resultado
  selected$group <- NULL
  return(selected)
}

# Simular el dataset basado en la descripción proporcionada
sample = c("GSM944831", "GSM944838", "GSM944845", "GSM944852", "GSM944859",
           "GSM944833", "GSM944840", "GSM944847", "GSM944854", "GSM944861",
           "GSM944834", "GSM944841", "GSM944848", "GSM944855", "GSM944862",
           "GSM944832", "GSM944839", "GSM944846", "GSM944853", "GSM944860",
           "GSM944835", "GSM944842", "GSM944849", "GSM944856", "GSM944863",
           "GSM944836", "GSM944843", "GSM944850", "GSM944857", "GSM944864",
           "GSM944837", "GSM944844", "GSM944851", "GSM944858", "GSM944865")

allTargets <- data.frame(
  sample = sample,
  infection = c(rep("uninfected", 15), rep("S. aureus USA300", 20)),
  time = c(rep("hour 0", 15), rep("hour 2", 5), rep("hour 24", 15)),
  agent = c(rep("untreated", 5), rep("linezolid", 5), rep("vancomycin", 5),
            rep("untreated", 5), rep("untreated", 5), rep("linezolid", 5), rep("vancomycin", 5)),
  filename = paste0(sample, ".CEL"),
  color = c(rep("lightblue", 5), rep("pink", 5), rep("lightgreen", 5),
            rep("grey", 5), rep("blue", 5), rep("red", 5), rep("green", 5))
)

# Aplicar la función (cambiar 123 por vuestro ID de la UOC u otro número que podáis escribir en el documento)
result <- filter_microarray(allTargets, seed=25178976)

result$grupo <- c(rep("INF.LIN",4),rep("NOINF.LIN",4),rep("INF.NOTTREAT",4),
                  rep("NOINF.NOTTREAT",4),rep("INF.VAN",4),rep("NOINF.VAN",4))
result

# Se crea el objeto AnnotatedDataFrame
targets <- AnnotatedDataFrame(result)

# Cargar todos los ficheros uno a uno con extensión CEL
celFiles <- result$filename
rawData <- read.celfiles(file.path(params$data_folder,celFiles),phenoData=targets)


rawData

##### Análisis exploratorio y control de calidad #####

#```{r fig.cap="Intensidad distribución muestras \\label{fig:g1}",fig.align='center', fig.width=10, fig.height=5}
#Diagrama de cajas
sampleNames <- as.character(result$grupo)
sampleColor <- as.character(result$color)

boxplot(rawData, which="all",las=2, main="Distribución intensidad datos en crudo", 
        cex.axis=0.6, col=sampleColor, names=sampleNames)
legend("topright", legend = c("US3000 - linezolid", "No infectado - linezolid", 
                              "US3000 - sin tratamiento","No infectado - Sin Tratamiento",
                              "US3000 - vancomycin","No infectado - vancomycin"),
       fill = c("red", "pink", "blue","lightblue","green","lightgreen"), title = "Grupos")
  #```

# Clustering jerarquico
clust.euclid.average <- hclust(dist(t(exprs(rawData))),method="average")
plot(clust.euclid.average, labels<-sampleNames,main="Clustering jerarquizado datos en crudo",cex=0.7,hang=-1 )

# Análisis de componentes principales

# Función para realizar una representacion de PCA

plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE, 
                      formapunts=NULL, myCex=0.8,...)
{
  pcX<-prcomp(t(X), scale=scale) # o prcomp(t(X))
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors, pch=formapunts, 
       xlim=c(min(pcX$x[,1])-100000, max(pcX$x[,1])+100000),ylim=c(min(pcX$x[,2])-100000, max(pcX$x[,2])+100000))
  text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=myCex)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
  legend("bottomright", legend = c("US3000 - linezolid", "No infectado - linezolid", 
                                "US3000 - sin tratamiento","No infectado - Sin Tratamiento",
                                "US3000 - vancomycin","No infectado - vancomycin"),
         col = c("red", "pink", "blue","lightblue","green","lightgreen"), 
         pch = c(15,22,16,21,17,24), title = "Grupos")
}

plotPCA(exprs(rawData), labels=sampleNames, dataDesc="raw data", colors=sampleColor,
        formapunts=c(rep(15,4),rep(22,4),rep(16,4),rep(21,4),rep(17,4),rep(24,4)), 
        myCex=0.5)

# Ejecutar el paquete para estudiar la calidad de los datos.
rerun <-FALSE
if(rerun){
  arrayQualityMetrics(rawData,  reporttitle="Anexo I. Calidad de los datos", force=TRUE)
}

# Normalización de los datos
eset<-rma(rawData)
#params$results_folder
write.exprs(eset, file.path(results_folder, "NormData.txt"))

##### Filtrado de los datos #####

annotation(eset)<-"mouse4302.db"
eset_filtered <- nsFilter(eset, var.func=IQR,
                          var.cutoff=0.9, var.filter=TRUE, require.entrez=FALSE,
                          filterByQuantile=TRUE)
#NUMBER OF GENES REMOVED
print(eset_filtered)

#NUMBER OF GENES IN
print(eset_filtered$eset)
dim(eset_filtered$eset)

# Matriz filtrada
filteredEset <- eset_filtered$eset
filteredData <- exprs(filteredEset)
colnames(filteredData) <- rownames(pData(eset_filtered$eset))

##### Matriz de diseño y contraste #####

result$grupo <- c(rep("INF.LIN",4),rep("NOINF.LIN",4),rep("INF.NOTTREAT",4),
                  rep("NOINF.NOTTREAT",4),rep("INF.VAN",4),rep("NOINF.VAN",4))
result

# Construcción de la matriz de diseño
lev <- factor(result$grupo, levels = unique(result$grupo))
design <-model.matrix(~0+lev,pData(filteredEset))
colnames(design) <- levels(lev)
rownames(design) <- sampleNames

design

# Construcción de la matriz de contrastes
mCont <- makeContrasts(LIN_INF_NOINF = INF.LIN - NOINF.LIN, 
                       NOTTREAT_INF_NOINF = INF.NOTTREAT - NOINF.NOTTREAT,
                       VAN_INF_NOINF = INF.VAN - NOINF.VAN,
                       levels = design)


##### Obtención de las listas de genes diferencialmente expresados #####

# Construcción del modelo lineal

fit <- lmFit(filteredEset,design)
# Aplicar los contrastes
fit.main <- contrasts.fit(fit,mCont)
fit.main <- eBayes(fit.main)
class(fit.main)

# Función para generar los mapas de calor
plotHeatmap <- function ( X, title="")
{
  selectedData <- filteredData[X,]
  
  #HEATMAP PLOT
  my_palette <- colorRampPalette(c("yellow", "red"))(n = 299)
  library(gplots)
  heatmap.2(selectedData,
            Rowv=TRUE,
            Colv=TRUE,
            main=title,
            scale="row",
            col=my_palette,
            sepcolor="white",
            sepwidth=c(0.05,0.05),
            cexRow=0.5,
            cexCol=0.9,
            key=TRUE,
            keysize=1.5,
            density.info="histogram",
            ColSideColors=sampleColor,
            tracecol=NULL,
            srtCol=30)
  
}

# Tratamiento Linezolid - diferencia entre Infectados y no infectados
topTabLIN <-  topTable (fit.main, number=nrow(fit.main), coef="LIN_INF_NOINF", 
                      adjust="fdr", lfc=2, p.value=0.05)
dim(topTabLIN)
head(topTabLIN)

selectedRows <- rownames(filteredData) %in% rownames(topTabLIN)
title <- paste("HeatMap LINEZOLID", "Infectados vs No Infectados FC>=2", sep="\n")
plotHeatmap(selectedRows,title)

# Sin tratamiento - diferencia entre Infectados y no infectados
topTabNOTTREAT <-  topTable (fit.main, number=nrow(fit.main), coef="NOTTREAT_INF_NOINF", 
                        adjust="fdr", lfc=2, p.value=0.05)
dim(topTabNOTTREAT)
head(topTabNOTTREAT)

selectedRows <- rownames(filteredData) %in% rownames(topTabNOTTREAT)
title <- paste("HeatMap Sin Tratamiento", "Infectados vs No Infectados FC>=2", sep="\n")
plotHeatmap(selectedRows,title)

# Vancomycin - diferencia entre Infectados y no infectados
topTabVAN <-  topTable (fit.main, number=nrow(fit.main), coef="VAN_INF_NOINF", 
                             adjust="fdr", lfc=2, p.value=0.05)
dim(topTabVAN)
head(topTabVAN)

selectedRows <- rownames(filteredData) %in% rownames(topTabVAN)
title <- paste("HeatMap VANCOMYCIN", "Infectados vs No Infectados FC>=2", sep="\n")
plotHeatmap(selectedRows,title)


# Gráfico de volcan para visulizar los genes más significativos
volcanoplot(fit.main, highlight=40, names= rownames(topTabLIN), 
            main = "Genes expresado diferencialmente")
abline(v = c(-2, 2))


rtdo_contrates <- decideTests(fit.main, method = "separate", adjust.method = "fdr", 
                       p.value = 0.05, lfc = 2)
# LINZENOLID
genes_up_LIN <- rownames(filteredEset)[rtdo_contrates[, "LIN_INF_NOINF"] == 1]
genes_notSig_LIN <- rownames(filteredEset)[rtdo_contrates[, "LIN_INF_NOINF"] == 0]
genes_down_LIN <- rownames(filteredEset)[rtdo_contrates[, "LIN_INF_NOINF"] == -1]
up_LIN <- length(genes_up_LIN)
notSign_LIN <- length(genes_notSig_LIN)
down_LIN <- length(genes_down_LIN)

# SIN TRATAMIENTO
genes_up_NT <- rownames(filteredEset)[rtdo_contrates[, "NOTTREAT_INF_NOINF"] == 1]
genes_notSig_NT <- rownames(filteredEset)[rtdo_contrates[, "NOTTREAT_INF_NOINF"] == 0]
genes_down_NT <- rownames(filteredEset)[rtdo_contrates[, "NOTTREAT_INF_NOINF"] == -1]
up_NT <- length(genes_up_NT)
notSign_NT <- length(genes_notSig_NT)
down_NT <- length(genes_down_NT)

# VANCOMICINA
genes_up_VAN <- rownames(filteredEset)[rtdo_contrates[, "VAN_INF_NOINF"] == 1]
genes_notSig_VAN <- rownames(filteredEset)[rtdo_contrates[, "VAN_INF_NOINF"] == 0]
genes_down_VAN <- rownames(filteredEset)[rtdo_contrates[, "VAN_INF_NOINF"] == -1]
up_VAN <- length(genes_up_VAN)
notSign_VAN <- length(genes_notSig_VAN)
down_VAN <- length(genes_down_VAN)


# Genes comunes entre contrastes
common_UP_LIN_NT <- intersect(genes_up_LIN, genes_up_NT)
common_UP_LIN_VAN <- intersect(genes_up_LIN, genes_up_VAN)
common_UP_NT_VAN <- intersect(genes_up_NT, genes_up_VAN)
common_UP_LIN_NT_VAN <- intersect(common_UP_LIN_NT,common_UP_LIN_VAN )


# venn diagram para visualizar la superposición entre contrastes
vennDiagram(rtdo_contrates, include = "up")

vennTable <- vennCounts(rtdo_contrates,include="up")

##### Anotación de Genes ####

# Estudio de contraste Linezolid
topTabLIN$PROBEID <- rownames(topTabLIN)
myProbes_LIN <- topTabLIN$PROBEID
# Se realiza la anotación de los genes
geneAnots_LIN <- AnnotationDbi::select(mouse4302.db, keys=myProbes_LIN,
                    columns = c("SYMBOL","GENENAME","ENSEMBL","ENTREZID"))

# Se mezclan ambas tablas por PROBEID, para ordenarlas se mantiene el valor de B
annotatedTopTab_LIN <- merge(x=geneAnots_LIN, y=topTabLIN, by.x="PROBEID", by.y="PROBEID")
# Se ordenan las claves en orden descendente pro la columna B
sortAnnotatedTopTabLIn <- annotatedTopTab_LIN[order(-topTabLIN$B),
                                 c("PROBEID","SYMBOL","GENENAME","ENSEMBL","ENTREZID",
                                   "B")]
sortAnnotatedTopTabLIn <- sortAnnotatedTopTabLIn[!is.na(sortAnnotatedTopTabLIn$SYMBOL),]
dim(sortAnnotatedTopTabLIn)

head(sortAnnotatedTopTabLIn[,c("PROBEID","SYMBOL","GENENAME","ENSEMBL","ENTREZID")],10)

# Estudio de contraste Sin tratamiento
topTabNOTTREAT$PROBEID <- rownames(topTabNOTTREAT)
myProbes_NOTTREAT <- topTabNOTTREAT$PROBEID
# Se realiza la anotación de los genes
geneAnots_NOTTREAT <- AnnotationDbi::select(mouse4302.db, keys=myProbes_NOTTREAT,
                        columns = c("SYMBOL","GENENAME","ENSEMBL","ENTREZID"))

# Se mezclan ambas tablas por PROBEID, para ordenarlas se mantiene el valor de B
annotatedTopTab_NOTTREAT <- merge(x=geneAnots_NOTTREAT, y=topTabNOTTREAT, by.x="PROBEID", by.y="PROBEID")
# Se ordenan las claves en orden descendente pro la columna B
sortAnnotatedTopTabNOTTREAT <- annotatedTopTab_NOTTREAT[order(-topTabNOTTREAT$B),
                                              c("PROBEID","SYMBOL","GENENAME","ENSEMBL","ENTREZID",
                                                "B")]
sortAnnotatedTopTabNOTTREAT <- sortAnnotatedTopTabNOTTREAT[!is.na(sortAnnotatedTopTabNOTTREAT$SYMBOL),]
dim(sortAnnotatedTopTabNOTTREAT)

head(sortAnnotatedTopTabNOTTREAT[,c("PROBEID","SYMBOL","GENENAME","ENSEMBL","ENTREZID")],10)

# Estudio de contraste Vancomicina
topTabVAN$PROBEID <- rownames(topTabVAN)
myProbes_VAN <- topTabVAN$PROBEID
# Se realiza la anotación de los genes
geneAnots_VAN <- AnnotationDbi::select(mouse4302.db, keys=myProbes_VAN,
                             columns = c("SYMBOL","GENENAME","ENSEMBL","ENTREZID"))

# Se mezclan ambas tablas por PROBEID, para ordenarlas se mantiene el valor de B
annotatedTopTab_VAN <- merge(x=geneAnots_VAN, y=topTabVAN, by.x="PROBEID", by.y="PROBEID")
# Se ordenan las claves en orden descendente pro la columna B
sortAnnotatedTopTabVAN <- annotatedTopTab_VAN[order(-topTabVAN$B),
                                                        c("PROBEID","SYMBOL","GENENAME","ENSEMBL","ENTREZID",
                                                          "B")]
sortAnnotatedTopTabVAN <- sortAnnotatedTopTabVAN[!is.na(sortAnnotatedTopTabVAN$SYMBOL),]
dim(sortAnnotatedTopTabVAN)

head(sortAnnotatedTopTabVAN[,c("PROBEID","SYMBOL","GENENAME","ENSEMBL","ENTREZID")],10)

# Listado de genes comunes
# Se realiza la anotación de los genes
geneAnots_common <- AnnotationDbi::select(mouse4302.db, keys=common_UP_LIN_NT_VAN,
                                       columns = c("SYMBOL","GENENAME","ENSEMBL","ENTREZID"))

head(geneAnots_common[,c("PROBEID","SYMBOL","GENENAME","ENSEMBL","ENTREZID")],10)



###### Análisis de significancia biológica #####

# LINEZOLID
# Lista de genes significativos (símbolos)
genes_ENTREZ_LIN <- sortAnnotatedTopTabLIn$ENTREZID

ego_LIN <- enrichGO(
  gene = genes_ENTREZ_LIN,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",                # Ontología: BP (Biological Process), MF (Molecular Function), CC (Cellular Component)
  pAdjustMethod = "fdr",      # Método de ajuste para múltiples pruebas
  pvalueCutoff = 0.05,       # Umbral de significancia
  qvalueCutoff = 0.2         # Umbral de Q-valor
)

# Ver los resultados
head(ego_LIN)

barplot(ego_LIN, showCategory = 10, title = "Top 10 GO Terms")

dotplot(ego_LIN, showCategory = 10, title = "GO Enrichment Dot Plot")
ego_LIN_df <- as.data.frame(ego_LIN)
head(ego_LIN_df)

cnetplot(ego_LIN, categorySize = "geneNum", schowCategory = 15, vertex.label.cex = 0.75)

# SIN TRATAMIENTO
# Lista de genes significativos (símbolos)
genes_ENTREZ_NT <- sortAnnotatedTopTabNOTTREAT$ENTREZID

ego_NT <- enrichGO(
  gene = genes_ENTREZ_NT,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",                # Ontología: BP (Biological Process), MF (Molecular Function), CC (Cellular Component)
  pAdjustMethod = "fdr",      # Método de ajuste para múltiples pruebas
  pvalueCutoff = 0.05,       # Umbral de significancia
  qvalueCutoff = 0.2         # Umbral de Q-valor
)

# Ver los resultados
head(ego_NT)

barplot(ego_NT, showCategory = 10, title = "Top 10 GO Terms")

dotplot(ego_NT, showCategory = 10, title = "GO Enrichment Dot Plot")
ego_NT_df <- as.data.frame(ego_NT)
head(ego_NT_df)

cnetplot(ego_NT, categorySize = "geneNum", schowCategory = 15, vertex.label.cex = 0.75)

# VANCOMICINA
# Lista de genes significativos (símbolos)
genes_ENTREZ_VAN <- sortAnnotatedTopTabVAN$ENTREZID

ego_VAN <- enrichGO(
  gene = genes_ENTREZ_VAN,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",                # Ontología: BP (Biological Process), MF (Molecular Function), CC (Cellular Component)
  pAdjustMethod = "fdr",      # Método de ajuste para múltiples pruebas
  pvalueCutoff = 0.05,       # Umbral de significancia
  qvalueCutoff = 0.2         # Umbral de Q-valor
)

# Ver los resultados
head(ego_VAN)

barplot(ego_VAN, showCategory = 10, title = "Top 10 GO Terms")

dotplot(ego_VAN, showCategory = 10, title = "GO Enrichment Dot Plot")
ego_VAN_df <- as.data.frame(ego_VAN)
head(ego_VAN_df)

cnetplot(ego_VAN, categorySize = "geneNum", schowCategory = 15, vertex.label.cex = 0.75)

# Genes comunes
# Lista de genes significativos (símbolos)
genes_ENTREZ_COMMON <- geneAnots_common$ENTREZID

ego_COMMON <- enrichGO(
  gene = genes_ENTREZ_COMMON,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",                # Ontología: BP (Biological Process), MF (Molecular Function), CC (Cellular Component)
  pAdjustMethod = "fdr",      # Método de ajuste para múltiples pruebas
  pvalueCutoff = 0.05,       # Umbral de significancia
  qvalueCutoff = 0.2         # Umbral de Q-valor
)

# Ver los resultados
head(ego_COMMON)

barplot(ego_COMMON, showCategory = 10, title = "Top 10 GO Terms")

dotplot(ego_COMMON, showCategory = 10, title = "GO Enrichment Dot Plot")
ego_COMMON_df <- as.data.frame(ego_COMMON)
head(ego_COMMON_df)

cnetplot(ego_COMMON, categorySize = "geneNum", schowCategory = 15, vertex.label.cex = 0.75)