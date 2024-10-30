#Script para armar una red de coexpresion con WGCNA y obtener correlaciones con rasgos fenotípicos
#instalar las librerias que hagan falta
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("WGCNA")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("lattice")
#cargar las librerias necesarias
library(WGCNA)
library(ggplot2)
library(dplyr)
library(lattice)
#Establecer estos parametros es necesario para el funcionamiento del paquete
options(stringsAsFactors = FALSE);
enableWGCNAThreads()
#IMPORTANTE: establecer el directorio donde se trabajara y donde se encuentra el archivo de entrada
getwd()
dir = "C:/Users/micha/OneDrive - Universidad de los andes/9no Semestre Final/Algoritmos en Biologia Computacional/Tarea 5/CoexpresiónYMedidasDeCentralidad/data" #cambiar a la carpeta donde descargaron los archivos de entrada.
setwd(dir)
#nombre de la matriz de expresion
file = "LiverFemale3600.csv"
#Leer la matriz de expresion y dejar solo los resultados importantes
df = read.csv(file, header = TRUE)
df_expr = as.data.frame(t(df[, -c(1:8)]))
names(df_expr) = df$substanceBXH;
rownames(df_expr) = names(df)[-c(1:8)]
#filtrar los genes y las muestras que tengan muy pocos valores
gsg = goodSamplesGenes(df_expr, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  df_expr = df_expr[gsg$goodSamples, gsg$goodGenes]
}

#Deteccion de muestras outliers y retirarlas en caso de que existan
sampleTree = hclust(dist(df_expr), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
abline(h = 15, col = "red")

#PUNTO I: cutHeight es el parametro clave, se debe observar en el dendrograma cual es el punto de corte 
# de los outliers de acuerdo a la altura en el arbol, modificarlo de acuerdo al umbral
cutHeight = 15 #cambiar de acuerdo a la gráfica
clust = cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 10)
table(clust)
keepSamples = (clust==1)
df_expr_filtrado = df_expr[keepSamples, ]
nGenes = ncol(df_expr_filtrado)
nSamples = nrow(df_expr_filtrado)

#Cargar los datos de fenotipo
fenotipos = read.csv("ClinicalTraits.csv")
dim(fenotipos)
names(fenotipos)

#Se remueven las columnas que contienen informacion innecesaria
fenotipos_filtrado = fenotipos[, -c(31, 16)]
fenotipos_filtrado = fenotipos_filtrado[, c(2, 11:36) ]
dim(fenotipos_filtrado)
names(fenotipos_filtrado)

#Se debe organizar en un dataframe que sea analogo al de expresion
muestrasFemeninas = rownames(df_expr_filtrado)
filasFenotipos = match(muestrasFemeninas, fenotipos_filtrado$Mice)
dfFenotipos = fenotipos_filtrado[filasFenotipos, -1]
rownames(dfFenotipos) = fenotipos_filtrado[filasFenotipos, 1]

#Se vuelve a generar un dendrograma, esta vez sin los outliers de expresion 
#y mostrando los rasgos fenotipicos 
sampleTree2 = hclust(dist(df_expr_filtrado), method = "average")
#Codigos de color para fenotipos: valor bajo -> blanco,valor alto -> rojo, gris -> Nulo
traitColors = numbers2colors(dfFenotipos, signed = FALSE);
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(dfFenotipos), 
                    main = "Sample dendrogram and trait heatmap")

#save(datExpr, datTraits, file = "FemaleLiver-01-dataInput.RData")

#Codigo para escoger los potenciadores de coexpresion para calcular adyacencia, como soft-tresholds,
#basado en la topologia de red
powers = c(c(1:10), seq(from = 12, to=24, by=2))
sft = pickSoftThreshold(df_expr_filtrado, dataIsExpr = TRUE, powerVector = powers,  corFnc = cor,networkType = "signed")
sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)
ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.80, col = "red") +
  ylim(c(min(sft_df$model_fit), 1.05)) +
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  theme_classic()

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.80,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#PUNTO 2: Escoger un potenciador soft treshold para la produccion de la matriz de adyacencias y la red,
#basado en el umbral de R^2 establecido en la grafica
potenciador = 20 #cambiar de acuerdo a la gráfica

#Con este codigo se corre la red y la asignacion de modulos como clusters de genes 
#basados en su conectividad, todo se guarda como archivos de R aparte para manejar 
#la complejidad de espacio y tiempo de la computacion 
red = blockwiseModules(df_expr_filtrado, power = potenciador,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)
#podemos observar cuantos modulos y cuantos genes por modulo quedaron en el proceso de clustering
table(red$colors)

#grafica que genera un dendrograma de los modulos con colores especificos para cada uno
sizeGrWindow(12, 9)
mergedColors = labels2colors(red$colors)
plotDendroAndColors(red$dendrograms[[1]], mergedColors[red$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


etiquetasModulos = red$colors
coloresModulos = labels2colors(red$colors)
MEs = red$MEs
geneTree = red$dendrograms[[1]]

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(df_expr_filtrado, coloresModulos)$eigengenes
MEs = orderMEs(MEs0)
CorModuloFenotipo = cor(MEs, dfFenotipos, use = "p")
ValorPModuloFenotipo = corPvalueStudent(CorModuloFenotipo, nSamples)


#Esta grafica presenta las correlaciones y valores p entre fenotipos y clusters,
#lo cual permitira escoger los mas relevantes
textMatrix =  paste(signif(CorModuloFenotipo, 2), "\n(",
                    signif(ValorPModuloFenotipo, 1), ")", sep = "");
dim(textMatrix) = dim(CorModuloFenotipo)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = CorModuloFenotipo,
               xLabels = names(dfFenotipos),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Relaciones fenotipo-cluster"))

#PUNTO 3: Basado en la gráfica anterior escoger un modulo de interes diferente al default
#de acuerdo a su relacion con ciertos fenotipos y trabajar con este en adelante; ademas
#escoger un rasgo fenotipico con el cual se quiera trabajar, diferente a peso en este caso,
#cambiar el color para que quede ese modulo
modulo = "turquoise" #cambiar el color del modulo por el escogido
rasgo = as.data.frame(dfFenotipos$Glucose_Insulin) #en este caso se escoge el atributo weight_g, escoger otro
names(rasgo) = "rasgo" #no cambiar este String, el rasgo que escogio tendra este nombre en los atributos
nombresModulos = substring(names(MEs), 3)

pertenenciaGenModulo = as.data.frame(cor(df_expr_filtrado, MEs, use = "p"));
valorMMP = as.data.frame(corPvalueStudent(as.matrix(pertenenciaGenModulo), nSamples));

names(pertenenciaGenModulo) = paste("MM", nombresModulos, sep="");
names(valorMMP) = paste("p.MM", nombresModulos, sep="");

significanciaGenFenotipo = as.data.frame(cor(df_expr_filtrado, rasgo, use = "p"));
valorGSP = as.data.frame(corPvalueStudent(as.matrix(significanciaGenFenotipo), nSamples));

names(significanciaGenFenotipo) = paste("GS.", names(rasgo), sep="");
names(valorGSP) = paste("p.GS.", names(rasgo), sep="");


columna = match(modulo, nombresModulos);
genesModulo = coloresModulos==modulo;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(pertenenciaGenModulo[genesModulo, columna]),
                   abs(significanciaGenFenotipo[genesModulo, 1]),
                   xlab = paste("Module Membership in", modulo, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = modulo)

names(df_expr_filtrado)[coloresModulos=="brown"] #cambiar por el color de modulo escogido

#Se lee el archivo de anotaciones de genes
annot = read.csv(file = "GeneAnnotation.csv");
dim(annot)
names(annot)
probes = names(df_expr_filtrado)
probes2annot = match(probes, annot$substanceBXH)
sum(is.na(probes2annot))
#esta suma debe dar 0 ya que son los genes sin anotacion

#Crear dataframe con informacion de anotacion y valores calculados
infoGenes0 = data.frame(substanceBXH = probes,
                       geneSymbol = annot$gene_symbol[probes2annot],
                       LocusLinkID = annot$LocusLinkID[probes2annot],
                       moduleColor = coloresModulos,
                       significanciaGenFenotipo,
                       valorGSP)
ordenModulos = order(-abs(cor(MEs, rasgo, use = "p")))
for (mod in 1:ncol(pertenenciaGenModulo))
{
  nombresAnteriores = names(infoGenes0)
  infoGenes0 = data.frame(infoGenes0, pertenenciaGenModulo[, ordenModulos[mod]], 
                         valorMMP[, ordenModulos[mod]])
  names(infoGenes0) = c(nombresAnteriores, paste("MM.", nombresModulos[ordenModulos[mod]], sep=""),
                       paste("p.MM.", nombresModulos[ordenModulos[mod]], sep=""))
}
ordenGenes = order(infoGenes0$moduleColor, -abs(infoGenes0$GS.rasgo)) #
infoGenes = infoGenes0[ordenGenes, ]

#se imprime un dataframe con la informacion de cada gen anotada, su pertenencia a cierto
#modulo y la informacion de correlacion gen fenotipo
write.csv(infoGenes, file = "infoGenes.csv")

#Exportar a Cytoscape para visualizar el modulo escogido
TOM = TOMsimilarityFromExpr(df_expr_filtrado, power = potenciador)
modulos = c(modulo)
probes = names(df_expr_filtrado)
enModulo = is.finite(match(coloresModulos, modulos))
modProbes = probes[enModulo]
modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)]
modTOM = TOM[enModulo, enModulo]
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-",
                                                paste(modulos, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-",
                                                paste(modulos, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = coloresModulos[enModulo])



