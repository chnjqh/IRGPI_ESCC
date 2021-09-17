library(WGCNA)
##WGCNA##
powers <- 1:20
sft <- pickSoftThreshold(rnaExpr_ESCC_immune, powerVector = powers, verbose = 5)

par(mfrow = c(1, 2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], type = 'n',
     xlab = 'Soft Threshold (power)', ylab = 'Scale Free Topology Model Fit,signed R^2',
     main = paste('Scale independence'))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, col = 'red');
abline(h = 0.90, col = 'red')

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = 'Soft Threshold (power)', ylab = 'Mean Connectivity', type = 'n',
     main = paste('Mean connectivity'))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = 'red')
power <- sft$powerEstimate
power

sampleTree = hclust(dist(rnaExpr_ESCC_immune), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

adjacency <- adjacency(rnaExpr_ESCC_immune, power = power)
tom_sim <- TOMsimilarity(adjacency)
rownames(tom_sim) <- rownames(adjacency)
colnames(tom_sim) <- colnames(adjacency)
tom_sim[1:6,1:6]

tom_dis  <- 1 - tom_sim

geneTree <- hclust(as.dist(tom_dis), method = 'average')
plot(geneTree, xlab = '', sub = '', main = 'Gene clustering on TOM-based dissimilarity',
     labels = FALSE, hang = 0.04)

minModuleSize <- 30  
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = tom_dis,deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
table(dynamicMods)

dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, 'Dynamic Tree Cut',
                    dendroLabels = FALSE, addGuide = TRUE, hang = 0.03, guideHang = 0.05,
                    main = 'Gene dendrogram and module colors')

plot_sim <- -(1-tom_sim)
#plot_sim <- log(tom_sim)
diag(plot_sim) <- NA
TOMplot(plot_sim, geneTree, dynamicColors,
        main = 'Network heatmap plot, selected genes')

MEList <- moduleEigengenes(rnaExpr_ESCC_immune, colors = dynamicColors)
MEs <- MEList$eigengenes
head(MEs)[1:6]

write.table(MEs, 'moduleEigengenes.txt', sep = '\t', col.names = NA, quote = FALSE)
trait <- read.delim('trait.txt', row.names = 1, check.names = FALSE)
module <- MEs

moduleTraitCor <- cor(module, trait, use = 'p')
moduleTraitCor[1:6,1:6]  
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(module))

textMatrix <- paste(signif(moduleTraitCor, 2), '\n(', signif(moduleTraitPvalue, 1), ')', sep = '')
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, main = paste('Module-trait relationships'),
               xLabels = names(trait), yLabels = names(module), ySymbols = names(module),
               colorLabels = FALSE, colors = blueWhiteRed(50), cex.text = 0.7, zlim = c(-1,1),
               textMatrix = textMatrix, setStdMargins = FALSE)

gene_module <- data.frame(gene_name = colnames(rnaExpr_ESCC_immune), module = dynamicColors, stringsAsFactors = FALSE)
head(gene_module)

gene_module_select_brown <- subset(gene_module, module == 'brown')$gene_name
gene_module_select_turquoise <- subset(gene_module, module == 'turquoise')$gene_name
gene_module_select<-c(gene_module_select_brown,gene_module_select_turquoise)
gene_select <- t(rnaExpr_ESCC_immune[,gene_module_select])


dir.create('cytoscape', recursive = TRUE)
for (mod in 1:nrow(table(dynamicColors))) {
  modules <- names(table(dynamicColors))[mod]
  probes <- colnames(rnaExpr_ESCC_immune)
  inModule <- (dynamicColors == modules)
  modProbes <- probes[inModule]
  modGenes <- modProbes
  modtom_sim <- tom_sim[inModule, inModule]
  dimnames(modtom_sim) <- list(modProbes, modProbes)
  outEdge <- paste('cytoscape/', modules, '.edge_list.txt',sep = '')
  outNode <- paste('cytoscape/', modules, '.node_list.txt', sep = '')
  exportNetworkToCytoscape(modtom_sim,
                           edgeFile = outEdge,
                           nodeFile = outNode,
                           weighted = TRUE,
                           threshold = 0.2,  
                           nodeNames = modProbes,
                           altNodeNames = modGenes,
                           nodeAttr = dynamicColors[inModule])
}