install.packages("devtools")
devtools::install_github("sqjin/Cellchat")
install.packages("NMF")
devtools::install_github("jokergoo/circlize")
devtools::install_github("jokergoo/ComplexHeatmap")
#python-pacakge: pip install umap-learn (installed correctly?)

remotes::install_version("igraph", version = "1.3.5") #1.4.1 is not compatible with cellchat

options(stringsAsFactors = FALSE)

#data

data = read.table("final.txt")
data


data1 = as.matrix(data)
data1

ncol(data1)
nrow(data1)
#rownames(data1)
colnames(data1)

#meta

meta = read.csv("mydata.csv")
meta

meta1 = as.data.frame(meta)
meta1

row.names(meta1)
colnames(meta1)

rownames(meta1) = c("B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "B_11", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", 
                    "C_8", "C_9", "M_1","M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "M_11")
row.names(meta1)

all(colnames(data1) == row.names(meta1))

setdiff(colnames(data1), rownames(meta1))

colnames(data1) = rownames(meta1)

#cellchat object 

cellchat = createCellChat(object = data1, meta = meta1, group.by = "condition")
cellchat

#set LR interaction database

CellChatDB = CellChatDB.human
showDatabaseCategory(CellChatDB) #show multpile categories available

#use a subset of Cellchat for cell-cell communication analysis
CellChatDB.use = subsetDB(CellChatDB, search = "Secreted Signaling") #use Secreted signaling

CellChatDB.use

cellchat@DB = CellChatDB.use

#Additional QC
#subset the expression data of signaling genes (Only those genes related to signaling, ignoring the rest)

cellchat = subsetData(cellchat)
future::plan("multisession", workers = 4) # of threads

#overexpressed genes are identified
cellchat = identifyOverExpressedGenes(cellchat)
cellchat

#overexpressed interactions have to be identified
cellchat = identifyOverExpressedInteractions(cellchat)
cellchat

#project gene expression data onto PPI (optional)
#cellchat = projectData(cellchat, PPI.human)

#calculate probabilities of all PPI 

cellchat = computeCommunProb(cellchat, raw.use = TRUE)

#let's filter bad-quality communications/less number of cells per group

#cellchat = filterCommunication(cellchat, min.cells = 10)

#returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors
df.net = subsetCommunication(cellchat, sources.use = c("B", 'M'), targets.use = c("C"))
df.net

df_final1 = write.csv(df.net, "df.Signaling.csv")

#calculate the pathways probability for CCC found above

cellchat = computeCommunProbPathway(cellchat)

#calculate the aggregated CCC network

cellchat = aggregateNet(cellchat, sources.use = c("Bcell", "Mono"), targets.use = c("CD4"), remove.isolate = FALSE)

#netwotk visualization

groupSize = as.numeric(table(cellchat@idents))
groupSize
par(mfrow = c(1,2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weight/strength")

#Visualization of each cell type againts other cell types 1 by 1

mat = cellchat@net$weight
mat
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#let's visualize all the pathways

cellchat@netP$pathways


p = saveRDS(cellchat@netP$pathways, "pathways_signaling.rds")
p1 = readRDS("pathways_signaling.rds")
p1

pathways.show = p1

#Hierarchy plot

#Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling 
(levels(cellchat@idents))

vertex.receiver = seq(1,2) #numeric vector
vertex.receiver
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()

# Chord diagram (total)
par(mfrow=c(1,1))
pdf(file = "signaling pathway_cord.pdf", width = 20, height = 16) #it takes time
netVisual_aggregate(cellchat, signaling = "IL6", layout = "chord", color.use = c("darkgreen", "darkorange"))
dev.off()

# Heatmap
par(mfrow=c(1,2))
netVisual_heatmap(cellchat, signaling = c("IL6"), color.heatmap = "Blues",color.use = c("darkorange", "darkgreen", "darkred"))
dev.off()

#Let's visualize highest number of interactions of FN1 pathway
#netAnalysis_contribution(cellchat, signaling = pathways.show = c("FN1"), title = "Contribution of each LR pair/all"
netAnalysis_contribution(cellchat, signaling = "IL6", title = "Contribution of each LR pair/IL6")

netAnalysis_contribution(cellchat, signaling = "IL6", title = "Contribution of each LR pair/IL6", return.data = T)

pairLR.THBS = extractEnrichedLR(cellchat_aqp4, signaling = "IL6", geneLR.return = FALSE)
pairLR.THBS
LR.show = pairLR.THBS[1:1, ]
LR.show

#Hierarchy plot
vertex.receiver = seq(1,2) # a numeric vector
netVisual_individual(cellchat, signaling = "IL6", pairLR.use = LR.show, vertex.receiver = vertex.receiver)
netVisual_individual(cellchat, signaling = "IL6", pairLR.use = LR.show, layout = "circle")

#show all significant interactions (LR pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = c("Bcell","Mono") , targets.use = c("CD4"), signaling = c("ANGPT", "CCL", "CD70", "COMPLEMENT", "CXCL", 
                                                                                                   "GALECTIN","IL1", "IL10", "IL12","IL2", "NPR2", "IL6", "TGFb", "SLURP", "FGF"),remove.isolate = FALSE, angle.x = 45)
dev.off()

#Chord figure (all interactions)
par(mfrow=c(1,2))
pdf(file ="cord_from bubble_new.pdf", width = 20, height =16)
netVisual_chord_gene(cellchat, sources.use = c("Bcell", "Mono"), targets.use = c("CD4"), 
                     signaling = c("IL2", "TGFb", "IL6", 
                                   "VEGI", "CXCL", "IL1","CCL", "IL10", "NRG","ACTIVIN", "FGF", "SEMA3", "IL12"), lab.cex = 1, legend.pos.y = 100, color.use = c("darkblue", "darkred", "darkgreen")) 
dev.off()

plotGeneExpression(cellchat, signaling = "IL6")

info = saveRDS(cellchat,"cellchat_signaling.rds")

#compute network centrality scores 
cellchat = netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = "TGFb", width = 8, height = 2.5, font.size = 10, color.use = c("darkgreen", "darkorange", "darkred"), color.heatmap = "Blues")

ht1 = netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", signaling =c("ANGPT", "CCL", "CD70", "COMPLEMENT", "CXCL", 
                                                                                     "GALECTIN","IL1", "IL10", "IL12","IL2", "NPR2", "IL6", "TGFb", "SLURP", "FGF"), color.use = c("darkgreen", "darkorange", "darkred"), color.heatmap = "Blues")
ht1

ht2 = netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", signaling =c("ANGPT", "CCL", "CD70", "COMPLEMENT", "CXCL", 
                                                                                     "GALECTIN","IL1", "IL10", "IL12","IL2", "NPR2", "IL6", "TGFb", "SLURP", "FGF"), color.use = c("darkgreen", "darkorange", "darkred"), color.heatmap = "Blues")
ht2

ht3 = computeAveExpr(cellchat, features = c("IL6"), type = "truncatedMean", trim = 0.1)
ht3
