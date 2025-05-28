library(BiocGenerics)
library("org.Hs.eg.db")
library("hgu133plus2.db")
library(jetset)
library(ggplot2)
library(dplyr)
library(icellnet)
library(gridExtra)

db=as.data.frame(read.csv(curl::curl(url="https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/data/ICELLNETdb.tsv"), sep="\t",header = T, check.names=FALSE, 
                          stringsAsFactors = FALSE, na.strings = ""))
db.name.couple=name.lr.couple(db, type="Family")
head(db.name.couple)

db = CellPhoneDB_convert(
  complex = utils::read.csv(curl::curl(url =
                                         "https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/complex_input.csv"),
                            sep = ",", header = T, check.names = FALSE, stringsAsFactors = FALSE, na.strings =
                              ""),
  interaction = utils::read.csv(curl::curl(url =
                                             "https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/interaction_input.csv"),
                                sep = ",", header = T, check.names = FALSE, stringsAsFactors = FALSE, na.strings =
                                  ""),
  gene_info = utils::read.csv(curl::curl(url =
                                           "https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/gene_input.csv"),
                              sep = ",", header = T, check.names = FALSE, stringsAsFactors = FALSE, na.strings =
                                ""),
  ppi = T
)


seurat <- readRDS(file = "NMO.rds")
# keep genes at least expressed by defined percentage of cell in their respective cluster (below 10%)
filter.perc=10
Idents(seurat)=seurat$singleR.labels

average.clean= sc.data.cleaning(object = seurat, db = db, filter.perc = filter.perc, save_file = F, force.file = F)
dim(average.clean)

data.icell=as.data.frame(gene.scaling(as.data.frame(average.clean), n=1, db= db))

senders = c("B", "C")
receivers = "A"

PC.data=as.data.frame(data.icell[,c(senders,"Symbol")], row.names = rownames(data.icell))
my.selection= senders

rownames(PC.data)
rownames(CC.data)

# Extract receiver data for CC.data
CC.data = as.data.frame(data.icell[, receivers, drop = FALSE], row.names = rownames(data.icell))

score.computation.1= icellnet.score(direction="in", PC.data=PC.data, 
                                    CC.data= CC.data,  
                                    PC=my.selection, CC.type = "RNAseq", 
                                    PC.type = "RNAseq",  db = db)
score.computation.1
str(score.computation.1)

score1=as.data.frame(score.computation.1[[1]])

lr1=(score.computation.1[[2]])

lr1

ymax=round(max(score1))+ 1 #to define the y axis range of the barplot

contribution = LR.family.score(lr=lr1, db.couple= db.name.couple, plot= NULL)
contribution

LR.family.score(lr=lr1, db.couple= db.name.couple, plot="heatmap", title = "CD4 T cells")

colnames(lr1)=c("B_to_C", "A_to_C")
var = LR.heatmap(lr = lr1, thresh = 0 , topn=50 , sort.by="var",  title="Most different interactions")  +
  theme(text = element_text(family = "Verdana", size = 10, hjust = 0.8)) 

colnames(lr1)=c("B_to_T", "A_to_C")
sum = LR.heatmap(lr = lr1, thresh = 0 , topn=50 , sort.by="sum",  title="Most contributing interactions") +
  theme(text = element_text(family = "Verdana", size = 10, hjust = 0.8)) 


lr1 = as.matrix(lr1)

colSums(lr1)

data = colSums(lr1[,1,drop = FALSE])

pheatmap(lr1)

contribution = matrix(nrow = (length(family) + 1), ncol = length(lr1[1,]), 0)

lr1 = matrix(lr1, byrow = T)

lr_contributions = LR.family.score(lr=lr1, db.couple=db.name.couple, plot=NULL)

data2 = lr_contributions[ ,1, drop = FALSE]

str(lr_contributions)
head(lr_contributions)
dim(data2)
LR.family.score(lr=data, db.couple=db.name.couple, plot="heatmap", title="CCC")
