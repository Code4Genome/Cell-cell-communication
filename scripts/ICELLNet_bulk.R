Icellent (BULK)

library(BiocGenerics)
library("org.Hs.eg.db")
library("hgu133plus2.db")
library(jetset)
library(ggplot2)
library(dplyr)
library(icellnet)
library(gridExtra)

# Load and select database
#db=as.data.frame(read.csv(curl::curl(url="https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/data/ICELLNETdb.tsv"), sep="\t",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""))

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
  ppi = F
)

db.name.couple=name.lr.couple(db)
head(db.name.couple)

# Load bulk RNAseq data
data = as.data.frame(read.table("final.txt", header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""))
rownames(data)
colnames(data)

PC.target <- data.frame(
  ID = c("B", "B", "B", "B","B", "B", "B", "B", "B",
         "C", "C", "C", "C", "C", "C", "C", "C", "C", "C",
        "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M"),
  Cell_type = c("B", "B", "B", "B","B", "B", "B", "B", "B",
                "C", "C", "C", "C", "C", "C", "C", "C", "C", "C",
                "M", "M", "M", "M", M", "M", "M", "M", M", "M", "M", "M")
  Class = "immune_cells"
)


# Data scaling
data.scaled=gene.scaling(data = data, db = db)

# data selection
CC.data.selection=PC.target$ID[which(PC.target$Cell_type=="C")] 
PC.data.selection=PC.target$ID[which(PC.target$Cell_type%in% c("B", "M"))]

my_Central_Cell_data=data.scaled[, CC.data.selection]
my_Partner_Cell_data = data.scaled[, PC.data.selection]

rm(target)
target = PC.target[-c(10:19),]
target

PC = unique(target$Cell_type)

PC = c("B", "M")

score.computation.1= icellnet.score(direction="in", PC.data=my_Partner_Cell_data, PC =PC, PC.target = target,
                                    CC.data= my_Central_Cell_data, CC.type = "RNAseq",  PC.type = "RNAseq",  db = db)
#undebug(icellnet.score)

score1=as.data.frame(score.computation.1[[1]])
lr1=score.computation.1[[2]] 

ymax=round(max(score1))+1

colnames(lr1)=c("B", "M")

LR.family.score(lr=lr1, db.couple=db, plot= "NULL")

LR.family.score(lr=lr1, db.couple=db.name.couple, plot="heatmap", title= "CD4 T cell")
  
LR.heatmap(lr = lr1, thresh = 0 , topn=60 , sort.by="sum",  title="Most contributing interactions_Bulk Data") +
  theme(text = element_text(family = "Verdana", size = 10, hjust = 0.8)) 


LR.heatmap(lr = lr1, thresh = 0 , topn=60 , sort.by="var",  title= "Most different interactions_Bulk Data")

pairs=LR.selection(lr = lr1, thresh = 0 , topn=10 , sort.by="var")
