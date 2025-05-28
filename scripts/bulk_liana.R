library(SingleCellExperiment)
library(liana)
library(dplyr)

#singlecellexperiment
#---------------------------------------------------------------------------------------------

ncells = 31

metadata = as.data.frame(SingleCellExperiment::colData(sce))

u = read.table("total_c.txt", fill = TRUE, header= T ,sep="\t",check.names= F, row.names = 1)
v = read.table("final.txt") 

pca <- matrix(runif(ncells*31), ncells)
tsne <- matrix(rnorm(ncells*31), ncells)

Cell_type = c("B", "B", "B", "B","B", "B", "B", "B", "B", "B",
              "C", "C", "C", "C", "C", "C", "C", "C", "C", "C",
              "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A")

celltype = as.data.frame(Cell_type)

colData(sce)$celltype = Cell_type

head(colData(sce))

colData(sce)

colData(sce)$celltype <- factor(colData(sce)$celltype)

liana = liana_wrap(sce, #method = c("cellphonedb", "sca", "natmi"),
                   resource = c("Consensus"), 
                   idents_col = "celltype")

liana %>% dplyr::glimpse()

liana <- liana %>%
  liana_aggregate()

dplyr::glimpse(liana)


cpdb_int <- liana %>%
  # only keep interactions with p-val <= 0.1
  filter(cellphonedb.pvalue <= 0.1) 
  # keep top 20 interactions (regardless of cell type)

cpdb_int %>%
  liana_dotplot(source_groups = c("Bcell"),
                target_groups = c("CD4"),
                #ntop = 70, 
                specificity = "aggregate_rank",
                size.label = "pvalue")

int = data.table(cpdb_int)

int2 = int[source == "B" & target == "C"]
