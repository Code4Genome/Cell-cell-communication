remotes::install_github('saezlab/liana')

path = system.file(package = "liana")

show_resources()
show_methods()

#check database------------------------------

#lr_db <- liana::select_resource("Consensus")
#data2 = lr_db$Consensus

#write.xlsx(data2, "database.xlsx")

#---------------------------------------------
data = readRDS("NMO_annotated.rds")

View(data@meta.data)
data@assays$RNA@data

#Run liana

Idents(data) = "singleR.labels"

colnames(data@meta.data)[colnames(data@meta.data) == "singleR.labels"] <- "colLabels"

liana = liana_wrap(data, #method = c("cellphonedb", "sca"),
                   resource = c("Consensus"),
                   idents_col = "singleR.labels",  min_cells = 5)

liana %>% dplyr::glimpse()

liana <- liana %>%
  liana_aggregate()

dplyr::glimpse(liana)

cpdb_int <- liana %>%
  # only keep interactions with p-val <= 0.05
  filter(cellphonedb.pvalue <= 0.05) # this reflects interactions `specificity`
  # then rank according to `magnitude` (lr_mean in this case)
 #rank_method(method_name = "cellphonedb",
              #mode = "magnitude") %>%
  # keep top interactions (regardless of cell type)
  #distinct_at(c("ligand.complex", "receptor.complex")) %>%
  #head(50)

cpdb_int %>%
  liana_dotplot(source_groups = c("Monocytes"),
                target_groups = c("CD8+ T cells"),
                ntop = 70, 
                specificity = "cellphonedb.pvalue",
                size.label = "pvalue")

int2 = data.table(cpdb_int)

int2 = int2[source == "M cells" & target == "C cells"]

int3 = data.table(cpdb_int)

int3 = int3[source == "A cells" & target == "C cells"]

df_unique = int2[!duplicated(int2$aggregate_rank)]

df_unique_2 = int3[!duplicated(int3$aggregate_rank)]
