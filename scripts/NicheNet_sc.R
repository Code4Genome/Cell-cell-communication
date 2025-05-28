#NicheNet network.

#read seurat object (NM0)
data2 = readRDS("NMO.rds")

data2@meta.data %>% head()
View(data2@meta.data)

data2@meta.data$singleR.labels %>% table() 

DimPlot(data2, reduction = "umap", group.by = "sample")

#Database (Human)

lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

lr_network = lr_network %>% distinct(from, to)
head(lr_network)

ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

ligand_target_matrix[1:5,1:5]

head(weighted_networks$lr_sig)
View(weighted_networks)

Idents(data2)

Idents(data2) = "singleR.labels"

#Perform NicheNet Analysis

DefaultAssay(data2) = "RNA"

DotPlot(data2, features = c("IL1B", "IL16", "CCL3"))
#Define a sender and a receiver cell population

senders = "B"
receiver = "C"

#get expressed genes from senders and receivers

list = senders %>% unique() %>%
  lapply(get_expressed_genes, data2, 0.10)

weighted_networks_lr = weighted_networks$lr_sig %>% 
  inner_join(lr_network, by = c('from', "to"))

expressed_genes = list %>% unlist() %>%  unique()

expressed_genes_receiver = get_expressed_genes(receiver, data2, pct = 0.10)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

#Define genes of interest
#these are the genes in the “receiver/target” cell population that are potentially affected by ligands expressed byinteracting cells (eg. DEGs)

geneset =  read.xlsx("upregulated_C.xlsx") %>% pull(gene)
geneset_up <- geneset %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.20) %>% pull(gene)

geneset <- geneset %>% .[. %in% rownames(ligand_target_matrix)]

#Define potential ligands

expressed_ligands = intersect(ligands, expressed_genes)
expressed_receptors = intersect(receptors, expressed_genes_receiver)

potential_ligands = lr_network %>%
  filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
  pull(from) %>% unique()

#perform NicheNet

options(timeout = 600)

ligand_activities = predict_ligand_activities(geneset = geneset, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix,
                                              potential_ligands = potential_ligands)

ligand_activities = arrange(ligand_activities, -aupr_corrected)

best_upstream_ligands = ligand_activities %>% top_n(30, aupr_corrected) %>%
  arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()

DotPlot(data2, features = best_upstream_ligands, cols = "RdYlBu") + RotatedAxis()

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() 

nrow(active_ligand_target_links_df)

head(active_ligand_target_links_df)

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

head(active_ligand_target_links)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])


make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

ordered_ligand_receptor_network <- vis_ligand_receptor_network[
  order(rownames(vis_ligand_receptor_network)),  # Sort ligands alphabetically
  order(colnames(vis_ligand_receptor_network))   # Sort receptors alphabetically
]

(make_heatmap_ggplot(t(ordered_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential_B cells"))


vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()))  
