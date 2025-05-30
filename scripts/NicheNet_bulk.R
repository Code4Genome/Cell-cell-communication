#Nichenet using bulk data

#Database (Human)

lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

lr_network <- lr_network %>% distinct(from, to)
head(lr_network)

ligand_target_matrix[1:5,1:5]

head(weighted_networks$lr_sig)

head(weighted_networks$gr)

data = read.table("final.txt", header = TRUE, row.names = 1)
meta = read.csv("mydata.csv")

converted <- convert_alias_to_symbols(rownames(data), "human", verbose = FALSE)

duplicated_genes <- converted[duplicated(converted)]
print(duplicated_genes)

unique_genes <- !duplicated(converted)
data <- data[unique_genes, ]
rownames(data) <- converted[unique_genes]

senders = c("B")
receiver = "C"

# Define sample groups
sender_cells <- c("B1", "B2","B3", "B4","B5"","B6","B7","B8")  # Modify based on column names
receiver_cells <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12")  # Modify based on column names

sender_cells_mono <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11")

sender_bulk <- data[, sender_cells]
receiver_bulk <- data[, receiver_cells]
mono_bulk <- data[, sender_cells_mono]

expressed_genes_sender = rownames(sender_bulk)[rowMeans(sender_bulk) > 1]

expressed_genes_receiver = rownames(receiver_bulk)[rowMeans(receiver_bulk) > 1]

expressed_genes_sender_mono = rownames(mono_bulk)[rowMeans(mono_bulk) > 1]

expressed_ligands <- expressed_genes_sender_mono [expressed_genes_sender_mono%in% lr_network$from]
expressed_receptors <- expressed_genes_receiver [expressed_genes_receiver%in% lr_network$to]

length(expressed_ligands)
length(expressed_receptors)

#senders = c("M")
#receiver = "C"

ligands <- lr_network %>% pull(from) %>% unique()
expressed_ligands <- intersect(ligands,expressed_ligands)

receptors <- lr_network %>% pull(to) %>% unique()
expressed_receptors <- intersect(receptors,expressed_receptors)

potential_ligands <-  lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% 
  pull(from) %>% unique()

head(potential_ligands)

geneset =  read.xlsx("up_CD4_AQP4.xlsx") %>% pull(gene_name)

geneset <- geneset %>% .[. %in% rownames(ligand_target_matrix)]

length(geneset)

background_expressed_genes <- expressed_receptors %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)

ligand_activities <- predict_ligand_activities(geneset = geneset,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

(ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>%
    mutate(rank = rank(desc(aupr_corrected))))

best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>%
  arrange(-aupr_corrected) %>% pull(test_ligand)

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset,
         ligand_target_matrix = ligand_target_matrix,
         n = 200) %>% bind_rows()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.25)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target_network <- make_heatmap_ggplot(vis_ligand_target, "Prioritized CAF-ligands", "p-EMT genes in malignant cells",
                                               color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

p_ligand_target_network

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
                     y_name = "Prioritized ligands", x_name = "Receptors expressed by C cells",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))
