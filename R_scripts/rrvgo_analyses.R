
# Load libraries
library(rrvgo)
library(dbplyr)
library(dplyr)
library(ggplot2)
library(piano)
library(snowfall)
library(clusterProfiler)
library(ggdendro)
library(org.Hs.eg.db)
library(tidyr)
library(pheatmap)

################################################################################################
# Set inputs
################################################################################################

###################### For TM protein data

# CC output - uncomment following 8 lines:
#set_ontology = "CC"
#groups <-  list("CC_celltype" = read.csv("../results/piano/piano_CC/piano_CC_celltype.csv", sep="\t")) # CC piano result
#DEGs_df <- read.csv("../results/DE/contrasts/DEGs_Vulko_vs_Expi.csv")
#group_of_interest = "CC_celltype" # CC rrvgo output
#path_name = "../results/rrvgo/CC_celltype_top50/" # CC rrvgo output
#save_file <- "TM_CC_heatmap.pdf"
#json_file <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c5.go.cc.v2023.2.Hs.json" # Corresponding json file - CC
#GO_XX <- "GOCC_"

# BP output - uncomment following 8 lines:
#set_ontology = "BP"
#groups <-  list("BP_celltype" = read.csv("../results/piano/piano_BP/piano_BP_celltype.csv", sep="\t")) # BP piano results
#DEGs_df <- read.csv("../results/DE/contrasts/DEGs_Vulko_vs_Expi.csv")
#group_of_interest = "BP_celltype" # BP rrvgo output
#path_name = "../results/rrvgo/BP_celltype_top50/" # BP rrvgo output
#save_file <- "TM_BP_heatmap.pdf"
#json_file <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c5.go.bp.v2023.2.Hs.json" # Corresponding json file - BP
#GO_XX <- "GOBP_"


###################### For CriptoFC data 

# Negavtive vs Negative (vulko vs expi) CC
#set_ontology = "CC"
#groups <-  list("CC_celltype" = read.csv("../CriptoFc/results/piano/piano_CC/piano_CC_Negative.csv", sep="\t")) # CC piano results
#DEGs_df <- read.csv("../CriptoFc/results/DE/contrasts/DEGs_Vulko293_neg_vs_Expi293_neg.csv")
#group_of_interest = "CC_celltype" # CC rrvgo output
#path_name = "../CriptoFc/results/rrvgo/CC_Negative_top50/"
#save_file <- "rrvgo_CC_Negative.pdf"
#json_file <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c5.go.cc.v2023.2.Hs.json" # Corresponding json file - CC
#GO_XX <- "GOCC_"

# Negavtive vs Negative (vulko vs expi) BP
#set_ontology = "BP"
#groups <-  list("BP_celltype" = read.csv("../CriptoFc/results/piano/piano_BP/piano_BP_Negative.csv", sep="\t")) # BP piano results
#DEGs_df <- read.csv("../CriptoFc/results/DE/contrasts/DEGs_Vulko293_neg_vs_Expi293_neg.csv")
#group_of_interest = "BP_celltype" # BP rrvgo output
#path_name = "../CriptoFc/results/rrvgo/BP_Negative_top50/"
#save_file <- "rrvgo_BP_Negative.pdf"
#json_file <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c5.go.bp.v2023.2.Hs.json" # Corresponding json file - BP
#GO_XX <- "GOBP_"


# CriptoFC vs CriptoFC (vulko vs expi) BP
set_ontology = "BP"
groups <-  list("BP_celltype" = read.csv("../CriptoFc/results/piano/piano_BP/piano_BP_CriptoFC.csv", sep="\t")) # BP piano results
DEGs_df <- read.csv("../CriptoFc/results/DE/contrasts/DEGs_Vulko293_CriptoFC_vs_Expi293_CriptoFC.csv")
group_of_interest = "BP_celltype" # BP rrvgo output
path_name = "../CriptoFc/results/rrvgo/BP_CriptoFC_top50/"
save_file <- "rrvgo_BP_CriptoFC.pdf"
json_file <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c5.go.bp.v2023.2.Hs.json" # Corresponding json file - BP
GO_XX <- "GOBP_"


################################################################################################
# A. Load PIANO results
################################################################################################


# A.1 Get all Piano results of interest
################################################################################################

datalist = list()

x <- 1
i <- 1
for (df in groups) {
  
  # Create one df for each direction for one comparison
  dn <- df %>% dplyr::select(Name, p.adj..dist.dir.dn., Genes..tot., Genes..up., Genes..down.) %>% dplyr::rename("padj" = "p.adj..dist.dir.dn.")
  up <- df %>% dplyr::select(Name, p.adj..dist.dir.up., Genes..tot., Genes..up., Genes..down.) %>% dplyr::rename("padj" = "p.adj..dist.dir.up.")
  mix_up <- df %>% dplyr::select(Name, p.adj..mix.dir.up., Genes..tot., Genes..up., Genes..down.) %>% dplyr::rename("padj" = "p.adj..mix.dir.up.")
  mix_dn <- df %>% dplyr::select(Name, p.adj..mix.dir.dn., Genes..tot., Genes..up., Genes..down.) %>% dplyr::rename("padj" = "p.adj..mix.dir.dn.")
  mix_nondir <- df %>% dplyr::select(Name, p.adj..non.dir.., Genes..tot., Genes..up., Genes..down.) %>% dplyr::rename("padj" = "p.adj..non.dir..")
  
  # Create named list of dataframes
  df_name <- names(groups[i])
  sets <- list()
  sets[[paste(df_name,"_up",sep="")]] <- up
  sets[[paste(df_name,"_dn",sep="")]] <- dn
  sets[[paste(df_name,"_mix_up",sep="")]] <- mix_up
  sets[[paste(df_name,"_mix_dn",sep="")]] <- mix_dn
  sets[[paste(df_name,"_mix_nondir",sep="")]] <- mix_nondir
  
  direction <- c("up", "dn", "mix_up", "mix_dn", "mix_nondir")
  
  k <- 1
  
  # Loop through set and add column with: 1. Group name (df_name_DIR),  2. direction, 3. Comparison (df_name)
  for (set in sets){
    set_name <- names(sets[k])
    df <- set
    df <- df %>%
      dplyr::select(Name, padj, Genes..tot., Genes..up., Genes..down.) %>%    
      arrange(padj) %>%
      mutate(Group = set_name)  %>%
      mutate(Direction = direction[k]) %>%
      mutate(Comparison = df_name)
    
    datalist[[x]] <- df
    
    k <- k + 1
    x <- x + 1
  }
  i <- i + 1
}

all_dfs = do.call(rbind, datalist)


# A.2 Use json file to map ontology "Standard name" to the "Exact source" -> meaning, get the GO:XXXX
################################################################################################

# Load corresponding json
json_file <- jsonlite::fromJSON(json_file, flatten = TRUE)

# Extract exactSource (GO:XXXXX) from json file 
data_list <- list()
for (group_name in names(json_file)) {
  exactSources <- json_file[[group_name]]$exactSource
  groups <- rep(group_name, length(exactSources))
  data_list[[group_name]] <- data.frame(group = groups, exactSource = exactSources)
}

df_json <- do.call(rbind, data_list)

# Add GO:XXXX to df
dfs_GO <- merge(all_dfs, df_json, by.x = "Name", by.y = "group", all.x = TRUE) # Keep all rows present in og df
dfs_GO <- dfs_GO[!is.na(dfs_GO$exactSource), ]



################################################################################################################################################################
# B. Get DEGs (%&average(abs(log2FC)) & sum (abs(log2FC))) for each ontology
################################################################################################################################################################

# B.1 Function: get DEGs for each ontology and corresponding: avg(abs(log2FC)), sum(abs(log2FC)) 
################################################################################################

get_go_genes <- function(go_term, DEGs_df) {
  
  ############# ENTREZ
  entrez_ids <- mget(go_term, org.Hs.egGO2ALLEGS, ifnotfound = NA)[[go_term]]
  if (is.null(entrez_ids)) {
    return(list("go_term" = go_term, "gene_list" = NA, "avg_log2_FC" = NA, "sum_log2FC" = NA))
  }
  
  entrez_ids <- as.character(entrez_ids) # Character vector
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  
  ############# GENE SYMBOLS
  gene_symbols <- tryCatch({
    mget(entrez_ids, org.Hs.egSYMBOL, ifnotfound = NA)
  }, error = function(e) {
    return(NA)
  })
  
  # If gene symbols only contain NA, return NA
  if (length(gene_symbols) == 1 && is.na(gene_symbols)) {
    return(list("go_term" = go_term, "gene_list" = NA, "avg_log2_FC" = NA, "sum_log2FC" = NA))
  }
  
  gene_symbols <- unlist(gene_symbols, use.names = FALSE)
  gene_symbols <- gene_symbols[!is.na(gene_symbols)]
  
  ############# Join gene_symbols and DEGS
  df <- DEGs_df %>% filter(X %in% gene_symbols)
  
  # If no overlap between gene_symbols and DEGS, return NA
  if (nrow(df) == 0) {
    return(list("go_term" = go_term, "gene_list" = NA, "avg_log2_FC" = NA, "sum_log2FC" = NA))
  }
  
  df$log2FoldChange <- abs(df$log2FoldChange) # Absolute log2 FC
  gene_list <- paste0('c("', paste(df$X, collapse = '", "'), '")')
  
  avg_log2_FC <- mean(df$log2FoldChange, na.rm = TRUE)
  sum_log2FC <- sum(df$log2FoldChange, na.rm = TRUE)
  
  return(list("go_term" = go_term, "gene_list" = gene_list, "avg_log2_FC" = avg_log2_FC, "sum_log2FC" = sum_log2FC))
}


# B.2 Get DEGs for specific comparison (ex for DEGs_Vulko_vs_Expi.csv)
################################################################################################

# Select comparison
GO_df <- subset(dfs_GO, Comparison == group_of_interest)

# Get DEGs for each ontology and corresponding info
unique_GOs <- unique(GO_df$exactSource)

# Get DEGs from GO
GO_DEGs_result <- lapply(unique_GOs, function(unique_GOs) get_go_genes(unique_GOs, DEGs_df))
GO_DEGs_result <- purrr::keep(GO_DEGs_result, ~ !is.na(.x$gene_list))

GO_DEGs_df <- do.call(rbind, lapply(GO_DEGs_result, function(x) {
  data.frame(
    go_term = x$go_term,
    avg_log2_FC = x$avg_log2_FC,
    sum_log2FC = x$sum_log2FC,
    gene_list = x$gene_list,
    stringsAsFactors = FALSE
  )
}))

# Join
GO_DEGs_df <- left_join(GO_df, GO_DEGs_df, by = c("exactSource" = "go_term"))




################################################################################################################################################################
# C. REVIGO & DENDROGRAMS FUNCTIONS
################################################################################################################################################################


# 1. Select top 50 GOs
# Select GOs > 5 genes in direction of interest
# Arrange after lowest padj and largest average log2FC of genes in each set

# 2. Run ReviGO
# Generate similarity matrix -> Resnik
# Scores -> here padj
# reduceSimMatrix -> threshold 0.9

# 3. Select representative of cluster
# Sort after cluster, pDadj, termUniqueness (lower, less unique)
# Get count of GOs in each cluster


# C.1 Run REVIGO & select representative cluster:
################################################################################################

revigo_cluster_terms <- function(df_in, group_of_interest, GO_DEGs_df, dir, path_name){
  
  # 1. Prep df for revigo
  revigo_df <- df_in %>%
    filter(Comparison == group_of_interest & Genes..tot. > 5 & Direction == dir) %>%
    arrange(padj, desc(avg_log2_FC)) %>%
    dplyr::select(exactSource, padj, Genes..tot., avg_log2_FC, sum_log2FC, gene_list, Direction) %>%
    head(50)
  
  # 2. REVIGO
  simMatrix <- calculateSimMatrix(revigo_df$exactSource,orgdb="org.Hs.eg.db",ont=set_ontology,method="Resnik")
  scores <- setNames(-log10(revigo_df$padj), revigo_df$exactSource)
  reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold=0.9, orgdb="org.Hs.eg.db")
  
  # Get info from selected reducedTerms
  df_reduced <- merge(reducedTerms, revigo_df, by.x = "go", by.y = "exactSource", all.x = TRUE) # Keep all rows present in og df
  
  # SAVE ALL ONTOLOGIES
  df_reduced_save <- df_reduced %>% arrange(cluster, padj, desc(size)) # SORT LFC
  write.csv(df_reduced_save, file = paste0(path_name, "reduced", group_of_interest, "_", dir, ".csv"), row.names = FALSE)
  
  # 3. Select representative GOs of cluster
  df_rep_GOs <- df_reduced %>%
    group_by(cluster, Direction) %>%
    mutate(gos_per_cluster = n()) %>%
    ungroup() %>%
    arrange(cluster, padj, desc(size)) %>% # SELECT GO FROM LOWEST PADJ & HIGHEST NR GENES
    distinct(cluster, .keep_all = TRUE) %>% # Keep first occurence of groups  (representative)
    dplyr::select(go, Direction, avg_log2_FC, gos_per_cluster, term, parentTerm)
  
  # CREATE DENDOGRAMS
  make_dendogram(df_in, df_rep_GOs, GO_DEGs_df, simMatrix, group_of_interest, dir, path_name)
  
  # Renmae
  names(df_rep_GOs)[names(df_rep_GOs) == 'Direction'] <- 'revigo_source'
  return(df_rep_GOs)
}


# C.2 - Dendogram plot function
################################################################################################

make_dendogram <- function(df_in, df_rep_GOs, GO_DEGs_df, simMatrix, group_of_interest, dir, path_name){
  
  label_df <- GO_DEGs_df %>% dplyr::select(Name, exactSource) %>% distinct(Name, exactSource) # Get GOs and terms
  label_df$Name <- label_df$Name %>% gsub(GO_XX, "", .) %>% gsub("_", " ", .) %>%  tolower() # Convert terms to lower key and remove GOBP_ or GOCC_
  
  # Get distance matrix and cluster complete linkage (as rrvgo does)
  distMatrix <- 1 - simMatrix
  hc <- hclust(as.dist(distMatrix), method = "complete")
  dendr <- ggdendro::dendro_data(hc)
  threshold_value = 0.9
  
  # Get clusters with threshold 0.9, used in rrvgo
  cluster <- cutree(hc, h = threshold_value)
  clust.df <- data.frame(label=names(cluster), cluster=factor(cluster))
  
  # Add terms instead of ontology ID as leaf labels
  dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")
  dendr[["labels"]] <- merge(dendr[["labels"]],label_df, by.x="label", by.y="exactSource", all.x = TRUE)
  
  # Add terms as cluster labels
  df_rep_GOs <- df_rep_GOs %>% left_join(dendr[["labels"]], by = c("go" = "label")) %>% dplyr::select(cluster, term) # create df linking term to cluster
  dendr[["labels"]] <- merge(dendr[["labels"]], df_rep_GOs, by = "cluster", all.x = TRUE) # add terms to df
  
  dendogram_plot <- ggplot() + 
    geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
    geom_text(data=label(dendr), aes(x, y, label=Name, hjust=0, color=term), size=3, key_glyph="point") +
    labs(color = "  Representative GO") +
    guides(shape = guide_legend(override.aes = list(size = 1))) +
    guides(color = guide_legend(override.aes = list(size = 1))) +
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.spacing.x = unit(1, 'mm')) +
    coord_flip() +
    scale_y_reverse(expand=c(0.8, 0)) + 
    geom_hline(yintercept = threshold_value, color = "red", linetype = "dashed") +
    theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(),axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank(),
          aspect.ratio=1,
          plot.margin = margin(5, 0, 5, 0))
    
  #plot.margin = margin(5, 0, 5, 0)) +
  #scale_x_continuous(expand = expansion(mult = c(0, 0)))
  
  pdf(file=paste0(path_name, "dendo_", group_of_interest, "_", dir, ".pdf"), width = 15, height = 10)
  print(dendogram_plot)
  dev.off()
  
}


##########################################################################################
#  D. RUN REVIGO AND DENDROGRAMS FOR COMPARISON OF INTEREST
##########################################################################################

# D.1 Run REVIGO for each direction for specified data set
################################################################################################

df = GO_DEGs_df
df_up <- revigo_cluster_terms(df, group_of_interest, GO_DEGs_df, "up", path_name)
df_mix_up <- revigo_cluster_terms(df, group_of_interest, GO_DEGs_df, "mix_up", path_name)
df_mix <- revigo_cluster_terms(df, group_of_interest, GO_DEGs_df, "mix_nondir", path_name)
df_mix_dn <- revigo_cluster_terms(df, group_of_interest, GO_DEGs_df, "mix_dn", path_name)
df_dn <- revigo_cluster_terms(df, group_of_interest, GO_DEGs_df, "dn", path_name)
df_revigo_clustered <- rbind(df_up, df_mix_up, df_mix, df_mix_dn, df_dn)



# D.2 Complement Revigo GOs with info from all_dfs_GO
################################################################################################
df <- GO_DEGs_df

# Get info from original all_dfs_GO df
df_info <- df %>%
  filter(Comparison == group_of_interest & exactSource %in% df_revigo_clustered$go) %>% 
  dplyr::select(Name, exactSource, padj, Direction)

# Merge df:info with revigo result df
df_plot <- merge(df_info, df_revigo_clustered, by.x = "exactSource", by.y = "go", all.x = TRUE) %>%
  pivot_wider(names_from = Direction, values_from = padj) %>%
  dplyr::select(exactSource, Name, parentTerm, revigo_source, up, mix_up, mix_nondir, mix_dn, dn)

# Reformat GO Name term
df_plot$Name <- df_plot$Name %>% gsub(GO_XX, "", .) %>% gsub("_", " ", .) %>%  tolower() # REFORMAT GO stuff (remove leading GO_XX) 

# Add number of GOs in each cluster info
df_plot <- merge(df_plot, df_revigo_clustered, by.x = c("exactSource", "revigo_source"), by.y = c("go", "revigo_source")) %>%
  dplyr::select(exactSource, Name, revigo_source, gos_per_cluster, up, mix_up, mix_nondir, mix_dn, dn)

# Match to ensure same order before
df_plot <- df_plot[match(paste(df_revigo_clustered$go, df_revigo_clustered$revigo_source), paste(df_plot$exactSource, df_plot$revigo_source)), ]



# D.3 Get distinct GOs (the first from first group -> the representative GO, with lowest padj or most genes)& set factors
################################################################################################
df_plot$Name <- paste(df_plot$revigo_source, " ", df_plot$Name, sep = "")
df_plot$Name <- as.factor(df_plot$Name) # SET FACTORS
df_plot$Name <- factor(df_plot$Name, levels(unique(df_plot$Name))) # SET FACTORS




# D.4 Barplot of nr GOs per clsuter
################################################################################################

barplot <- ggplot(df_plot, aes(x=Name, y=gos_per_cluster, fill=revigo_source, level=Name)) + 
  geom_bar(stat = "identity", aes(x = factor(Name, level = Name))) +
  scale_fill_manual(values = c("dn" = "#6E7CB9FF", "mix_dn" = "#7BBCD5FF", "mix_nondir" = "#F8CD9CFF", "mix_up" = "#86C767FF", "up" = "#608F3DFF")) +
  theme_minimal() +
  ylab("# GOs") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))


pdf(file=paste0(path_name,"barplot.pdf"), width = 8, height = 1.5)
barplot
dev.off()


# D.5 Barplot of nr GOs per clsuter
################################################################################################

# Heatmap of representing GOs per cluster
plot_matrix <- df_plot %>% dplyr::select(-exactSource & -Name & -revigo_source & -gos_per_cluster) %>% as.matrix()
row.names(plot_matrix) <- df_plot$Name
plot_matrix <- t(plot_matrix)

d <- plot_matrix
d[d <= 0.01] <- as.character('**')
d[d > 0.01 & d <= 0.05] <- as.character('*')
d[d > 0.05] <- as.character(' ')


# Plot
plot_matrix_log <- apply(plot_matrix, 1:2, function(x) ifelse(is.na(x), NA, -log10(x)))

# Colors
color_palette <- colorRampPalette(c("white", "#E1595CFF"))(10)

annot_groups = data.frame(Groups=as.factor(df_plot$revigo_source))
row.names(annot_groups) <- colnames(plot_matrix_log)
ann_colors = list(Groups = c(dn="#6E7CB9FF", mix_dn="#7BBCD5FF", mix_nondir="#F8CD9CFF", mix_up="#86C767FF", up="#608F3DFF"))

#row.names(annot_groups) <- colnames(plot_matrix_log)

clusterTerms <- t(as.data.frame(plot_matrix_log))
clusterTerms <- data.frame(
  heatmap_nr = seq_len(nrow(clusterTerms)),
  term_name = row.names(clusterTerms))

#clusterTerms$heatmap_nr <- seq_len(nrow(clusterTerms))

write.csv(clusterTerms, file = paste0(path_name, "clusterTerms", group_of_interest, ".csv"), row.names = FALSE)


genes_heatmap <- pheatmap(plot_matrix_log,
                          display_numbers=d, # padj matrix (not logged)
                          number_color="black", # color of text in matrix
                          fontsize_number = 4, # size of text in matrix
                          color = color_palette, 
                          na_col = "white",
                          cluster_rows=F, cluster_cols = F,
                          annotation_col = annot_groups,
                          annotation_colors = ann_colors,
                          border_color = "grey92",
                          cex = 1,
                          #legend = TRUE,   # Show the legend
                          #labels_col = as.character(colnames(plot_matrix_log)),
                          angle_col = 90,
                          labels_col = 1:ncol(plot_matrix_log),
                          legend_width = 0.2,  # Adjust legend width (default is 0.15)
                          fontsize = 8,
                          fontsize_col = 6,
                          fontsize_row = 8,
                          cellwidth = 8, cellheight = 14)
#cellwidth = 8, cellheight = 14)

pdf(file=paste0(path_name,save_file), width = 16, height=5)
genes_heatmap
dev.off()



