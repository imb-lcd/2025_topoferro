library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(readxl)
library(cluster)


library(mgcv)


setwd("D:/Spatiotemporal_analysis/scRNA-seq")

anno_data <- read.table('jh_figure/gene_ferroptosis_label.txt', 
                  header = TRUE, sep = '\t', row.names = 1, stringsAsFactors = TRUE, na.string="")
anno_data <- anno_data[, 1:3]
anno_data$ox <- as.factor(anno_data$ox)

ferrdb_genelist <- rownames(anno_data[anno_data$db %in% c("driver", "suppressor"), ])

ann_colors <- list(
  reg = c("NRF2" = "#000000", "TEAD1" = "darkorange", "NA" = "#FFFFFF"),
  db = c("driver" = "darkorange", "suppressor" = "grey15", "related" = "#FFFFFF", "unknown" = "#FFFFFF"),
  ox = c("1" = "darkorange", "0" = "#FFFFFF")
)

padj_cut = 0.005
abc_cut_list = seq(50, 55, 1) #seq(20, 30, 1)
fold_cut_list = c(0.05, 0.1)
end_diff_cut_list = seq(0.5, 1.5, 0.1) #seq(0.3, 0.33, 0.01)


# selected
padj_cut = 0.005
abc_cut_list = 45 #seq(20, 30, 1)
fold_cut_list = 0.1
end_diff_cut_list = 0.65 #seq(0.3, 0.33, 0.01)

# k = 6;
degfree = 5;

sil_summary <- data.frame(
  file = character(), 
  mean_sil = numeric(),
  stringsAsFactors = FALSE
)


for (abc_cut in abc_cut_list) {
  for (fold_cut_range in fold_cut_list) {
    fold_cut = c(1-fold_cut_range, 1+fold_cut_range)
    for (end_diff_cut in end_diff_cut_list) {
      # for (genelist_sel in c(1, 2)) {
      for (genelist_sel in c(1, 1)) {
        de_gene = read_excel('jh_data/DE_statistics_all_quant.xlsx')
        
        abc = paste0('ABC', degfree)
        foldauc = paste0('foldAUC', degfree)
        end_diff = paste0('end_diff', degfree)
        
        de_gene <- de_gene[de_gene$ts_btwn_padj < padj_cut, ]
        de_gene <- de_gene[de_gene$`m_qval-65` == 1, ]
        
        de_gene <- de_gene[de_gene[[abc]] > abc_cut, ]
        de_gene <- de_gene[de_gene[[foldauc]] < fold_cut[1] | de_gene[[foldauc]] > fold_cut[2], ]
        de_gene <- de_gene[de_gene[[end_diff]] > end_diff_cut, ]
        
        de_gene_list = de_gene$Gene
        
        # ferrgene_list <- c('FTH1', 'XIST', 'CDKN1A', 'GPX4', 'SAT1', 'FTL', 'CYP1B1', 'TXNIP', 'BEX1', 'PRDX2', 'CCDC6', 'TGFB2', 'SMURF2', 'DNAJB6', 'CFL1', 'ACSL4', 'HSPA5', 'HSPB1', 'FSCN1')
        
        if (genelist_sel == 1) {
          genelist = de_gene_list
          outtype = paste0('df', degfree, '_m3_de', padj_cut, '_abc', abc_cut, '_auc', fold_cut_range, '_endiff', end_diff_cut, "_degene")
          k = 5
          print(paste(length(genelist), outtype))
        } else {
          print(paste(length(genelist), outtype))
          # genelist = ferrgene_list
          genelist= intersect(de_gene_list, ferrdb_genelist)
          outtype = paste0('df', degfree, '_m3_de', padj_cut, '_abc', abc_cut, '_auc', fold_cut_range, '_endiff', end_diff_cut, "_ferrgene")
          k = 2
        }

        # # load gene expression from upr and lwr traj
        upr <- readRDS(file = paste0('jh_data/all_upper_model_all_df', degfree, '_all.rds'))
        lwr <- readRDS(file = paste0('jh_data/all_lower_model_all_df', degfree, '_all.rds'))
      
        # upr <- readRDS(file = 'jh_data/upr_model_821_specificknot.rds')
        # lwr <- readRDS(file = 'jh_data/lwr_model_821_specificknot.rds')
        
        anno_data_curr <- anno_data
        
        upr <- upr[genelist, ]
        lwr <- lwr[genelist, ]
        anno_data_curr <- anno_data[genelist, ]
      
        # print(dim(upr))
        # print(dim(lwr))
  
        # set up pred data for pheatmap
        colnames(lwr) <- c(101:200)
        
        pred <- cbind(upr[, 100:1], lwr)
        
      
        pred <- (pred - rowMeans(pred)) / apply(pred, 1, sd)
        row_dist <- as.dist((1-cor(Matrix::t(pred)))/2)
        pred[pred>3] <- 3
        pred[pred<(-3)] <- (-3)
      
        # clust_method = list('ward.D', 'ward.D2') #, 'average', 'complete') #'single', 'complete', 'average', 'mcquitty', 'median', 'centroid')
        clust_method = list('ward.D') 
        
        for (cm in clust_method) {
          # print(cm)
          
          labrow <- rownames(lwr)
          labrow[ !(labrow %in% ferrdb_genelist) ] <- ''
      
          ph <- pheatmap::pheatmap(
            pred, 
            file = paste0('jh_figure/cluster_df', degfree, '_sel2/all_branched_',outtype,'_',cm,'_cl', k, '_heatmap.pdf'),
            cluster_cols = FALSE, 
            cluster_rows = TRUE, 
            show_colnames = F, 
            labels_row = labrow,
            clustering_distance_rows = row_dist, 
            clustering_method = cm,
            cutree_rows = k,
            silent = TRUE,
            annotation_row = anno_data,  # Add annotation
            annotation_colors = ann_colors, # Define colors for annotation
            color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(300),
            breaks = seq(-3, 3, length.out = 301)
          )
          
          tree <- ph$tree_row
          clusters <- cutree(tree, k = k)
          
          sil <- silhouette(clusters, row_dist)
          mean_sil <- mean(sil[, "sil_width"])
          
          # sil_summary <- rbind(sil_summary, data.frame(
          #   file = paste0(outtype,'_',cm),
          #   mean_sil = mean_sil
          # ))
          
          # show gene names in order of the pheatmap
          row_order <- ph$tree_row$order
          labrow_order <- labrow[row_order]
          
          # show the cluster it belongs
          treecuts <- cutree(ph$tree_row, k = k)
          treecuts <- treecuts[row_order]
          
          # reorder the anno_data and add cluster
          anno_data_reorder = data.frame(Cluster = treecuts, label = anno_data_curr[row_order, ])
          
          write.csv(anno_data_reorder, paste0('jh_figure/cluster_df', degfree, '_sel2/all_branched_',outtype,'_',cm,'_cl', k, '.csv'))
          
          group_labels <- c(rep("upper", 100), rep("lower", 100))
          names(group_labels) <- colnames(pred)

          cluster_df <- data.frame(Gene = names(treecuts), Cluster = treecuts)

          long_cluster_df <- as.data.frame(pred) %>%
            rownames_to_column("Gene") %>%
            pivot_longer(-Gene, names_to = "Cell", values_to = "Expression") %>%
            mutate(
              Group = group_labels[Cell],
              Cell = as.numeric(Cell)
            ) %>%
            left_join(cluster_df, by = "Gene")

          plot_data <- long_cluster_df %>%
            mutate(Cell_aligned = if_else(as.numeric(Cell) > 100,
                                          as.numeric(Cell) - 100,
                                          as.numeric(Cell))) %>%
            group_by(Cluster, Group, Cell_aligned) %>%
            summarise(
              mean_expr = mean(Expression, na.rm = TRUE),
              sd_expr = sd(Expression, na.rm = TRUE),
              .groups = "drop"
            ) %>%
            mutate(
              ymin = mean_expr - sd_expr,
              ymax = mean_expr + sd_expr
            )

          for (cl in seq(1,k)) {
            df <- filter(plot_data, Cluster == cl)

            p <- ggplot(df, aes(x = Cell_aligned, y = mean_expr, color = Group, fill = Group)) +
              geom_line(size = 2) +
              geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, linetype = 0) +
              labs(title = paste("Cluster", cl),
                   x = "Pseudotime",
                   y = "Mean Expression ± SD") +
              ylim(-3, 3) + 
              theme_minimal() +
              scale_color_manual(values = c(upper = "gray10", lower = "darkorange")) +
              scale_fill_manual(values = c(upper = "gray10", lower = "darkorange")) +
              theme(legend.position="none")

            print(p)

            ggsave(paste0("jh_figure/cluster_lineplots_sel2/cluster_", cl, "_upper_lower_plot.pdf"),
                   plot = p, width = 8, height = 5)
          }
          
        }
      }
    }
  }
}

write.csv(sil_summary, 'jh_figure/cluster_df5_tuning/sil_summary.csv')

## plot the mean and sd for each cluster

    # show gene names in order of the pheatmap
    row_order <- ph$tree_row$order
    labrow_order <- labrow[row_order]
    
    # show the cluster it belongs
    treecuts <- cutree(ph$tree_row, k = k)
    treecuts <- treecuts[row_order]
    
    # reorder the anno_data and add cluster
    anno_data_reorder = data.frame(Cluster = treecuts, label = anno_data_curr[row_order, ])
    
    write.csv(anno_data_reorder, paste0('jh_figure/all_branched_',outtype,'_',cm,'_cl', k, '.csv'))    
    

## for MSigDB

library(msigdbr)
library(clusterProfiler)
library(dplyr)

# ctype = 'C5_GOBPCCMF'
ctype = 'GOBP'

outpath = paste0("jh_figure/msigdb_df5_cl", k, "_", ctype, "_sel2")
dir.create(outpath, showWarnings = FALSE)

# Get C2 gene sets for Homo sapiens
msig <- msigdbr(species = "Homo sapiens", collection = 'C5')

# msig <- msig %>%
#   filter(gs_subcollection %in% c("GO:CC", "GO:BP", "GO:MF"))
# msig <- msig %>%
#   filter(gs_subcollection %in% c("GO:BP", "GO:MF"))
msig <- msig %>%
  filter(gs_subcollection %in% c("GO:BP"))

# Create TERM2GENE data frame
term2gene <- msig %>% select(gs_name, gene_symbol)

cluster_df <- data.frame(
  Gene = names(treecuts),
  Cluster = treecuts
)

# Initialize a list to store enrichment results
enrich_results <- list()

# Loop over each cluster
for (cl in sort(unique(cluster_df$Cluster))) {
# for (cl in c(1)) {
  
  genes_cl <- cluster_df %>% 
    filter(Cluster == cl) %>% 
    pull(Gene)
  
  # Run overrepresentation analysis
  enrich_res <- enricher(genes_cl, TERM2GENE = term2gene)
  
  # Filter for adjusted p-value < 0.05
  # enrich_res_filtered <- enrich_res@result %>%
  #   filter(p.adjust < 0.5) %>%
  #   arrange(p.adjust)
  enrich_res_filtered <- enrich_res@result %>%
    filter(qvalue < 0.05) %>%
    arrange(qvalue)
  
  # Save as .txt file
  write.table(enrich_res_filtered,
              file = paste0(outpath, "/cluster_", cl, "_enrichment.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Store in list
  # enrich_results[[paste0("cluster_", cl)]] <- enrich_res_filtered
  
  # Create horizontal barplot
  # if (nrow(enrich_res_filtered) > 0) {
    top_terms <- head(enrich_res_filtered, 20)  # Top 20 terms
    top_terms$Description <- factor(top_terms$Description,
                                    levels = rev(top_terms$Description))  # Reverse for plotting
    
    p <- ggplot(top_terms, aes(x = Description, y = -log10(p.adjust))) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(title = paste("MSigDB C2 Enrichment - Cluster", cl),
           x = "Enriched Terms", y = "-log10(q-value)") +
      theme_minimal(base_size = 12)
    
    ggsave(p, filename = paste0(outpath, "/cluster_", cl, "_barplot_top5.jpg"),
           width = 10, height = 6)
  # }
}

# supperessor in cluster 5
# 5 44 12 393
k = 5;
K = 12;
n = 44;
N = 393;
nrf2_pval <- phyper(q = k - 1, m = K, n = N - K, k = n, lower.tail = FALSE)
nrf2_pval