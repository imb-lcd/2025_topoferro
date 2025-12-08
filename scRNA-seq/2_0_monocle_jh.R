# calculate differential expression between two trajectories
if (FALSE)
{
  library(monocle3)
  cds <- readRDS('R_data/monocle3/cds_all_unorder_2')
  cds <- order_cells(cds, root_pr_nodes = c('Y_25'))
  
  # load upper trajectory
  cds_sub <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_25'), ending_pr_nodes = c('Y_64','Y_78'))
  cds_sub2 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_64'), ending_pr_nodes = c(paste0('Y_78')))
  final <- union(colnames(cds_sub),colnames(cds_sub2))
  print(head(final))
  cds_upr <- cds[,final]
  
  # load lower trajectory
  cds_sub <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_25'), ending_pr_nodes = c('Y_58','Y_82'))
  cds_sub2 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_58'), ending_pr_nodes = c(paste0('Y_82')))
  final <- union(colnames(cds_sub),colnames(cds_sub2))
  cds_lwr <- cds[,final]
  
  # extract the cells belonging to each branch
  upr_cells <- colnames(cds_upr)
  lwr_cells <- colnames(cds_lwr)
  
  # add a branch column in cds
  colData(cds)$branch <- NA
  
  colData(cds)$branch[colnames(cds) %in% upr_cells] <- "upr"
  colData(cds)$branch[colnames(cds) %in% lwr_cells] <- "lwr"
  
  # subset the cds into the two branches only
  cds_sub <- cds[, !is.na(colData(cds)$branch)]
  
  # perform differential epression
  # Fit the model
  gene_fits <- fit_models(cds_sub, model_formula_str = "~ branch")
  # saveRDS(gene_fits, file = 'figure_jh/gene_fits_uprlwr.rds')
  
  # Extract coefficient estimates and p-values
  fit_coefs <- coefficient_table(gene_fits)
  # saveRDS(fit_coefs, file = 'figure_jh/fit_coefs_uprlwr.rds')
  
  # Filter for genes with a significant difference between branch2 and branch1
  DE_genes <- subset(
    fit_coefs,
    term == "branchupr" & q_value < 0.05
  )
  
  head(DE_genes)
  saveRDS(DE_genes, file = 'figure_jh/DE_genes_uprlwr.rds')
  
  DE_genes <- apply(DE_genes, 2, as.character)
  write.csv(DE_genes,  'figure_jh/DE_genes_uprlwr.csv', row.names=FALSE)
}

# obtain cublic spline for the two trajecotires
if (FALSE)
{
  library(monocle3)
  cds <- readRDS('R_data/monocle3/cds_all_unorder_2')
  cds <- order_cells(cds, root_pr_nodes = c('Y_25'))
  
  # load upper trajectory
  cds_sub <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_25'), ending_pr_nodes = c('Y_64','Y_78'))
  cds_sub2 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_64'), ending_pr_nodes = c(paste0('Y_78')))
  final <- union(colnames(cds_sub),colnames(cds_sub2))
  print(head(final))
  cds_upr <- cds[,final]
  
  # load lower trajectory
  cds_sub <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_25'), ending_pr_nodes = c('Y_58','Y_82'))
  cds_sub2 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_58'), ending_pr_nodes = c(paste0('Y_82')))
  final <- union(colnames(cds_sub),colnames(cds_sub2))
  cds_lwr <- cds[,final]
  
  # load differential expression gene list
  DE_genes <- readRDS('figure_jh/DE_genes_uprlwr.rds')
  
  
  # extract gm smooth spline data for cds upr differential genes
  cds_upr_de <- cds_upr[DE_genes$gene_id, ]
  
  
  trend_formula = "~ splines::ns(pseudotime, df=5)"
  model_tbl <- fit_models(cds_upr_de, cores=1, model_formula_str = trend_formula)
  colData(cds_upr_de)$pseudotime <- monocle3::pseudotime(cds_upr_de)
  newdata <- data.frame(pseudotime = seq(0, max(colData(cds_upr_de)$pseudotime), length.out=100))
  model_expectation <- monocle3:::model_predictions(model_tbl, new_data = newdata)
  
  saveRDS(model_expectation, file = 'figure_jh/cds_upr_DE_genes_uprlwr.rds')
  
  
  # extract gm smooth spline data for cds lwr differential genes
  cds_lwr_de <- cds_lwr[DE_genes$gene_id, ]
  
  
  trend_formula = "~ splines::ns(pseudotime, df=5)"
  model_tbl <- fit_models(cds_lwr_de, cores=1, model_formula_str = trend_formula)
  colData(cds_lwr_de)$pseudotime <- monocle3::pseudotime(cds_lwr_de)
  newdata <- data.frame(pseudotime = seq(0, max(colData(cds_lwr_de)$pseudotime), length.out=100))
  model_expectation <- monocle3:::model_predictions(model_tbl, new_data = newdata)
  
  saveRDS(model_expectation, file = 'figure_jh/cds_lwr_DE_genes_uprlwr.rds')
}
