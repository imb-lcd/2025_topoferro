#adding lower RNA velocity
if(FALSE)
{
    gene <- readRDS(file = 'R_data/monocle3/pr_test_all_sub_upper_2')
    gene_nrf2 <- row.names(subset(gene, q_value < 10^(-40)))

    gene <- readRDS(file = 'R_data/monocle3/pr_test_all_sub_lower_2')
    gene_tead <- row.names(subset(gene, q_value < 10^(-40)))

    gene <- union(gene_nrf2,gene_tead)
    df <- as.data.frame(matrix(0,nrow = length(gene), ncol = 5))
    colnames(df) <- c('name', 'upper', 'lower', 'both', 'anno')
    df[,'name'] <- gene
    df[df[,'name'] %in% gene_nrf2,'upper'] <- 1
    df[df[,'name'] %in% gene_tead,'lower'] <- 1
    df[ (df[,'upper'] == 1)&(df[,'lower'] == 1) ,'both'] <- 1
    df[df[,'upper'] == 1,'anno'] <- 'upper'
    df[df[,'lower'] == 1,'anno'] <- 'lower'
    df[df[,'both'] == 1,'anno'] <- 'both'
    write.table( df, file = 'R_data/monocle3/gene_annotation.csv',quote = F,row.names = T,col.names = T,sep = ',')
    saveRDS(df, file = 'R_data/monocle3/gene_annotation')
}


#adding lower RNA velocity
if(FALSE)
{
    df <- readRDS(df, file = 'R_data/monocle3/gene_annotation')
    #divide <- '50'
    #all_exp <- read.csv( paste0('R_data/python_heatmap_diff_gene_',divide,'.csv'), header = TRUE, row.names = 1, sep = ',', check.names = FALSE)
    #print(dim(all_exp))
    regulon_link <- read.table(file = 'output/Step2_regulonTargetsInfo.tsv', sep = '\t', header = TRUE)
    regulon_link <- regulon_link[regulon_link[,'gene'] %in% df[,'name'],]
    print(dim(regulon_link))
    #print(head(regulon_link))
    library(dplyr)
    library(tibble)
    library(tidygraph)
    library(tidyr)
    library(ggraph)
    library(uwot)
    library(ggrepel)

    reg_mat <- regulon_link %>%
        as_tibble() %>%
        select(TF,gene,CoexWeight) %>%
        pivot_wider(names_from = gene, values_from = CoexWeight, values_fill = 0)%>%
        column_to_rownames('TF') %>% as.matrix()
    corr <- readRDS(file = 'R_data/monocle3/gene_corr_matrix')
    corr <- corr[rownames(reg_mat),colnames(reg_mat)] * (reg_mat*100+1)
    print(dim(reg_mat))
    #print(head(reg_mat))
    #if (ncol(reg_mat)>100){
    reg_mat[is.na(reg_mat)] <- 0
    if (TRUE){
        pca_mat <- prcomp(reg_mat, n=10)$x
        rownames(pca_mat) <- rownames(reg_mat)
        reg_mat <- as.matrix(pca_mat)
    }
    umap_tbl <- umap(reg_mat)
    #print(umap_tbl)
    umap_tbl <- as.data.frame(umap_tbl)
    colnames(umap_tbl) <- c('UMAP_1', 'UMAP_2')

    umap_tbl[,'name'] = rownames(umap_tbl)

    library(ggplot2)
	pdf("figure/TF_umap.pdf", width = 8, height = 7, onefile = TRUE)
    graph <- ggplot(umap_tbl, aes(y = UMAP_2, x = UMAP_1)) +
        geom_point( alpha = 0.6, size = 2) +
        geom_label_repel( aes(label = name))+
        #scale_fill_brewer(palette = 'Dark2', guide = FALSE)+
        labs()
    #graph <- graph+
    #    theme_classic()+
    #    theme(strip.background = element_blank(), 
    #          strip.text.x = element_blank(),
    #          axis.ticks.x = element_blank())
    graph <- graph+
        theme_classic()+
        theme(strip.background = element_blank())
    print(graph)
	dev.off()
}

if(FALSE)
{
    df <- readRDS(df, file = 'R_data/monocle3/gene_annotation')
    #divide <- '50'
    #all_exp <- read.csv( paste0('R_data/python_heatmap_diff_gene_',divide,'.csv'), header = TRUE, row.names = 1, sep = ',', check.names = FALSE)
    #print(dim(all_exp))
    regulon_link <- read.table(file = 'output/Step2_regulonTargetsInfo.tsv', sep = '\t', header = TRUE)
    regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% df[,'name']) | (regulon_link[,'gene'] %in% regulon_link[,'TF']),]
    print(dim(regulon_link))
    gene <- union( regulon_link[,'TF'], regulon_link[,'gene'])
    corrMatrix <- readRDS(file = 'int/1.2_corrMat.Rds')
    print(gene)
    corrMatrix <- corrMatrix[gene,gene]
    saveRDS(corrMatrix, file = 'R_data/monocle3/gene_corr_matrix')
    print(dim(corrMatrix))
}

if(FALSE)
{
    df <- readRDS(df, file = 'R_data/monocle3/gene_annotation')
    #divide <- '50'
    #all_exp <- read.csv( paste0('R_data/python_heatmap_diff_gene_',divide,'.csv'), header = TRUE, row.names = 1, sep = ',', check.names = FALSE)
    #print(dim(all_exp))
    regulon_link <- read.table(file = 'output/Step2_regulonTargetsInfo.tsv', sep = '\t', header = TRUE)
    regulon_link <- regulon_link[regulon_link[,'gene'] %in% df[,'name'],]
    print(dim(regulon_link))
    #print(head(regulon_link))
    library(dplyr)
    library(tibble)
    library(tidygraph)
    library(tidyr)
    library(ggraph)
    library(uwot)
    library(ggrepel)
    
    corr <- readRDS(file = 'R_data/monocle3/gene_corr_matrix')

    reg_mat <- regulon_link %>%
        as_tibble() %>%
        select(TF,gene,CoexWeight) %>%
        pivot_wider(names_from = TF, values_from = CoexWeight, values_fill = 0)%>%
        column_to_rownames('gene') %>% as.matrix()
    #reg_mat <- corr[rownames(reg_mat),colnames(reg_mat)] * (reg_mat*1000+1)
    reg_mat <- corr[rownames(reg_mat),colnames(reg_mat)] * (reg_mat*10+1)
    
    print(dim(reg_mat))

    reg_mat[is.na(reg_mat)] <- 0
    if (TRUE){
        pca_mat <- prcomp(reg_mat)$x[,1:10]
        rownames(pca_mat) <- rownames(reg_mat)
        reg_mat <- as.matrix(pca_mat)
    }

    set.seed(123)
    umap_tbl <- umap(reg_mat)
    #print(umap_tbl)
    umap_tbl <- as.data.frame(umap_tbl)
    colnames(umap_tbl) <- c('UMAP_1', 'UMAP_2')

    umap_tbl <- umap_tbl %>%
        as_tibble(rownames='gene')

    library(ggplot2)
	pdf("figure/TF_umap.pdf", width = 8, height = 7, onefile = TRUE)
    graph <- ggplot(umap_tbl, aes(y = UMAP_2, x = UMAP_1)) +
        geom_point( alpha = 0.6, size = 2) +
        #geom_label_repel( aes(label = name))+
        #scale_fill_brewer(palette = 'Dark2', guide = FALSE)+
        labs()
    #graph <- graph+
    #    theme_classic()+
    #    theme(strip.background = element_blank(), 
    #          strip.text.x = element_blank(),
    #          axis.ticks.x = element_blank())
    graph <- graph+
        theme_classic()+
        theme(strip.background = element_blank())
    print(graph)
	dev.off()
}



if(FALSE)
{
    .libPaths( c('~/.conda/pkgs/r-igraph-1.5.1-r43hbd8eb98_0/lib/R/',.libPaths()))
    df <- readRDS(df, file = 'R_data/monocle3/gene_annotation')
    #print(df)

    corr <- read.table(file = 'R_data/regulon_cor_python.csv', sep = ',', header = TRUE, row.names = 1)
    #print(head(regulon_link))
    regulon_link <- read.table(file = 'output/Step2_regulonTargetsInfo.tsv', sep = '\t', header = TRUE)
    #regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% df[,'name']) | (regulon_link[,'gene'] %in% regulon_link[,'TF']),]
    regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% colnames(corr)),]
    regulon_link <- regulon_link[ (regulon_link[,'TF'] %in% colnames(corr)),]
    print(dim(regulon_link))

    library(Matrix, lib='~/.conda/pkgs/r-matrix-1.5_4-r42he1ae0d6_0/lib/R/library')
    library(igraph, lib="~/.conda/pkgs/r-igraph-1.5.1-r43hbd8eb98_0/lib/R/library/")
    library(irlba, lib="~/.conda/pkgs/r-irlba-2.3.5.1-r43h316c678_1/lib/R/library/")
    library(tidygraph)
    library(dplyr)
    library(tibble)
    library(tidyr)
    library(ggraph)
    library(uwot)
    library(ggrepel)
    library(ggplot2, lib = '~/.conda/pkgs/r-ggplot2-3.4.0-r42hc72bb7e_1/lib/R/library')
    pdf("figure/TF_umap_test_0.pdf", width = 8, height = 7, onefile = TRUE)
    
    #for( i in c(100,105,110,115,116,117,118,119,120,121,122,123,124,125,126))
    #for( i in c(2,5,10,20,25,30,35,40,45,50,60,70,75,80,100,125,150))
    for( i in c(30,31))
    {
        #corr <- readRDS( file = 'R_data/monocle3/gene_corr_matrix')

        reg_mat <- regulon_link %>%
            as_tibble() %>%
            select(TF,gene,CoexWeight) %>%
            pivot_wider(names_from = TF, values_from = CoexWeight, values_fill = 0)%>%
            column_to_rownames('gene') %>% as.matrix()
        reg_mat[is.na(reg_mat)] <- 0
        #reg_mat <- corr[rownames(reg_mat),colnames(reg_mat)] * (reg_mat*0+1)
        reg_mat <- corr[rownames(reg_mat),colnames(reg_mat)]
        
        print(dim(reg_mat))

        
        if (TRUE){
            pca_mat <- prcomp(reg_mat)$x[,1:i]
            rownames(pca_mat) <- rownames(reg_mat)
            reg_mat <- as.matrix(pca_mat)
        }

        set.seed(123)
        umap_tbl <- umap(reg_mat)
        umap_tbl <- as.data.frame(umap_tbl)
        colnames(umap_tbl) <- c('UMAP_1', 'UMAP_2')

        umap_tbl <- umap_tbl %>%
            as_tibble(rownames='gene')
        
        if( FALSE)
        {
        gene_graph <- as_tbl_graph(gene_net) %>%
            activate(edges) %>%
            mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
            activate(nodes) %>%
            mutate(centrality=centrality_pagerank()) %>%
            inner_join(umap_tbl, by=c('name'='gene'))

        print(gene_graph)
        }

        umap_tbl$name <- umap_tbl$gene
        genelist <- c('TEAD1','NFE2L2','GPX4','SQSTM1','ACSL4','FTL','FTH1','E2F1')
        #umap_tbl$name[ !(umap_tbl$gene %in% regulon_link[,'TF'])] <- ''
        umap_tbl$name[ !(umap_tbl$gene %in% genelist)] <- ''

        print(i)
        graph <- ggplot(umap_tbl , aes(y = UMAP_2, x = -UMAP_1)) +
            geom_point( alpha = 0.6, size = 2) +
            geom_text_repel( aes(label = name), max.overlaps=99999)+
            #scale_fill_brewer(palette = 'Dark2', guide = FALSE)+
            labs(title = toString(i))
        #graph <- graph+
        #    theme_classic()+
        #    theme(strip.background = element_blank(), 
        #          strip.text.x = element_blank(),
        #          axis.ticks.x = element_blank())
        graph <- graph+
            theme_classic()+
            theme(strip.background = element_blank())
        print(graph)
    }
	dev.off()
}

if(FALSE)
{
    df <- readRDS(df, file = 'R_data/monocle3/gene_annotation')
    regulon_link <- read.table(file = 'output/Step2_regulonTargetsInfo.tsv', sep = '\t', header = TRUE)
    regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% df[,'name']) | (regulon_link[,'gene'] %in% regulon_link[,'TF']),]
    print(dim(regulon_link))

    library(dplyr)
    library(tibble)
    library(tidygraph)
    library(tidyr)
    library(ggraph)
    library(uwot)
    library(ggrepel)
    
    gene_net <- regulon_link %>%
            select(TF, gene, everything()) %>%
            group_by(gene)
    print(gene_net)
    
    #corr <- readRDS(file = 'int/1.2_corrMat.Rds')
    corr <- readRDS( file = 'R_data/monocle3/gene_corr_matrix')

    reg_mat <- regulon_link %>%
        as_tibble() %>%
        select(TF,gene,CoexWeight) %>%
        pivot_wider(names_from = TF, values_from = CoexWeight, values_fill = 0)%>%
        column_to_rownames('gene') %>% as.matrix()
    #reg_mat <- corr[rownames(reg_mat),colnames(reg_mat)] * (reg_mat*1000+1)
    reg_mat[is.na(reg_mat)] <- 0
    print(dim(corr))
    print(dim(reg_mat))
    reg_mat <- corr[rownames(reg_mat),colnames(reg_mat)] * (reg_mat*120+1)
    
    print(dim(reg_mat))

    
    if (TRUE){
        pca_mat <- prcomp(reg_mat)$x[,1:10]
        rownames(pca_mat) <- rownames(reg_mat)
        reg_mat <- as.matrix(pca_mat)
    }

    set.seed(123)
    umap_tbl <- umap(reg_mat)
    umap_tbl <- as.data.frame(umap_tbl)
    colnames(umap_tbl) <- c('UMAP_1', 'UMAP_2')

    umap_tbl <- umap_tbl %>%
        as_tibble(rownames='gene')
    
    if( FALSE)
    {
    gene_graph <- as_tbl_graph(gene_net) %>%
        activate(edges) %>%
        mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
        activate(nodes) %>%
        mutate(centrality=centrality_pagerank()) %>%
        inner_join(umap_tbl, by=c('name'='gene'))

    print(gene_graph)
    }

    umap_tbl$name <- umap_tbl$gene
    genelist <- c('TEAD1','NFE2L2','GPX4','SQSTM1','ACSL4','FTL','FTH1','E2F1')
    #umap_tbl$name[ !(umap_tbl$gene %in% regulon_link[,'TF'])] <- ''
    umap_tbl$name[ !(umap_tbl$gene %in% genelist)] <- ''

    library(ggplot2)
	pdf("figure/TF_umap.pdf", width = 8, height = 7, onefile = TRUE)
    graph <- ggplot(umap_tbl , aes(y = UMAP_2, x = UMAP_1)) +
        geom_point( alpha = 0.6, size = 2) +
        geom_text_repel( aes(label = name), max.overlaps=99999)+
        #scale_fill_brewer(palette = 'Dark2', guide = FALSE)+
        labs()
    #graph <- graph+
    #    theme_classic()+
    #    theme(strip.background = element_blank(), 
    #          strip.text.x = element_blank(),
    #          axis.ticks.x = element_blank())
    graph <- graph+
        theme_classic()+
        theme(strip.background = element_blank())
    print(graph)
	dev.off()
}

if(FALSE)
{
    df <- readRDS(df, file = 'R_data/monocle3/gene_annotation')
    regulon_link <- read.table(file = 'output/Step2_regulonTargetsInfo.tsv', sep = '\t', header = TRUE)
    regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% df[,'name']) | (regulon_link[,'gene'] %in% regulon_link[,'TF']),]
    print(dim(regulon_link))

    library(dplyr)
    library(tibble)
    library(tidygraph)
    library(tidyr)
    library(ggraph)
    library(uwot)
    library(ggrepel)
    library(ggplot2)
    pdf("figure/TF_umap.pdf", width = 8, height = 7, onefile = TRUE)
    
    #for( i in c(100,105,110,115,116,117,118,119,120,121,122,123,124,125,126))
    for( i in c(1:200))
    {
        corr <- readRDS( file = 'R_data/monocle3/gene_corr_matrix')

        reg_mat <- regulon_link %>%
            as_tibble() %>%
            select(TF,gene,CoexWeight) %>%
            pivot_wider(names_from = TF, values_from = CoexWeight, values_fill = 0)%>%
            column_to_rownames('gene') %>% as.matrix()
        reg_mat[is.na(reg_mat)] <- 0
        reg_mat <- corr[rownames(reg_mat),colnames(reg_mat)] * (reg_mat*i+1)
        
        print(dim(reg_mat))

        
        if (TRUE){
            pca_mat <- prcomp(reg_mat)$x[,1:10]
            rownames(pca_mat) <- rownames(reg_mat)
            reg_mat <- as.matrix(pca_mat)
        }

        set.seed(123)
        umap_tbl <- umap(reg_mat)
        umap_tbl <- as.data.frame(umap_tbl)
        colnames(umap_tbl) <- c('UMAP_1', 'UMAP_2')

        umap_tbl <- umap_tbl %>%
            as_tibble(rownames='gene')
        
        if( FALSE)
        {
        gene_graph <- as_tbl_graph(gene_net) %>%
            activate(edges) %>%
            mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
            activate(nodes) %>%
            mutate(centrality=centrality_pagerank()) %>%
            inner_join(umap_tbl, by=c('name'='gene'))

        print(gene_graph)
        }

        umap_tbl$name <- umap_tbl$gene
        genelist <- c('TEAD1','NFE2L2','GPX4','SQSTM1','ACSL4','FTL','FTH1','E2F1')
        #umap_tbl$name[ !(umap_tbl$gene %in% regulon_link[,'TF'])] <- ''
        umap_tbl$name[ !(umap_tbl$gene %in% genelist)] <- ''

        print(i)
        graph <- ggplot(umap_tbl , aes(y = -UMAP_2, x = -UMAP_1)) +
            geom_point( alpha = 0.6, size = 2) +
            geom_text_repel( aes(label = name), max.overlaps=99999)+
            #scale_fill_brewer(palette = 'Dark2', guide = FALSE)+
            labs(title = toString(i))
        #graph <- graph+
        #    theme_classic()+
        #    theme(strip.background = element_blank(), 
        #          strip.text.x = element_blank(),
        #          axis.ticks.x = element_blank())
        graph <- graph+
            theme_classic()+
            theme(strip.background = element_blank())
        print(graph)
    }
	dev.off()
}

if(FALSE)
{
    df <- readRDS(df, file = 'R_data/monocle3/gene_annotation')

    write.table( df,file = 'R_data/monocle3/gene_annotation.csv',quote = F,row.names = T,col.names = T,sep = '\t')
}


if(FALSE)
{
    df <- readRDS(df, file = 'R_data/monocle3/gene_annotation')
    regulon_link <- read.table(file = 'output/Step2_regulonTargetsInfo.tsv', sep = '\t', header = TRUE)
    #regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% df[,'name']) | (regulon_link[,'gene'] %in% regulon_link[,'TF']),]
    #regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% regulon_link[,'TF']),]
    regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% df[,'name']) ,]
    print(dim(regulon_link))

	tf_ratio <- read.table('python_data/tf_ratio.csv',sep = ',', header = TRUE, row.names = 1)
    print(head(tf_ratio))
    library(dplyr)
    library(tibble)
    library(tidygraph)
    library(tidyr)
    library(ggraph)
    library(uwot)
    library(ggrepel)
    library(ggplot2)
    pdf("figure/TF_umap_only_TF_only_TF_noPCA.pdf", width = 8, height = 7, onefile = TRUE)
    
    #for( i in c(100,105,110,115,116,117,118,119,120,121,122,123,124,125,126))
    #for( i in c(0,1,25,50,100,200,500,1000,2000))
    for( i in c(0))
    {
        corr <- readRDS( file = 'R_data/monocle3/gene_corr_matrix')

        reg_mat <- regulon_link %>%
            as_tibble() %>%
            select(TF,gene,CoexWeight) %>%
            pivot_wider(names_from = gene, values_from = CoexWeight, values_fill = 0)%>%
            column_to_rownames('TF') %>% as.matrix()
        reg_mat[is.na(reg_mat)] <- 0
        #reg_mat <- corr[rownames(reg_mat),colnames(reg_mat)] * (reg_mat*i+1)
        
        print(dim(reg_mat))

        
        if (FALSE){
            pca_mat <- prcomp(reg_mat)$x[,1:10]
            rownames(pca_mat) <- rownames(reg_mat)
            reg_mat <- as.matrix(pca_mat)
        }

        set.seed(123)
        umap_tbl <- umap(reg_mat)
        umap_tbl <- as.data.frame(umap_tbl)
        colnames(umap_tbl) <- c('UMAP_1', 'UMAP_2')

        umap_tbl <- umap_tbl %>%
            as_tibble(rownames='TF')
        
        if( FALSE)
        {
        gene_graph <- as_tbl_graph(gene_net) %>%
            activate(edges) %>%
            mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
            activate(nodes) %>%
            mutate(centrality=centrality_pagerank()) %>%
            inner_join(umap_tbl, by=c('name'='gene'))

        print(gene_graph)
        }

        umap_tbl$name <- umap_tbl$TF
        rownames(umap_tbl) <- umap_tbl$TF
        genelist <- c('TEAD1','NFE2L2','E2F1')
        #umap_tbl$name[ !(umap_tbl$gene %in% regulon_link[,'TF'])] <- ''
        umap_tbl$name[ !(umap_tbl$TF %in% genelist)] <- ''
        umap_tbl <- umap_tbl[umap_tbl$TF %in% rownames(tf_ratio),]
        tf_ratio[,'TF'] <- rownames(tf_ratio)
        tf_ratio <- tf_ratio[,c('TF','ratio')]
        umap_tbl <- umap_tbl %>%
            select(TF, everything()) %>%
            group_by(TF) %>%
            left_join(tf_ratio, by=c('TF'='TF'))
        print(umap_tbl)

        print(i)
        graph <- ggplot(umap_tbl , aes(y = UMAP_2, x = UMAP_1)) +
            geom_point( alpha = 0.6, aes(size = ratio)) +
            geom_text_repel( aes(label = name), max.overlaps=99999)+
            #scale_fill_brewer(palette = 'Dark2', guide = FALSE)+
            labs(title = toString(i))
        #graph <- graph+
        #    theme_classic()+
        #    theme(strip.background = element_blank(), 
        #          strip.text.x = element_blank(),
        #          axis.ticks.x = element_blank())
        graph <- graph+
            theme_classic()+
            theme(strip.background = element_blank())
        print(graph)
    }
	dev.off()
}


if(FALSE)
{
    df <- readRDS(df, file = 'R_data/monocle3/gene_annotation')
    #print(df)

    corr <- read.table(file = 'R_data/regulon_cor_python.csv', sep = ',', header = TRUE, row.names = 1)
    #print(head(regulon_link))
    regulon_link <- read.table(file = 'output/Step2_regulonTargetsInfo.tsv', sep = '\t', header = TRUE)
    #regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% df[,'name']) | (regulon_link[,'gene'] %in% regulon_link[,'TF']),]

    regulon_link_tmp <- regulon_link[ (regulon_link[,'gene'] %in% df[,'name']) ,]
    regulon_link_TF_tmp <- regulon_link_tmp[ (regulon_link_tmp[,'gene'] %in% regulon_link_tmp[,'TF']) ,]

    regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% colnames(corr)),]
    regulon_link <- regulon_link[ (regulon_link[,'TF'] %in% colnames(corr)),]
    regulon_link <- regulon_link[ (regulon_link[,'TF'] %in% regulon_link_TF_tmp[,'TF']),]
    regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% regulon_link_TF_tmp[,'TF']),]
    regulon_link_TF <- regulon_link
    print(dim(regulon_link))

    library(Matrix, lib='~/.conda/pkgs/r-matrix-1.5_4-r42he1ae0d6_0/lib/R/library')
    library(igraph, lib="~/.conda/pkgs/r-igraph-1.5.1-r43hbd8eb98_0/lib/R/library/")
    library(irlba, lib="~/.conda/pkgs/r-irlba-2.3.5.1-r43h316c678_1/lib/R/library/")
    library(tidygraph)
    library(dplyr)
    library(tibble)
    library(tidyr)
    library(ggraph)
    library(uwot)
    library(ggrepel)
    library(ggplot2, lib = '~/.conda/pkgs/r-ggplot2-3.4.0-r42hc72bb7e_1/lib/R/library')
    pdf("figure/TF_umap_test_0_edge.pdf", width = 8, height = 7, onefile = TRUE)
    
    #for( i in c(100,105,110,115,116,117,118,119,120,121,122,123,124,125,126))
    #for( i in c(2,5,10,20,25,30,35,40,45,50,60,70,75,80,100,125,150))
    #for( i in c(30,31))
    for( i in c(2,5,6,10,11,20,21,25,26,30,31,35,36,40,41,45,46,50,51,60,61,70,71,75))
    {
        #corr <- readRDS( file = 'R_data/monocle3/gene_corr_matrix')
        gene_net <- regulon_link_TF %>%
                select(TF, gene, everything()) %>%
                group_by(TF)
        print(gene_net)
        reg_mat <- regulon_link %>%
            as_tibble() %>%
            select(TF,gene,CoexWeight) %>%
            pivot_wider(names_from = TF, values_from = CoexWeight, values_fill = 0)%>%
            column_to_rownames('gene') %>% as.matrix()
        reg_mat[is.na(reg_mat)] <- 0
        #reg_mat <- corr[rownames(reg_mat),colnames(reg_mat)] * (reg_mat*0.1+1)
        reg_mat <- corr[rownames(reg_mat),colnames(reg_mat)]
        
        print(dim(reg_mat))

        
        if (TRUE){
            pca_mat <- prcomp(reg_mat)$x[,1:i]
            rownames(pca_mat) <- rownames(reg_mat)
            reg_mat <- as.matrix(pca_mat)
        }

        set.seed(123)
        umap_tbl <- umap(reg_mat)
        umap_tbl <- as.data.frame(umap_tbl)
        colnames(umap_tbl) <- c('UMAP_1', 'UMAP_2')

        umap_tbl <- umap_tbl %>%
            as_tibble(rownames='TF')
        
        tf_ratio <- read.table('python_data/tf_ratio_anno.csv',sep = ',', header = TRUE, row.names = 1)
        print(head(tf_ratio))
        umap_tbl$name <- umap_tbl$TF
        rownames(umap_tbl) <- umap_tbl$TF
        umap_tbl <- umap_tbl[umap_tbl$TF %in% rownames(tf_ratio),]
        tf_ratio[,'TF'] <- rownames(tf_ratio)
        print('setting umap')
        tf_ratio <- tf_ratio[,c('TF','ratio','ratio_u','ratio_l','weight_time','alpha_u','alpha_l','u_name','l_name')]
        umap_tbl <- umap_tbl %>%
            select(TF, everything()) %>%
            group_by(TF) %>%
            left_join(tf_ratio, by=c('TF'='TF'))
        print('incorporating ratio info')
        if( TRUE)
        {
            gene_graph <- as_tbl_graph(gene_net) %>%
                activate(edges) %>%
                mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
                activate(nodes) %>%
                mutate(centrality=centrality_pagerank()) %>%
                inner_join(umap_tbl, by=c('name'='TF'))

            print(gene_graph)
        }
        p <- ggraph(gene_graph, x=UMAP_1, y=UMAP_2)
        p <- p + geom_edge_diagonal(alpha = 0.05)
        #p <- p + geom_node_point(aes(fill=centrality, size=ratio), color='darkgrey', shape=21)
        p <- p + geom_node_point(aes(fill=weight_time, size=ratio), color='darkgrey', shape=21)
        p <- p + scale_fill_viridis(option = "B", end = 0.85)
        p <- p + geom_node_text(
            aes(label=name),
            repel=T, max.overlaps=99999
        )
        p <- p +
            labs(title = toString(i))+
            theme_void() 
        print(p)
        p <- ggraph(gene_graph, x=UMAP_1, y=UMAP_2)
        p <- p + geom_edge_diagonal(alpha = 0.05)
        #p <- p + geom_node_point(aes(fill=centrality, size=ratio), color='darkgrey', shape=21)
        p <- p + geom_node_point(aes(alpha=alpha_u, size=ratio), fill = 'darkorange', color='darkgrey', shape=21)
        #p <- p + scale_fill_gradient(low = '#FFD580', high = 'orange')
        #p <- p + scale_fill_viridis(option = "B", end = 0.85)
        p <- p + geom_node_text(
            aes(label=u_name),
            repel=T, max.overlaps=99999
        )
        p <- p +
            theme_void() 
        print(p)
        p <- ggraph(gene_graph, x=UMAP_1, y=UMAP_2)
        p <- p + geom_edge_diagonal(alpha = 0.05)
        #p <- p + geom_node_point(aes(fill=centrality, size=ratio), color='darkgrey', shape=21)
        p <- p + geom_node_point(aes(alpha=alpha_l, size=ratio), fill = 'purple', color='darkgrey', shape=21)
        #p <- p + scale_fill_gradient(low = 'light purple', high = 'purple')
        #p <- p + scale_fill_viridis(option = "B", end = 0.85)
        p <- p + geom_node_text(
            aes(label=l_name),
            repel=T, max.overlaps=99999
        )
        p <- p +
            theme_void() 
        print(p)
        if( FALSE)
        {
        gene_graph <- as_tbl_graph(gene_net) %>%
            activate(edges) %>%
            mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
            activate(nodes) %>%
            mutate(centrality=centrality_pagerank()) %>%
            inner_join(umap_tbl, by=c('name'='gene'))

        print(gene_graph)
        }

    }
	dev.off()
}


if(FALSE)
{
    df <- readRDS(df, file = 'R_data/monocle3/gene_annotation')
    regulon_link <- read.table(file = 'output/Step2_regulonTargetsInfo.tsv', sep = '\t', header = TRUE)
    #regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% df[,'name']) | (regulon_link[,'gene'] %in% regulon_link[,'TF']),]
    #regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% regulon_link[,'TF']),]
    regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% df[,'name']) ,]
    regulon_link_TF <- regulon_link[ (regulon_link[,'gene'] %in% regulon_link[,'TF']) ,]
    print(dim(regulon_link))

	tf_ratio <- read.table('python_data/tf_ratio_anno.csv',sep = ',', header = TRUE, row.names = 1)
    print(head(tf_ratio))
    library(dplyr)
    library(tibble)
    library(tidygraph)
    library(tidyr)
    library(ggraph)
    library(uwot)
    library(ggrepel)
    library(ggplot2)
    library(viridis)
    pdf("figure/TF_umap_only_TF_only_TF_noPCA_graph_sign.pdf", width = 8, height = 7, onefile = TRUE)
    
    #for( i in c(100,105,110,115,116,117,118,119,120,121,122,123,124,125,126))
    #for( i in c(0,1,25,50,100,200,500,1000,2000))
    for( i in c(10,15,20,25,50,100,150))
    {
        corr <- readRDS( file = 'R_data/monocle3/gene_corr_matrix')

        #corr <- readRDS(file = 'int/1.2_corrMat.Rds')
        gene_net <- regulon_link_TF %>%
                select(TF, gene, everything()) %>%
                group_by(TF)
        print(gene_net)

        reg_mat <- regulon_link %>%
            as_tibble() %>%
            select(TF,gene,CoexWeight) %>%
            pivot_wider(names_from = gene, values_from = CoexWeight, values_fill = 0)%>%
            column_to_rownames('TF') %>% as.matrix()
        reg_mat[is.na(reg_mat)] <- 0
        #print(head(sign(corr[rownames(reg_mat),colnames(reg_mat)])))
        #reg_mat <- sign(corr[rownames(reg_mat),colnames(reg_mat)]) * (reg_mat)
        
        print(dim(reg_mat))

        
        if (TRUE){
            print(i)
            pca_mat <- prcomp(reg_mat)$x[,1:i]
            rownames(pca_mat) <- rownames(reg_mat)
            reg_mat <- as.matrix(pca_mat)
        }

        set.seed(123)
        umap_tbl <- umap(reg_mat)
        umap_tbl <- as.data.frame(umap_tbl)
        colnames(umap_tbl) <- c('UMAP_1', 'UMAP_2')

        umap_tbl <- umap_tbl %>%
            as_tibble(rownames='TF')
        
        

        umap_tbl$name <- umap_tbl$TF
        rownames(umap_tbl) <- umap_tbl$TF
        genelist <- c('TEAD1','NFE2L2','E2F1')
        #umap_tbl$name[ !(umap_tbl$gene %in% regulon_link[,'TF'])] <- ''
        umap_tbl$name[ !(umap_tbl$TF %in% genelist)] <- ''
        umap_tbl <- umap_tbl[umap_tbl$TF %in% rownames(tf_ratio),]
        tf_ratio[,'TF'] <- rownames(tf_ratio)
        tf_ratio <- tf_ratio[,c('TF','ratio','ratio_u','ratio_l','weight_time','alpha_u','alpha_l','u_name','l_name')]
        umap_tbl <- umap_tbl %>%
            select(TF, everything()) %>%
            group_by(TF) %>%
            left_join(tf_ratio, by=c('TF'='TF'))
        print(umap_tbl)
        if( TRUE)
        {
            gene_graph <- as_tbl_graph(gene_net) %>%
                activate(edges) %>%
                mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
                activate(nodes) %>%
                mutate(centrality=centrality_pagerank()) %>%
                inner_join(umap_tbl, by=c('name'='TF'))

            print(gene_graph)
        }
        p <- ggraph(gene_graph, x=UMAP_1, y=UMAP_2)
        p <- p + geom_edge_diagonal(alpha = 0.05)
        #p <- p + geom_node_point(aes(fill=centrality, size=ratio), color='darkgrey', shape=21)
        p <- p + geom_node_point(aes(fill=weight_time, size=ratio), color='darkgrey', shape=21)
        p <- p + scale_fill_viridis(option = "B", end = 0.85)
        p <- p + geom_node_text(
            aes(label=name),
            repel=T, max.overlaps=99999
        )
        p <- p +
            theme_void() 
        print(p)
        p <- ggraph(gene_graph, x=UMAP_1, y=UMAP_2)
        p <- p + geom_edge_diagonal(alpha = 0.05)
        #p <- p + geom_node_point(aes(fill=centrality, size=ratio), color='darkgrey', shape=21)
        p <- p + geom_node_point(aes(alpha=alpha_u, size=ratio), fill = 'darkorange', color='darkgrey', shape=21)
        #p <- p + scale_fill_gradient(low = '#FFD580', high = 'orange')
        #p <- p + scale_fill_viridis(option = "B", end = 0.85)
        p <- p + geom_node_text(
            aes(label=u_name),
            repel=T, max.overlaps=99999
        )
        p <- p +
            theme_void() 
        print(p)
        p <- ggraph(gene_graph, x=UMAP_1, y=UMAP_2)
        p <- p + geom_edge_diagonal(alpha = 0.05)
        #p <- p + geom_node_point(aes(fill=centrality, size=ratio), color='darkgrey', shape=21)
        p <- p + geom_node_point(aes(alpha=alpha_l, size=ratio), fill = 'purple', color='darkgrey', shape=21)
        #p <- p + scale_fill_gradient(low = 'light purple', high = 'purple')
        #p <- p + scale_fill_viridis(option = "B", end = 0.85)
        p <- p + geom_node_text(
            aes(label=l_name),
            repel=T, max.overlaps=99999
        )
        p <- p +
            theme_void() 
        print(p)
    }
	dev.off()
}
if(FALSE)
{
    auc <- readRDS('working_data/regulon_auc.rds')
    write.table( auc,file = 'R_data/regulon_auc_for_python.csv',quote = F,row.names = T,col.names = T,sep = '\t')
}

if(FALSE)
{
    df <- readRDS(df, file = 'R_data/monocle3/gene_annotation')
    #print(df)

    corr <- read.table(file = 'R_data/regulon_cor_python.csv', sep = ',', header = TRUE, row.names = 1)
    #print(head(regulon_link))
    regulon_link <- read.table(file = 'output/Step2_regulonTargetsInfo.tsv', sep = '\t', header = TRUE)
    #regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% df[,'name']) | (regulon_link[,'gene'] %in% regulon_link[,'TF']),]

    regulon_link_tmp <- regulon_link[ (regulon_link[,'gene'] %in% df[,'name']) ,]
    regulon_link_TF_tmp <- regulon_link_tmp[ (regulon_link_tmp[,'gene'] %in% regulon_link_tmp[,'TF']) ,]

    regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% colnames(corr)),]
    regulon_link <- regulon_link[ (regulon_link[,'TF'] %in% colnames(corr)),]
    regulon_link <- regulon_link[ (regulon_link[,'TF'] %in% regulon_link_TF_tmp[,'TF']),]
    regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% regulon_link_TF_tmp[,'TF']),]
    regulon_link_TF <- regulon_link
    print(dim(regulon_link))

    library(Matrix, lib='~/.conda/pkgs/r-matrix-1.5_4-r42he1ae0d6_0/lib/R/library')
    library(igraph, lib="~/.conda/pkgs/r-igraph-1.5.1-r43hbd8eb98_0/lib/R/library/")
    library(irlba, lib="~/.conda/pkgs/r-irlba-2.3.5.1-r43h316c678_1/lib/R/library/")
    library(tidygraph)
    library(dplyr)
    library(tibble)
    library(tidyr)
    library(ggraph)
    library(uwot)
    library(ggrepel)
    library(ggplot2, lib = '~/.conda/pkgs/r-ggplot2-3.4.0-r42hc72bb7e_1/lib/R/library')
    pdf("figure/TF_umap_day2_3_final_opt_1.pdf", width = 8, height = 7, onefile = TRUE)
    
    #for( i in c(100,105,110,115,116,117,118,119,120,121,122,123,124,125,126))
    #for( i in c(2,5,10,20,25,30,35,40,45,50,60,70,75,80,100,125,150))
    #for( i in c(30,31))
    for( i in c(35))
    {
        #corr <- readRDS( file = 'R_data/monocle3/gene_corr_matrix')
        gene_net <- regulon_link_TF %>%
                select(TF, gene, everything()) %>%
                group_by(TF)
        print(gene_net)
        reg_mat <- regulon_link %>%
            as_tibble() %>%
            select(TF,gene,CoexWeight) %>%
            pivot_wider(names_from = TF, values_from = CoexWeight, values_fill = 0)%>%
            column_to_rownames('gene') %>% as.matrix()
        reg_mat[is.na(reg_mat)] <- 0
        #reg_mat <- corr[rownames(reg_mat),colnames(reg_mat)] * (reg_mat*0.1+1)
        reg_mat <- corr[rownames(reg_mat),colnames(reg_mat)]
        
        print(dim(reg_mat))

        
        if (TRUE){
            pca_mat <- prcomp(reg_mat)$x[,1:i]
            rownames(pca_mat) <- rownames(reg_mat)
            reg_mat <- as.matrix(pca_mat)
        }

        set.seed(123)
        umap_tbl <- umap(reg_mat)
        umap_tbl <- as.data.frame(umap_tbl)
        colnames(umap_tbl) <- c('UMAP_1', 'UMAP_2')

        umap_tbl <- umap_tbl %>%
            as_tibble(rownames='TF')
        
        tf_ratio <- read.table('python_data/tf_ratio_anno.csv',sep = ',', header = TRUE, row.names = 1)
        print(head(tf_ratio))
        umap_tbl$name <- umap_tbl$TF
        rownames(umap_tbl) <- umap_tbl$TF
        umap_tbl <- umap_tbl[umap_tbl$TF %in% rownames(tf_ratio),]
        tf_ratio[,'TF'] <- rownames(tf_ratio)
        print('setting umap')
        tf_ratio <- tf_ratio[,c('TF','ratio','ratio_u','ratio_l','weight_time','alpha_u','alpha_l','u_name','l_name')]
        umap_tbl <- umap_tbl %>%
            select(TF, everything()) %>%
            group_by(TF) %>%
            left_join(tf_ratio, by=c('TF'='TF'))
        print('incorporating ratio info')
        if( TRUE)
        {
            gene_graph <- as_tbl_graph(gene_net) %>%
                activate(edges) %>%
                mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
                activate(nodes) %>%
                mutate(centrality=centrality_pagerank()) %>%
                inner_join(umap_tbl, by=c('name'='TF'))

            print(gene_graph)
        }
        p <- ggraph(gene_graph, x=UMAP_1, y=-UMAP_2)
        p <- p + geom_edge_diagonal(alpha = 0.05)
        #p <- p + geom_node_point(aes(fill=centrality, size=ratio), color='darkgrey', shape=21)
        p <- p + geom_node_point(aes(fill=weight_time, size=ratio), color='darkgrey', shape=21)
        p <- p + scale_fill_viridis(option = "B", end = 0.85)
        p <- p + geom_node_text(
            aes(label=name),
            repel=T, max.overlaps=99999
        )
        p <- p +
            labs(title = toString(i))+
            theme_void() 
        print(p)
        p <- ggraph(gene_graph, x=UMAP_1, y=-UMAP_2)
        p <- p + geom_edge_diagonal(alpha = 0.05)
        #p <- p + geom_node_point(aes(fill=centrality, size=ratio), color='darkgrey', shape=21)
        p <- p + geom_node_point(aes(alpha=alpha_u, size=ratio), fill = 'darkorange', color='darkgrey', shape=21)
        #p <- p + scale_fill_gradient(low = '#FFD580', high = 'orange')
        #p <- p + scale_fill_viridis(option = "B", end = 0.85)
        p <- p + geom_node_text(
            aes(label=u_name),
            repel=T, max.overlaps=99999
        )
        p <- p +
            theme_void() 
        print(p)
        p <- ggraph(gene_graph, x=UMAP_1, y=-UMAP_2)
        p <- p + geom_edge_diagonal(alpha = 0.05)
        #p <- p + geom_node_point(aes(fill=centrality, size=ratio), color='darkgrey', shape=21)
        p <- p + geom_node_point(aes(alpha=alpha_l, size=ratio), fill = 'purple', color='darkgrey', shape=21)
        #p <- p + scale_fill_gradient(low = 'light purple', high = 'purple')
        #p <- p + scale_fill_viridis(option = "B", end = 0.85)
        p <- p + geom_node_text(
            aes(label=l_name),
            repel=T, max.overlaps=99999
        )
        p <- p +
            theme_void() 
        print(p)
        if( FALSE)
        {
        gene_graph <- as_tbl_graph(gene_net) %>%
            activate(edges) %>%
            mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
            activate(nodes) %>%
            mutate(centrality=centrality_pagerank()) %>%
            inner_join(umap_tbl, by=c('name'='gene'))

        print(gene_graph)
        }

    }
	dev.off()
}


if(FALSE)
{
    df <- readRDS(df, file = 'R_data/monocle3/gene_annotation')
    #print(df)

    corr <- read.table(file = 'R_data/regulon_cor_python.csv', sep = ',', header = TRUE, row.names = 1)
    #print(head(regulon_link))
    regulon_link <- read.table(file = 'output/Step2_regulonTargetsInfo.tsv', sep = '\t', header = TRUE)
    #regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% df[,'name']) | (regulon_link[,'gene'] %in% regulon_link[,'TF']),]

    regulon_link_tmp <- regulon_link[ (regulon_link[,'gene'] %in% df[,'name']) ,]
    regulon_link_TF_tmp <- regulon_link_tmp[ (regulon_link_tmp[,'gene'] %in% regulon_link_tmp[,'TF']) ,]

    regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% colnames(corr)),]
    regulon_link <- regulon_link[ (regulon_link[,'TF'] %in% colnames(corr)),]
    regulon_link <- regulon_link[ (regulon_link[,'TF'] %in% regulon_link_TF_tmp[,'TF']),]
    regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% regulon_link_TF_tmp[,'TF']),]
    regulon_link_TF <- regulon_link
    print(dim(regulon_link))

    library(Matrix, lib='~/.conda/pkgs/r-matrix-1.5_4-r42he1ae0d6_0/lib/R/library')
    library(igraph, lib="~/.conda/pkgs/r-igraph-1.5.1-r43hbd8eb98_0/lib/R/library/")
    library(irlba, lib="~/.conda/pkgs/r-irlba-2.3.5.1-r43h316c678_1/lib/R/library/")
    library(tidygraph)
    library(dplyr)
    library(tibble)
    library(tidyr)
    library(ggraph)
    library(uwot)
    library(ggrepel)
    library(ggplot2, lib = '~/.conda/pkgs/r-ggplot2-3.4.0-r42hc72bb7e_1/lib/R/library')
    
    #for( i in c(100,105,110,115,116,117,118,119,120,121,122,123,124,125,126))
    #for( i in c(2,5,10,20,25,30,35,40,45,50,60,70,75,80,100,125,150))
    #for( i in c(30,31))
#for( i in c(15,20,25,30,35,40,45,50))
#corr <- readRDS( file = 'R_data/monocle3/gene_corr_matrix')
    gene_net <- regulon_link_TF %>%
            select(TF, gene, everything()) %>%
            group_by(TF)
    print(gene_net)
    reg_mat <- regulon_link %>%
        as_tibble() %>%
        select(TF,gene,CoexWeight) %>%
        pivot_wider(names_from = TF, values_from = CoexWeight, values_fill = 0)%>%
        column_to_rownames('gene') %>% as.matrix()
    reg_mat[is.na(reg_mat)] <- 0
    #reg_mat <- corr[rownames(reg_mat),colnames(reg_mat)] * (reg_mat*0.1+1)
    reg_mat <- corr[rownames(reg_mat),colnames(reg_mat)]
	write.table(rownames(reg_mat), file = paste0('R_data/regulon_name.csv'), quote = F, sep = ',')
}        

if(TRUE)
{
    df <- readRDS(df, file = 'R_data/monocle3/gene_annotation')
    #print(df)

    corr <- read.table(file = 'R_data/regulon_cor_python.csv', sep = ',', header = TRUE, row.names = 1)
    #print(head(regulon_link))
    regulon_link <- read.table(file = 'output/Step2_regulonTargetsInfo.tsv', sep = '\t', header = TRUE)
    #regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% df[,'name']) | (regulon_link[,'gene'] %in% regulon_link[,'TF']),]

    regulon_link_tmp <- regulon_link[ (regulon_link[,'gene'] %in% df[,'name']) ,]
    regulon_link_TF_tmp <- regulon_link_tmp[ (regulon_link_tmp[,'gene'] %in% regulon_link_tmp[,'TF']) ,]

    regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% colnames(corr)),]
    regulon_link <- regulon_link[ (regulon_link[,'TF'] %in% colnames(corr)),]
    regulon_link <- regulon_link[ (regulon_link[,'TF'] %in% regulon_link_TF_tmp[,'TF']),]
    regulon_link <- regulon_link[ (regulon_link[,'gene'] %in% regulon_link_TF_tmp[,'TF']),]
    regulon_link_TF <- regulon_link
    print(dim(regulon_link))

    library(Matrix, lib='~/.conda/pkgs/r-matrix-1.5_4-r42he1ae0d6_0/lib/R/library')
    library(igraph, lib="~/.conda/pkgs/r-igraph-1.5.1-r43hbd8eb98_0/lib/R/library/")
    library(irlba, lib="~/.conda/pkgs/r-irlba-2.3.5.1-r43h316c678_1/lib/R/library/")
    library(tidygraph)
    library(dplyr)
    library(tibble)
    library(tidyr)
    library(ggraph)
    library(uwot)
    library(ggrepel)
    library(ggplot2, lib = '~/.conda/pkgs/r-ggplot2-3.4.0-r42hc72bb7e_1/lib/R/library')
    # pdf("figure/TF_umap_day2_3_final_opt_sel.pdf", width = 8, height = 7, onefile = TRUE)
    pdf("figure_jh/TF_umap_day2_3_final_opt_sel_color.pdf", width = 8, height = 7, onefile = TRUE)
    
    #for( i in c(100,105,110,115,116,117,118,119,120,121,122,123,124,125,126))
    #for( i in c(2,5,10,20,25,30,35,40,45,50,60,70,75,80,100,125,150))
    #for( i in c(30,31))
    #for( i in c(15,20,25,30,35,40,45,50))
    for( i in c(10))
    {
        #corr <- readRDS( file = 'R_data/monocle3/gene_corr_matrix')
        gene_net <- regulon_link_TF %>%
                select(TF, gene, everything()) %>%
                group_by(TF)
        print(gene_net)
        reg_mat <- regulon_link %>%
            as_tibble() %>%
            select(TF,gene,CoexWeight) %>%
            pivot_wider(names_from = TF, values_from = CoexWeight, values_fill = 0)%>%
            column_to_rownames('gene') %>% as.matrix()
        reg_mat[is.na(reg_mat)] <- 0
        #reg_mat <- corr[rownames(reg_mat),colnames(reg_mat)] * (reg_mat*0.1+1)
        reg_mat <- corr[rownames(reg_mat),colnames(reg_mat)]
        
        print(dim(reg_mat))

        
        if (TRUE){
            pca_mat <- prcomp(reg_mat)$x[,1:35]
            rownames(pca_mat) <- rownames(reg_mat)
            reg_mat <- as.matrix(pca_mat)
        }

        set.seed(123)
        #umap_tbl <- umap(reg_mat, n_neighbors = i)
        umap_tbl <- umap(reg_mat, min_dist = 0.01*20)
        umap_tbl <- as.data.frame(umap_tbl)
        colnames(umap_tbl) <- c('UMAP_1', 'UMAP_2')

        umap_tbl <- umap_tbl %>%
            as_tibble(rownames='TF')
        
        tf_ratio <- read.table(paste0('python_data/tf_ratio_anno_',toString(i),'_sel.csv'),sep = ',', header = TRUE, row.names = 1)
        print(head(tf_ratio))
        umap_tbl$name <- umap_tbl$TF
        rownames(umap_tbl) <- umap_tbl$TF
        umap_tbl <- umap_tbl[umap_tbl$TF %in% rownames(tf_ratio),]
        tf_ratio[,'TF'] <- rownames(tf_ratio)
        print('setting umap')
        tf_ratio <- tf_ratio[,c('TF','ratio','ratio_u','ratio_l','weight_time','alpha_u','alpha_l','u_name','l_name','all_name')]
        umap_tbl <- umap_tbl %>%
            select(TF, everything()) %>%
            group_by(TF) %>%
            left_join(tf_ratio, by=c('TF'='TF'))
        print('incorporating ratio info')
        if( TRUE)
        {
            gene_graph <- as_tbl_graph(gene_net) %>%
                activate(edges) %>%
                mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
                activate(nodes) %>%
                mutate(centrality=centrality_pagerank()) %>%
                inner_join(umap_tbl, by=c('name'='TF'))

            print(gene_graph)
        }
        p <- ggraph(gene_graph, x=UMAP_1, y=UMAP_2)
        p <- p + geom_edge_diagonal(alpha = 0.02)
        #p <- p + geom_node_point(aes(fill=centrality, size=ratio), color='darkgrey', shape=21)
        # p <- p + geom_node_point(aes(fill=weight_time, size=ratio), color='darkgrey', shape=21)
        p <- p + geom_node_point(aes(fill=weight_time, size=ratio), color='darkgrey', shape=21)
        # p <- p + scale_fill_viridis(option = "B", end = 0.85)
        # p <- p + scale_fill_viridis(option = "D", direction =-1, end = 0.85)
        p <- p + scale_colour_gradient(low="white", high="darkblue")
        p <- p + geom_node_text(
            aes(label=all_name),
            repel=T, max.overlaps=99999
        )
        p <- p +
            labs(title = toString(i))+
            theme_void() 
        print(p)
        p <- ggraph(gene_graph, x=UMAP_1, y=UMAP_2)
        p <- p + geom_edge_diagonal(alpha = 0.02)
        #p <- p + geom_node_point(aes(fill=centrality, size=ratio), color='darkgrey', shape=21)
        p <- p + geom_node_point(aes(alpha=alpha_u, size=ratio), fill = 'darkorange', shape=21)
        # p <- p + geom_node_point(aes(alpha=alpha_u, size=ratio), fill = 'green4', color='darkgrey', shape=21)
        #p <- p + scale_fill_gradient(low = '#FFD580', high = 'orange')
        #p <- p + scale_fill_viridis(option = "B", end = 0.85)
        p <- p + geom_node_text(
            aes(label=u_name),
            repel=T, max.overlaps=99999
        )
        p <- p +
            theme_void() 
        print(p)
        p <- ggraph(gene_graph, x=UMAP_1, y=UMAP_2)
        p <- p + geom_edge_diagonal(alpha = 0.02)
        #p <- p + geom_node_point(aes(fill=centrality, size=ratio), color='darkgrey', shape=21)
        # p <- p + geom_node_point(aes(alpha=alpha_l, size=ratio), fill = 'purple', color='darkgrey', shape=21)
        p <- p + geom_node_point(aes(alpha=alpha_l, size=ratio), fill = 'gray15', shape=21)
        #p <- p + scale_fill_gradient(low = 'light purple', high = 'purple')
        #p <- p + scale_fill_viridis(option = "B", end = 0.85)
        p <- p + geom_node_text(
            aes(label=l_name),
            repel=T, max.overlaps=99999
        )
        p <- p +
            theme_void() 
        print(p)
        if( FALSE)
        {
        gene_graph <- as_tbl_graph(gene_net) %>%
            activate(edges) %>%
            mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
            activate(nodes) %>%
            mutate(centrality=centrality_pagerank()) %>%
            inner_join(umap_tbl, by=c('name'='gene'))

        print(gene_graph)
        }

    }
	dev.off()
}
