if(FALSE)
{
	#BiocManager::install('DropletUtils', lib = '/tmp2/b04401068/R/R/4.0/')
	install.packages('Matrix', lib = '/tmp2/b04401068/lib/R/')
}

#seurat preprocessing, normalizing
if(FALSE)
{
	library(Seurat)
    
    #read file
	sc.data <- Read10X(data.dir = './R_data/aggr/')
	sc <- CreateSeuratObject( counts = sc.data, project = 'sc')
	print('Finish creating')
    
    #save intermediate file for possible future use
	saveRDS(sc, file = './R_data/Seurat/sc_create.rds')
	
    #quality control
    sc[['RNA']]@counts
	sc[['percent.mt']] <- PercentageFeatureSet(sc, pattern = '^MT-')
	sc <- subset(sc, subset = percent.mt < 10 & nCount_RNA > 5000)
	
    #normalize data
    sc <- NormalizeData(sc, normalization.method = 'LogNormalize', scale.factor = 10000)
	sc[['RNA']]@data
	
    #save normalize data
    saveRDS(sc, file = './R_data/Seurat/sc_normal.rds')
	print('Finish normalization')
}

#seurat 1-3d count
if(FALSE)
{
	library(Seurat)
	sc <- readRDS(file = 'R_data/Seurat/sc_create.rds')
	sc[['RNA']]@counts
	sc[['percent.mt']] <- PercentageFeatureSet(sc, pattern = '^MT-')
	sc <- subset(sc, subset = percent.mt < 10 & nCount_RNA > 5000)

	name = colnames(x = sc)
	name = grep("-2|-3|-4",x = name)
	print(head(name))
	sc <- subset(sc, cells = name)

	motif <- as.data.frame(sc@assays$RNA@counts)
	rownames(motif) <- rownames(x = sc)

	print(dim(motif))
	print(head(motif))
	saveRDS(motif, file = 'R_data/Seurat/sc_count_1-3d_with_doublets.rds')
}

#removing doublets
if( FALSE )
{
	exp_tmp <- readRDS(file = 'R_data/Seurat/sc_count_1-3d_with_doublets.rds')
	exp <- as.matrix(exp_tmp)
	rm(exp_tmp)
	cell <- read.csv( paste0('R_data/Seurat/Singlet.csv'), header = TRUE, row.names = 1, sep = ',')
	cell <- rownames(cell)
	exp <- exp[,cell]
	saveRDS(exp, file = 'R_data/Seurat/sc_count_1-3d.rds')
}



