#scenic initiation
if( FALSE )
{
	exp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	exp <- as.matrix(exp)
	print('finish loading')

	#dir.create('./int')
	library(SCENIC)
	org <- 'hgnc'
    #need to download the transcription factor database
	dbDir <- '/tmp2/b04401068/working_data/hg19'
	scenicOptions <- initializeScenic(org = org, dbDir = dbDir, nCores = 20)
	saveRDS(scenicOptions, file = 'int/scenicOptions.Rds')
    
    #use the scenic filtering for keeping genes 
	genesKept <- geneFiltering( exp, scenicOptions = scenicOptions)

    #checking some genes
	interestingGenes <- c('YAP1','NFE2L2','ACSL4','TFRC')
    print(interestingGenes[which(!interestingGenes %in% genesKept)])
}

#exporting geneskept
if(FALSE)
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.rds')
	genesKept <- loadInt( scenicOptions, 'genesKept')
	print(head(genesKept))
	write.table( genesKept,file = 'working_data/genesKept.csv',quote = F,row.names = T,col.names = T,sep = ',')
}

#output correlation matrix, to only keep positive correlated interaction
if( FALSE )
{
	exp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	print(dim(exp))
	print('finish loading')

	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	genesKept <- loadInt( scenicOptions, 'genesKept')
	exp <- exp[genesKept,]
	exp <- t(exp)
	exp <- as.data.frame(exp)
	print(dim(exp))
	library(feather)
	write_feather( exp, paste0('working_data/exp_1-3d_python.feather'))
}

#retrive spearman correlation from python input
if( FALSE )
{
	library(SCENIC)
	cor <- read.csv( paste0('working_data/exp_spear_cor_tmp.csv'), header = TRUE, row.names = 1, sep = ',')
	print(head(rownames(cor)))
	print(dim(cor))
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	saveRDS( cor, file = getIntName(scenicOptions, "corrMat"))
}

#fix gene name problem of the spearman correlation matrix
if(FALSE)
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	cor <- readRDS(file = getIntName(scenicOptions, "corrMat"))
	old <- c('IQCJ.SCHIP1','NKX2.5','HLA.F','HLA.A','HLA.E','HLA.C','HLA.B','NKX3.1','BDNF.AS')
	real <- c('IQCJ-SCHIP1','NKX2-5','HLA-F','HLA-A','HLA-E','HLA-C','HLA-B','NKX3-1','BDNF-AS')
	for( i in c(1:9))
	{
		names(cor)[names(cor) == old[i]] <- real[i]
	}
	cor <- as.matrix(cor)
	print(dim(cor))
	print(head(rownames(cor)))
	saveRDS( cor, file = getIntName(scenicOptions, "corrMat"))
}

#scenic runGenie3, this doesen't work, program will crash
if( FALSE )
{
	exp_tmp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	exp <- as.matrix(exp_tmp)
	rm(exp_tmp)
	print(head(colnames(exp)))
	print(head(rownames(exp)))
	print(dim(exp))
	print('finish loading')

	#dir.create('./int')
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	#tmp <- getIntName(scenicOptions, 'genie3wm')
	#print(tmp)
	genesKept <- loadInt( scenicOptions, 'genesKept')
	exp_filter_tmp <- exp[genesKept,]
	exp_filter <- log2(exp_filter_tmp+1)
	rm(exp)
	rm(exp_filter_tmp)
	runGenie3(exp_filter, scenicOptions, nParts = 10)
}

#runGenie3 check parrellel on 10, create data
if( FALSE )
{
	dir.create('./working_data/para/')
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	regulator <- getDbTfs(scenicOptions)
	genekept <- readRDS(file = './int/1.1_genesKept.Rds')
	exp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	exp <- exp[genekept,]
	exp <- log2(exp+1)
	regulator <- regulator[regulator %in% rownames(exp)]
	print(length(regulator))
	exp <- t(exp)
	exp <- as.data.frame(exp)
	print(dim(exp))
	write.table( regulator, file = paste0('working_data/para/regulator.csv'), quote = F,row.names = T,col.names = T,sep = '\t')
	library(feather)
	increment <- 950
	for( i in c(1:10) )
	{
		start <- (i-1)*increment + 1
		end <- start + increment - 1
		if( end > length(genekept) )
			end <- length(genekept)
		target <- genekept[c(start:end)]
		print(i)
		print(length(target))
		output <- union( regulator, target)
		write.table( target, file = paste0('working_data/para/target_', toString(i),'.csv'), quote = F,row.names = T,col.names = T,sep = '\t')
		write_feather( exp[,output], paste0('working_data/para/exp_', toString(i),'.feather'))
	}
}

#runGenie3 check parrellel on 10, running
if( FALSE )
{
	library(GENIE3)
	library(feather)
	i <- toString(2)
	regulator <- read.csv( paste0('working_data/para/regulator.csv'), header = TRUE, sep = '\t')
	target <- read.csv( paste0('working_data/para/target_', toString(i),'.csv'), header = TRUE, sep = '\t')
	exp <- read_feather( paste0('working_data/para/exp_', toString(i),'.feather'))
	regulator <- as.vector(regulator[,1])
	target <- as.vector(target[,1])
	exp <- as.matrix(exp)
	exp <- t(exp)
	print(head(rownames(exp)))
	print(head(target))
	print(head(regulator))
	set.seed(123)
	print('finish reading')
	weightMat <- GENIE3(exp, regulators = regulator, targets = target, verbose = TRUE, nCores = 20)
	write.table( weightMat,file = paste0('genie3_1-3d_weightmatrix_', i,'.csv'),quote = F,row.names = T,col.names = T,sep = ',')
}

#retrive parallel
if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	for( i in c(1,2,3,4,5,6,7,8,9,10))
	{
		mat <- read.csv( paste0('working_data/genie3_1-3d_weightmatrix_',toString(i),'.csv'), header = TRUE, row.names = 1, sep = ',')
		print(head(rownames(mat)))
		print(dim(mat))
		saveRDS( mat, file = paste0("1.3_GENIE3_weightMatrix_part_",toString(i),'.Rds'))
	}
}

#fit the combined data into scenic pipeline
if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	genesKept <- loadInt( scenicOptions, 'genesKept')
	final <- list()
	if(TRUE)
	{
		final <- readRDS( file = paste0("int/1.3_GENIE3_weightMatrix_part_",toString(1),'.Rds'))
		genesDone <- colnames(final)
		if( FALSE)
		{
            for( i in c(2,3,4,5,6,7,8,9,10))
            {
                mat <- readRDS( file = paste0("int/1.3_GENIE3_weightMatrix_part_",toString(i),'.Rds'))
                final <- cbind(final, mat)
            }
		}
	}
	print(length(genesDone))
	print(length(genesKept))
	genes <- setdiff( genesKept, genesDone)
	saveRDS(final, file = '1.3_GENIE3_weightMatrix_part_1.Rds')
	print(genes)
}

#copying code for runGenie3
if(FALSE)
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	print(getSettings(scenicOptions, 'modules/weightThreshold'))
	weight <- readRDS(file = 'int/1.3_GENIE3_weightMatrix_part_1.Rds')
	print(head(colnames(weight)))
	print(head(rownames(weight)))
	print(dim(weight))
	weight <- as.matrix(weight)
	linklist <- GENIE3::getLinkList(weight, threshold=getSettings(scenicOptions, 'modules/weightThreshold'))
	print('finish link list')
	rm(weight)
	colnames(linklist) <- c('TF','Target','weight')
	linklist <- linklist[order(linklist[,'weight'],decreasing=TRUE),]
	linklist <- unique(linklist)
	rownames(linklist) <- NULL
	saveRDS(linklist, file=getIntName(scenicOptions,'genie3ll'))
}

#run coexnetwork
if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	scenicOptions@settings$verbose <- TRUE
	scenicOptions@settings$seed <- 123
	print('starting')
	start_time <- Sys.time()
	scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
	end_time <- Sys.time()
	print('finishing')
	print(start_time - end_time)
	saveRDS(scenicOptions, file = 'int/scenicOptions.Rds')
}

#run regulon
if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	scenicOptions@settings$verbose <- TRUE
	scenicOptions@settings$seed <- 123
	scenicOptions@settings$nCores <- 20
	print('starting')
	start_time <- Sys.time()
	scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
	end_time <- Sys.time()
	print('finishing')
	print(start_time - end_time)
	saveRDS(scenicOptions, file = 'int/scenicOptions.Rds')
}

#run cell
if( FALSE )
{
	exp_tmp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	exp <- as.matrix(exp_tmp)
	rm(exp_tmp)
	print(dim(exp))
	print('finish loading')

	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	genesKept <- loadInt( scenicOptions, 'genesKept')
	exp_filter_tmp <- exp[genesKept,]
	exp_filter <- log2(exp_filter_tmp+1)
	rm(exp)
	rm(exp_filter_tmp)

	scenicOptions@settings$verbose <- TRUE
	scenicOptions@settings$seed <- 123
	scenicOptions@settings$nCores <- 5
	print('starting')
	start_time <- Sys.time()
	scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exp_filter)
	end_time <- Sys.time()
	print('finishing')
	print(start_time - end_time)
	saveRDS(scenicOptions, file = 'int/scenicOptions.Rds')
}

#run binary, this doesn't work, the shiny app is not available for remote access
if( FALSE )
{
	exp_tmp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	exp <- as.matrix(exp_tmp)
	rm(exp_tmp)
	print(dim(exp))
	print('finish loading')

	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	genesKept <- loadInt( scenicOptions, 'genesKept')
	exp_filter_tmp <- exp[genesKept,]
	exp_filter <- log2(exp_filter_tmp+1)
	rm(exp)
	rm(exp_filter_tmp)
	
	scenicOptions@settings$verbose <- TRUE
	scenicOptions@settings$seed <- 123
	scenicOptions@settings$nCores <- 5

	aucellApp <- plotTsne_AUCellApp(scenicOptions, exp_filter)
	savedSelections <- shiny::runApp(aucellApp, port = 4694)

	newthres <- savedSelections$thresholds
	scenicOptions@fileNames$int['aucell_thresholds',1] <- 'int/newThresholds.Rds'
	#saveRDS(newthres, file=getIntName(scenicOptions, 'aucell_thresholds'))
	#saveRDS(scenicOptions, file = 'int/scenicOptions.Rds')
}

#run binary
if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	thres = readRDS(file=getIntName(scenicOptions, 'aucell_thresholds'))
	print(thres)
}

#output for python
if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	#regulon <- loadInt(scenicOptions, 'aucell_regulonAUC')
	regulon <- readRDS('int/4.1_binaryRegulonActivity.Rds')
	print(dim(regulon))
	print(head(rownames(regulon)))
	print(head(colnames(regulon)))
	write.table( rownames(regulon), file = paste0('working_data/regulon_auc_bin_row.csv'), quote = F,row.names = F,col.names = F,sep = '\t')
	print('write name')
	library(feather)
	#regulon <- as.data.frame(regulon@assays@data@listData$AUC)
	regulon <- as.data.frame(regulon)
	print(dim(regulon))
	write_feather( regulon, paste0('working_data/regulon_auc_bin.feather'))
}

#plot heatmap optimizing
if( FALSE )
{
	regulon <- readRDS('int/4.1_binaryRegulonActivity.Rds')
	
	anno <- read.table('R_data/Regulon/rna_velocity_latent_time.csv', sep = ',', check.names = FALSE, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
	anno <- anno[,colnames(regulon)]
	anno <- t(anno)
	anno <- as.data.frame(anno)
	#regulon <- t(t(regulon) * anno[,'day'])
	cindex <- 4-anno$day
	anno$day <- factor(anno$day)
	anno[,'latent'] <- as.numeric(anno[,'latent'])
	print(dim(anno))
	#regulon[regulon == 0] <- -1
	#regulon <- t(t(regulon) * anno[,'latent'])
	#regulon[regulon < 0] <- -1
	#print('finish transformation')
	
	
	#ann_colors <- list(day = c('darkgreen','navyblue','oldlace'), latent = colorRampPalette(c('navy','white','firebrick3'))(50))
	ann_colors <- list( day = c('darkgoldenrod1','wheat','tomato4'), latent = colorRampPalette(c('navy','white','firebrick3'))(50),cluster = c('#F8766D','#7CAE00','#00BFC4','#C77CFF'))
	#ann_colors <- list( day = c('darkgoldenrod1','seashell','tomato4'), latent = colorRampPalette(c('navy','white','firebrick3'))(50),cluster = c('#F8766D','#7CAE00','#00BFC4','#C77CFF'))
	#ann_colors <- list( day = c('#1B9E77','#7570B3','#D95F02'), latent = colorRampPalette(c('navy','white','firebrick3'))(50),cluster = c('#F8766D','#7CAE00','#00BFC4','#C77CFF'))
	#ann_colors <- list( day = c('darkgoldenrod1','lightcyan1','tomato4'), latent = colorRampPalette(c('navy','white','firebrick3'))(50))
	print('color annotation')
	
	#change threshold for some regulons
	if( TRUE )
	{
        temp <- readRDS('int/3.4_regulonAUC.Rds')
		temp <- as.data.frame(temp@assays@data@listData$AUC)
		temp <- temp[,colnames(regulon)]
		regulon['TEAD1 (106g)',] <- ifelse(temp['TEAD1 (106g)',]>0.29, 1, 0)
		#regulon['TEAD1_extended (107g)',] <- ifelse(temp['TEAD1_extended (107g)',]>0.3, 1, 0)
		#regulon['ZEB1_extended (401g)',] <- ifelse(temp['ZEB1_extended (401g)',]>0.1, 1, 0)
		#regulon['HIF1A_extended (166g)',] <- ifelse(temp['HIF1A_extended (166g)',]>0.21, 1, 0)
		print('adj thres')
	}

	chrc <- hclust(dist(t(regulon)))
	#print('finish clustering')

	treecuts <- cutree(chrc, k = 4)
	treecuts<-data.frame(treecuts, row.names = colnames(regulon))
	anno$cluster <- factor(treecuts$treecuts)

	dendro <- as.dendrogram(chrc)
	dendro <- reorder(dendro, cindex)
	#regulon[regulon==0] <- NA
	#regulon <- regulon[rhcr$order, chrc$order]
	#anno <- anno[chrc$order,]
	print(dim(regulon))

	rhcr <- hclust(dist(regulon))
	dendro_r <- as.dendrogram(rhcr)
	dendro_r <- reorder(dendro_r, rowMeans(regulon))
	
    #pdf('cluster_output/R_regulon_bin_all_time_wide_color_0.29_adj.pdf', width = 15, height = 7, onefile = TRUE)
	#NMF::aheatmap(regulon, scale = 'none', revC = TRUE, Rowv = T, Colv = T, color = c('white','black'), filename = 'cluster_output/R_regulon_bin_all_adj.pdf')
	#NMF::aheatmap(regulon, scale = 'none', annCol = anno, revC = TRUE, Colv = T, Rowv = T, color = c('white','black'), filename = 'cluster_output/R_regulon_bin_all_time.pdf')
	#NMF::aheatmap(regulon, scale = 'none', annCol = anno, annColors = ann_colors, revC = T, Colv = cindex, Rowv = T, color = c('white','black'), filename = 'R_closer_optimize.png')
	#NMF::aheatmap(regulon, scale = 'none', annCol = anno, annColors = ann_colors, revC = T, Colv = cindex, Rowv = T, color = c('white','black'), filename = 'R_closer_optimize.png')
	reorderfun = function(d,w){d}
	#NMF::aheatmap(regulon, scale = 'none', annCol = anno, annColors = ann_colors, revC = T, Colv = dendro, reorderfun = reorderfun, Rowv = dendro_r, color = c('white','black'), filename = 'R_closer_optimize.png')
	#NMF::aheatmap(regulon, scale = 'none', annCol = anno, annColors = ann_colors, revC = T, Colv = dendro, reorderfun = reorderfun, Rowv = dendro_r, color = c('white','black'), filename = 'figure/Regulon_heatmap_final.pdf')
	NMF::aheatmap(regulon, scale = 'none', annCol = anno, annColors = ann_colors, revC = T, Colv = dendro, reorderfun = reorderfun, Rowv = dendro_r, color = c('white','black'), filename = 'figure/Regulon_heatmap_final_all_regulons.pdf')
	#NMF::aheatmap(regulon, scale = 'none', annCol = anno, annColors = ann_colors, revC = T, Colv = dendro, reorderfun = reorderfun, Rowv = T, color = c('white','black'), filename = 'R_closer_optimize.png')
	#library(NMF)
	#NMF::aheatmap(regulon, scale = 'none', annCol = anno, annColors = ann_colors, revC = TRUE, Colv = T, Rowv = T, color = colorRampPalette(c('white','navy','grey80','firebrick3'))(4))
	#NMF::aheatmap(regulon, scale = 'none', annCol = anno, annColors = ann_colors, revC = TRUE, Colv = NA, Rowv = NA, color = colorRampPalette(c('white','navy','grey80','firebrick3'))(4))
	#breaks <- seq( -1,1,length.out = 101)
	#NMF::aheatmap(regulon, scale = 'none', breaks = breaks, annCol = anno, annColors = ann_colors, revC = TRUE, Colv = T, Rowv = T, color = colorRampPalette(c('white','yellow','navy','grey80','firebrick3'))(200))
	#NMF::aheatmap(regulon, scale = 'none', breaks = breaks, annCol = anno, annColors = ann_colors, Colv = NA, Rowv = NA, color = colorRampPalette(c('white','navy','grey80','firebrick3'))(200))
	#NMF::aheatmap(regulon, scale = 'none', annCol = anno, annColors = ann_colors, Colv = NA, Rowv = NA, color = colorRampPalette(c('white','navy','grey80','firebrick3'))(200))
	#dev.off()
}

getTF <- function(regulonName, sep="\\s")
{
  gsub(paste0(sep,"\\(\\d+g)"), "",regulonName)
}

onlyNonDuplicatedExtended <- function(regulonNames)
{
  regulonNames <- unname(regulonNames)
  tfs <- getTF(regulonNames)
  tfs <- gsub("_extended", "", tfs)
  splitRegulons <- split(regulonNames, tfs)[unique(tfs)]
  
  ret <- sapply(splitRegulons, function(x) 
  {
    split(x, grepl("_extended", x))[[1]] # False (direct) will be first
  })
  
  return(ret)
}

if( FALSE )
{
    regulon <- readRDS('int/3.4_regulonAUC.Rds')
	regulon <- regulon[which(rownames(regulon)%in% onlyNonDuplicatedExtended(rownames(regulon))),]
    print(head(regulon))
}

if( FALSE )
{
	regulon <- readRDS('int/4.1_binaryRegulonActivity.Rds')
	
	anno <- read.table('R_data/Regulon/rna_velocity_latent_time.csv', sep = ',', check.names = FALSE, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
	anno <- anno[,colnames(regulon)]
	anno <- t(anno)
	anno <- as.data.frame(anno)
	cindex <- 4-anno$day
	anno$day <- factor(anno$day)
	anno[,'latent'] <- as.numeric(anno[,'latent'])
	print(dim(anno))
	
	#change threshold for some regulons
	if( TRUE )
	{
        temp <- readRDS('int/3.4_regulonAUC.Rds')
		temp <- as.data.frame(temp@assays@data@listData$AUC)
		temp <- temp[,colnames(regulon)]
		regulon['TEAD1 (106g)',] <- ifelse(temp['TEAD1 (106g)',]>0.29, 1, 0)
		print('adj thres')
	}

	chrc <- hclust(dist(t(regulon)))
	dendro <- as.dendrogram(chrc)
	dendro <- reorder(dendro, cindex)
	reorderfun = function(d,w){d}
	print('finish clustering')
	
	regulon <- regulon[which(rownames(regulon)%in% onlyNonDuplicatedExtended(rownames(regulon))),]
	regMinCells <- names(which(rowSums(regulon)> ncol(regulon)*0.01))
	reguCor <- cor( t(regulon[regMinCells,]) )
	reguCor[which(is.na(reguCor))] <- 0
	diag(reguCor) <- 0
	corrRegs <- names(which(rowSums(abs(reguCor)>0.3)>0))
	regulon <- regulon[corrRegs,]

	rhcr <- hclust(dist(regulon))
	dendro_r <- as.dendrogram(rhcr)
	dendro_r <- reorder(dendro_r, -rowMeans(regulon))
	
	treecuts <- cutree(chrc, k = 4)
	treecuts<-data.frame(treecuts, row.names = colnames(regulon))
	print(dim(regulon))

	anno$cluster <- factor(treecuts$treecuts)
	#ann_colors <- list( day = c('#1B9E77','#7570B3','#D95F02'), latent = colorRampPalette(c('navy','white','firebrick3'))(50),cluster = c('#F8766D','#7CAE00','#00BFC4','#C77CFF'))
	ann_colors <- list( day = c('darkgoldenrod1','wheat','tomato4'), latent = colorRampPalette(c('navy','white','firebrick3'))(50),cluster = c('#F8766D','#7CAE00','#00BFC4','#C77CFF'))
	
	NMF::aheatmap(regulon, scale = 'none', annCol = anno, annColors = ann_colors, Colv = dendro, reorderfun = reorderfun, Rowv = dendro_r, color = c('white','black'), filename = 'figure/Regulon_heatmap_final_fewer_regulons.pdf')
}

#getting regulon cluster information, setting correlation threshold at 0.3 to filter relevant regulon
if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	regulon <- readRDS('int/4.1_binaryRegulonActivity.Rds')
	
	#change threshold for some regulons
	if( TRUE )
	{
		temp <- loadInt(scenicOptions, 'aucell_regulonAUC')
		temp <- as.data.frame(temp@assays@data@listData$AUC)
		temp <- temp[,colnames(regulon)]
		regulon['TEAD1 (106g)',] <- ifelse(temp['TEAD1 (106g)',]>0.29, 1, 0)
		print('adj thres')
	}

	regulon <- regulon[which(rownames(regulon)%in% onlyNonDuplicatedExtended(rownames(regulon))),]
	regMinCells <- names(which(rowSums(regulon)> ncol(regulon)*0.01))
	reguCor <- cor( t(regulon[regMinCells,]) )
	reguCor[which(is.na(reguCor))] <- 0
	diag(reguCor) <- 0
	corrRegs <- names(which(rowSums(abs(reguCor)>0.3)>0))
	regulon <- regulon[corrRegs,]
	rhcr <- hclust(dist(regulon))
	treecuts <- cutree(rhcr, k = 7)
	treecuts<-data.frame(treecuts, row.names = rownames(regulon))
	write.table( treecuts, file = paste0('working_data/regulon_sel_cluster.csv'), quote = F,row.names = T,col.names = T,sep = ',')
}

#extracting cell for UMAP presentation
if( FALSE )
{
	regulon <- readRDS('int/4.1_binaryRegulonActivity.Rds')
	set.seed(123456)

	if( TRUE )
	{
        temp <- readRDS('int/3.4_regulonAUC.Rds')
		temp <- as.data.frame(temp@assays@data@listData$AUC)
		temp <- temp[,colnames(regulon)]
		regulon['TEAD1 (106g)',] <- ifelse(temp['TEAD1 (106g)',]>0.29, 1, 0)
		print('adj thres')
	}

	chrc <- hclust(dist(t(regulon)))
	print('finish clustering')

	treecuts <- cutree(chrc, k = 4)
	treecuts<-data.frame(treecuts, row.names = colnames(regulon))
	print(head(treecuts))
	write.table( treecuts, file = paste0('R_data/Regulon/cellcluster.csv'), quote = F,sep = '\t')

}

#export regulon
if( FALSE )
{
    temp <- readRDS('int/3.4_regulonAUC.Rds')
    temp <- as.data.frame(temp@assays@data@listData$AUC)
	saveRDS(temp, file = 'R_data/Regulon/regulon_auc.rds')
}



#monocle3 code
#adding day imformation
if(FALSE)
{
    library(monocle3)
    library(limma)
    cds <- readRDS('working_data/monocle3/cds_500_reduce')
    day <- strsplit2(colnames(cds),split = '-')[,2]
    day <- data.frame(day)
    rownames(day) <- colnames(cds)
    colnames(day) <- c('day')
    print(head(day))
    saveRDS(day, 'working_data/cell_meta_data')
}

#combining day information
if(FALSE)
{
	exp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	gene <- read.csv( paste0('working_data/genesKept_trim.csv'), header = TRUE, row.names = 1, sep = ',')
    print(dim(exp))
    gene <- gene[,'x']
    keep <- intersect(gene, rownames(exp))
    exp <- exp[keep,]
	#exp <- log2(exp+1)
    library(monocle3)
    cell_meta <- readRDS('working_data/cell_meta_data')
    cds <- new_cell_data_set(exp, cell_metadata = cell_meta)
    print(cds)
	gene <- read.csv( paste0('working_data/sc_var_1000.csv'), header = TRUE, row.names = 1, sep = ',')
    #gene <- gene[,'x']
	#gene <- read.csv( paste0('data/ferroptosis.txt'), header = TRUE, sep = '\t')
    #gene <- gene[,'gene']
    keep <- intersect(gene, rownames(exp))
    cds <- preprocess_cds(cds, use_genes = keep)
    #cds <- preprocess_cds(cds)
	set.seed(123)
    cds <- reduce_dimension(cds, umap.n_neighbors = 35L)
    #cds <- reduce_dimension(cds)
	pdf("umap_all_2.pdf", onefile = TRUE)
    p <- plot_cells(cds, color_cells_by = 'day')
    print(p)
    dev.off()
    saveRDS(cds, file = 'working_data/monocle3/cds_1000_reduce_2')
}


if(FALSE)
{
    library(monocle3)
	cds <- readRDS('working_data/monocle3/cds_1000_reduce')
    pdf("umap_1000_day.pdf", onefile = TRUE)
    p <- plot_cells(cds, color_cells_by = 'day')
    print(p)
    dev.off()
}

#tune clustering
if(FALSE)
{
    library(monocle3)
	cds <- readRDS('working_data/monocle3/cds_1000_reduce')
    #cds <- cluster_cells(cds, random_seed = 10000, k = 50)
    #cds <- cluster_cells(cds, resolution = 2)
	set.seed(1)
    cds <- cluster_cells(cds)
    #cds <- learn_graph(cds)
    cds <- learn_graph(cds, use_partition = FALSE)
	#saveRDS(cds, file = 'working_data/monocle3/cds_1000_unorder')
	saveRDS(cds, file = 'working_data/monocle3/cds_1000_unorder_no_partition')
	#pdf("umap_all_prin_2_3.pdf", onefile = TRUE)
    p <- plot_cells(cds, label_groups_by_cluster = TRUE, label_principal_points = TRUE)
    print(p)
    dev.off()
}

if(FALSE)
{
    library(monocle3)
	cds <- readRDS('working_data/monocle3/cds_all_unorder')
    cds <- order_cells(cds, root_pr_nodes = c('Y_101'))
    pdf("umap_all_time.pdf", onefile = TRUE)
    p <- plot_cells(cds, color_cells_by = 'pseudotime')
    print(p)
    dev.off()
}
#choosing subsegments
if(FALSE)
{
    library(monocle3)
	cds <- readRDS('working_data/monocle3/cds_1000_unorder_no_partition')
    cds <- order_cells(cds, root_pr_nodes = c('Y_76'))
    colData(cds)$closest_vertex <- cds@principal_graph_aux[['UMAP']]$pr_graph_cell_proj_closest_vertex[,1]
    pdf("umap_1000_pick_vertex.pdf", onefile = TRUE)
    p <- plot_cells(cds, color_cells_by = 'closest_vertex', label_cell_groups = FALSE)
    print(p)
    dev.off()
}
if(FALSE)
{
    library(monocle3)
	cds <- readRDS('working_data/monocle3/cds_1000_unorder_no_partition')
    cds <- order_cells(cds, root_pr_nodes = c('Y_76'))
    pdf("umap_1000_sub.pdf", onefile = TRUE)
    p <- plot_cells(cds, label_groups_by_cluster = TRUE, label_principal_points = TRUE)
    print(p)
    for( i in 200:210)
    {
        cds_sub <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_76'), ending_pr_nodes = c(paste0('Y_',toString(i))))
        p <- plot_cells(cds_sub)
        print(p)
    }
    dev.off()
}

#lower trajectory
if(FALSE)
{
    library(monocle3)
	cds <- readRDS('working_data/monocle3/cds_1000_unorder_no_partition')
    cds <- order_cells(cds, root_pr_nodes = c('Y_76'))
    cds_sub <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_76'), ending_pr_nodes = c('Y_206','Y_12','Y_4'))
    cds_sub2 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_206'), ending_pr_nodes = c(paste0('Y_189')))
    cds_sub3 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_206'), ending_pr_nodes = c(paste0('Y_201')))
    cds_sub4 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_210'), ending_pr_nodes = c(paste0('Y_189')))
    cds_sub5 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_189'), ending_pr_nodes = c(paste0('Y_201')))
    final <- union(colnames(cds_sub),colnames(cds_sub2))
    final <- union(final,colnames(cds_sub3))
    final <- union(final,colnames(cds_sub4))
    final <- union(final,colnames(cds_sub5))
    #final <- union(final,colnames(cds_sub3))
    print(head(final))
    cds_sub <- cds[,final]
    print(cds_sub)
    print('subset')
    #final <- c(colnames(cds_sub),colnames(cds_sub2),colnames(cds_sub3),colnames(cds_sub4))
    pdf("umap_1000_sub_2.pdf", onefile = TRUE)
    p <- plot_cells(cds_sub)
    print(p)
    dev.off()
}

#higher trajectory
if(FALSE)
{
    library(monocle3)
	cds <- readRDS('working_data/monocle3/cds_1000_unorder_no_partition')
    cds <- order_cells(cds, root_pr_nodes = c('Y_76'))
    cds_sub <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_76'), ending_pr_nodes = c('Y_208','Y_12','Y_4'))
    cds_sub2 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_208'), ending_pr_nodes = c(paste0('Y_202')))
    final <- union(colnames(cds_sub),colnames(cds_sub2))
    #final <- union(final,colnames(cds_sub3))
    print(head(final))
    cds_sub <- cds[,final]
    print(cds_sub)
    print('subset')
    #final <- c(colnames(cds_sub),colnames(cds_sub2),colnames(cds_sub3),colnames(cds_sub4))
    pdf("umap_1000_sub_3.pdf", onefile = TRUE)
    p <- plot_cells(cds_sub)
    print(p)
    dev.off()
}

#ploting pseudotime
if(FALSE)
{
    library(monocle3)
	cds <- readRDS('working_data/monocle3/cds_1000_unorder_k50')
    cds <- order_cells(cds, root_pr_nodes = c('Y_30','Y_107'))
    #cds <- order_cells(cds, root_pr_nodes = c('Y_76'))
	pdf("umap_1000_order_k50.pdf", onefile = TRUE)
    p <- plot_cells(cds, color_cells_by = 'pseudotime')
    print(p)
    dev.off()
}

#finding differential genes
if(FALSE)
{
    library(monocle3)
	cds <- readRDS('working_data/monocle3/cds_1000_unorder_no_partition')
    cds <- order_cells(cds, root_pr_nodes = c('Y_76'))
    pr_test_res <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 16)
    saveRDS(pr_test_res, file = 'working_data/monocle3/pr_test_1000')
    print('test')
    pr_deg_ids <- row.names(subset(pr_test_res, q_value < 0.05))
    print(head(pr_deg_ids))
	pdf("umap_1000_order_no_partition.pdf", onefile = TRUE)
    p <- plot_cells(cds, color_cells_by = 'pseudotime')
    print(p)
    dev.off()
}

#finding differential genes tead (actually, this part is nrf2...)
if(FALSE)
{
    library(monocle3)
	cds <- readRDS('working_data/monocle3/cds_1000_unorder_no_partition')
    cds <- order_cells(cds, root_pr_nodes = c('Y_76'))
    cds_sub <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_76'), ending_pr_nodes = c('Y_208','Y_12','Y_4'))
    cds_sub2 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_208'), ending_pr_nodes = c(paste0('Y_202')))
    final <- union(colnames(cds_sub),colnames(cds_sub2))
    #final <- union(final,colnames(cds_sub3))
    print(head(final))
    cds_sub <- cds[,final]

    pr_test_res <- graph_test(cds_sub, neighbor_graph = 'principal_graph', cores = 16)
    saveRDS(pr_test_res, file = 'working_data/monocle3/pr_test_1000_sub_tead')
    print('test')
	pdf("umap_1000_order_no_partition_tead.pdf", onefile = TRUE)
    p <- plot_cells(cds_sub, color_cells_by = 'pseudotime')
    print(p)
    dev.off()
}

#finding differential genes nrf2 (actually, this part is tead1)
if(FALSE)
{
    library(monocle3)
	cds <- readRDS('working_data/monocle3/cds_1000_unorder_no_partition')
    cds <- order_cells(cds, root_pr_nodes = c('Y_76'))
    cds_sub <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_76'), ending_pr_nodes = c('Y_206','Y_12','Y_4'))
    cds_sub2 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_206'), ending_pr_nodes = c(paste0('Y_189')))
    cds_sub3 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_206'), ending_pr_nodes = c(paste0('Y_201')))
    cds_sub4 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_210'), ending_pr_nodes = c(paste0('Y_189')))
    cds_sub5 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_189'), ending_pr_nodes = c(paste0('Y_201')))
    final <- union(colnames(cds_sub),colnames(cds_sub2))
    final <- union(final,colnames(cds_sub3))
    final <- union(final,colnames(cds_sub4))
    final <- union(final,colnames(cds_sub5))
    #final <- union(final,colnames(cds_sub3))
    print(head(final))
    cds_sub <- cds[,final]

    pr_test_res <- graph_test(cds_sub, neighbor_graph = 'principal_graph', cores = 16)
    saveRDS(pr_test_res, file = 'working_data/monocle3/pr_test_1000_sub_nrf2')
    print('test')
	pdf("umap_1000_order_no_partition_nrf2.pdf", onefile = TRUE)
    p <- plot_cells(cds_sub, color_cells_by = 'pseudotime')
    print(p)
    dev.off()
}


#adding NRF2 and TEAD1 into graph
if(FALSE)
{
    library(monocle3)
	cds <- readRDS('working_data/monocle3/cds_1000_unorder_no_partition')
    cds <- order_cells(cds, root_pr_nodes = c('Y_76'))
	regulon <- readRDS('working_data/regulon_auc.rds')
	#regulon <- readRDS('int/4.1_binaryRegulonActivity.Rds')
    #print(dim(regulon))
    #print(head(regulon))
    #regulon <- as.numeric(unlist(regulon[,colnames(cds)]))
    regulon <- regulon[,colnames(cds)]
    #print(head(regulon['NFE2L2_extended (53g)',]))
    #colData(cds)$nrf2 <- regulon['NFE2L2_extended (53g)',]
    colData(cds)$NRF2 <- as.numeric(unlist(regulon['NFE2L2_extended (53g)',]))
    #colData(cds)$TEAD1 <- as.numeric(unlist(regulon['TEAD1 (106g)',]))
    #colData(cds)$ZEB1 <- as.numeric(unlist(regulon['ZEB1_extended (401g)',]))
    #print(head(colData(cds)$nrf2))
    #print(typeof(colData(cds)$nrf2))
    #colData(cds)$nrf2 <- regulon['TEAD1 (106g)',]
    #colData(cds)$nrf2 <- regulon['XBP1 (377g)',]
    pdf("umap_1000_pick_nrf2.pdf", onefile = TRUE)
    p <- plot_cells(cds, color_cells_by = 'NRF2', show_trajectory_graph = FALSE, label_cell_groups = FALSE)
    print(p)
    dev.off()
}


#gemSmoooth for two trajectory
if(FALSE)
{
    library(monocle3)
	cds <- readRDS('working_data/monocle3/cds_1000_unorder_no_partition')
    cds_sub <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_76'), ending_pr_nodes = c('Y_206','Y_12','Y_4'))
    cds_sub2 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_206'), ending_pr_nodes = c(paste0('Y_189')))
    cds_sub3 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_206'), ending_pr_nodes = c(paste0('Y_201')))
    cds_sub4 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_210'), ending_pr_nodes = c(paste0('Y_189')))
    cds_sub5 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_189'), ending_pr_nodes = c(paste0('Y_201')))
    final <- union(colnames(cds_sub),colnames(cds_sub2))
    final <- union(final,colnames(cds_sub3))
    final <- union(final,colnames(cds_sub4))
    final <- union(final,colnames(cds_sub5))
    #final <- union(final,colnames(cds_sub3))
    #print(head(final))
    cds_sub <- cds[,final]
    cds_nrf2 <- order_cells(cds_sub, root_pr_nodes = c('Y_76'))
    
    gene <- readRDS(file = 'working_data/monocle3/pr_test_1000_sub_nrf2')
    #gene_nrf2 <- row.names(subset(gene, q_value < 10^(-100)))
    gene_nrf2 <- row.names(subset(gene, q_value < 10^(-35)))
    gene <- readRDS(file = 'working_data/monocle3/pr_test_1000_sub_tead')
    gene_tead <- row.names(subset(gene, q_value < 10^(-35)))

    gene <- union(gene_nrf2,gene_tead)
    cds_nrf2 <- cds_nrf2[gene,]
    print(cds_nrf2)
    print('finish data preparation')

    trend_formula = "~ splines::ns(pseudotime, df=3)"
    model_tbl <- monocle3:::fit_models(cds_nrf2, cores = 1, model_formula_str = trend_formula)
    colData(cds_nrf2)$pseudotime <- monocle3::pseudotime(cds_nrf2)
    newdata <- data.frame(pseudotime = seq(0, max(colData(cds_nrf2)$pseudotime), length.out=100))
    model_expectation <- monocle3:::model_predictions(model_tbl, new_data = newdata)
    saveRDS(model_expectation, file = 'working_data/monocle3/nrf2_model_2')
    print('nrf2 finish')
}
#adding NRF2 and TEAD1 into graph
if(FALSE)
{
    library(monocle3)
	cds <- readRDS('working_data/monocle3/cds_1000_unorder_no_partition')
    cds <- order_cells(cds, root_pr_nodes = c('Y_76'))
    cds_sub <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_76'), ending_pr_nodes = c('Y_206','Y_12','Y_4'))
    cds_sub2 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_206'), ending_pr_nodes = c(paste0('Y_189')))
    cds_sub3 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_206'), ending_pr_nodes = c(paste0('Y_201')))
    cds_sub4 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_210'), ending_pr_nodes = c(paste0('Y_189')))
    cds_sub5 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_189'), ending_pr_nodes = c(paste0('Y_201')))
    final <- union(colnames(cds_sub),colnames(cds_sub2))
    final <- union(final,colnames(cds_sub3))
    final <- union(final,colnames(cds_sub4))
    final <- union(final,colnames(cds_sub5))
    #final <- union(final,colnames(cds_sub3))
    #print(head(final))
    cds_sub <- cds[,final]
    #cds_nrf2 <- order_cells(cds_sub, root_pr_nodes = c('Y_76'))
    
	time <- read.table('working_data/rna_velo_time.csv',sep = ',', header = TRUE, row.names = 1)
    time <- time[colnames(cds_sub),]
    colData(cds_sub)$time <- time
    cds_nrf2 <- cds_sub
    
    gene <- readRDS(file = 'working_data/monocle3/pr_test_1000_sub_nrf2')
    gene_nrf2 <- row.names(subset(gene, q_value < 10^(-35)))
    gene <- readRDS(file = 'working_data/monocle3/pr_test_1000_sub_tead')
    gene_tead <- row.names(subset(gene, q_value < 10^(-35)))

    gene <- union(gene_nrf2,gene_tead)
    cds_nrf2 <- cds_nrf2[gene,]
    print(cds_nrf2)
    print('finish data preparation')

    trend_formula = "~ splines::ns(time, df=3)"
    model_tbl <- monocle3:::fit_models(cds_nrf2, cores = 1, model_formula_str = trend_formula)
    #print(colData(cds_nrf2)$pseudotime)
    #colData(cds_nrf2)$pseudotime <- monocle3::pseudotime(cds_nrf2)
    newdata <- data.frame(pseudotime = seq(0, max(colData(cds_nrf2)$time), length.out=100))
    model_expectation <- monocle3:::model_predictions(model_tbl, new_data = newdata)
    saveRDS(model_expectation, file = 'working_data/monocle3/nrf2_model_rna_velocity')
    print('nrf2 finish')
}

if(FALSE)
{
    library(monocle3)
	cds <- readRDS('working_data/monocle3/cds_1000_unorder_no_partition')
    gene <- readRDS(file = 'working_data/monocle3/pr_test_1000_sub_nrf2')
    gene_nrf2 <- row.names(subset(gene, q_value < 10^(-35)))

    cds_sub <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_76'), ending_pr_nodes = c('Y_208','Y_12','Y_4'))
    cds_sub2 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_208'), ending_pr_nodes = c(paste0('Y_202')))
    final <- union(colnames(cds_sub),colnames(cds_sub2))
    cds_tead <- cds[,final]
    cds_tead <- order_cells(cds_tead, root_pr_nodes = c('Y_76'))
    
    gene <- readRDS(file = 'working_data/monocle3/pr_test_1000_sub_tead')
    gene_tead <- row.names(subset(gene, q_value < 10^(-35)))

    gene <- union(gene_nrf2,gene_tead)
    cds_tead <- cds_tead[gene,]
    print(cds_tead)
    print('finish data preparation')

    trend_formula = "~ splines::ns(pseudotime, df=3)"
    model_tbl <- monocle3:::fit_models(cds_tead, cores = 1, model_formula_str = trend_formula)
    colData(cds_tead)$pseudotime <- monocle3::pseudotime(cds_tead)
    newdata <- data.frame(pseudotime = seq(0, max(colData(cds_tead)$pseudotime), length.out=100))
    #print(head(newdata))
    model_expectation <- monocle3:::model_predictions(model_tbl, new_data = newdata)
    print(head(model_expectation))
    saveRDS(model_expectation, file = 'working_data/monocle3/tead_model_2')
}
if(FALSE)
{
    library(monocle3)
	cds <- readRDS('working_data/monocle3/cds_1000_unorder_no_partition')
    gene <- readRDS(file = 'working_data/monocle3/pr_test_1000_sub_nrf2')
    gene_nrf2 <- row.names(subset(gene, q_value < 10^(-35)))

    cds_sub <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_76'), ending_pr_nodes = c('Y_208','Y_12','Y_4'))
    cds_sub2 <- choose_graph_segments(cds, clear_cds = FALSE, starting_pr_node = c('Y_208'), ending_pr_nodes = c(paste0('Y_202')))
    final <- union(colnames(cds_sub),colnames(cds_sub2))
    cds_tead <- cds[,final]
    cds_tead <- order_cells(cds_tead, root_pr_nodes = c('Y_76'))
    
	time <- read.table('working_data/rna_velo_time.csv',sep = ',', header = TRUE, row.names = 1)
    time <- time[colnames(cds_tead),]
    colData(cds_tead)$time <- time
    names(time) <- colnames(cds_tead)
    print(head(time))
    cds_tead@principal_graph_aux[['UMAP']]$pseudotime <- time
    #print(head(cds_nrf2@principal_graph_aux[['UMAP']]$pseudotime))
    print(head(cds_tead@principal_graph_aux[['UMAP']]$pseudotime))
    
    gene <- readRDS(file = 'working_data/monocle3/pr_test_1000_sub_tead')
    gene_tead <- row.names(subset(gene, q_value < 10^(-35)))

    gene <- union(gene_nrf2,gene_tead)
    cds_tead <- cds_tead[gene,]
    #cds_tead <- cds_tead['GPX4',]
    #print(cds_tead)
    print('finish data preparation')
    

    trend_formula = "~ splines::ns(pseudotime, df=3)"

    #print(colData(cds_tead)$time)
    #model_tbl <- monocle3:::fit_models(cds_tead, cores = 1, model_formula_str = trend_formula)
    model_tbl <- fit_models(cds_tead, cores = 1, model_formula_str = trend_formula)
    #colData(cds_tead)$pseudotime <- monocle3::pseudotime(cds_tead)
    newdata <- data.frame(pseudotime = seq(0, max(colData(cds_tead)$time), length.out=100))
    print(head(newdata))
    model_expectation <- monocle3:::model_predictions(model_tbl, new_data = newdata)
    print(head(model_expectation))
    saveRDS(model_expectation, file = 'working_data/monocle3/tead_model_rna_velocity')
}


#plotting heatmap
if(FALSE)
{
    library(monocle3)
    nrf2 <- readRDS(file = 'working_data/monocle3/nrf2_model_2')
    tead <- readRDS(file = 'working_data/monocle3/tead_model_2')
    print(dim(tead))
    print(dim(nrf2))
    colnames(tead) <- c(101:200)
    pred <- cbind( nrf2[,100:1], tead)
    #saveRDS(pred, file = 'time')
    #pred <- Matrix::t(scale(Matrix::t(pred),center=TRUE))
    #print(head(pred))
    #print(pred)
    row_dist <- as.dist((1-cor(Matrix::t(pred)))/2)
    #row_dist[is.na(row_dist)] <- 1
    pred <- (pred - rowMeans(pred)) / apply(pred,1,sd)
    pred[pred>3] <- 3
    pred[pred<(-3)] <- (-3)
    labrow <- rownames(tead)
    labrow[ !(labrow %in% c('GPX4','ACSL4','TFRC', 'FTH1', 'FTL', 'AKT3')) ] <- ''
    pheatmap::pheatmap(
                pred, 
                cluster_cols = FALSE, 
                #color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(300),
                #breaks= seq(-3,3,length.out=301),
                #scale = 'row',
                file = 'branched_1000.pdf',
                cluster_rows = TRUE, 
                #show_rownames = T, 
                show_colnames = F, 
                #fontsize_row = 0.5,
                labels_row = labrow,
                clustering_distance_rows = row_dist, 
                hclust_method = 'ward.D2', cutree_rows = 5, silent = TRUE)
}

#plotting heatmap
if(FALSE)
{
    library(monocle3)
    nrf2 <- readRDS(file = 'working_data/monocle3/nrf2_model_rna_velocity')
    tead <- readRDS(file = 'working_data/monocle3/tead_model_rna_velocity')
    print(dim(tead))
    print(dim(nrf2))
    colnames(tead) <- c(101:200)
    pred <- cbind( nrf2[,100:1], tead)
    #saveRDS(pred, file = 'time')
    #pred <- Matrix::t(scale(Matrix::t(pred),center=TRUE))
    #print(head(pred))
    #print(pred)
    row_dist <- as.dist((1-cor(Matrix::t(pred)))/2)
    #row_dist[is.na(row_dist)] <- 1
    pred <- (pred - rowMeans(pred)) / apply(pred,1,sd)
    pred[pred>3] <- 3
    pred[pred<(-3)] <- (-3)
    labrow <- rownames(tead)
    labrow[ !(labrow %in% c('GPX4','ACSL4','TFRC', 'FTH1', 'FTL')) ] <- ''
    pheatmap::pheatmap(
                pred, 
                cluster_cols = FALSE, 
                #color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(300),
                #breaks= seq(-3,3,length.out=301),
                #scale = 'row',
                file = 'branched_rna_velocity.pdf',
                cluster_rows = TRUE, 
                #show_rownames = T, 
                show_colnames = F, 
                #fontsize_row = 0.5,
                labels_row = labrow,
                clustering_distance_rows = row_dist, 
                hclust_method = 'ward.D2', cutree_rows = 5, silent = TRUE)
}

#codes may be useful
#exporting results
#export regulon
if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.rds')
	temp <- loadInt(scenicOptions, 'aucell_regulonAUC')
	temp <- as.data.frame(temp@assays@data@listData$AUC)
	saveRDS(temp, file = 'regulon_auc.rds')
}

#calculate ferro sensitivity using AUC
if( FALSE )
{
	exp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	print(dim(exp))
	print('finish loading')
	set.seed(123)
	library(SCENIC)
	library(AUCell)
	scenicOptions <- readRDS(file = 'int/scenicOptions.rds')
	genesKept <- loadInt( scenicOptions, 'genesKept')
	exp_filter_tmp <- exp[genesKept,]
	exp_filter <- log2(exp_filter_tmp+1)
	rm(exp)
	rm(exp_filter_tmp)
	cell_rankings <- AUCell_buildRankings(exp_filter, nCores = 10)
	saveRDS(cell_rankings, file = 'working_data/cell_rank')
}

if(FALSE)
{
	library(AUCell)
	anti <- read.table('data/antiferroptosis_3.txt',sep = '\t',header = TRUE,stringsAsFactors = FALSE)[,1]
	pro <- read.table('data/proferroptosis_3.txt',sep = '\t',header = TRUE,stringsAsFactors = FALSE)[,1]
	gene <- list(anti = anti, pro = pro)
	cell_rankings <- readRDS('working_data/cell_rank')
	#cells_AUC <- AUCell_calcAUC(gene, cell_rankings, aucMaxRank = 1000, nCores = 10)
	cells_AUC <- AUCell_calcAUC(gene, cell_rankings, aucMaxRank = ceiling(0.15*nrow(cell_rankings)), nCores = 10)
	saveRDS(cells_AUC, file = 'working_data/auc_ferro_test_3')
}

if(FALSE)
{
	library(AUCell)
	obj <- readRDS('working_data/auc_ferro_test_3')
    	mat <- getAUC(obj)[,]
	mat <- t(mat)
	mat <- data.frame(mat)
	mat$all <- mat$anti - mat$pro
	print(head(mat))
	write.table( mat,file = 'auc_ferro_3.csv',quote = F,row.names = T,col.names = T,sep = ',')
	#saveRDS(mat, file = 'working_data/auc_ferro_mat')
}

#plot branched on NMF
if(FALSE)
{

	time <- readRDS(file = 'time')
	print(dim(time))
	print(head(rownames(time)))
	row_dist <- as.dist((1-cor(Matrix::t(time)))/2)
	time <- (time- rowMeans(time)) / apply(time,1,sd)
	print(time[,100])
	time[time>3] <- 3
	time[time<(-3)] <- (-3)
	#time[,1:100] <- (time[,1:100] - rowMeans(time[,1:100])) / apply(time[,1:100],1,sd)
	#time[,101:200] <- (time[,101:200] - rowMeans(time[,101:200])) / apply(time[,101:200],1,sd)
	#NMF::aheatmap(time, breaks = seq( -3,3,length.out = 300), 
	NMF::aheatmap(time, 
		      Rowv = row_dist, 
		      hclustfun = 'ward', Colv = NA, filename = 'test_3.pdf')
}


#SCENIC internal codes for modification
#this part of code is from scenic internal
#plot heatmap
.openDev <- function(fileName, devType, ...)
{
  if(devType=="pdf")
    pdf(paste0(fileName, ".pdf"), ...)
  
  if(devType=="png")
    png(paste0(fileName, ".png", type="cairo"), ...)
  
  if(devType=="cairo_pfd") # similar to Cairo::CairoPDF?
    grDevices::cairo_pdf(paste0(fileName, ".pdf"), ...)
}

.openDevHeatmap <- function(fileName, devType)
{
  if(devType!="pdf") 
  {
    if(devType=="png") .openDev(fileName=fileName, devType=devType, width=1200,height=1200)
    if(devType!="png") .openDev(fileName=fileName, devType=devType)
    fileName <- NA
  }else{
    fileName <- paste0(fileName,".pdf")
  }
  return(fileName)
}

.closeDevHeatmap <- function(devType)
{
  if(devType!="pdf") 
  {
    dev.off()
  }
}

if(FALSE)
{
	library(NMF)
	library(SCENIC)
	library(AUCell)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	skipHeatmap <- FALSE
	regulonAUC <- readRDS(file = getIntName(scenicOptions, 'aucell_regulonAUC'))
msg <- paste0(format(Sys.time(), "%H:%M"), "\tFinished running AUCell.")
  if(getSettings(scenicOptions, "verbose")) message(msg)
  
  if(!skipHeatmap){
    msg <- paste0(format(Sys.time(), "%H:%M"), "\tPlotting heatmap...")
    if(getSettings(scenicOptions, "verbose")) message(msg)
    
    nCellsHeatmap <- min(500, ncol(regulonAUC))
    cells2plot <- sample(colnames(regulonAUC), nCellsHeatmap)
    
    cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")   #TODO check if exists, if not... create/ignore?
    if(!is.null(cellInfo)) cellInfo <- data.frame(cellInfo)[cells2plot,,drop=F]
    colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")
    
    fileName <- getOutName(scenicOptions, "s3_AUCheatmap")
    
    fileName <- .openDevHeatmap(fileName=fileName, devType=getSettings(scenicOptions, "devType"))
    NMF::aheatmap(getAUC(regulonAUC)[,cells2plot],
                  annCol=cellInfo,
                  annColor=colVars,
                  main="AUC",
                  sub=paste("Subset of",nCellsHeatmap," random cells"),
                  filename=fileName)
    .closeDevHeatmap(devType=getSettings(scenicOptions, "devType"))
  }
}

#plot Tsne
if( FALSE)
{
	library(Rtsne)
	library(R2HTML)
	library(SCENIC)
	library(AUCell)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	skipTsne <- FALSE
	regulonAUC <- readRDS(file = getIntName(scenicOptions, 'aucell_regulonAUC'))

	exp_tmp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	exp <- as.matrix(exp_tmp)
	rm(exp_tmp)
	print(dim(exp))
	print('finish loading')

	genesKept <- loadInt( scenicOptions, 'genesKept')
	exp_filter_tmp <- exp[genesKept,]
	exprMat <- log2(exp_filter_tmp+1)
	rm(exp)
	rm(exp_filter_tmp)
if(!skipTsne){
    msg <- paste0(format(Sys.time(), "%H:%M"), "\tPlotting t-SNEs...")
    if(getSettings(scenicOptions, "verbose")) message(msg)
    
    tSNE_fileName <- tsneAUC(scenicOptions, aucType="AUC", onlyHighConf=FALSE) # default: nPcs, perpl, seed, tsne prefix
    tSNE <- readRDS(tSNE_fileName)
    
    # AUCell (activity) plots with the default tsne, as html: 
    fileName <- getOutName(scenicOptions, "s3_AUCtSNE_colAct")
    plotTsne_AUCellHtml(scenicOptions, exprMat, fileName, tSNE) #open the resulting html locally

    # Plot cell properties:
    sub <- ""; if("type" %in% names(tSNE)) sub <- paste0("t-SNE on ", tSNE$type)
    cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null") 
    colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")
    pdf(paste0(getOutName(scenicOptions, "s3_AUCtSNE_colProps"),".pdf"))
    plotTsne_cellProps(tSNE$Y, cellInfo=cellInfo, colVars=colVars, cex=1, sub=sub)
    dev.off()
  }
}
###end scenic internal

########this part of code is just for debugging
####SCENIC checking
#runGenie3 check parrellel on 10
if( FALSE )
{
	library(feather)
	exp <- readRDS(file = './int/1.3_GENIE3_weightMatrix_part_103.Rds')
	print(length(colnames(exp)))
	print(length(rownames(exp)))
	i <- 2
	regulator <- read.csv( paste0('working_data/para/regulator.csv'), header = TRUE, sep = '\t')
	target <- read.csv( paste0('working_data/para/target_', toString(i),'.csv'), header = TRUE, sep = '\t')
	exp <- read_feather( paste0('working_data/para/exp_', toString(i),'.feather'))
	exp <- t(exp)
	print(head(regulator))
	print(head(target))
	print(dim(exp))
	print(head(rownames(exp)))
	print(head(rownames(exp[target[,'x'],])))
}

#retrive parallel
if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	for( i in c(1,2,3,4,5,6,7,8,9,10))
	{
		mat <- read.csv( paste0('working_data/genie3_1-3d_weightmatrix_',toString(i),'.csv'), header = TRUE, row.names = 1, sep = ',')
		print(head(rownames(mat)))
		print(dim(mat))
		saveRDS( mat, file = paste0("1.3_GENIE3_weightMatrix_part_",toString(i),'.Rds'))
	}
}

#rerun 
if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	genesKept <- loadInt( scenicOptions, 'genesKept')
	final <- list()
	if(TRUE)
	{
		final <- readRDS( file = paste0("int/1.3_GENIE3_weightMatrix_part_",toString(1),'.Rds'))
		genesDone <- colnames(final)
		if( FALSE)
		{
		for( i in c(2,3,4,5,6,7,8,9,10))
		{
			mat <- readRDS( file = paste0("int/1.3_GENIE3_weightMatrix_part_",toString(i),'.Rds'))
			final <- cbind(final, mat)
			#print(intersect(genesDone,colnames(mat)))
			#genesDone <- union(colnames(mat), genesDone)
		}
		}
	}
	print(length(genesDone))
	print(length(genesKept))
	genes <- setdiff( genesKept, genesDone)
	#genes <- intersect( genesKept, colnames(final))
	#final <- final[,genes]
	#print(dim(final))
	#saveRDS(final, file = '1.3_GENIE3_weightMatrix_part_1.Rds')
	print(genes)
}

if( FALSE )
{
	final <- readRDS(file = '1.3_GENIE3_weightMatrix_part_1.Rds')
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	regulator <- getDbTfs(scenicOptions)
	genesKept <- loadInt( scenicOptions, 'genesKept')
	exp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	exp <- exp[genesKept,]
	exp <- log2(exp+1)
	regulator <- regulator[regulator %in% rownames(exp)]
	print(length(regulator))
	exp <- t(exp)
	exp <- as.data.frame(exp)
	print(dim(exp))
	write.table( regulator, file = paste0('regulator.csv'), quote = F,row.names = T,col.names = T,sep = '\t')
	target <- setdiff( colnames(exp), colnames(final))
	print(length(target))
	output <- union( regulator, target)
	print(length(output))
	i<-1
	write.table( target, file = paste0('target_test_', toString(i),'.csv'), quote = F,row.names = T,col.names = T,sep = '\t')
	library(feather)
	write_feather( exp[,output], paste0('exp_test_', toString(i),'.feather'))
	if( FALSE )
	{
		increment <- 20
		for( i in c(1:1) )
		{
			start <- (i-1)*increment + 1
			end <- start + increment - 1
			if( end > length(genekept) )
				end <- length(genekept)
			target <- genekept[c(start:end)]
			print(i)
			print(length(target))
			output <- union( regulator, target)
			write.table( target, file = paste0('target_test_', toString(i),'.csv'), quote = F,row.names = T,col.names = T,sep = '\t')
			write_feather( exp[,output], paste0('exp_test_', toString(i),'.feather'))
		}
	}
}

if( FALSE )
{
	library(SCENIC)
	final <- readRDS(file = '1.3_GENIE3_weightMatrix_part_1.Rds')
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	genesKept <- loadInt( scenicOptions, 'genesKept')
	mat <- read.csv( paste0('genie3_1-3d_weightmatrix_test_1.csv'), header = TRUE, row.names = 1, sep = ',')
	colnames(mat) <- c('IQCJ-SCHIP1','NKX2-5','HLA-F','HLA-A','HLA-E','HLA-C','HLA-B','NKX3-1','BDNF-AS')
	print(intersect(genesKept,colnames(mat)))
	final <- cbind(final, mat)
	saveRDS( final, file = paste0("int/1.3_GENIE3_weightMatrix_part_1.Rds"))
}
######debuggin part finish

