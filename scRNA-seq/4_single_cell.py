import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import numpy as np

#import sys
#sys.path.insert(0, './utility/')

def save_obj( data, name ):
    with open(name,'wb') as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL);


def load_obj( name ):
    with open(name,'rb') as f:
        return pickle.load(f);


#remove doublets
'''doub = pd.read_csv('python_data/DoubletFinder.csv', index_col = 0)
doub = doub.loc[doub.loc[:,'DF.classifications_0.25_0.09_2540'] == 'Singlet',:]
doub = doub.loc[~doub.index.str.contains('-1'),:]
doub.index = doub.index.set_names('Cell')
print(doub)
doub.to_csv('R_data/Seurat/Singlet.csv')
#'''

#output RNA velocity time info
'''velo = pd.read_pickle('python_data/rna_all_once.pkl')
print(velo)
timeinfo = pd.DataFrame( velo.loc[:,'latent_time'].values, index = velo.index, columns = ['latent'])
timeinfo.loc[timeinfo.index.str.contains('-2'),'day'] = 1
timeinfo.loc[timeinfo.index.str.contains('-3'),'day'] = 2
timeinfo.loc[timeinfo.index.str.contains('-4'),'day'] = 3
print(timeinfo)
timeinfo.T.to_csv('R_data/Regulon/rna_velocity_latent_time.csv')
#'''

#removing MT RPL RPS gene
'''x = pd.read_csv('python_data/ribosomal_gene.txt', sep = '\t')
print(x)
gene = pd.read_csv('R_data/genesKept.csv', sep = ',', index_col = 0)
print(gene)
gene = gene.loc[~gene.iloc[:,0].isin(x.loc[:,'Approved symbol']),:]
gene = gene.loc[~gene.iloc[:,0].str.contains('MT-'),:]
print(gene)
gene.to_csv('R_data/genesKept_trim_2.csv')
for i in ['1000']:
    gene = pd.read_csv('R_data/monocle3/sc_vargene_'+i+'.csv', sep = '\t', index_col = 0)
    gene = gene.loc[~gene.iloc[:,0].isin(x.loc[:,'Approved symbol']),:]
    gene = gene.loc[~gene.iloc[:,0].str.contains('MT-'),:]
    print(gene)
    gene.to_csv('R_data/sc_var_'+i+'_2.csv')
#'''
'''gene = pd.read_csv('R_data/Regulon/genesKept.csv', sep = ',', index_col = 0)
gene = gene.loc[~gene.iloc[:,0].str.contains('MT-'),:]
gene = gene.loc[~gene.iloc[:,0].str.contains('RPL'),:]
gene = gene.loc[~gene.iloc[:,0].str.contains('RPS'),:]
print(gene)
gene.to_csv('R_data/genesKept_trim.csv')
for i in ['1000']:
    gene = pd.read_csv('R_data/monocle3/sc_vargene_'+i+'.csv', sep = '\t', index_col = 0)
    gene = gene.loc[~gene.iloc[:,0].str.contains('MT-'),:]
    gene = gene.loc[~gene.iloc[:,0].str.contains('RPL'),:]
    gene = gene.loc[~gene.iloc[:,0].str.contains('RPS'),:]
    print(gene)
    gene.to_csv('R_data/sc_var_'+i+'.csv')
#'''

#scvelo
'''import scvelo as scv
print('importing scvelo')
adata = scv.read('../Jen-Hao/RNA_velocity/M4340-43_rmDoub_dynamical.h5ad', cache = True)
print(adata)
#'''

###subsetting data and integrate UMAP
#clustdf = pd.read_csv('../Jen-Hao/RNA_velocity/N417_cell_cluster.tsv', index_col=0, sep='\t')
#print(clustdf)
'''genes = pd.read_csv( 'R_data/genesKept_trim_2.csv', index_col = 0, sep = ',')
genes = genes.loc[:,'x'].tolist()
#print(genes)
monoumap = pd.read_csv('R_data/monocle3/monocle_umap.csv', sep = '\t', index_col = 0)
clustdf = monoumap
newname = []
cat_list = ['M4340', 'M4341', 'M4342', 'M4343']
for i in clustdf.index:
    barcode = i.split('-')[0]
    sample = i.split('-')[1]
    if sample == '2':
        sample = 'M4341'
    if sample == '3':
        sample = 'M4342'
    if sample == '4':
        sample = 'M4343'
    
    key = ''
    
    if sample in cat_list:
        if len(cat_list) == 1:
            key = sample+'_data:'+barcode+'x'
        elif cat_list.index(sample) == len(cat_list)-1:
            key = sample+'_data:'+barcode+'x-1'
        elif cat_list.index(sample) == 0:
            key = sample+'_data:'+barcode+'x'
            for j in range(len(cat_list)-1):
                key = key + '-0'
        else:
            key = sample+'_data:'+barcode+'x-1'
            for j in range(len(cat_list)-cat_list.index(sample)-1):
                key = key + '-0'
    newname.append(key)

clustdf.index = newname
print(clustdf)
clustdf.to_pickle('python_data/monocle_umap.pkl')
if True:
    import scvelo as scv
    adata = scv.read('../../../../Jen-Hao/RNA_velocity/M4340-43_rmDoub_dynamical.h5ad', cache = True)
    #scv.pp.filter_genes(adata, retain_genes = genes)
    #print('retain genes')
    adata = adata[clustdf.index]
    print(adata)
    adata.obsm['X_umap'] = clustdf.loc[adata.obs_names].values
    print('umap')
    adata.write('python_data/new_rna_velo_model')
#'''

#trim to 9000 genes
'''import scvelo as scv
genes = pd.read_csv( 'R_data/genesKept_trim_2.csv', index_col = 0, sep = ',')
genes = genes.loc[:,'x'].tolist()
print(genes[0:10])
adata = scv.read('python_data/new_rna_velo_model.h5ad', cache = True)
#print(adata.var)
x = adata.var.index.tolist()
x = set(genes) & set(x)
adata = adata[:,list(x)]
print(adata)
adata.write('python_data/new_rna_velo_model_trim.h5ad')
#'''

#rerun rna velocity from the beginning
'''import scvelo as scv
adata = scv.read('python_data/new_rna_velo_model_trim.h5ad', cache = True)
scv.utils.cleanup(adata, clean = 'layers')
scv.utils.cleanup(adata, clean = 'obs')
scv.utils.cleanup(adata, clean = 'var')
scv.utils.cleanup(adata, clean = 'uns')
del getattr(adata, 'uns')['neighbors']
del getattr(adata, 'obsm')['X_pca']
del getattr(adata, 'obsm')['X_umap']
del getattr(adata, 'obsm')['velocity_umap']
del getattr(adata, 'varm')['PCs']
del getattr(adata, 'varm')['loss']
del getattr(adata, 'layers')['Ms']
del getattr(adata, 'layers')['Mu']
del getattr(adata, 'obsp')['connectivities']
del getattr(adata, 'obsp')['distances']
#genes = pd.read_csv('working_data/sc_var_1000.csv', index_col = 0)
#genes = genes.loc[:,'x'].tolist()
#print(genes)
#x = adata.var.index.tolist()
#x = set(genes) & set(x)
#adata = adata[:,list(x)]
#print(adata)
scv.pp.moments(adata, n_pcs = 50, n_neighbors = 30)
print('m')
adata.write('python_data/small_rna_velo_model_m.h5ad')
#'''

#construct different velocity parameters
'''import scvelo as scv
adata = scv.read('python_data/small_rna_velo_model_m.h5ad', cache = True)
scv.tl.recover_dynamics(adata, n_jobs = 16)
print('d')
adata.write('python_data/small_rna_velo_model_d.h5ad')
scv.tl.velocity(adata, mode = 'dynamical')
scv.tl.velocity_graph(adata)
print('g')
adata.write('python_data/small_rna_velo_model_g.h5ad')
#'''
'''import scvelo as scv
adata = scv.read('python_data/small_rna_velo_model_g.h5ad', cache = True)
scv.tl.velocity_graph(adata, n_jobs = 15)
print(adata)
#scv.pl.velocity_embedding_stream(adata, basis = 'umap', save = 'figure/velo/newumap_stream_recalculate.pdf')
#scv.pl.velocity_embedding(adata, density = 0.2, basis = 'umap', save = 'figure/velo/newumap_embedding_recalculate.pdf')
#scv.pl.velocity_grid(adata, basis = 'umap', save = 'figure/velo/newumap_embedding_grid_recalculate.pdf')
#print(adata[:,'GPX4'].layers['Ms'])
#adata.obsm['X_umap_new'] = adata.obsm['X_umap']
#scv.pl.velocity_embedding_grid(adata, density = 1, basis = 'umap_new', save = 'figure/velo/newumap_6.pdf')
#scv.pl.velocity_embedding_stream(adata, basis = 'umap_new', save = 'figure/velo/newumap_stream_recalculate.pdf')
#scv.pl.velocity_embedding_grid(adata, basis = 'umap', save = 'figure/velo/newumap_grid.pdf')
#'''

#try out stochastic
'''import scvelo as scv
adata = scv.read('python_data/small_rna_velo_model_d.h5ad', cache = True)
clustdf = pd.read_pickle('python_data/monocle_umap.pkl')
adata.obsm['X_umap_sm'] = clustdf.loc[adata.obs_names].values
#scv.tl.velocity(adata, mode = 'stochastic')
scv.tl.velocity(adata, n_jobs = 15)
adata.write('python_data/small_rna_velo_model_g_1000_umap.h5ad')
#scv.tl.velocity(adata, mode = 'deterministic')
#scv.tl.velocity_graph(adata, basis = 'umap_sm')
#scv.tl.velocity_graph(adata)
#scv.pl.velocity_embedding_stream(adata, basis = 'umap_sm', save = 'newumap_stream_recalculate_diff_bas_2.png')
#scv.pl.velocity_embedding_grid(adata, basis = 'umap_sm', save = 'newumap_embedding_grid_recalculate_diff_bas_2.pdf')
#'''

'''import scvelo as scv
adata = scv.read('python_data/small_rna_velo_model_g.h5ad', cache = True)
clustdf = pd.read_pickle('python_data/monocle_umap.pkl')
adata.obsm['X_umap_sm'] = clustdf.loc[adata.obs_names].values
scv.tl.velocity_graph(adata, n_jobs = 15)
scv.pl.velocity_embedding_stream(adata, basis = 'umap_sm', save = 'newumap_stream_sm_recalculate_1000.png')
scv.pl.velocity_embedding(adata, density = 0.2, basis = 'umap_sm', save = 'newumap_embedding_sm_recalculate_1000.pdf')
scv.pl.velocity_embedding_grid(adata, basis = 'umap_sm', save = 'newumap_embedding_grid_sm_recalculate_1000.pdf')
#'''

'''import scvelo as scv
adata = scv.read('python_data/small_rna_velo_model_g.h5ad', cache = True)
clustdf = pd.read_pickle('python_data/monocle_umap.pkl')
adata.obsm['X_umap_sm'] = clustdf.loc[adata.obs_names].values
#scv.tl.velocity(adata, mode = 'stochastic')
scv.tl.velocity(adata, n_jobs = 15)
#scv.tl.velocity(adata, mode = 'deterministic')
scv.tl.velocity_graph(adata, basis = 'umap_sm')
#scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis = 'umap_sm', save = 'newumap_stream_recalculate_diff_bas_2.png')
scv.pl.velocity_embedding_grid(adata, basis = 'umap_sm', save = 'newumap_embedding_grid_recalculate_diff_bas_2.pdf')
#'''

#try out different clusters
'''import scvelo as scv
import scanpy as sc
adata = scv.read('data/small_rna_velo_model_g_1000.h5ad', cache = True)
#sc.tl.leiden(adata)
scv.tl.louvain(adata, resolution = 2, random_state = 10)
clustdf = pd.read_pickle('data/monocle_umap.pkl')
adata.obsm['X_umap_sm'] = clustdf.loc[adata.obs_names].values
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
#scv.tl.paga(adata, groups = 'leiden')
scv.tl.paga(adata, groups = 'louvain')
scv.pl.paga(adata, basis = 'umap_sm', size =50, alpha = .1, min_edge_width=2, node_size_scale = 1.5, save = 'trajectory_paga_louvain.pdf')
#scv.tl.velocity(adata, mode = 'deterministic')
#scv.tl.velocity_graph(adata, basis = 'umap_sm')
#scv.tl.velocity_graph(adata)
#scv.pl.velocity_embedding_stream(adata, basis = 'umap_sm', save = 'newumap_stream_recalculate_diff_bas_2.png')
#scv.pl.velocity_embedding_grid(adata, basis = 'umap_sm', save = 'newumap_embedding_grid_recalculate_diff_bas_2.pdf')
#'''

#using different root region
'''import scvelo as scv
adata = scv.read('python_data/small_rna_velo_model_g.h5ad', cache = True)
clustdf = pd.read_pickle('python_data/monocle_umap.pkl')
#print(clustdf)
#root = ((clustdf.loc[:,'UMAP1'] > -6) & (clustdf.loc[:,'UMAP1'] <-5) & ((clustdf.loc[:,'UMAP2'] > -1)|(clustdf.loc[:,'UMAP2']<-4))).astype('double')
root = ((clustdf.loc[:,'UMAP1'] > -7) & (clustdf.loc[:,'UMAP1'] <-5) & ((clustdf.loc[:,'UMAP2'] > 0))).astype('double')
#root = ((clustdf.loc[:,'UMAP1'] <-5) ).astype('double')
#root = ((clustdf.loc[:,'UMAP1'] > -7) & (clustdf.loc[:,'UMAP1'] <-5)).astype('double')
#print(root)
#clustdf
adata.obsm['X_umap_sm'] = clustdf.loc[adata.obs_names].values
#day = [int(i.split('-')[1])-1 for i in adata.obs_names]
day = []
for i in adata.obs_names:
    sample = i.split('_')[0]
    if sample == 'M4341':
        day.append('1')
    if sample == 'M4342':
        day.append('2')
    if sample == 'M4343':
        day.append('3')
#print(day)
adata.obs['day'] = day
#scv.tl.terminal_states(adata)
scv.tl.terminal_states(adata, **{'basis':'umap_sm'})
tmp = adata.obs['root_cells']
tmp = tmp * root
#tmp[tmp.index.str.contains('M4343')] = 0
#tmp[tmp.index.str.contains('M4342')] = 0
#adata.obs['new_root'] = root.values
adata.obs['new_root'] = tmp.values
#adata.obs['new_end'] = 1-tmp
#print(adata.obs['new_root'])
#scv.tl.latent_time(adata, root_key = 'new_root', end_key = 'new_end')
scv.tl.latent_time(adata, root_key = 'new_root')
#scv.tl.latent_time(adata)
scv.pl.scatter(adata, color = 'latent_time', basis = 'umap_sm', save = 'latent_time_all.pdf')

scv.tl.velocity(adata, n_jobs = 15)
scv.tl.velocity_graph(adata, n_jobs = 15)
#scv.pl.velocity_embedding_stream(adata, basis = 'umap_sm', save = 'newumap_stream_sm_recalculate_1000.png')
#scv.pl.velocity_embedding(adata, density = 0.2, basis = 'umap_sm', save = 'newumap_embedding_sm_recalculate_1000.pdf')
#scv.pl.velocity_embedding_grid(adata, basis = 'umap_sm', save = 'newumap_embedding_grid_sm_recalculate_1000.pdf')
scv.pl.velocity_embedding_grid(adata, basis = 'umap_sm', density = 0.5, color = 'latent_time', scale = 0.7, save = 'newumap_embedding_grid_sm_recalculate_1000_time_3.pdf')
scv.pl.velocity_embedding_grid(adata, basis = 'umap_sm', density = 0.5, color = 'day', scale = 0.7, save = 'newumap_embedding_grid_sm_recalculate_1000_day.pdf')

#scv.pl.scatter(adata, c = ['root_cells','end_points'], basis = 'umap_sm', save = 'root_end_cells_all.pdf')
#scv.tl.latent_time(adata, root_key = 200)
#scv.tools.velocity_embedding(adata, basis = 'umap_sm', vkey = 'velocity')
#tmpx = clustdf.loc[cell.index].values
#vemb = np.array(adata.obsm['velocity_umap_sm'][adata.obs.index.isin(cell.index),:2])
#print(tmpx.shape)
#print(vemb.shape)
#scv.pl.velocity_embedding_grid(adata, basis = 'umap_sm', groups = ['1'], color = 'regu_c', density = 0.6, scale = 0.4, save = 'newumap_embedding_grid_sm_recalculate_1000_time.pdf')
#scv.pl.velocity_embedding_grid(adata, basis = 'umap_sm', X = tmpx, V = vemb, density = 0.5, color = 'lat_time', scale = 0.3, save = 'newumap_embedding_grid_sm_recalculate_1000_time.pdf')
#
#print(adata[:,'GPX4'].layers['Ms'])
#adata.obsm['X_umap_new'] = adata.obsm['X_umap']
#scv.pl.velocity_embedding_grid(adata, density = 1, basis = 'umap_new', save = 'newumap_6.pdf')
#scv.pl.velocity_embedding_stream(adata, basis = 'umap_new', save = 'newumap_stream_recalculate.pdf')
#scv.pl.velocity_embedding_grid(adata, basis = 'umap', save = 'newumap_grid.pdf')
#'''

#getting the latent time
'''import scvelo as scv
adata = scv.read('python_data/small_rna_velo_model_g.h5ad', cache = True)
clustdf = pd.read_pickle('python_data/monocle_umap.pkl')
print(clustdf)
root = ((clustdf.loc[:,'UMAP1'] > -7) & (clustdf.loc[:,'UMAP1'] <-5) & ((clustdf.loc[:,'UMAP2'] > 0))).astype('double')
adata.obsm['X_umap_sm'] = clustdf.loc[adata.obs_names].values
scv.tl.terminal_states(adata, **{'basis':'umap_sm'})
tmp = adata.obs['root_cells']
tmp = tmp * root
#tmp[tmp.index.str.contains('M4343')] = 0
#tmp[tmp.index.str.contains('M4342')] = 0
#adata.obs['new_root'] = root.values
adata.obs['new_root'] = tmp.values
#adata.obs['new_end'] = 1-tmp
#print(adata.obs['new_root'])
#scv.tl.latent_time(adata, root_key = 'new_root', end_key = 'new_end')
scv.tl.latent_time(adata, root_key = 'new_root')
#scv.tl.latent_time(adata)
adata.obs['latent_time'].to_pickle('python_data/latent_time.pkl')
#'''

#help plot R pseudo
'''c = pd.read_csv('R_data/monocle3/branched_gene_cluster_all.csv', sep = '\t', index_col = 0)
nrf2 = pd.read_csv('R_data/monocle3/upper_model_all.csv', sep = ',', index_col = 0)
tead = pd.read_csv('R_data/monocle3/lower_model_all.csv', sep = ',', index_col = 0)
nrf2.columns = range(1,101)
tead.columns = range(101,201)
nrf2 = pd.concat([nrf2,tead], axis = 1)
nrf2 = nrf2.T
nrf2 = (nrf2 - nrf2.mean())/nrf2.std()
nrf2 = nrf2.T
nrf2.loc[:,'cluster'] = c.iloc[:,0]
nrf2 = nrf2.groupby('cluster').mean()
nrf2 = nrf2.T
x = [i for i in range(1,101)]
x_2 = [i for i in range(101,201)]
print(nrf2)
nrf2.loc[x,'time'] = list(range(1,101))
nrf2.loc[x_2,'time'] = list(range(1,101))
nrf2.loc[x,'group'] = 'S'
nrf2.loc[x_2,'group'] = 'R'
print(nrf2)
nrf2.to_csv('R_data/monocle3/module_pseudo_drawing.csv')
#'''

#help plot R pseudo #2
'''c = pd.read_csv('R_data/monocle3/branched_gene_cluster_all.csv', sep = '\t', index_col = 0)
nrf2 = pd.read_csv('R_data/monocle3/upper_model_all.csv', sep = ',', index_col = 0)
tead = pd.read_csv('R_data/monocle3/lower_model_all.csv', sep = ',', index_col = 0)
nrf2.columns = range(1,101)
tead.columns = range(101,201)
nrf2 = pd.concat([nrf2,tead], axis = 1)
gene = nrf2.index
nrf2 = nrf2.T
nrf2 = (nrf2 - nrf2.mean())/nrf2.std()
nrf2 = nrf2.T
nrf2.loc[:,'cluster'] = c.iloc[:,0]
x = [i for i in range(1,101)]
x_2 = [i for i in range(101,201)]
nrf2 = nrf2.T
print(nrf2)
nrf2.loc[x,'time'] = list(range(1,101))
nrf2.loc[x_2,'time'] = list(range(1,101))
nrf2.loc[x,'group'] = 'S'
nrf2.loc[x_2,'group'] = 'R'
print(nrf2)
output = []
for i in gene:
    tmp = nrf2.loc[list(range(1,201)),[i,'time','group']]
    tmp.columns = ['exp','time','group']
    tmp.loc[:,'cluster'] = nrf2.loc['cluster',i]
    tmp.loc[:,'gene'] = i
    output.append(tmp)
output = pd.concat(output, axis = 0)
output.index = list(range(len(output.index)))
print(output)
output.to_csv('R_data/monocle3/module_pseudo_drawing_se.csv')
#'''

#regulon network construction
'''import pyreadr
ncoef = pyreadr.read_r('python_data/lower_trajectory_count.rds')
for key, value in ncoef.items():
    print(value)
    value.to_pickle('python_data/lower_trajectory_count.pkl')
ncoef = pyreadr.read_r('python_data/lower_trajectory_time.rds')
for key, value in ncoef.items():
    print(value)
    print(key)
    value.to_pickle('python_data/lower_trajectory_time.pkl')
#'''

'''import pyreadr
ncoef = pyreadr.read_r('python_data/higher_trajectory_count.rds')
for key, value in ncoef.items():
    print(value)
    value.to_pickle('python_data/higher_trajectory_count.pkl')
ncoef = pyreadr.read_r('python_data/higher_trajectory_time.rds')
for key, value in ncoef.items():
    print(value)
    print(key)
    value.to_pickle('python_data/higher_trajectory_time.pkl')
#'''

'''target = pd.read_csv('output/Step2_regulonTargetsInfo.tsv', sep = '\t')
print(target)
count = pd.read_pickle('python_data/higher_trajectory_count.pkl')
print(count)
tf = pd.DataFrame('No',index = count.index, columns = ['NFE2L2','ATF4','TEAD1','YAP1'])
for i in tf.columns:
    tmp = target.loc[target.loc[:,'TF'] == i,]
    tf.loc[tf.index.isin(tmp.loc[:,'gene']),i] = 'Yes'
    print(tf.loc[:,i].value_counts())
#print(target.loc[target.loc[:,'TF'] == 'NFE2L2',:])
#gene = pd.read_csv('python_data/YAP1_targets_plos.txt')
gene = pd.read_csv('python_data/TF_target_YAP.csv')
tf.loc[ tf.index.isin(gene.iloc[:,0]),'YAP1'] = 'Yes'
print(tf)
tf.to_csv('python_data/TF_info.csv')
#'''

#derive branch specific TF to put in tf_ratio
'''target = pd.read_csv('output/Step2_regulonTargetsInfo.tsv', sep = '\t')
diff = pd.read_csv('R_data/monocle3/gene_annotation.csv', sep = ',')
print(diff)
dtf = target.loc[target.loc[:,'gene'].isin(diff.loc[:,'name']),:]
utf = target.loc[target.loc[:,'gene'].isin(diff.loc[diff.loc[:,'anno']=='upper','name']),:]
ltf = target.loc[target.loc[:,'gene'].isin(diff.loc[diff.loc[:,'anno']=='lower','name']),:]

atf = target.loc[:,'TF'].value_counts()
dtf = dtf.loc[:,'TF'].value_counts()
utf = utf.loc[:,'TF'].value_counts()
ltf = ltf.loc[:,'TF'].value_counts()

atf = pd.DataFrame(atf.values, index = atf.index, columns = ['cnt_2'])
dtf = pd.DataFrame(dtf.values, index = dtf.index, columns = ['cnt'])
utf = pd.DataFrame(utf.values, index = utf.index, columns = ['cnt_u'])
ltf = pd.DataFrame(ltf.values, index = ltf.index, columns = ['cnt_l'])

tf = pd.concat([atf, dtf, utf, ltf], axis = 1)
tf.loc[:,'ratio'] = tf.loc[:,'cnt']/tf.loc[:,'cnt_2']
tf.loc[:,'ratio_u'] = tf.loc[:,'cnt_u']/tf.loc[:,'cnt_2']
tf.loc[:,'ratio_l'] = tf.loc[:,'cnt_l']/tf.loc[:,'cnt_2']
tf = tf.loc[tf.loc[:,'cnt_2'] > 10,:]
#print(tf)
#tf = tf.sort_values(by = ['ratio'], ascending = False)
#print(tf.nlargest(50,columns = ['ratio']))
#print(tf.index.get_loc('NFE2L2'))
tf = tf.sort_values(by = ['ratio_u'], ascending = False)
print(tf.nlargest(50,columns = ['ratio_u']))
print(tf.index.get_loc('NFE2L2'))
tf = tf.sort_values(by = ['ratio_l'], ascending = False)
print(tf.nlargest(50,columns = ['ratio_l']))
#print(tf.index.get_loc('NFE2L2'))
tf.fillna(0, inplace = True)
tf.to_csv('python_data/tf_ratio.csv')
#'''

#heatmap
'''import pyreadr
ncoef = pyreadr.read_r('R_data/Regulon/regulon_auc.rds')
for key, value in ncoef.items():
    regulon = value
regulon = regulon.loc[['NFE2L2_extended (53g)','ATF4_extended (276g)','TEAD1_extended (107g)'],:]
print(regulon)
regulon.index = ['NRF2','ATF4','TEAD1']
regulon_high = regulon.copy()
regulon_low = regulon.copy()
divide = 50
cut = 30
count = pd.read_pickle('python_data/higher_trajectory_count.pkl')
time = pd.read_pickle('python_data/higher_trajectory_time.pkl')
regulon_high = regulon_high.reindex(columns = count.columns)
time.index = count.columns
time.columns = ['time']
timemin = time.loc[:,'time'].min()
timemax = time.loc[:,'time'].max()
step = (timemin - timemax) / divide
time.loc[:,'group'] = (time.loc[:,'time'] - timemin)//step-1
timecount = time.loc[:,'group'].value_counts()
timecount = timecount[timecount > cut]
count = count.T
count = pd.concat([count, time.loc[:,'group']], axis = 1)
regulon_high = regulon_high.T
regulon_high = pd.concat([regulon_high,time.loc[:,'group']], axis = 1)
count = count.loc[count.loc[:,'group'].isin(timecount.index),:]
regulon_high = regulon_high.loc[regulon_high.loc[:,'group'].isin(timecount.index),:]
count_2 = count.groupby('group').mean()
regulon_high = regulon_high.groupby('group').mean()
#print(count_2)

count = pd.read_pickle('python_data/lower_trajectory_count.pkl')
time = pd.read_pickle('python_data/lower_trajectory_time.pkl')
regulon_low = regulon_low.reindex(columns = count.columns)
time.index = count.columns
time.columns = ['time']
timemin = time.loc[:,'time'].min()
timemax = time.loc[:,'time'].max()
step = (timemax - timemin) / divide
time.loc[:,'group'] = (time.loc[:,'time'] - timemin)//step
timecount = time.loc[:,'group'].value_counts()
timecount = timecount[timecount > cut]
count = count.T
count = pd.concat([count, time.loc[:,'group']], axis = 1)
regulon_low = regulon_low.T
regulon_low = pd.concat([regulon_low,time.loc[:,'group']], axis = 1)
count = count.loc[count.loc[:,'group'].isin(timecount.index),:]
regulon_low = regulon_low.loc[regulon_low.loc[:,'group'].isin(timecount.index),:]
count = count.groupby('group').mean()
regulon_low = regulon_low.groupby('group').mean()
#print(count)
count = pd.concat([count_2,count], axis = 0)
#print(count)
count = (count - count.mean()) / count.std()
count = count.T
print(count)
count.to_csv('R_data/python_heatmap_diff_gene_'+str(divide)+'.csv')
regulon = pd.concat([regulon_high,regulon_low], axis = 0)
regulon = (regulon - regulon.mean()) / regulon.std()
regulon = regulon.T
print(regulon)
regulon.to_csv('R_data/python_heatmap_TF_'+str(divide)+'.csv')
#'''

#regulon network annotation
'''
clustdf = pd.read_csv('R_data/regulon_auc_for_python.csv', sep = '\t', index_col = 0)
clustdf = clustdf.T
print(clustdf)
newname = []
cat_list = ['M4340', 'M4341', 'M4342', 'M4343']
for i in clustdf.index:
    barcode = i.split('-')[0]
    sample = i.split('-')[1]
    if sample == '2':
        sample = 'M4341'
    if sample == '3':
        sample = 'M4342'
    if sample == '4':
        sample = 'M4343'
    
    key = ''
    
    if sample in cat_list:
        if len(cat_list) == 1:
            key = sample+'_data:'+barcode+'x'
        elif cat_list.index(sample) == len(cat_list)-1:
            key = sample+'_data:'+barcode+'x-1'
        elif cat_list.index(sample) == 0:
            key = sample+'_data:'+barcode+'x'
            for j in range(len(cat_list)-1):
                key = key + '-0'
        else:
            key = sample+'_data:'+barcode+'x-1'
            for j in range(len(cat_list)-cat_list.index(sample)-1):
                key = key + '-0'
    newname.append(key)

clustdf.index = newname
clustdfall = clustdf.loc[:,~clustdf.columns.str.contains('extended')]
clustdf = clustdf.loc[:,clustdf.columns.str.contains('extended')]
colname = []
for i in clustdf.columns:
    tmp = i.split('_')[0]
    colname.append(tmp)
clustdf.columns = colname
extend = clustdf.copy()
clustdf = clustdfall.copy()
colname = []
for i in clustdf.columns:
    tmp = i.split(' ')[0]
    colname.append(tmp)
clustdf.columns = colname
clustdf = clustdf.loc[:,~clustdf.columns.isin(extend.columns)]
clustdf = pd.concat([extend,clustdf], axis = 1)
clustdf = clustdf.loc[:,~clustdf.columns.duplicated()]

print(clustdf)
#clustdf = (clustdf - clustdf.min())/(clustdf.max() - clustdf.min())
clustdf = (clustdf - clustdf.mean()) / clustdf.std()
print(clustdf)
clustdf = (clustdf - clustdf.min())/(clustdf.max() - clustdf.min())
print(clustdf)
clustdf = clustdf/clustdf.sum()
latent = pd.read_pickle('python_data/latent_time.pkl')
latent = latent.reindex(clustdf.index)
latent = latent.values.tolist()
clustdf = clustdf.T
clustdf = clustdf*latent
clustdf = clustdf.T
weight_time = clustdf.sum()

tf = pd.read_csv('python_data/tf_ratio.csv', index_col = 0)
tfsel = pd.read_csv('R_data/regulon_name.csv', index_col = 0)
tf = tf.loc[tf.index.isin(tfsel.iloc[:,0]),:]
print(tf)
tf.loc[:,'weight_time'] = weight_time[tf.index]
#tf.loc[:,'u_z'] = (tf.loc[:,'ratio_u'] - tf.loc[:,'ratio_u'].mean()) / tf.loc[:,'ratio_u'].std()
#tf.loc[:,'l_z'] = (tf.loc[:,'ratio_l'] - tf.loc[:,'ratio_l'].mean()) / tf.loc[:,'ratio_l'].std()
#tf.loc[:,'ratio_diff_u'] = tf.loc[:,'u_z'] - tf.loc[:,'l_z']
#tf.loc[:,'ratio_diff_l'] = -tf.loc[:,'ratio_diff_u']
#tf.loc[tf.loc[:,'ratio_diff_u'] < -2, 'ratio_diff_u'] = -2
#tf.loc[tf.loc[:,'ratio_diff_u'] > 2, 'ratio_diff_u'] = 3.5
#tf.loc[tf.loc[:,'ratio_diff_l'] < -2, 'ratio_diff_l'] = -2
#tf.loc[tf.loc[:,'ratio_diff_u'] > 2, 'ratio_diff_u'] = 3.5
#tf.loc[:,'alpha_u'] = tf.loc[:,'ratio_u']
#tf.loc[:,'alpha_l'] = tf.loc[:,'ratio_l']
tf.loc[:,'u_name'] = tf.index
tf.loc[:,'l_name'] = tf.index
tf.loc[:,'all_name'] = tf.index
#tf.loc[ tf.loc[:,'ratio_u'] > 0.1,'alpha_u'] = 0.3
#tf.loc[ tf.loc[:,'ratio_u'] > 0.1,'alpha_u'] = 0.3
tf.loc[:,'alpha_u'] = tf.loc[:,'ratio_l']
tf.loc[:,'alpha_l'] = tf.loc[:,'ratio_u']
tf.loc[tf.loc[:,'ratio']<0.2,'all_name'] = ''

large_index = tf.loc[:,'alpha_u'].nlargest(20).index
tf.loc[large_index,'alpha_u'] = 0.5
tf.loc[ ~tf.index.isin(large_index),'u_name'] = ''

large_index = tf.loc[:,'alpha_l'].nlargest(20).index
tf.loc[large_index,'alpha_l'] = 0.5
tf.loc[ ~tf.index.isin(large_index),'l_name'] = ''

tf.loc[:,'alpha_u'] = (tf.loc[:,'alpha_u'] - tf.loc[:,'alpha_u'].min()) / (tf.loc[:,'alpha_u'].max() - tf.loc[:,'alpha_u'].min())
tf.loc[:,'alpha_l'] = (tf.loc[:,'alpha_l'] - tf.loc[:,'alpha_l'].min()) / (tf.loc[:,'alpha_l'].max() - tf.loc[:,'alpha_l'].min())
print(tf)
tf.to_csv('python_data/tf_ratio_anno_20_sel.csv')
#print(weight_time.iloc[0:50])
#print(weight_time.iloc[51:100])
#print(weight_time.iloc[100:150])
#print(weight_time.iloc[150:200])
#'''



#pathway enrichment
'''c = pd.read_csv('R_data/monocle3/branched_gene_cluster_all_raw.csv', sep = '\t', index_col = 0)
print(c)
#c = c.loc[c.iloc[:,0] != 4,:]
c = c - 1
from gprofiler import GProfiler
gp = GProfiler(return_dataframe=True)
goanno = {}
for i in range(6):
    tmp = c[c.iloc[:,0] == i]
    if 'GPX4' in tmp.index:
        print(i, 'GPX4')
    elif 'ACSL4' in tmp.index:
        print(i, 'ACSL4')
    elif 'AKT3' in tmp.index:
        print(i, 'AKT3')
    else:
        print(i)
    print(len(tmp))
    #result = gp.profile(organism = 'hsapiens',  query = tmp.index.tolist(), sources = ['GO:BP','GO:CC','GO:MF'], no_evidences = False)
    result = gp.profile(organism = 'hsapiens',  query = tmp.index.tolist(), sources = ['KEGG'], no_evidences = False)
    #result = gp.profile(organism = 'hsapiens',  query = tmp.index.tolist(), sources = ['GO:MF'], no_evidences = False)
    #result = result.loc[result.loc[:,'term_size']<300,:]
    result = result.loc[:,['name','intersections']]
    if True:
        if len(result.index) > 0:
            j = 0
            while j < len(result.index):
                print(result.iloc[j:(j+50),:])
                j = j + 50
            print()
#'''
