#!/usr/bin/env python
# coding: utf-8

# In[1]:



from node2vec import Node2Vec
import networkx as nx
import matplotlib as plt
import csv
import pandas as pd
from sklearn.neighbors import kneighbors_graph
import numpy as np
import scanpy as sc
from annoy import AnnoyIndex
from sklearn.manifold import TSNE
import seaborn as sns
import umap.umap_ as umap
from sklearn import preprocessing
from sklearn.decomposition import PCA
import subprocess
import os
import shutil
import sys
import math
from sklearn.metrics.cluster import adjusted_rand_score
# # DEFAULT PARAMETER SETTING


sc.settings.verbosity = 2
np.set_printoptions(precision=2)
reducer = umap.UMAP(random_state=42)
G_normalized = nx.Graph()
adata = None;

# # READ 10X DATA


def read_data(inp_path,out_path,type):
    import anndata
    import math
    if type=="csv":
        #data_df = pd.read_csv(inp_path+"GSE115469_Data.csv") #liver
        data_df = pd.read_csv(inp_path+"expression.csv")   #melanoma
        if 'Unnamed: 0' in list(data_df.columns):
            data_df.set_index('Unnamed: 0',inplace=True)
        global adata;
        obs=pd.DataFrame(index=data_df.index)
        var=pd.DataFrame(index=data_df.columns)
        X=data_df.values
        adata=anndata.AnnData(X=X,obs=obs,var=var)
    elif type=="10x":
        adata = sc.read_10x_mtx(
        inp_path,  # the directory with the `.mtx` file
        var_names='gene_symbols',                  # use gene symbols for the variable names (variables-axis index)
        cache=True)


# # FILTER CELL AND GENES
def compute(pca_n,cN,gN,out_path,n_pass):
    global adata;
    import math
    sc.pp.filter_cells(adata, min_genes=0)
    val = math.ceil(0.3*adata.obs['n_genes'].mean())
    sc.pp.filter_cells(adata, min_genes=val)
    sc.pp.filter_genes(adata,min_cells=3)
    data_df = adata.to_df()
    orig_data_df = data_df
    if n_pass==1:
        #data_df.to_csv(out_path+'/heat_map_df.csv')
        data_df.to_hdf(out_path+"heat_map_df.h5",key='df')

    import pickle
    with open(output_path+'/gene_list.txt', 'wb') as fp:
            pickle.dump(list(data_df.columns), fp)
    with open(output_path+'/cell_list.txt', 'wb') as fp:
            pickle.dump(list(data_df.index), fp)

    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata,n_top_genes=1500)

    variable_genes=adata.var['highly_variable']
    variable_gene_list = []
    for ind in variable_genes.index:
        if variable_genes[ind]==True:
            variable_gene_list.append(ind)
    with open(output_path+'/variable_genes.txt', 'wb') as fp:
            pickle.dump(variable_gene_list, fp)

    if pca_n !=0:
        adata = adata[:, adata.var.highly_variable]

# # NORMALIZATION
    sc.pp.normalize_total(adata, target_sum=None, inplace=True)
    data_df = adata.to_df()
    read_data_df = adata.to_df()
    column_names = data_df.columns


# # LOUVAIN CLUSTERING
    print("Starting louvain clustering....")
    if pca_n !=0:
        sc.pp.neighbors(adata,n_neighbors=25)
        sc.tl.umap(adata)
        import pickle
        umap_coord = adata.obsm['X_umap']
        with open(output_path+'/umap_coord.txt', 'wb') as fp:
            pickle.dump(umap_coord, fp)

        sc.tl.louvain(adata)
        p=np.array(adata.uns['neighbors']['distances'].todense())
        neighbor_df = pd.DataFrame(p,columns=list(read_data_df.index),index=list(read_data_df.index))
        louvain_cluster_df=sc.get.obs_df(adata,keys=["louvain"])
        louvain_cluster_df.to_csv(out_path+'/louvain_cluster_df.csv')
        n_cluster=louvain_cluster_df['louvain'].nunique()
        size_per_cluster = louvain_cluster_df['louvain'].value_counts()


# # SUBSAMPLING
    print("Starting subsampling....")
    if pca_n != 0:
        pl=0.1
        pu=0.9
        k=500
        prob_per_cluster=[]
        for i in range(0,n_cluster):
            prob_per_cluster.append(pl-np.exp(-(size_per_cluster[i]/k))*(pl-pu))
        sample_per_cluster=[]
        total_sample = []
        for i in range(0,n_cluster):
            sample_per_cluster.append(int(prob_per_cluster[i]*size_per_cluster[i]))
        sampled_df = pd.DataFrame(columns=data_df.columns)
        left_out_cells=[]
        for i in range(0,n_cluster):
            cells=list(louvain_cluster_df.loc[louvain_cluster_df['louvain'] == str(i)].index)
            sample = np.random.choice(cells,sample_per_cluster[i],replace=False)
            remainder = list(set(cells)-set(sample))
            left_out_cells=left_out_cells+remainder
            for smp in sample:
                sampled_df=sampled_df.append(data_df.loc[smp])
                total_sample.append(smp)
        import pickle
        with open(output_path+'/sampled_cells.txt', 'wb') as fp:
            pickle.dump(total_sample, fp)

    else:
        import pickle
        with open ('Stardust_results/build_output/1_pass/sampled_cells.txt', 'rb') as fp:
            sampled_cells = pickle.load(fp)
        sampled_df = data_df.loc[sampled_cells]


    data_df=sampled_df
    column_names = data_df.columns

    ## TOP PCA componnets############
    if pca_n!=0:
         import anndata
         obs = pd.DataFrame(index=data_df.index)
         var = pd.DataFrame(index=data_df.columns)
         X = data_df.values
         adata = anndata.AnnData(X=X,obs=obs,var=var)
         sc.tl.pca(adata, n_comps=500,svd_solver='arpack')
         d2 = adata.to_df()
         orig_pca_data = pd.DataFrame(index=list(d2.index),data=adata.obsm['X_pca'])
         #orig_pca_data.to_csv(out_path+'/pca_df.csv')
         orig_pca_data.to_hdf(out_path+"pca_df.h5",key='df')
    else:
         #orig_pca_data = pd.read_csv("Stardust_results/build_output/1_pass/pca_df.csv")
         #orig_pca_data.set_index("Unnamed: 0",inplace=True)
         orig_pca_data = pd.read_hdf('Stardust_results/build_output/1_pass/pca_df.h5', 'df')
         if "Unnamed: 0" in list(orig_pca_data.columns):
            orig_pca_data.set_index("Unnamed: 0",inplace=True)


# # TOP PCA GENES
    print("PCA gene selection....")
    if pca_n!=0:
    	pca = PCA(n_components=pca_n+200, svd_solver='full')
    	pca.fit(data_df)
    	rotation_mat = np.transpose(abs(pca.components_))
    	rotation_matrix_df = pd.DataFrame(rotation_mat,index=column_names)
    	gene_max = rotation_matrix_df.max(axis=1)
    	gene_max=gene_max.sort_values(ascending=False)
    	pca_genes = list(gene_max[:pca_n].index)
    else:
        import pickle
        if n_pass ==2:
            filename='Stardust_results/visualization_output/1_pass/de_genes.txt'
        if n_pass ==3:
            filename='Stardust_results/visualization_output/2_pass/de_genes.txt'
        if n_pass ==4:
            filename='Stardust_results/visualization_output/3_pass/de_genes.txt'
        with open (filename, 'rb') as fp:
            pca_genes = pickle.load(fp)
        orig_data_df = orig_data_df[pca_genes]

    pca_genes_df = data_df[pca_genes]
    cell_names = list(pca_genes_df.index.values)
    column = list(pca_genes_df.columns)
    pca_genes_df.to_csv(out_path+'/pca_genes_df.csv')


# # CELL CELL GRAPH BY ANNOY
    print("Network building ....")
    if n_pass ==1:
        f = len(orig_pca_data.columns)
        t = AnnoyIndex(f, 'angular')  # Length of item vector that will be indexed
        for i in range(len(orig_pca_data.index)):
            v  = orig_pca_data.iloc[i]
            t.add_item(i, v)

        t.build(30)
        print("graph build....")
        l=[]
        upper_bound = cN-1
        for i in range(upper_bound):
            l.append(i)
        cell_cell_df = pd.DataFrame(index=l)
        for i in range(len(orig_pca_data.index)):
               knn_values = list(t.get_nns_by_item(i,cN,include_distances=True))
               cell_cell_df["item."+str(i)]=knn_values[0][1:]
               cell_cell_df["distance."+str(i)] = knn_values[1][1:]
        cell_cell_df.to_hdf(out_path+"cell_cell_df.h5",key='df')
    else:
        upper_bound = cN-1
        cell_cell_df = pd.read_hdf('Stardust_results/build_output/1_pass/cell_cell_df.h5', 'df')
    cell_graph_column = list(cell_cell_df.columns)
    print("NN done....")


# # CELL-GENE CONNECTIONS
    print("cell gene connection building....")
    c1=0
    c2=0
    if n_pass ==4:
        con_cells = []
        connections = 5
        cell_gene_edge_list = []
        import pickle
        with open ('Stardust_results/visualization_output/3_pass/de_genes_cluster.txt', 'rb') as fp:
            de_gene_cluster = pickle.load(fp)
        data_df = pd.read_csv("Stardust_results/visualization_output/3_pass/data.csv",delimiter=",",index_col=False)
        data_df.set_index('data',inplace=True)
        cell_data_df = data_df[data_df['type']=='cell']
        for rank in range(len(de_gene_cluster[:20])):
            connections = 5-math.ceil(rank/10)*1
            n_genes = de_gene_cluster[rank]
            for gene in n_genes:
                if gene == "NULL":
                    continue
                c1=c1+1
                cluster = n_genes.index(gene)
                n_cells = list(cell_data_df[cell_data_df['cluster']==cluster].index)
                cell_type_df = pd.DataFrame(columns=["expression_value"],index=n_cells,dtype=float)
                for cell in n_cells:
                    cell_type_df.loc[cell]["expression_value"] = read_data_df.loc[cell][gene]
                n_cells = cell_type_df.nlargest(connections,"expression_value")["expression_value"].index.values
                for cell in n_cells:
                    c2=c2+1
                    temp=[]
                    temp.append(gene)
                    temp.append(cell)
                    temp.append(float(1))
                    cell_gene_edge_list.append(tuple(temp))

        G_normalized.add_weighted_edges_from(cell_gene_edge_list)


    else:
        cell_gene_edge_list = []
        for col in column:
    	    n_cells = pca_genes_df.nlargest(gN,col)[col].index.values

    	    for cell in n_cells:
        	    temp=[]
        	    temp.append(col)
        	    #temp.append(cell_names.index(cell)+1)
        	    temp.append(cell)
        	    temp.append(float(1))
        	    cell_gene_edge_list.append(tuple(temp))
        	    i=5
        	    while(i>=1):
            		    distance="distance."+str(cell_names.index(cell))
            		    if cell_cell_df.iloc[i][distance]!=999:
                		    cell_cell_df.at[i,distance]=999
                		    break
            		    else:
                		    i=i-1
        G_normalized.add_weighted_edges_from(cell_gene_edge_list)

# # CELL -CELL WITHOUT DROPPING CELL CONNECTIONS
    print("cell cell connection building ....")
    cell_cell_edge_list=[]
    c=0
    g=0
    for col in cell_graph_column:
        if "distance" in col or "Unnamed" in col:
            continue
        else:
            if col=="item":
                distance="distance"
            else:
                distance="distance."+str(col[5:])
            for i in range(1,upper_bound):
                    temp_edge=[]
                    #temp = list(cell_cell_df["distance."+str(int(cell_cell_df.iloc[i][col]))])
                    temp_edge.append(pca_genes_df.index[int(col[5:])])
                    temp_edge.append(pca_genes_df.index[int(cell_cell_df.iloc[i][col])])
                    temp_edge.append(float(1))
                    cell_cell_edge_list.append(tuple(temp_edge))
                    c=c+1
    G_normalized.add_weighted_edges_from(cell_cell_edge_list)


# # WRITE GRAPH

    nx.write_gexf(G_normalized, out_path+"/graph_normalized_own_cc.gexf")
    nx.write_gexf(G_normalized, "graph_normalized_own_cc.gexf")
    if pca_n !=0:
        neighbor_df=neighbor_df.drop(list(sampled_df.index),axis=0)
        neighbor_df=neighbor_df.drop(left_out_cells,axis=1)
        #neighbor_df.to_csv(out_path+'/neighbor.csv')
        neighbor_df.to_hdf(out_path+"neighbor.h5",key='df')

    subprocess.call(['java', '-jar', 'stardust/run_stardust/openOrd.jar',out_path+'graph_normalized_own_cc.gexf',out_path+'output_normalized_own_cc.csv'])


input_path=""
output_path=""
file=""
n_pass=1
#def build_start(data_path,out_path,file_type,N_pass,pca_n=300,cn=20,gn=10):
#    input_path = data_path
#    output_path = out_path
#    type = file_type
#    pca_n = pca_n
#    cN=cn
#    gN=gn
#    n_pass=N_pass
#    read_data(input_path,output_path,type);
#    compute(pca_n,cN,gN,output_path,n_pass);


input_path=""
output_path=""
file=""
n_pass=1
arguments=sys.argv[1:]
if len(arguments)>=8 and arguments[0]=="-i" and arguments[2]=="-t" and arguments[4]=="-o" and arguments[6]=="-n_pass" :
    input_path=arguments[1]
    type = arguments[3]
    output_path=arguments[5]
    temp = arguments[8:]
    n_pass =int( arguments[7])
    pca_n = 300
    cN = 20   #20
    gN = 10   #10

    if len(temp)>0:
        if temp[0] == '-pca_n':
            pca_n = int(temp[1])
        if temp[0] == '-cN':
            cN = int(temp[1])
        elif len(temp)>2 and temp[2] == '-cN':
            cN = int(temp[3])
        if temp[0] == '-gN':
            gN = int(temp[1])
        elif len(temp)>2 and temp[2] == '-gN':
            gN = int(temp[3])
        elif len(temp)>4 and temp[4] == '-gN':
            gN = int(temp[5])
    read_data(input_path,output_path,type);
    compute(pca_n,cN,gN,output_path,n_pass);
else:
        print ("\nPlease follow the command: <file>.py -i inputpath -o outputpath -pca_n integer -cN integer -gN integer\n")
        print('inputpath : path for read 10x data\n')
        print('outputpath : path where the output is to be stored\n')
        print('pca_n : number of pca components\n');
        print('cN : number of connections per cell\n');
        print('gN : number of connections per gene\n');
        sys.exit()







