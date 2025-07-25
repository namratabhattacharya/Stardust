import matplotlib as plt
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import random
import sys
import os

def silhouette():
        if not os.path.exists("Stardust_results"):
            print("The directory structure Stardust_results doest not exist. Please run run_stardust first")
            sys.exit()
        if not os.path.exists("Stardust_results/analysis"):
                os.mkdir("Stardust_results/analysis")
        output_path="Stardust_results/analysis/"
        from sklearn.metrics import silhouette_samples, silhouette_score
        data_df = pd.read_csv('Stardust_results/visualization_output/3_pass/data.csv',delimiter=",",index_col=False)
        data_df.set_index('data',inplace=True)
        silhouette_avg = silhouette_score(data_df[['x','y']], data_df['cluster'])
        sample_silhouette_values = silhouette_samples(data_df[['x','y']], data_df['cluster'])
        print("silhouette score ",silhouette_avg)

        y_lower = 10
        import matplotlib.cm as cm
        fig=plt.figure(figsize=(4,7))
        n_clusters = len(list(data_df['cluster'].unique()))
        for i in range(n_clusters):
                # Aggregate the silhouette scores for samples belonging to
                # cluster i, and sort them
                ith_cluster_silhouette_values = \
                    sample_silhouette_values[data_df['cluster'] == i]

                ith_cluster_silhouette_values.sort()

                size_cluster_i = ith_cluster_silhouette_values.shape[0]
                y_upper = y_lower + size_cluster_i

                color = cm.nipy_spectral(float(i) / n_clusters)
                plt.fill_betweenx(np.arange(y_lower, y_upper),
                                  0, ith_cluster_silhouette_values,
                                  facecolor=color, edgecolor=color, alpha=0.7)

                # Label the silhouette plots with their cluster numbers at the middle
                plt.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

                # Compute the new y_lower for next plot
                y_lower = y_upper + 10  # 10 for the 0 samples

        plt.title("The silhouette plot for the various clusters.")
        plt.xlabel("silhouette coefficient",fontsize=20)
        plt.ylabel("Cluster label",fontsize=20)
        plt.axvline(x=silhouette_avg, color="red", linestyle="--")

        plt.yticks([])  # Clear the yaxis labels / ticks
        plt.xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])
        sns.despine(bottom = False, left = False)
        fig.savefig(output_path+"/silhouette.pdf",bbox_inches='tight',dpi=600)
        fig.savefig(output_path+"/silhouette.png",bbox_inches='tight',dpi=600)



def alluvial():
        if not os.path.exists("Stardust_results"):
            print("The directory structure Stardust_results doest not exist. Please run run_stardust first")
            sys.exit()
        if not os.path.exists("Stardust_results/analysis/"):
                os.mkdir("Stardust_results/analysis")
        output_path="Stardust_results/analysis/"
        from .alu import  plot
        import matplotlib.pyplot as plt
        import matplotlib.cm
        import matplotlib.colors as mcolors
        data_df = pd.read_csv('Stardust_results/visualization_output/3_pass/data.csv',delimiter=",",index_col=False)
        data_df.set_index('data',inplace=True)
        n_clusters = len(list(data_df['cluster'].unique()))
        louvain_cluster_df = pd.read_csv('Stardust_results/build_output/1_pass/louvain_cluster_df.csv')
        louvain_cluster_df.set_index('Unnamed: 0',inplace=True)
        louvain_cluster_df.index = louvain_cluster_df.index.map(str)
        meenet_cluster_df = data_df[data_df['type']=='cell']
        meenet_cluster_df = meenet_cluster_df[['cluster']]
        meenet_index = list(meenet_cluster_df.index)
        louvain_cluster_df = louvain_cluster_df.loc[meenet_index]
        louvain_cluster_df = louvain_cluster_df.sort_values(by=['louvain'])
        input_data=[]
        for ind in list(louvain_cluster_df.index):
            clust=[]
            clust.append("scanpy_"+str(louvain_cluster_df.loc[ind]['louvain']))
            clust.append("sd_"+str(meenet_cluster_df.loc[ind]['cluster']))
            input_data.append(clust)

        seed=7
        np.random.seed(seed)
        css_color = mcolors.CSS4_COLORS
        css_colors =[]
        for key in css_color.keys():
            css_colors.append(key)
        color = random.sample(css_colors,n_clusters)
        cmap = matplotlib.cm.get_cmap('jet')
        ax = plot(
            input_data,  alpha=0.4, color_side=1, rand_seed=seed, figsize=(7,7),
            disp_width=True, wdisp_sep=' '*2,colors=color,fontname='Monospace',
            labels=('Scanpy clusters','Stardust cluster'),label_shift=0)
        ax.set_title('Alluvial', fontsize=14, fontname='Monospace')
        plt.savefig(output_path+"alluvial_plot.png",bbox_inches='tight',dpi=600);
        plt.savefig(output_path+"alluvial_plot.pdf",bbox_inches='tight',dpi=600);


def heatmap():
        if not os.path.exists("Stardust_results"):
            print("The directory structure Stardust_results doest not exist. Please run run_stardust first")
            sys.exit()
        if not os.path.exists("Stardust_results/analysis/"):
                os.mkdir("Stardust_results/analysis")
        output_path="Stardust_results/analysis/"

        import matplotlib.colors as mcolors
        import anndata

        heat_map_df = pd.read_hdf('Stardust_results/build_output/1_pass/heat_map_df.h5', 'df')
        heat_map_df.fillna(0)
        if "Unnamed: 0" in list(heat_map_df.columns):
            heat_map_df.set_index(["Unnamed: 0"],inplace=True)
        data_df = pd.read_csv("Stardust_results/visualization_output/3_pass/data.csv",delimiter=",",index_col=False)
        data_df.set_index('data',inplace=True)
        data_df.drop('Unnamed: 0',axis=1,inplace=True)
        cells = list(data_df[data_df['type']=='cell'].index)
        genes = list(data_df[data_df['type']=='gene'].index)
        cluster = []
        gene_order = []
        cell_order = []
        heat_map_df = heat_map_df.loc[cells]
        heat_map_df = heat_map_df[genes]
        data_df = data_df.sort_values(by=['cluster'])
        index_list = data_df.index
        gene_cluster =[]
        for ind in index_list:
            if data_df.loc[ind]['type'] == 'cell':
                cell_order.append(ind)
                cluster.append(data_df.loc[ind]['cluster'])
            else:
                gene_order.append(ind)
        heat_map_df = heat_map_df.reindex(cell_order)
        heat_map_df = heat_map_df.reindex(columns = gene_order)
        cluster = [str(s) for s in cluster]
        heat_map_df['cluster'] = cluster
        clustering = heat_map_df['cluster']
        clustering  = clustering.astype('category')
        heat_map_df.drop('cluster',axis=1,inplace=True)
        heat_map_df = heat_map_df[heat_map_df.columns].astype(float)

        obs = pd.DataFrame(index=heat_map_df.index)
        var = pd.DataFrame(index=heat_map_df.columns)
        X = heat_map_df.values
        adata1 = anndata.AnnData(X=X,obs=obs,var=var)
        adata1.obs['cluster'] = clustering
        import scanpy as sc
        sc.pp.log1p(adata1)
        sc.tl.rank_genes_groups(adata1, 'cluster', method='t-test_overestim_var',key_added="t-test")
        sc.pl.rank_genes_groups_heatmap(adata1, n_genes=3, key="t-test", groupby="cluster", show_gene_labels=True,figsize=(10,5),save='.pdf')


