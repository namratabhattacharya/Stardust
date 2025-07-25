#!/usr/bin/env python
# coding: utf-8

# In[17]:

import json
import matplotlib as plt
import csv
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from decimal import Decimal
import seaborn as sns
import pandas as pd
import networkx as nx
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans
import operator
import numpy as np
import random
import sys
from sklearn.metrics.cluster import adjusted_rand_score
import pickle
# In[18]:


import seaborn
seaborn_colors=[]
for k in seaborn.xkcd_rgb.keys():
    seaborn_colors.append(k)


def compute(inp_dataset,input_path,output_path,de_analysis,n_pass):

    print("Current pass ",n_pass)
    import json
    import matplotlib as plt
    import csv
    from sklearn.manifold import TSNE
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from decimal import Decimal
    import seaborn as sns
    import pandas as pd
    import networkx as nx
    from sklearn.cluster import DBSCAN
    from sklearn.cluster import KMeans
    import operator
    import numpy as np
    import random
    import sys


    #csvData=[['data','x','y','type']]
    print("Processing the input data into datafames....")
    csvData=[]
    count=0
    #filename = "G:/Thesis/Dropclust/plots/output_normalized_own_cc.csv" filename = "G:/Thesis/Dropclust/plots/PCA_GENES/output_normalized_own_cc.csv" filename =
    #"G:/Thesis/Dropclust/output_normalized_zscore_cc1.csv" filename = "C:/Users/Swagatam/IdeaProjects/openOrd/output_normalized_own_cc.csv"
    filename = input_path+"/output_normalized_own_cc.csv"
    coord_data = pd.read_csv(filename,names=['data','x','y'])
    coord_data.set_index('data',inplace=True)
    data=[]
    data_outlier=[]
    with open(filename, 'r') as csvfile:
         csvreader = csv.reader(csvfile)
         for row in csvreader:
            #f=0
             #row=[float(i) for i in row]
            data.append(row)
            temp_outlier=[]
            temp_outlier.append(row[1])
            temp_outlier.append(row[2])
            data_outlier.append(temp_outlier)
            temp=row
            #if row[0].isnumeric():
            #    temp.append('cell')
            if len(row[0]) >=16:
                temp.append('cell')
            else:
                temp.append('gene')
                count=count+1
            csvData.append(temp)


    # # DB SCAN

    # In[20]:

    if n_pass != 4:
        noise =[]
        print("Performing clustering....")
        db = DBSCAN(eps=180,min_samples=55).fit_predict(data_outlier)
        final_data=[]
        csvData=[['data','x','y','type']]
        for i in range(0,len(list(db))):
            if db[i]!=-1:
                final_data.append(data[i])
                csvData.append(data[i])
            if db[i]==-1:
                noise.append(data[i][0])
        data=final_data

        n_clusters = len(set(db)) - (1 if -1 in list(db) else 0)
        print("Clustering done. the number of obtained clusters: ",n_clusters)
    else:
        remove_data = []

        prev_df = pd.read_csv("Stardust_results/visualization_output/3_pass/data.csv",delimiter=",",index_col=False)
        prev_df.set_index('data',inplace=True)
        clusters_info = []
        for i in range(0,len(csvData)):
            if csvData[i][3] == 'cell':
                if csvData[i][0] in (prev_df.index):
                    clusters_info.append(prev_df.loc[csvData[i][0]]['cluster'])
                else:
                    remove_data.append(csvData[i])
            else:
                f=0
                import pickle
                with open ('Stardust_results/visualization_output/3_pass/de_genes_cluster.txt', 'rb') as fp:
                    de_gene_cluster = pickle.load(fp)
                for rank in range(0,len(de_gene_cluster)):
                    if csvData[i][0] in de_gene_cluster[rank]:
                        f=1
                        clusters_info.append(de_gene_cluster[rank].index(csvData[i][0]))
                        break
                if f==0:
                    remove_data.append(csvData[i])
        for r in remove_data:
            csvData.remove(r)
        temp=[['data','x','y','type']]
        temp.extend(csvData)
        csvData = temp



    # In[13]:



    # # OUTLIER VISUALIZATION

    # In[21]:
    if n_pass !=4:
        print("Starting outlier detection....")
        data_type = []
        c=0
        g=0
        for i in range(0,len(coord_data)):
            if db[i]!=-1:
                data_type.append("data")
            else:
                if len(coord_data.index[i])>=16:
                    data_type.append("cell_outliers")
                else:
                    g=g+1
                    data_type.append("gene_outliers")
        coord_data["data_type"] = data_type
        data_colors = ["lightblue"]
        if g>0:
            noise_colors = ['blue','red']
        else:
            noise_colors = ['blue']
        coord_data["alpha"] = np.where(coord_data['data_type'] == 'data', 0.5,1.0)
        plt.figure(figsize=(6,4.5))
        #ax = sns.scatterplot(x="x", y="y", data=coord_data[coord_data['alpha']==0.5],hue="data_type",palette=sns.xkcd_palette(data_colors),sizes=(50,100),size="data_type",alpha=0.3)
        #sns.scatterplot(x="x", y="y", data=coord_data[coord_data['alpha']==1.0],hue="data_type",palette=sns.xkcd_palette(noise_colors),sizes=(50,100),size="data_type",marker="^",alpha=1.0,ax=ax)
        marker = {"gene_outliers": "^", "cell_outliers": "^"}
        ax = sns.scatterplot(x="x", y="y", data=coord_data[coord_data['alpha']==0.5],hue="data_type",palette=sns.xkcd_palette(data_colors),sizes=(50,100),size="data_type",linewidth=0.0,s=10,alpha=0.3)
        sns.scatterplot(x="x", y="y", data=coord_data[coord_data['alpha']==1.0],hue="data_type",palette=sns.xkcd_palette(noise_colors),sizes=(100,50),size="data_type",style="data_type",markers=marker,alpha=1.0,linewidth=0.0,s=10,legend='brief',ax=ax)
        #plt.legend(title=='')
        ax.legend(bbox_to_anchor=(1.1, 1.05),frameon=False)
        sns.despine(bottom = False, left = False)
        plt.xlabel("dim1")
        plt.ylabel("dim2")
        plt.savefig(output_path+'outliers_visualization.png',bbox_inches='tight')
        print("Outliers removed from the dataset....")

    # # POST-HOC CLUSTER ASSIGNMENT

    # In[23]:

    print("Starting post hoc clustering....")
    neighbor_df = pd.read_hdf('Stardust_results/build_output/1_pass/neighbor.h5', 'df')
    if 'Unnamed: 0' in list(neighbor_df.columns):
        neighbor_df.set_index('Unnamed: 0',inplace=True)
    p=0
    col = list(neighbor_df.columns)
    index = list(neighbor_df.index)
    cell_dict = dict()
    column_dict = dict()
    for i in range(len(col)):
        column_dict[i]=col[i]
    for i in range(len(list(neighbor_df.index))):
        row = neighbor_df.iloc[i]
        col_ind = list(row.to_numpy().nonzero())[0]
        for ind in col_ind:
            if index[i] in cell_dict.keys():
                cell_dict[index[i]].append(column_dict[ind])
            else:
                temp=[]
                temp.append(column_dict[ind])
                cell_dict[index[i]]=temp
    cluster_assign = []
    for key_cell in cell_dict.keys():
        clust = dict()
        cells = cell_dict[key_cell]
        for cell in cells:
                if n_pass == 4 :
                    if cell in list(prev_df.index):
                        cluster = prev_df.loc[cell]['cluster']
                    else:
                        cluster = -1
                else:
                    cluster = db[list(coord_data.index).index(cell)]
                if cluster not in clust.keys():
                    clust[cluster]=1
                else:
                    clust[cluster]=clust[cluster]+1
        max_cluster=max(clust.items(), key=operator.itemgetter(1))[0]
        if max_cluster == -1:
            continue
        cluster_assign.append(max_cluster)
        x_total = 0
        y_total = 0
        count = 0
        for cell in cells:
            if (n_pass!=4 and db[list(coord_data.index).index(cell)] == max_cluster) or (n_pass==4 and cell in list(prev_df.index) and prev_df.loc[cell]['cluster'] == max_cluster):
                count=count+1
                x_total=x_total+coord_data.loc[cell]['x']
                y_total=y_total+coord_data.loc[cell]['y']
        temp=[]
        temp.append(key_cell)
        temp.append(x_total/count)
        temp.append(y_total/count)
        temp.append('cell')
        p=p+1
        csvData.append(temp)
    print("Post hoc clustering done....")

    # In[24]:



    with open(output_path+'data.csv', 'w') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(csvData)
    csvFile.close()
    data_df = pd.read_csv(output_path+"data.csv",delimiter=",",index_col=False)
    if n_pass !=4 :
        clusters_info = [x for x in db if x!=-1]
        clusters_info = clusters_info + cluster_assign
    else:
        clusters_info = clusters_info + cluster_assign
        data_df['cluster'] = clusters_info
    data_df.to_csv(output_path+'data.csv')
    n_clusters = len(list(set(clusters_info)))
    print("cluster saved ....")






    n_clusters = len(data_df['cluster'].unique())
    colors = random.sample(seaborn_colors,n_clusters)


    colors = random.sample(seaborn_colors,n_clusters)
    plt.figure(figsize=(5,5))
    #cmap = sns.cubehelix_palette(dark=.3, light=.8, as_cmap=True)
    ax = sns.scatterplot(x="x", y="y", data=data_df,hue="cluster",palette=sns.xkcd_palette(colors),linewidth=0.0,s=2)
    ax.legend(bbox_to_anchor=(1.0, 1.00),frameon=False)
    for cl in range(n_clusters):
        plt.annotate(cl,data_df.loc[data_df['cluster']==cl,['x','y']].mean(),
                horizontalalignment='center',
                verticalalignment='center',
                size=10, weight='bold',
                color="black")
    sns.despine(bottom = False, left = False)
    plt.xlabel("sd1",fontsize=20)
    plt.ylabel("sd2",fontsize=20)
    plt.setp(ax.spines.values(), linewidth=2)
    plt.yticks([],linewidth=20)
    plt.xticks([])
    plt.savefig(output_path+"cluster_visualization.png",bbox_inches='tight',dpi=600);
    plt.savefig(output_path+"cluster_visualization.pdf",bbox_inches='tight',dpi=600);

    if n_pass ==3:
        from sklearn.datasets import make_blobs
        from sklearn.metrics import silhouette_samples, silhouette_score
        silhouette_avg = silhouette_score(data_df[['x','y']], data_df['cluster'])
        sample_silhouette_values = silhouette_samples(data_df[['x','y']], data_df['cluster'])
        print(silhouette_avg)

        y_lower = 10
        import matplotlib.cm as cm
        #fig, (ax1, ax2) = plt.subplots(1, 2)
        fig=plt.figure(figsize=(4,7))
        #fig.set_size_inches(18, 7)
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

    #  #  MARKER FINDING
    data_df = pd.read_csv(output_path+"data.csv",delimiter=",",index_col=False)
    data_df.set_index('data',inplace=True)
    import pickle
    if n_pass == 2:
        path = 'Stardust_results/visualization_output/1_pass'
    if n_pass == 3:
        path = 'Stardust_results/visualization_output/2_pass'
    if n_pass == 4:
        path = 'Stardust_results/visualization_output/3_pass'
    if n_pass != 1:
        with open (path+'/de_genes_cluster.txt', 'rb') as fp:
                de_gene_cluster = pickle.load(fp)

        marker =[]
        disp_marker = []
        for cl in range(n_clusters):
            cls = data_df[data_df['cluster']==cl]
            gene_df = cls[cls['type']=='gene']
            f=0
            for rank in range(len(de_gene_cluster)):
                if f==1:
                    break
                for gene in de_gene_cluster[rank]:
                    if gene in list(gene_df.index):
                        disp_marker.append(gene)
                    #print(cl)
                        f=1
                        break
        marker = disp_marker

        #sys.exit(0)

    # # CELL GENE MARKER

    # In[28]:
    from sklearn.neighbors import KNeighborsRegressor
    prev_pass_data =  pd.read_csv('Stardust_results/visualization_output/3_pass/data_openOrd.csv')
    prev_pass_data.set_index('data',inplace=True)
    data_df = pd.read_csv(output_path+'/data.csv')
    data_df.set_index('data',inplace=True)
    gene_df = data_df[data_df['type']=='gene']
    x_gene_fit = list(gene_df['x'])
    y_gene_fit = list(gene_df['y'])
    cells = list(prev_pass_data.index)
    cell_list=[]
    x_coord = []
    y_coord = []

    for i in range(len(cells)):
             if cells[i] in list(data_df.index):
                 cell_list.append(cells[i])
                 x_coord.append(prev_pass_data.iloc[i]['x'])
                 y_coord.append(prev_pass_data.iloc[i]['y'])

    prev_df = pd.DataFrame(index=cell_list)
    prev_df['x'] = x_coord
    prev_df['y'] = y_coord

    import numpy as np
    from sklearn.linear_model import Lasso
    from sklearn.neighbors import KNeighborsRegressor
    import pickle
    cells=[]
    genes=[]
    gene_coord_x = []
    gene_coord_y = []

    for i in range(n_clusters):
             clust_data = data_df[data_df['cluster']==i]
             clust_cells = clust_data[clust_data['type']=='cell']
             clust_genes = clust_data[clust_data['type']=='gene']
             cells.extend(list(clust_cells.index))
             genes.extend(list(clust_genes.index))
             if len(list(clust_genes.index))==0:
                 continue
             model1 = KNeighborsRegressor(n_neighbors=4)

             model2 = KNeighborsRegressor(n_neighbors=4)
             temp =[]
             for cell in list(clust_cells.index):
                 if cell in list(prev_df.index):
                     temp.append(cell)
             clust_cells = clust_cells.loc[temp]
             model1.fit(np.array(list(clust_cells['x'])).reshape((-1, 1)),np.array(list(prev_df.loc[list(clust_cells.index)]['x'])).reshape((-1, 1)))

             filename = output_path+'/sd_x_KNN_model.sav'
             pickle.dump(model1, open(filename, 'wb'))
             #model1 = pickle.load(open(filename, 'rb'))
             x_gene_pred = model1.predict(np.array(list(clust_genes['x'])).reshape((-1, 1)))
            #  gene_coord_x.extend(x_gene_pred)  #old
             gene_coord_x.extend(np.ravel(x_gene_pred).tolist())
             model2.fit(np.array(list(clust_cells['y'])).reshape((-1, 1)),np.array(list(prev_df.loc[list(clust_cells.index)]['y'])).reshape((-1, 1)))

             filename = output_path+'/sd_y_KNN_model.sav'
             pickle.dump(model2, open(filename, 'wb'))
             #model2 = pickle.load(open(filename, 'rb'))
             y_gene_pred = model2.predict(np.array(list(clust_genes['y'])).reshape((-1, 1)))
            #  gene_coord_y.extend(y_gene_pred) #old
             gene_coord_y.extend(np.ravel(y_gene_pred).tolist())

    with open(output_path+"/sd_gene_coord_x.txt", 'wb') as fp:
            pickle.dump(gene_coord_x,fp)
    with open(output_path+"/sd_gene_coord_y.txt", 'wb') as fp:
            pickle.dump(gene_coord_y,fp)

    #with open (output_path+"/sd_gene_coord_x.txt", 'rb') as fp:
    #        gene_coord_x = pickle.load(fp)
    #with open (output_path+"/sd_gene_coord_y.txt", 'rb') as fp:
    #        gene_coord_y = pickle.load(fp)

    import matplotlib.pyplot as plt, mpld3
    from scipy.spatial import ConvexHull, convex_hull_plot_2d
    prev_pass_data =  pd.read_csv('Stardust_results/visualization_output/3_pass/data_openOrd.csv')
    prev_pass_data["alpha"] = np.where(prev_pass_data['type'] == 'gene', 1.0, 0.5)
    color_gene = ["light blue"]
    color_cell = ["red"]
    #fig,ax1 = plt.subplots()
    plt.figure(figsize=(6,6))
    ax = sns.scatterplot(x="x", y="y", data=prev_pass_data[prev_pass_data['alpha']==0.5],hue="type",palette=sns.xkcd_palette(color_gene),sizes=(10,5),size="type",alpha=0.3,s=10)
    #sns.scatterplot(x="x", y="y", data=data_df[data_df['alpha']==1.0],hue="type",palette=sns.xkcd_palette(color_cell),sizes=(20,5),size="type",marker="^",alpha=1.0,ax=ax,s=10)
    sns.scatterplot(x=gene_coord_x, y=gene_coord_y,palette=sns.xkcd_palette(color_cell),sizes=(20,5),marker="^",alpha=1.0,ax=ax,s=10)
    for c in range(n_clusters):
        p=data_df[data_df["cluster"]==c]
        p=p[['x','y']]
        points = p.values
        hull = ConvexHull(points)
        #for simplex in hull.simplices:
       #    sns.lineplot(points[simplex, 0], points[simplex, 1])
    x_list = []
    y_list = []
    if n_pass !=1:
        for m in marker:
            #x_list.append(data_df.loc[m]['x'])
            x_list.append(gene_coord_x[genes.index(m)])
            #y_list.append(data_df.loc[m]['y'])
            y_list.append(gene_coord_y[genes.index(m)])
        for label, x, y in zip(marker, x_list, y_list):
            plt.annotate(
                 label,
                 xy=(x, y), xytext=(-20, 20),
                 textcoords='offset points', ha='right', va='bottom',
                 #bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                 arrowprops=dict(arrowstyle = '-', connectionstyle='arc3,rad=0'))
    ax.legend(bbox_to_anchor=(1.0, 1.00),frameon=False)
    sns.despine(bottom = False, left = False)
    plt.xlabel("sd1",fontsize=20)
    plt.ylabel("sd2",fontsize=20)
    plt.setp(ax.spines.values(), linewidth=2)
    plt.yticks([],linewidth=20)
    plt.xticks([])
    plt.savefig(output_path+"sd_embedding.png",bbox_inches='tight',dpi=600);
    plt.savefig(output_path+"sd_embedding.pdf",bbox_inches='tight',dpi=600);


    import matplotlib.pyplot as plt, mpld3
    from scipy.spatial import ConvexHull, convex_hull_plot_2d
    #data_df["alpha"] = np.where(data_df['type'] == 'gene', 1.0, 0.5)
    prev_pass_data.set_index('data',inplace=True)
    temp_data = prev_pass_data[prev_pass_data['type']=='cell']
    temp_genes = data_df[data_df['type']=='gene']
    for pos in range(0,len(genes)):
        temp_genes.at[genes[pos],'x'] = gene_coord_x[pos]
        temp_genes.at[genes[pos],'y'] = gene_coord_y[pos]
    temp_data.append(temp_genes)
    color_gene = ["light blue"]
    color_cell = ["red"]
    n_clusters = len(data_df['cluster'].unique())
    colors = random.sample(seaborn_colors,n_clusters)
    #fig,ax1 = plt.subplots()
    plt.figure(figsize=(6,6))
    ax = sns.scatterplot(x="x", y="y", data=temp_data,hue="cluster",palette=sns.xkcd_palette(colors),s=2,linewidth=0.0)
    #sns.scatterplot(x="x", y="y", data=data_df[data_df['alpha']==1.0],hue="type",palette=sns.xkcd_palette(color_cell),sizes=(20,5),size="type",marker="^",alpha=1.0,ax=ax,s=10)
    #sns.scatterplot(x=gene_coord_x, y=gene_coord_y,palette=sns.xkcd_palette(color_cell),sizes=(20,5),marker="^",alpha=1.0,ax=ax,s=20)
    for c in range(n_clusters):
        p=data_df[data_df["cluster"]==c]
        p=p[['x','y']]
        points = p.values
        hull = ConvexHull(points)
        #for simplex in hull.simplices:
       #    sns.lineplot(points[simplex, 0], points[simplex, 1])
    x_list = []
    y_list = []
    d1 = prev_pass_data[prev_pass_data['alpha']==0.5]
    for cl in range(n_clusters):
            plt.annotate(cl,d1.loc[d1['cluster']==cl,['x','y']].mean(),
                        horizontalalignment='center',
                        verticalalignment='center',
                        size=10, weight='bold',
                        color="black")
    if n_pass !=1:
        for m in marker:
            #x_list.append(data_df.loc[m]['x'])
            x_list.append(gene_coord_x[genes.index(m)])
            #y_list.append(data_df.loc[m]['y'])
            y_list.append(gene_coord_y[genes.index(m)])
        for label, x, y in zip(marker, x_list, y_list):
            plt.annotate(
                 label,
                 xy=(x, y), xytext=(-20, 20),
                 textcoords='offset points', ha='right', va='bottom',
                 #bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                 arrowprops=dict(arrowstyle = '-', connectionstyle='arc3,rad=0'))
    ax.legend(bbox_to_anchor=(1.0, 1.00),frameon=False)
    sns.despine(bottom = False, left = False)
    plt.xlabel("sd1",fontsize=20)
    plt.ylabel("sd2",fontsize=20)
    plt.setp(ax.spines.values(), linewidth=2)
    plt.yticks([],linewidth=20)
    plt.xticks([])
    plt.savefig(output_path+"sd_color_embedding.png",bbox_inches='tight',dpi=600);
    plt.savefig(output_path+"sd_color_embedding.pdf",bbox_inches='tight',dpi=600);
    #sys.exit(0)
# # UMAP CELL GENE MARKER # #

    if n_pass ==4:

        import pickle
        with open ('Stardust_results/build_output/1_pass/umap_coord.txt', 'rb') as fp:
            umap_coord = pickle.load(fp)
        louvain_df = pd.read_csv('Stardust_results/build_output/1_pass/louvain_cluster_df.csv')
        louvain_df.set_index('Unnamed: 0',inplace=True)
        #data_df = pd.read_csv('F:/output/output_visualize_melanoma_pca/3rd_pass/data.csv')
        data_df = pd.read_csv(output_path+'/data.csv')
        data_df.set_index('data',inplace=True)
        gene_df = data_df[data_df['type']=='gene']
        x_gene_fit = list(gene_df['x'])
        y_gene_fit = list(gene_df['y'])
        cells = list(louvain_df.index)
        cell_list=[]
        x_coord = []
        y_coord = []
        for i in range(len(cells)):
            if cells[i] in list(data_df.index):
                cell_list.append(cells[i])
                x_coord.append(umap_coord[i][0])
                y_coord.append(umap_coord[i][1])
        umap_df = pd.DataFrame(index=cell_list)
        umap_df['x'] = x_coord
        umap_df['y'] = y_coord

        import numpy as np
        from sklearn.linear_model import Lasso
        from sklearn.neighbors import KNeighborsRegressor
        import pickle
        cells=[]
        genes=[]
        gene_coord_x = []
        gene_coord_y = []
        for i in range(n_clusters):
            clust_data = data_df[data_df['cluster']==i]
            clust_cells = clust_data[clust_data['type']=='cell']
            clust_genes = clust_data[clust_data['type']=='gene']
            cells.extend(list(clust_cells.index))
            genes.extend(list(clust_genes.index))
            if len(list(clust_genes.index))==0:
                 continue
            model1 = KNeighborsRegressor(n_neighbors=5)

            model2 = KNeighborsRegressor(n_neighbors=5)

            model1.fit(np.array(list(clust_cells['x'])).reshape((-1, 1)),np.array(list(umap_df.loc[list(clust_cells.index)]['x'])).reshape((-1, 1)))

            filename = output_path+'/scanpy_x_KNN_model.sav'
            pickle.dump(model1, open(filename, 'wb'))
            #model1 = pickle.load(open(filename, 'rb'))
            x_gene_pred = model1.predict(np.array(list(clust_genes['x'])).reshape((-1, 1)))
#             gene_coord_x.extend(x_gene_pred)  #old
            gene_coord_x.extend(np.ravel(x_gene_pred).tolist())
            model2.fit(np.array(list(clust_cells['y'])).reshape((-1, 1)),np.array(list(umap_df.loc[list(clust_cells.index)]['y'])).reshape((-1, 1)))

            filename = output_path+'/scanpy_y_KNN_model.sav'
            pickle.dump(model2, open(filename, 'wb'))
            #model2 = pickle.load(open(filename, 'rb'))
            y_gene_pred = model2.predict(np.array(list(clust_genes['y'])).reshape((-1, 1)))
#             gene_coord_y.extend(y_gene_pred)  #old
            gene_coord_y.extend(np.ravel(y_gene_pred).tolist())

            



        with open(output_path+"/scanpy_gene_coord_x.txt", 'wb') as fp:
                pickle.dump(gene_coord_x,fp)
        with open(output_path+"/scanpy_gene_coord_y.txt", 'wb') as fp:
                pickle.dump(gene_coord_y,fp)

        #with open (output_path+"/scanpy_gene_coord_x.txt", 'rb') as fp:
        #    gene_coord_x = pickle.load(fp)
        #with open (output_path+"/scanpy_gene_coord_y.txt", 'rb') as fp:
        #    gene_coord_y = pickle.load(fp)

        #n_clusters = len(list(data_df['cluster'].unique()))

        u_map_x =[]
        u_map_y =[]
        for ind in list(data_df.index):
            if ind in list(louvain_df.index):

                u_map_x.append(umap_coord[list(louvain_df.index).index(ind)][0])
                u_map_y.append(umap_coord[list(louvain_df.index).index(ind)][1])
            else:
                u_map_x.append(gene_coord_x[genes.index(ind)])
                u_map_y.append(gene_coord_y[genes.index(ind)])
        data_df['umap_x'] = u_map_x
        data_df['umap_y'] = u_map_y


#        colors = random.sample(seaborn_colors,n_clusters)
        #colors = colors3
        plt.figure(figsize=(5,5))
        #cmap = sns.cubehelix_palette(dark=.3, light=.8, as_cmap=True)
        ax = sns.scatterplot(x="umap_x", y="umap_y", data=data_df,hue="cluster",palette=sns.xkcd_palette(colors),linewidth=0.0,s=2)
        ax.legend(bbox_to_anchor=(1.0, 1.00),frameon=False)
        for cl in range(n_clusters):
            plt.annotate(cl,data_df.loc[data_df['cluster']==cl,['umap_x','umap_y']].mean(),
                        horizontalalignment='center',
                        verticalalignment='center',
                        size=10, weight='bold',
                        color="black")
        sns.despine(bottom = False, left = False)
        plt.xlabel("umap1",fontsize=20)
        plt.ylabel("umap2",fontsize=20)
        plt.setp(ax.spines.values(), linewidth=2)
        plt.yticks([],linewidth=20)
        plt.xticks([])
        plt.savefig(output_path+'umap_clustering.png',bbox_inches='tight',dpi=600)
        plt.savefig(output_path+'umap_clustering.pdf',bbox_inches='tight',dpi=600)

        import matplotlib.pyplot as plt, mpld3
        from scipy.spatial import ConvexHull, convex_hull_plot_2d
        data_df["alpha"] = np.where(data_df['type'] == 'gene', 1.0, 0.5)
        color_gene = ["light grey"]
        color_cell = ["red"]
        #fig,ax1 = plt.subplots()
        plt.figure(figsize=(6,6))

        ax = sns.scatterplot(x="umap_x", y="umap_y", data=data_df[data_df['alpha']==0.5],hue="type",palette=sns.xkcd_palette(color_gene),sizes=(10,5),size="type",alpha=0.3,s=10)
        sns.scatterplot(x="umap_x", y="umap_y", data=data_df[data_df['alpha']==1.0],hue="type",palette=sns.xkcd_palette(color_cell),sizes=(20,5),size="type",marker="^",alpha=1.0,ax=ax,s=10)
        for c in range(n_clusters):
            p=data_df[data_df["cluster"]==c]
            p=p[['umap_x','umap_y']]
            points = p.values
            hull = ConvexHull(points)
            #for simplex in hull.simplices:
            #    sns.lineplot(points[simplex, 0], points[simplex, 1])
        x_list = []
        y_list = []
        for m in marker:
            x_list.append(data_df.loc[m]['umap_x'])
            #x_list.append(gene_coord_x[genes.index(m)])
            y_list.append(data_df.loc[m]['umap_y'])
            #y_list.append(gene_coord_y[genes.index(m)])
        for cl in range(n_clusters):
            plt.annotate(cl,data_df.loc[data_df['cluster']==cl,['umap_x','umap_y']].mean(),
                        horizontalalignment='center',
                        verticalalignment='center',
                        size=10, weight='bold',
                        color="black")
        for label, x, y in zip(marker, x_list, y_list):
            plt.annotate(
                label,
                xy=(x, y), xytext=(-20, 20),
                textcoords='offset points', ha='right', va='bottom',
                #bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                arrowprops=dict(arrowstyle = '-', connectionstyle='arc3,rad=0'))
        ax.legend(bbox_to_anchor=(1.0, 1.00),frameon=False)
        sns.despine(bottom = False, left = False)
        plt.xlabel("umap1",fontsize=20)
        plt.ylabel("umap2",fontsize=20)
        plt.setp(ax.spines.values(), linewidth=2)
        plt.yticks([],linewidth=20)
        plt.xticks([])
        plt.savefig(output_path+'umap_embedding.png',bbox_inches='tight',dpi=600)
        plt.savefig(output_path+'umap_embedding.pdf',bbox_inches='tight',dpi=600)


        import matplotlib.pyplot as plt, mpld3
        from scipy.spatial import ConvexHull, convex_hull_plot_2d
        data_df["alpha"] = np.where(data_df['type'] == 'gene', 1.0, 0.5)
        color_gene = ["light grey"]
        color_cell = ["red"]
        #fig,ax1 = plt.subplots()
        plt.figure(figsize=(6,6))
 #       colors = color
        ax = sns.scatterplot(x="umap_x", y="umap_y", data=data_df[data_df['alpha']==0.5],hue="cluster",linewidth=0.0,sizes=(2,5),size="type",palette=sns.xkcd_palette(colors),s=2)
        sns.scatterplot(x="umap_x", y="umap_y", data=data_df[data_df['alpha']==1.0],hue="type",palette=sns.xkcd_palette(color_cell),linewidth=0.1,marker="^",ax=ax,alpha=1.0,s=10)
        for c in range(n_clusters):
            p=data_df[data_df["cluster"]==c]
            p=p[['umap_x','umap_y']]
            points = p.values
            hull = ConvexHull(points)
            #for simplex in hull.simplices:
            #    sns.lineplot(points[simplex, 0], points[simplex, 1])
        x_list = []
        y_list = []
        for m in marker:
            x_list.append(data_df.loc[m]['umap_x'])
            y_list.append(data_df.loc[m]['umap_y'])
        for cl in range(n_clusters):
            plt.annotate(cl,data_df.loc[data_df['cluster']==cl,['umap_x','umap_y']].mean(),
                        horizontalalignment='center',
                        verticalalignment='center',
                        size=10, weight='bold',
                        color="black")
        for label, x, y in zip(marker, x_list, y_list):
            plt.annotate(
                label,
                xy=(x, y), xytext=(-20, 20),
                textcoords='offset points', ha='right', va='bottom',
                #bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                arrowprops=dict(arrowstyle = '-', connectionstyle='arc3,rad=0'))
        ax.legend(bbox_to_anchor=(1.0, 1.00),frameon=False)
        sns.despine(bottom = False, left = False)
        plt.xlabel("umap1",fontsize=20)
        plt.ylabel("umap2",fontsize=20)
        plt.setp(ax.spines.values(), linewidth=2)
        plt.yticks([],linewidth=20)
        plt.xticks([])
        plt.savefig(output_path+'umap_color_embedding.png',bbox_inches='tight',dpi=600)
        plt.savefig(output_path+'umap_color_embedding.pdf',bbox_inches='tight',dpi=600)







# In[44]:


input_path=""
output_path=""
de_analysis= "False"
file=""
n_pass = 1
arguments=sys.argv[1:]

if len(arguments)>=8 and arguments[0]=="-inp_data" and arguments[2]=="-i" and arguments[4]=="-o" and arguments[6]=="-n_pass":
    input_dataset=arguments[1]
    input_path=arguments[3]
    output_path=arguments[5]
    n_pass = int(arguments[7])
    if len(arguments) >8 and arguments[8]=="-de_analysis":
        de_analysis = arguments[9]
    #file=open(output_file,"w")

    compute(input_dataset,input_path,output_path,de_analysis,n_pass);
else:
        print ("\nPlease follow the command: <file>.py -inp_data inp_dataset -i inputpath -o outputpath \n")
        print('inputpath : path for read 10x data\n')
        print('outputpath : path where the output is to be stored\n')
        sys.exit()

