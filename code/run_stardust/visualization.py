#!/usr/bin/env python
# coding: utf-8

# In[17]:

import sys
# In[18]:


import seaborn
seaborn_colors=[]
for k in seaborn.xkcd_rgb.keys():
    seaborn_colors.append(k)


def compute(inp_dataset,input_path,output_path,de_analysis,n_pass):

    print("Current pass ",n_pass)
    import csv
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    from sklearn.cluster import DBSCAN
    import operator
    import numpy as np
    import random

    #csvData=[['data','x','y','type']]
    print("Processing the input data into datafames....")
    csvData=[]
    count=0
    filename = input_path+"/output_normalized_own_cc.csv"
    coord_data = pd.read_csv(filename,names=['data','x','y'])
    coord_data.set_index('data',inplace=True)
    data=[]
    data_outlier=[]
    import pickle
    with open ('Stardust_results/build_output/1_pass/cell_list.txt', 'rb') as fp:
            cell_list = pickle.load(fp)
    with open ('Stardust_results/build_output/1_pass/gene_list.txt', 'rb') as fp:
            gene_list = pickle.load(fp)

    with open(filename, 'r') as csvfile:
         csvreader = csv.reader(csvfile)
         for row in csvreader:
            data.append(row)
            temp_outlier=[]
            temp_outlier.append(row[1])
            temp_outlier.append(row[2])
            data_outlier.append(temp_outlier)
            temp=row
            if row[0] in cell_list:
                temp.append('cell')
            else:
                temp.append('gene')
                count=count+1
            csvData.append(temp)


    # # DB SCAN

    # In[20]:


    if n_pass !=4:
        noise =[]
        print("Performing clustering....")
        #36
        #db = DBSCAN(eps=35,min_samples=40).fit_predict(data_outlier) #liver
        db = DBSCAN(eps=180,min_samples=55).fit_predict(data_outlier) #melanoma
        final_data=[]
        csvData=[['data','x','y','type']]
        for i in range(0,len(list(db))):
            if db[i]!=-1:
                final_data.append(data[i])
                csvData.append(data[i])
            if db[i]==-1:
                noise.append(data[i][0])
        data=final_data
        print("outliers count ",len(noise))
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
                import pickle
                with open ('stardust/Stardust_results/visualization_output/3_pass/de_genes_cluster.txt', 'rb') as fp:
                    de_gene_cluster = pickle.load(fp)
                for rank in range(0,len(de_gene_cluster)):
                    if csvData[i][0] in de_gene_cluster[rank]:
                        clusters_info.append(de_gene_cluster[rank].index(csvData[i][0]))
                        break
        for r in remove_data:
            csvData.remove(r)
        temp=[['data','x','y','type']]
        temp.extend(csvData)
        csvData = temp


    # In[13]:



    # # POST-HOC CLUSTER ASSIGNMENT

    # In[23]:


    print("Starting post hoc clustering....")
    neighbor_df = pd.read_hdf('Stardust_results/build_output/1_pass/neighbor.h5', 'df')
    if 'Unnamed: 0' in list(neighbor_df.columns):
        neighbor_df.set_index('Unnamed: 0',inplace=True)
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
            if db[list(coord_data.index).index(cell)] == max_cluster:
                count=count+1
                x_total=x_total+coord_data.loc[cell]['x']
                y_total=y_total+coord_data.loc[cell]['y']
        temp=[]
        temp.append(key_cell)
        temp.append(x_total/count)
        temp.append(y_total/count)
        temp.append('cell')
        csvData.append(temp)
    print("Post hoc clustering done....")


    # In[24]:



    with open(output_path+'data.csv', 'w') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(csvData)
    csvFile.close()
    data_df = pd.read_csv(output_path+"data.csv",delimiter=",",index_col=False)
    clusters_info = [x for x in db if x!=-1]
    clusters_info = clusters_info + cluster_assign
    data_df['cluster'] = clusters_info
    data_df.to_csv(output_path+'data.csv')
    data_df.to_csv(output_path+'data_openOrd.csv')
    print("cluster saved ....")





    # In[26]:


    colors = random.sample(seaborn_colors,n_clusters)
    fig=plt.figure(figsize=(5,5))
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
    fig.savefig(output_path+"cluster_visualization_with_umap.png",bbox_inches='tight',dpi=600);
    fig.savefig(output_path+"cluster_visualization_with_umap.pdf",bbox_inches='tight',dpi=600);




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




    import matplotlib.pyplot as plt
    from scipy.spatial import ConvexHull
    data_df["alpha"] = np.where(data_df['type'] == 'gene', 1.0, 0.5)
    color_gene = ["light blue"]
    color_cell = ["red"]
    plt.figure(figsize=(5,5))
    ax = sns.scatterplot(x="x", y="y", data=data_df[data_df['alpha']==0.5],hue="type",palette=sns.xkcd_palette(color_gene),sizes=(10,5),size="type",alpha=0.3,s=10)
    sns.scatterplot(x="x", y="y", data=data_df[data_df['alpha']==1.0],hue="type",palette=sns.xkcd_palette(color_cell),sizes=(20,5),size="type",marker="^",alpha=1.0,ax=ax,s=10)
    for c in range(n_clusters):
        p=data_df[data_df["cluster"]==c]
        p=p[['x','y']]
        points = p.values
        hull = ConvexHull(points)
    x_list = []
    y_list = []
    if n_pass !=1:
        for m in marker:
            x_list.append(data_df.loc[m]['x'])
            y_list.append(data_df.loc[m]['y'])
            #x_list.append(gene_coord_x[genes.index(m)])
            #y_list.append(gene_coord_y[genes.index(m)])
        for label, x, y in zip(marker, x_list, y_list):
            plt.annotate(
                 label,
                 xy=(x, y), xytext=(-20, 20),
                 textcoords='offset points', ha='right', va='bottom',
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


    import matplotlib.pyplot as plt
    from scipy.spatial import ConvexHull
    data_df["alpha"] = np.where(data_df['type'] == 'gene', 1.0, 0.5)
    color_gene = ["light blue"]
    color_cell = ["red"]
    plt.figure(figsize=(5,5))
    ax = sns.scatterplot(x="x", y="y", data=data_df[data_df['alpha']==0.5],hue="cluster",palette=sns.xkcd_palette(colors),sizes=(10,5),size="type",alpha=0.3,s=10)
    sns.scatterplot(x="x", y="y", data=data_df[data_df['alpha']==1.0],hue="type",palette=sns.xkcd_palette(color_cell),sizes=(20,5),size="type",marker="^",alpha=1.0,ax=ax,s=10)
    for c in range(n_clusters):
        p=data_df[data_df["cluster"]==c]
        p=p[['x','y']]
        points = p.values
        hull = ConvexHull(points)
    x_list = []
    y_list = []
    if n_pass !=1:
        for m in marker:
            x_list.append(data_df.loc[m]['x'])
            y_list.append(data_df.loc[m]['y'])
            #x_list.append(gene_coord_x[genes.index(m)])
            #y_list.append(gene_coord_y[genes.index(m)])
        for label, x, y in zip(marker, x_list, y_list):
            plt.annotate(
                 label,
                 xy=(x, y), xytext=(-20, 20),
                 textcoords='offset points', ha='right', va='bottom',
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
# # UMAP CELL GENE MARKER # #
    #sys.exit(0)
    if n_pass ==4:

        import pickle
        with open ('stardust/Stardust_results/build_output/1_pass/umap_coord.txt', 'rb') as fp:
            umap_coord = pickle.load(fp)
        louvain_df = pd.read_csv('Stardust_results/build_output/1_pass/louvain_cluster_df.csv')
        louvain_df.set_index('Unnamed: 0',inplace=True)
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
            model1 = KNeighborsRegressor(n_neighbors=5)

            model2 = KNeighborsRegressor(n_neighbors=5)

            model1.fit(np.array(list(clust_cells['x'])).reshape((-1, 1)),np.array(list(umap_df.loc[list(clust_cells.index)]['x'])).reshape((-1, 1)))

            filename = output_path+'/x_KNN_model.sav'
            pickle.dump(model1, open(filename, 'wb'))
            #model1 = pickle.load(open(filename, 'rb'))
            x_gene_pred = model1.predict(np.array(list(clust_genes['x'])).reshape((-1, 1)))
            gene_coord_x.extend(x_gene_pred)
            model2.fit(np.array(list(clust_cells['y'])).reshape((-1, 1)),np.array(list(umap_df.loc[list(clust_cells.index)]['y'])).reshape((-1, 1)))

            filename = output_path+'/y_KNN_model.sav'
            pickle.dump(model2, open(filename, 'wb'))
            #model2 = pickle.load(open(filename, 'rb'))
            y_gene_pred = model2.predict(np.array(list(clust_genes['y'])).reshape((-1, 1)))
            gene_coord_y.extend(y_gene_pred)

        #n_clusters = len(list(data_df['cluster'].unique()))

        u_map_x =[]
        u_map_y =[]
        for ind in list(data_df.index):
            if ind in list(louvain_df.index):
            #print(ind)
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

        import matplotlib.pyplot as plt
        from scipy.spatial import ConvexHull
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
        plt.savefig(output_path+'umap_embedding.png',bbox_inches='tight',dpi=600)
        plt.savefig(output_path+'umap_embedding.pdf',bbox_inches='tight',dpi=600)


        import matplotlib.pyplot as plt
        from scipy.spatial import ConvexHull
        data_df["alpha"] = np.where(data_df['type'] == 'gene', 1.0, 0.5)
        color_gene = ["light grey"]
        color_cell = ["red"]
        plt.figure(figsize=(6,6))
        ax = sns.scatterplot(x="umap_x", y="umap_y", data=data_df[data_df['alpha']==0.5],hue="cluster",palette=sns.xkcd_palette(colors),sizes=(10,5),size="type",alpha=0.2,s=10)
        sns.scatterplot(x="umap_x", y="umap_y", data=data_df[data_df['alpha']==1.0],hue="type",palette=sns.xkcd_palette(color_cell),sizes=(20,5),size="type",marker="^",alpha=1.0,ax=ax,s=20)
        for c in range(n_clusters):
            p=data_df[data_df["cluster"]==c]
            p=p[['umap_x','umap_y']]
            points = p.values
            hull = ConvexHull(points)
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


    # In[50]:





     # # DE GENES FINDING
    import anndata
    import scanpy as sc
    de_analysis = True
    f=True
    if f == True and n_pass !=4:

        print("Finding top DE genes....");
        #heat_map_df = pd.read_csv("Stardust_results/build_output/1_pass/heat_map_df.csv",delimiter=",",index_col=False)
        heat_map_df = pd.read_hdf('Stardust_results/build_output/1_pass/heat_map_df.h5', 'df')
        heat_map_df.fillna(0)
        if "Unnamed: 0" in list(heat_map_df.columns):
            heat_map_df.set_index(["Unnamed: 0"],inplace=True)
        data_df = pd.read_csv(output_path+"/data.csv",delimiter=",",index_col=False)
        data_df.set_index('data',inplace=True)
        cells = list(data_df[data_df['type']=='cell'].index)
        genes = list(data_df[data_df['type']=='gene'].index)
        cluster = []
        cell_order =[]
        heat_map_df = heat_map_df.loc[cells]
        #data_df = data_df.sort_values(by=['cluster'])
        #index_list = data_df.index
        #for ind in index_list:
        #    if data_df.loc[ind]['type'] == 'cell':
        #         cell_order.append(ind)
        #         cluster.append(data_df.loc[ind]['cluster'])
        #heat_map_df = heat_map_df.reindex(cell_order)

        index_list = heat_map_df.index
        for ind in index_list:
                 cluster.append(data_df.loc[ind]['cluster'])

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
        #sc.pp.log1p(adata1)
        #sc.pp.normalize_total(adata1, target_sum=None, inplace=True)
        #sc.pp.log1p(adata1)
        sc.tl.rank_genes_groups(adata1,'cluster',n_genes=len(list(heat_map_df.columns)),  method='t-test_overestim_var',rankby_abs=True)

        genes_df = pd.DataFrame(adata1.uns['rank_genes_groups']['names'])
        p_val_df = pd.DataFrame(adata1.uns['rank_genes_groups']['pvals_adj'])
        log_fold_df = pd.DataFrame(adata1.uns['rank_genes_groups']['logfoldchanges'])
        DE_genes = []
        col = list(genes_df.columns)
        row = list(genes_df.index)
        de_gene_cluster = []
        for r in row:
            temp = []
            for c in col:
                if p_val_df.loc[r][c]<=0.05 and log_fold_df.loc[r][c]>=1.2:
                    if genes_df.loc[r][c] not in DE_genes and len(DE_genes)<500:
                        DE_genes.append(genes_df.loc[r][c])
                    temp.append(genes_df.loc[r][c])
                else:
                    temp.append("NULL")
            de_gene_cluster.append(temp)
        de_genes = list(set(DE_genes))
        print("No. of de genes is ",len(de_genes))
        import pickle
        with open(output_path+'/de_genes.txt', 'wb') as fp:
            pickle.dump(de_genes, fp)
        with open(output_path+'/de_genes_cluster.txt', 'wb') as fp:
            pickle.dump(de_gene_cluster, fp)







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

