
import pandas as pd
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score

def load_clusters():
    stardust_path = "Stardust_results/visualization_output/3_pass/data.csv"
    louvain_path = "Stardust_results/build_output/1_pass/louvain_cluster_df.csv"

    # Load Stardust clusters
    stardust_df = pd.read_csv(stardust_path).set_index('data')
    stardust_clusters = stardust_df[stardust_df['type'] == 'cell']['cluster']

    # Load Louvain clusters
    louvain_df = pd.read_csv(louvain_path).set_index('Unnamed: 0')
    louvain_df.index = louvain_df.index.map(str)
    louvain_clusters = louvain_df['louvain']

    # Get common cell indices
    common_cells = stardust_clusters.index.intersection(louvain_clusters.index)

    return stardust_clusters.loc[common_cells], louvain_clusters.loc[common_cells]

def compute_nmi(stardust_clusters, louvain_clusters):
    nmi_score = normalized_mutual_info_score(stardust_clusters, louvain_clusters)
    print(f"NMI between Stardust and Louvain: {nmi_score:.4f}")

def compute_ari(stardust_clusters, louvain_clusters):
    ari_score = adjusted_rand_score(stardust_clusters, louvain_clusters)
    print(f"ARI between Stardust and Louvain: {ari_score:.4f}")

if __name__ == "__main__":
    stardust_clusters, louvain_clusters = load_clusters()
    compute_nmi(stardust_clusters, louvain_clusters)
    compute_ari(stardust_clusters, louvain_clusters)





