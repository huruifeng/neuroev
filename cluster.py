import streamlit as st
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import scipy
import scipy.spatial as sp, scipy.cluster.hierarchy as hc

import pickle

from matplotlib import cm
from func.plotly_3dlines import plot_cluster3D

cluster_data_file = "data/mfuzz_cluster_centers_100.txt"
EV_data_file = "data/EV_gene_expr_avg_8points_gt001.tsv"
membership_file = "data/clusters/mfuzz_cluster_membership_100.txt"

@st.cache
def load_data(data_file,rows=[],cols=[],nrows=0):
    if nrows!=0:
        data = pd.read_csv(data_file,nrows=nrows, index_col=0, header=0, sep="\t")
        if rows != []:
            data = data.loc[rows, :]
        if cols != []:
            data = data.loc[:, cols]
    else:
        data = pd.read_csv(data_file, index_col=0,header=0,sep="\t")
        if rows != []:
            data = data.loc[rows, :]
        if cols!=[]:
            data = data.loc[:,cols]

    return data


def render_cluster_page():
    # ===============================================================
    ## main content
    st.header('Gene expression pattern' )
    cluster_df = load_data(cluster_data_file)
    if st.checkbox('Show gene expression pattern cluster table'):
        st.subheader('Gene expression cluster data')
        st.write(cluster_df)

    title_str = '<span style="font-family:sans-serif; color:Black; font-size: 20px;">Cluster options</span>'
    st.markdown(title_str, unsafe_allow_html=True)
    cols = st.columns([4,4,3,5])
    with cols[0]:
        cluster_distance = st.selectbox('Distance',( "DTW", 'euclidean','hamming', 'cosine', 'cityblock', 'correlation'))
    with cols[1]:
        cluster_method = st.selectbox('Method', ('centroid','average', 'single', 'complete', 'median','ward'))
    with cols[2]:
        cluster_N = st.selectbox('Pick a pattern to show', list(range(1,101)))

    EV_expr = load_data(EV_data_file)
    membership = pd.read_csv(membership_file, sep="\t", index_col=0, header=0)

    col1,col2,col3 = st.columns([12,12,3])
    with col1:
        cluster_scaled = cluster_df.apply(lambda x: 2 * (x - x.min()) / (x.max() - x.min()) - 1, axis=1)
        if cluster_distance == "DTW":
            DTW_distance =  pickle.load(open("data/DTW_distance.pkl", "rb" ))
            linkage = hc.linkage(sp.distance.squareform(DTW_distance), method=cluster_method, optimal_ordering=False)
            g = sns.clustermap(cluster_scaled, cmap="vlag", center=0, figsize=(6, 14),
                               col_cluster=False, row_linkage=linkage, yticklabels=True,
                               dendrogram_ratio = (0.2,0.01))
        else:
            g = sns.clustermap(cluster_scaled, cmap="vlag", center=0, figsize=(6, 14),
                               method = cluster_method, metric = cluster_distance,
                               col_cluster=False, yticklabels=True,
                               dendrogram_ratio = (0.2,0.01))
        ax = g.ax_heatmap
        ax.set_ylabel("")
        ax.set_xticklabels(labels=ax.get_xticklabels(), fontsize=8)
        ax.set_yticklabels(labels=ax.get_yticklabels(), fontsize=7)
        g.fig.subplots_adjust(right=0.78)
        g.ax_cbar.set_position((0.85, .4, .03, .4))

        st.pyplot(g)
    with col2:
        cluster_i_genes = pd.read_csv("data/clusters/gene_list_" + str(cluster_N) + ".txt", sep="\t",
            index_col=0, header=0)

        cluster_i_expr = EV_expr.loc[cluster_i_genes.index, :]
        cluster_i_expr = np.log10(cluster_i_expr + 0.01)
        cluster_i_expr = cluster_i_expr.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
        gene_ls = list(cluster_i_genes.index)

        color_map = cm.get_cmap("Reds")(np.linspace(0, 1, 101))
        cluster_i_member = membership.loc[gene_ls, str(cluster_N)]

        fig = plt.figure(figsize=(4, 4))
        ax = fig.add_subplot(111)

        x_data = list(int(i.split("_")[1][1:]) for i in cluster_i_expr.columns)
        for gene_i in gene_ls:
            z_data = cluster_i_expr.loc[gene_i, :]

            x_space = np.linspace(np.min(x_data), np.max(x_data), 100)
            a_BSpline = scipy.interpolate.make_interp_spline(x_data, z_data)
            z_space = a_BSpline(x_space)

            ax.plot(x_space, z_space,lw=1, c=color_map[int(cluster_i_member[gene_i] * 100)])
        ax.tick_params(labelsize=6)
        ax.set_ylim(-4, 4)
        ax.set_xlabel("Day",size=6)
        ax.set_ylabel("Scaled expression",size=6)
        ax.set_title("Gene expression pattern in cluster " + str(cluster_N),fontsize=8)

        st.pyplot(fig)

    title_str = '<span style="font-family:sans-serif; color:Black; font-size: 20px;">3D cluster</span>'
    st.markdown(title_str, unsafe_allow_html=True)
    cols = st.columns([3, 3, 3, 3])
    with cols[0]:
        cluster_1 = st.selectbox('Pick 1st cluster', list(range(0, 101)))
    with cols[1]:
        cluster_2 = st.selectbox('Pick 2nd cluster', list(range(0, 101)))
    with cols[2]:
        cluster_3 = st.selectbox('Pick 3rd cluster', list(range(0, 101)))
    with cols[3]:
        # cluster_4 = st.selectbox('Pick 4th cluster', list(range(0, 101)))
        cluster_4 = 0

    if cluster_1 !=0 or cluster_3 !=0 or cluster_3 !=0 or cluster_4 !=0:
        clusters_tmp = [cluster_1,cluster_2,cluster_3,cluster_4]
        clusters = [cluster_i for cluster_i in clusters_tmp if cluster_i !=0]
        if len(clusters)==3:
            fig = plot_cluster3D(clusters,EV_expr,membership)
            st.plotly_chart(fig,use_container_width=True)


def main():
    render_cluster_page()
