import pandas as pd
import numpy as np
from matplotlib import cm
from scipy.interpolate import make_interp_spline
from matplotlib.colors import rgb2hex

import plotly.express as px

import streamlit as st

@st.cache
def plot_cluster3D(clusters = [],gene_expr=[], membership=[]):
    color_ls=["Reds","Blues","Greens","Purples"]
    insert_n = 50
    range_y = 1
    cluster_dict = {}
    i = 0
    for cluster_i in clusters:
        if cluster_i !=0 and cluster_i not in cluster_dict:
            cluster_dict [cluster_i]=[color_ls[i],i*2+1]
            i += 1

    print(cluster_dict)
    y_max = len(cluster_dict)*2+1

    data_df = pd.DataFrame(columns=["gene_i","cluster_i","X","Y","Z","C"])
    for cluster_i in cluster_dict:
        cluster_i_genes = pd.read_csv("data/clusters/gene_list_"+str(cluster_i)+".txt",sep="\t",index_col=0,header=0)

        cluster_i_expr = gene_expr.loc[cluster_i_genes.index,:]
        cluster_i_expr = np.log10(cluster_i_expr + 0.01)
        cluster_i_expr = cluster_i_expr.apply(lambda x: (x-x.mean())/x.std(),axis=1)

        gene_ls = list(cluster_i_genes.index)

        ##
        color_map = cm.get_cmap(cluster_dict[cluster_i][0])(np.linspace(0,1,100))
        cluster_i_member = membership.loc[gene_ls,str(cluster_i)]

        x_data = list(int(i.split("_")[1][1:]) for i in cluster_i_expr.columns)
        for gene_i in gene_ls:
            z_data = cluster_i_expr.loc[gene_i,:]
            y_data = np.random.random_sample((len(x_data),)) * range_y + cluster_dict[cluster_i][1]

            x_space = np.linspace(np.min(x_data), np.max(x_data), insert_n)

            a_BSpline = make_interp_spline(x_data, z_data)
            z_space = a_BSpline(x_space)

            # b_BSpline = make_interp_spline(x_data, y_data)
            # y_space = b_BSpline(x_space)
            y_space = [np.random.random(1)[0] * range_y +cluster_dict[cluster_i][1]] * insert_n

            c = color_map[int(cluster_i_member[gene_i] * insert_n)]
            c_hex = rgb2hex(c)

            for i in range(insert_n):
                data_df.loc[len(data_df)] = [gene_i,cluster_i, x_space[i], y_space[i], z_space[i], c_hex]

    data_df.index = data_df["gene_i"]
    gene_ls= list(data_df.index.unique())
    fig = px.line_3d()
    for gene_i in gene_ls:
        # print(gene_i)
        fig.add_scatter3d(x=data_df.loc[gene_i,"X"], y=data_df.loc[gene_i,"Y"], z=data_df.loc[gene_i,"Z"],
                          line={"color":data_df.loc[gene_i,"C"].unique()[0]},
                          name = gene_i,
                          mode='lines')

    fig.update_layout(
        scene = dict(xaxis = dict(nticks=7, range=[0,30]),
                     yaxis = dict(nticks=6, range=[0,y_max]),
                     zaxis = dict(nticks=7, range=[-4,4]),
                     xaxis_title='Days',
                     yaxis_title='Cluster',
                     zaxis_title='Scaled exprression'),
        scene_aspectmode='manual',
        scene_aspectratio=dict(x=2, y=1, z=1)
    )
    return fig


# clusters = [76,15,0,4]
# EV_expr = data = pd.read_csv("../../data/EV_gene_expr_avg_8points_gt001.tsv", index_col=0,header=0,sep="\t")
# membership = pd.read_csv("../../data/clusters/mfuzz_cluster_membership_100.txt",
#                              sep="\t", index_col=0, header=0)
# fig = plot_cluster3D(clusters,EV_expr,membership)
# fig.write_html("3D_clusters_Ynochange_test.html")



