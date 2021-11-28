import streamlit as st

import pandas as pd
import numpy as np

import plotly.express as px

EV_data_file = "data/EV_gene_expr_avg_8points_gt001.tsv"
cell_data_file = "data/Cell_gene_expr_avg_3points_gt001.tsv"

@st.cache
def load_data(data_file,rows=[],cols=[],nrows=0):
    if nrows!=0:
        data = pd.read_csv(data_file,nrows=nrows, index_col=0, header=0, sep="\t")
        if rows != []:
            data = data.loc[data.index.isin(rows), :]
        if cols != []:
            data = data.loc[:, data.columns.isin(cols)]
    else:
        data = pd.read_csv(data_file, index_col=0,header=0,sep="\t")
        if rows != []:
            data = data.loc[data.index.isin(rows), :]
        if cols!=[]:
            data = data.loc[:,data.columns.isin(cols)]

    return data


def render_expr_page():
    # ===============================================================
    ## main content
    st.header('Expression data table' )

    # show table
    if st.checkbox('Show EV expression data'):
        EV_expr = load_data(EV_data_file)
        st.subheader('EV expression data')
        st.write(EV_expr)
    if st.checkbox('Show Cell expression data'):
        cell_expr = load_data(cell_data_file)
        st.subheader('Cell expression data')
        st.write(cell_expr)

    ## gene boxplot
    gene_id = "ENSG00000228794.4"
    gene_id = st.text_input(label="Query a gene",key="gene_id")

    if len(gene_id) >=2:
        EV_gene_expr = load_data(EV_data_file, rows=[gene_id])
        cell_gene_expr = load_data(cell_data_file, rows=[gene_id])

        if EV_gene_expr.shape[0] >=1:
            st.write(EV_gene_expr)
            EV =1
        else:
            EV = 0
            st.info("No EV expression record matches " + gene_id)

        if cell_gene_expr.shape[0] >=1:
            st.write(cell_gene_expr)

            cell = 1
        else:
            cell = 0
            st.info("No Cell expression record matches " + gene_id)

        fig = px.line()
        if EV:
            EV_x = [0, 1, 3, 5, 7, 9, 11, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]
            EV_y =  EV_gene_expr.loc[gene_id,:].to_numpy()
            fig.add_scatter(x=EV_x, y=EV_y,  line={"color": "red"},
                              name="EV",  mode='lines')
        if cell:
            cell_x = [0,11,30]
            cell_y = cell_gene_expr.loc[gene_id, :].to_numpy()
            fig.add_scatter(x=cell_x, y=cell_y, line={"color": "blue"},
                            name="cell", mode='lines')
        st.plotly_chart(fig,use_container_width=True)

    else:
        st.warning("Input at least 2 letters.")



def main():
    render_expr_page()