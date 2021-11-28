import streamlit as st
import pandas as pd
import numpy as np

EV_data_file = "data/EV_gene_expr_avg_8points_gt001.tsv"
cell_data_file = "data/Cell_gene_expr_avg_3points_gt001.tsv"

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

# Load 10,000 rows of data into the dataframe.
EV_expr = load_data(EV_data_file,nrows=1000)
cell_expr = load_data(cell_data_file,nrows=1000)
# Notify the reader that the data was successfully loaded.
# st.text("Done! (using st.cache)")

def render_expr_page():
    # ===============================================================
    ## main content
    st.header('Expression data table' )

    if st.checkbox('Show EV expression data'):
        st.subheader('EV expression data')
        st.write(EV_expr)
    if st.checkbox('Show Cell expression data'):
        st.subheader('Cell expression data')
        st.write(cell_expr)

def main():
    with st.spinner("Loading page..."):
        render_expr_page()