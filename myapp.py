import base64

import streamlit as st
import pandas as pd
import numpy as np

st.set_page_config(
    page_title="NeuroEV",
    # page_icon="favicon.ico",
    layout = 'wide',
    initial_sidebar_state="expanded",
)
############################################

from .pages import home
from .pages import expression
from .pages import cluster
from .pages import de


## Functins / Style
st.markdown(
    '''
        <style>
             img.logo {
                    margin: 0px auto;
                    margin-top: -60px;
                    margin-bottom: 25px;
                    width: 60px;
                }
        </style>
    ''',
    unsafe_allow_html=True
)

def local_image(url):
    file_ = open(url, "rb")
    contents = file_.read()
    data_url = base64.b64encode(contents).decode("utf-8")
    file_.close()
    return data_url

PAGES = {
    "Home": home,
    "Expression Data": expression,
    "Gene Clusters": cluster,
    "DE Analysis": de
}

##================================================================
def main():
    ## Sidebar
    col1, col2 = st.sidebar.columns([3,9])
    with col1:
        st.image('images/side_img.png', width=70)
    with col2:
        title_str='<span style="font-family:sans-serif; color:Black; font-size: 50px;">NeuroEV</span>'
        st.markdown(title_str,unsafe_allow_html=True)

    st.sidebar.title("Navigation")
    select_page = st.sidebar.radio("Go to page", list(PAGES.keys()),key="page_radio")

    page = PAGES[select_page]

    st.sidebar.title("About")
    st.sidebar.info(
        """
        This app is used for exploring the gene expression ptterns in Human iPSC-derived Exosomes 
        in Dopamine Neuron Differentiation.
        """
    )

    page.main()

if __name__ == "__main__":
    main()







