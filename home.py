import streamlit as st

def render_home_page():
    # ===============================================================
    ## main content
    st.header('Transcriptional Dynamics of Human iPSC-derived Exosomes in Dopamine Neuron Differentiation')
    cols = st.columns(4)
    with cols[0]:
        st.image('images/ipsc1.png', width=200)
    with cols[1]:
        st.image('images/ipsc2.png', width=200)
    with cols[2]:
        st.image('images/ipsc3.png', width=200)
    with cols[3]:
        st.image('images/ipsc4.png', width=200)

def main():
    render_home_page()