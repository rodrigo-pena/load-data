# -*- coding: utf-8 -*-
"""
Functions for loading and displaying the Higgs-Twitter Dataset

If used by itself, the syntax is
    $ python higgs-twitter.py HIGGS_DIR

Author: Rodrigo Pena
E-mail: rodrigo.pena@ephl.ch
Created on 30 Oct 2015

"""
import os
import sys
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


def higgs_social_network(HIGGS_DIR, directed=1):
    '''Assemble graph from higgs-social_network.edgelist

    :parameters:
        - HIGGS_DIR : string
          Directory where the file 'higgs-social_network.edgelist' is located

        - directed : int = 0, 1
          Assemble directed (1) or undirected (0) graph

    :returns:
        - G : networkx graph
          Unweighted graph whose node connections are specified in
          'higgs-social_network.edgelist'

    '''
    path = os.path.join(HIGGS_DIR, 'higgs-social_network.edgelist')
    G = nx.read_edgelist(path)
    return G

if __name__ == '__main__':
    HIGGS_DIR = sys.argv[1]
    G = higgs_social_network(HIGGS_DIR, 1)
    nx.draw(G, pos=nx.spring_layout(G))  # Use spring layout
    plt.show()
