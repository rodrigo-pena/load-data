# load-data
Repository for maintaining functions and scripts to load and display data. The main focus is on data that can be represented as graphs (or networks).

## spectral-layout
Contains functions for assigning new coordinates to the nodes of a graph. It does so by computing (or approximating) the eigenvectors of the graph Laplacian (L), or the generalized eigenvectors of (L,D), where D is the diagonal degree matrix of the graph.

###References
1. Y. Koren, "Drawing Graphs by Eigenvectors," Computers and Mathematics with Applications, vol. 49, pp. 1867-1888, 2005.