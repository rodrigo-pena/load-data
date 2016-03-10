# load-data
Repository for maintaining functions and scripts to load and display data. The main focus is on data that can be represented as graphs (or networks).

## etex
Contains functions for assembling and plotting a graph from the data the the European Tracer Experiment (ETEX), gathered in the middle of the 1990's. ETEX essentially consisted of two releases to atmosphere of tracers (perfluorocarbons) sampled for three days after the beginning of the emission using a sampling network spread over a large part of Europe [1].

In order to use these functions, you'll need to install the [GSPBox](https://lts2.epfl.ch/gsp/).

###References
1. "European Tracer Experiment (ETEX)," [Online](https://rem.jrc.ec.europa.eu/RemWeb/etex/).

## snow-gis
Contains functions for assembling and plotting a graph from the data of a 1854 cholera outbreak in Soho, London, UK. The data was gathered by physician John Snow, considered the father of epidemiology.

In order to use these functions, you'll need to install the [GSPBox](https://lts2.epfl.ch/gsp/) and [MatlabBGL](http://dgleich.github.io/matlab-bgl/).

###References
1. Robin Wilson "John Snow’s famous cholera analysis data in modern GIS formats," [Online](http://blog.rtwilson.com/john-snows-famous-cholera-analysis-data-in-modern-gis-formats/).

## spectral-layout
Contains functions for assigning new coordinates to the nodes of a graph. It does so by computing (or approximating) the eigenvectors of the graph Laplacian (L), or the generalized eigenvectors of (L,D), where D is the diagonal degree matrix of the graph.

###References
1. Y. Koren, "Drawing Graphs by Eigenvectors," Computers and Mathematics with Applications, vol. 49, pp. 1867-1888, 2005.

