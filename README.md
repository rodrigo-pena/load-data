# load-data
Repository for maintaining functions and scripts to load and display data. The main focus is on data that can be represented as graphs (or networks).

For MATLAB users, by running the script [load_data_start.m](https://github.com/rodrigo-pena/load-data/tree/master/load_data_start.m), you automatically add all the folders of the **load-data** repository to your path.

We provide below some detailed information about each subfolder in the **load-data** repository.

## [chernobyl](https://github.com/rodrigo-pena/load-data/tree/master/chernobyl)
Contains functions for assembling and plotting a graph from the [Chernobyl incident](https://en.wikipedia.org/wiki/Chernobyl_disaster) particle deposition data.

### Requirements
In order to use these functions, you'll need to install the [GSPBox](https://lts2.epfl.ch/gsp/).

To get your local Chernobyl data folder to be properly visible to the functions, follow these steps:
1. Go to this [link][chernobyl_link] and download the files *CUMDEP.DAT*, *NORWAY.UPD*, *POLAND.UPD*, and *RUMANIA.UPD*.
2. Save those files into the same directory, e.g., "chernobyl-data/". 
3. Provide the full path to this directory when running the Chernobyl functions.


### References
1. "Chernobyl deposition data," [Online][chernobyl_link]. Institute for Transuranium Elements, a European Comission Joint Research Centre.

[chernobyl_link]: https://rem.jrc.ec.europa.eu/RemWeb/Browse.aspx?path=\Chernobyl%20Data\Deposition


## [etex](https://github.com/rodrigo-pena/load-data/tree/master/etex)
Contains functions for assembling and plotting a graph from the data the the European Tracer Experiment (ETEX), gathered in the middle of the 1990's. ETEX essentially consisted of two releases to atmosphere of tracers (perfluorocarbons) sampled for three days after the beginning of the emission using a sampling network spread over a large part of Europe [1].

### Requirements
In order to use these functions, you'll need to install the [GSPBox](https://lts2.epfl.ch/gsp/).

To get your local ETEX data folder to be properly visible to the functions, follow these steps:
1. Go to this [link][etex_link] and scroll down on the left of the page until you see the link "Datasets". Click on it. 
2. Download all the files in "Tracer data for release 1" into a folder that you should name "release-1". 
3. Similarly, save the files in "Tracer data for release 2" into folder "release-2". 
4. Save both folders on a parent directory, e.g., "etex-data/"
5. Provide the full path to this directory to the ETEX functions.

### References
1. "European Tracer Experiment (ETEX)," [Online][etex_link].

[etex_link]: https://rem.jrc.ec.europa.eu/RemWeb/etex/.

## [misc](https://github.com/rodrigo-pena/load-data/tree/master/misc)
Contains some miscellaneous functions dealing with graph data.

## [snow-gis](https://github.com/rodrigo-pena/load-data/tree/master/snow-gis)
Contains functions for assembling and plotting a graph from the data of a 1854 cholera outbreak in Soho, London, UK [1]. The data was gathered by physician John Snow, considered the father of epidemiology.

### Requirements
In order to use these functions, you'll need to install:
* [GSPBox](https://lts2.epfl.ch/gsp/)
* [MatlabBGL](http://dgleich.github.io/matlab-bgl/).

To get your local Snow GIS data folder to be properly visible to the functions, follow these steps:
1. Go to this [link][snow_link] and download the file *SnowGIS.zip*.
2. Unzip the file and save its contents into a directory, e.g., "snow-gis-data/". 
3. Provide the full path to this directory when running the Snow GIS functions.

### References
1. Robin Wilson "John Snow’s famous cholera analysis data in modern GIS formats," [Online][snow_link].

[snow_link]: http://blog.rtwilson.com/john-snows-famous-cholera-analysis-data-in-modern-gis-formats/

## spectral-layout
Contains functions for drawing graphs with spectral-based layouts [1]. It does so by computing (or approximating) the eigenvectors of the graph Laplacian (L), or the generalized eigenvectors of (L,D), where D is the diagonal degree matrix of the graph.

### References
1. Y. Koren, ["Drawing Graphs by Eigenvectors"](https://www.math.ucdavis.edu/~saito/data/acha.read.w12/koren-graph-drawing.pdf) Computers and Mathematics with Applications, vol. 49, pp. 1867-1888, 2005.

