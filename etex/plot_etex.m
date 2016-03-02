function plot_etex(G, x)
%PLOT_ETEX plots ETEX graph G and signal x on a European map
%
%   Usage:
%       plot_etex(G, x)
%
%   Input:
%       G   : A Matlab structure encoding graph information.
%       x   : The signal on the graph.
%
%   Output:
%
%   Example:
%       [G, B] = etex_graph();
%       plot_snow_gis(G, B(:, 15));
%
%   Requires:
%
%   Reference: https://rem.jrc.ec.europa.eu/RemWeb/etex/
%              http://www.arcgis.com/home/item.html?id=6d611f8d87d54227b494d4c3becef6a0
%
% Author: Rodrigo Pena (rodrigo.pena@epfl.ch)
% Date: 2 Mar 2016

%% Parse input
assert(isfield(G, 'coords'), 'G.coords doens''t exist');
assert(isfield(G, 'N'), 'G.N doens''t exist');
assert(sum(size(x) ~= [G.N, 1]) == 0, 'x must be a G.N-by-1 vector');

%% Initialization
% TODO: Check if there's a better shapefile online
% Read Europe map shape file 
S = shaperead('MyEurope.shp'); 

% Define colormap
cmap = jet(256);

% Normalize of color (signal) for plotting
color  = (x - min(x))./(max(x) - min(x)); 
color = floor(255.*color + 1);

%% Display graph & signals
mapshow(S, 'FaceColor', [1 1 1], 'EdgeColor', 'black')
colormap(cmap)
hold on
scatter3(G.coords(:,1), G.coords(:,2), color, 200, cmap(color, :), '.');
if isfield(G, 'idx_release_site')
    scatter3(G.coords(G.idx_release_site, 1), ...
        G.coords(G.idx_release_site, 2), color(G.idx_release_site), ...
        100, 'r', 'o');
end
hold off
axis([-11, 30, 35, 70]);
axis off
h = colorbar;
set(h, 'TickLabels', linspace(min(x), max(x), length(h.Ticks))');
set(get(h,'Title'),'String','PMCH ng/m^3');
%TODO: find a way to refresh the plot when we call plot_etex successively

end