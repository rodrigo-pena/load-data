function plot_etex(G, x)
%PLOT_ETEX plots ETEX graph G and signal x on a European map
%
%   Usage:
%       plot_etex(G, x)
%
%   Input:
%       G   : A Matlab structure encoding graph information.
%       x   : A G.N-by-1 vector representing the signal on the graph.
%
%   Output:
%
%   Example:
%       [G, B] = etex_graph();
%       plot_snow_gis(G, B(:, 15));
%
%   Requires:  Files ne_110m_land.* on path (see reference below for 
%              the download link)
%
%   Reference: https://rem.jrc.ec.europa.eu/RemWeb/etex/
%              http://www.naturalearthdata.com/downloads/
%
% Author: Rodrigo Pena (rodrigo.pena@epfl.ch)
% Date: 2 Mar 2016

%% Parse input
assert(isfield(G, 'coords'), 'G.coords doens''t exist');
assert(isfield(G, 'N'), 'G.N doens''t exist');
assert(sum(size(x) ~= [G.N, 1]) == 0, 'x must be a G.N-by-1 vector');

%% Initialization
% Read World Map shape file 
S = shaperead('ne_110m_land.shp'); 

% Define colormap
cmap = jet(256);

% Normalize of color (signal) for plotting
if nnz(x) > 0
    color  = (x - min(x))./(max(x) - min(x));
else
    color = x;
end
color = floor(255.*color + 1);

% Set edge color
G.plotting.edge_color = [0.6 0.6 0.6];

%% Display graph & signals
clf(gcf); % To avoid having mapshow add another layer to the current plot

mapshow(S, 'FaceColor', [1 1 1], 'EdgeColor', 'black')
colormap(cmap)
hold on

scatter3(G.coords(:,1), G.coords(:,2), color./10, 500, cmap(color, :), '.');
if isfield(G, 'idx_release_site')
    scatter3(G.coords(G.idx_release_site, 1), ...
        G.coords(G.idx_release_site, 2), color(G.idx_release_site)./10, ...
        100, 'r', 'o');
end

if G.Ne <= 1000
    gsp_plot_edges(G);
end

hold off
axis([-11, 30, 35, 70]); % Focus on Europe
axis off
h = colorbar;
set(h, 'TickLabels', linspace(min(x), max(x), length(h.Ticks))');
set(get(h,'Title'),'String','PMCH ng/m^3');

end