function plot_chernobyl(G, x, unit)
%PLOT_ETEX plots Chernobyl graph G and signal x on a European map
%
%   Usage:
%       plot_chernobyl(G, x)
%
%   Input:
%       G   : A Matlab structure encoding graph information.
%       x   : A G.N-by-1 vector representing the signal on the graph.
%       unit: (Optional) A string specifying the unit of measurement of the
%             signal on the graph.
%
%   Output:
%
%   Example:
%       [G, measurements] = chernobyl();
%       plot_chernobyl(G, measurements.Cs137_fallout);
%
%   Requires:  Files ne_110m_land.* on path (see reference below for 
%              the download link)
%
%   Reference: https://rem.jrc.ec.europa.eu/RemWeb/Browse.aspx?path=Chernobyl%20Data
%              http://www.naturalearthdata.com/downloads/
%
% Author: Rodrigo Pena (rodrigo.pena@epfl.ch)
% Date: 11 Mar 2016

%% Parse input
assert(isfield(G, 'coords'), 'G.coords doens''t exist');
assert(isfield(G, 'N'), 'G.N doens''t exist');
assert(sum(size(x) ~= [G.N, 1]) == 0, 'x must be a G.N-by-1 vector');
if nargin < 3; unit = []; end    

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

if G.Ne <= 1000
    gsp_plot_edges(G);
end

hold off
axis([-15, 40, 35, 75]); % Focus on Europe
axis off
h = colorbar;
set(h, 'TickLabels', linspace(min(x), max(x), length(h.Ticks))');
if ~isempty(unit); set(get(h,'Title'),'String', unit); end
set(gcf, 'Position', [500 500 960 720]);

end