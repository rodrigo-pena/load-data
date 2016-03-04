function plot_snow_gis(G, x, SNOW_GIS_DIR)
%PLOT_SNOW_GIS plots graph and signal on Soho map
%
%   Usage:
%       plot_snow_gis(G, x, SNOW_GIS_DIR)
%
%   Input:
%       G               : A Matlab structure encoding graph information.
%       x               : The signal on the graph.
%       SNOW_GIS_DIR    : A string specifying the directory where the GIS
%                       dataset is located (see reference below).
%                       (DEFAULT: '~/data/snow_gis/');
%
%   Output:
%
%   Example:
%       [G, x, b] = snow_gis();
%       plot_snow_gis(G, b)
%
%   Requires:
%
%   Reference: http://blog.rtwilson.com/john-snows-famous-cholera-analysis-data-in-modern-gis-formats/
%
% Author: Rodrigo Pena (rodrigo.pena@epfl.ch)
% Date: 8 Feb 2016

%% Parse input
if nargin < 3 || isempty(SNOW_GIS_DIR)
    SNOW_GIS_DIR = '~/data/snow_gis/';
end
assert(isa(SNOW_GIS_DIR, 'char'), 'SNOW_GIS_DIR must be a string');

%% Load Soho map
[soho, cmap_gs, R] = geotiffread([SNOW_GIS_DIR, 'OSMap_Grayscale.tif']);

%% Initialization
min_c = 136;
max_c = 256;
cmap_gs(min_c:max_c, :) = colormap(jet(max_c - min_c + 1));

% Normalization for plotting
if nnz(x) > 0
    color  = (x - min(x))./(max(x) - min(x));
else
    color = x;
end
color = floor((max_c - min_c - 1).*color + min_c);

% Pixel coordinates of points in the map
[XG, YG] = world2pixel(G.coords(:,1), G.coords(:,2), R);

% Set edge color
G.plotting.edge_color = [0.1 0.1 0.1];

%% Display graph & signals
image(soho);
colormap(cmap_gs)
hold on

scatter3(XG, YG, color - min_c, 200, cmap_gs(color, :), '.');
if isfield(G, 'idx_pump')
    scatter3(XG(G.idx_pump), YG(G.idx_pump), color(G.idx_pump) - min_c, 100, 'r', 'o');
end

if G.Ne <= 1000
    G.coords = [XG, YG];
    gsp_plot_edges(G);
end

hold off
axis([min(XG) - 45, max(XG) + 45, min(YG) - 45, max(YG) + 45]);
axis off
h = colorbar;
set(h, 'Limits', [min_c, max_c])
set(h, 'Ticks', linspace(min_c, max_c, 6)')
set(h, 'TickLabels', linspace(min(x), max(x), 6)')

end

%% Auxiliary functions
function [Xp, Yp] = world2pixel(Xw, Yw, R)
%Converts world coordinates to pixel coordinates
    Xconv = sum(R.XIntrinsicLimits)./R.RasterExtentInWorldX;
    Xp = ( ( (Xw - R.XWorldLimits(1)) .* R.RasterExtentInWorldX ./ diff(R.XWorldLimits) ) - ...
        R.XIntrinsicLimits(1) )./Xconv + 1;
    Xp = floor(Xp);
    
    Yconv = sum(R.YIntrinsicLimits)./R.RasterExtentInWorldY;
    Yp = -( ( (Yw - R.YWorldLimits(1)) .* R.RasterExtentInWorldY ./ diff(R.YWorldLimits) ) - ...
        R.YIntrinsicLimits(1) )./Yconv + R.RasterExtentInWorldY;
    Yp = floor(Yp);
end