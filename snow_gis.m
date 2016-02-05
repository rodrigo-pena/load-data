function [G, x, b] = snow_gis(SNOW_GIS_DIR)
%SNOW_GIS creates a GSPBox-compatible graph from GIS gathered by John Snow
%about a cholera outbreak in Soho, London, in the 19th century
%
%   Usage:
%       [G, x, b] = snow_gis(SNOW_GIS_DIR)
%
%   Input:
%       SNOW_GIS_DIR    : A string specifying the directory where the GIS
%                       dataset is located (see reference below).
%                       (DEFAULT: '~/data/snow_gis/');
%
%   Output:
%       G   : A Matlab structure encoding graph information.
%       x   : A vector with the locations of the original cholera deaths.
%       b   : A vector with the final accounts of cholera deaths.
%
%   Example:
%       G = snow_gis();
%
%   Requires: GSPBox (https://lts2.epfl.ch/gsp/)
%
%   Reference: http://blog.rtwilson.com/john-snows-famous-cholera-analysis-data-in-modern-gis-formats/
%
% Author: Rodrigo Pena (rodrigo.pena@epfl.ch)
% Date: 27 Jan 2016

%% Parse input
if nargin < 1 || isempty(SNOW_GIS_DIR)
    SNOW_GIS_DIR = '~/data/snow_gis/';
end
assert(isa(SNOW_GIS_DIR, 'char'), ...
    'SNOW_GIS_DIR must be a string');

%% Initialization
gsp_reset_seed(10); % Set the seed for random processes

%% Load Soho map
[soho, cmap, R] = geotiffread([SNOW_GIS_DIR, 'OSMap.tif']);
[soho_gs, cmap_gs] = geotiffread([SNOW_GIS_DIR, 'OSMap_Grayscale.tif']);

%% Compute a road mask
road_mask = (soho ~= 204);
crop_mask = false(size(road_mask));
crop_mask(143:885, 335:1064) = true;
road_mask = road_mask .* crop_mask;
road_mask = bwmorph(road_mask, 'close', Inf);
road_mask = bwmorph(road_mask, 'skel', Inf);
road_mask = bwmorph(road_mask, 'spur', 20);
road_mask = bwmorph(road_mask, 'clean');
endp = bwmorph(road_mask, 'endpoints');
branchp = bwmorph(road_mask, 'branchpoints');
road_mask = road_mask.*(rand(size(road_mask)) < 0.25) | endp | branchp;

%% Compute world coordinates and pick the points defined by the road mask
[X, Y] = pixel2world(1:R.RasterExtentInWorldX, 1:R.RasterExtentInWorldY, R); 

% Plane coordinates as complex numbers
coords_map = repmat(X, [R.RasterExtentInWorldY, 1]) + ...
    1i.*repmat(Y', [1, R.RasterExtentInWorldX]);
coords_map = coords_map .* road_mask;
coords_map = coords_map(:);
coords_map = coords_map(coords_map(:) ~= 0);
coords_map = [real(coords_map), imag(coords_map)];
nb_pts_map = size(coords_map, 1);

%% Load cholera deaths info
cholera_deaths = shaperead([SNOW_GIS_DIR, 'Cholera_Deaths.shp']);
coords_chol = [ [cholera_deaths.X]' , [cholera_deaths.Y]' ];
nb_pts_chol = length(cholera_deaths);
death_count = [cholera_deaths.Count]';

%% Load pump info
pumps = shaperead([SNOW_GIS_DIR, 'Pumps.shp']);
coords_pump = [ [pumps.X]' , [pumps.Y]' ];
nb_pts_pump = length(pumps);

%% Assemble graph
coords = [coords_map; coords_chol; coords_pump];
nb_pts = nb_pts_map + nb_pts_chol + nb_pts_pump;

% Treat differently the weighting of the edges on the road and the egdes
% from death points
param = struct('type', 'knn', 'use_flann', 0, 'k', 5);
G_standard = gsp_nn_graph(coords, param);

% Change sigma of the exponential kernel based on in the subgraph
% containing the death points (and the egdes in the frontier)
[spi, spj, dist] = gsp_nn_distanz(coords', coords', param);
param.sigma = mean(dist((spi > nb_pts_map) | (spj > nb_pts_map)).^2);

% Merge the two graphs
G = gsp_nn_graph(coords, param);
G.W(1:nb_pts_map, 1:nb_pts_map) = G_standard.W(1:nb_pts_map, 1:nb_pts_map);

G = gsp_graph_default_parameters(G);

%% Retain only the biggest connected component
[G_components, node_indexes] = connected_subgraphs(G);
G = G_components{1};
nodes = node_indexes{1};

%% Generate signal vectors
% Infected pump
x = zeros(nb_pts, 1);
x(nb_pts_map + nb_pts_chol + 1) = 1;
x = x(nodes);

% Observed death count
b = zeros(nb_pts, 1);
b((nb_pts_map + 1):(nb_pts_map + nb_pts_chol)) = death_count;
b = b(nodes);

%% Display graph & signals
% TODO: make function to plot the graph on the map.
[XG, YG] = world2pixel(G.coords(:,1), G.coords(:,2), R);

cmap_gs(217:232, :) = colormap(jet(16));

image(soho_gs);
colormap(cmap_gs)
hold on
scatter(XG, YG, 200, cmap_gs(15*x + 217, :), '.');
hold off

figure
image(soho_gs);
colormap(cmap_gs)
hold on
scatter(XG, YG, 200, cmap_gs(b + 217, :), '.');
hold off

nnz(b)
length(b(nodes>nb_pts_map+nb_pts_chol))

figure
mapshow(soho, cmap, R);
xlabel('(m)')
ylabel('(m)')
mapshow(cholera_deaths);

figure
mapshow(soho, cmap, R);
xlabel('(m)')
ylabel('(m)')
mapshow(pumps);

%% Cleanup
rng default

end

%% Auxiliary functions

function [Xw, Yw] = pixel2world(Xp, Yp, R)
%Converts pixel coordinates to world coordinates
    Xconv = sum(R.XIntrinsicLimits)./R.RasterExtentInWorldX;
    Xw = R.XWorldLimits(1) + (Xconv.*(Xp - 1) + R.XIntrinsicLimits(1)).*...
        diff(R.XWorldLimits)./R.RasterExtentInWorldX;
    
    Yconv = sum(R.YIntrinsicLimits)./R.RasterExtentInWorldY;
    Yw = R.YWorldLimits(1) + (Yconv.*(R.RasterExtentInWorldY - Yp) + R.YIntrinsicLimits(1)).*...
        diff(R.YWorldLimits)./R.RasterExtentInWorldY;
end

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