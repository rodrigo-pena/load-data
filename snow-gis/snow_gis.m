function [G, x, b, M] = snow_gis(SNOW_GIS_DIR, downsample)
%SNOW_GIS creates a GSPBox-compatible graph from GIS gathered by John Snow
%about a cholera outbreak in Soho, London, in the 19th century
%
%   Usage:
%       [G, x, b, M] = snow_gis(SNOW_GIS_DIR)
%
%   Input:
%       SNOW_GIS_DIR    : A string specifying the directory where the GIS
%                       dataset is located (see reference below).
%                       (DEFAULT: '~/data/snow_gis/');
%       downsample      : (Optional) A number in [0, 1]. Percentage of 
%                         points in the road to keep.
%                         (DEFAULT: 0.1)
%
%   Output:
%       G   : A Matlab structure encoding graph information.
%       x   : A vector with non-zero entry at the location of the infected
%             water pump.
%       b   : A vector with whose entries represent the observed death 
%             count by cholera at each point.
%       M   : A diagonal matrix encoding the observation mask (nodes where 
%             the deaths were observed)
%
%   Example:
%       G = snow_gis();
%
%   See also: snow_gis_reduce.m
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
assert(isa(SNOW_GIS_DIR, 'char'), 'SNOW_GIS_DIR must be a string');

if nargin < 2 || isempty(downsample)
    downsample = 0.1;
end
assert(isnumeric(downsample), 'downsample must be a number');
assert((downsample >= 0) && (downsample <= 1), ...
    'downsample must be between 0 and 1');

%% Initialization
gsp_reset_seed(2); % Set the seed for random processes

%% Load Soho map
[soho, ~, R] = geotiffread([SNOW_GIS_DIR, 'OSMap.tif']);

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
road_mask = road_mask.*(rand(size(road_mask)) < downsample) | endp | branchp;

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

%% Assemble k-NN graph
coords = [coords_map; coords_chol; coords_pump];
nb_pts = nb_pts_map + nb_pts_chol + nb_pts_pump;

param = struct('type', 'knn', 'use_flann', 0, 'k', 6, 'center', 0, ...
    'rescale', 0, 'epsilon', 0.01, 'use_l1', 0, 'target_degree', 0, ...
    'symmetrize_type', 'average', 'light', 0);

[spi, spj, dist] = gsp_nn_distanz(coords', coords', param);
Dist = sparse(spi, spj, dist, nb_pts, nb_pts);
Dist = gsp_symmetrize(Dist, param.symmetrize_type);

% Treat differently the weighting of the edges on the road and the egdes
% from death points
sigma_roads = mean(dist)^2;
W_roads = sparse(spi, spj, double(exp(-dist.^2/sigma_roads)), nb_pts, nb_pts);
W_roads(1:(nb_pts + 1):end) = 0;  % Remove values in the main diagonal
W_roads = gsp_symmetrize(W_roads, param.symmetrize_type);

% Change sigma of the exponential kernel based on in the subgraph
% containing the death points (and the egdes in the frontier)
d_nr = tril(Dist);
d_nr = d_nr(nb_pts_map + 1:end, :);
sigma_non_roads = max(mean(d_nr(d_nr ~= 0).^2, 2));
W_non_roads = sparse(spi, spj, double(exp(-dist.^2/sigma_non_roads)), nb_pts, nb_pts);
W_non_roads(1:(nb_pts + 1):end) = 0;  % Remove values in the main diagonal
W_non_roads = gsp_symmetrize(W_non_roads, param.symmetrize_type);

% Merge the two graphs
G = struct('N', nb_pts, 'W', W_non_roads, 'coords', coords, ...
    'type', 'nearest neighbors', 'sigma', [sigma_roads, sigma_non_roads]);
G.W(1:nb_pts_map, 1:nb_pts_map) = W_roads(1:nb_pts_map, 1:nb_pts_map);
G = gsp_graph_default_parameters(G);

%% Retain only the biggest connected component
[G_components, node_indexes] = connected_subgraphs(G);
G = G_components{1};
nodes = node_indexes{1};
Dist = Dist(nodes, nodes);
G.Dist = Dist;
G.idx_road = find(nodes <= nb_pts_map);
G.idx_cholera = find((nodes > nb_pts_map)&(nodes <= nb_pts - nb_pts_pump));
G.idx_pump = find((nodes > nb_pts - nb_pts_pump));
G.data_source = 'snow_gis';

%% Generate signal vectors
% Infected pump
x = zeros(nb_pts, 1);
x(nb_pts_map + nb_pts_chol + 1) = 1;
x = x(nodes);

% Observed death count
b = zeros(nb_pts, 1);
b((nb_pts_map + 1):(nb_pts_map + nb_pts_chol)) = death_count;
b = b(nodes);

% Observation mask
M = double(b ~= 0);

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