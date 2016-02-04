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

%% Load Soho map
[soho, cmap, R] = geotiffread([SNOW_GIS_DIR, 'OSMap.tif']);
[soho_gs, cmap_gs] = geotiffread([SNOW_GIS_DIR, 'OSMap_Grayscale.tif']);
% info = geotiffinfo([SNOW_GIS_DIR, 'OSMap.tif']);
% mapshow(soho, cmap, R);
% xlabel('(m)')
% ylabel('(m)')

%% Compute world coordinates of the pixels on the map
[X, Y] = pixel2world(1:R.RasterExtentInWorldX, 1:R.RasterExtentInWorldY, R); 

% Plane coordinates as complex numbers
coords_map = repmat(X, [R.RasterExtentInWorldY, 1]) + ...
    1i.*repmat(Y', [1, R.RasterExtentInWorldX]);

%% Get road mask & apply it to world coordinates
road_mask = (soho ~= 204);
road_mask = bwmorph(road_mask, 'close', Inf);
road_mask = bwmorph(road_mask, 'skel', Inf);
road_mask = bwmorph(road_mask, 'spur', 20);
road_mask = bwmorph(road_mask, 'clean');

% figure; 
% image(soho);
% colormap(cmap)
% hold on
% h = image(3*double(road_mask));
% set(h, 'AlphaData', road_mask);
% hold off

coords_map = coords_map .* road_mask;
coords_map = coords_map(:);
coords_map = coords_map(coords_map(:) ~= 0);
coords_map = [real(coords_map), imag(coords_map)];
nb_pts_map = size(coords_map, 1);

%% Load cholera deaths info
cholera_deaths = shaperead([SNOW_GIS_DIR, 'Cholera_Deaths.shp']);
% mapshow(cholera_deaths);

coords_chol = [ [cholera_deaths.X]' , [cholera_deaths.Y]' ];
nb_pts_chol = length(cholera_deaths);
death_count = [cholera_deaths.Count]';

%% Load pump info
pumps = shaperead([SNOW_GIS_DIR, 'Pumps.shp']);
% mapshow(pumps);

coords_pump = [ [pumps.X]' , [pumps.Y]' ];
nb_pts_pump = length(pumps);

%% Assemble graph
coords = [coords_map; coords_chol; coords_pump];
param = struct('type', 'knn', 'use_flann', 1, 'k', 10);
G = gsp_nn_graph(coords, param);
W = 0.*G.W;
W(G.W > 1e-1) = G.W(G.W > 1e-1);
G.W = W;
G = gsp_graph_default_parameters(G);

%% Generate signal vectors
% Infected pump
x = zeros(G.N, 1);
x(nb_pts_map + nb_pts_chol + 1) = 1; 

% Observed death count
b = zeros(G.N, 1);
b((nb_pts_map + 1):(nb_pts_map + nb_pts_chol)) = death_count;

%% Plot signals on the graph
[XG, YG] = world2pixel(G.coords(:,1), G.coords(:,2), R);

cmap_gs(217:232, :) = colormap(jet(16));

image(soho_gs);
colormap(cmap_gs)
hold on
scatter(XG, YG, 200, cmap_gs(b + 217, :), '.');
hold off

end

function [Xw, Yw] = pixel2world(Xp, Yp, R)
    Xconv = sum(R.XIntrinsicLimits)./R.RasterExtentInWorldX;
    Xw = R.XWorldLimits(1) + (Xconv.*(Xp - 1) + R.XIntrinsicLimits(1)).*...
        diff(R.XWorldLimits)./R.RasterExtentInWorldX;
    
    Yconv = sum(R.YIntrinsicLimits)./R.RasterExtentInWorldY;
    Yw = R.YWorldLimits(1) + (Yconv.*(R.RasterExtentInWorldY - Yp) + R.YIntrinsicLimits(1)).*...
        diff(R.YWorldLimits)./R.RasterExtentInWorldY;
end

function [Xp, Yp] = world2pixel(Xw, Yw, R)
    Xconv = sum(R.XIntrinsicLimits)./R.RasterExtentInWorldX;
    Xp = ( ( (Xw - R.XWorldLimits(1)) .* R.RasterExtentInWorldX ./ diff(R.XWorldLimits) ) - ...
        R.XIntrinsicLimits(1) )./Xconv + 1;
    Xp = floor(Xp);
    
    Yconv = sum(R.YIntrinsicLimits)./R.RasterExtentInWorldY;
    Yp = -( ( (Yw - R.YWorldLimits(1)) .* R.RasterExtentInWorldY ./ diff(R.YWorldLimits) ) - ...
        R.YIntrinsicLimits(1) )./Yconv + R.RasterExtentInWorldY;
    Yp = floor(Yp);
end