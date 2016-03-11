function [G, B, M] = etex_graph(ETEX_DIR, release, k)
%ETEX_GRAPH creates a GSPBox-compatible graph from ETEX tracer release data
%
%   Usage:
%       [G, B, M] = etex_graph(ETEX_DIR, release, k)
%
%   Input:
%       ETEX_DIR    : A string specifying the directory where the ETEX
%                     dataset is located (see reference below).
%                     (DEFAULT: '~/data/etex/');
%       release     : (Optional) A number in {1, 2}, specifying with
%                     release data to use.
%                     (DEFAULT: 1)
%       k           : (Optional) A number in [2, 167]. Number of neighbors
%                     of each station when constructing the k-NN graph.
%                     (DEFAULT: 6)
%
%   Output:
%       G   : A Matlab structure encoding graph information.
%       B   : A G.N-by-30 matrix with the observed tracer concentration at
%             each station, at each of the experiment's time steps.
%       M   : A G.N-by-30 matrix encoding the observation mask (nodes where
%             a valid tracer concentration was observed).
%
%   Example:
%       [G, B, M] = etex_graph();
%
%   See also: plot_etex.m
%
%   Requires: GSPBox (https://lts2.epfl.ch/gsp/)
%
%   Reference: https://rem.jrc.ec.europa.eu/RemWeb/etex/
%
% Author: Rodrigo Pena (rodrigo.pena@epfl.ch)
% Date: 02 Mar 2016

%% Parse input
if nargin < 1 || isempty(ETEX_DIR)
    ETEX_DIR = '~/data/etex/';
end
assert(isa(ETEX_DIR, 'char'), 'ETEX_DIR must be a string');

if nargin < 2 || isempty(release)
    release = 1;
end
assert(isnumeric(release), 'release must be a number');
assert(release == 1 || release == 2, 'release must be equal to 1 or 2');

if nargin < 3 || isempty(k); k = 6; end
assert((k >= 2) || (k <= 167), 'k must be between 2 and 167');

%% Read text files
switch release
    case 1
        ETEX_DIR = [ETEX_DIR, 'release-1/'];
        
        % Quality flags file path
        quality_file_path = [ETEX_DIR, 'pmch.cod'];
        
        % PMCH concentration file path
        concentration_file_path = [ETEX_DIR, 'pmch.dat'];
        
    case 2
        ETEX_DIR = [ETEX_DIR, 'release-2/'];
        
        % Quality flags file path
        quality_file_path = [ETEX_DIR, 'pmcp.codv1.2'];
        
        % PMCH concentration file path
        concentration_file_path = [ETEX_DIR, 'pmcp.datv1.2'];
end

% Read quality flags file
fid = fopen(quality_file_path, 'r');
assert(fid ~= -1, ['File ', quality_file_path, ' could not be opened']);
header_lines = textscan(fid, '%s', 2, 'delimiter', '\n');
cols = textscan(header_lines{1}{2}, '%d');
num_cols = length(cols{1});
format_string = repmat('%f', 1, num_cols);
format_string = ['%d%s', format_string];
input_text = textscan(fid,format_string);
Quality_Flags = cell2mat(input_text(3:end));
fclose(fid);

% Read PMCH concentration file
fid = fopen(concentration_file_path, 'r');
assert(fid ~= -1, ['File ', concentration_file_path, ' could not be opened']);
header_lines = textscan(fid, '%s', 2, 'delimiter', '\n');
cols = textscan(header_lines{1}{2}, '%d');
num_cols = length(cols{1});
format_string = repmat('%f', 1, num_cols);
format_string = ['%d%s', format_string];
input_text = textscan(fid,format_string);
PMCH_Concentration = cell2mat(input_text(3:end));
fclose(fid);

% Read station list file
fid = fopen([ETEX_DIR, 'stationlist.950130'], 'r');
assert(fid ~= -1, 'File stationlist.950130 could not be opened');
header_lines = textscan(fid, '%s', 5, 'delimiter', '\n');
format_string = '%s%s%f%f%d%s%s';
input_text = textscan(fid,format_string, 'delimiter', ',');
station_labels = input_text{1}(1:168);
station_names = input_text{2}(1:168);
lat = input_text{3}(1:168); % Given in degree.minutes
long = input_text{4}(1:168); % Given in degree.minutes
alt = input_text{5}(1:168);
coords = [long, lat];
fclose(fid);

%% Assemble k-NN graph
param = struct('type', 'knn', 'use_flann', 0, 'k', k, 'center', 0, ...
    'rescale', 0, 'epsilon', 0.01, 'use_l1', 0, 'target_degree', 0, ...
    'symmetrize_type', 'average', 'light', 0);

G = gsp_nn_graph(coords, param);
G.station_labels = station_labels;
G.station_names = station_names;
G.alt = alt;
G.data_source = 'etex';
G.idx_release_site = find(strcmp(station_names, 'Rennes'));

%% Generate signal matrices
% Observation mask
M = double(Quality_Flags == 1 | Quality_Flags == 11 | ...
    Quality_Flags == 21 | Quality_Flags == 31);

% Observed PMCH concentration
B = M .* PMCH_Concentration;

end