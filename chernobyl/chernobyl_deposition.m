function [G, measurements] = chernobyl_deposition(CHERNOBYL_DIR, k)
%CHERNOBYL creates a GSPBox-compatible graph from the Chernobyl incident
%Cs deposition data
%
%   Usage:
%       [G, measurements] = chernobyl(CHERNOBYL_DIR, k)
%
%   Input:
%       CHERNOBYL_DIR   : A string specifying the directory where the
%                         CHERNOBYL cumulative deposition data are located
%                         (see reference below).
%                         (DEFAULT: '~/data/chernobyl/deposition/');
%       k               : Number of neighbors to use when assembling the
%                         k-NN graph.
%
%   Output:
%       G               : A Matlab structure encoding graph information.
%       measurements    : A Matlab structure encoding measurement
%                         information.
%
%   Example:
%       [G, measurements] = chernobyl(CHERNOBYL_DIR, k);
%
%   See also: plot_chernobyl.m
%
%   Requires: GSPBox (https://lts2.epfl.ch/gsp/)
%
%   Reference: https://rem.jrc.ec.europa.eu/RemWeb/Browse.aspx?path=Chernobyl%20Data
%
%   Remark: Your Chernobyl deposition data folder should contain the files
%           CUMDEP.DAT, NORWAY.UPD, POLAND.UPD, and RUMANIA.UPD.
%
% Author: Rodrigo Pena (rodrigo.pena@epfl.ch)
% Date: 11 Mar 2016

%% Parse input
% CHERNOBYL_DIR
if nargin < 1 || isempty(CHERNOBYL_DIR)
    CHERNOBYL_DIR = '~/data/chernobyl/deposition/';
end
assert(isa(CHERNOBYL_DIR, 'char'), 'CHERNOBYL_DIR must be a string');

% k
if nargin < 4 || isempty(k); k = 6; end
assert(isnumeric(k) && length(k) == 1, 'k must be a number');

%% Initialization
data_filepaths = cell(1);
data_filepaths{1,1} = [CHERNOBYL_DIR, 'CUMDEP.DAT'];
data_filepaths{2,1} = [CHERNOBYL_DIR, 'NORWAY.UPD'];
data_filepaths{3,1} = [CHERNOBYL_DIR, 'POLAND.UPD'];
data_filepaths{4,1} = [CHERNOBYL_DIR, 'RUMANIA.UPD'];

ptr_i = 1;
data = cell(4514, 11);

%% Process data files
for i = 1:length(data_filepaths)
    % Open file
    fid = fopen(data_filepaths{i}, 'r');
    assert(fid ~= -1, ...
        ['File ', data_filepaths{i}, ' could not be opened']);
    
    % Establish the format_spec
    format_spec = '%21s%9s%9s%8s%4s%6s%5s%8s%4s%5s%[^\n\r]';
    num_cols = [2, 3, 6, 8];
    
    % Read columns of data according to format string.
    data_array = textscan(fid, format_spec, 'Delimiter', '', ...
        'WhiteSpace', '',  'ReturnOnError', false);
    
    % Remove white space around all cell columns.
    for j = 1:length(data_array); 
        data_array{j} = strtrim(data_array{j}); 
    end
    
    % Close file
    fclose(fid);
    
    % Reshape data_array
    raw = repmat({''}, length(data_array{1}), length(data_array));
    for col = 1:length(data_array)
        raw(1:length(data_array{col}), col) = data_array{col};
    end
    
    % Convert strings in the numerical columns to numbers.
    for col = num_cols
        raw_data = data_array{col};
        
        for row = 1:size(raw_data, 1);
            % Create a regular expression to detect and remove non-numeric
            % prefixes and suffixes.
            regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
            try
                result = regexp(raw_data{row}, regexstr, 'names');
                numbers = result.numbers;
                
                % Detected commas in non-thousand locations.
                invalid_thousands_separator = false;
                if any(numbers==',');
                    thousands_reg_exp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                    if isempty(regexp(thousands_reg_exp, ',', 'once'));
                        numbers = NaN;
                        invalid_thousands_separator = true;
                    end
                end
                
                % Convert numeric strings to numbers.
                if ~invalid_thousands_separator;
                    numbers = textscan(strrep(numbers, ',', ''), '%f');
                    raw{row, col} = numbers{1};
                end
                
            catch
            end
        end
    end
    
    % Replace non-numeric strings in these numeric columns with NaN.
    raw_numeric_columns = raw(:, num_cols);
    R = cellfun(@(x) isempty(x) || (ischar(x) && all(x==' ')), ...
        raw_numeric_columns);
    raw_numeric_columns(R) = {NaN};
    raw(:, num_cols) = raw_numeric_columns;
    
    % Remove entries for which we don't have the geographical coordinates
    geo_coords = cell2mat(raw(:, num_cols(2:-1:1)));
    raw(isnan(sum(geo_coords,2)), :) = [];
    
    % Concatenate data
    [r, ~] = size(raw);
    data(ptr_i:(ptr_i + r - 1), :) = raw;
    ptr_i = ptr_i + r;
    
end
%% Assemble k-NN graph
coords = cell2mat(data(:, num_cols(2:-1:1)));
param = struct('type', 'knn', 'use_flann', 0, 'k', k, 'center', 0, ...
    'rescale', 0, 'epsilon', 0.01, 'use_l1', 0, 'target_degree', 0, ...
    'symmetrize_type', 'average', 'light', 0);
G = gsp_nn_graph(coords, param);

%% Record measurement information
G.locality_name = data(:, 1);
measurements.date = data(:, 4);
measurements.Cs134_fallout = cell2mat(data(:, 6));
measurements.Cs137_fallout = cell2mat(data(:, 8));
measurements.thickness = data(:, 9);
measurements.sample_type = data(:, 10);
measurements.reference = data(:, 11);
measurements.unit = 'kBq/m2';

end