function [G, measurements] = chernobyl(CHERNOBYL_DIR, data_type, ...
    data_file, k)
%CHERNOBYL creates a GSPBox-compatible graph from the Chernobyl incident
%data
%
%   Usage:
%       [G, measurements] = chernobyl(CHERNOBYL_DIR, data_type, ...
%       data_file, k)
%
%   Input:
%       CHERNOBYL_DIR   : A string specifying the directory where the
%                         CHERNOBYL dataset is located (see reference
%                         below).
%                         (DEFAULT: '~/data/chernobyl/');
%       data_type       : A string specifying which type of data to read.
%           'air_concentration' : Read air concentration data.
%           'deposition'        : Read deposition data.
%           (DEFAULT: 'deposition')
%       data_file       : A string specifying which data file to read.
%           If data_type = 'air_concentration':
%               'CHERNAIR.TXT'  : Concentration of I-131, Cs-134 and Cs-137
%                                 (aerosol particles)
%               'CHERNIOD.TXT'  : Concentration of total I-131 (gas +
%                                 aerosol particles)
%           If data_type = 'deposition':
%               'CUMDEP.DAT'    : Cumulative deposition measurements of
%                                 Cs-134 and Cs-137 (across Europe).
%               'NORWAY.UPD'    : Cumulative deposition measurements in
%                                 Norway.
%               'POLAND.UPD'    : Cumulative deposition measurements in
%                                 Poland.
%               'RUMANIA.UPD'   : Cumulative deposition measurements in
%                                 Romania.
%           (DEFAULT: 'CUMDEP.DAT', if data_type = 'deposition';
%           'CHERNAIR.TXT', otherwise)
%       k               : Number of neighbors to use when assembling the
%                         k-NN graph.
%
%   Output:
%       G               : A Matlab structure encoding graph information.
%       measurements    : A Matlab structure encoding measurement
%                         information.
%
%   Example:
%       [G, measurements] = chernobyl(CHERNOBYL_DIR, data_type, ...
%       data_file, k);
%
%   See also: plot_chernobyl.m
%
%   Requires: GSPBox (https://lts2.epfl.ch/gsp/)
%
%   Reference: https://rem.jrc.ec.europa.eu/RemWeb/Browse.aspx?path=Chernobyl%20Data
%
%   Remark: Your Chernobyl data folder should consist of two directories,
%           namely, "air_concentration/" and "deposition/". The directory
%           air_concentration/ should contain files CHERNAIR.TXT and
%           CHERNIOD.TXT. The directory deposition/ should contain files
%           CUMDEP.DAT, NORWAY.UPD, POLAND.UPD, and
%           RUMANIA.UPD.
%
% Author: Rodrigo Pena (rodrigo.pena@epfl.ch)
% Date: 11 Mar 2016

%% Parse input
% CHERNOBYL_DIR
if nargin < 1 || isempty(CHERNOBYL_DIR)
    CHERNOBYL_DIR = '~/data/chernobyl/';
end
assert(isa(CHERNOBYL_DIR, 'char'), 'CHERNOBYL_DIR must be a string');

% data_type
if nargin < 2 || isempty(data_type); data_type = 'deposition'; end
assert(isa(data_type, 'char'), 'data_type must be a number');
assert(strcmp(data_type, 'air_concentration') || ...
    strcmp(data_type, 'deposition'), ...
    'data_type must be ''air_concentration'' or ''deposition''.');

% data_file
switch data_type
    case 'deposition'
        if nargin < 3 || isempty(data_file); data_file = 'CUMDEP.DAT'; end
        assert(strcmp(data_file, 'CUMDEP.DAT') || ...
            strcmp(data_file, 'NORWAY.UPD') || ...
            strcmp(data_file, 'POLAND.UPD') || ...
            strcmp(data_file, 'RUMANIA.UPD'), ...
            ['data_file must be', ' ''CUMDEP.DAT''', ...
            ' ''NORWAY.UPD''' , ' ''POLAND.UPD''', ...
            ' or', ' ''RUMANIA.UPD''', ...
            ' when data_type = ''deposition''.']);
        
    case 'air_concentration'
        if nargin < 3 || isempty(data_file); data_file = 'CHERNAIR.TXT'; end
        assert(strcmp(data_file, 'CHERNAIR.TXT') || ...
            strcmp(data_file, 'CHERNIOD.TXT'), ...
            ['data_file must be', ' ''CHERNAIR.TXT''', ...
            ' or', ' ''CHERNIOD.TXT''', ...
            ' when data_type = ''air_concentration''.']);
end

% k
if nargin < 4 || isempty(k); k = 6; end
assert(isnumeric(k) && length(k) == 1, 'k must be a number');

%% Read data file
file_path = [CHERNOBYL_DIR, data_type, '/', data_file];
fid = fopen(file_path, 'r');
assert(fid ~= -1, ['File ', file_path, ' could not be opened']);

% Establish the format_spec
switch data_file
    case 'CHERNAIR.TXT'
        format_spec = '%4s%33s%9s%11s%9s%6s%9s%5s%9s%5s%9s%5s%[^\n\r]';
        num_cols = [3, 4, 9, 11, 13];
        
    case 'CHERNIOD.TXT'
        format_spec = '%4s%33s%9s%11s%9s%6s%9s%4s%10s%[^\n\r]';
        num_cols = [3, 4, 9];
        
    otherwise
        %TODO: Norway, Poland, and Romania's data are not in CUMDEP.DAT.
        %Restric this function to read only cumulative deposition
        %measurements, and create a graph with the data from all the
        %deposition files.
        format_spec = '%21s%9s%9s%8s%4s%6s%5s%8s%4s%5s%[^\n\r]';
        num_cols = [2, 3, 6, 8];
end

% Read columns of data according to format string.
data_array = textscan(fid, format_spec, 'Delimiter', '', ...
    'WhiteSpace', '',  'ReturnOnError', false);

% Remove white space around all cell columns.
for i = 1:length(data_array); data_array{i} = strtrim(data_array{i}); end

fclose(fid);

%% Process numeric fields of the data array
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
coords = cell2mat(raw(:, num_cols(2:-1:1)));
raw(isnan(sum(coords,2)), :) = [];
coords(isnan(sum(coords,2)), :) = [];

%% Assemble k-NN graph
param = struct('type', 'knn', 'use_flann', 0, 'k', k, 'center', 0, ...
    'rescale', 0, 'epsilon', 0.01, 'use_l1', 0, 'target_degree', 0, ...
    'symmetrize_type', 'average', 'light', 0);
G = gsp_nn_graph(coords, param);

%% Record measurement information
switch data_file
    case 'CHERNAIR.TXT'
        %TODO: Data has repeated nodes. Convert repeated node information 
        %into a time series of measurements!
        G.country_code = raw(:, 1);
        G.locality_name = raw(:, 2);
        measurements.date = raw(:, 5);
        measurements.hour_end_sampling = raw(:, 6);
        measurements.duration = raw(:, 7);
        measurements.I131_concentration = cell2mat(raw(:, 9));
        measurements.Cs134_concentration = cell2mat(raw(:, 11));
        measurements.Cs137_concentration = cell2mat(raw(:, 13));
        measurements.unit = 'Bq/m3';
        
    case 'CHERNIOD.TXT'
        %TODO: Data has repeated nodes. Convert repeated node information 
        %into a time series of measurements!
        G.country_code = raw(:, 1);
        G.locality_name = raw(:, 2);
        measurements.date = raw(:, 5);
        measurements.hour_end_sampling = raw(:, 6);
        measurements.duration = raw(:, 7);
        measurements.I131_concentration = cell2mat(raw(:, 9));
        measurements.unit = 'Bq/m3';
        
    otherwise
        G.locality_name = raw(:, 1);
        measurements.date = raw(:, 4);
        measurements.Cs134_fallout = cell2mat(raw(:, 6));
        measurements.Cs137_fallout = cell2mat(raw(:, 8));
        measurements.thickness = raw(:, 9);
        measurements.sample_type = raw(:, 10);
        measurements.reference = raw(:, 11);
        measurements.unit = 'kBq/m2';
end

end