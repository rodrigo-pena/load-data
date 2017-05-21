function load_data_start()
%SRC_LOC_START Initialize the toolbox
%   Details:
%       Initialization script for the load-data toolbox. This
%       script adds to the path the folders with the files needed to run
%       the toolbox.
%
% Author: Rodrigo Pena (rodrigo.pena@epfl.ch)
% Date: 18 March 2016

%% Add dependencies
global GLOBAL_loadpath;
GLOBAL_loadpath = fileparts(mfilename('fullpath'));
addpath([GLOBAL_loadpath, ':', ...
    GLOBAL_loadpath, '/chernobyl:', ...
    GLOBAL_loadpath, '/etex:', ...
    GLOBAL_loadpath, '/snow-gis:', ...
    GLOBAL_loadpath, '/spectral-layout:']);

end
