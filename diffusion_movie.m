function diffusion_movie(G, g, x, nframes, tstep, SAVEPATH, plot_fun)
%DIFFUSION_MOVIE creates a movie out of the diffusion of a signal x on a
%graph G. The diffusion is determined by the diffusion kernel g.
%
%   Usage:
%       diffusion_movie(G, g, x)
%
%   Input:
%       G       : A Matlab structure encoding graph information.
%                   G.N is the number of nodes in the graph
%                   G.W is the graph's weighted adjacency matrix
%       g       : Diffusion kernel. Can be a function handle acting on the
%                 eigenvalues of the graph Laplacian, or a symmetric diffusion
%                 matrix acting on x.
%       x       : A G.N-by-1 vector with the initial signal on the nodes of
%                 the graph G.
%       nframes : Number of frames to record.
%                 (DEFAULT: 125).
%       tstep   : Timestep between frame captures.
%                 (DEFAULT: 0.04).
%       SAVEPATH: (Optional) A string indicating the path (with filename)
%                 in which save the video.
%                 (DEFAULT: [pwd, '/diffusion_movie.avi']).
%       plot_fun: Function handle specifying how to plot a signal on a
%                 graph. It should take G as the first input, and x as the
%                 second input.
%                 (DEFAULT: @(G, x) gsp_plot_signal(G, x)).
%
%   Output:
%
%   Example:
%
%   Requires: GSPBox (https://lts2.epfl.ch/gsp/)
%
%   Reference:
%
% Author: Rodrigo Pena (rodrigo.pena@epfl.ch)
% Date: 22 Feb 2016

%% Parse input
% G
assert(isfield(G, 'W'), ...
    'G should have field ''W'' specifying the adjacency matrix');
assert(isfield(G, 'N'), ...
    'G should have field ''N'' specifying the number of nodes');
assert(sum(size(G.W) ~= [G.N, G.N]) == 0, ...
    'G.W should be a G.N-by-G.N matrix');

% g
if isa(g, 'function_handle')
    handle_flag = 1;
elseif isnumeric(g)
    assert(sum(size(g) ~= [G.N, G.N]) == 0, ...
        'If g is numeric, it should be a G.N-by-G.N matrix');
    handle_flag = 0;
else
    error(['g should be a function handle acting on the eigenvalues ', ...
        'of the graph Laplacian, or a symmetric diffusion matrix ',...
        'acting on x']);
end

% x
assert(sum(size(x) ~= [G.N, 1]) == 0, 'x should be a G.N-by-1 vector');

% nframes
if nargin < 4 || isempty(nframes); nframes = 125; end
assert(isnumeric(nframes) && sum(size(nframes) ~= 1) == 0, ...
    'nframes should be a number')

% tstep
if nargin < 5 || isempty(tstep); tstep = 0.04; end
assert(isnumeric(tstep) && sum(size(tstep) ~= 1) == 0, ...
    'tstep should be a number')

% SAVEPATH
if nargin < 6 || isempty(SAVEPATH)
    SAVEPATH = [pwd, '/diffusion_movie.avi'];
end
assert(isa(SAVEPATH, 'char'), 'SAVEPATH should be a string');

% plot_fun
if nargin < 7 || isempty(plot_fun)
    plot_fun = @(G, x) gsp_plot_signal(G, x);
end
assert(isa(plot_fun, 'function_handle'), ['plot_fun should be a ', ...
    'function handle on G and x']);


%% Initialization
writerObj = VideoWriter(SAVEPATH);
writerObj.FrameRate = 1./tstep;
open(writerObj);

figure; clf;
set(gca,'nextplot','replacechildren');
set(gcf,'color', 'white');
axis tight

%% Make movie
plot_fun(G, x);
frame = getframe(gcf);
writeVideo(writerObj,frame);

for i = 2:nframes
    
    % TODO: check if the following forward process is correct
    if handle_flag
        x = x + tstep .* gsp_filter(G, g, x);
    else
        x = x + tstep .* g * x;
    end

    plot_fun(G, x);
    
    frame = getframe(gcf);
    
    writeVideo(writerObj,frame);
    
end

%% Cleanup
close(writerObj);

end