function diffusion_movie(G, x, a, nframes, framerate, SAVEPATH, plot_fun)
%DIFFUSION_MOVIE creates a movie out of the heat diffusion of a signal x on
%a graph G. The diffusion is parametrized by the diffusion time t.
%
%   Usage:
%       diffusion_movie(G, x, t, nframes, framerate, SAVEPATH, plot_fun)
%
%   Input:
%       G           : A Matlab structure encoding graph information.
%                       G.N is the number of nodes in the graph
%                       G.W is the graph's weighted adjacency matrix
%                       G.L is the graph Laplacian
%       a           : A number representing the thermal diffusivity, i.e., 
%                     the a in dx/dt = -a * G.L * x.
%                     (DEFAULT: 1./normest(G.L)).
%       x           : A G.N-by-1 vector with the initial signal on the 
%                     nodes of the graph G.
%       nframes     : Number of frames to record.
%                     (DEFAULT: 125).
%       framerate   : Number of frames per second.
%                     (DEFAULT: 25).
%       SAVEPATH    : (Optional) A string indicating the path (with 
%                     filename) in which save the video.
%                     (DEFAULT: [pwd, '/diffusion_movie.avi']).
%       plot_fun    : Function handle specifying how to plot a signal on a
%                     graph. It should take G as the first input, and x as 
%                     the second input.
%                     (DEFAULT: @(G, x) gsp_plot_signal(G, x)).
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
assert(isfield(G, 'N'), ...
    'G should have field ''N'' specifying the number of nodes');
assert(isfield(G, 'W'), ...
    'G should have field ''W'' specifying the adjacency matrix');
assert(isfield(G, 'W'), ...
    'G should have field ''L'' specifying the graph Laplacian');
assert(sum(size(G.W) ~= [G.N, G.N]) == 0, ...
    'G.W should be a G.N-by-G.N matrix');
assert(sum(size(G.L) ~= size(G.W)) == 0, ...
    'G.L should be the same size as G.W');

% x
assert(sum(size(x) ~= [G.N, 1]) == 0, 'x should be a G.N-by-1 vector');

% a
if nargin < 3 || isempty(a); a = 1./normest(G.L); end
assert(isnumeric(a) && sum(size(a) ~= 1) == 0, 'a should be a number');


% nframes
if nargin < 4 || isempty(nframes); nframes = 125; end
assert(isnumeric(nframes) && sum(size(nframes) ~= 1) == 0, ...
    'nframes should be a number');

% framerate
if nargin < 5 || isempty(framerate); framerate = 25; end
assert(isnumeric(framerate) && sum(size(framerate) ~= 1) == 0, ...
    'framerate should be a number')

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
    'function handle on G and x (in this order)']);


%% Initialization
writerObj = VideoWriter(SAVEPATH);
writerObj.FrameRate = framerate;
writerObj.Quality = 100;
open(writerObj);

fig_handle = figure;
set(fig_handle, 'Position', [100, 100, 1080, 720]);
set(gca,'nextplot','replacechildren');
set(gcf,'color', 'white');
axis tight

%% Make movie
plot_fun(G, x);
frame = getframe(gcf);
writeVideo(writerObj,frame);

for n = 2:nframes
    
    x = x - a * (G.L * x);
    
    plot_fun(G, x);
    drawnow;
    
    frame = getframe(gcf);
    
    writeVideo(writerObj,frame);
    
end

%% Cleanup
close(writerObj);

end