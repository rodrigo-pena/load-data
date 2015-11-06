function [G, x, b] = digg_graph(DIGG_DIR, id, directed, pinitial, pend)
% Outputs a GSP Toolbox-compatible graph from the Digg 2009 dataset.
%   Inputs:
%       DIGG_DIR:   Directory where the mutual.mat, friend_date.mat, 
%                   user_id.mat, friend_id.mat, vote_date.mat, 
%                   voter_id.mat, and story_id.mat files are located.
%       id:         ID of the Digg story to consider as signal on the 
%                   graph. The ID is an integer between 1 (default) and 97.
%       directed:   0(default) or 1: assemble and undirected 
%                   graph, or a directed graph, respectively. 
%       pinitial:   A number between 0 and 1 encoding the percentage of
%                   timesteps for the story to consider when determining 
%                   the original voters for this story. If set to 0, only 
%                   the very first recorded vote is considered as the 
%                   original vote. The closer this number is to 1, more 
%                   votes are considered as original votes.
%       pend:       A number between 0 and 1 encoding the percentage of
%                   timesteps for the story to consider when determining 
%                   the end of the information diffusion. If set to 0, then 
%                   b = x, and we don't observe the diffusion. If set to 1, 
%                   all the recorded timeteps for the story are considered,
%                   and b encodes all the votes observed until the end of 
%                   the recording for such story.
%
%   Outputs:
%       G:          A GSP Toolbox-compatible graph assembled from
%                   user friendship information.
%       x:          A vector with ones on the indices corresponding to the
%                   first voters of the Digg story being considered.
%       b:          A vector with ones on the indices corresponding to the
%                   voters observed after the information has been diffused
%
%   Example:
%       home = getenv('HOME');
%       DIGG_DIR = strcat(home, 'data/digg');
%       [G, x, b] = digg_graph(DIGG_DIR, 1);
% 
% Author: Rodrigo Pena (rodrigo.pena@epfl.ch)
% Date: 30 Oct 2015

%% Process input
if nargin < 2 || isempty(id) || (id < 1) || (id > 97)
    id = 1;
end

if nargin < 3 || isempty(directed) || (directed ~= 1 && directed ~= 0)
    directed = 0;
end

if nargin < 4 || isempty(pinitial) || (pinitial > 1.0 && pinitial < 0.0)
    pinitial = 0.0;
end

if nargin < 5 || isempty(pend) || (pend > 1.0 && pend < 0.0)
    pend = 1.0;
end

%% Load pre-assembled .mat files
mutual = []; friend_date = []; user_id = []; friend_id = [];
load(strcat(DIGG_DIR, 'mutual.mat'));
load(strcat(DIGG_DIR, 'friend_date.mat'));
load(strcat(DIGG_DIR, 'user_id.mat'));
load(strcat(DIGG_DIR, 'friend_id.mat'));

vote_date = []; voter_id = []; story_id = [];
load(strcat(DIGG_DIR, 'vote_date.mat'));
load(strcat(DIGG_DIR, 'voter_id.mat'));
load(strcat(DIGG_DIR, 'story_id.mat'));

%% Keep info only for the story we want
vote_date = vote_date(story_id == id);
voter_id = voter_id(story_id == id);

%% Establish the start and end times for the information diffusion
t_first = min(vote_date);
t_last = max(vote_date);
t_start = floor(pinitial * (t_last - t_first) + t_first);
t_end = floor(pend * (t_last - t_first) + t_first);

initial_voters = voter_id(vote_date <= t_start);
voter_id = voter_id(vote_date <= t_end);

%% Keep friendship connections formed previous to t_end
mutual = mutual(friend_date <= t_end);
user_id = user_id(friend_date <= t_end);
friend_id = friend_id(friend_date <= t_end);

%% Decide on mutual friendships
if ~directed % Keep only mutual friendships
    user_id = user_id(mutual == 1);
    friend_id = friend_id(mutual == 1);
else % Account for all friendships
    user_id = [user_id; friend_id(mutual == 1)];
    friend_id = [friend_id; user_id(mutual == 1)];
end

%% Identify unique user IDs
users = unique([user_id; friend_id; voter_id]);
N = max(users);

%% Create sparse adjacency matrix for the friend network
spi = user_id;
spj = friend_id;
spv = ones(length(user_id), 1);
A = sparse(spi,spj,spv,N,N);
A = spdiags(zeros(N,1), 0, A); % Remove self-friendships
if ~directed
    A = max(A, A');
end

%% Create graph structure
G.N = N;
G.W = A;
G.coords = [randperm(N); randperm(N)]';
G.A = logical(A);
G.type = 'unknown';
G.directed = directed;
G.d = sum(G.W,2);
G.Ne = nnz(G.W)*(directed + 0.5*~directed);
G.L = diag(G.d) - G.W;
G.lap_type = 'combinatorial';

try
    G = gsp_graph_default_plotting_parameters(G);
catch
    warning('GSPBox not found. Could not set default plotting parameters.')
end

%% Create the signals on the graph
x = zeros(N, 1);
x(initial_voters) = 1;
b = zeros(N, 1);
b(voter_id) = 1;

end