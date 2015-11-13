%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test spectral layout functions 
%
% Requires: - GSPBox (https://lts2.epfl.ch/gsp/)
%           - MatlabBGL (http://dgleich.github.io/matlab-bgl/)
% 
% Author: Rodrigo Pena
% E-mail: rodrigo.pena@epfl.ch
% Date: 13 Nov 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup
addpath('../');
method = 'standard';
m = 100;
p = 4;
indices = 2:p;
TOL = 1e-8;
MAX_ITER = 1000;

%% Create graph
G = gsp_logo();

%% Test spectral_layout.m
tic;
G1 = spectral_layout(G, indices, method);
t1 = toc;
fprintf('Time to execute spectral_layout.m: %1.2f\n', t1);

%% Test spectral_layout_approx.m
tic;
G2 = spectral_layout_approx(G, p, m, method, TOL, MAX_ITER);
t2 = toc;
fprintf('Time to execute spectral_layout.m: %1.2f\n', t2);

%% Display figures
figure();
subplot(221)
gsp_plot_graph(G);
title('Original Coordinates')

subplot(222)
gsp_plot_graph(G1);
title('Spectral Layout')

subplot(223)
gsp_plot_graph(G2);
title('Spectral Layout Approx.')

%% Cleanup
rmpath('../');