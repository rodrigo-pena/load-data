% Tonnetz as a torus

G = gsp_torus(12,4);
pts = G.coords;
pts_norm = sqrt(sum(pts(:,1:2).^2, 2));
m1 = min(pts_norm); m2 = max(pts_norm);
pts(:,1:2) = pts(:,1:2) .* repmat(1./pts_norm, [1, 2]);
pts_norm = 1 + (1./(m2 - m1)).*(pts_norm - m1);
pts(:,1:2) = pts(:,1:2) .* repmat(pts_norm, [1, 2]);
pts(:,3) = pts(:,3)./(range(pts(:,3)));

cross_mat = @(u) [0 -u(3), u(2); u(3), 0, -u(1); -u(2), u(1), 0];
R = @(n, t) cos(t).*eye(3) + sin(t).*cross_mat(n) + (1 - cos(t)).*kron(n, n');
for i = 2:12
    n = cross(pts(i,:), pts(i+24,:));
    n = n./norm(n);
    T = (pts(i+12,:) + pts(i+36,:))./2;
    pts(i:12:end, :) = pts(i:12:end, :) - repmat(T, [4, 1]);
    pts(i:12:end, :) = (R(n, -(i-1)*pi/6)*pts(i:12:end, :)')';
    pts(i:12:end, :) = pts(i:12:end, :) + repmat(T, [4, 1]);
end
G.coords = pts;
G.plotting = [];
G = gsp_graph_default_plotting_parameters(G);

% Connect the 6 nearest neighbors
param = struct('type', 'knn', 'use_flann', 0, 'k', 6);
G = gsp_nn_graph(G.coords, param);
figure;
daspect([1 1 1])
gsp_plot_graph(G);

C1 = {'C';'G';'D';'A';'E';'B';'F#/Gb';'C#/Db';'G#/Ab';'D#/Eb';'A#/Bb';'F'};
C2 = circshift(C1, 5);
C3 = circshift(C2, 5);
C4 = circshift(C3, 5);
C = vertcat(C1, C2, C3, C4);
disp = 0.05; % Text displacement parameter
text(pts(:,1) + disp, pts(:,2) + disp, pts(:,3) + disp, C);
