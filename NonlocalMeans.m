% % Nonlocal Means - uses nonlocal-means to construct a fully-connected
%     undirected graph, whose vertices are pixels, which are each 
%     represented as NxN neighbourhoods of pixels about that pixel. The
%     eigenvectors of the symmetric graph Laplacian are determined, and the
%     top 10 visualized. 

%% Clear workspace
clear all; close all; clc;

%%
A = getImage('img/BaboonEye.jpg');
szA = size(A);

tooBig = prod(szA) > 5000;

N = 5; % nbhd size


%% compute "feature vectors"
FtrArr = getFeatureArr(A, N); % Neumann by default; LoD by default.

%% Compute Laplacian
W = getWeights(FtrArr, @gaussNorm, 1, 1); % weight matrix; computationally intensive
clear FtrArr N;

%% Temporary interruption
if tooBig
    warning('I hope you have enough RAM to compute this');
    GL = getGraphLaplacian(W); % use default settings to return symmetric graph Laplacian
else
    [GL, D] = getGraphLaplacian(W); % use default settings to return symmetric graph Laplacian
end
clear W;

%% Get EigenInformation of Graph Laplacian
[V,Lam] = eig(GL);

%% Plot the top 10 eigenvectors
for j = 1:10
    imshow(imgnmz(reshape(V(:,j), szA)).');
    pause;
end