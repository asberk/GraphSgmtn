% % PredImgSgmtn - Predictive image segmentation using Diffuse Interface models on Graphs.
%    Ginzburg-Landau-type functionals with fidelity terms are used to group
%    objects in an image by patterns, where patterns are determined by
%    considering each pixel in the image as an NxN feature vector/matrix
%    and treating the image as fully-connected, undirected graph with
%    corresponding adjacency matrix. 
%
%    The segmentation is computed by solving the gradient descent equation 
%      (u^{n+1} - u^n)/dt 
%           = - \del{E_vex}{u}(u^{n+1}) + \del{E_cave}{u}(u^n),
%    where
%      E_vex(u) = \varepsilon/2 \int |\nabla u(x)|^2 \d x 
%               + c/2 \int |u(x)|^2 \d x
%
%      E_cave = -1/4\varepsilon \int(u(x)^2-1)^2 \d x 
%             + c/2 \int|u(x)|^2 \d x 
%             - \int \lambda(x)/2 (u(x) - u_0(x))^2 \d x
%
%    This file makes use of the Nyström Extension to approximate the
%    elements of the adjacency matrix, so that this script can be run even
%    for large images. 
%
% This script is distinct from PredUSHouseVoting.
% This script depends on the file(s):
% * gaussNorm.m
% * getFeatureArr.m
% * getImage.m
% * getNormalizedSampleWeights.m
% * getWeights.m
% * getWeightsXY.m
% * imgnmz.m
% * widen.m


%% Clear Workspace
clear all; close all; clc;

%% Import initial data and set computational parameter values
inputImg = getImage('img/tiger.gif');
szu0 = size(inputImg);

detectObject = 'tr'; % tiger, ti, tr
L = 500; % Number of random samples to take for X, where Z = X\cup Y
M = 100; % Number of iterations for convergence
N = 3; % neighbourhood radius
tau = .2; % weight scaling parameter
c = 1; % arbitrary constant to adjust region of convexity/concavity (how to determine this???)
dt = .1; % time stepping value
epsilon = 2; % interface width

u0 = zeros(size(inputImg));
switch detectObject
    case {'tiger', 'ti'}
        u0(30:52, 59:84) = 1; % face
        u0(53:69, 66:84) = 1; % neck
        u0(70:99, 66:146) = 1; % body
        u0(7:41, 104:147) = -1; % background
    case {'trees', 'tr'}
        % left tree
        u0(5:16, 21:33) = 1; % top-left of tree
        u0(41:55, 31:43) = 1; % bottom-right of tree
        % right tree
        u0(10:17, 95:109) = 1;
        % not right tree
        u0(6:15, 81:88) = -1;
        % tiger
        u0(73:95, 81:134) = -1; % tiger body
    otherwise
        error('whoops');
end
u0 = u0(:).'; % row vector

%% Get features, and choose L at random

FtrArr = getFeatureArr(inputImg, N); % Neumann, 2D matrix
szFtr = size(FtrArr);

% Get L unique pseudo-random integers
display('obtaining pseudo-random vector of integers...');
tic;
ranL = unique(randi(szFtr(2), [1, L]));
while size(unique(ranL), 2) < L 
    ranL = unique([ ranL, randi(szFtr(2), [1,L])]);
end
ranL = ranL(1:L); % trim extra
toc;

% FtrArr = Z = X \cup Y
ranLC = setdiff(1:szFtr(2), ranL); % the complement of ranL
X = FtrArr(:, ranL);
Y = FtrArr(:, ranLC);
u0 = u0([ranL ranLC]);

%% Compute normalized weight matrices Wxx and Wxy
% where, e.g., W_XY = [ w(x1, y1), ..., w(x1, y{n-L}) ] 
%                     [            ...                ]
%                     [ w(xL, y1), ..., w(xL, y{n-L]) ]

display('Approximate graph Laplacian using Nyström Extension.');
display('Computing partial adjacency matrix...')
[Wxx, Wxy] = getNormalizedSampleWeights(X,Y, @gaussNorm, tau);

%% More stuff
display('Computing eigenvectors of sampled subspace, components of orthognal approximation...');
tic;
[Bx, Gamma] = eig(Wxx);
sqrtGamma = Gamma^(1/2);
S = Bx*(sqrtGamma\Bx.');
Q = Wxx + S*(Wxy*Wxy.')*S;
[A Xi] = eig(Q);
toc;

%%
display('Approximating leading eigenvectors of graph Laplacian on orthogonal subspace...');
tic;
Phi = [Bx*sqrtGamma; Wxy.'*Bx/sqrtGamma]*Bx.'*(A/(Xi^(1/2)));
Lambda = 1- diag(Xi);
toc;

%% Visualize eigenvectors
% szU0 = size(unaught);
% for k = 1:15
%     imshow(imgnmz(reshape(Phi(:,k), szU0([2, 1]))).');
%     pause;
% end

%% Convex splitting for Graph Laplacian
display('Running convex splitting scheme for graph Laplacian...');

lenu0 = length(u0);
% initialize (all are 1-by-L vectors)
a = u0*Phi./lenu0; 
b = (u0.^3)*Phi./lenu0; 
d = zeros(1, L);
D = 1 + dt*(epsilon*Lambda.' + c);

eta = abs(u0); %fidelity term; known information

tic;
for m = 1:M
    a = ((1+dt/epsilon + c*dt)*a - dt/epsilon*b - dt*d)./D;
    u = a*Phi.';
    b = (u.^3)*Phi;
    d = (eta.*(u-u0))*Phi;

end
toc;

%% Post convergence rearrangement
% after we run the algorithm, we need to un-sort the elements
uM = zeros(size(u));
uM(ranL) = u(1:L);
uM(ranLC) = u(L+1:end);
uM = reshape(uM, [szu0(2), szu0(1)]).';

%% Plot initial data and final image
figure(1);
subplot(2,2,1);
imshow(inputImg);
switch detectObject
    case {'tiger', 'ti'}
        rectangle('Position', [30, 59, 52-30, 84-59]) = 1; % face
        rectangle('Position', [66, 53, 84-66, 69-53], 'EdgeColor', [1 0 0]); % neck
        rectangle('Position', [66, 70, 146-66, 99-70], 'EdgeColor', [1 0 0]); % body
        rectangle('Position', [104, 7, 147-104, 41-7], 'EdgeColor', [0 0 1]); % background
    case {'trees', 'tr'}
        % left tree
        rectangle('Position', [21, 5, 33-21, 16-5], 'EdgeColor', [1 0 0]); % top-left of tree
        rectangle('Position', [31, 41, 43-31, 55-41], 'EdgeColor', [1 0 0]); % bottom-right of tree
        % right tree
        rectangle('Position', [95, 10, 109-95, 17-10], 'EdgeColor', [1 0 0]); 
        % not right tree
        rectangle('Position', [81, 6, 88-81, 15-6], 'EdgeColor', [0 0 1]); 
        % tiger
        rectangle('Position', [81, 73, 134-81, 95-73], 'EdgeColor', [0 0 1]);  % tiger body
    otherwise
        error('whoops');
end
subplot(2,2,2);
imshow(imgnmz(uM));
subplot(2,2,3);
hist(uM(:), 20);

%% Threshold images
figure(2);
for j = 1:4
    subplot(2, 2, j);
    imagesc(uM > .4*j-1);
    colormap gray;
end
