% % PredWvltSgmtn - Predictive image segmentation using Diffuse Interface models on Graphs.
%    Ginzburg-Landau-type functionals with fidelity terms are used to group
%    objects in an image by patterns, where patterns are determined by
%    considering each pixel in the image as an NxN feature vector/matrix
%    and treating the image as fully-connected, undirected graph with
%    corresponding adjacency matrix, the weights of which are computed using
%    the wavelet coefficients of the blocks. 
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
runDebug = 0;
runBlur = 0;
if runBlur
    blurType = 'gaussian';
    blurParms{1} = 5; % HSIZE
    blurParms{2} = 1; % SIGMA
end

%% Import initial data 
%imgstr = '~/Downloads/cows/cows1.png';
imgstr = 'img/tiger.gif';
inputImg = getImage(imgstr);
szu0 = size(inputImg);
if strcmp(imgstr, 'img/tiger.gif');
    % set detect object to RHS
    detectObject = 'tr'; % tiger, ti, tr
else
    % a different image was chosen
    detectObject = 'different';
end

if runBlur
    psf = fspecial(blurType, blurParms{:});
    inputImg = imfilter(inputImg, psf, 'symmetric', 'conv');
    
end


%% Set computational parameter values and initial formatting
L = 100; % Number of random samples to take for X, where Z = X\cup Y
M = 100; % Number of iterations for convergence
N = 3; % neighbourhood radius
tau = 2; % weight scaling parameter
c = .25; % arbitrary constant to adjust region of convexity/concavity (how to determine this???)
dt = .1; % time stepping value
epsilon = 2; % interface width
wname = 'db1'; % haar

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
    case {'different', 'diff'}
        message = sprintf('Manually selecting SSL regions.\nLeft click and hold to begin drawing.\nSimply lift the mouse button to finish');
        message1 = sprintf('Select portion of foreground region.');
        message2 = sprintf('Select protion of background region.');
        figure(1);
        imshow(inputImg, []);
        axis on;
        set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
        title('Select portion of foreground region', 'FontSize', 16);
        uiwait(msgbox(message));
        %uiwait(msgbox(message1));
        hFHfg = imfreehand(gca, 'Closed', 1);
        binaryImageObject = hFHfg.createMask();
        xyObject = hFHfg.getPosition;
        %axis on;
        %set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
        title('Select portion of background region', 'FontSize', 16);
        %uiwait(msgbox(message2));
        hFHbg = imfreehand(gca, 'Closed', 1);
        binaryImageBG = hFHbg.createMask();
        xyBG = hFHbg.getPosition;
        close 1;
        
        u0 = binaryImageObject - binaryImageBG;
        clear message message1 message2 hFHfg hFHbg
    otherwise
        error('whoops');
end
u0 = u0(:).'; % row vector

%% Get features, and choose L at random
myFFT = @(x) abs(fftshift(fft2(x)));
FtrArr = getFeatureArr(inputImg, N, myFFT); % Neumann, 2D matrix
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

%% Nyström extension
[Phi, Lambda] = nystrom(Wxx, Wxy);

%% Visualize eigenvectors
if runDebug
figure(1);
sampledPoints = zeros(size(inputImg));
sampledPoints(ranL) = 1;
rChannel = min(1, inputImg + sampledPoints);
gbChannel = max(0, inputImg - sampledPoints);

[Y I] = sort(Lambda, 'descend');
for k = 0:floor(size(Phi, 2))
    Jvec = (2:16) + 15*k;
    Jvec = Jvec(Jvec <= size(Phi,2));
    subplot(4,4,1);
    imshow(cat(3, rChannel, gbChannel, gbChannel), []);
    title(['$j = ', num2str(1 + 15*k), '\ldots', num2str(15*(k+1)), '$'], 'Interpreter', 'latex');
for j = Jvec
    subplot(4,4,j-15*k);
    imshow(reOrder(Phi(:,I(j-1)), szu0, ranL), []);
    title(['$\lambda = ', num2str(Y(j-1), 3), '$'], 'Interpreter', 'latex');
    colorbar;
end
pause;
end
pause;
close;
end

%% Convex splitting for Graph Laplacian
display('Running convex splitting scheme for graph Laplacian...');

Parms.c = c;
Parms.dt = dt;
Parms.epsilon = epsilon;
Parms.eta = abs(u0);
Parms.M = M;
Parms.L = L;

u = pdeGraphSgmtn(u0, Phi, Lambda, Parms);

%% Post convergence rearrangement
% after we run the algorithm, we need to un-sort the elements
uM = zeros(size(u));
uM(ranL) = u(1:L);
uM(ranLC) = u(L+1:end);
uM = reshape(uM, szu0);

%% Plot initial data and final image
figure(1);
subplot(2,2,1);
imshow(inputImg);
switch detectObject
    case {'tiger', 'ti'}
        rectangle('Position', [59, 30, 84-59, 52-30], 'EdgeColor', [1 0 0]); % face
        rectangle('Position', [66, 53, 84-66, 69-53], 'EdgeColor', [1 0 0]); % neck
        rectangle('Position', [66, 70, 146-66, 99-70], 'EdgeColor', [1 0 0]); % body
        rectangle('Position', [104, 7, 147-104, 41-7], 'EdgeColor', [0 1 1]); % background
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
        warning('whoops');
end
subplot(2,2,2);
imshow(uM, []);
colorbar;
subplot(2,2,3);
hist(uM(:), 20);

%% Threshold images
figure(2);
for j = 1:4
    subplot(2, 2, j);
    imagesc(uM > .4*j-1);
    colormap gray;
end

%%

figure;
imshow(imadd(double(uM <= 0.2), inputImg), [])
colormap gray;
