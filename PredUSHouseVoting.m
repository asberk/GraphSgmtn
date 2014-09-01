% % PredUSHouseVoting - Designed to test Bertozzi and Flenner's example, US House of Representatives Voting Records from 1984.
%     This script file loads a data set used by Bertozzi and Flenner to
%     verify that our implementation of the diffuse interface method for
%     arbitrary graphs yields the same results as that of Bertozzi and
%     Flenner, section 4.1. See the help file for PredImgSgmtn.m for more 
%     detail. 
%
% This script depends on the file(s):
% * gaussNorm.m
% * getGraphLaplacian.m
% * getImage.m
% * getWeights.m
% * imgnmz.m

%% clear the workspace
clear all; close all; clc;

%% Load and format dataset
fid = fopen('data/house-votes-84.data');
tmp = textscan(fid, '%s', 'Delimiter', ',');
tmp = reshape(tmp{1}, [17, 435]);
fclose(fid);
hv84 = strcmp('y', tmp) - strcmp('n', tmp);
hv84 = hv84(2:end, :);
DRlabels = strcmp('democrat', tmp(1,:)) - strcmp('republican',tmp(1,:));
clear tmp fid

%% Initialize parameters, variables
c = 1; 
dt = .1; 
epsilon = 2; 

M = 500; % Number of iterations for convergence
tau = 3; % use 3 for 1:5, 2 for 4:8; tau makes a difference! (probably in part due to numerical error)

%% Get weights
W = getWeights(hv84, @gaussNorm, 1, 1, tau);

%% Check if W is positive semi-definite
[R, p] = chol(W);

if p ~= 0 % then W is not positive semi definite:
    warning('W is not positive semi-definite.');
end

%% Compute graph Laplacian
[GL, Deg] = getGraphLaplacian(W); % use default settings to return symmetric graph Laplacian

%% Get EigenInformation of Graph Laplacian
[Phi,Lambda] = eig(GL);
%Phi = real(Phi);
Lambda = diag(Lambda);
%Lambda = real(Lambda)

%% Convex splitting for Graph Laplacian
% initialize (all are 1-by-L vectors)
u0 = zeros(size(DRlabels));
u0(1:5) = DRlabels(1:5);

display('Running convex splitting scheme for graph Laplacian...');
Parms.eta = abs(u0);
Parms.M = M;
Parms.L = 435;
Parms.c = c;
Parms.dt = dt;
Parms.epsilon = epsilon; 

u = pdeGraphSgmtn(u0, Phi, Lambda, Parms);

%% Display prediction accuracy
DRpred = sign(u);

predacc = sum(DRpred == DRlabels)/length(DRlabels);
fprintf('The algorithm determined party affiliation with %04.02f percent accuracy.\n', 100*predacc);

%% Compare results graphically
imagesc([DRlabels; DRpred]); colormap gray; colorbar; 
set(gcf, 'Color', [1 1 1]);
figure;
imagesc(abs(DRlabels - DRpred)<=0); colormap gray; colorbar;
set(gcf, 'Color', [1 1 1]);
