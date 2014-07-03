function [Wxx, Wxy] = getNormalizedSampleWeights(X, Y, weightFunc, varargin)
% Returns normalized sample weights Wxx, Wxy as given by (3.9) of Bertozzi and Flenner.

%% Compute unnormalized matrices
display('computing X-weights...');
W_XX = getWeights(X, weightFunc, 1, 1, varargin);
display('computing co-weights...');
W_XY = getWeightsXY(X, Y, weightFunc, 1, 1, varargin);

%% Normalize weight matrices
d_X = sum(W_XX, 2) + sum(W_XY, 2);
d_Y = sum(W_XY.', 2) + sum(W_XY.'/W_XX*W_XY, 2);

if any(d_Y<0)
    warning([num2str(sum(d_Y(:)<0)), ' value(s) < 0; use their absolute value(s) instead?']);
    absVal = input('answer [y]es (default) or [n]o >> ', 's');
    if sum(strcmp(absVal, {'yes', 'y', ''})) || isempty(absVal) 
        d_Y = abs(d_Y);
    else
        warning('using d_Y as is.');
    end
end

s_X = sqrt(d_X);
s_Y = sqrt(d_Y);

% s_X = sqrt(sum(W_XX, 2) + sum(W_XY, 2));
% s_Y = sqrt(sum(W_XY.', 2) + sum(W_XY.'/W_XX*W_XY, 2));

Wxx = W_XX./(s_X*s_X.'); % s_X*s_X.' is an L-by-L outer product matrix
Wxy = W_XY./(s_X*s_Y.'); % s_X*s_Y.' is an L-by-(n-L) outer product matrix