function WXY = getWeightsXY(X, Y, weightFunc, varargin)
% Computes co-weights for feature matrices X and Y
% % % Input
% X is (2N+1)^2-by-L; Y is (2N+1)^2-by-(n-L), 
% where N is the neighbourhood radius, n is the number of pixels in the
% original image, and L is the size of the random sample from the feature
% vectors/matrices.
% weightFunc: (2N+1)^2-by-(n-L) -> 1-by-(n-L) is a function (e.g., gaussNorm),
% which computes an adjacency weight for edges nu_x and mu_y represented by 
% the columns of X, Y respectively. 
% varargin{1} = displayProgres. Displays text progress bar if 1 (default 0)
% varargin{2} = displayElapsedTime. Displays elapsed time if 1 (default 0)
% % % Output
% WXY is L-by-(n-L) and represents the adjacency submatrix or the submatrix
% of the weight function that gives the co-weights between features from X
% and Y.

szX = size(X);
szY = size(Y);

WXY = zeros(szX(2), szY(2));

% housekeeping
nAuxParam = nargin(weightFunc) - 1; % number of auxillary parameters of weight function
opt = cell(1, nAuxParam);

if nargin > 5 + nAuxParam
    error('too many input arguments.');
end
switch nargin
    case 4 + nAuxParam
        displayProgress = varargin{1};
        for j = 1:nAuxParam
            opt{j} = varargin{j+1};
        end
    case 5 + nAuxParam
        displayProgress = varargin{1};
        displayElapsedTime = varargin{2};
        for j = 1:nAuxParam
            opt{j} = varargin{j+2};
        end        
    otherwise
        displayProgress = 0;
        displayElapsedTime = 0;
end

if displayProgress
    reverseStr = '';
end
if displayElapsedTime
    tic;
end

for j = 1:szX(2)
    % compute vector giving weight function values for feature j of X and
    % features of Y
    WXY(j,:) = weightFunc(bsxfun(@minus, X(:,j), Y), opt{:});
    
    % Display the progress
    if (displayProgress) && (mod(j,10^(round(log10(szX(2)))-2))==0)
        percentDone = 100 * j / szX(2);
        msg = sprintf('Percent done: %3.1f', percentDone);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
end

% more housekeeping
if displayProgress
    fprintf('\n');
end
if displayElapsedTime
    toc;
end
