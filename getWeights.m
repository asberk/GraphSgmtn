function W = getWeights(FtrArr, weightFunc, varargin)
% GETWEIGHTS computes the adjacency matrix correpsonding to a given feature array. 
%   FtrArr : a vector of feature vectors corresponding to each vertex in
%            the graph
%   weightFunc : the function used to compute the weight matrix; must
%                return a PSD matrix
%   varargin{1} : display progress? (default 0)
%   varargin{2} : display elapsed time? (default 0)
%   varargin{j:end} : auxiliary parameters to weightFunc

%% set-up
szFtr = size(FtrArr);
W = zeros(szFtr(2)); % nxn matrix where n = #(pixels of A)
% W(i,j) = (weight of edge/ adjacency strength) between pixel i and pixel j

%% Housekeeping
nAuxParam = nargin(weightFunc) - 1; % number of auxillary parameters of weight function
opt = cell(1, nAuxParam);
if nargin > 4 + nAuxParam
    error('too many input arguments.');
else
    displayProgress = 0;
    displayElapsedTime = 0;
end
switch nargin
    case 3 + nAuxParam
        displayProgress = varargin{1};
        for j = 1:nAuxParam
            opt{j} = varargin{j+1};
        end
    case 4 + nAuxParam
        displayProgress = varargin{1};
        displayElapsedTime = varargin{2};
        for j = 1:nAuxParam
            opt{j} = varargin{j+2};
        end
end
if displayProgress
    reverseStr = '';
end
if displayElapsedTime
    tic;
end

%% Main computation
for j = 1:szFtr(2)
    W(j,:) = weightFunc(bsxfun(@minus, FtrArr(:,j), FtrArr), opt{:});
%    wvec(j,:) = exp(-sum(bsxfun(@minus, ftrmat(:,j), ftrmat).^2));
    
    % Display the progress
    if (displayProgress) && (mod(j,10^(round(log10(szFtr(2)))-2))==0)
        percentDone = 100 * j / szFtr(2);
        msg = sprintf('Percent done: %3.1f', percentDone);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
end

%% More housekeeping
if displayProgress
    fprintf('\n');
end
if displayElapsedTime
    toc;
end

end