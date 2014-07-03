function W = getWeights(FtrArr, weightFunc, varargin)
% % % code for high-dim feature array
% wvec = zeros([szA, szA]);
% 
% for j1 = 1:szA(1)
%     display(j1);
%     for k1 = 1:szA(2)
%         display(k1);
%         for j2 = j1:szA(1) % suffices to compute for upper-triangular matrix.
%             if j2 == j1
%                 k2vec = k1:szA(2);
%             else
%                 k2vec = 1:szA(2);
%             end
%             for k2 = k2vec
%                 wvec(j1, k1, j2, k2) = exp(-norm(shiftdim(ftrvec(j1,k1,:,:) - ftrvec(j2,k2,:,:),2)));
%             end
%         end
%     end
% end

% % % lo-dim feature array
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

if 0 % doing this makes W non-positive definite!
W = W - diag(diag(W)); % make diagonal zero; assume no vertex has an edge to itself. 
end
%% More housekeeping
if displayProgress
    fprintf('\n');
end
if displayElapsedTime
    toc;
end

end