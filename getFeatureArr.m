function FtrArr = getFeatureArr(A, N, varargin)
% Feature Vectors/Matrices stored as matrix or high-dimensional array
% varargin{1} : aux param to widen: Neumann, Dirichlet, etc.
% varargin{2} : specify whether to return high-dim feature array ('highD'), or 2-dim
% feature array ('loD')

% ***
%     N.B. in loD, columns of FtrArr must be set in the order 
%     (k-1)*szA(1) + j and NOT in the order (j-1)*szA(2) +k
%     This is because Matlab column-wise vectorizes matrices --- i.e., 
%     Matlab turns [ 1 2 3 ; 4 5 6 ; 7 8 9 ] into [ 1 4 7 2 5 8 3 6 9 ].'
%     and the code needs to be consistent with that paradigm.
% ***

szA = size(A);
Awide = widen(A, N, 'Neumann'); % Widen A for convolution

if (nargin >= 4) && strcmp(varargin{2}, 'highD')
    FtrArr = zeros([szA, 2*N+1, 2*N+1]); % empty array to store result

    for j = (1:szA(1))+N % there's probably a loop-free way to do this...
        for k = (1:szA(2))+N
            FtrArr(j,k, :, :) = Awide((-N:N)+j, (-N:N)+k);
        end
    end
else % => loD
    FtrArr = zeros((2*N+1)^2, prod(szA));
    szftr = size(FtrArr);
    
    for k = 1:szA(2)
        for j = 1:szA(1)
            tmp = Awide((-N:N)+j+N, (-N:N)+k+N);
            FtrArr(:, (k-1)*szA(1) + j) = tmp(:);  % (j-1)*szA(2) + k
        end
    end
end



end