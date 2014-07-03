function FtrArr = getFeatureArr(A, N, varargin)
% Feature Vectors/Matrices stored as matrix or high-dimensional array
% varargin{1} : aux param to widen: Neumann, Dirichlet, etc.
% varargin{2} : specify whether to return high-dim feature array ('highD'), or 2-dim
% feature array ('loD')

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
    
    for j = 1:szA(1)
        for k = 1:szA(2)
            tmp = Awide((-N:N)+j+N, (-N:N)+k+N);
            FtrArr(:, (j-1)*szA(2) + k) = tmp(:);  % (j-1)*szA(2) + k
        end
    end
end



end