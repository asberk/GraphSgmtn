function Awide = widen(A, N, type)
% WIDEN(A,N,type) widens a matrix A by N pixels on each border by type TYPE, so that it can be convolved with padding. 
%   type : Dirichlet, Neumann (default)
szA = size(A);

Awide = zeros(szA+2*N);
Awide((1:szA(1))+N, (1:szA(2))+N) = A;
switch type
    case {'Dirichlet'}
        % done.
    otherwise % {'Neumann', 'Neu', 'N'}
        % % % edges
        % top edge
        Awide(1:N, (1:szA(2))+N) = bsxfun(@plus, Awide(1:N, (1:szA(2))+N), Awide(N+1, (1:szA(2))+N));
        % left edge
        Awide((1:szA(1))+N, 1:N) = bsxfun(@plus, Awide((1:szA(1))+N, 1:N), Awide((1:szA(1))+N, N+1));
        % right edge
        Awide((1:szA(1))+N, szA(2)+N+1:end) = bsxfun(@plus, Awide((1:szA(1))+N, szA(2)+N+1:end), Awide((1:szA(1))+N, szA(2)+N));
        % bottom edge
        Awide(szA(1)+N+1:end, (1:szA(2))+N) = bsxfun(@plus, Awide(szA(1)+N+1:end, (1:szA(2))+N), Awide(szA(1)+N, (1:szA(2))+N));
        
        % % % corners
        % topleft corner
        Awide(1:N+1, 1:N+1) = Awide(N+1,N+1);
        % topright corner
        Awide(1:N+1, end-N:end) = Awide(N+1, end-N);
        % bottom right corner
        Awide(end-N:end, end-N:end) = Awide(end-N, end-N);
        % bottom left corner
        Awide(end-N:end, 1:N+1) = Awide(end-N, N+1);
end
end