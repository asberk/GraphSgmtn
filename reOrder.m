function uNew = reOrder(u, szU, ranL)
%% reOrder reformats the vector as the original U.
%    Given a permutation vector, corresponding to an image, U, this 
%    function unpermutes the elements of the vector u and returns u as an 
%    image of size szU. 
%    
%    This function is especially useful for debugging the loop sequence of
%    the functional minimization in Pred*.m

lenRanL = length(ranL);
tmp = zeros(size(u));
tmp(ranL) = u(1:lenRanL);
tmp(setdiff(1:length(u), ranL)) = u(lenRanL+1:end);
uNew = imgnmz(reshape(tmp, szU([2 1]))).';