function uNew = reOrder(u, szU, ranL)
%REORDER reformats the permuted vector u as the original image/matrix U.
%   reOrder(u, szU, ranL) takes a permutation vector ranL, which 
%   corresponds to an image U. It then unpermutes the elemtns of the vector
%   u and returns u as an image of size szU. Where u is a permutation of 
%   U(:).', this process returns the original image U.
%   
%   This function is especially useful for debugging the loop sequence of
%   the functional minimization in Pred*.m

lenRanL = length(ranL);
tmp = zeros(size(u));
tmp(ranL) = u(1:lenRanL);
tmp(setdiff(1:length(u), ranL)) = u(lenRanL+1:end);
uNew = imgnmz(reshape(tmp, szU));