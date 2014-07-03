function [A, sizeA] = getImage(str)
%% getImage retrieves the image filename.ext from pathname
%  where str = [pathname filename.ext]

B = imread(str);

if min(B(:)) < 0 
    error('cannot have image values less than zero');
end

if isinteger(B) 
    if length(size(B)) == 3 % colour image
        B = rgb2gray(B);
    end
    if max(B(:)) > 1
        % probably rgb image
        B = double(B)./255;
    % else
    %   binary image; don't touch
    end        
else % B is a floating point array
    % don't know how to deal with colour images
    if max(B(:)) > 1 % not binary image
        warning('are you sure you want a floating image taking values greater than 1?');
        if isfinite(max(B(:)))
            B = B./max(B);
        else
            error('image cannot have infinite values');
        end
    % else
    %   already in the format we want
    end
end

A = B;
sizeA = size(A);

