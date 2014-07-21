function arrExplore(W, j,k)
% ARREXPLORE Use arrow keys and left mouse button to explore 2D cross sections of a 4D array. Useful for exploring a weight matrix whose elements represent edge weights of an undirected graph, whose vertices are pixels in an image. 
% e.g., given W = getWeights(FtrArr, @gaussNorm, 1, 1, tau), 
% set W_arr = reshape(W, szu0, szu0) and execute arrExplore(W_arr, 1, 1),
% where szu0 = size(u0) and u0 is the original image (or a matrix having
% size equal to the original image). 

H = imshow(squeeze(W(j,k,:,:)));
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
drawnow;
title('(q)uit');
[x,y, btn] = ginput(1);

szW = size(W);

while btn ~= 113 % press 'q' to quit. 
    switch btn
        case 1 % left mouse button
            k = min(szW(2), round(x));
            j = min(szW(1), round(y));
            if k < 1
                k = 1;
            end
            if j < 1
                j = 1;
            end
        case 3 % right mouse button
            % so far, do nothing
        case 28 % left arrow key
            if k > 1
                k = k - 1;
            end
        case 29 % right arrow key
            if k < szW(2)
                k = k + 1;
            end
        case 30 % up arrow key
            if j > 1
                j = j - 1;
            end
        case 31 % down arrow key
            if j < szW(1)
                j = j + 1;
            end
    end
    set(H, 'CData', squeeze(W(j,k, :, :)));
    %set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    %drawnow;
    %title('(q)uit');
    [x,y, btn] = ginput(1);
end