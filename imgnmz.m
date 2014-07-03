function B = imgnmz(A)
A = A - min(A(:));
B = A./max(A(:));