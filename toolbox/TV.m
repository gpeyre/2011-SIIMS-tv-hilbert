function tv = TV(M)
d = sqrt(sum(grad(M).^2,3));
tv = sum(d(:));