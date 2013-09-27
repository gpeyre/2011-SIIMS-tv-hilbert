function Im_atten = attenue_bord(Im, marge)
n = size(Im,1);
Im_atten = Im;
if (nargin<2)
    marge = n/50;
end
for i=1:n
    for j=1:n
        if (i<=marge || i>(n-marge) || j<=marge || j>(n-marge))
            d = min([i-1 j-1 n-i n-j]);
            Im_atten(i,j) = sin(pi*d/(2*(marge-1)))^2 * Im(i,j);
        end
    end
end

