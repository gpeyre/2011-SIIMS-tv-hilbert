function M = masque(n,rho,siz)
% for inpainting
%siz: siz of holes
%rho amount of remove pixels
%n size of the image

M = zeros(n,n);
tot = 0;
for i=0:n/siz-1
    for j=0:n/siz-1
        if (rand < rho)
            for x=0:siz-1
                for y=0:siz-1
                    M(1+siz*i+x,1+siz*j+y)=1;
                end
            end
        end
    end
end
        

