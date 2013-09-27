function G = gabor_weights(Ori,q,sigma)

if (~exist('sigma','var'))
    sigma=1;
end
G = ones(q,q,size(Ori,1),size(Ori,2));

for i=1:size(Ori,1)
    for j=1:size(Ori,1)
        x = Ori(i,j,1);
        y = Ori(i,j,2);
        if (x~=0 || y ~=0)
            G(:,:,i,j) = masque_gabor(q,sigma,x,-y);
        end
    end
end
%G(G<.4)=0;