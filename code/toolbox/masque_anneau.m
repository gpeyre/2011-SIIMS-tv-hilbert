function G=masque_anneau(n,sigma,freq)

G = zeros(n,n);
for f1=0:n-1
    for f2=0:n-1
        G(f1+1,f2+1) =            exp(-(sqrt((f1-n/2)^2+(f2-n/2)^2)-freq)^2/(2*sigma*sigma));
    end
end 
G = abs(1-G);
%G = 1 - G/(2*sigma*sqrt(2*pi));