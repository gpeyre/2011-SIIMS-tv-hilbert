function G=masque_gabor(n,sigma,freq1,freq2)

G = zeros(n,n);
for f1=0:n-1
    for f2=0:n-1
        G(f1+1,f2+1) =  (1-exp(-((f1-n/2 - freq1)^2  +  (f2-n/2 - freq2)^2) / (2*sigma*sigma)));
        G(f1+1,f2+1) = G(f1+1,f2+1)* (1-exp(-((f1-n/2 + freq1)^2  +  (f2-n/2 + freq2)^2) / (2*sigma*sigma)));
    end
end 
G(G<0.4)=0;