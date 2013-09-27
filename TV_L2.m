% decomposition TV-L2

path(path, 'images/cartoons');
path(path, 'images/');
path(path, 'toolbox');
n = 512;

 M0 = rescale(load_image('barb',n));
 R = .05 * randn(n,n);
 M = M0 + R;


mu=0.1;


tau = 1/4.;
nrj=zeros(1);

v=zeros(n,n);

     p = zeros(n,n,2);

     evolution=1;k=1;
     while(evolution>0.0001)
         p2 = p;
         p = p + tau*grad( div(p) - M/mu );
         d = repmat( sqrt(sum(p.^2,3)), [1 1 2] );
         p = p ./ max(d,ones(n,n,2));

        
        v2 = v;
        v = mu*div(p);
        if (mod(k,10)==0)
          imageplot({M-v v})
          pause(0.001);
        end

        n_tmp = norm(v,'fro');
        nn = TV(M-v) + n_tmp*n_tmp/(2*mu);
        nrj(k) = nn; 
        evolution = max(abs(v2(:)-v(:)))
        k=k+1;
     end
    imageplot(M-v)
%end