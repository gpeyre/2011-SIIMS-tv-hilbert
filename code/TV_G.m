% decomposition TV-G Meyer

path(path, 'images/cartoons');
path(path, 'images/');
path(path, 'toolbox');

n=256;close all;

 M0 = rescale(load_image('barb',n));
 R = .05 * randn(n,n);
 M = M0 + R;

% step size
tau = 1/4.;

nrj = zeros(1);
u=zeros(n,n);v=zeros(n,n);
p = zeros(n,n,2);q = zeros(n,n,2);

alpha = 1/300;
mu = 0.05;

epsilon = 0.00001;

evolution = 1.;

i=1;
nn = (TV(u) + norm(M-u-v,'fro')/(2*alpha))
nrj(i) = nn;
while(evolution>epsilon)
    u2 = u;
    v2 = v;
	for k=1:50
        p = p + tau*grad( div(p) - (M-u)/mu );
        d = repmat( sqrt(sum(p.^2,3)), [1 1 2] );
        p = p ./ max(d,ones(n,n,2));
    end
    v = mu*div(p);
    for k=1:50
         q = q + tau*grad( div(q) - (M-v)/alpha );
         d = repmat( sqrt(sum(q.^2,3)), [1 1 2] );
         q = q ./ max(d,ones(n,n,2));
    end
    u = M - v - alpha*div(q);
    
    imageplot({u v});
    pause(0.01);

    nn = (TV(u) + norm(M-u-v,'fro')/(2*alpha));
    nrj(i) = nn;
    if (i~=1)
        diff=nrj(i)-nrj(i-1)
    end
    evolution = max(abs(u2(:)-u(:)))
    i = i + 1;
end

imageplot(u)