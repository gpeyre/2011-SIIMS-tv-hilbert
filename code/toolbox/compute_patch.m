function P = compute_patch(M, w, q, options)

options.null = 0;
bound = getoptions(options, 'bound', 'sym');
n = size(M,1);

% perform sampling
t = 1:q:n+1;
[Y,X] = meshgrid(t,t);
p = size(X,1);

if mod(w,2)==1
% w = ceil((w-1)/2)*2+1;
    w1 = (w-1)/2;
    t = -w1:w1;
else
    t = -w/2:w/2-1;
end
[dY,dX] = meshgrid(t,t);

X = reshape(X,[1 1 p p]);
Y = reshape(Y,[1 1 p p]);
X = repmat( X, [w w 1 1] );
Y = repmat( Y, [w w 1 1] );
dX = repmat( dX, [1 1 p p] );
dY = repmat( dY, [1 1 p p] );

X1 = X+dX;
Y1 = Y+dY;

switch lower(bound)
    case 'sym'
        X1(X1<1) = 1-X1(X1<1);
        X1(X1>n) = 2*n+1-X1(X1>n);
        Y1(Y1<1) = 1-Y1(Y1<1);
        Y1(Y1>n) = 2*n+1-Y1(Y1>n);
    case 'per'
        X1 = mod(X1-1,n)+1;
        Y1 = mod(Y1-1,n)+1;
end


I = X1 + (Y1-1)*n;

P = M(I);