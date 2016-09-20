% Inpainting TV

path(path, 'toolbox/');
path(path, 'images/');
path(path, 'images/cartoons');
n = 256;
mu = 1;
lambda = 0;

M = rescale(load_image('bugsbunny_carreNB',n));
mask = load_image('mask2',n);mask=squeeze(mask(:,:,1));mask(mask<128)=0;mask(mask~=0)=1;mask = 1-mask;

% regularization parameter for the TV norm
epsilon = 1e-2;
% gradient descent step size
tau = .005;
% initialization


% remove pixels
y = (1-mask).*M;
% R = rand(n,n);
% y(mask==1) = R(mask==1);
% display


Mtv = y;

nrj = zeros(1);
evolution = 1;
i=1;k=1;
while(evolution>0.000001)
    Mtv2 = Mtv;
    Gr = grad(Mtv);
    d = sqrt( epsilon^2 + sum(Gr.^2,3) );
    G = -div( Gr./repmat(d, [1 1 2]));
    %G = mu*G + 2*lambda*(mask==0) .* (Mtv-y);
    Mtv = Mtv - tau*G;
    Mtv(mask==0) = y(mask==0);
    k=k+1;
    if (mod(k,20)==0)
        n_tmp = norm(Mtv(mask==0)-y(mask==0),'fro');
        nn = mu*TV(Mtv) + lambda*n_tmp*n_tmp;
        nrj(i) = nn;
        evolution = max(abs(Mtv(:)-Mtv2(:)))
        if (i>1)
            diff = nrj(i-1)-nrj(i)
%             if (diff<0)
%                 evolution=0;
%             end
        end
        i=i+1;
        imageplot({y Mtv});
        pause(0.001);
        disp([num2str(i)]);
    end
end

clf;

