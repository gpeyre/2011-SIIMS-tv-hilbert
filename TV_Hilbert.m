% decomposition TV-Hilbert JF 
% masque_anneau enleve un anneau autour d'une frequence donnée (ttes
% orientations)

path(path, 'images/cartoons');
path(path, 'images/');
path(path, 'toolbox');

n=256;close all;


% pour energie: poids sur TV
lambda = .3;
q=16;
Text = rescale(load_image('fingerprint',n));
Text = Text - mean(Text(:));
Struct = rescale(load_image('bugsbunny_carreNB_ret',n));
M = 2*Struct+attenue_bord(Text);
% M = rescale(load_image('barb',n));
% M = image_test(n);

Text = image_texture(n,100/n,100);
Text = attenue_bord(Text,q);

Struct = load_image('bugsbunny_carreNB_ret',n);
Struct = rescale(squeeze(Struct(:,:,1)));

M = Struct+Text;

 operateur_psi = @psi_fourier
 operateur_psi_star = @psi_star_fourier

%Gabor
gamma = ifftshift(masque_anneau(n,18, 75));
%step size ( gradient descent)
mu = .9; 

%H^-1
%  x = -n/2+1:n/2;
%  [P,Q] = meshgrid(x,x);
%  gamma = 1./sqrt(2*(2-cos(2*pi*P/n)-cos(2*pi*Q/n)));
%  gamma(n/2,n/2)=0;
%  gamma = ifftshift(gamma);
%  mu = .0005;

% step size (Chambolle)
tau = 1/5.;

% pr affichage
affich = 5;

% init
y = gamma .* operateur_psi(M);
clear nrj;
nrj = zeros(1);
u=M;w=zeros(n,n,2);

evolution=1;i=1;
while(evolution~=0)
    % gradient descent step
    u2=u;
    u = u + real(mu * operateur_psi_star(gamma.*( y - gamma .* operateur_psi(u))));

    % TV denoising
    for k=1:10
        dw = grad( u/(lambda*mu) + div(w));
        w = w + tau * dw;
        d = repmat( sqrt(sum(w.^2,3)), [1 1 2] );
        w = w ./ max(d,ones(n,n,2));
    end
    u = u + lambda * mu * div(w);
    evolution = max(abs(u2(:)-u(:)));
    if (mod(i,affich)==0)
        imageplot({u M-u});
        pause(0.001);
        n_tmp = norm(gamma.*(operateur_psi(u))-y,'fro');
        nn = .5*n_tmp*n_tmp + lambda * TV(u);
        nrj(i/affich) = nn;
        if ((i/affich)~=1)
            diff = nrj(i/affich)-nrj(i/affich-1)
            if (diff>0)
                evolution=0;
            end
        end
    end
    i=i+1;
end