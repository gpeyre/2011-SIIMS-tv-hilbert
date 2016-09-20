% inpainting avec notre norme convexe

path(path, 'toolbox/');
path(path, 'images/');
path(path, 'images/cartoons');
n = 128;close all;
% pour energie: poids sur TV
lambda = .1;
% pour energie: poids sur norme de texture
mu = 1;
%sigma for gaussian
sigma = 1.4;

% size of the window
q = 16;
% spacing between window, redundancy is q/Dx
Dx = 4;


Text = image_texture(n,50/n,100);
%     Text(Text>0)=1;
%     Text(Text<0)=-1;

Struct = load_image('bugsbunny_carreNB',n);
%Struct = load_image('simple',n);
Struct = rescale(squeeze(Struct(:,:,1)));
% Text = rescale(load_image('fingerprint',n));
% Text = Text-mean(Text(:));

Text = attenue_bord(Text,q/2);
M0 = 2*Struct+Text;
M = M0;


% amount of removed pixels and size of holes
rho = .3; siz = 15;
% random mask, mask==1 for removed pixels
mask = masque(n,rho,siz);

% remove pixels
y = (1-mask) .* M;
% boundary handling
options.bound = 'sym';
% energy conservation
options.normalization = 'tightframe';
% windowing function
options.window_type = 'sin';

% regularization parameter for the TV norm
epsilon = 1e-2;
% gradient descent step size
tau = .01;
tolinit=25;

clear nrj;
nrj = zeros(1);
nrj_glob = zeros(1);
MF = perform_windowed_fourier_transform(M,q,Dx,n, options);
Orientationsss = zeros(1,size(MF,3),size(MF,4),2);
Imgs_u = zeros(1,n,n);
Imgs_v = zeros(1,n,n);
u=zeros(n,n);v=y;
w=zeros(n,n,2);
evolution=1;i=0;
evolution_glob=1;i_glob=1;

while(evolution_glob>-0.001)
    u2_glob=u;
    v2_glob=v;
    FM = perform_windowed_fourier_transform(v,q,Dx,n, options);

    Orientations = estimate_orientations(FM,mask,q,Dx,max(10,tolinit-i_glob));
    Gamma = gabor_weights(Orientations,q,sigma);

    Orientationsss(i_glob,:,:,:) = Orientations;
    Imgs_u(i_glob,:,:) = u;
    Imgs_v(i_glob,:,:) = v;
    evolution=1;
     stop = (0.5/(5*i_glob));
    while(evolution>stop)
%         mu =mu*.9;
%         lambda =lambda*.9;
        i = i + 1;
        u2=u;v2=v;
        %energy
%         normT = perform_windowed_fourier_transform(v,q,Dx,n, options);
%         normT = Gamma .* normT;
%         normT = norm(normT(:),2);
%         normT = normT*normT;
%         normTV = TV(u);
%         normL2 = norm(y(:)-(1-mask(:)).*(u(:)+v(:)));
%         normL2 = normL2*normL2;
%         nn = mu*normT + lambda*normTV + .5 * normL2;
%         nrj(i) = nn;

        % v is fixed: TV denoising
        y_bar = y - (1-mask).*v;
        for k=1:min(400,50+5*i_glob)
            Gr = grad(u);
            d = sqrt( epsilon^2 + sum(Gr.^2,3) );
            G = -div( Gr./repmat(d, [1 1 2]));
            G = lambda * G + (1-mask).*u-y_bar;
            u = u - tau*G;
        end

        % u is fixed: conjugate gradient
        y_bar = y - (1-mask).*u;
        v = gradconj(y_bar,Gamma,mu,q,Dx,n,options,stop/(50+i*5),mask);

        evolution = max(max(abs(u2(:)-u(:))),max(abs(v2(:)-v(:))));
        disp(['i_glob: ' num2str(i_glob) ' - evol: ' num2str(evolution) '>stop: ' num2str(stop)]);
        imageplot({u v u+v},{'u' 'v' 'u+v'});
        pause(0.1);
        if (i>(3+i_glob))
            evolution=0;
        end
    end
    evolution_glob = max(max(abs(u2_glob(:)-u(:))),max(abs(v2_glob(:)-v(:))));
    normT = perform_windowed_fourier_transform(v,q,Dx,n, options);
    normT = Gamma .* normT;
    normT = norm(normT(:),2);
    normT = normT*normT;
    normTV = TV(u);
    normL2 = norm(y(:)-(1-mask(:)).*(u(:)+v(:)));
    normL2 = normL2*normL2;
    nn = mu*normT + lambda*normTV + .5 * normL2;
    nrj_glob(i_glob)=nn;
     if (i_glob~=1)
        diff = nrj_glob(i_glob)-nrj_glob(i_glob-1);
        disp(['i_glob=' num2str(i_glob) ' diff_nrj=' num2str(diff) ' evolution_glob=' num2str(evolution_glob)]);
        if (diff>0)
            evolution_glob=0;
        end
    end
    i_glob = i_glob+1;
end
