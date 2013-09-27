% Decomposition Struct/Text avec notre norme convexe avec bruit

path(path, 'toolbox/');
path(path, 'images/');
path(path, 'images/cartoons/');
path(path, 'images/textures/');
close all;
n = 512;
% pour energie: poids sur TV
lambda = .15;
% pour energie: poids sur norme de texture
mu = 5;
alpha = 1;

lambda = alpha*lambda;
mu = alpha*mu;

% size of the window
q = 32;
% spacing between window, redundancy is q/Dx
Dx = 16;

%sigma for gaussian
sigma = 2.2;
% tolerance pour orientations
tol = 6;


% Text = image_texture(n,100/n,100);
% Text = attenue_bord(Text,q);
% Struct = load_image('bugsbunny_carreNB',n);
% Struct = rescale(squeeze(Struct(:,:,1)));
%M0 = 2*Struct+Text;

M0 = load_image('desert_sable4_carre',n); M0 =  rescale(squeeze(M0(:,:,1)));
M0 = rescale(load_image('barb',n));
%R0 = .05*randn(n,n);
M = M0 + R0;


% boundary handling
options.bound = 'sym';
% energy conservation
options.normalization = 'tightframe';
% windowing function
options.window_type = 'sin';

% step size (Chambolle)
tau = 1/5.;

clear nrj;
nrj = zeros(1);
nrj_glob = zeros(1);
MF = perform_windowed_fourier_transform(M,q,Dx,n, options);
Orientationsss = zeros(1,size(MF,3),size(MF,4),2);
Imgs_u = zeros(1,n,n);dif=0;
Imgs_v = zeros(1,n,n);
u=zeros(n,n);v=M;
w=zeros(n,n,2);
evolution=1;i=0;
evolution_glob=1;i_glob=1;

temps = zeros(1);
temps(1)=0;
tic;
Zero = zeros(n,n);
while(evolution_glob>0.001)
    u2_glob=u;
    v2_glob=v;
    Estim = v;
    F_estim = perform_windowed_fourier_transform(Estim,q,Dx,n, options);
    Orientations = estimate_orientations(F_estim,Zero,q,Dx,tol);
    Gamma = gabor_weights(Orientations,q,sigma);
    Orientationsss(i_glob,:,:,:) = Orientations;
    Imgs_u(i_glob,:,:) = u;
    Imgs_v(i_glob,:,:) = v;
    evolution=1;
    stop = 0.1/i_glob;
    while(evolution>stop)
        i = i + 1;
        u2=u;v2=v;
%         %energy
%         normT = perform_windowed_fourier_transform(v,q,Dx,n, options);
%         normT = Gamma .* normT;
%         normT = norm(normT(:),2);
%         normT = normT*normT;
%         normTV = TV(u);
%         normL2 = norm(M(:)-u(:)-v(:));
%         normL2 = normL2*normL2;
%         nn = mu*normT + lambda*normTV + .5 * normL2;
%         nrj(i) = nn;

        % v is fixed: TV denoising
        for k=1:min(50,(10*i_glob))
            dw = grad( (M-v)/lambda + div(w));
            w = w + tau * dw;
            d = repmat( sqrt(sum(w.^2,3)), [1 1 2] );
            w = w ./ max(d,ones(n,n,2));
        end
        u = M-v + lambda * div(w);

        % u is fixed: conjugate gradient
%        v = gradconj(M-u,Gamma,mu,q,Dx,n,options,50*stop);
        %synthese
         Mx = perform_windowed_fourier_transform(M-u,q,Dx,n, options);
         Tmp = (1./(2*mu*Gamma.*Gamma+1)).*Mx;
         v = perform_windowed_fourier_transform(Tmp,q,Dx,n, options);
        
        evolution = max(max(abs(u2(:)-u(:))),max(abs(v2(:)-v(:))));

        imageplot({u v M-u-v});
        pause(0.001);
    end
    evolution_glob = max(max(abs(u2_glob(:)-u(:))),max(abs(v2_glob(:)-v(:))));
    %energy
    normT = perform_windowed_fourier_transform(v,q,Dx,n, options);
    normT = Gamma .* normT;
    normT = norm(normT(:),2);
    normT = normT*normT;
    normTV = TV(u);
    normL2 = norm(M(:)-u(:)-v(:));
    normL2 = normL2*normL2;
    nn = mu*normT + lambda*normTV + .5 * normL2;
    nrj_glob(i_glob)=nn;
    if (i_glob~=1)
        dif = nrj_glob(i_glob)-nrj_glob(i_glob-1);
        if (dif>0)
            evolution_glob=0;
        end
    end
   i_glob = i_glob+1;
   disp(['i_glob=' num2str(i_glob) ' nrj_glob=' num2str(nrj_glob(i_glob-1)) ' diff=' num2str(dif) ' evolution_glob=' num2str(evolution_glob)]);
   temps(i_glob) = toc;
end

i=i_glob-1;
close all;
u = squeeze(Imgs_u(i,:,:));
v = squeeze(Imgs_v(i,:,:));
figure;imageplot({M u+v},{'M' 'u+v'})
figure;imageplot({u v M-u-v},{'u' 'v' 'w'})
snr(M0,u+v)