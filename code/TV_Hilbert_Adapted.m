% Decomposition Struct/Text avec notre norme convexe sans bruit

path(path, 'toolbox/');
path(path, 'images/');
path(path, 'images/cartoons');
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parametres:
%taille image
n = 256;
% pour energie: poids sur TV 
lambda = .3;
% size of the window
q = 24;
% spacing between window, redundancy is q/Dx
Dx = 6;
%sigma for gaussian
sigma = 1.4;
%pour T&J_fingerprint:
%n = 512;lambda = .6;q = 48;Dx = 16;sigma = 1.2;
%pour Bugsbun_fingerprint:
n = 512;lambda = .3;q = 48;Dx = 16;sigma = 2.2;
%n = 512;lambda = .35;q = 36;Dx = 12;sigma = 2.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Text = image_texture(n,100/n,100);
 Text(Text>.2)=1;
 Text(Text<-.2)=-1;
Text = attenue_bord(Text,q);

Text = rescale(load_image('fingerprint',n));
Text = Text-mean(Text(:));
Text = attenue_bord(Text,q/2);

Struct = load_image('bugsbunny_carreNB',n);
Struct = rescale(squeeze(Struct(:,:,1)));

M = 2*Struct+Text;
%M = image_test(n);

%P = compute_patch(M, q, Dx, options);

% ici notation du papiwer donc differente du code de
% perform_windowed_transform (q est la taille de la fenetre)

% boundary handling
options.bound = 'sym';
% energy conservation
options.normalization = 'tightframe';
% windowing function
options.window_type = 'sin';
options.window_type = 'constant';

% step size (Chambolle)
tau = 1/4.;
% step size ( gradient descent)
mu = 1;


nrj = zeros(1);
nrj_global = zeros(1);
[MF,Weight] = perform_windowed_fourier_transform(M,q,Dx,n, options);
Orientationsss = zeros(1,size(MF,3),size(MF,4),2);
Imgs = zeros(1,n,n);
u=zeros(n,n);w=zeros(n,n,2);
evolution=1;i=1;
evolution_glob=1;i_glob=1;

while(evolution_glob>0.001)
    % init
    u2_glob=u;
    MF = perform_windowed_fourier_transform(M,q,Dx,n, options);
    MTest = perform_windowed_fourier_transform(M-u,q,Dx,n, options);
    
    Orientations = estimate_orientations(MTest,zeros(n,n),q,Dx);
    Gamma = gabor_weights(Orientations,q,sigma);

%    imageplot({M M});
%    display_orientations(Orientations,Dx,n,1);
    Orientationsss(i_glob,:,:,:) = Orientations;
    Imgs(i_glob,:,:) = u;
    y = Gamma .* MF;
    evolution=1;
    i0=i;
%    while(evolution>(0.05/i_glob))
    while((i-i0)<(i_glob+1))
        %energy
        u2=u;
         Mu = perform_windowed_fourier_transform(u,q,Dx,n, options);
         Tmp = Gamma.*Mu-y;
%         n_tmp = norm(Tmp(:),2);
%         nn = .5*n_tmp*n_tmp + lambda * TV(u);
%         nrj(i) = nn;
        i = i + 1;
        % gradient descent step
        MTmp = -Gamma.*Tmp;
        Tmp = perform_windowed_fourier_transform(MTmp,q,Dx,n, options);
        u = u + real(mu * Tmp);

        % TV denoising
        for k=1:(10+i_glob)
            dw = grad( u/(lambda*mu) + div(w));
            w = w + tau * dw;
            d = repmat( sqrt(sum(w.^2,3)), [1 1 2] );
            w = w ./ max(d,ones(n,n,2));
        end
        u = u + lambda * mu * div(w);
        evolution = max(abs(u2(:)-u(:)))
        if (mod(i,1)==0)
            imageplot({u M-u});
            pause(0.001);
        end
    end

    evolution_glob = max(abs(u2_glob(:)-u(:)))
    Mu = perform_windowed_fourier_transform(u,q,Dx,n, options);
    Tmp = Gamma.*Mu-y;
    n_tmp = norm(Tmp(:),2);
    nn = .5*n_tmp*n_tmp + lambda * TV(u);
    nrj_global(i_glob) = nn;
    i_glob = i_glob+1
%    evolution_glob=0.0000001;
end

imageplot({u M-u});