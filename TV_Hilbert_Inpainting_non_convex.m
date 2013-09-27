% inpainting avec notre norme non convexe

path(path, 'toolbox/');
path(path, 'images/');
path(path, 'images/cartoons');
path(path, 'images/textures');
close all;
n = 256;
% pour energie: poids sur TV
lambda0 = .06;
% pour energie: poids sur norme de texture
mu0 = .1;
% attache aux données (1/alpha)
alpha = 1;

lambda = alpha*lambda0;
mu = alpha*mu0;


%ptr_h = @test_h;
%ptr_diff_h = @test_diff_h;
%ptr_h = @h_identity;
%ptr_diff_h = @diff_h_identity;
ptr_h = @test_h2;
ptr_diff_h = @test_diff_h2;
a=4;b=0;

% size of the window
q = 32;
% spacing between window, redundancy is q/Dx
Dx = 8;


% gradient descent step size
mu_tv = 1/2;
tau_tv = 1/4;
tau = 0.01;
tolerance = 5;


% boundary handling
options.bound = 'sym';
% energy conservation
options.normalization = 'tightframe';
% windowing function
options.window_type = 'sin';

% mask = zeros(n,n);
% taille0 = 25;
% taille=taille0;
% mask(n/2-taille:n/2+taille,n/2-taille:n/2+taille) = 1;
% mask = load_image('mask_desert_sable4',n);mask = squeeze(mask(:,:,1));mask(mask<128)=0;mask(mask~=0)=1;mask = 1-mask;
mask = load_image('mask2_bis',n);mask = squeeze(mask(:,:,1));mask(mask<128)=0;mask(mask~=0)=1;mask = 1-mask;
%mask=zeros(n,n);

% mask = zeros(n,n);
% taille = 2*floor(q/4);
% xc = 3*n/4;
% yc = 3*n/4;
% for i=1:n
%     for j=1:n
%         if (sqrt((i-xc)*(i-xc)+(j-yc)*(j-yc))<taille)
%             mask(i,j)=1;
%         end
%     end
% end
% xc = n/4;
% yc = n/4; 
% mask(xc-taille:xc+taille,yc-taille:yc+taille) = 1;

P = compute_patch(mask, q, Dx, options);
mask_amplitudes = squeeze(sum(sum(P))>.2*q*q);

%
% Text = image_texture(n,50/n,100);
%      Text(Text>0)=1;
%      Text(Text<0)=-1
%Text = test_h(Text,a,b);

Struct = load_image('bugsbunny_carreNB',n);
%Struct = load_image('simple',n);
Struct = rescale(squeeze(Struct(:,:,1)));
Text = rescale(load_image('fingerprint',n),.1,.9);
%Text = Text-mean(Text(:));
% load('texture_sinusoide')
% Sinusoid = v(1:n,1:n);
 Text = image_texture(n,70/n,100); 
 Text = feval(ptr_h,Text,a,b);

Text = attenue_bord(Text,q/2);
M0 = 2*Struct+Text;
%M0 = load_image('desert_sable4_carre',n); M0 =  rescale(squeeze(M0(:,:,1)));
R = randn(n,n);
M = M0+.05*R;

y = (1-mask) .* M;


Imgs_u = zeros(1,n,n);
Imgs_v = zeros(1,n,n);


%%% INIT
%v=(1-mask).*Text;u=(1-mask).*Struct;
v = y/2; u = y/2;

load('images\resultats\BB_inpainting\gros-trou-non-cvxe\bb_inpainting_mask2bis_noncvxe_h_id.mat');

v = y/2; u = y/2;
lambda = .3;

w=zeros(n,n,2);
evolution=1;i=0;
evolution_glob=1;i_glob=1;

while(evolution_glob>-0.001)
    u2_glob=u;
    v2_glob=v;
%              alpha = alpha*.95;
%              lambda = alpha*lambda0;
%               mu = alpha*mu0;

    Imgs_u(i_glob,:,:) = u;
    Imgs_v(i_glob,:,:) = v;

    i = i + 1;
    u2=u;v2=v;

    % v is fixed: TV denoising
    y_bar = y - .15*(1-mask).*feval(ptr_h,v/.15,a,b);

    for k1=1:20
        u2=u;
        u = u + mu_tv*(1-mask).*(y_bar-u);
        % TV denoising
        for k2=1:10
            dw = grad( u/(lambda*mu_tv) + div(w));
            w = w + tau_tv * dw;
            d = repmat( sqrt(sum(w.^2,3)), [1 1 2] );
            w = w ./ max(d,ones(n,n,2));
        end
        u = u + lambda * mu_tv * div(w);
    end


    % u is fixed: conjugate gradient
    y_bar = y - (1-mask).*u;

    mask_erod = erosion(mask,floor(((i_glob-10)/2)));

    Estim = (1-mask_erod).*v;
    [F_estim, Weight] = perform_windowed_fourier_transform(Estim,q,Dx,n, options);
    Orientations = estimate_orientations(F_estim,mask_erod,q,Dx,tolerance);
    %Orientations = estimate_orientations_zeropadding(Estim,mask_erod,q,Dx,5);
    Amplitudes = ones(size(Orientations,1),size(Orientations,2));
    
        Chi = non_convex_weight(Orientations, Amplitudes, Weight, q);
        Amplitudes = estimate_amplitudes(F_estim, Chi);%, mask_erod, q, Dx, options);
        diff=1;
        while(diff>1e-4)
            Amplitudes2= perform_blurring(Amplitudes,3);
            Amplitudes2(mask_amplitudes==0) = Amplitudes(mask_amplitudes==0);
            diff = norm(Amplitudes2-Amplitudes,'fro');
            Amplitudes = Amplitudes2;
        end
    
     amplit_moy = mean(Amplitudes(:));
    Chi = non_convex_weight(Orientations, Amplitudes, Weight, q);

    %    v = gradconj(y_bar,Gamma,mu,q,Dx,n,options,stop/(50+i*5),mask);
    [v,niter] = gradient_descent_armijo(v,mu,alpha,1/(2*i_glob),Chi,y_bar,mask,q,Dx,n,options,ptr_h,ptr_diff_h,a,b,10+max(0,floor(i_glob/5)),amplit_moy);
%    [v,niter] = gradient_descent_armijo((1-mask).*attenue_bord(Sinusoid,q/2)+0.0*randn(n,n),mu,alpha,-1/(2*i_glob),Chi,y_bar,mask,q,Dx,n,options,ptr_h,ptr_diff_h,a,b,30+max(0,floor(i_glob/5)-25));
%   [v,niter] = gradient_descent_armijo((1-mask).*Text,mu,alpha,1/(2*i_glob),Chi,(1-mask).*Text,mask,q,Dx,n,options,ptr_h,ptr_diff_h,a,b,30+max(0,floor(i_glob/5)-25));
    evolution = max(max(abs(u2(:)-u(:))),max(abs(v2(:)-v(:))));
    disp(['(i_glob: ' num2str(i_glob) ' - evol: ' num2str(evolution)]);
    imageplot(u,'u',2,3,1);    imageplot(v,'v',2,3,2);    imageplot(amplit_moy*feval(ptr_h,v/amplit_moy,a,b),'h(v)',2,3,3);
    imageplot(u+amplit_moy*feval(ptr_h,v/amplit_moy,a,b),'u+h(v)',2,3,4);  imageplot(M,'M',2,3,5);  imageplot(y-(1-mask).*(u+amplit_moy*feval(ptr_h,v/amplit_moy,a,b)),'w',2,3,6); 

    pause(0.001);
    if (i>(3+i_glob))
        evolution=0;
    end
    evolution_glob = max(max(abs(u2_glob(:)-u(:))),max(abs(v2_glob(:)-v(:))));
    i_glob = i_glob+1;
end

i=i_glob-1;
close all;
u = squeeze(Imgs_u(i,:,:));
v = squeeze(Imgs_v(i,:,:));
w = (1-mask).*(y-u-amplit_moy*feval(ptr_h,v/amplit_moy,a,b));
figure;imageplot({y u+amplit_moy*feval(ptr_h,v/amplit_moy,a,b) amplit_moy*feval(ptr_h,v/amplit_moy,a,b)} ,{'y' 'u+h(v)' 'h(v)'})
figure;imageplot({u v w},{'u' 'v' 'w'})
snr(M0,u+amplit_moy*feval(ptr_h,v/amplit_moy,a,b))

figure;
imageplot(u+amplit_moy*feval(ptr_h,v/amplit_moy,a,b),'u+h(v)',2,3,2)
imageplot(y,'y',2,3,1)
imageplot(M,'M',2,3,3)
imageplot(v,'v',2,3,4)
imageplot(amplit_moy*feval(ptr_h,v/amplit_moy,a,b),'h(v)',2,3,5)
imageplot(Text,'Text',2,3,6)