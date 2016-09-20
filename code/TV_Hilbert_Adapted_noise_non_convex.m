% Decomposition Struct/Text avec notre norme non convexe

path(path, 'toolbox/');
path(path, 'images/');
path(path, 'images/cartoons/');
close all;
n = 256;
% pour energie: poids sur TV
lambda0 = .2;
% pour energie: poids sur norme de texture
mu0 = .2;
 alpha = 1;
 
 lambda = alpha*lambda0;
 mu = alpha*mu0;

% size of the window
q = 32;
% spacing between window, redundancy is q/Dx
Dx = 4;

Text = image_texture(n,50/n,100);
% Text(Text>0)=1;
% Text(Text<0)=-1;

ptr_h = @test_h;
ptr_diff_h = @test_diff_h;
a=10;b=0;
Text = test_h(Text,a,b);

Struct = load_image('bugsbunny_carreNB',n);
%Struct = load_image('simple',n);
Struct = rescale(squeeze(Struct(:,:,1)));
Text = rescale(load_image('fingerprint',n));
Text = Text-mean(Text(:));

Text = attenue_bord(Text,q/2);
M0 = 2*Struct+Text;
%R0 = .2*randn(n,n);
%M = M0 + R0;
M = M0;


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
[MF,weight] = perform_windowed_fourier_transform(M,q,Dx,n, options);
Orientationsss = zeros(1,size(MF,3),size(MF,4),2);
Imgs_u = zeros(1,n,n);dif=0;
Imgs_v = zeros(1,n,n);
u=zeros(n,n);v=M;
%v=zeros(n,n);u=M;
w=zeros(n,n,2);
evolution=1;i=0;
evolution_glob=1;i_glob=1;

temps = zeros(1);
temps(1)=0;
tic;
while(evolution_glob>-0.000001)
%     alpha = alpha*.9;
%     lambda = alpha*lambda0;
%     mu = alpha*mu0;

    u2_glob=u;
    v2_glob=v;
    Estim = v;
    [F_estim, Weight] = perform_windowed_fourier_transform(Estim,q,Dx,n, options);
    Orientations = estimate_orientations(F_estim,zeros(n,n),q,Dx);
    Amplitudes = ones(size(Orientations,1),size(Orientations,2));
    Chi = non_convex_weight(Orientations, Amplitudes, Weight, q);
    Amplitudes = estimate_amplitudes(F_estim, Chi);%, mask_erod, q, Dx, options);
    Chi = non_convex_weight(Orientations, Amplitudes, Weight, q);

    Orientationsss(i_glob,:,:,:) = Orientations;
    Imgs_u(i_glob,:,:) = u;
    Imgs_v(i_glob,:,:) = v;
    evolution=1;
    stop = -0.05/i_glob;
    i=0;
    while(evolution>stop && i<(1+i_glob))
        i = i + 1;
        u2=u;v2=v;

        % v is fixed: TV denoising
        for k=1:min(10,10*i_glob)
            dw = grad( (M-test_h(v,a,b))/lambda + div(w));
            w = w + tau * dw;
            d = repmat( sqrt(sum(w.^2,3)), [1 1 2] );
            w = w ./ max(d,ones(n,n,2));
        end
        u = M-test_h(v,a,b) + lambda * div(w);

        % u is fixed: gradient descent
        [v,niter] = gradient_descent_armijo(v,mu,1,.0001/i_glob,Chi,M-u,zeros(n,n),q,Dx,n,options,ptr_h,ptr_diff_h,a,b,5);


        evolution = max(max(abs(u2(:)-u(:))),max(abs(v2(:)-v(:))));
        %         Fv = perform_windowed_fourier_transform(v,q,Dx,n, options);
        %         normT = non_convex_norm(Fv,Chi);
        %         normTV = TV(u);
        %         normL2 = norm(M(:)-u(:)-v(:));
        %         normL2 = normL2*normL2;
        %         nn = mu*normT + lambda*normTV + .5 * normL2;
        %         disp([num2str(i_glob) ': ' num2str(i) ' evolution=' num2str(evolution) ' stop=' num2str(stop) ' nrj: ' num2str(nn) ' niter=' num2str(niter)]);
        disp([num2str(i_glob) ': ' num2str(i) ' evolution=' num2str(evolution) ' lambda=' num2str(lambda)]);
        imageplot({u v test_h(v,a,b) M-u-test_h(v,a,b)});
        pause(0.001);
    end
    evolution_glob = max(max(abs(u2_glob(:)-u(:))),max(abs(v2_glob(:)-v(:))));
    %energy
    Fv = perform_windowed_fourier_transform(v,q,Dx,n, options);
    normT = non_convex_norm(Fv,Chi);
    normTV = TV(u);
    normL2 = norm(M(:)-u(:)-test_h(v(:),a,b));
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
    disp(['i_glob=' num2str(i_glob) ' nrj_glob=' num2str(nrj_glob(end)) ' diff=' num2str(dif) ' evolution_glob=' num2str(evolution_glob)]);
    temps(i_glob) = toc;
end
i=i_glob-1;
close all;
u = squeeze(Imgs_u(i,:,:));
v = squeeze(Imgs_v(i,:,:));
w = M-u-test_h(v,a,b);
figure;imageplot({M u+test_h(v,a,b)},{'M' 'u+h(v)'})
figure;imageplot({u v w},{'u' 'v' 'w'})
figure;imageplot({u v test_h(v,a,b)},{'u' 'v' 'h(v)'})
snr(M0,u+test_h(v,a,b))
