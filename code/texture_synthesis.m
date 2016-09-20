% synthetise une "texture" à partir d'un champ
% d'orientations/frequence/amplitudes

path(path, 'toolbox/');
path(path, 'images/');
path(path, 'images/cartoons');
options.bound = 'sym';
options.normalization = 'tightframe';
options.window_type = 'sin';

ptr_h = @h_identity;
ptr_diff_h = @diff_h_identity;
a=1;b=0;

close all;
%param gabor
n=256;
mu = 1;
lambda = 0; %0 = enlever attache donnee

q = 16;
Dx = 4;
v = zeros(n,n);v(n/2,n/2)=.01;
Tmp = perform_windowed_fourier_transform(v,q,Dx,n, options);
zzz = size(Tmp,4);
Orientations = zeros(zzz,zzz,2);
Amplitudes = ones(zzz,zzz);
for i=2:zzz-1
    for j=2:zzz-1
        Orientations(i,j,2)= 1.5;
        Orientations(i,j,1)= 1.5*sin(3*pi*i/zzz);
        Amplitudes(i,j) = 30 + i;
    end
end
niter=0;
k_glob=0;
Estim = v;
    [F_estim, Weight] = perform_windowed_fourier_transform(Estim,q,Dx,n, options);
    Chi = non_convex_weight(Orientations, Amplitudes, Weight, q);
while(1)
    k_glob = k_glob+1;
    mask_erod=zeros(n,n);

    v2 = v;
    
           
           % non convexe
   [v,niter] = gradient_descent_armijo(v,mu,lambda,-1,Chi,v,zeros(n,n),q,Dx,n,options,ptr_h,ptr_diff_h,a,b,3,mean(Amplitudes(:)));
    imageplot(v);
    pause(0.0001);
    disp(['    ' num2str(k_glob) ' : ' num2str(norm(v2-v,'fro')) ' min_ampli=' num2str(min(min(Amplitudes(3:end-3,3:end-3)))) ' (niter=' num2str(niter) ')']);
end