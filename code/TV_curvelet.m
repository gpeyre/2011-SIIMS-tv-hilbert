% Inpainting TV+Curvelet (pour cartoon)
path(path, 'toolbox/');
path(path, 'images/');
path(path, 'images/cartoons');
close all;

n=256;
%M0 = rescale(load_image('lena',n));
M0 = rescale(load_image('bugsbunny_carreNB',n));
%M0 = rescale(load_image('test_tvcurv',n)); M0 = squeeze(M0(:,:,1));

sigma = 0;
M = M0 + sigma*randn(n,n);

%mask = load_image('mask_text_512x512',n);mask=squeeze(mask(:,:,1));mask(mask<128)=0;mask(mask~=0)=1;%mask = 1-mask;
%mask = load_image('mask_lena',n);mask=squeeze(mask(:,:,1));mask(mask<128)=0;mask(mask~=0)=1;mask = 1-mask;
mask = load_image('mask2',n);mask=squeeze(mask(:,:,1));mask(mask<128)=0;mask(mask~=0)=1;mask = 1-mask;
%mask = load_image('mask_test_tvcurv',n);mask=squeeze(mask(:,:,1));mask(mask<128)=0;mask(mask~=0)=1;mask = 1-mask;
mask=mmm;
y = (1-mask).*M;
y(y==0)=.5;
%y = uuu;

%attache aux données:
g = y;


u = g;
lambda_max = .1;
lambda_min = .01;
niter = 100;
lambda_list = linspace(lambda_max,lambda_min, niter);
gam =  .75; % TV:gam=1 Curvelet:gam=0

mu = 1/2;
tau = 1/4;

options.n=n;
options.is_real = 1;
options.finest = 1;
coarsest_scale = 2;
options.nbscales = log2(n)-coarsest_scale;

p1=zeros(n,n,2);Qstar_p=zeros(n,n);
p2 = perform_curvelet_transform(y,options);
nrjs = zeros(1);k_glob=0;
tmps = zeros(1);
tic;
while(1)
    u2=u;
    k_glob=k_glob+1;
    if (k_glob<niter)
        lambda = lambda_list(k_glob);
    else
        %break;
    end
    v = u + mu*(1-mask).*(g-u);
    k_qstar=0;
    while(k_qstar==0 || (norm(Qstar_p2(:)-Qstar_p(:))>1))
        Qstar_p2 = Qstar_p;
        Qstar_p = -gam*div(p1);
        if (gam~=1)
            Qstar_p = Qstar_p + (1-gam)*perform_curvelet_transform(p2,options);
        end
%         norm(Qstar_p(:)-Qstar_p2(:))
%         imageplot(Qstar_p);pause(0.01);

        q1 = gam*grad(Qstar_p - v/(lambda*mu));
        p1 = p1 - tau * q1;
        d = repmat( sqrt(sum(p1.^2,3)), [1 1 2] );
        p1 = p1 ./ max(d,ones(n,n,2));

        if (gam~=1)
            q2 = perform_curvelet_transform(Qstar_p- v/(lambda*mu),options);
            for i=1:size(q2,1)
                for j=1:size(q2,2)
                    for k=1:size(q2{i,j},1)
                        for l=1:size(q2{i,j},2)
                            p2{i,j}{k,l} = p2{i,j}{k,l} - tau * (1-gam) * q2{i,j}{k,l};
                            p2{i,j}{k,l} = p2{i,j}{k,l}./max(ones(size(p2{i,j}{k,l})),abs(p2{i,j}{k,l}));
                        end
                    end
                end
            end
        end
        k_qstar = k_qstar+1;
    end
    Qstar_p = -gam*div(p1);
    if (gam~=1)
        Qstar_p = Qstar_p + (1-gam)*perform_curvelet_transform(p2,options);
    end
    u = v - lambda * mu * Qstar_p;


    nrj = norm(g-(1-mask).*u,'fro');
    nrj = 1/(2*lambda) * nrj*nrj + gam*TV(u);
    if (gam~=1)
        nrj = nrj + (1-gam)*L1_curv(u,options);
    end
    nrjs(k_glob)=nrj;
    tmps(k_glob)=toc;
    disp(['k_glob:' num2str(k_glob) ' nrj:' num2str(nrj) ' : ' num2str(norm(u2-u,'fro')) ' lambda:' num2str(lambda) '(' num2str(k_qstar) ')']);
    %     if (k_glob~=1 && nrjs(k_glob)>nrjs(k_glob-1))
    %         disp('----------------------------------');
    %         break;
    %     end
    imageplot({g u M0},{'M' num2str(snr(M0,u)) 'M0'});pause(0.001);
end

