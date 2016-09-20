function [v2,niter]=gradient_descent_armijo(v_init,mu,lambda,stop,Chi,y,mask,q,Dx,n,options,ptr_h,ptr_diff_h,a,b,maxiter,amplit_moy)

if (~exist('maxiter','var'))
    maxiter=100;
end
v = v_init;
v2 = v+10*stop;
rho0 = 1;
epsilon = 1e-5;
k=0;
accu=0;

while(norm(v2-v,'fro')>stop && k<maxiter)
    Fv = perform_windowed_fourier_transform(v,q,Dx,n, options);
    nrj = mu*non_convex_norm(Fv,Chi)+lambda*.5*norm((1-mask).*(y-amplit_moy*feval(ptr_h,v/amplit_moy,a,b)),'fro')^2;
    C = abs(Fv);
    C = (1 - Chi./(max(C,epsilon))) .* Fv; 
    gr = 2*perform_windowed_fourier_transform(C,q,Dx,n, options);
    gr = mu*gr - lambda*(1-mask).*(y-amplit_moy*feval(ptr_h,v/amplit_moy,a,b)).*feval(ptr_diff_h,v,a,b);

    rho = rho0;
    v_tmp = v - rho*gr;
    Fv = perform_windowed_fourier_transform(v_tmp,q,Dx,n, options);
    nrj_tmp = mu*non_convex_norm(Fv,Chi)+lambda*.5*norm((1-mask).*(y-amplit_moy*feval(ptr_h,v_tmp/amplit_moy,a,b)),'fro')^2;
    while(nrj_tmp>(nrj-1e-4*rho*sum(gr(:).^2)))
        rho = rho/2;
        v_tmp = v - rho*gr;
        Fv = perform_windowed_fourier_transform(v_tmp,q,Dx,n, options);
        nrj_tmp = mu*non_convex_norm(Fv,Chi)+lambda*.5*norm((1-mask).*(y-amplit_moy*feval(ptr_h,v_tmp/amplit_moy,a,b)),'fro')^2;
    end
    if (rho==rho0)
        accu = accu+1;
        if (accu==5)
            rho0 = 2*rho0;
            accu=0;
        end
    else
        rho0 = rho;
        accu=0;
    end
    k=k+1;
    nrj = nrj_tmp;
    v2 = v;
    v = v_tmp;
    imageplot(v);
    disp([ num2str(k) ' : ' num2str(nrj) ' , ' num2str(norm(v2-v,'fro'))]);
    pause(0.001);
end

niter=k;
v2 = v;