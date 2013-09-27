function y=Psi_star_Gamma2_Psi(x,Gamma,q,Dx,n,options)

Mx = perform_windowed_fourier_transform(x,q,Dx,n, options);

Tmp = Gamma.*Gamma.*Mx;

y = perform_windowed_fourier_transform(Tmp,q,Dx,n, options);

