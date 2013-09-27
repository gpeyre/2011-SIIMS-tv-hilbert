function x=gradconj(b,Gamma,mu,q,Dx,n,options,tol,mask)

if (nargin<9)
    mask = zeros(n,n);
end
if (nargin<8)
    tol = 0.05;
end

x0=ones(n,n);


% initialisation
r=b-((1-mask).*x0 + 2*mu*Psi_star_Gamma2_Psi(x0,Gamma,q,Dx,n,options));
k=0;

% programme principal

xprec=x0;
rprec=r;
pprec=r;
normerprec=norm(r(:));
nor(1)=normerprec;

while normerprec>tol
    prod = (1-mask).*pprec + 2*mu*Psi_star_Gamma2_Psi(pprec,Gamma,q,Dx,n,options);
	alpha=(normerprec)^2/(sum(pprec(:).*prod(:)));
	xsuiv=xprec+alpha*pprec;
	rsuiv=rprec-alpha*prod;
	beta=(norm(rsuiv(:)))^2/(normerprec)^2;
	psuiv=rsuiv+beta*pprec;
	
	xprec=xsuiv;
	rprec=rsuiv;
	normerprec=norm(rprec(:));
	nor(k+2)=normerprec;
	pprec=psuiv;
	k=k+1;
end
x=xprec;