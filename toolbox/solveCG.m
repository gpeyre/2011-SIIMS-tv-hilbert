function [u, niter] = solveCG(ptr_A, f, s, tol, Gamma, mu, q, Dx, n, options, mask)
% SOLVECG   Conjugate Gradients method.
%
%    Input parameters: 
%           A : Symmetric, positive definite NxN matrix 
%           f : Right-hand side Nx1 column vector 
%           s : Nx1 start vector (the initial guess)
%         tol : relative residual error tolerance for break
%               condition 
%     maxiter : Maximum number of iterations to perform
%
%    Output parameters:
%           u : Nx1 solution vector
%       niter : Number of iterations performed

% Author : Andreas Klimke, Universität Stuttgart
% Version: 1.0
% Date   : May 13, 2003
	

if (nargin<11)
    mask = zeros(n,n);
end

u = s;         % Set u_0 to the start vector s
r = f - feval(ptr_A, s, Gamma, mu, q, Dx, n, options, mask);   % Compute first residuum
p = r;         
rho = r(:)'*r(:);
niter = 0;     % Init counter for number of iterations

% Compute norm of right-hand side to take relative residuum as
% break condition.
normf = norm(f(:));
if normf < eps  % if the norm is very close to zero, take the
                % absolute residuum instead as break condition
                % ( norm(r) > tol ), since the relative
                % residuum will not work (division by zero).
  warning(['norm(f) is very close to zero, taking absolute residuum' ... 
					 ' as break condition.']);
	normf = 1;
end

while (norm(r(:))/normf > tol)   % Test break condition
	a = feval(ptr_A, p, Gamma, mu, q, Dx, n, options, mask);
	alpha = rho/(a(:)'*p(:));
	u = u + alpha*p;
	r = r - alpha*a;
	rho_new = r(:)'*r(:);
	p = r + rho_new/rho * p;
	rho = rho_new;
	niter = niter + 1;
end