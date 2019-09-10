%rouwenhorst.m
%
%[zgrid, P] = rouwenhorst(rho, sigma_eps, n)
%
% rho is the 1st order autocorrelation
% sigma_eps is the standard deviation of the error term
% n is the number of points in the discrete approximation
%
function [zgrid, P, pi] = KLS2_rouwenhorst(rho,sigma_eps,mu,n)

options = optimset('Display','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',1e8,'MaxFunEval',1e8); %,'Algorithm','trust-region-reflective');
x = fsolve(@(x)KLS2_rouwenhorst_fn(x,rho,sigma_eps,mu,n),[0.5 0.5 0.5],options);

p = x(1); 
q = x(2); 
nu = x(3); 

mu_eps = 0;

P = [p 1-p;1-q q];

for i=2:n-1
   P = p*[P zeros(i,1);zeros(1,i+1)] + (1-p)*[zeros(i,1) P;zeros(1,i+1)] + ...
       (1-q)*[zeros(1,i+1); P zeros(i,1)] + q*[zeros(1,i+1); zeros(i,1) P];
   P(2:i,:) = P(2:i,:)/2;
end

zgrid = linspace(mu_eps/(1-rho)-nu,mu_eps/(1-rho)+nu,n);

evalc('[V,D]=eigs(P'',1)');
pi=V/sum(V);



