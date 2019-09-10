%rouwenhorst.m
%
%[zgrid, P] = rouwenhorst(rho, sigma_eps, n)
%
% rho is the 1st order autocorrelation
% sigma_eps is the standard deviation of the error term
% n is the number of points in the discrete approximation
%
function f = KLS2_rouwenhorst_fn(x,rho,sigma_eps,mu,N)


p = x(1); 
q = x(2); 
nu = x(3); 

sigma_z = sigma_eps/sqrt(1-rho^2);
s = (1-q)/(2-(p+q));

f(1) = (p + q - 1) - rho;
f(2) = (q-p)*nu/(2-(p+q)) - mu;
f(3) = nu^2*(1-4*s*(1-s)+4*s*(1-s)/(N-1)) - (sigma_z^2 + mu^2);


end



