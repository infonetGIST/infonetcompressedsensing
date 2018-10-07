function [x] = cgls(A,b,supp)
%This function performns the conjugate gradient least square method with
%the estimated support set,i.e., supp.
%It should be noted that it solves an overdetermined system as long as the
%support set is correctly estimated.
supp = sign(abs(supp));
k = sum(supp); n = length(supp); m = length(b);
maxit = min([m,n,20]); tol = 1e-6; x = zeros(n,1);

if isnumeric(A)    
    r = b - A*(supp.*x);
    tmp = A'*(r);    
else
    r = b - A.times(supp.*x);
    tmp = A.trans(r);
end

s = zeros(k,1); s = tmp(supp==1);


% Initialize
p      = s;
norms0 = norm(s);
gamma  = norms0^2;
t      = 0; 
flag   = 0;
supp_ind = supp==1;

while (t < maxit) && (flag == 0)    
    t = t+1;
    
    tmp = zeros(n,1);
    tmp(supp_ind) = p;
    if isnumeric(A)
        q = A*(supp.*tmp);
    else
        q = A.times(supp.*tmp);
    end
           
    delta = norm(q)^2;
    if delta == 0
        delta = eps; 
    end
    alpha = gamma / delta;    
    r     = r - alpha*q;
    
    x(supp_ind) = x(supp_ind) + alpha*p;          
    
    if isnumeric(A)
        tmp = A'*r;
    else
        tmp = A.trans(r);
    end        
    s = zeros(k,1);
    s = tmp(supp_ind);
   
    norms  = norm(s);
    gamma1 = gamma;
    gamma  = norms^2;
    beta   = gamma / gamma1;
    p      = s + beta*p;
    
    % Convergence
    normx = norm(x);   
    flag  = (norms <= norms0 * tol) || (normx * tol >= 1);    
end % while