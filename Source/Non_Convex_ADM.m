%Copyright The original package is available at http://infonet.gist.ac.kr/
%COPYRIGHT (c) 2018 Heung-No Lee, and Sangjun Park
%E-mail: heungno@gist.ac.kr, sjpark1@gist.ac.kr

function [sol] = Non_Convex_ADM(F_ori,b_ori,opts)
[tau, scale,maxit, tol, n, rho, beta, b] = Set_Parameter(opts,b_ori);
[F,Ft] = my_linear_operators(F_ori);
Ftb = Ft(b);

dn = 2*n; it = 1;

if ~isfield(opts,'x_init'); opts.x_init = zeros(n,1) ; end
if ~isfield(opts,'u_init'); opts.u_init = zeros(n,1) ; end
if ~isfield(opts,'z_init'); opts.z_init = zeros(dn,1); end
if ~isfield(opts,'l_init'); opts.l_init = zeros(dn,1); end
if opts.orth == 0 && opts.dct == 0
    D = inv(F_ori*F_ori' + 2*rho*eye(opts.m));
end

x_cur = opts.x_init; u_cur = opts.u_init;
z_cur = opts.z_init; l_cur = opts.l_init;

zero_vec = zeros(n,1);
zero_vecz = zeros(dn,1);

q = [-Ftb ;tau.*ones(n,1)];
rho1 = rho;
r_beta_square = -(rho1*beta^2).*ones(n,1);
t0 = cputime;
while it <= maxit     
    tmp = beta.*u_cur;
    tmp1 = [x_cur-tmp ; -x_cur-tmp ];    
    z_cur = max(zero_vecz,l_cur./rho - tmp1);
       
    
    x_prv = x_cur;
    tmp = l_cur-rho.*z_cur ;
    q1 = q - [ tmp(1:n)-tmp(n+1:dn) ; -beta.*(tmp(1:n)+tmp(n+1:dn))];
    q1_rho = q1(1:n)./(2*rho);

    if opts.dct == 0 && opts.orth == 0
        x_cur = Ft(D*F(q1_rho)) - q1_rho;
    else
        x_cur = (-q1_rho + Ft(F(q1_rho))./(1+2*rho));  
    end
    
    
    u_prv = u_cur;
    u_cur = zero_vec;    
    u_cur(q1(n+1:dn)<=r_beta_square) = 1;     
    
    tmp = beta.*u_cur;    
    tmp1 = [x_cur-tmp ; -x_cur-tmp ];    
    l_cur = l_cur - rho.*(tmp1+z_cur);           
    if (norm([x_cur;u_cur]-[x_prv;u_prv])/norm([x_cur;u_cur]) <= tol ) && it >= 50  
        break
    end    
    it = it + 1;               
end
sol.x1 = x_cur*scale;
sol.u = u_cur;   
sol.x = cgls(F_ori,b_ori,u_cur);  
sol.iter = min(it,maxit);
sol.time = cputime-t0;
end