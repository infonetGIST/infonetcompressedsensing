function [sol] = MDAL(F_ori,b,opts)
[tau, scale, maxit, tol, n, ~, beta, b] = Set_Parameter(opts,b);
[F,Ft] = my_linear_operators(F_ori);
rho  = tau*10;

z_cur = zeros(n,1); l_cur = zeros(n,1);
x_sol = zeros(n,1); z_sol = zeros(n,1);
if ~isfield(opts,'pt'); opts.pt = 0; end
it = 1;
t0 = cputime;
while it <= maxit
    
    tmp = Ft(b) + rho * (z_cur - l_cur);
    tmp_rho = tmp./(rho);
    
    x_cur = tmp_rho - Ft(F(tmp_rho))./(1+rho);
    x_cur = (max(-beta,min(x_cur,beta)));
    
    z_cur = hard_thresholding(x_cur+l_cur,tau,rho);
           
    l_cur = l_cur + x_cur - z_cur;

    x_sol_prv = x_sol;
    x_sol = (it+1)/(it+2)*x_sol + 1/(it+2)*x_cur;
    z_sol = (it+1)/(it+2)*z_sol + 1/(it+2)*z_cur;
    %% Termination    
    
    if opts.pt == 1
        a = norm(x_sol*scale-opts.x)/norm(opts.x);
        if a^2 <= 1e-4
            break
        end
    else
        if (norm(x_sol_prv-x_sol)/norm(x_sol) <= tol ) && it >= 150            
            break
        end
    end
    it = it + 1;
end
sol.x = x_sol*scale;
sol.x1 = x_cur*scale;
sol.iter = min(it,opts.maxit);
sol.time = cputime-t0;

end



%% tau => lambda in [ref], i.e., regularization value
%% rho => mu in [Re], i.e., penatly value
%% gamma => gamma in [Re], i.e., regularization value for inexact approahces
function sol = hard_thresholding(x,tau,rho)
std = sqrt(2*tau./(rho));
tmp = (rho.*x)./(rho);
ind1 = abs(tmp) < std;
sol = tmp;
sol(ind1) = 0;
end

