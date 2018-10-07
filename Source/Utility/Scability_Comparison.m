%Copyright The original package is available at http://infonet.gist.ac.kr/
%COPYRIGHT (c) 2018 Heung-No Lee, and Sangjun Park
%E-mail: heungno@gist.ac.kr, sjpark1@gist.ac.kr

function Scability_Comparison(snr,opts,nrun)
addpath('./Utility');
addpath('./Other Algorithms');
aa = sprintf('N = %d, M = %d, K = %d, SNR = %d [db] tol = %d iter = %d',opts.n,opts.m,opts.k,snr,opts.tol*1e5,opts.maxit);

my_sse  = 0; md_sse  = 0; ya_sse  = 0;             
my_mse  = 0; md_mse  = 0; ya_mse  = 0; or_mse  = 0; 
my_rerr = 0; md_rerr = 0; ya_rerr = 0; or_rerr = 0; 
my_time = 0; md_time = 0; ya_time = 0;           
my_iter = 0; md_iter = 0; ya_iter = 0;
for j = 1 : nrun
    n = opts.n;
    [x,true_support,pp] = Signal_Generation(opts);
    [b,sigma,F] = Get_Measurement_Vector(x,snr,opts,0);
    opts.xmax = max(abs(x));
    if snr == inf
        opts.tau = 1e-4;
    else
        opts.tau = sigma * sqrt(log(2*opts.n));
        opts.rho = opts.tau;
    end
    
    x_norm = norm(x);
    thres_value = min(abs(x(pp)))*0.8;
    
    sol_my   = Non_Convex_ADM(F,b,opts);
    sol_md   = MDAL(F,b,opts);
    sol_ya   = yall1(F,b,opts);
    sol_or   = oracle(F,b,true_support,opts);
    
    sol_md.u = zeros(n,1);
    tmp = zeros(n,1);
    tmp2 = abs(sol_md.x);
    ind = find(tmp2~=0 & tmp2>=thres_value);
    tmp(ind) = ones(length(ind),1);
    sol_md.u = tmp;
    
    
    sol_ya.u = zeros(n,1);
    tmp = zeros(n,1);
    tmp2 = abs(sol_ya.x);
    ind = find(tmp2~=0 & tmp2>=thres_value);
    tmp(ind) = ones(length(ind),1);
    sol_ya.u = tmp;
    
    a = norm(sol_my.x-x);
    b = norm(sol_or.x-x);
    c = norm(sol_md.x-x);
    d = norm(sol_ya.x-x);
    
    my_mse  =   my_mse   + a^2/opts.n;
    my_rerr =   my_rerr  + a/x_norm;
    
    my_time =   my_time  + sol_my.time;
    my_iter =   my_iter  + sol_my.iter;
    my_sse  =   my_sse   + sum(norm(true_support-sol_my.u,1))/opts.k;
    
    or_mse  = or_mse  + b^2/opts.n;
    or_rerr = or_rerr + b/x_norm;
    
    md_mse  = md_mse  + c^2/opts.n;
    md_rerr = md_rerr + c/x_norm;
    md_time = md_time + sol_md.time;
    md_iter = md_iter + sol_md.iter;
    md_sse  = md_sse  + sum(norm(true_support-sol_md.u,1))/opts.k;
    
    ya_mse  = ya_mse  + d^2/opts.n;
    ya_rerr = ya_rerr + d/x_norm;
    ya_time = ya_time + sol_ya.time;
    ya_iter = ya_iter + sol_ya.iter;
    ya_sse  = ya_sse  + sum(norm(true_support-sol_ya.u,1))/opts.k;
end
my.mse = 0; my.rerr = 0; my.time = 0; my.iter = 0; my.sse = 0; 
md.mse = 0; md.rerr = 0; md.time = 0; md.iter = 0; md.sse = 0;
ya.mse = 0; ya.rerr = 0; ya.time = 0; ya.iter = 0; ya.sse = 0;
or.mse = 0; or.rerr = 0; 

my.sse  = my_sse;  my.mse  = my_mse;  my.rerr = my_rerr;   my.time = my_time;  my.iter = my_iter;
md.sse  = md_sse;  md.mse  = md_mse;  md.rerr = md_rerr;   md.time = md_time;  md.iter = md_iter;
ya.sse  = ya_sse;  ya.mse  = ya_mse;  ya.rerr = ya_rerr;   ya.time = ya_time;  ya.iter = ya_iter;
or.mse  = or_mse;  or.rerr = or_rerr;

disp(aa);
disp(sprintf('my_sse = %0.9f   my_mse = %0.9f   my_rerr = %0.9f    my_time = %f   my_iter = %f,   ',my.sse/nrun,  my.mse/nrun,   my.rerr/nrun,  my.time/nrun,  my.iter/nrun));
disp(sprintf('md_sse = %0.9f   md_mse = %0.9f   md_rerr = %0.9f    md_time = %f   md_iter = %f,   ',md.sse/nrun,  md.mse/nrun,   md.rerr/nrun,  md.time/nrun,  md.iter/nrun));
disp(sprintf('ya_sse = %0.9f   ya_mse = %0.9f   ya_rerr = %0.9f    ya_time = %f   ya_iter = %f,   ',ya.sse/nrun,  ya.mse/nrun,   ya.rerr/nrun,  ya.time/nrun,  ya.iter/nrun));
disp(sprintf('\t\t\t\t\t   or_mse = %0.9f   or_rerr = %0.9f\n',                                                          or.mse/nrun,   or.rerr/nrun                              ));
save(aa,'opts','my','md','ya','or','nrun');
end
