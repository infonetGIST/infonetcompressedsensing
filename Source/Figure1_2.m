%Copyright The original package is available at http://infonet.gist.ac.kr/
%COPYRIGHT (c) 2018 Heung-No Lee, and Sangjun Park
%E-mail: heungno@gist.ac.kr, sjpark1@gist.ac.kr

clear; clc; close all;
addpath('./Utility'); addpath('./Other Algorithms');

opts.x0 = [];  opts.e = 10;
opts.orth = 0; opts.dct = 1;
iter = 2000; nrun = 30;
opts.maxit = iter;
snr = 45;
opts.tol = 0; opts.n = 1024; opts.m = 307; opts.k = 30;
my_MSE= zeros(2000,1);
my_SSE= zeros(2000,1);
for j = 1 : nrun
    [x,true_support,pp] = Signal_Generation(opts);
    [b,sigma,F] = Get_Measurement_Vector(x,snr,opts,0);
    opts.xmax = max(abs(x)); opts.x0 = x;
    if snr == inf
        opts.tau = 1e-4;
    else
        opts.tau = sigma * sqrt(log(2*opts.n));
    end
    my_MSE  = my_MSE + Non_Convex_ADM_MSE(F,b,opts,true_support,pp);
    my_SSE  = my_SSE + Non_Convex_ADM_SSE(F,b,opts,true_support,pp); 
end

mse = my_MSE/nrun;
figure(1);
semilogy(mse,'-o'); hold on; legend('ADM-MIQP');
xlabel('The number of iterations'); ylabel('Averaged mean square error (MSE)'); grid on;

sse = my_SSE/nrun;
figure(2);
semilogy(sse,'-o'); hold on; legend('ADM-MIQP');
xlabel('The number of iterations'); ylabel('Averaged support set error (MSE)'); grid on;
