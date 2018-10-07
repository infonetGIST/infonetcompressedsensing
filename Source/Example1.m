%Copyright The original package is available at http://infonet.gist.ac.kr/
%COPYRIGHT (c) 2018 Heung-No Lee, and Sangjun Park
%E-mail: heungno@gist.ac.kr, sjpark1@gist.ac.kr

addpath('./Utility');
addpath('./Other Algorithms');
clear; close all;

opts.n = 512*4;
opts.m = floor(0.3*opts.n);
opts.k = floor(0.35*opts.m);
snr = 60;
opts.maxit = 1000; opts.x0  = [];
opts.tol   = 1e-4; opts.e   = 10;
opts.orth  = 0; opts.dct = 0;

[x,true_support] = Signal_Generation(opts);
[b,sigma,F] = Get_Measurement_Vector(x,snr,opts,0);
opts.tau = sigma * sqrt(log(2*opts.n));
opts.xmax = 1.5*max(abs(x));
sol_my   = Non_Convex_ADM(F,b,opts);

a    = norm(sol_my.x-x);
mse  = a^2/opts.n;
sse  = sum(norm(true_support-sol_my.u,1))/opts.k;
rerr = a/norm(x);
time = sol_my.time;

figure(1); 
subplot(2,1,1); plot(x,'-o'); grid on; hold on; xlabel('index'); ylabel('coefficients'), plot(sol_my.x,'-*'); legend('original x','x estimated by ADM-MIQP');
subplot(2,1,2); plot(x-sol_my.x,'-'); grid on; xlabel('index'); ylabel('x-x_{sol}');
fprintf('\nN = %d, M = %d, K = %d and SNR [dB] = %d\n',opts.n,opts.m,opts.k,snr);
fprintf('MSE (mean square error) = %0.8f\t\tSSE (support set error) = %0.8f\t\tRerr (relative error) = %0.8f\n',mse,sse,rerr);
fprintf('Reconstruction time = %0.2f (seconds)\n\n',time);
