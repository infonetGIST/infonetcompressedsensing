%Copyright The original package is available at http://infonet.gist.ac.kr/
%COPYRIGHT (c) 2018 Heung-No Lee, and Sangjun Park
%E-mail: heungno@gist.ac.kr, sjpark1@gist.ac.kr

clc; clear; close all;
addpath('./Utility');
addpath('./Other Algorithms');

opts.x0 = [];  opts.e = 10;
opts.orth = 0; opts.dct = 1;
k_list = [30 40 50 60 70 80 90 100 110];
iter = 2000; nrun = 30;
opts.maxit = iter; snr = 45; opts.tol= 1e-4; 
opts.n = 1024;
opts.m = 307;

for i = 1 : length(k_list)
    opts.k = k_list(i);
    Scability_Comparison(snr,opts,nrun);
end

mse  = zeros(length(k_list),5); sse  = zeros(length(k_list),5);
time = zeros(length(k_list),5); iter = zeros(length(k_list),5);
rerr = zeros(length(k_list),5);

for i = 1 : length(k_list)
    k = k_list(i);
    aa = sprintf('N = %d, M = %d, K = %d, SNR = %d [db] tol = %d iter = %d',opts.n,opts.m,k,snr,opts.tol*1e5,opts.maxit);
    load(aa);
    rerr(i,1) = my.rerr; mse(i,1) = my.mse;  sse(i,1) = my.sse;  time(i,1) = my.time;  iter(i,1) = my.iter;
    rerr(i,2) = md.rerr; mse(i,2) = md.mse;  sse(i,2) = md.sse;  time(i,2) = md.time;  iter(i,2) = md.iter;
    rerr(i,3) = ya.rerr; mse(i,3) = ya.mse;  sse(i,3) = ya.sse;  time(i,3) = ya.time;  iter(i,3) = ya.iter;
    rerr(i,4) = or.rerr; mse(i,4) = or.mse;
end

sse  =  sse/nrun; mse  =  mse/nrun; rerr = rerr/nrun; time = time/nrun; iter = iter/nrun;
figure(4);
semilogy(k_list,mse(:,1),'-o'); hold on; semilogy(k_list,mse(:,2),'-*'); semilogy(k_list,mse(:,3),'-+'); semilogy(k_list,mse(:,4),'-');
legend('ADM-MIQP','MDAL','YALL1','ORACLE');
xlabel('sparsity level k'); ylabel('Averaged mean square error (MSE)'); grid on;

figure(5);
semilogy(k_list,sse(:,1),'-o'); hold on; semilogy(k_list,sse(:,2),'-*'); semilogy(k_list,sse(:,3),'-+');  
legend('ADM-MIQP','MDAL','YALL1'); 
xlabel('sparsity level k'); ylabel('Averaged support set error (SSE)'); grid on;
