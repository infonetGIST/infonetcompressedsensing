%Copyright The original package is available at http://infonet.gist.ac.kr/
%COPYRIGHT (c) 2018 Heung-No Lee, and Sangjun Park
%E-mail: heungno@gist.ac.kr, sjpark1@gist.ac.kr

clc; clear; close all;
addpath('./Utility');
addpath('./Other Algorithms');

%% the parameters for our algorithms
opts.x0 = [];  opts.e = 10;
opts.orth = 0; opts.dct = 1;
snr = 45;

n_list = [2^10 2^11 2^12 2^13 2^14];
iter = 1000; nrun = 50;
opts.maxit = iter; opts.tol= 1e-4; 

for i = 1 : length(n_list)
    opts.n = n_list(i);
    opts.m = floor(opts.n*0.3);
    opts.k = floor(opts.m*0.3);
    Scability_Comparison(snr,opts,nrun)
end

mse  = zeros(length(n_list),5);
sse  = zeros(length(n_list),5);
time = zeros(length(n_list),5);
iter = zeros(length(n_list),5);
rerr = zeros(length(n_list),5);

for i = 1 : length(n_list)
    n = n_list(i);
    m = floor(n*0.3);
    k = floor(m*0.3);
    tol = opts.tol;
    maxit = opts.maxit;
    aa = sprintf('N = %d, M = %d, K = %d, SNR = %d [db] tol = %d iter = %d',n,m,k,snr,tol*1e5,maxit);
    load(aa);
    rerr(i,1) = my.rerr; mse(i,1) = my.mse;  sse(i,1) = my.sse;  time(i,1) = my.time;  iter(i,1) = my.iter;
    rerr(i,2) = md.rerr; mse(i,2) = md.mse;  sse(i,2) = md.sse;  time(i,2) = md.time;  iter(i,2) = md.iter;
    rerr(i,3) = ya.rerr; mse(i,3) = ya.mse;  sse(i,3) = ya.sse;  time(i,3) = ya.time;  iter(i,3) = ya.iter;
    rerr(i,4) = or.rerr; mse(i,4) = or.mse;
end
iter = iter/nrun;
time = time/nrun;
rerr = rerr/nrun;
mse  = mse/nrun;
sse  = sse/nrun;

figure(6);
semilogy(n_list,time(:,1),'-o');
hold on; grid on;
semilogy(n_list,time(:,2),'-+');
semilogy(n_list,time(:,3),'-^');
legend('ADM-MIQP','MDAL','YALL1');
ylabel('Averaged running times (seconds)');
xlabel('The problem dimension n');
