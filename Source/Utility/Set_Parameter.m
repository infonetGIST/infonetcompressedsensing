%Copyright The original package is available at http://infonet.gist.ac.kr/
%COPYRIGHT (c) 2018 Heung-No Lee, and Sangjun Park
%E-mail: heungno@gist.ac.kr, sjpark1@gist.ac.kr

function [tau, scale, maxit, tol, n, rho, beta, b] = Set_Parameter(opts,b_ori)
maxit = opts.maxit; tol = opts.tol;
scale = norm(b_ori,inf);
n = opts.n; m = opts.m;
b = b_ori/scale;
beta = (opts.xmax)/scale;
tau = opts.tau/scale;
if opts.dct == 0 && opts.orth == 0
    rho = 1;       
else    
    rho = tau/beta;
end
end

