%Copyright The original package is available at http://infonet.gist.ac.kr/
%COPYRIGHT (c) 2018 Heung-No Lee, and Sangjun Park
%E-mail: heungno@gist.ac.kr, sjpark1@gist.ac.kr

function [x,true_support,pp] = Signal_Generation(opts)
p = randperm(opts.n);
x = zeros(opts.n,1);
x(p(1:opts.k)) = (randn(opts.k,1));
x = (opts.e*x/norm(x));   
true_support = zeros(opts.n,1);
true_support(p(1:opts.k)) = 1;
pp = p(1:opts.k);
end