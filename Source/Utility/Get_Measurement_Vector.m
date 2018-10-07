%Copyright The original package is available at http://infonet.gist.ac.kr/
%COPYRIGHT (c) 2018 Heung-No Lee, and Sangjun Park
%E-mail: heungno@gist.ac.kr, sjpark1@gist.ac.kr

function [b,sigma,F] = Get_Measurement_Vector(x,var,opts,mode)
if opts.dct == 1
    op_handle = ['@pdct_operator'];
    p = randperm(opts.n);  picks = sort(p(1:opts.m)); 
    perm = sort(randperm(opts.n));
    F = feval(eval(op_handle),picks,perm);
    b_noise_free = F.times(x);
elseif opts.orth == 1
    F = randn(opts.m,opts.n);
    F = orth(F')';
    b_noise_free = F*x;
else
    F = randn(opts.m,opts.n);
    b_noise_free = F*x;
end
if mode == 1
    sigma = var;
else
    sigma = sqrt((norm(b_noise_free,2)^2/(opts.m*10^(var/10))));    
end
b = b_noise_free + sigma*randn(opts.m,1);
end
