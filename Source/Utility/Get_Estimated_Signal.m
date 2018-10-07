%Copyright The original package is available at http://infonet.gist.ac.kr/
%COPYRIGHT (c) 2018 Heung-No Lee, and Sangjun Park
%E-mail: heungno@gist.ac.kr, sjpark1@gist.ac.kr

function [tmp]=Get_Estimated_Signal(u_cur,F_ori,b,opts)
n = opts.n; m = opts.m;
ind = find(u_cur>0);
tmp = zeros(n,1);
if ~isempty(ind)
     if opts.dct
         sub_F = zeros(m,length(ind));
         x = zeros(n,1);
         [F,~] = my_linear_operators(F_ori);
         for i = 1 : size(ind)
             x(ind(i)) = 1;
             sub_F(:,i) = F(x);
             x(ind(i)) = 0;
         end
     else
         sub_F = F_ori(:,ind);
         tmp = zeros(n,1);
    end
    tmp(ind) = pinv(sub_F)*b;
end
end