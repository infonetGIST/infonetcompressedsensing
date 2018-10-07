%Copyright The original package is available at http://infonet.gist.ac.kr/
%COPYRIGHT (c) 2018 Heung-No Lee, and Sangjun Park
%E-mail: heungno@gist.ac.kr, sjpark1@gist.ac.kr

addpath('./Utility');
addpath('./Other Algorithms');
clc; clear; close all;
opts.n = 1024;
compression_ratio =  0.15:0.025:0.35;
sparsity_ratio =  0.15:0.025:0.35;

snr = inf;
opts.maxit = 10000; opts.x0  = [];
opts.tol   = 1e-4; opts.e   = 10;
opts.orth  = 0;    opts.dct = 1;
nrun = 100;

Number = zeros(length(sparsity_ratio), length(compression_ratio),3);
n = opts.n;
for cnt = 1 : nrun
    for i = 1 : length(compression_ratio)
        m = floor(n*compression_ratio(i));
        for j = 1 : length(sparsity_ratio)
            k = floor(m*sparsity_ratio(j));
            a = zeros(3,1);
            opts.n = n; opts.m = m; opts.k = k; opts.pt = 0;
            [x,true_support] = Signal_Generation(opts);
            [b,sigma,F] = Get_Measurement_Vector(x,snr,opts,0);
            if snr == inf
                opts.tau = 1e-4;      % YALL1 solves the basis pursuit problem
            else
                opts.tau = sigma * sqrt(log(2*opts.n));
                opts.rho = opts.tau;  % the regularization parameter for YALL! is referred to as rho
            end
            opts.xmax = max(abs(x));   opts.x = x;
            sol_my   = Non_Convex_ADM(F,b,opts);
            
            opts.x = x;
            opts.norm_x0 = norm(opts.x);
            opts.pt = 1;
            sol_md   = MDAL(F,b,opts);
            sol_ya   = yall1(F,b,opts);
            a(1) = norm(sol_my.x-x)/norm(x);
            a(2) = norm(sol_md.x-x)/norm(x);
            a(3) = norm(sol_ya.x-x)/norm(x);
            for l = 1 : 3
                if a(l)^2 <= 1e-4
                    Number(j,i,l) = Number(j,i,l) + 1;
                end
            end
        end
    end
    aa = sprintf('PT_data');
    save(aa,'cnt','nrun','compression_ratio','sparsity_ratio','snr','Number');
end
A1 = Number(:,:,1)/nrun;
A2 = Number(:,:,2)/nrun;
A3 = Number(:,:,3)/nrun;

figure1 = figure();
axes1 = axes('Parent',figure1); hold(axes1,'on');
imagesc(compression_ratio,sparsity_ratio,A1,'Parent',axes1,...
    'CDataMapping','scaled');
grid on; colorbar('peer',axes1);
xlabel('Under sampling ration m/n'); ylabel('Over sampling ratio k/m');
title('Phase transition of ADM-MIQP');

figure1 = figure();
axes1 = axes('Parent',figure1); hold(axes1,'on');
imagesc(compression_ratio,sparsity_ratio,A2,'Parent',axes1,...
    'CDataMapping','scaled');
box(axes1,'on'); grid on; colorbar('peer',axes1);
xlabel('Under sampling ration m/n'); ylabel('Over sampling ratio k/m');
title('Phase transition of MDAL');

figure1 = figure();
axes1 = axes('Parent',figure1); hold(axes1,'on');
imagesc(compression_ratio,sparsity_ratio,A3,'Parent',axes1,...
    'CDataMapping','scaled');
box(axes1,'on'); grid on; colorbar('peer',axes1);
xlabel('Under sampling ration m/n'); ylabel('Over sampling ratio k/m');
title('Phase transition of YALL1');
