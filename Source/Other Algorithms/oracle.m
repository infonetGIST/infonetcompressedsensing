function sol = oracle(F_ori,b_ori,true_support,opts)
if opts.dct || opts.n >= 1e4
    sol.x = cgls(F_ori,b_ori,true_support);  
else    
    sol.x = Get_Estimated_Signal(true_support,F_ori,b_ori,opts);
end
end
