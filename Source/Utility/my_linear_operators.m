function [A,At] = my_linear_operators(A0)
if isnumeric(A0);
    A  = @(x) A0*x;
    At = @(y) (y'*A0)';
elseif isstruct(A0) && isfield(A0,'times') && isfield(A0,'trans');
    A  = A0.times;
    At = A0.trans;
elseif isa(A0,'function_handle')
    A  = @(x) A0(x,1);
    At = @(x) A0(x,2);
else
    A = @(x) x;
    At = @(x) x;
end