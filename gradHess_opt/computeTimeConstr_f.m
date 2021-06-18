function [ f_time] = computeTimeConstr_f( x, functionParams , params)
%Compute function value, and gradient of time constraint
%t_i^(j) - t_i^(j-1) - b_i^(j) = 0

%get parameters
N = params.n_well;
T = params.n_period;

%evaluate constraint function value
f_time = zeros(N*T,1);
for n = 1:N
    t_j = x((N*T+(n-1)*T + 1):(N*T+n*T));
    t_j_1 = x((N*T+(n-1)*T):(N*T+n*T-1));
    t_j_1(1) = 0;
    b_j = x((2*N*T+(n-1)*T + 1):(2*N*T+n*T));
    f_time(((n-1)*T + 1): n*T) = t_j-t_j_1 - b_j;
end
end