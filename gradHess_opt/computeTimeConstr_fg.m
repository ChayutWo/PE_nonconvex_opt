function [ f_time, g_time] = computeTimeConstr_fg( x, functionParams , params )
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

%calculate constraint gradient AT
g_time = zeros(length(x),N*T);
for i = 1:N
    for j = 1:T
        if j == 1
            g_time(N*T+T*(i-1) + j,((i-1)*T+j)) = 1; % wrt t_i^(j)
            g_time(2*N*T+T*(i-1) + j,((i-1)*T+j)) = -1; % wrt b
        else
            g_time(N*T+T*(i-1) + j-1,((i-1)*T+j)) = -1; % wrt t_i^(j-1)
            g_time(N*T+T*(i-1) + j,((i-1)*T+j)) = 1; % wrt t_i^(j)
            g_time(2*N*T+T*(i-1) + j,((i-1)*T+j)) = -1; % wrt b
        end
    end
end
end