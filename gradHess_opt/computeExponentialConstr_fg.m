function [ f_exp, g_exp] = computeExponentialConstr_fg( x, functionParams , params)
%Compute function value, gradient, and Hessian of the Exponential constraint
%q_i^(j)*deltat - (q_poti/d_i)(exp(-d_i*t_i^(j-1)) - exp(-d_i*t_i^(j))) = 0

%parameters specific to wells with exponential decline
n_well_type = params.n_expo;
index = params.expo_index;

%get parameters
q_pot = params.q_pot(index);
d = params.decline(index);
deltaT = params.deltaT;
N = params.n_well;
T = params.n_period;
hyperbolic_const = params.hyperbolic_const;

q_g = x(index); % q_g_i^(j)
t_ij = x(N*T+index); % t_i^(j)
t_ij_prev = x(N*T+index-1); % t_i^(j-1)
for n = 1:n_well_type
    t_ij_prev(1+(n-1)*T) = 0;
end

%evaluate constraint function value
f_exp = deltaT*q_g - (q_pot./d).*(exp(-d.*t_ij_prev) - exp(-d.*t_ij));

%calculate constraint gradient AT
g_exp = zeros(length(x),n_well_type*T);

start = min(index) - 1; %shift from wells with other types
for i = 1:n_well_type
    for j = 1:T
        cur_index = (i-1)*T+j;
        q_pot_i = q_pot(cur_index);
        d_i = d(cur_index);
        t_cur = t_ij(cur_index);
        g_exp(start+cur_index,cur_index) = deltaT;
        if j == 1
            g_exp(N*T + start + cur_index,cur_index) = -q_pot_i*exp(-d_i*t_cur);
        else
            t_prev = t_ij(cur_index-1);
            g_exp(N*T + start + cur_index-1,cur_index) = q_pot_i*exp(-d_i*t_prev);
            g_exp(N*T + start + cur_index,cur_index) = -q_pot_i*exp(-d_i*t_cur);
        end
    end
end
end