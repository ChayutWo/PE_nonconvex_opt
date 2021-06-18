function [ f_exp] = computeExponentialConstr_f(x, functionParams , params)
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
end