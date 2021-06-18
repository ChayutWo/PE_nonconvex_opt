function [ f_hyp] = computeHyperbolicConstr_f( x, functionParams , params)
%Compute function value, gradient, and Hessian of the Hyperbolic constraint
%q_i^(j)*deltat - (q_poti/(b-1)d_i)((1+b*d_i*t_i^(j))^(1-1/b)-(1+b*d_i*t_i^(j-1))^(1-1/b))

%parameters specific to wells with hyperbolic decline
n_well_type = params.n_hyper;
index = params.hyper_index;

%get parameters
q_pot = params.q_pot(index);
d = params.decline(index);
deltaT = params.deltaT;
N = params.n_well;
T = params.n_period;
hyperbolic_const = params.hyperbolic_const;
b = hyperbolic_const;

q_g = x(index); % q_g_i^(j)
t_ij = x(N*T+index); % t_i^(j)
t_ij_prev = x(N*T+index-1); % t_i^(j-1)
for n = 1:n_well_type
    t_ij_prev(1+(n-1)*T) = 0;
end

%evaluate constraint function value
f_hyp = deltaT*q_g - (q_pot./(d*(b-1))).*((1+b*d.*t_ij).^(1-1/b)-(1+b*d.*t_ij_prev).^(1-1/b));
end