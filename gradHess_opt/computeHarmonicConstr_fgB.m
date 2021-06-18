function [ f_har, g_har, B_har] = computeHarmonicConstr_fgB( x, functionParams ,params )
%Compute function value, gradient, and Hessian of the Harmonic constraint
%q_i^(j)*deltat - (q_poti/d_i)ln((1+d_i*t_i^(j))/(1+d_i*t_i^(j-1))) = 0
%B_exp here is (lambda(i) - mu * f(i)) * hess_c(:,:,i)

% get mu and lambda for augmented lagrangian function
mu = functionParams.penalty;
lambda = functionParams.lambda;

%parameters specific to wells with harmonic decline
n_well_type = params.n_harm;
index = params.harm_index;

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
f_har = deltaT*q_g - (q_pot./d).*log((1+d.*t_ij)./(1+d.*t_ij_prev));

%calculate constraint gradient AT
g_har = zeros(length(x),n_well_type*T);
%calculate hessian: 3NT+T x 3NT+T x NT
B_har = zeros(length(x),length(x));

start = min(index) - 1; %shift from wells with other types
for i = 1:n_well_type
    for j = 1:T
        cur_index = (i-1)*T+j;
        ind = N*T + start + cur_index;
        factor = lambda(start+cur_index) - mu*f_har(cur_index);
        q_pot_i = q_pot(cur_index);
        d_i = d(cur_index);
        t_cur = t_ij(cur_index);
        g_har(start+cur_index,cur_index) = deltaT;
        if j == 1
            g_har(ind,cur_index) = -q_pot_i/(1+d_i*t_cur);
            B_har(ind,ind) = B_har(ind,ind)+factor*q_pot_i*d_i/(1+d_i*t_cur)^2;
        else
            t_prev = t_ij(cur_index-1);
            g_har(ind-1,cur_index) = q_pot_i/(1+d_i*t_prev);
            g_har(ind,cur_index) = -q_pot_i/(1+d_i*t_cur);
            B_har(ind-1,ind-1) = B_har(ind-1,ind-1)-factor*q_pot_i*d_i/(1+d_i*t_prev)^2;
            B_har(ind,ind) = B_har(ind,ind)+factor*q_pot_i*d_i/(1+d_i*t_cur)^2;
        end
    end
end
end