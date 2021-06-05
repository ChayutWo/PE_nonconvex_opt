function [ f_har, g_har, B_har] = computeHarmonicConstr( x, params )
%Compute function value, gradient, and Hessian of the Harmonic constraint
%q_i^(j)*deltat - (q_poti/d_i)ln((1+d_i*t_i^(j))/(1+d_i*t_i^(j-1))) = 0

%get parameters
q_pot = params.q_pot;
deltaT = params.deltaT;
d = params.decline;
N = params.n_well;
T = params.n_period;
q_g = x(1:N*T);
t_ij = x((N*T+1):2*N*T);
t_ij_prev = x((N*T):2*N*T-1);
for n = 1:N
    t_ij_prev(1+(n-1)*T) = 0;
end

%evaluate constraint function value
f_har = deltaT*q_g - (q_pot./d).*log((1+d.*t_ij)./(1+d.*t_ij_prev));

%calculate constraint gradient AT
g_har = zeros(length(x),N*T);
%calculate hessian: 3NT+T x 3NT+T x NT
B_har = zeros(length(x),length(x), N*T);
for i = 1:N
    for j = 1:T
        cur_index = (i-1)*T+j;
        q_pot_i = q_pot(cur_index);
        d_i = d(cur_index);
        t_cur = t_ij(cur_index);
        g_har(cur_index,cur_index) = deltaT;
        if j == 1
            g_har(N*T + cur_index,cur_index) = -q_pot_i/(1+d_i*t_cur);
            B_har(N*T + cur_index,N*T +cur_index,cur_index) = q_pot_i*d_i/(1+d_i*t_cur)^2;
        else
            t_prev = t_ij(cur_index-1);
            g_har(N*T + cur_index-1,cur_index) = q_pot_i/(1+d_i*t_prev);
            g_har(N*T + cur_index,cur_index) = -q_pot_i/(1+d_i*t_cur);
            B_har(N*T + cur_index-1,N*T +cur_index-1,cur_index) = -q_pot_i*d_i/(1+d_i*t_prev)^2;
            B_har(N*T + cur_index,N*T +cur_index,cur_index) = q_pot_i*d_i/(1+d_i*t_cur)^2;
        end
    end
end
end