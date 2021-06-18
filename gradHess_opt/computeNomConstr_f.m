function [ f_nom] = computeNomConstr_f( x, functionParams , params)
%Compute function value, and gradient of nomination constraint
%Q_nom^(j) - deltat*Sum_i q_i^(j) - s^(j) = 0

%get parameters
q_nom = params.q_nom;
deltaT = params.deltaT;
N = params.n_well;
T = params.n_period;

%compute function value: Tx1
q_g = x(1:N*T);
s = x((3*N*T+1):end);
sum_q_g = sum(reshape(q_g,T,N),2);
f_nom = q_nom - sum_q_g - s;
end

