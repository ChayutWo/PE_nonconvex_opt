function [ f_nom, g_nom] = computeNomConstr_fg( x, functionParams , params )
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

%compute gradient: 3NT+TxT (there are T constraints)
g_nom = zeros(length(x),T);
for j = 1:T
    grad = zeros(T,1);
    grad(j) = -deltaT; %partial derivative with respect to q_g_i^(j)
    grad = repmat(grad, N, 1);
    g_nom(:,j) = vertcat(grad, zeros(length(x)-N*T,1));
    g_nom(3*N*T+j,j) = -1; %partial derivative with respect to s
end
end

