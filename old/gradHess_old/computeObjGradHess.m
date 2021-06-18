function [ f_obj, g_obj, B_obj] = computeObjGradHess( x, params )
%Compute function value, and gradient of the objective
%Sum_j Sum_i {CGR_i*q_i^(j)*deltat^(j)*Price^(j)}

%get parameters
price = params.price;
CGR = params.CGR;
deltaT = params.deltaT;
N = params.n_well;
T = params.n_period;
%compute objective
f_obj = -sum(deltaT*CGR.*price.*x(1:N*T));

%compute gradient of f at x: 3NT+T x 1
g_obj = zeros(length(x),1);
g_obj(1:N*T) = -deltaT*CGR.*price;

%0 hessian: 3NT+T x 3NT+T
B_obj = zeros(length(x),length(x));
end