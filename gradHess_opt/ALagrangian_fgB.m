function [ f_AL, grad_AL, hess_AL ] = ALagrangian_fgB( x, functionParams , params)
%UNTITLED3 Summary of this function goes here
%   Augmented Lagrangian function value, gradient and hessian

% get mu and lambda for augmented lagrangian function
mu = functionParams.penalty;
lambda = functionParams.lambda;

%compute value, gradient, hessian for objective function
[ f_obj, dfdx, hess_f] = computeObjGradHess( x, functionParams , params );
%compute value, gradient, hessian for constraint function
[c, gradC, hess_c] = combineConst_fgB( x, functionParams , params );
%compute value of augmented lagrangian function
f_AL = f_obj - lambda'*c + 1/2*mu*(c'*c);

%compute gradient of augmented lagrangian function
grad_AL = dfdx - gradC*(lambda - mu*c);

%compute hessian of augmented lagrangian function
hess_AL = hess_f +  mu * (gradC*gradC') - hess_c;
end

