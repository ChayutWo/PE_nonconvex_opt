function [ f_AL] = ALagrangian_f( x, functionParams , params)
%UNTITLED3 Summary of this function goes here
%   Augmented Lagrangian function value, gradient and hessian

% get mu and lambda for augmented lagrangian function
mu = functionParams.penalty;
lambda = functionParams.lambda;

%compute value, gradient, hessian for objective function
[ f_obj, ~, ~] = computeObjGradHess( x, functionParams , params);
%compute value, gradient, hessian for constraint function
[c] = combineConst_f( x, functionParams , params );
%compute value of augmented lagrangian function
f_AL = f_obj - lambda'*c + 1/2*mu*(c'*c);
end

