function [KKT_error] = computeKKT_AL(x,functionParams,params,l,u)
%Check KKT condition for AL method with bound constraints
%   ||x - P(x - gradAL,l,u)||

%compute gradient of AL function at x
[ ~, grad_AL, ~ ] = ALagrangian( x, functionParams , params);
%compute KKT error
KKT_error = norm(x - project(x-grad_AL,l,u));
end

