function [f_comb] = combineConst_f( x, functionParams , params )
%Get grad and hess for all constraints

%decline constraints for different decline patterns
[ f_exp] = computeExponentialConstr_f(x, functionParams , params );
[ f_hyp] = computeHyperbolicConstr_f( x, functionParams , params );
[ f_har] = computeHarmonicConstr_f( x, functionParams , params );
%time constraints
[f_time] = computeTimeConstr_f( x, functionParams , params );
%nomination constraints
[f_nom] = computeNomConstr_f( x, functionParams , params );
%combine them into constraint functions, gradients, and hessians
f_comb = vertcat(f_exp, f_hyp, f_har, f_time, f_nom);
end

