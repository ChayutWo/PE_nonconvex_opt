function [f_comb, g_comb] = combineConst_fg( x, functionParams , params )
%Get grad and hess for all constraints

%decline constraints for different decline patterns
[ f_exp, g_exp] = computeExponentialConstr_fg(x, functionParams , params);
[ f_hyp, g_hyp] = computeHyperbolicConstr_fg(x, functionParams , params );
[ f_har, g_har] = computeHarmonicConstr_fg( x, functionParams , params );
%time constraints
[ f_time, g_time] = computeTimeConstr_fg( x, functionParams , params );
%nomination constraints
[ f_nom, g_nom] = computeNomConstr_fg( x, functionParams , params);
%combine them into constraint functions, gradients, and hessians
f_comb = vertcat(f_exp, f_hyp, f_har, f_time, f_nom);
g_comb = horzcat(g_exp, g_hyp, g_har, g_time, g_nom);
end

