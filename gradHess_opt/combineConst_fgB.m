function [f_comb, g_comb, B_comb] = combineConst_fgB( x, functionParams , params )
%Get grad and hess for all constraints

%decline constraints for different decline patterns
[ f_exp, g_exp, B_exp] = computeExponentialConstr_fgB( x, functionParams , params );
[ f_hyp, g_hyp, B_hyp] = computeHyperbolicConstr_fgB( x, functionParams , params );
[ f_har, g_har, B_har] = computeHarmonicConstr_fgB( x, functionParams , params );
%time constraints
[ f_time, g_time, B_time] = computeTimeConstr_fgB( x, functionParams , params);
%nomination constraints
[ f_nom, g_nom, B_nom] = computeNomConstr_fgB( x, functionParams , params );
%combine them into constraint functions, gradients, and hessians
f_comb = vertcat(f_exp, f_hyp, f_har, f_time, f_nom);
g_comb = horzcat(g_exp, g_hyp, g_har, g_time, g_nom);
B_comb = B_exp+ B_hyp+ B_har+ B_time+ B_nom;
end

