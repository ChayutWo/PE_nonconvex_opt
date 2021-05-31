function [f_comb, g_comb, B_comb] = combineConst( x, params )
%Get grad and hess for all constraints
[ f_har, g_har, B_har] = computeHarmonicConstr( x, params );
[ f_time, g_time, B_time] = computeTimeConstr( x, params );
[ f_nom, g_nom, B_nom] = computeNomConstr( x, params );
f_comb = vertcat(f_har, f_time, f_nom);
g_comb = horzcat(g_har, g_time, g_nom);
B_comb = cat(3,B_har, B_time, B_nom);

end

