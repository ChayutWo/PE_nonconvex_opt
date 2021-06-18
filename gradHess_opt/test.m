clc
clear all
rng(3) %set random seed
option = 1; %the case that we want to test

%generate a simple setting and randomly set initialization
if option == 1
    %case1: Infeasible start in a small field to check the optimization
    N_exponential = 10;
    N_hyperbolic = 10;
    N_harmonic = 10;
    N = N_exponential + N_hyperbolic + N_harmonic; %number of wells
    T = 50; %number of time steps
    [x,functionParams,params,l,u] = gen_case_1(N_exponential, N_hyperbolic, N_harmonic, T);
end


%eval = @computeTimeConstr;
%eval = @computeNomConstr;
%eval = @computeExponentialConstr;
%eval = @computeHyperbolicConstr;
%eval = @computeHarmonicConstr;
%eval = @computeObjGradHess;
eval = @combineConst;

iter = 1;
T_eval = zeros(iter,1);
T_comb_f = zeros(iter,1);
T_comb_g = zeros(iter,1);
T_comb_B = zeros(iter,1);
for i = 1:100
    tic
    [ f_exp, g_exp, B_exp] = computeExponentialConstr_fgB( x, functionParams , params );
    [ f_hyp, g_hyp, B_hyp] = computeHyperbolicConstr_fgB( x, functionParams , params );
    [ f_har, g_har, B_har] = computeHarmonicConstr_fgB( x, functionParams , params );
    %time constraints
    [ f_time, g_time, B_time] = computeTimeConstr_fgB( x, functionParams , params);
    %nomination constraints
    [ f_nom, g_nom, B_nom] = computeNomConstr_fgB( x, functionParams , params );
    T_eval(i) = toc;

    %combine them into constraint functions, gradients, and hessians
    tic
    f_comb = vertcat(f_exp, f_hyp, f_har, f_time, f_nom);
    T_comb_f(i) = toc;
    tic
    g_comb = horzcat(g_exp, g_hyp, g_har, g_time, g_nom);
    T_comb_g(i) = toc;
    tic
    B_comb = B_exp+ B_hyp+ B_har+ B_time+ B_nom;
    T_comb_B(i) = toc;
end
sum(T_eval)
sum(T_comb_f)
sum(T_comb_g)
sum(T_comb_B)
    
    
