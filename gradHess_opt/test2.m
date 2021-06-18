clc
clear all
%rng(2) %set random seed
option = 1; %the case that we want to test

%generate a simple setting and randomly set initialization
if option == 1
    %case1: Infeasible start in a small field to check the optimization
    N_exponential = 2;
    N_hyperbolic = 2;
    N_harmonic = 2;
    N = N_exponential + N_hyperbolic + N_harmonic; %number of wells
    T = 5; %number of time steps
    [x,functionParams,params,l,u] = gen_case_1(N_exponential, N_hyperbolic, N_harmonic, T);
end

err = [];
count = 0;
for trial = 1:1000
    [x,functionParams,params,l,u] = gen_case_1(N_exponential, N_hyperbolic, N_harmonic, T);
    x = 1000*x;
    [ f, g, B] = ALagrangian( x, functionParams , params);
    [ f_AL, grad_AL, hess_AL ] = ALagrangian_fgB( x, functionParams , params);
    if norm(B - hess_AL)>0
        count = count +1;
        err(trial) =  norm(B - hess_AL);
    end
end

