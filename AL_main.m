clc
clear all
rng(5) %set random seed
option = 1; %the case that we want to test

%generate a simple setting and randomly set initialization
if option == 1
    %case1: Feasible start
    N = 3; %number of wells
    T = 5; %number of time steps
    [x,functionParams,params,l,u] = gen_case_1(N,T);
elseif option == 2
    %case2: Infeasible start
    N = 3; %number of wells
    T = 5; %number of time steps
    [x,functionParams,params,l,u] = gen_case_2(N,T);
elseif option == 3
    %case3: long production and check decline pattern
    N = 3; %number of wells
    T = 10; %number of time steps
    [x,functionParams,params,l,u] = gen_case_3(N,T);
elseif option == 4
    %case4: Tradeoff between production and decline
    N = 3; %number of wells
    T = 3; %number of time steps
    %T = 15;
    [x,functionParams,params,l,u] = gen_case_4(N,T);
    %params.CGR(7:9) = 0.1;
end

%case2: Infeasible start
%[x,functionParams,params,l,u] = gen_case_2(N,T);
x_old = x;
%set AL parameters: lambda and mu
lambda=functionParams.lambda;
lambda_prev=lambda;
lambda_error=0;

mu=functionParams.penalty;
mu_max = 1e5;
omega = 1/mu;
eta = 1/mu^(0.1);

%set optimization parameters
maxIter = 2000; %max iteration for TR to solve each subproblem
kkt_tol=1e-3;
feas_tol = 1e-4;
iter=1;

%need to change this function
tr_options = struct('maxIterations', 2000,...,
    'tolerance', omega,...,
    'delta_max', 100,...,
    'delta_init', 10,...,
    'eta',1e-3); 

functionParams =  struct('penalty',mu,...,
                         'lambda',lambda);


fprintf('Iteration \t |c|\t Lambda_update \t x_update \t   norm_gAL \t    mu\n');
lambda_track = [];
KKT_error_0 = computeKKT_AL(x,functionParams,params,l,u);

while(iter < maxIter)
    tr_options.tolerance = omega;
    functionParams.lambda=lambda;
    functionParams.penalty=mu;
    
    %need to change this function
    [x, k, error, delta, rho] = solveWithTR(x, @ALagrangian, tr_options, functionParams, params,l,u);
    change_x = norm(x-x_old)/norm(x_old);
    x_old = x;
    
    %evaluate current infeasibility
    [c, g_comb, B_comb] = combineConst( x, params );
    %evaluate KKT error
    KKT_error = computeKKT_AL(x,functionParams,params,l,u);

    fprintf('   %d  \t  %10.3e \t %10.3e \t %10.3e \t %10.3e \t %10.3e\n', iter,norm(c),...
        lambda_error(end), change_x, KKT_error, mu);
     
    if (norm(c) < eta)
        if norm(c)<= feas_tol && KKT_error/KKT_error_0 < kkt_tol
            break;
        end
        
        lambda = lambda - mu*c;
        lambda_error(iter+1) = norm(lambda-lambda_prev)/norm(lambda_prev);
        lambda_prev = lambda;
        eta = eta/mu^(0.9);
        omega = omega/mu;
    else
        lambda_error(iter+1) = 0;
        mu = 100*mu;
        eta = 1/mu^(0.1);
        omega = 1/mu;
    end
    if mu > mu_max
        mu = mu_max;
    end
    lambda_track(iter) = lambda(1);
    iter =iter+1;
end

fprintf('\n Solution is:\n');
disp(x);
disp(lambda);
q_tab = array2table(reshape(x(1:N*T), T, N), 'VariableNames',{'well_1','well_2','well_3'});
figure;
plot_t(x, params);
figure;
plot_q(x, params);
