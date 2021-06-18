clc
clear all
rng(0) %set random seed
option = 1; %the case that we want to test

%generate a simple setting and randomly set initialization
if option == 1
    %case1: Infeasible start in a small field to check the optimization
    N_exponential = 1;
    N_hyperbolic = 1;
    N_harmonic = 1;
    N = N_exponential + N_hyperbolic + N_harmonic; %number of wells
    T = 5; %number of time steps
    [x,functionParams,params,l,u] = gen_case_1(N_exponential, N_hyperbolic, N_harmonic, T);
end

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
maxIter = 1000; %max iteration for TR to solve each subproblem
kkt_tol=1e-3;
feas_tol = 1e-3;
iter=1;

%save optimization parameters into a variable: tr_options
tr_options = struct('maxIterations', 1000,...,
    'tolerance', omega,...,
    'delta_max', 100,...,
    'delta_init', 10,...,
    'eta',1e-3); 

functionParams =  struct('penalty',mu,...,
                         'lambda',lambda);


fprintf('Iteration \t |c| \t Lambda_update \t x_update \t norm_gAL \t mu\n');
lambda_track = [];
KKT_error_0 = computeKKT_AL(x,functionParams,params,l,u);

time_TR = [];
time_update_1 = [];
time_update_2 = [];
while(iter < maxIter)
    tr_options.tolerance = omega;
    tr_options.eta = eta;
    functionParams.lambda=lambda;
    functionParams.penalty=mu;
    tic
    %solve subproblem using Trust-Region method with gradient projection
    [x, k, error, delta, rho] = solveWithTR(x, @ALagrangian_f, @ALagrangian_fgB, tr_options, functionParams, params,l,u);
    time_TR(iter) = toc;
    tic
    change_x = norm(x-x_old)/norm(x_old); % keep track how x change over iterations
    x_old = x;
    
    %evaluate current infeasibility
    [c] = combineConst_f( x, functionParams , params );
    %evaluate KKT error
    KKT_error = computeKKT_AL(x,functionParams,params,l,u);

    fprintf('   %d  \t  %10.3e \t %10.3e \t %10.3e \t %10.3e \t %10.3e\n', iter, norm(c),...
        lambda_error(end), change_x, KKT_error, mu);
    time_update_1(iter) = toc;
    tic
    if (norm(c) < eta)
        %we have good feasibility
        if norm(c)<= feas_tol && KKT_error/KKT_error_0 < kkt_tol
            %we have good feasibility and small KKT error, finish
            break;
        end
        %update optimization parameters
        lambda = lambda - mu*c;
        %track how multipliers change over iterations
        lambda_error(iter+1) = norm(lambda-lambda_prev)/norm(lambda_prev);
        lambda_prev = lambda;
        eta = eta/mu^(0.9);
        omega = omega/mu;
    else
        %we have bad feasibility, not update lambda and ease up criteria
        lambda_error(iter+1) = 0;
        mu = 100*mu;
        %mu = 1.2*mu;
        eta = 1/mu^(0.1);
        omega = 1/mu;
    end
    if mu > mu_max
        mu = mu_max;
    end
    lambda_track(iter) = lambda(1);
    iter =iter+1;
    time_update_2(iter) = toc;
end

%plotting results using plot_t and plot_q functions
fprintf('\n Solution is:\n');
%disp(x);
%disp(lambda);
figure;
plot_t(x, params);
figure;
plot_q(x, params);
