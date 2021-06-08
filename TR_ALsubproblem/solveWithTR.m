function [x, k, error, delta, rho] = solveWithTR(x, evalAL, options, functionParams, params,l,u)
%solve AL subproblem using TR method
%input: x - current iterate, evalAL - function to compute obj, grad, hess
%of AL, options - parameters for TR, functionParams - parameters for AL,
%params - parameters for the problem, l - lower bound, u - upper bound

%output: x - new iterate, k - number of TR iterations, error - KKT
%condition error, delta - records of TR size used

%General
max_TR_iterations = options.maxIterations;
tol = options.tolerance; %omega as a final tolerance for each TR solver

% Trust region parameters
delta(1)=options.delta_init;
delta_max = options.delta_max;                              
eta = options.eta; % accept reject parameter for TR                              

rho(1:max_TR_iterations)=0; %agreement
delta(2:max_TR_iterations)=0; %TR size
error(1:max_TR_iterations)=0; %store error from KKT condition

k=1;
[fprev, ~, ~ ]= evalAL( x, functionParams , params);

while ( k < max_TR_iterations)
    [~, g, B ]= evalAL( x, functionParams , params);
    
    %Check error
    [KKT_error] = computeKKT_AL(x,functionParams,params,l,u);
    error(k)=KKT_error;
    if (KKT_error<=tol)
        return;
    end
        
    % Solve TR subproblem
    % calculate lower bound and upper bound taken into account TR size
    % (l_infinity norm)
    lower_bound = max(l,x-delta(k)*ones(length(x),1));
    upper_bound = min(u, x+delta(k)*ones(length(x),1));
    x_new = getCauchypoint(x,lower_bound,upper_bound,B,g-B*x);
    x_new = CG_subproblem(x_new,lower_bound,upper_bound, B,g-B*x);
    p = x_new-x; % TR step   
    
    % Size the trust region appropriately
    m1= g'*p + 1/2*p'*B*p; %expected reduction
    [f, ~, ~ ]= evalAL(x+p, functionParams , params); %obj at x + p
    rho(k) = -(fprev-f)/m1;
    
    if rho(k) < 1/4
        % poor approximation within TR, reduce TR size
        delta(k+1) = 1/4 * delta(k);
    else
        if (rho(k) > 3/4 && max(p) == delta(k))
            % good quadratic approximation and hit the bound
            delta(k+1) = min(2*delta(k), delta_max);
        else
            % otherwise, keep TR size unchanged
            delta(k+1)=delta(k);
        end
    end
    
    %Accept or reject step
    if (rho(k) > eta)
        x=x+p;
        fprev=f;
    else
        %disp('reject');
    end
    k=k+1;
end
%disp('exit due to max iteration reaches')
[KKT_error] = computeKKT_AL(x,functionParams,params,l,u);
error(k)=KKT_error;

