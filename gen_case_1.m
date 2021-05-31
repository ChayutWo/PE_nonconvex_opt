function [x,functionParams,params,l,u] = gen_case_1(N,T)
%Generate data to be used for optimization
%Input: N - number of wells, T - number of time steps
%Output: x - vector of variables to be optimized
%        functionParams - mu and lambda for AL method
%        params - array of values for constants
%        l - lower bound of x
%        u - upper bound of x

% Case_1: feasible start case. One well with low potential but high CGR
%1. Initial Potential
q_pot_i = ones(N,1)*20;
q_pot = repelem(q_pot_i,T);

%2. Nomination rate
q_nom = ones(T,1)*30;

%3. deltaT stepsize
deltaT = 1;

%4. decline rate of each well
decline_i = ones(N,1)*0.2;
decline = repelem(decline_i,T);

%5. CGR
CGR_i = linspace(20,2,N);
CGR = repelem(CGR_i,T)';

%6. Price
depreciation_rate = 0.9;
price_j = 10*depreciation_rate.^((1:T)-1)';
price = repmat(price_j, N,1);

% Constants in optimization problem
params = struct('q_pot', q_pot, 'q_nom', q_nom, 'deltaT', ...
    deltaT, 'decline', decline, 'CGR', CGR, 'price', price, 'n_well', N,...
    'n_period', T);

% Initial values for optimization
t_mat = rand(T, N)/3;
t = reshape(cumsum(t_mat,1), T*N,1); % open/close time
t_prev = t(1:N*T-1);
t_prev = vertcat(0,t_prev);
for n = 1:N
    t_prev(1+(n-1)*T) = 0;
end
% average production rate per time step to satisfy harmonic decline
q_g = (1/deltaT)*(q_pot./decline).*log((1+decline.*t)./(1+decline.*t_prev));
b = reshape(t_mat,T*N,1); %auxiliary variable for time constraints
s = q_nom - sum(reshape(q_g,T,N),2); %auxiliary variable for nomination constraints
x = vertcat(q_g, t,b,s);

% initial value for mu and lambda
mu=10;
lambda = rand(2*N*T+T,1);
functionParams =  struct('penalty',mu,...,
                         'lambda',lambda);
                     
l = ones(length(x),1)*(-inf);
l((2*N*T+1):end) = 0;
l = zeros(length(x),1);
u = ones(length(x),1)*(inf);
u((2*N*T+1):(3*N*T)) = deltaT; %max bound for b
u(1:N*T) = q_pot; %max bound of q_g
u(N*T+1:2*N*T) = repmat(deltaT*(1:T)',N,1);
end

