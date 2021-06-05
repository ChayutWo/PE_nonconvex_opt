function [x,functionParams,params,l,u] = gen_case_1(N_exponential, N_hyperbolic, N_harmonic, T)
%Generate data to be used for optimization
%Input: 
% N_exponential, N_hyperbolic, N_harmonic - number of wells, 
% T - number of time steps
%Output: x - vector of variables to be optimized
%        functionParams - mu and lambda for AL method
%        params - array of values for constants
%        l - lower bound of x
%        u - upper bound of x

% index
N = N_exponential + N_hyperbolic + N_harmonic; % total number of wells
expo_index = (1:N_exponential*T);
hyper_index = ((N_exponential*T+1):(N_exponential*T+N_hyperbolic*T));
harm_index = ((N_exponential*T+N_hyperbolic*T+1):(N*T));


%1. Initial Potential: size N*T
q_pot_i = ones(N,1)*20;
q_pot = repelem(q_pot_i,T);

%2. Nomination rate: size T
q_nom = ones(T,1)*30;

%3. deltaT stepsize
deltaT = 1;

%4. decline rate of each well: size N*T
decline_i = ones(N,1)*0.2;
decline = repelem(decline_i,T);
hyperbolic_const = 0.5; % same b parameter for hyperbolic decline

%5. CGR: size N*T
CGR_i = linspace(20,2,N);
CGR = repelem(CGR_i,T)';

%6. Price: size N*T
depreciation_rate = 0.9;
price_j = 10*depreciation_rate.^((1:T)-1)';
price = repmat(price_j, N,1);

%7. other parameters
q_abandonment = 1; % abandonment rate

% Calculate maximum time until reach abandonment rate
t_max_expo = log(q_pot./q_abandonment)./decline;
t_max_hyper = (((q_pot./q_abandonment).^hyperbolic_const)-1)./(hyperbolic_const*decline);
t_max_harm = ((q_pot./q_abandonment)-1)./decline;

% Calculate Gp or reserves for each well
Gp_expo = (q_pot./decline).*(exp(-decline.*0) - exp(-decline.*t_max_expo));
Gp_hyper = (q_pot./((hyperbolic_const-1)*decline)).*((1+hyperbolic_const*decline.*t_max_hyper).^(1-1/hyperbolic_const)-(1+hyperbolic_const*decline.*0).^(1-1/hyperbolic_const));
Gp_harm = (q_pot./decline).*log((1+decline.*t_max_harm)./(1+decline.*0));
Gp = [Gp_expo(expo_index); Gp_hyper(hyper_index); Gp_harm(harm_index)];

% Constants in optimization problem
params = struct('q_pot', q_pot, 'q_nom', q_nom, 'deltaT', ...
    deltaT, 'decline', decline, 'CGR', CGR, 'price', price, 'n_well', N,...
    'n_expo', N_exponential, 'n_hyper', N_hyperbolic, 'n_harm', N_harmonic,... 
    'expo_index', expo_index, 'hyper_index', hyper_index, 'harm_index', harm_index,...
    'hyperbolic_const', hyperbolic_const ,'n_period', T, 'Gp', Gp);

% Initial values for optimization
t_mat = rand(T, N)/3;
t = reshape(cumsum(t_mat,1), T*N,1); % open/close time
t_prev = t(1:N*T-1);
t_prev = vertcat(0,t_prev);
for n = 1:N
    t_prev(1+(n-1)*T) = 0;
end
% average production rate per time step to satisfy different decline
q_expo = (1/deltaT)*(q_pot./decline).*(exp(-decline.*t_prev) - exp(-decline.*t));
q_hyper = (1/deltaT)*(q_pot./((hyperbolic_const-1)*decline)).*((1+hyperbolic_const*decline.*t).^(1-1/hyperbolic_const)-(1+hyperbolic_const*decline.*t_prev).^(1-1/hyperbolic_const));
q_harm = (1/deltaT)*(q_pot./decline).*log((1+decline.*t)./(1+decline.*t_prev));
q_g = q_expo;
q_g(hyper_index) =  q_hyper(hyper_index);
q_g(harm_index) = q_harm(harm_index);

b = reshape(t_mat,T*N,1); %auxiliary variable for time constraints
s = q_nom - sum(reshape(q_g,T,N),2); %auxiliary variable for nomination constraints
x = vertcat(q_g, t,b,s);
x = rand(length(x),1); % infeasible start
% initial value for mu and lambda
mu=10;
lambda = rand(2*N*T+T,1);
functionParams =  struct('penalty',mu,...,
                         'lambda',lambda);
                     
l = zeros(length(x),1);
u = ones(length(x),1)*(inf);
u(N*T+expo_index) = t_max_expo(expo_index); %max t for wells with exponential decline
u(N*T+hyper_index) = t_max_hyper(hyper_index); %max t for wells with hyperbolic decline
u(N*T+harm_index) = t_max_harm(harm_index); %max t for wells with harmonic decline
u((2*N*T+1):(3*N*T)) = deltaT; %max bound for b
u(1:N*T) = q_pot; %max bound of q_g
end

