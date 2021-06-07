% perform gradient and hessian checking on various constraints

clc
clear all
N_exponential = 1;
N_hyperbolic = 1;
N_harmonic = 1;
N = N_exponential + N_hyperbolic + N_harmonic; %number of wells
T = 5; %number of time steps

%generate a simple setting and randomly set initialization
[x,functionParams,params,l,u] = gen_case_1(N_exponential, N_hyperbolic, N_harmonic, T);
x = x + rand(length(x),1); %perturb x so that it is not feasible

%generate a unit vector of direction h
h = rand(length(x),1)-0.5;
h = h/norm(h);


%eval = @computeTimeConstr;
%eval = @computeNomConstr;
%eval = @computeExponentialConstr;
eval = @computeHyperbolicConstr;
%eval = @computeHarmonicConstr;
%eval = @combineConst;
%eval = @computeObjGradHess;

%perform numerical gradient checking
error = [];
for n = 1:21
    epsilon = 10^(-n+11); % step size from 1e-10 to 1e10
    x_new = x + epsilon*h; % perturb x from starting value: x_new
    [ f, g, B] = eval( x, params );
    size(g)
    size(B)
    [ f_new, g_new, B_new] = eval( x_new, params );
    %error(n) = norm(g'*h*epsilon - (f_new-f));
    error(n) = norm(g'*h - (f_new-f)/epsilon);
end
figure
plot(linspace(-10,10,21),flip(log10(error)),'LineWidth',2)
xlabel('\epsilon size (10^n)');
ylabel('log ||\nabla f(x)^Th - (f(x+\epsilon h) - f(x))/\epsilon||');
title('First order approximation error');
ax = gca; 
ax.FontSize = 11;

%perform numerical Hessian checking
error = [];
for n =  1:21
    epsilon = 10^(-n+11);
    x_new = x + epsilon*h;
    [ f, g, B] = eval( x, params );
    [ f_new, g_new, B_new] = eval( x_new, params );
    %error(n) = norm(epsilon*B(:,:,1)*h - (g_new(:,1)-g(:,1)));
    error(n) = norm(B(:,:,1)*h - (g_new(:,1)-g(:,1))/epsilon);
end
figure
plot(linspace(-10,10,21),flip(log10(error)),'LineWidth',2)
xlabel('\epsilon size (10^n)');
ylabel('log ||\nabla^2 f(x)h - (\nabla f(x + \epsilon h) - \nabla f(x))/\epsilon||');
title('Gradient approximation error');
ax = gca; 
ax.FontSize = 11;