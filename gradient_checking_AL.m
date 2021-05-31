% perform gradient and hessian checking on various constraints

clc
clear all

N = 3; %number of wells
T = 5; %number of time steps

%generate a simple setting and randomly set initialization
[x,functionParams,params,l,u] = gen_case_1(N,T);
x = x + rand(length(x),1); %perturb x so that it is not feasible
h = rand(length(x),1)-0.5;
h = h/norm(h);
error = [];

eval = @ALagrangian;
for n = 1:21
    epsilon = 10^(-n+11);
    x_new = x + epsilon*h;
    [ f, g, B] = eval( x, functionParams , params );
    [ f_new, g_new, B_new] = eval( x_new, functionParams , params );
    error(n) = norm(g'*h*epsilon - (f_new-f));
    error(n) = norm(g'*h - (f_new-f)/epsilon);
end
figure
plot(linspace(-10,10,21),flip(log10(error)),'LineWidth',2)
xlabel('\epsilon size (10^n)');
ylabel('log ||\nabla f(x)^Th - (f(x+\epsilon h) - f(x))/\epsilon||');
title('First order approximation error');
ax = gca; 
ax.FontSize = 11;

for n = 1:21
    epsilon = 10^(-n+11);
    x_new = x + epsilon*h;
    [ f, g, B] = eval( x, functionParams , params );
    [ f_new, g_new, B_new] = eval( x_new, functionParams , params );
    error(n) = norm(epsilon*B(:,:,1)*h - (g_new(:,1)-g(:,1)));
    error(n) = norm(B(:,:,1)*h - (g_new(:,1)-g(:,1))/epsilon);
end
figure
plot(linspace(-10,10,21),flip(log10(error)),'LineWidth',2)
xlabel('\epsilon size (10^n)');
ylabel('log ||\nabla^2 f(x)h - (\nabla f(x + \epsilon h) - \nabla f(x))/\epsilon||');
title('Gradient approximation error');
ax = gca; 
ax.FontSize = 11;