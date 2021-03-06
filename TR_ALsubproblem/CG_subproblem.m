function [output] = CG_subproblem(x,l,u, G, c)
%Use CG to solve subproblem
%min q(x) = 1/2xTGx+xTc st Ax = b

%get current active set
Ax =  getActiveSet(x,l,u);
if length(Ax) == length(x)
    %all constraints are active. Don't have to solve this step
    output = x;
    return;
end   
%Construct matrix A for equality constraint
if isempty(Ax)
    %unconstrain optimization. Solve with normal CG
    P = eye(length(x), length(x));
else
    %equality constrained optimization
    %construction projection matrix P
    A = ones(length(x),1);
    A(Ax) = 0;
    P = diag(A);
end
% x is already feasible so we can start CG iteration right away
r = G*x+c; % residual
% g = P*r
g = r;
g(Ax) = 0;
d = -g;
tol = 1e-6;
max_iter = 250;
count = 1;
output = x;
while abs(r'*g) > tol && count < max_iter
    Gd = G*d;
    curvature = d'*Gd;
    rTg = r'*g;
    alpha = rTg/curvature; %step size in CG iteration
    [~,t_sorted] = calculate_t_bound(x,l,u,-d);
    t_min = t_sorted(1);
    if (curvature <=0) || (alpha > t_min)
        %nonpositive curvature is found
        %move along d until hit some boundary
        %if the update would hit the boundary, stop at the boundary
        output = project(x+t_min*d,l,u);
        return
    end
    % update parameter for that iteration
    x = x + alpha*d; % move x to the new location
    r_new = r+alpha*Gd; % new residual
    % g_new = P*r_new
    g_new = r_new;
    g_new(Ax) = 0;
    beta = r_new'*g_new/rTg; % new step size to update CG direction
    d = -g_new+beta*d; % update CG direction
    g = g_new;
    r = r_new;
    count = count + 1;
end
end

