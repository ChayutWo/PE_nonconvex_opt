function [xc] = getCauchypoint(x,l,u,G,c)
%Calculate Cauchy point associated with box constraint
%inpute: x - current iterate, l - lower bound, u - upper bound
% G - hessian, c = gradient
% problem f = 1/2xTGx+cTx

%calculate gradient of the function f = 1/2xTGx+cTx: g
g = G*x+c;
[t_bound,t_sorted] = calculate_t_bound(x,l,u,g);
if isempty(t_sorted)
    %there are nothing to search
    xc = x;
    return;
end
t_old = 0;
x_old = project(x-t_old*g,l,u);
for j = 1:length(t_sorted)
    % search in each interval for the local minima
    t_new = t_sorted(j);
    % get new s(j)
    s_j = -g;
    s_j(t_old>=t_bound) = 0; %if we hit the bound already, don't move further
    
    %calculate coefficient
    f_grad = c'*s_j+x_old'*G*s_j;
    f_hess = s_j'*G*s_j;
    
    %check whether a local minima is in this interval or not
    if f_grad > 0
        % function is increasing
        delta = 0;
    elseif f_grad == 0
        % function doesn't change for first order approximation
        if f_hess <=0
            delta = t_new - t_old;
        else
            delta = 0;
        end
    else
        %function is decreasing
        if f_hess <=0
            delta = t_new - t_old;
        else
            delta = min(t_new - t_old, -f_grad/f_hess);
        end
    end
    
    %check whether we need to search for the next interval or not
    if delta < t_new-t_old
        t_opt = t_old + delta;
        xc = project(x-t_opt*g,l,u);
        return;
    end
    %continue to search
    t_old = t_new;
    x_old = project(x-t_old*g,l,u);
end

%move to the end
t_opt = max(t_sorted);
xc = project(x-t_opt*g,l,u);
end

