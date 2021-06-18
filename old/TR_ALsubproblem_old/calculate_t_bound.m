function [t_bound,t_sorted] = calculate_t_bound(x,l,u,g)
%Calculate t and sorted t for projected gradient path given x,u,l,g
%Input: x - current iterate, u - upper bound, l - lower bound, g - gradient
%output: t_bound - t for each boundary, t_sorted - sorted version of
%t_bound
t_bound = inf(length(x),1);
for i = 1:length(x)
    if g(i) < 0 && u(i) < inf
        t_bound(i) = (x(i)-u(i))/g(i);
    elseif g(i) > 0 && l(i) > -inf
        t_bound(i) = (x(i)-l(i))/g(i);
    end
end
t_sorted = unique(sort(t_bound));

if t_sorted(1) == 0
    % no need to keep 0 as it hits the boundary right away
    t_sorted = t_sorted(2:end);
end
end

