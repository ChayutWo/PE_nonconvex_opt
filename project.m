function [output] = project(x,l,u)
%project x onto the box constraint [l, u]
output = x;
for i = 1:length(x)
    if x(i) < l(i)
        %less than lower bound
        output(i) = l(i);
    elseif x(i) > u(i)
        %above upper bound
        output(i) = u(i);
    end
end
end