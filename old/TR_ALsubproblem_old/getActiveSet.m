function [Ax] = getActiveSet(x,l,u)
%Get active set at current x(t)
Ax = [];
count = 1;
for i = 1:length(x)
    if x(i) == l(i) || x(i) == u(i)
        Ax(count) = i;
        count = count +1;
    end
end
end

