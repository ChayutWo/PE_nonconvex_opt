function [output] = project(x,l,u)
%project x onto the box constraint [l, u]
output = max(x,l);
output = min(output, u);
end