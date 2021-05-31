function plot_q(x, params)
%PLOT_Q Summary of this function goes here
%   Detailed explanation goes here
%get parameters
N = params.n_well;
T = params.n_period;

q_tab = array2table(reshape(x(1:N*T), T, N), 'VariableNames',{'well_1','well_2','well_3'});
q_tab = reshape(x(1:N*T), T, N);
name = {};
for n = 1:N
    name{end+1} = strcat('well ',string(n));
end

area(q_tab)
grid on
colormap summer
set(gca,'Layer','top')
ylim([0,35])
xlim([1,T])
xlabel('timestep')
ylabel('production volume')
title('Production over time')
legend(name)
end

