function plot_q(x, params)
% Perform stacked area plot of production rate from each well over time

% get parameters
N = params.n_well;
T = params.n_period;
deltaT = params.deltaT;
nomination = max(params.q_nom);

% get production volumne over time for each well from x
q_tab = reshape(x(1:N*T), T, N);
q_tab = q_tab*deltaT; %take into account the length of each time step

% create well name to be used as a plot legend
name = {};
for n = 1:N
    name{end+1} = strcat('well-',string(n));
end
%perform area plot
area(q_tab)
grid on
colormap summer
set(gca,'Layer','top')
ylim([0,nomination+10])
xlim([1,T])
xlabel('timestep')
ylabel('production volume')
title('Production over time')
legend(name)
end

