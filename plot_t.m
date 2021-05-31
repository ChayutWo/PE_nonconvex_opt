function plot_t(x, params)
%PLOT_T Summary of this function goes here
%   %plot open duration for each well

%get parameters
N = params.n_well;
T = params.n_period;


t_tab = array2table(reshape(x(2*N*T+1:3*N*T), T, N));
name = {};
for n = 1:N
    name{end+1} = strcat('well_',string(n));
end
t_tab.Properties.VariableNames = cellstr(name);
s = stackedplot(t_tab, '--o');
s.LineWidth = 2;
s.MarkerEdgeColor = 'black';
s.Color = [0.4660 0.6740 0.1880];
s.AxesProperties(1).YLimits = [0,1.1];
s.AxesProperties(2).YLimits = [0,1.1];
s.AxesProperties(3).YLimits = [0,1.1];
xlabel('timestep')
title('Open duration in each time step')
ax = gca; 
ax.FontSize = 10;
end

