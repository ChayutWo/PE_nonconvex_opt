function plot_t(x, params)
% plot open duration for each well

% get parameters
N = params.n_well;
T = params.n_period;
deltaT = params.deltaT;

% get t_i^(j) - t_i^(j-1)
t_tab = array2table(reshape(x(2*N*T+1:3*N*T), T, N));

% create well names to be used as a legend
name = {};
for n = 1:N
    name{end+1} = strcat('well-',string(n));
end
t_tab.Properties.VariableNames = cellstr(name);

% create a stacked plot for production time at each period for each well
s = stackedplot(t_tab, '--o');
s.LineWidth = 2;
s.MarkerEdgeColor = 'black';
s.Color = [0.4660 0.6740 0.1880];
% change ylimit
for n = 1:N
    s.AxesProperties(n).YLimits = [0,deltaT*1.2];
end
% add reference horizontal line to indicate maximum production time
ax = findobj(s.NodeChildren, 'Type','Axes');
arrayfun(@(s)yline(s,deltaT,'LineWidth',1, 'color', uint8([17 17 17]),...
    'LineStyle','-.'),ax);
xlabel('timestep')
title('Open duration in each time step')
ax = gca; 
ax.FontSize = 9;
end

