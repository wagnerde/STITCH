function noLegend(h)
% Excludes a given curve on a plot from having a legend entry
% h = a handle to the plot.
% Example:   h = plot(x,y); noLegend(h);

set(get(get(h,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off'); % Exclude line from legend

end