function semilogyNumbMarkers(x, y, numbMarkers, markerStyle, color, makerSize, LineWidth)
% similar to semilogy but plot with defined # markers
% Using the idea in http://www.mathworks.com/matlabcentral/answers/2165-too-many-markers

%%
semilogy(x(1), y(1),markerStyle, 'color', color, 'MarkerSize', makerSize,'Linewidth',LineWidth);hold on;
semilogy(x, y, markerStyle(2),'color', color, 'HandleVisibility', 'off', 'MarkerSize', makerSize,'Linewidth',LineWidth);hold on;
semilogy(x(1:numbMarkers:end), y(1:numbMarkers:end), markerStyle(1), 'color', color, 'MarkerSize', makerSize, 'HandleVisibility', 'off');

end
