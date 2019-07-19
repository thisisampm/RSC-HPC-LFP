hold on

%maze
rectangle('Position', [ 0 0.4692 1 0.0617 ])
axis( [ -0.2 1.2 0.025 1.175 ] )
set(gca, 'TickLength', [0 0])

%bowl
plot(0.10 * cos(0:pi/50:2*pi) + 0.50, 0.10 * sin(0:pi/50:2*pi) + 0.70, 'k');