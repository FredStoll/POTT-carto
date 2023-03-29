%- individiual boxplots!

function boxplot_ind(X,yl,wdth,colors)


% info for making boxplot
quartiles   = quantile(X, [0.25 0.75 0.5]);
iqr         = quartiles(2) - quartiles(1);
Xs          = sort(X);
whiskers(1) = min(Xs(Xs > (quartiles(1) - (1.5 * iqr))));
whiskers(2) = max(Xs(Xs < (quartiles(2) + (1.5 * iqr))));
Y           = [quartiles whiskers];

% jitter for raindrops
jit = (rand(size(X)) - 0.5) * wdth;
drops_pos = jit + yl(1) ;

hold on

h{2} = scatter(X, drops_pos);
h{2}.SizeData = 10;
h{2}.MarkerFaceColor = colors(1,:);
h{2}.MarkerEdgeColor = 'none';

box_pos = [Y(1) yl(1)-(wdth * 0.5) Y(2)-Y(1) wdth];
% mean line
h{4} = line([Y(3) Y(3)], [yl(1)-(wdth * 0.5) yl(1) + (wdth * 0.5)], 'col', colors(2,:), 'LineWidth', 2);

% whiskers
h{5} = line([Y(2) Y(5)], [yl(1) yl(1)], 'col', colors(2,:), 'LineWidth', 1);
h{6} = line([Y(1) Y(4)], [yl(1) yl(1)], 'col', colors(2,:), 'LineWidth', 1);

h{3} = rectangle('Position', box_pos);
set(h{3}, 'EdgeColor', colors(2,:))
set(h{3}, 'LineWidth', 1.5);

