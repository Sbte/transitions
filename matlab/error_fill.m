function error_fill(x, y, below, above, col)

d1 = y - below;
d2 = fliplr(y + above);
patch([x, fliplr(x)], [d1, d2], 'red', 'FaceColor', col, 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off')