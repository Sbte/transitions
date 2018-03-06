function y=dist_fun(x)
y = (x - 1) / 2;
z = (x + 1) / 2;
f1 = 0.25;
f2 = 0.75;
y = f1 - f1 * exp(-0.5 * (z / 0.25) .^2) + f2 * exp(-0.5 * (y / 0.25) .^2);
end