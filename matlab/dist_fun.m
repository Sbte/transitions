function y=dist_fun(x, a, b)
    y = ones(1,size(x,1)) * (x-b).^2;
    z = ones(1,size(x,1)) * (x-a).^2;
    f1 = 0.5;
    f2 = 0.5;
    y = f1 - f1 * exp(-8 * z) + f2 * exp(-8 * y);
end