function [trans_prob] = transitions_direct(F, B, dt, tmax, N, rho)
% Compute the transition probability directly
    ntrans = 0;
    tsteps = tmax / dt;

    t = 0;
    z = ones(N,1) * -1;
    dW = sqrt(dt) * randn(N,tsteps);

    % Loop until tmax and see if a transition happens
    bstr = '';
    for i=1:tsteps
        fprintf([bstr, 't=%10f'], t);
        bstr = repmat('\b', 1, 12);

        t = t + dt;
        z = z + dt * F(z) + B * dW(:,i);
    end

    ntrans = sum(dist_fun(z) > 1-rho);
    trans_prob = ntrans / N;
end