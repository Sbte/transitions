function [trans_prob] = transitions_direct(F, B, dt, tmax, N, rho)
% Compute the transition probability directly
    tsteps = tmax / dt;

    t = 0;
    z = ones(N,1) * -1;

    % Loop until tmax and see if a transition happens
    for i=1:tsteps
        dW = sqrt(dt) * randn(N,1);

        t = t + dt;
        z = z + dt * F(z) + B * dW;
    end

    ntrans = sum(dist_fun(z) > 1-rho);
    trans_prob = ntrans / N;
end