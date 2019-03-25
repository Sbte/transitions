function [trans_prob, time_steps] = transitions_direct(F, B, z0, phi, dt, tmax, N, rho)
% Compute the transition probability directly
    tsteps = tmax / dt;

    t = 0;
    time_steps = 0;
    z = z0 * ones(1,N);

    % Loop until tmax and see if a transition happens
    for i=1:tsteps
        dW = sqrt(dt) * randn(size(z));

        t = t + dt;
        z = z + dt * F(z) + B * dW;

        time_steps = time_steps + size(z, 2);
        converged = phi(z) > 1-rho;
        z = z(:,~converged);
    end

    ntrans = N - size(z,2);
    trans_prob = ntrans / N;
end