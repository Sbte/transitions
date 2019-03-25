function [trans_prob, time_steps, mfpt] = transitions_mfpt(F, B, z0, phi, dt, tmax, N, rho)
% Compute mean first passage time and from that the transition probability
    mfpt = 0;

    t = 0;
    time_steps = 0;
    z = z0 * ones(1,N);

    % Loop until all samples have transitioned
    while ~isempty(z)
        dW = sqrt(dt) * randn(size(z));

        t = t + dt;
        z = z + dt * F(z) + B * dW;
        time_steps = time_steps + size(z,2);

        converged = phi(z) > 1-rho;
        z = z(:, ~converged);

        mfpt = mfpt + t * sum(converged);
    end

    mfpt = mfpt / N;
    trans_prob = 1 - exp(-1 / mfpt * tmax);
end