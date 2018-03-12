function [trans_prob, mfpt] = transitions_mfpt(F, B, dt, tmax, N, rho)
% Compute mean first passage time and from that the transition probability
    mfpt = 0;

    t = 0;
    z = ones(N,1) * -1;

    % Loop until all samples have transitioned
    while ~isempty(z)
        dW = sqrt(dt) * randn(length(z),1);

        t = t + dt;
        z = z + dt * F(z) + B * dW;

        converged = dist_fun(z) > 1-rho;
        z = z(find(~converged));

        mfpt = mfpt + t * sum(converged);
    end

    mfpt = mfpt / N;
    trans_prob = 1 - exp(-1 / mfpt * tmax);
end