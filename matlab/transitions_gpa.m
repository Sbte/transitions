function [trans_prob] = transitions_gpa(F, B, z0, phi, dt, tmax, N, rho)
% Compute the transition probability with GPA

    tstep = dt * 100;
    z = z0 * ones(1,N);
    converged = zeros(1, N);
    eta = 1;
    beta = 1;
    Y = ones(1, N);

    for ti=1:tmax/tstep
        weights = exp(beta * V(z));

        eta = eta * 1 / N * sum(weights);

        samples = datasample(1:N, N, 'Weights', weights);

        z = z(:,samples);
        Y = Y(samples);
        converged = converged(samples);

        Y = Y .* exp(-beta * V(z));

        for j=1:tstep / dt
            dW = randn(size(z)) * sqrt(dt);
            z = z + dt * F(z) + B * dW;
            converged = converged | (phi(z) > 1-rho);
        end
    end

    trans_prob = 0;

    for i=1:N
        if converged(i)
            trans_prob = trans_prob + Y(i);
        end
    end

    trans_prob = trans_prob * eta;
    trans_prob = trans_prob / N;

    function y=V(x)
        y = phi(x);
    end
end