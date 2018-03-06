function [trans_prob] = transitions_gpa(F, B, dt, tmax, N, rho)
% Compute the transition probability with GPA

    N = 5000;
    tstep = dt * 100;
    z = ones(1, N) * -1;
    converged = zeros(1, N);
    eta = 1;
    beta = 1;
    Y = ones(1, N);

    tmax = 100;

    bstr = '';
    for ti=1:tmax/tstep
        fprintf([bstr, 't=%10f'], ti * tstep);
        bstr = repmat('\b', 1, 12);

        weights = exp(beta * V(z));

        eta = eta * 1 / N * sum(weights);

        samples = datasample(1:N, N, 'Weights', weights);

        z = z(samples);
        Y = Y(samples);
        converged = converged(samples);

        Y = Y .* exp(-beta * V(z));

        for i = 1:N
            dW = randn(1,tstep/dt) * sqrt(dt);
            for j=1:tstep / dt
                z(i) = z(i) + dt * F(z(i)) + B * dW(j);
                if dist_fun(z(i)) > 0.95
                    converged(i) = true;
                end
            end
        end
    end
    fprintf('\n');

    trans_prob = 0;

    for i=1:N
        if converged(i)
            trans_prob = trans_prob + Y(i);
        end
    end

    trans_prob = trans_prob * eta;
    trans_prob = trans_prob / N;
end

function y=V(x)
    y = (x + 1) / 2;
end