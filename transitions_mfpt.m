function [trans_prob, mfpt] = transitions_mfpt(F, B, dt, tmax, N, rho)
% Compute mean first passage time and from that the transition probability
    mfpt = 0;
    M = 1000;

    bstr = '';
    for n = 1:N
        fprintf([bstr, 'n=%d'], n);
        bstr = repmat('\b', 1, length(num2str(n))+2);

        t = 0;
        z = -1;

        % Loop with 1000 random numbers at a time for performance
        converged = false;
        while ~converged
            dW = sqrt(dt) * randn(1,M);
            for i=1:M
                t = t + dt;
                z = z + dt * F(z) + B * dW(i);
                if dist_fun(z) > 1-rho
                    mfpt = mfpt + t;
                    converged = true;
                    break;
                end
            end
        end
    end
    fprintf('\n');

    mfpt = mfpt / N;
    trans_prob = 1 - exp(-1 / mfpt * tmax);
end