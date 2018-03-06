function [trans_prob] = transitions_direct(F, B, dt, tmax, N, rho)
% Compute the transition probability directly
    ntrans = 0;
    tsteps = tmax / dt;
    bstr = '';
    for n = 1:N
        fprintf([bstr, 'n=%d'], n);
        bstr = repmat('\b', 1, length(num2str(n))+2);

        t = 0;
        z = -1;
        dW = sqrt(dt) * randn(1,tsteps);

        % Loop until tmax and see if a transition happens
        for i=1:tsteps
            t = t + dt;
            z = z + dt * F(z) + B * dW(i);
            if dist_fun(z) > 1-rho
                ntrans = ntrans + 1;
                break;
            end
        end
    end
    fprintf('\n');

    trans_prob = ntrans / N;
end