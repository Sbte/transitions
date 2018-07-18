function [trans_prob] = transitions_tams(F, B, z0, phi, dt, tmax, N, N3, rho)
% Compute the transition probability with TAMS

    experiments = {};

    for i=1:N
        t = 0;
        z = z0;

        M = ceil(1 / dt * (tmax - t));
        dW = randn(size(z,1),M) * sqrt(dt);

        clear experiment;
        experiment.max_dist = 0;
        experiment.x = [z];
        experiment.t = [t];
        experiment.d = [0];
        for j=1:M
            t = t + dt;
            z = z + dt * F(z) + B * dW(:,j);
            dist = phi(z);
            if dist > experiment.max_dist
                experiment.x = [experiment.x, z];
                experiment.t = [experiment.t, t];
                experiment.d = [experiment.d, dist];
                experiment.max_dist = dist;
            end
        end
        experiments{i} = experiment;
    end

    its = 0;
    for i=1:N3
        min_val = 1;
        min_idx = 0;
        for j=1:N
            if experiments{j}.max_dist < min_val
                min_val = experiments{j}.max_dist;
                min_idx = j;
            end
        end

        if min_val > 1-rho
            break;
        end

        idx = min_idx;
        while idx == min_idx
            idx = randi(N,1);
        end

        same_dist_idx = 1;
        while experiments{idx}.d(same_dist_idx) < min_val
            same_dist_idx = same_dist_idx + 1;
        end

        experiments{min_idx}.x = experiments{idx}.x(:,same_dist_idx);
        experiments{min_idx}.t = experiments{idx}.t(same_dist_idx);
        experiments{min_idx}.d = experiments{idx}.d(same_dist_idx);
        experiments{min_idx}.max_dist = experiments{idx}.d(same_dist_idx);
        t = experiments{min_idx}.t(end);
        z = experiments{min_idx}.x(:,end);
        M = ceil(1/dt*(tmax-t));
        dW = randn(size(z,1),M) * sqrt(dt);
        for j=1:M
            t = t + dt;
            z = z + dt * F(z) + B * dW(:,j);
            dist = phi(z);
            if dist > experiments{min_idx}.max_dist
                experiments{min_idx}.x = [experiments{min_idx}.x, z];
                experiments{min_idx}.t = [experiments{min_idx}.t, t];
                experiments{min_idx}.d = [experiments{min_idx}.d, dist];
                experiments{min_idx}.max_dist = dist;
                if dist > 1-rho
                    break;
                end
            end
        end
        its = its + 1;
    end

    converged = 0;
    for i = 1:N
        if experiments{i}.max_dist > 1-rho
            converged = converged + 1;
        end
    end

    W = N * (1.0 - 1.0 / N) ^ its;
    for i=1:(its-1)
        W = W + (1.0 - 1.0 / N) ^ i;
    end
    trans_prob = (converged * (1.0 - 1.0 / N) ^ its) / W;
end