function [trans_prob] = transitions_tams(F, B, z0, phi, dt, tmax, N, N3, rho)
% Compute the transition probability with TAMS

    experiments = {};

    % Compute all initial trajectories
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
    l = [];
    w = [1];
    for i=1:N3
        % Determine the trajectories with the smallest maximum
        % value of the reaction coordinate
        min_val = 1;
        min_idx_list = [0];
        for j=1:N
            if experiments{j}.max_dist < min_val
                min_val = experiments{j}.max_dist;
                min_idx_list = [j];
            elseif experiments{j}.max_dist == min_val
                min_idx_list = [min_idx_list, j];
            end
        end

        if min_val > 1-rho
            % All converged
            break;
        end

        % Number of trajectories that will be branched and the
        % weights associated with them
        l = [l, length(min_idx_list)];
        w = [w, w(end) * (1 - l(end) / N)];

        if l(end) == N
            % Extinction
            break;
        end

        % Branch trajectories
        for min_idx = min_idx_list
            idx = min_idx;
            while ~isempty(find(idx == min_idx_list))
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
        end
        its = its + 1;
    end

    % Compute the probability based on the weights
    W = N * w(end);
    for i = 1:its
        W = W + l(i) * w(i);
    end

    converged = 0;
    for i = 1:N
        if experiments{i}.max_dist > 1-rho
            converged = converged + 1;
        end
    end
    trans_prob = converged * w(end) / W;
end