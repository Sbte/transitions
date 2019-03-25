function [trans_prob, time_steps, mfpt] = transitions_ams(F, B, z0, phi, dt, tmax, N, rho)
% Compute the mean first passage time and transition probability with AMS

    N2 = N * 10;
    N3 = 1000000;
    M = 1000;

    experiments = {};

    % Compute trajectories until they reach C
    for i=1:N2
        t = 0;
        z = z0;
        steps = 0;
        converged = false;
        while ~converged
            dW = randn(size(z,1),M) * sqrt(dt);
            for j=1:M
                t = t + dt;
                z = z + dt * F(z) + B * dW(:,j);
                dist = phi(z);
                if dist > 0.1
                    clear experiment;
                    experiment.start_time = t;
                    experiment.x = [z];
                    experiment.t = [0];
                    experiment.d = [dist];
                    experiment.return_time = 0;
                    experiment.max_dist = dist;
                    experiment.steps = steps + j;
                    experiments{i} = experiment;
                    converged = true;
                    break;
                end
            end
            steps = steps + M;
        end
    end

    % Compute trajectories until they reach A or B from C
    for i=1:N2
        t = 0;
        z = experiments{i}.x(:,end);
        steps = 0;
        converged = false;
        while ~converged
            dW = randn(size(z,1),M) * sqrt(dt);
            for j=1:M
                t = t + dt;
                z = z + dt * F(z) + B * dW(:,j);
                dist = phi(z);
                if dist > experiments{i}.max_dist
                    experiments{i}.x = [experiments{i}.x, z];
                    experiments{i}.t = [experiments{i}.t, t];
                    experiments{i}.d = [experiments{i}.d, dist];
                    experiments{i}.steps = experiments{i}.steps + steps + j;
                    experiments{i}.max_dist = dist;
                    if dist > 1-rho
                        converged = true;
                        break;
                    end
                elseif dist < rho
                    if experiments{i}.return_time == 0
                        experiments{i}.return_time = t;
                    end
                    experiments{i}.steps = experiments{i}.steps + steps + j;
                    converged = true;
                    break;
                end
            end
            steps = steps + M;
        end
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
            while idx == min_idx
                idx = randi(N,1);
            end

            same_dist_idx = 1;
            while experiments{idx}.d(same_dist_idx) < min_val
                same_dist_idx = same_dist_idx + 1;
            end

            experiments{min_idx}.x = experiments{idx}.x(:,1:same_dist_idx);
            experiments{min_idx}.t = experiments{idx}.t(1:same_dist_idx);
            experiments{min_idx}.d = experiments{idx}.d(1:same_dist_idx);
            experiments{min_idx}.max_dist = experiments{idx}.d(same_dist_idx);
            t = experiments{min_idx}.t(end);
            z = experiments{min_idx}.x(:,end);
            steps = 0;
            converged = false;
            while ~converged
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
                            experiments{min_idx}.steps = experiments{min_idx}.steps + steps + j;
                            converged = true;
                            break;
                        end
                    elseif dist < rho
                        experiments{min_idx}.steps = experiments{min_idx}.steps + steps + j;
                        converged = true;
                        break;
                    end
                end
                steps = steps + M;
            end
        end
        its = its + 1;
    end

    total_tr = 0;
    total_t1 = 0;
    total_t2 = 0;
    num_t1 = N2;
    num_t2 = 0;

    for i=1:N2
        total_t1 = total_t1 + experiments{i}.start_time;
        total_t2 = total_t2 + experiments{i}.return_time;
        if (experiments{i}.return_time > 0)
            num_t2 = num_t2 + 1;
        end
    end

    % Compute the probability of reaching B before A
    W = N * w(end);
    for i = 1:its
        W = W + l(i) * w(i);
    end

    converged = 0;
    time_steps = 0;
    for i = 1:N
        if experiments{i}.max_dist > 1-rho
            converged = converged + 1;
            total_tr = total_tr + experiments{i}.t(end);
        end
        time_steps = time_steps + experiments{i}.steps;
    end
    alpha = converged * w(end) / W;

    % Compute the mean first passage time
    meann = 1.0 / alpha - 1.0;
    mfpt = meann * (total_t1 / num_t1 + total_t2 / num_t2) + total_t1 / num_t1 + total_tr / converged;

    % Compute an estimate of the transition probability
    trans_prob = 1.0 -  exp(-1.0 / mfpt * tmax);
end