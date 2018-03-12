function [trans_prob, mfpt] = transitions_ams(F, B, dt, tmax, N, rho)
% Compute the mean first passage time and transition probability with AMS

    N2 = 1000000;
    M = 1e7;

    experiments = {};

    for i=1:N
        t = 0;
        z = -1;
        for j=1:M
            t = t + dt;
            dW = randn(1,1) * sqrt(dt);
            z = z + dt * F(z) + B * dW;
            dist = dist_fun(z);
            if dist > 0.1
                clear experiment;
                experiment.start_time = t;
                experiment.x = [z];
                experiment.t = [0];
                experiment.d = [dist];
                experiment.return_time = 0;
                experiment.max_dist = dist;
                experiments{i} = experiment;
                break;
            end
        end
    end

    for i=1:N
        t = 0;
        z = experiments{i}.x(end);
        for j=1:M
            t = t + dt;
            dW = randn(1,1) * sqrt(dt);
            z = z + dt * F(z) + B * dW;
            dist = dist_fun(z);
            if dist > experiments{i}.max_dist
                experiments{i}.x = [experiments{i}.x, z];
                experiments{i}.t = [experiments{i}.t, t];
                experiments{i}.d = [experiments{i}.d, dist];
                experiments{i}.max_dist = dist;
            elseif dist < rho
                if experiments{i}.return_time == 0
                    experiments{i}.return_time = t;
                end
                break;
            end
            if dist > 1-rho
                break;
            end
        end
    end

    its = 0;
    for i=1:N2
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

        experiments{min_idx}.x = experiments{idx}.x(1:same_dist_idx);
        experiments{min_idx}.t = experiments{idx}.t(1:same_dist_idx);
        experiments{min_idx}.d = experiments{idx}.d(1:same_dist_idx);
        experiments{min_idx}.max_dist = experiments{idx}.d(same_dist_idx);
        t = experiments{min_idx}.t(end);
        z = experiments{min_idx}.x(end);
        for j=1:M
            t = t + dt;
            dW = randn(1,1) * sqrt(dt);
            z = z + dt * F(z) + B * dW;
            dist = dist_fun(z);
            if dist > experiments{min_idx}.max_dist
                experiments{min_idx}.x = [experiments{min_idx}.x, z];
                experiments{min_idx}.t = [experiments{min_idx}.t, t];
                experiments{min_idx}.d = [experiments{min_idx}.d, dist];
                experiments{min_idx}.max_dist = dist;
            elseif dist < rho
                break;
            end
            if dist > 1-rho
                break;
            end
        end
        its = its + 1;
    end

    total_tr = 0;
    total_t1 = 0;
    total_t2 = 0;
    num_t1 = N;
    num_t2 = 0;
    for i=1:N
        total_tr = total_tr + experiments{i}.t(end);
        total_t1 = total_t1 + experiments{i}.start_time;
        total_t2 = total_t2 + experiments{i}.return_time;
        if (experiments{i}.return_time > 0)
            num_t2 = num_t2 + 1;
        end
    end
    alpha = (1.0 - 1.0 / N) ^ its;

    meann = 1.0 / alpha - 1.0;
    mfpt = meann * (total_t1 / num_t1 + total_t2 / num_t2) + total_t1 / num_t1 + total_tr / N;

    trans_prob = 1.0 -  exp(-1.0 / mfpt * tmax);
end