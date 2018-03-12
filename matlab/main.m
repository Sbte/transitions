F = @(x) x-x.^3;
B = sqrt(0.08);

dt = 0.001;
tmax = 2;
rho = 0.05;

samples = 10;

trans_prob_list = [];
trans_prob_list2 = [];
trans_prob_list3 = [];
trans_prob_list4 = [];

trans_time = transitions_fastest(F, B, dt, 10000, rho)

mfpt2 = 0;
for i=1:samples
    [trans_prob, mfpt] = transitions_mfpt(F, B, 0.01, tmax, 1000, rho);
    mfpt2 = mfpt2 + mfpt / samples;
end

mfpt4 = 0;
for i=1:samples
    [trans_prob, mfpt] = transitions_ams(F, B, 0.01, tmax, 1000, rho);
    mfpt4 = mfpt4 + mfpt / samples;
end

for tmax=0.5:0.1:2
    fprintf('T=%d\n', tmax);

    trans_prob2 = 0;
    for i=1:samples
        trans_prob = transitions_direct(F, B, dt, tmax, 200000, rho);
        trans_prob2 = trans_prob2 + trans_prob / samples;
    end
    trans_prob_list = [trans_prob_list, trans_prob2];

    trans_prob_list2 = [trans_prob_list2, 1 - exp(-1 / mfpt2 * tmax)];

    trans_prob2 = 0;
    for i=1:samples
        trans_prob = transitions_gpa(F, B, dt, tmax, 100000, rho);
        trans_prob2 = trans_prob2 + trans_prob / samples;
    end
    trans_prob_list3 = [trans_prob_list3, trans_prob2];

    trans_prob_list4 = [trans_prob_list4, 1 - exp(-1 / mfpt4 * tmax)];
end

plot(trans_prob_list);
hold on
plot(trans_prob_list2);
plot(trans_prob_list3);
plot(trans_prob_list4);
hold off
legend('Direct', 'MFPT', 'GPA', 'AMS')
xlabel('Tmax')
ylabel('Transition probability')