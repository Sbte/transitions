F = @(x) x-x.^3;
B = sqrt(0.1);

dt = 0.005;
tmax = 2;
rho = 0.05;

samples = 10;

trans_prob_list = [];
trans_prob_list2 = [];
trans_prob_list3 = [];

mfpt2 = 0;
for i=1:samples
    [trans_prob, mfpt] = transitions_mfpt(F, B, 0.02, tmax, 2000, rho);
    mfpt2 = mfpt2 + mfpt / samples;
end

for tmax=1:20
    fprintf('T=%d\n', tmax);

    trans_prob2 = 0;
    for i=1:samples
        trans_prob = transitions_direct(F, B, dt, tmax, 20000, rho);
        trans_prob2 = trans_prob2 + trans_prob / samples;
    end
    trans_prob_list = [trans_prob_list, trans_prob2];

    trans_prob_list2 = [trans_prob_list2, 1 - exp(-1 / mfpt2 * tmax)];

    trans_prob2 = 0;
    for i=1:samples
        trans_prob = transitions_gpa(F, B, 0.001, tmax, 10000, rho);
        trans_prob2 = trans_prob2 + trans_prob / samples;
    end
    trans_prob_list3 = [trans_prob_list3, trans_prob2];
end

plot(trans_prob_list);
hold on
plot(trans_prob_list2);
plot(trans_prob_list3);
hold off
legend('Direct', 'MFPT', 'GPA')
xlabel('Tmax')
ylabel('Transition probability')