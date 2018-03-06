F = @(x) x-x.^3;
B = sqrt(0.1);

dt = 0.005;
tmax = 2;
rho = 0.05;

tic
trans_prob = transitions_direct(F, B, dt, tmax, 200000, rho)
toc
% [trans_prob2, mfpt2] = transitions_mfpt(F, B, dt, tmax, 2000, rho)
% trans_prob3 = transitions_gpa(F, B, dt, tmax, 1000, rho)