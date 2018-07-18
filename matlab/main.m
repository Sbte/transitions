% 1D problem
% F = @(x) x-x.^3;
% sigma = sqrt(0.1);
% B = sigma;

% z0 = -1;
% zA = z0;
% zB = 1;
% zC = 0;

% V = @(x) (1/4)*x.^4-(1/2)*x.^2;
% Vx = @(x) x.^3-x;
% Vxx = @(x) 3*x.^2-1;

% 2D problem
F = @(x) [x(1,:)-x(1,:).^3;-2*x(2,:)];
sigma = sqrt(0.1);
B = sigma;

z0 = [-1;0];
zA = z0;
zB = [1;0];
zC = [0;0];

V = @(x) (1/4)*x(1).^4-(1/2)*x(1).^2+x(2)^2;
Vx = @(x) [x(1).^3-x(1);2*x(2)];
Vxx = @(x) [3*x(1).^2-1,0;0,2];

% General parameters
dt = 0.01;
rho = 0.05;
Trange = 1:50;
Brange = [B];

samples = 100;

Nmfpt = 1000;
Ndirect = 1000;
Ntams = 1000;

% Generic part
phi = @(x) dist_fun(x, zA, zB);
VxxEv = eig(Vxx(zC));

trans_prob_list = {};
trans_prob_list2 = {};
trans_prob_list3 = {};
trans_prob_list4 = {};
trans_prob_list5 = {};
trans_prob_list6 = {};

error_list = {};
error_list2 = {};
error_list3 = {};
error_list4 = {};
error_list5 = {};
error_list6 = {};

mfpt_list = [];
mfpt_list2 = [];
mfpt_list4 = [];

mfpt_error_list = [];
mfpt_error_list2 = [];
mfpt_error_list4 = [];

for Bi=1:length(Brange)
    B = Brange(Bi);

    trans_prob_list{Bi} = [];
    trans_prob_list2{Bi} = [];
    trans_prob_list3{Bi} = [];
    trans_prob_list4{Bi} = [];
    trans_prob_list5{Bi} = [];
    trans_prob_list6{Bi} = [];

    error_list{Bi} = [];
    error_list2{Bi} = [];
    error_list3{Bi} = [];
    error_list4{Bi} = [];
    error_list5{Bi} = [];
    error_list6{Bi} = [];

    mfptt = 2*pi / -min(VxxEv) * sqrt(abs(det(Vxx(zC)))/det(Vxx(zA))) * exp((V(zC) - V(zA)) / (sigma^2/2));
    mfpt_list = [mfpt_list, mfptt];
    mfpt_error_list = [mfpt_error_list, 0];

    mfpt = 0;
    error = 0;
    if mfptt < 1e5
        [error, trans_prob, mfpt] = make_samples(@transitions_mfpt, samples, F, B, z0, phi, dt, 1, Nmfpt, rho);
    end
    mfpt_list2 = [mfpt_list2, mfpt];
    mfpt_error_list2 = [mfpt_error_list2, error];

    [error, trans_prob, mfpt] = make_samples(@transitions_ams, samples, F, B, z0, phi, dt, 1, Nmfpt, rho);
    mfpt_list4 = [mfpt_list4, mfpt];
    mfpt_error_list4 = [mfpt_error_list4, error];

    for tmax=Trange
        fprintf('T=%d\n', tmax);

        [error, trans_prob] = make_samples(@transitions_direct, samples, F, B, z0, phi, dt, tmax, Ndirect, rho);
        trans_prob_list{Bi} = [trans_prob_list{Bi}, trans_prob];
        error_list{Bi} = [error_list{Bi}, error * trans_prob];

        trans_prob_list2{Bi} = [trans_prob_list2{Bi}, 1 - exp(-1 / mfpt_list2(Bi) * tmax)];
        error_list2{Bi} = [error_list2{Bi}, mfpt_error_list2(Bi) * trans_prob_list2{Bi}(end)];

        [error, trans_prob] = make_samples(@transitions_gpa, samples, F, B, z0, phi, dt, tmax, Ndirect, rho);
        trans_prob_list3{Bi} = [trans_prob_list3{Bi}, trans_prob];
        error_list3{Bi} = [error_list3{Bi}, error * trans_prob];

        trans_prob_list4{Bi} = [trans_prob_list4{Bi}, 1 - exp(-1 / mfpt_list4(Bi) * tmax)];
        error_list4{Bi} = [error_list4{Bi}, mfpt_error_list4(Bi) * trans_prob_list4{Bi}(end)];

        trans_prob_list5{Bi} = [trans_prob_list5{Bi}, 1 - exp(-1 / mfpt_list(Bi) * tmax)];
        error_list5{Bi} = [error_list5{Bi}, mfpt_error_list(Bi) * trans_prob_list5{Bi}(end)];

        [error, trans_prob] = make_samples(@transitions_tams, samples, F, B, z0, phi, dt, tmax, Nmfpt, Ntams, rho);
        trans_prob_list6{Bi} = [trans_prob_list6{Bi}, trans_prob];
        error_list6{Bi} = [error_list6{Bi}, error * trans_prob];
    end
end

plot(cell2mat(trans_prob_list5))
hold on
errorbar(cell2mat(trans_prob_list), cell2mat(error_list))
errorbar(cell2mat(trans_prob_list2), cell2mat(error_list2))
errorbar(cell2mat(trans_prob_list3), cell2mat(error_list3))
errorbar(cell2mat(trans_prob_list4), cell2mat(error_list4))
errorbar(cell2mat(trans_prob_list6), cell2mat(error_list6))
legend('Theory', 'Direct', 'MFPT', 'GPA', 'AMS', 'TAMS', 'Location', 'northwest')
xlabel('Tmax')
ylabel('Transition probability')
hold off