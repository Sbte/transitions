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

Nmfpt = 10000;
Ndirect = 10000;
Ntams = 10000;

% Generic part
phi = @(x) dist_fun(x, zA, zB);
VxxEv = eig(Vxx(zC));

trans_prob_list1 = {};
trans_prob_list2 = {};
trans_prob_list3 = {};
trans_prob_list4 = {};
trans_prob_list5 = {};
trans_prob_list6 = {};

mfpt_list1 = [];
mfpt_list3 = [];
mfpt_list5 = [];

data_list1 = [];
data_list2 = {};
data_list3 = [];
data_list4 = {};
data_list5 = [];
data_list6 = {};

error_list1 = {};
error_list2 = {};
error_list3 = {};
error_list4 = {};
error_list5 = {};
error_list6 = {};

normalized_error_list1 = {};
normalized_error_list2 = {};
normalized_error_list3 = {};
normalized_error_list4 = {};
normalized_error_list5 = {};
normalized_error_list6 = {};

for Bi=1:length(Brange)
    B = Brange(Bi);

    trans_prob_list1{Bi} = [];
    trans_prob_list2{Bi} = [];
    trans_prob_list3{Bi} = [];
    trans_prob_list4{Bi} = [];
    trans_prob_list5{Bi} = [];
    trans_prob_list6{Bi} = [];

    data_list2{Bi} = [];
    data_list4{Bi} = [];
    data_list6{Bi} = [];

    error_list1{Bi} = {};
    error_list2{Bi} = {};
    error_list3{Bi} = {};
    error_list4{Bi} = {};
    error_list5{Bi} = {};
    error_list6{Bi} = {};

    normalized_error_list1{Bi} = [];
    normalized_error_list2{Bi} = [];
    normalized_error_list3{Bi} = [];
    normalized_error_list4{Bi} = [];
    normalized_error_list5{Bi} = [];
    normalized_error_list6{Bi} = [];

    mfptt = 2*pi / -min(VxxEv) * sqrt(abs(det(Vxx(zC)))/det(Vxx(zA))) * exp((V(zC) - V(zA)) / (sigma^2/2));
    mfpt_list1 = [mfpt_list1, mfptt];
    data_list1 = [data_list1, 0];

    mfpt = 0;
    data = 0;
    if mfptt < 1e5
        [data, trans_prob, mfpt] = make_samples(...
            @transitions_mfpt, samples, F, B, z0, phi, dt, 1, Nmfpt, rho);
    end
    mfpt_list3 = [mfpt_list3, mfpt];
    data_list3 = [data_list3, data];

    [data, trans_prob, mfpt] = make_samples(...
        @transitions_ams, samples, F, B, z0, phi, dt, 1, Nmfpt, rho);
    mfpt_list5 = [mfpt_list5, mfpt];
    data_list5 = [data_list5, data];

    Ti = 0;
    for tmax=Trange
        fprintf('T=%d\n', tmax);
        Ti = Ti + 1;

        trans_prob_list1{Bi} = [trans_prob_list1{Bi}, 1 - exp(-1 / mfpt_list1(Bi) * tmax)];

        [data, trans_prob] = make_samples(...
            @transitions_direct, samples, F, B, z0, phi, dt, tmax, Ndirect, rho);
        trans_prob_list2{Bi} = [trans_prob_list2{Bi}, trans_prob];
        data_list2{Bi} = [data_list2{Bi}, data];
        error_list2{Bi}{Ti} = [trans_prob - data.Q1, data.Q3 - trans_prob];
        normalized_error_list2{Bi} = [normalized_error_list2{Bi}, sqrt(dt) * data.normalized_error];

        trans_prob = 1 - exp(-1 / mfpt_list3(Bi) * tmax);
        data = data_list3(Bi);
        trans_prob_list3{Bi} = [trans_prob_list3{Bi}, trans_prob];
        error_list3{Bi}{Ti} = [trans_prob - data.Q1 / data.mu * trans_prob,
                            data.Q3 / data.mu * trans_prob - trans_prob];
        normalized_error_list3{Bi} = [normalized_error_list3{Bi}, sqrt(dt) * data.normalized_error / data.mu * trans_prob];

        [data, trans_prob] = make_samples(...
            @transitions_gpa, samples, F, B, z0, phi, dt, tmax, Ndirect, rho);
        trans_prob_list4{Bi} = [trans_prob_list4{Bi}, trans_prob];
        data_list4{Bi} = [data_list4{Bi}, data];
        error_list4{Bi}{Ti} = [trans_prob - data.Q1, data.Q3 - trans_prob];
        normalized_error_list4{Bi} = [normalized_error_list4{Bi}, sqrt(dt) * data.normalized_error];

        trans_prob = 1 - exp(-1 / mfpt_list5(Bi) * tmax);
        data = data_list5(Bi);
        trans_prob_list5{Bi} = [trans_prob_list5{Bi}, trans_prob];
        error_list5{Bi}{Ti} = [trans_prob - data.Q1 / data.mu * trans_prob,
                            data.Q3 / data.mu * trans_prob - trans_prob];
        normalized_error_list5{Bi} = [normalized_error_list5{Bi}, sqrt(dt) * data.normalized_error / data.mu * trans_prob];

        [data, trans_prob] = make_samples(...
            @transitions_tams, samples, F, B, z0, phi, dt, tmax, Nmfpt, Ntams, rho);
        trans_prob_list6{Bi} = [trans_prob_list6{Bi}, trans_prob];
        data_list6{Bi} = [data_list6{Bi}, data];
        error_list6{Bi}{Ti} = [trans_prob - data.Q1, data.Q3 - trans_prob];
        normalized_error_list6{Bi} = [normalized_error_list6{Bi}, sqrt(dt) * data.normalized_error];
    end
end

generate_plots