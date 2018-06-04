F = @(x) x-x.^3;
B = sqrt(0.1);

dt = 0.01;
rho = 0.05;

samples = 1;

V = @(x) (1/4)*x.^4-(1/2)*x.^2;
Vx = @(x) x.^3-x;
Vxx = @(x) 3*x.^2-1;

trans_prob_list = {};
trans_prob_list2 = {};
trans_prob_list3 = {};
trans_prob_list4 = {};
trans_prob_list5 = {};
trans_prob_list6 = {};

mfpt_list = [];
mfpt_list2 = [];
mfpt_list4 = [];

N1=2000;
N2=10000;
N3=10000;

for Bi=1:15
    B = sqrt(0.8^(Bi-1))

    mfptt = 2*pi / -Vxx(0) * 1 * exp((V(0) - V(-1)) / (B^2/2))
    mfpt_list = [mfpt_list, mfptt];

    mfpt2 = 0;
    if mfptt < 1e5
        for i=1:samples
            [trans_prob, mfpt] = transitions_mfpt(F, B, dt, 1, N1, rho);
            mfpt2 = mfpt2 + mfpt / samples;
        end
    end
    mfpt_list2 = [mfpt_list2, mfpt2];

    mfpt4 = 0;
    for i=1:samples
        [trans_prob, mfpt] = transitions_ams(F, B, dt, 1, N1, rho);
        mfpt4 = mfpt4 + mfpt / samples;
    end
    mfpt_list4 = [mfpt_list4, mfpt4];

     trans_prob_list{Bi} = [];
    trans_prob_list2{Bi} = [];
    trans_prob_list3{Bi} = [];
    trans_prob_list4{Bi} = [];
    trans_prob_list5{Bi} = [];
    trans_prob_list6{Bi} = [];

    for tmax=1:10
        fprintf('T=%d\n', tmax);

        trans_prob2 = 0;
        for i=1:samples
            trans_prob = transitions_direct(F, B, dt, tmax, N2, rho);
            trans_prob2 = trans_prob2 + trans_prob / samples;
        end
        trans_prob_list{Bi} = [trans_prob_list{Bi}, trans_prob2];

        trans_prob_list2{Bi} = [trans_prob_list2{Bi}, 1 - exp(-1 / mfpt2 * tmax)];

        trans_prob2 = 0;
        for i=1:samples
            trans_prob = transitions_gpa(F, B, dt, tmax, N2, rho);
            trans_prob2 = trans_prob2 + trans_prob / samples;
        end
        trans_prob_list3{Bi} = [trans_prob_list3{Bi}, trans_prob2];

        trans_prob_list4{Bi} = [trans_prob_list4{Bi}, 1 - exp(-1 / mfpt4 * tmax)];

        trans_prob_list5{Bi} = [trans_prob_list5{Bi}, 1 - exp(-1 / mfptt * tmax)];

        trans_prob2 = 0;
        for i=1:samples
            trans_prob = transitions_tams(F, B, dt, tmax, N1, N3, rho);
            trans_prob2 = trans_prob2 + trans_prob / samples;
        end
        trans_prob_list6{Bi} = [trans_prob_list6{Bi}, trans_prob2];
    end

    if Bi > 1
        x = sqrt(0.8.^((1:length(mfpt_list))-1));
        t = 1:10

        figure(1)
        plot(x, mfpt_list);
        hold on
        plot(x, mfpt_list2);
        plot(x, mfpt_list4);
        hold off
        legend('Theory', 'MFPT', 'AMS')
        set(gca, 'YScale', 'log')
        drawnow;

        cols = colormap(lines);
        figure(2)
        surf(x, t, cell2mat(trans_prob_list5')','FaceAlpha', 0.8, 'FaceColor', cols(1,:));
        hold on
        surf(x, t, cell2mat(trans_prob_list')' ,'FaceAlpha', 0.8, 'FaceColor', cols(2,:));
        surf(x, t, cell2mat(trans_prob_list2')','FaceAlpha', 0.8, 'FaceColor', cols(3,:));
        surf(x, t, cell2mat(trans_prob_list3')','FaceAlpha', 0.8, 'FaceColor', cols(4,:));
        surf(x, t, cell2mat(trans_prob_list4')','FaceAlpha', 0.8, 'FaceColor', cols(5,:));
        surf(x, t, cell2mat(trans_prob_list6')','FaceAlpha', 0.8, 'FaceColor', cols(6,:));
        hold off
        legend('Theory', 'Direct', 'MFPT', 'GPA', 'AMS', 'TAMS')
        xlabel('Noise')
        ylabel('Tmax')
        zlabel('Transition probability')
        set(gca, 'ZScale', 'log')
        drawnow;
    end
end