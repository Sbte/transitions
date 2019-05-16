c = lines(6);

figure(1)
plot(Trange, cell2mat(trans_prob_list1), 'Color', c(1,:), 'LineWidth', 0.7)
hold on

error_fill(Trange, cell2mat(trans_prob_list2), getQ(error_list2, 1), getQ(error_list2, 2), c(2,:))
error_fill(Trange, cell2mat(trans_prob_list3), getQ(error_list3, 1), getQ(error_list3, 2), c(3,:))
error_fill(Trange, cell2mat(trans_prob_list4), getQ(error_list4, 1), getQ(error_list4, 2), c(4,:))
error_fill(Trange, cell2mat(trans_prob_list5), getQ(error_list5, 1), getQ(error_list5, 2), c(5,:))
error_fill(Trange, cell2mat(trans_prob_list6), getQ(error_list6, 1), getQ(error_list6, 2), c(6,:))

plot(Trange, cell2mat(trans_prob_list2), 'Color', c(2,:), 'LineWidth', 0.7, 'LineStyle', '-', 'Marker', 'o')
plot(Trange, cell2mat(trans_prob_list3), 'Color', c(3,:), 'LineWidth', 0.7, 'LineStyle', '-', 'Marker', '*')
plot(Trange, cell2mat(trans_prob_list4), 'Color', c(4,:), 'LineWidth', 0.7, 'LineStyle', '-', 'Marker', '+')
plot(Trange, cell2mat(trans_prob_list5), 'Color', c(5,:), 'LineWidth', 0.7, 'LineStyle', '-', 'Marker', 'd')
plot(Trange, cell2mat(trans_prob_list6), 'Color', c(6,:), 'LineWidth', 0.7, 'LineStyle', '-', 'Marker', 'x')

legend('Theory', 'Direct', 'Direct MFPT', 'GPA', 'AMS', 'TAMS', 'Location', 'northwest')
xlabel('T_{max}')
ylabel('Transition probability')
hold off

figure(2)
loglog(Trange, cell2mat(normalized_error_list2), 'Color', c(2,:))
hold on
loglog(Trange, cell2mat(normalized_error_list3), 'Color', c(3,:))
loglog(Trange, cell2mat(normalized_error_list4), 'Color', c(4,:))
loglog(Trange, cell2mat(normalized_error_list5), 'Color', c(5,:))
loglog(Trange, cell2mat(normalized_error_list6), 'Color', c(6,:))
legend('Direct', 'Direct MFPT', 'GPA', 'AMS', 'TAMS', 'Location', 'west')
xlabel('T_{max}')
ylabel('Time-normalized relative error')
hold off

function out = getQ(z, num)
out = cell2mat(cellfun(@(x) cellfun(@(y) y(num), x), z, 'UniformOutput', false));
end