plot(Trange, cell2mat(trans_prob_list1))
hold on
errorbar(Trange, cell2mat(trans_prob_list2), getQ(error_list2, 1), getQ(error_list2, 2))
errorbar(Trange, cell2mat(trans_prob_list3), getQ(error_list3, 1), getQ(error_list3, 2))
errorbar(Trange, cell2mat(trans_prob_list4), getQ(error_list4, 1), getQ(error_list4, 2))
errorbar(Trange, cell2mat(trans_prob_list5), getQ(error_list5, 1), getQ(error_list5, 2))
errorbar(Trange, cell2mat(trans_prob_list6), getQ(error_list6, 1), getQ(error_list6, 2))
legend('Theory', 'Direct', 'MFPT', 'GPA', 'AMS', 'TAMS', 'Location', 'northwest')
xlabel('Tmax')
ylabel('Transition probability')
hold off

function out = getQ(z, num)
out = cell2mat(cellfun(@(x) cellfun(@(y) y(num), x), z, 'UniformOutput', false));
end