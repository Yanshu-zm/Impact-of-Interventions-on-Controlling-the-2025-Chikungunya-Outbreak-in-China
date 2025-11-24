function CI_plot(mean_val, lower_bound, upper_bound)
time_index = 1:numel(mean_val);
fillyy(time_index, lower_bound, upper_bound, [0.9 0.9 0.9]);
plot(mean_val,'k');
end
