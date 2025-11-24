function CI_plot(mean_val, lower_bound, upper_bound,line_col,area_col,face_alpha)
if nargin <= 3
    line_col = 'k';
    area_col = [0.9 0.9 0.9];
    face_alpha = 1;
end
time_index = 1:numel(mean_val);
fillyy(time_index, lower_bound, upper_bound, area_col,face_alpha);
plot(mean_val,line_col);
end
