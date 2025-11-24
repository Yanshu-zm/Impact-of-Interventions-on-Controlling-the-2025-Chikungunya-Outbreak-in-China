function [xticks_positions, xticks_labels] = get_month_xticks_week()
    weeks_in_month  = [4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5];
    months = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
    
    xticks_positions = cumsum(weeks_in_month ); % Calculate the positions for month boundaries
    xticks_labels = cell(1, 12);
    for i = 1:12
        xticks_labels{i} = sprintf('%s %d', months{i});
    end
end
