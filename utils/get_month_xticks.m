function [xticks_positions, xticks_labels] = get_month_xticks()
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    months = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
    
    xticks_positions = cumsum(days_in_month); % Calculate the positions for month boundaries
    xticks_labels = cell(1, 12);
    for i = 1:12
        xticks_labels{i} = sprintf('%s %d', months{i});
    end
end
