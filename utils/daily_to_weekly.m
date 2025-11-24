function weekly_value = daily_to_weekly(daily_value, reduction,interval)
    if nargin < 2
        reduction = 'sum';
        interval = 7;
    end
    if nargin < 3
        interval = 7;
    end
    num_days = length(daily_value);
    num_weeks = ceil(num_days / interval);
    assert(mod(num_days,interval) == 0)
    weekly_matrix = reshape(daily_value, interval, num_weeks);
    if strcmp(reduction, 'mean')
        weekly_value = mean(weekly_matrix, 1);
    elseif strcmp(reduction, 'sum')
        weekly_value = sum(weekly_matrix, 1);
    elseif strcmp(reduction, 'outbreak_prob')
        weekly_value = outbreak_reduce(weekly_matrix);
    end
end

function out_prob = outbreak_reduce(weekly_matrix)
    % weekly_matrix = reshape(daily_value, interval, num_weeks);
    [interval, ~] = size(weekly_matrix);
    out_prob = 1 - weekly_matrix(1,:);
    for j = 2 : interval
        out_prob = out_prob .*  (1 - weekly_matrix(j,:));
    end
    out_prob = 1- out_prob;
end