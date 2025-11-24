function weekly_value = daily_to_weekly_2D(daily_value, reduction)
    if nargin < 2
        reduction = 'sum';
    end
    num_days = size(daily_value,1);
    ndims = size(daily_value,2);
    num_weeks = ceil(num_days / 7);
    assert(mod(num_days,7) == 0)
    weekly_matrix = reshape(daily_value, 7, num_weeks,ndims);
    if strcmp(reduction, 'mean')
        weekly_value = squeeze(mean(weekly_matrix, 1));
    elseif strcmp(reduction, 'sum')
        weekly_value = squeeze(sum(weekly_matrix, 1));
    end
end
