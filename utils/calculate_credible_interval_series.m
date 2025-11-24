function [lower_bound, upper_bound] = calculate_credible_interval_series(samples, confidence_level)
    [T, ~] = size(samples);
    lower_bound = zeros(T,1);
    upper_bound = zeros(T,1);

    for t = 1: T
        [l,u] = calculate_credible_interval(samples(t,:), confidence_level);
        lower_bound(t) = l;
        upper_bound(t) = u;
    end
end
