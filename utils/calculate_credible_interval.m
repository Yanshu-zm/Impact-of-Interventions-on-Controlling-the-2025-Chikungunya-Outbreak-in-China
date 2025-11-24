function [lower_bound, upper_bound] = calculate_credible_interval(samples, confidence_level)
    % Calculate the lower and upper bounds based on the specified confidence level
    lower_percentile = (100 - confidence_level) / 2;
    upper_percentile = 100 - lower_percentile;
    
    % Calculate the lower and upper bounds for the credible interval
    lower_bound = prctile(samples, lower_percentile);
    upper_bound = prctile(samples, upper_percentile);
    
    % Display the credible interval
    %fprintf('The %d%% credible interval is [%.4f, %.4f]\n', confidence_level, lower_bound, upper_bound);
end
