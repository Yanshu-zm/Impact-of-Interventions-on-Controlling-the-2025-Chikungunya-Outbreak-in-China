function CIs = possiCI(lambda, alpha)
    %alpha = 0.05; % Significance level for 95% credible interval
    lower_bound = icdf('Poisson', alpha/2, lambda);
    upper_bound = icdf('Poisson', 1 - alpha/2, lambda);

    lower_bound = reshape(lower_bound, [1, numel(lower_bound)]);
    upper_bound = reshape(upper_bound, [1, numel(upper_bound)]);
    CIs = [lower_bound;upper_bound];
end

