function fine_daily_mean = weekly_to_daily_mean(weekly_mean)
    length_v = size(weekly_mean,1)*7;
    mid = 3;
    fine_daily_mean = interp1(mid+1:7:mid+length_v-1, ...
                              weekly_mean, 1:length_v, 'spline'); 
                              % Interpolation method can be adjusted
end
