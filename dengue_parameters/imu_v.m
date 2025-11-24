function x = imu_v(T, param)
    x = Quadratic_function(T, param.imv_c, param.imv_tmin, param.imv_tmax);
    x = max(x,20);
end
