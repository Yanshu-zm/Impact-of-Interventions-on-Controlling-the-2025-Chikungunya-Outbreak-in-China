function x = carrying_capacity_Tpart(T,param)
    T0 = 29.0;
    EA = 0.05;
    theta = 8.617e-5;
    Nv_max = 2; % ratio or abandunce?
    
    rate = EFD(T0,param) .* pEA(T0,param) .* MDR(T0,param) .* imu_v(T0,param);
    exponetial_term = (- EA*(T-T0).^2)./(theta*(T+273).*(T0+273));
    
    x = Nv_max*(rate - 1./imu_v(T0,param))./rate .*exp(exponetial_term);
end
