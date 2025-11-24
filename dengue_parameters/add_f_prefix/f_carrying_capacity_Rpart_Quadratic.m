function x = carrying_capacity_Rpart_Quadratic(R,Rmin,Rmax,z)
c = -59.90e-4;
% z in [0.015, 0.025]
x = Quadratic_function(R,c,Rmin,Rmax)*z;
end
