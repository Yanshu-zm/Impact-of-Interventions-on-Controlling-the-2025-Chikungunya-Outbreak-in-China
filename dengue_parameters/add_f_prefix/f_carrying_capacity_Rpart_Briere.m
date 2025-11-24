function x = carrying_capacity_Rpart_Briere(R,Rmin,Rmax,z)
c = 0.79e-4;
% z in [0.15, 0.25]
x = Briere_function(R,c,Rmin,Rmax)*z;
end
