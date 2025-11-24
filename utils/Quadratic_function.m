function x = Quadratic_function(v, c, vmin, vmax)
    x = zeros(size(v));

    x(v < vmin ) = 0.1;
    x(v > vmax ) = 0.1;
    between_indexing = v >= vmin & v < vmax;
    vt = v(between_indexing);
    x(between_indexing) =  c*(vt-vmax).*(vt-vmin);
    
end
% FROM: https://github.com/jms5151/SEI-SEIR_Arboviruses/blob/master/models/SEI-SEIR_model_with_trait_variation.R
%     mu_th <- function(temp, hum, timestep){
%   if (hum <= 1){
%     inverted_quadratic(temp, mu_th_c, mu_th_T0, mu_th_Tm, timestep)+(1-(0.01256 + 2.00893*hum))*0.005
%   } else {
%     inverted_quadratic(temp, mu_th_c, mu_th_T0, mu_th_Tm, timestep)+(1-(1.2248 + 0.2673*hum))*0.01
%   }
% }

% inverted_quadratic <- function(x, c, T0, Tm, timestep){
%   if((x < T0) | (x > Tm))
%     1.0/timestep
%   else
%     1.0/(c*(x-T0)*(x-Tm))
% }