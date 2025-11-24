function x = Briere_function(v, c, vmin, vmax)
    x = zeros(size(v));

    x(v < vmin ) = 0;
    x(v > vmax ) = 0;
    between_indexing = v >= vmin & v < vmax;
    vt = v(between_indexing);
    x(between_indexing) =  c*vt.*(vt-vmin).*sqrt(vmax-vt);
end
