function z = compute_actual_z(p, T)
% p in Pa, T in K. Must have size(p,last) == size(T,last); Output z in
% meters

integrand = 287*T./(p*9.8);
geopot = cumtrapz(p, -integrand, ndims(T));
z = geopot - min(geopot(:));

end