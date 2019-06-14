function z = compute_standard_height(p)

i_trop = p >= 226.321;
i_strat = ~i_trop;

z_trop = 288.15 * (1 - (p/1013.25).^(6.5E-3*287.06/9.81)) / 6.5; % km
z_strat = 11 - 1E-3*287.06*216.65/9.81 * log(p/226.321);

z = zeros(size(p));
z(i_trop) = z_trop(i_trop);
z(i_strat) = z_strat(i_strat);

end