function [streamfun_max] = computeStreamfunMax (Q0, L, H, N, r, eps)
% [Q0] = m/s^3, H is the geometric height

H = N*H./sqrt(r*eps);
streamfun_max = 2*Q0.*L.*H./(3*N^2 * sqrt(2*L.^2+3*L.*H+H.^2));

end
