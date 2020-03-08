function [] = w_test (div)

if isrow(div)
  div = div';
end

w_half = [0; cumsum(div)];
w_trapz = 0.5*(w_half(1:end-1) + w_half(2:end));
w_simpson = compute_w(div, 1);
%w = w - w(end);

p = (0:length(div)-1) + 0.5;
figure;
subplot(2,1,1)
plot(p, div);
title('div')
subplot(2,1,2)
plot(p, w_trapz, p, w_simpson, p, min(w_half(2:end), w_half(1:end-1)));
legend('trapz','simpson','min')


end

function [w, w_surf] = compute_w(div, dp)

w = zeros(length(div),1);
w(1) = dp*div(1)/2;
w(2) = w(1) + (5*div(1)/12 + 2*div(2)/3 - div(3)/12) * dp;
for k = 3:length(div)
  w(k) = w(k-2) + 1*dp*(div(k-2) + 4*div(k-1) + div(k))/3;
end

w_surf = w(end) + dp*div(end)/2;
%w = w - w_surf;

end
