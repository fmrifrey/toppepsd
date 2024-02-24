function [y,h,x_ramp,x_plat] = opttrap(A,ymax,dydxmax,dx)
% A = area under trapezoid y(x)
% ymax = maximum allowed amplitude of trapezoid
% dydxmax = maximum first derivative of trapezoid
%
% y = optimal trapezoid waveform

% Calculate trapezoid height
h = sqrt(A*dydxmax);
if (h > ymax), h = ymax; end

% Calculate ramp and plateau time
x_ramp = h / dydxmax;
x_plat = A / h - x_ramp;

% Round to nearest sampling interval & correct height based on error
x_ramp = dx * ceil(x_ramp / dx);
x_plat = dx * ceil(x_plat / dx);
h = A / (x_ramp + x_plat);

% Formulate the trapezoid
x = 0:dx:2*x_ramp+x_plat;
x1 = x_ramp;
x2 = x_ramp + x_plat;
x3 = 2*x_ramp + x_plat;
y = ...
    (x <= 0) .* 0 + ...
    (x < x1) .* h/x_ramp .* x + ...
    (x >= x1) .* (x < x2) .* h + ...
    (x >= x2) .* (x < x3) .* h .* (1 - 1/x_ramp * (x - x2)) + ...
    (x >= x3) .* 0;

end