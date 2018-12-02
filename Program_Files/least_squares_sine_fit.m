function [amplitude,phase] = least_squares_sine_fit(omega,x,y);
% LEAST_SQUARES_SINE_FIT : function to perform a least squares sine fit.
% Given a vector x and corresponding vector y (both of which are either row
% or column vectors, and the x_i are not necessarily equally spaced) and
% given omega, this function applies a least squares sine fit to find the
% "best" sine function through the data (i.e. it finds best amplitude and
% phase):
%
% y = amplitude*sin(omega*x + phase)
%
% ******************************** Author *********************************
%
% Dr David Forehand
% The University of Edinburgh
%
% This version - 01/12/2018
%
% *************************************************************************

sin_omega_x = sin(omega*x);
cos_omega_x = cos(omega*x);
A = sum(sin_omega_x.*cos_omega_x);
B = sum(cos_omega_x.*cos_omega_x);
C = sum(sin_omega_x.*sin_omega_x);
D = sum(y.*sin_omega_x);
E = sum(y.*cos_omega_x);

phi_1 = atan((C*E-A*D)/(B*D-A*E));

if phi_1 < 0
    phi_2 = phi_1 + pi;
else
    phi_2 = phi_1 - pi;
end

cos_phi_1 = cos(phi_1);
sin_phi_1 = sin(phi_1);

denominator_1 = C*(cos_phi_1^2) + 2*A*sin_phi_1*cos_phi_1 ...
    + B*(sin_phi_1^2);
denominator_2 = A*((cos_phi_1^2)-(sin_phi_1^2)) ...
    + (B-C)*sin_phi_1*cos_phi_1;

if abs(denominator_1) > abs(denominator_2)
    a_1 = (D*cos_phi_1+E*sin_phi_1)/denominator_1;
else
    a_1 = (E*cos_phi_1-D*sin_phi_1)/denominator_2;
end

if a_1 > 0
    amplitude = a_1;
    phase = phi_1;
else
    amplitude = -a_1;
    phase = phi_2;
end