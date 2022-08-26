function out = rho_R2StarMEGRE_fitting(x,y,TEs)
% x(1): real part of rho
% x(2): imaginary part of rho
% x(3): R2 (1/msec)

real_rho = x(1);
%imag_rho = x(2);
R2      = x(2);

R_Ax = real_rho .* exp(-TEs*R2);
%I_Ax = imag_rho .* exp(-TEs*R2);

%Ax = [R_Ax; I_Ax]; % A_wave * u_wave 

Ax = R_Ax;

out = Ax - y;
end

