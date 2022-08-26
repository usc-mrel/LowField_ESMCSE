function out = rho_R2bSSFP_fitting(x,TEs,y)
% x(1): real part of rho
% x(2): imaginary part of rho
% x(3): R2 (1/msec)

real_rho = x(1);
% imag_rho = x(2);
T2      = x(2);

R_Ax = real_rho .* exp(-TEs*T2);
% I_Ax = imag_rho .* exp(-TEs*T2);

% Ax = [R_Ax; I_Ax]; % A_wave * u_wave 

% Ax = cat(1, R_Ax, I_Ax); % A_wave * u_wave 
Ax = R_Ax;

out = Ax - y;
end

