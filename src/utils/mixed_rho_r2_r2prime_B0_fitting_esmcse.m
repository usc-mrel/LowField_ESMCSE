function out = mixed_rho_r2_r2prime_B0_fitting_esmcse(x,y,t_shift,t_se)
% -----------------------
% Input:
%   x(1)   :   real part of rho
%   x(2)   :   imaginary part of rho
%   x(3)   :   R2  [1/msec]
%   x(4)   :   R2' [1/msec]
%   x(5)   :   off-reasoance [Hz]
%   y      :   obeserved data (2xN,1)
%   t_shift:   vector of shift time (N,1)
%   t_SE   :   vecotr of effective spin-echo time (N,1)

real_rho = x(1);
imag_rho = x(2);
R2       = x(3);
R2prime  = x(4);
delta_f  = x(5)/1000 ;

Ax = zeros(length(y),1);

for nt = 1:length(t_shift) % construct vector of signal at echa time from signal equations
    tshift = t_shift(nt);
    tse    = t_se(nt);
     if  tshift < 0
        Ax(nt) = ( real_rho * cos(2*pi*delta_f*tshift) + imag_rho * sin(2*pi*delta_f*tshift) ) * exp(-R2*(tse+tshift)) * exp(R2prime*tshift);
        Ax(nt+length(t_shift)) = ( imag_rho * cos(2*pi*delta_f*tshift) - real_rho *sin(2*pi*delta_f*tshift) ) * exp(-R2*(tse+tshift)) * exp(R2prime*tshift);
    elseif tshift >= 0
        Ax(nt) =  ( real_rho * cos(2*pi*delta_f*tshift) + imag_rho * sin(2*pi*delta_f*tshift) ) * exp(-R2*(tse+tshift)) * exp(-R2prime*tshift);
        Ax(nt+length(t_shift)) =  ( imag_rho * cos(2*pi*delta_f*tshift) - real_rho *sin(2*pi*delta_f*tshift) ) * exp(-R2*(tse+tshift)) * exp(-R2prime*tshift);
    end
end

out = Ax - y;
end
