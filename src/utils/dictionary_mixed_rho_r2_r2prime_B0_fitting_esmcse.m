function out = dictionary_mixed_rho_r2_r2prime_B0_fitting_esmcse(x,y,t_shift,t_se)
% -----------------------
% Input:
%   x(1)   :   S
%   x(2)   :   R2  [kHz]
%   x(3)   :   R2' [Khz]
%   x(4)   :   off-reasoance [Hz]
%   y      :   obeserved data (N,1)
%   t_shift:   vector of shift time (N,1)
%   t_SE   :   vecotr of effective spin-echo time (N,1)

abs_rho = x(1);
R2       = x(2);
R2prime  = x(3);
delta_f  = x(4) ;

Ax = zeros(length(y),1);

for nt = 1:length(t_shift) % construct vector of signal at echa time from signal equations
    tshift = t_shift(nt);
    tse    = t_se(nt);
     if  tshift < 0
        Ax(nt) = abs_rho * exp(-2*pi*delta_f*tshift) * exp(-R2*(tse+tshift)) * exp(R2prime*tshift);
    elseif tshift >= 0
        Ax(nt) = abs_rho * exp(-2*pi*delta_f*tshift) * exp(-R2*(tse+tshift)) * exp(-R2prime*tshift);
    end
end

out = Ax - y;
end
