function R = constructMeasNoise(p, elev_rad, sys_mark, dt)
% DGNSS measurement noise model from RTKLIB (rtkpos.c -> varerr)
% 2.0*(a*a+b*b/sinel/sinel+c*c)+d*d;
% a and b are both = code/carrier ratio (default: 300) * carrier-phase error STD (default: 0.003m)
% c is base length STD whose default is 0 for < 10km.
% d = speed of light * click stability (default: 5E-12) * delta time
R = zeros(length(elev_rad), length(elev_rad));
sigma_a = 300*0.003;
sigma_b = 400*0.003;
sigma_a_glo = 600*0.003;
sigma_b_glo = 1000*0.003;
d = p.c * 5e-12 * dt;
for i = 1:length(elev_rad)
    if sys_mark(i) == p.glo.sys_num
        R(i,i) = 4 * (sigma_a_glo^2+sigma_b_glo^2 / (sin(elev_rad(i)))^2) + d^2;
        continue;
    end
    % R(i,i) = 2 * (sigma_a^2 + sigma_b^2) + d^2;
    R(i,i) = 4 * (sigma_a^2 + sigma_b^2 / (sin(elev_rad(i)))^2 ) + d^2;
end
end