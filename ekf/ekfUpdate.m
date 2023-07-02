function [estState,res] = ekfUpdate(p, cpt, dt)

%-------------------%
% Initialize
estState.isb_dict = dictionary;
estState.isb_dict(p.glo.sys_num) = NaN;
estState.isb_dict(p.gal.sys_num) = NaN;
estState.isb_dict(p.bds.sys_num) = NaN;
x_minus = p.state0;
[H_isb,x_isb] = formIsbStatesAndH(cpt.num_sv);
if length(x_isb) + 5 ~= length(p.state0)
    error('current No. of sys does not match the previous epoch');
end
%------------------%
y = cpt.corr_range;
num = length(y); % The number of measurement
H = zeros(num,5);
Range = zeros(num,1);
r = zeros(num,1);
off = zeros(num,1);
if p.post_mode == 1
    if p.IGS_enable == 1
        s_pos_ecef = cpt.s_pos_prc;
    else
        s_pos_ecef = cpt.s_pos_ecef;
    end
else
    s_pos_ecef = cpt.s_pos_ecef;
end
for j=1:num
    Range(j)=norm(s_pos_ecef(:,j)-x_minus(1:3));
    V= (x_minus(1:3)-s_pos_ecef(:,j))'/Range(j)+...
        [-s_pos_ecef(2,j)*p.omge/p.c s_pos_ecef(1,j)*p.omge/p.c 0];
    H(j,:)=[V, 1, 0];
    r(j) = Range(j)+sagnac(p,s_pos_ecef(:,j),x_minus(1:3));
    if ~isempty(H_isb)
        ind = find(H_isb(j,:)==1);
        if ~isempty(ind)
            off(j) = x_minus(5+ind);
        end
    end
end
H_os = [H,H_isb];
R = constructMeasNoise(p.c, cpt.elev, dt);
% measurement residual
res = y - r - x_minus(4)-off;
% Innovation (or residual) covariance
S = H_os * p.state_cov * H_os' + R;
% Near-optimal Kalman Gain
Kk = p.state_cov * H_os' * S^(-1);
delta_x = Kk * res;
x_plus = x_minus + delta_x;
I = eye(size(p.state_cov));
cov_plus = (I - Kk*H_os)*p.state_cov*(I - Kk*H_os)' + Kk*R*Kk';

p.state0 = x_plus;
p.state_cov = cov_plus;

estState.pos = x_plus(1:3);
estState.clock_bias = x_plus(4);
estState.clock_drift = x_plus(5);

if ~isempty(H_isb)
    isb_est = x_plus(6:end);
    j = 1;
    for i = 2:length(cpt.num_sv)
        if cpt.num_sv(i) ~= 0
            estState.isb_dict(i) = isb_est(j);
            j = j+1;
        end
    end
end


