function [p,estState,res] = stateUpdate(p, cpt, dt)

%-------------------%
% Initialize
estState.isb_dict = dictionary;
estState.isb_dict(p.glo.sys_num) = NaN;
estState.isb_dict(p.gal.sys_num) = NaN;
estState.isb_dict(p.bds.sys_num) = NaN;
x_minus = p.state0;
num_user_states = p.modeToNumUserStates(p.state_mode);
[H_isb,x_isb] = formIsbStatesAndH(p, cpt.num_sv, p.state_mode);
if length(x_isb) + num_user_states + 2 ~= length(p.state0)
    error('current No. of sys does not match the previous epoch');
end
%------------------%
y = cpt.corr_range;
num = length(y); % The number of measurement
H = zeros(num,num_user_states+2);
H(:,num_user_states+1) = ones(num,1);
Range = zeros(num,1);
r = zeros(num,1);
off = zeros(num,1);
if p.post_mode == p.mode_ppp
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
    H(j,1:3)=V;
    r(j) = Range(j)+sagnac(p,s_pos_ecef(:,j),x_minus(1:3));
    if ~isempty(H_isb)
        ind = find(H_isb(j,:)==1);
        if ~isempty(ind)
            off(j) = x_minus(num_user_states+2+ind);
        end
    end
end
H_os = [H,H_isb];
R = constructMeasNoise(p, cpt.elev, cpt.svprn_mark, dt);
% measurement residual
res = y - r - x_minus(num_user_states+1)-off;
[x_minus, p.state_cov, flag] = checkClockReset(p, x_minus, p.state_cov, ...
    num_user_states+1, res, cpt);
if flag == true
    res = y - r - x_minus(num_user_states+1)-off;
end
% y - f(x0) = H (x - x0);
zk = res + H_os * x_minus;
switch p.est_mode
    case p.ekf_est
        [x_plus, cov_plus] = ekfUpdate(x_minus, p.state_cov, res, H_os, R);
    case p.map_est
        [x_plus, cov_plus] = mapUpdate(ones(length(y),1), x_minus, p.state_cov, res, H_os, R);
    otherwise
        error('Incorrect state estimation mode configuration');
end

p.state0 = x_plus;
p.state_cov = cov_plus;

estState.pos = x_plus(1:3);
estState.clock_bias = x_plus(num_user_states+1);
estState.clock_drift = x_plus(num_user_states+2);

if ~isempty(H_isb)
    isb_est = x_plus(num_user_states+2+1:end);
    j = 1;
    for i = 2:length(cpt.num_sv)
        if cpt.num_sv(i) ~= 0
            estState.isb_dict(i) = isb_est(j);
            j = j+1;
        end
    end
end


