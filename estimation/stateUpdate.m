function [p,estState,res] = stateUpdate(p, cpt, dt)

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
% y - f(x0) = H (x - x0);
zk = res + H_os * x_minus;
switch p.est_mode
    case p.ekf_est
        [x_plus, cov_plus] = ekfUpdate(x_minus, p.state_cov, res, H_os, R);
    case p.map_est
        [x_plus, cov_plus] = mapUpdate(ones(length(y),1), x_minus, p.state_cov, res, H_os, R);
    case p.raps_est   
        lla = ecef2lla(x_minus(1:3)', 'WGS84');
        % Horizontal: alpha = 1, beta = 0.15;
        % Vertical: alpha = 3, beta = 0.2.
        cov_spec_ecef = specInNedToEcef(lla', diag([p.raps.hor_cov_spec; ...
            p.raps.hor_cov_spec; p.raps.ver_cov_spec]));
        p_clk = diag([p.raps.clk_cov_spec;p.raps.dclk_cov_spec;...
            p.raps.isb_cov_spec*ones(length(x_minus)-5, 1)]);
        p_u = [cov_spec_ecef, zeros(3, length(x_minus)-3);
               zeros(length(x_minus)-3, 3), p_clk];
        J_l = p_u^(-1);
        p.raps_spec_xyz = diag(J_l(1:3, 1:3));
        [x_plus,cov_plus,b,augcost,J_out,exitflag,num_iter,comp_t] = ...
            mapRiskAverse(zk,H_os,p.state_cov,R,diag(diag(J_l)),x_minus);
        p.raps_J = J_out;
    otherwise
        error('Incorrect state estimation mode configuration');
end

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


