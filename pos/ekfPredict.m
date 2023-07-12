function [prior_state, prior_cov] = ekfPredict(p, state, cov, dt)

% x = [x,y,z,clk,clk_drift,isb..]
phi = diag([ones(1, 5), ...
    exp(p.ekf_para.u_isb * dt) * ones(1, length(state) - 5)]);
phi(4, 5) = dt; % clk(k) = clk(k-1) + clk_drift * dt

Qc = p.ekf_para.q_clkDrift * [dt^3 / 3, dt^2 / 2;
                   dt^2 / 2, dt];
Qp = p.ekf_para.q_pos * eye(3,3);
Q_pc = [Qp, zeros(3, 2);
        zeros(2, 3), Qc];
Lbias = length(state) - 5;
if Lbias ~= 0
    Qd = [Q_pc, zeros(5, Lbias);
          zeros(Lbias, 5), p.ekf_para.q_isb * eye(Lbias, Lbias)];
else
    Qd = Q_pc;
end

prior_state = phi * state;
prior_cov = phi * cov * phi' + Qd;